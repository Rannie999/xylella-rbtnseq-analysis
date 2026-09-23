#!/usr/bin/perl -w
#######################################################
## seqservice.cgi
##
## Copyright (c) 2015 University of California
##
## Authors:
## Morgan Price
#######################################################
#
# Required parameters: seq -- the protein sequence to look for matches to.
# Optional parameters: maxhits -- default is 20. Cannot be raised above 50.
# html -- return HTML code (with no wrappers) for the result, instead of a table
# wrap -- return HTML code with HTML tags (for testing). Must specify with html.
#
# Relies on a LAST database in ../cgi_data/ldb.*
#
# By default, returns a tab-delimited table of FastBLAST hits, or, a
# simple two-line table with Error and the error (or "no hits").
#
# For an example of a page that uses this service see
# ../images/fitblast_example.html
# ../images/fitblast.js

use strict;
use CGI qw(:standard Vars);
use CGI::Carp qw(warningsToBrowser fatalsToBrowser);
use DBI;

use lib "../lib";
use Utils;

my $cgi=CGI->new;
if ($cgi->param('wrap') && $cgi->param('html')) {
    print $cgi->header('text/html');
} else {
    print $cgi->header(-type => 'text/plain',  
                       # allow cross-site scripts to use this
                       -access_control_allow_origin => '*');
}

my $seq = $cgi->param('seq');
$seq = "" if !defined $seq;
$seq =~ s/[ \r\n]//g;
if (! $seq) {
    print "Error\nNo sequence specified\n";
    exit(0);
}
# Allow * as internal stop codon
unless ($seq =~ m/^[A-Z*]+$/) {
    print "Error\nIllegal sequence\n";
    exit(0);
}

if (length($seq) > 5000) {
  print "Error\nSorry, the sequence is too long\n";
  exit(0);
}

my $maxhits = $cgi->param('maxhits');
if (defined $maxhits && $maxhits > 0) {
    $maxhits = $maxhits + 0; # convert to integer
} else {
    $maxhits = 20;
}
$maxhits = 50 if $maxhits > 50;

my $lastal = "../bin/lastal";
die "No such executable: $lastal" unless -x $lastal;
my $ldb = "../cgi_data/ldb";
die "No such file: $ldb.suf" unless -e "$ldb.suf";
my $len = length($seq);
my $tmpPre = "/tmp/seqservice.$$";
open(my $fh, ">", "$tmpPre.faa") || die "Cannot write to $tmpPre.faa";
print $fh ">seq\n$seq\n";
close($fh) || die "Error writing to $tmpPre.faa";

my $cmd = "$lastal -f BlastTab -m $maxhits $ldb $tmpPre.faa > $tmpPre.hits";
system($cmd) ==0 || die "lastal failed\n$cmd\n$!";
my @hits = ();
open(my $fhResults, "<", "$tmpPre.hits") || die "No $tmpPre.hits file";
while(<$fhResults>) {
  next if m/^#/; # comment lines form lsatal
  chomp;
  my @F = split /\t/, $_;
  push @hits, \@F;
  last if @hits >= $maxhits;
}
close($fhResults) || die "Error reading $tmpPre.hits";

unlink("$tmpPre.faa");
unlink("$tmpPre.hits");

my $dbh = Utils::get_dbh() || die "Cannot access database";
my $base_url = url();
$base_url =~ s!/[a-zA-Z_.]+$!!;
my $blastURL = "$base_url/mySeqSearch.cgi?query=$seq";

my $html = $cgi->param('html');
# d3.tsv treats an empty table as an error, so instead, return an actual useful error code
if (@hits == 0 ) {
    if ($html) {
        print qq{No close hits (<A HREF="$blastURL" TITLE="View more hits">more</A>)\n};
    } else {
        print "Error\nNo hits\n";
    }
    exit(0);
}

print join("\t", qw{orgId organism locusId sysName name description identity coverage evalue bits minfit maxfit minT maxT maxcofit})."\n" unless $html;
my $orginfo = Utils::orginfo($dbh);
my @save = ();
foreach my $row (@hits) {
    my ($query,$locusspec,$identity,$alnlen,$mm,$gap,$qBeg,$qEnd,$lBeg,$lEnd,$evalue,$bits) = @$row;
    my ($orgId,$locusId) = split /:/, $locusspec;
    die "Invalid locus $locusspec" unless defined $orgId && defined $locusId;
    my ($sysName,$name,$desc) = $dbh->selectrow_array("SELECT sysName, gene, desc FROM Gene WHERE orgId = ? AND locusId = ?",
					  {}, $orgId, $locusId);
    $sysName = "" if !defined $sysName;
    $name = "" if !defined $name;
    my ($minFit,$maxFit,$minT,$maxT) = $dbh->selectrow_array(qq{ SELECT min(fit), max(fit), min(t), max(t)
                                                                 FROM GeneFitness WHERE orgId = ? AND locusId = ? ; },
							     {}, $orgId, $locusId);
    $minFit = "" unless defined $minFit;
    $maxFit = "" unless defined $maxFit;
    $minT = "" unless defined $minT;
    $maxT = "" unless defined $maxT;
    my $maxCofit = "";
    if (defined $minFit) {
	$maxCofit = $dbh->selectrow_array(qq{ SELECT cofit FROM Cofit WHERE orgId = ? AND locusId = ? AND rank = 1 LIMIT 1; },
					  {}, $orgId, $locusId);
        $maxCofit = "" if !defined $maxCofit; # possible if there is no cofitness for this bug
    }
    my $coverage = $alnlen/length($seq);
    if ($html) {
        push @save, [ $orgId, $orginfo->{$orgId}{genome}, $locusId, $sysName, $name, $desc,
                      $identity, $coverage,
                      $minFit, $maxFit, $minT, $maxT, $maxCofit ];
    } else {
        print join("\t",$orgId, $orginfo->{$orgId}{genome},
                   $locusId, $sysName, $name, $desc,
                   $identity,$coverage,$evalue,$bits,
                   $minFit,$maxFit,$minT,$maxT,$maxCofit)."\n";
    }
}

$dbh->disconnect();

if ($html) {
    # mimic fitblast_short() from fitblast.js
    my $mincoverage = 0.75; # hits below this coverage are ignored
    my $mincloseid = 80; # %identity to be considered a close hit
    my $minabsstrong = 2.0;
    # a hit is considered likely to be useful if either
    # cofitness is above threshold or (both absfit and absT are above thresholds)
    my $mincofit = 0.75;
    my $minabsfit = 1.0;
    my $minabsT = 4.0;
    my @pieces = (); # pieces of HTML, 1 for item shown so far

    foreach my $row (@save) {
        my ($orgId, $orgName, $locusId, $sysName, $name, $desc,
            $identity, $coverage, $minFit, $maxFit, $minT, $maxT, $maxCofit) = @$row;
        next unless $coverage >= $mincoverage;
        my $close = $identity >= $mincloseid;
        my $hascofit = $maxCofit >= $mincofit;
        my $hassig = ($minFit <= -$minabsfit && $minT <= $minabsT)
            || ($maxFit >= $minabsfit && $maxT >= $minabsT);
        my $hasstrong = $hassig && ($minFit <= -$minabsstrong || $maxFit >= $minabsstrong);
        my $useful = $hascofit || $hassig;
        if (($close && scalar(@pieces) == 0) || $useful) {
            my $sigstring = $hasstrong ? "strong phenotype" :
                ($hassig ? "has phenotype" : ($minFit ne "" ? "has data" : "no data"));
            if ($minFit ne "") {
                $sigstring = sprintf(qq{<A TITLE="fitness %+.1f to %+.1f" HREF="$base_url/singleFit.cgi?orgId=%s&locusId=%s">%s</A>},
                                     $minFit, $maxFit, $orgId, $locusId, $sigstring);
            }
            if ($hascofit) {
                $sigstring = sprintf(qq{%s, <A TITLE="top cofitness %.2f" HREF="$base_url/cofit.cgi?orgId=%s&locusId=%s">cofit</A>},
                                     $sigstring, $maxCofit, $orgId, $locusId);
            }
            push @pieces, sprintf(qq{%d%% id. to <A HREF="$base_url/singleFit.cgi?orgId=%s&locusId=%s">%s</A> from %s: %s\n},
                                  int($identity+0.5), $orgId, $locusId,
                                  $name eq "" ? ($sysName eq "" ? $locusId : $sysName) : $name,
                                  $orgName, $sigstring);
        }
        last if $useful;
    }
    print join("<BR>",@pieces) 
        . qq{ (<A HREF="$blastURL" TITLE="view all hits">more</A>)} . "\n";
}


exit(0);

