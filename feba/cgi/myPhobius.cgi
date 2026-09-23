#!/usr/bin/perl -w
# myPhobius.cgi -- given a single sequence, run Phobius to identify
# trans-membrane helices and signal peptides, and show the results.
#
# Required CGI parameters:
# seq -- the sequence
# name -- a name or description of the query (such as the definition line from a fasta file)

use strict;
use CGI qw{:standard Vars};
use CGI::Carp qw(warningsToBrowser fatalsToBrowser);
use HTML::Entities;
use Time::HiRes qw{gettimeofday};
use FindBin qw{$RealBin};
use lib "$RealBin/../lib";
use Utils qw{formatFASTA};

my $cgi=CGI->new;
my $seq = $cgi->param('seq') || "";
$seq =~ s/[ *\r\n]//g;
$seq = uc($seq);
die "No sequence specified" unless $seq =~ m/^[A-Z]+$/;
die "Maximum sequence length is 100K" if length($seq) > 100*1000;

my $name = $cgi->param('name');
die "No name specified" unless defined $name && $name ne "";
$name =~ s/[\r\n]+$//;
die "Invalid name" if $name =~ m/[\r\n]/;

my $procId = $$;
my $timestamp = int (gettimeofday * 1000);
my $tmpDir = 
my $tmpPre = Utils::tmp_dir() . "/$procId$timestamp";

open(my $fhFaa, ">", "$tmpPre.faa")
  || die "Cannot write to $tmpPre.faa";
print $fhFaa Utils::formatFASTA($name, $seq);
close($fhFaa) || die "Error writing to $tmpPre.faa";

my $nameSafe = HTML::Entities::encode($name);
print $cgi->header;
print start_html( -title => "Phobius results for $nameSafe"),
  "\n";

system("../bin/myPhobius.pl", "-in", "$tmpPre.faa", "-out", $tmpPre) == 0
  || die "myPhobius.pl failed on $tmpPre.faa -- $!";

system("cat", "$tmpPre.svg");
print "\n<PRE>\n\nTable:\n";
system("cat $tmpPre.tsv");
print "\n</PRE>\n";

print "\n<PRE>\nInput:\n";
system("cat $tmpPre.faa");
print "\n</PRE>\n";

foreach my $suffix (qw{faa plp R tsv svg}) {
  unlink("$tmpPre.$suffix");
}

print p("Or try the official",
        a({ -href => "https://phobius.sbc.su.se/" }, "Phobius web server,"),
        "which has a different display");

print p("Reference:",
        a( { -href => "http://www.ncbi.nlm.nih.gov/pmc/articles/pmc1933244/" },
           "Advantages of combined transmembrane topology and signal peptide prediction - the Phobius web server"));

print $cgi->end_html;
