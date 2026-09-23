#!/usr/bin/perl -w
#######################################################
## expCor.cgi -- heatmap of similarities of experiments
##
## Copyright (c) 2021 University of California
##
## Authors:
## Morgan Price
#######################################################
#
# Required parameters:
# orgId, the organism
#
# Optional parameters:
# e -- the experiment ids to show
# add -- search for matching experiments to add

use strict;
use CGI qw(:standard Vars);
use CGI::Carp qw(warningsToBrowser fatalsToBrowser);
use DBI;
use Time::HiRes qw(gettimeofday);

use lib "../lib";
use Utils;
use FEBA_Utils qw{ReadTable};

sub CorToCell($$$); # exp1, exp2, and correlation value => table element (including td tag)
my $maxExp = 50;

my $cgi=CGI->new;
my $orgId = $cgi->param('orgId');
die "Must specify orgId" unless defined $orgId && $orgId ne "";
my @expNames = $cgi->param('e');

my $dbh = Utils::get_dbh();
my $orginfo = Utils::orginfo($dbh);
die "Unknown organism" unless exists $orginfo->{$orgId};

splice(@expNames, $maxExp) if @expNames > $maxExp;
my %exps = (); # expName to experiment object
my @exps = ();
foreach my $expName (@expNames) {
  next if exists $exps{$expName};
  my $exp = $dbh->selectrow_hashref("SELECT * from Experiment WHERE orgId = ? AND expName = ?",
                                          {}, $orgId, $expName)
    || die "Unknown expName";
  $exps{$expName} = $exp;
  push @exps, $exp;
}

my @errors = ();
my $add = $cgi->param('add') || "";
$add =~ s/[ \t]+$//;
$add =~ s/^[ \t]+//;
if ($add) {
  my $exps = Utils::matching_exps($dbh,$orgId,$add);
  if (@$exps == 0) {
    push @errors, qq{No experiments matching "$add"};
  } else {
    my $tooMany = 0;
    my $nIgnore = 0;
    foreach my $exp (@$exps) {
      my $expName = $exp->{expName};
      if (exists $exps{$expName}) {
        $nIgnore++;
      } elsif (scalar(@exps) < $maxExp) {
        push @exps, $exp;
        $exps{$expName} = $exp;
      } else {
        $tooMany = 1;
      }
    }
    push @errors, qq{Not all matching experiments are shown: the limit is $maxExp.}
      if $tooMany;
    push @errors, qq{$nIgnore matching experiment(s) are already shown.}
      if $nIgnore > 0;
  }
}

print
  header,
  Utils::start_page("Experiment Correlations for $orginfo->{$orgId}{genome}"),
  '<div id="ntcontent">',
  h2("Experiment Correlations for",
     a({ -href => "org.cgi?orgId=$orgId" }, $orginfo->{$orgId}{genome}));

foreach my $error (@errors) {
  print h3($error);
}

my @hidden = ();
push @hidden, hidden('orgId', $orgId);
foreach my $exp (@exps) {
  push @hidden, hidden(-name => 'e', -default => $exp->{expName}, -override => 1);
}
print
  # corner compare box
  qq[<div style="position: relative;"><div class="floatbox">],
  start_form('-name' => 'select', -method => 'GET', -action => 'expCor.cgi'),
  join("", @hidden),
  "Add: ",
  textfield( -name => "add", -default => "", -override => 1, -size => 10, -maxLength => 500 ),
    "<button type='submit'>Go</button>",
  end_form,
  qq{</div></div>};

if (@exps > 0) {
  print h3("Linear (Pearson) correlations for", scalar(@exps), "experiments");
  print  '</div>'; # end ntcontent
  my %expFit = (); # expName => locusId => fit
  my %locusSeen = ();
  foreach my $expName (keys %exps) {
    my $fit = $dbh->selectall_arrayref("SELECT locusId,fit FROM 'FitByExp_${orgId}' WHERE expName = ?",
                                       {}, $expName);
    die "No fitness data for $expName" unless @$fit > 0;
    foreach my $row (@$fit) {
      my ($locusId, $fit) = @$row;
      $expFit{$expName}{$locusId} = $fit;
      $locusSeen{$locusId} = 1;
    }
  }
  die "No fitness data" unless keys %locusSeen > 0;

  # Compute pairwise correlations using R
  my $tmpDir = Utils::tmp_dir();
  my $procId = $$;
  my $timestamp = int (gettimeofday * 1000);
  my $prefix = $procId . $timestamp;
  my $dataFile = "$tmpDir/expCor.$prefix.tsv";
  open(my $fhData, ">", $dataFile) || die "Cannot write to $dataFile";
  my @expNameSorted = sort keys %exps;
  print $fhData join("\t", @expNameSorted)."\n";
  foreach my $locusId (sort keys %locusSeen) {
    my @out = ();
    foreach my $expName (@expNameSorted) {
      push @out, exists $expFit{$expName}{$locusId} ? $expFit{$expName}{$locusId} : "";
    }
    print $fhData join("\t", @out)."\n";
  }
  close($fhData) || die "Error writing to $dataFile";
  my $corFile = "$tmpDir/expCor.$prefix.cor";
  my $Rcmd = qq{ write.table(cor(read.delim("$dataFile"), use="pairwise"), "$corFile", sep="\\t", row.names=FALSE, quote=FALSE) };
  system("Rscript", "-e", $Rcmd) == 0 || die "Rscript failed on $Rcmd -- $!";
  my @cors = ReadTable($corFile, \@expNameSorted);
  die "R failed" unless scalar(@cors) == scalar(@exps);
  unlink($dataFile);
  unlink($corFile);

  my %cors = (); # exp name to exp name to similarity
  foreach my $i (0..(scalar(@expNameSorted)-1)) {
    my $expName1 = $expNameSorted[$i];
    my $cors = $cors[$i];
    foreach my $expName2 (@expNameSorted) {
      $cors{$expName1}{$expName2} = $cors->{$expName2};
      die "Self correlation for $expName1 is not 1"
        if $expName1 eq $expName2 && $cors->{$expName2} != 1;
    }
  }

  my @th = ();
  push @th, "&nbsp;";
  foreach my $exp (@exps) {
    my $desc = $exp->{expDesc};
    $desc =~ s/_/ /g;
    push @th, small(a({ -href => "exp.cgi?orgId=$orgId&expName=$exp->{expName}", -title => $exp->{expName} },
                      $desc));
  }
  my @trows = ();
  push @trows, td(\@th);
  my @remove = ();
  push @remove, "&nbsp;";
  my $rmmark = "&#10799;"; # Unicode Character 'VECTOR OR CROSS PRODUCT' (U+2A2F)
  foreach my $exp (@exps) {
    my $URL = "expCor.cgi?orgId=$orgId";
    foreach my $exp2 (@exps) {
      $URL .= "&e=$exp2->{expName}" unless $exp2->{expName} eq $exp->{expName};
    }
    push @remove, a({ -href => $URL, -title => "remove $exp->{expName}"}, $rmmark);
  }
  push @trows, Tr(td({-align => 'CENTER', -valign => 'TOP'}, \@remove));
  foreach my $i (0..(scalar(@exps)-1)) {
    my $exp = $exps[$i];
    my $cors = $cors{ $exp->{expName} };
    my @trow = td(a({ -href => "exp.cgi?orgId=$orgId&expName=$exp->{expName}", -title => $exp->{expName} },
                    $exp->{expDesc} ));
    push @trow, map CorToCell($exp, $_, $cors->{$_->{expName}}), @exps;
    push @trows, Tr({ -align => 'CENTER', -valign => 'TOP'}, @trow);
  }
  print table( { cellspacing => 0, cellpadding => 1 }, @trows);
} else {
  print p("No experiments have been selected. Please search for experiments to add.");
  print  '</div>'; # end ntcontent
}


sub CorToCell($$$) {
  my ($exp1, $exp2, $cor) = @_;
  my $orgId = $exp1->{orgId} || die;
  my $expName1 = $exp1->{expName};
  my $expName2 = $exp2->{expName};
  my $URL = $exp1->{expName} eq $exp2->{expName} ? "exp.cgi?orgId=$orgId&expName=$exp1->{expName}"
    : "compareExps.cgi?orgId=$orgId&expName1=$exp1->{expName}&expName2=$exp2->{expName}";
  my $shade = 0; # 0 for poorly correlated to 127 for perfectly correlated
  if ($cor ne "NA" && $cor ne "" && $cor > 0) {
    $shade = int(127 * $cor + 0.5);
  }
  # grey to green
  my $color = sprintf("#%2.2X%2.2X%2.2X", 127 - $shade, 128 + $shade, 127 - $shade);
  return td({ -bgcolor => $color },
            a({ -href => $URL, -title => "$expName1 vs. $expName2: $exp1->{expDesc} vs. $exp2->{expDesc}",
                -style => "color: black;" },
              sprintf("%.2f", $cor)));
}
