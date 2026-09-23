#!/usr/bin/perl -w
# cofitNegative.cgi -- given a query gene, find genes with the most negative cofitness
# Required parameters: orgId and locusId, for organism and which gene

use strict;
use CGI qw(:standard Vars);
use CGI::Carp qw(warningsToBrowser fatalsToBrowser);
use DBI;
use List::Util qw{sum};
use lib "../lib";
use Utils;
use IO::Handle; # for autoflush
sub rescale($); # fitness values to centered/rescaled with mean 0, sd 1

my $cgi=CGI->new;
print $cgi->header;
my $orgId = $cgi->param('orgId') || die "no species identifier\n";
my $locusId = $cgi->param('locusId') || die "no gene identifier\n";

Utils::fail($cgi, "$orgId is invalid. Please enter correct species name!") unless ($orgId =~ m/^[A-Za-z0-9_-]*$/);
Utils::fail($cgi, "$locusId is invalid. Please enter correct locusId!") unless ($locusId =~ m/^[A-Za-z0-9_.-]*$/);

my $dbh = Utils::get_dbh();
my $gene = $dbh->selectrow_hashref("SELECT * FROM Gene WHERE orgId=? AND locusId=?",
                                   undef, $orgId, $locusId);

if (!defined $gene) {
  my $xref = $dbh->selectall_arrayref("SELECT * from LocusXref WHERE orgId=? AND xrefId = ?",
                                      { Slice => {} }, $orgId, $locusId);
  if (@$xref == 1) {
    $locusId = $xref->[0]{locusId};
    $gene = $dbh->selectrow_hashref("SELECT * FROM Gene WHERE orgId=? AND locusId=?",
                                    undef, $orgId, $locusId);
  }
}
Utils::fail($cgi, "No locus $locusId in species $orgId") unless defined $gene;
my ($genus,$species,$strain) = $dbh->selectrow_array("SELECT genus,species,strain FROM Organism WHERE orgId=?", undef, $orgId);
Utils::fail($cgi, "No species information for $orgId") unless defined $species;

my $showId = $gene->{sysName} || $locusId;
my $showName = $gene->{gene} || $showId;
my $start = Utils::start_page("Negative cofitness for $showName ($genus $species $strain");
my $tabs = Utils::tabsGene($dbh,$cgi,$orgId,$locusId,0,$gene->{type},"cofit");

print
  $start, $tabs,
  h2("Negative cofitness for",
     a({-href => "singleFit.cgi?orgId=$orgId&locusId=$locusId"}, $showId),
     "from",
     $cgi->a({href => "org.cgi?orgId=$orgId"}, "$genus $species $strain")),
  p(Utils::gene_link($dbh, $gene, "lines"));

my $hasFit = $dbh->selectall_arrayref("SELECT * from GeneFitness WHERE orgId = ? AND locusId = ?",
                                     {}, $orgId, $locusId);
if (@$hasFit == 0) {
  print p("Sorry, there is no fitness data for",
          a({-href => "geneOverview.cgi?orgId=$orgId&gene=$locusId"}, $showId));
  Utils::endHtml($cgi);
  exit(0);
}
if (scalar(@$hasFit) < 10) {
  print p("Sorry, there are not enough fitness experiments for this organism");
  Utils::endHtml($cgi);
  exit(0);
}
#else
my $nExp = scalar(@$hasFit);
print p("Computing cofitness values with $nExp experiments...\n");
autoflush STDOUT 1; # show preliminary results
print "\n";
my $fits = $dbh->selectall_arrayref("SELECT locusId,expName,fit FROM GeneFitness WHERE orgId = ? ORDER BY locusId,expName",
                                    {}, $orgId);
my %fit = (); # locusId => list of fitness values
foreach my $row (@$fits) {
  my ($l, $expName, $fit) = @$row;
  push @{ $fit{$l} }, $fit;
}
die "No fitness data for $locusId" unless exists $fit{$locusId};
my $fitFocal = $fit{$locusId};
die unless $nExp == scalar(@$fitFocal);
my $fitFocalScaled = rescale($fitFocal);

my %cofit;
foreach my $l (sort keys %fit) {
  my $rescaled = rescale($fit{$l});
  my $cov;
  foreach my $i (0..($nExp-1)) {
    $cov += $rescaled->[$i] * $fitFocalScaled->[$i];
  }
  my $r = $cov/$nExp;
  $cofit{$l} = $r;
}

my @sortedLoci = sort { $cofit{$a} <=> $cofit{$b} } keys %fit;
my $maxHits = 20;
splice @sortedLoci, $maxHits if @sortedLoci > $maxHits;
@sortedLoci = grep $cofit{$_} < 0, @sortedLoci;
print p("Genes with the strongest negative cofitness:");

my @headRow = map { $cgi->th($cgi->b($_)) } qw{Rank Hit Name Description Cofit &nbsp;};
my @trows = ( $cgi->Tr(@headRow) );
my @colors = ('#FFFFDD', '#FFFFFF');
my $iRow = 0;

foreach my $hitId (@sortedLoci) {
  my $hitGene = $dbh->selectrow_hashref("SELECT * from Gene WHERE orgId=? AND locusId=?",
                                        {}, $orgId, $hitId);
  die $hitId unless defined $hitGene;
  my $rowcolor = $colors[ $iRow++ % scalar(@colors) ];
  my $cofit = sprintf("%.2f",$cofit{$hitId});
  my $url = "compareGenes.cgi?orgId=$orgId&locus1=$locusId&locus2=$hitId";
  push @trows,
    $cgi->Tr( {bgcolor => $rowcolor, align => 'left', valign => 'top' },
              $cgi->td($iRow),
              $cgi->td(Utils::gene_link($dbh, $hitGene, "name", "myFitShow.cgi")),
              $cgi->td($hitGene->{gene} || ""),
              $cgi->td(Utils::gene_link($dbh, $hitGene, "desc", "domains.cgi")),
              $cgi->td( $cgi->a({-href => $url},$cofit )),
              $cgi->td(checkbox('locusId',0,$hitId,''))
            );
}
print
  start_form(-name => 'input', -method => 'GET', -action => 'genesFit.cgi'),
  hidden('orgId', $orgId),
  hidden('locusId', $locusId),
  table( {cellpadding => 3, cellspacing => 0 }, @trows ),
	p(submit(-class=>"heatmap", -name=>"Heatmap of $showId with selected genes")),
  end_form;
print p("Or look for",
        a({-href => "cofit.cgi?orgId=$orgId&locusId=$locusId"}, "positive cofitness"));
Utils::endHtml($cgi);
exit(0);

# Not s.d. is dividing by n, not by n-1
sub rescale($) {
  my ($in) = @_;
  my $len = scalar(@$in);
  die if $len < 2;
  my $mean = sum(@$in) / $len;
  my @diff = map $_ - $mean, @$in;
  my $sumsq = sum(map $_ * $_, @diff);
  my $sd = sqrt($sumsq / $len);
  my @out = map $_ / $sd, @diff;
  return \@out;
}
