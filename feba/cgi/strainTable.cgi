#!/usr/bin/perl -w
#######################################################
## strainTable.cgi
##
## Copyright (c) 2015 University of California
##
## Authors:
## Morgan Price
#######################################################
#
# Required CGI parameters:
# orgId -- which organism
# and either
#	scaffoldId, begin, end -- the range of interest (limited to 50 kb)
# or
#	locusId -- which gene to show (by default, range is exactly the width of the gene)
# or
#	scaffoldId, object (see below). [Object can also be specified together with begin/end]
# Optional CGI parameters:
# expName -- which experiment(s) to show. (Can be more than one, or none.)
# addexp -- additional experiments (i.e. a setname or a condition)
# zoom -- in or out
# pan -- left or right
# object -- an additional object to show, specified by a colon-delimited set of arguments such as
#	b:1000:e:2000:s:1:n:name:d:popup
#	begin, end, and name must be specified, and begin should be less than end
#	object is shown as unstranded if strand is absent
#	use 1 for + strand, not "+", because of CGI encoding issues
#	d is shown as popup text and is optional
#	(unlike genomeBrowse.cgi, just one object is supported)
# format -- svg to download svg, tsv to download user-friendly table (default empty)

use strict;
use CGI qw(-nosticky :standard Vars);
use CGI::Carp qw(warningsToBrowser fatalsToBrowser);
use DBI;
use List::Util qw{sum min max};
use POSIX qw{floor ceil};

use lib "../lib";
use Utils;
use StrainFitness;
use svgPlot;

my $cgi=CGI->new;
my $style = Utils::get_style();
my $dbh = Utils::get_dbh();
my $orginfo = Utils::orginfo($dbh);
my $orgId = $cgi->param('orgId') || die "No orgId found";
my $genome = $orginfo->{$orgId}{genome} || die "Invalid orgId $orgId";

my @expNames = $cgi->param('expName');
my $scaffoldId = $cgi->param('scaffoldId');
my $begin = $cgi->param('begin');
my $end = $cgi->param('end');
my $locusSpec = $cgi->param('locusId');
my $locusSpecShow;
my $expName = $cgi->param('expName') || "";
my $debug = $cgi->param('debug') || "";
my $help = $cgi->param('help') || "";
my $objspec = $cgi->param('object') || "";
my $format = $cgi->param('format') || ""; # default is html page

# Set up begin/end from locusSpec
if (defined $locusSpec && $locusSpec ne "") {
  my $sysName;
  ($scaffoldId,$begin,$end,$sysName) = $dbh->selectrow_array(
                                                             "SELECT scaffoldId,begin,end,sysName FROM Gene WHERE orgId = ? AND locusId = ?",
                                                             {
                                                             }, $orgId, $locusSpec);
  die "Unknown locus $locusSpec in $orgId" if !defined $end;
  $locusSpecShow = $sysName || $locusSpec;
  # my $widen = int(1 + 0.2 * ($end-$begin+1));
  my $widen = 1000;
  $begin -= $widen;
  $end += $widen;
} elsif (defined $scaffoldId && defined $begin && defined $end) {
  die "Invalid scaffold parameter" if $scaffoldId eq "";
  die "Invalid begin parameter" unless $begin =~ m/^-?\d+$/;
  die "Invalid end parameter" unless $end =~ m/^\d+$/;
}

# Add the user-defined object, if any, and maybe modify begin/end
my %object = ();
if ($objspec) {
  die "object is used by scaffold is not set" unless $scaffoldId ne "";
  my %keynames = ("b" => "begin", "e" => "end", "s" => "strand", "n" => "name", "d" => "desc");
  my @parts = split /:/, $objspec;
  while (@parts > 0) {
    my $key = shift @parts;
    die "Unknown key $key in object argument\n" unless exists $keynames{$key};
    my $value = shift @parts;
    die "Wrong number of fields in $objspec\n" unless defined $value;
    $object{$keynames{$key}} = $value;
  }
  $object{object} = 1 ;   # so that geneArrows() knows it is an object
  die "Invalid object $objspec\n"
    unless exists $object{begin} && exists $object{end} && exists $object{name};
  die "begin must be less than end" unless $object{begin} < $object{end};
  # To build the URL for getNtSeq.cgi, may need to reverse begin and end
  my ($beginL, $endL) = ($object{begin}, $object{end});
  ($beginL,$endL) = ($endL,$beginL)
    if exists $object{strand} && ($object{strand} eq "-" || $object{strand} eq "-1");
  # this link is currently always on the + strand, should fix
  $object{URL} = "getNtSeq.cgi?orgId=$orgId&scaffoldId=$scaffoldId&begin=$beginL&end=$endL";

  if (!defined $begin || !defined $end) {
    $begin = $object{begin} - 500;
    $end = $object{begin} + 500;
    if ($begin - $end + 1 < 4000) {
      my $center = int(($begin+$end)/2);
      $begin = $center - 2000;
      $end = $center + 2000;
    }
  }
}

# Handle zoom/pan options
my $zoom = $cgi->param('zoom') || "";
my $initwidth = $end - $begin + 1;
if ($zoom eq "in") {
  $begin += 0.2 * $initwidth;
  $end -= 0.2 * $initwidth;
} elsif ($zoom eq "out") {
  $begin -= 0.4 * $initwidth;
  $end += 0.4 * $initwidth;
}
my $pan = $cgi->param('pan') || "";
if ($pan eq "left") {
  $begin -= 0.4 * $initwidth;
  $end -= 0.4 * $initwidth;
} elsif ($pan eq "right") {
  $begin += 0.4 * $initwidth;
  $end += 0.4 * $initwidth;
}
$begin = int($begin);
$end = int($end);

$end = $begin + 1 if $begin eq $end;
unless ($begin < $end) {
  print $cgi->header;
  Utils::fail($cgi, "Invalid begin/end $begin $end")
  }
my $maxWidth = 40*1000;
if ($end - $begin >= $maxWidth) {
  print $cgi->header;
  Utils::fail($cgi, "Too wide, the limit is to show a range of " . Utils::commify($maxWidth) . " nucleotides");
}

# Add additional experiments?
my $addexp = $cgi->param('addexp');
if (defined $addexp && $addexp ne "") {
  my @expsNew = @{ Utils::matching_exps($dbh,$orgId,$addexp) };
  if (@expsNew == 0) {
    print header;
    Utils::fail($cgi, qq{No experiment matching "$addexp"});
  }
  # else
  my %oldExps = map { $_ => 1 } @expNames;
  push @expNames, grep {!exists $oldExps{$_} } map { $_->{expName} } @expsNew;
}

# Fetch genes, experiments, and strain fitness values
my $expinfo = Utils::expinfo($dbh,$orgId);
foreach my $expName (@expNames) {
  die "No such experiment: $expName" unless exists $expinfo->{$expName};
}

my $genes = $dbh->selectall_arrayref("SELECT * FROM Gene WHERE orgId = ? AND scaffoldId = ? AND Gene.end >= ? AND Gene.begin <= ? ORDER by Gene.begin",
                                     {
                                      Slice => {} }, $orgId, $scaffoldId, $begin, $end);
my %genesh = map {$_->{locusId} => $_} @$genes;

my $strains = StrainFitness::GetStrainFitness("../cgi_data", $dbh, $orgId, $scaffoldId, $begin, $end);
# leave out gene if not a used strain
foreach my $row (@$strains) {
  $row->{locusId} = "" unless $row->{used} eq "TRUE";
}

if ($format eq "tsv") {
  # Table of strains and fitness values
  print "Content-Type:text/tab-separated-values\n";
  print "Content-Disposition: attachment; filename=strainTable_${orgId}.tsv\n\n";
  my @expShow = map { $_ . " " . $expinfo->{$_}{expDesc} } @expNames;
  print join("\t", qw{scaffoldId position strand locusId gene fraction}, @expShow)."\n";
  foreach my $row (@$strains) {
    my $locusId = $row->{locusId};
    my $displayName = $locusId eq "" ? ""
      : $genesh{$locusId}{gene} || $genesh{$locusId}{sysName} || $locusId;
    my $frac = "";
    $frac = sprintf("%.2f", ($row->{pos} - $genesh{$locusId}{begin})
                    / ( $genesh{$locusId}{end} -  $genesh{$locusId}{begin} + 1))
      if $locusId ne "";
    print join("\t", $scaffoldId, $row->{pos}, $row->{strand}, $locusId, $displayName,
               $frac,
               map { $row->{$_} } @expNames)."\n";
  }
  exit 0;
}

# Make the SVG
my $defs = "";
my $svg = "";
my $svgWidth = 800;
my @showObj = @$genes;
push @showObj, \%object if keys(%object) > 0;
@showObj = sort { $a->{begin} <=> $b->{begin} } @showObj;
my $plotTop = 0;
if (@showObj > 0) {
  ($defs, $svg, $svgWidth) = Utils::geneArrowsParts($dbh, \@showObj, $locusSpec, $begin, $end);
  $defs .= qq{<clipPath id="arrowsRegion"><rect x="0" y="-30" width="800" height="130" /></clipPath>};

  # add left and right borders to illustrate truncation
  foreach my $x (0, 800) {
    $svg .= qq{<line x1="$x" x2="$x" y1="-30" y2="100" stroke="black" />};
  }
  # for x, instead of 0 to 800, now at 80 to 880
  # for y, instead of -30 to 100, now at 0 to 130
  $svg = qq{<g transform="translate(80 30)" clip-path="url(#arrowsRegion)">$svg</g>};
  $plotTop = 130;
}
my $plotBottom = $plotTop;
if (@$strains > 0 && @expNames > 0) {
  # Plot the strain fitness values
  foreach my $row (@$strains) {
    my @values = map $row->{$_}, @expNames;
    $row->{avg} = sum(@values) / scalar(@values);
  }
  my @avgs = map $_->{avg}, @$strains;
  my @yRange = ( min(-1, min(@avgs) - 0.25),
                 max(1, max(@avgs) + 0.25) );
  my $plot = svgPlot->new('width' => 800+20*8, height => '350',
                          'lineSize' => 20,
                          'margin' => [3,4,0.5,4],
                          'xLeft' => 0,
                          'yTop' => $plotTop,
                          'xRange' => [$begin/1000,$end/1000],
                          'yRange' => \@yRange);
  $plotBottom = $plot->{yBottom};
  $svg .= $plot->axes();
  $svg .= $plot->marginText("Position (kb)", "bottom", 'lineAt' => 2);
  my @atX = ceil($begin/1000) .. floor($end/1000);
  if (scalar(@atX) < 2 || $atX[0] == $atX[1]) {
    @atX = map $_/10, ceil($begin/100) .. floor($end/100);
  }
  if (scalar(@atX) > 20) {
    @atX = grep $_ % 2 == 0, @atX;
  }
  $svg .= $plot->axisTicks("x", \@atX);

  $svg .= $plot->marginText("Strain fitness (log2 ratio)", "left", 'lineAt' => 1.75);
  my @atY = grep { $_ >= $yRange[0] && $_ <= $yRange[1] } (-20..20);
  $svg .= $plot->axisTicks("y", \@atY);
  my $y0 = $plot->convertY(0);
  $svg .= qq{<line x1="80" x2="880" y1="$y0" y2="$y0" stroke="grey" />};
  foreach my $row (@$strains) {
    my $locusId = $row->{locusId};
    my $fill = $locusId eq "" ? "grey" :
      ($row->{strand} eq "+" ? "green" : "red");
    my $label = sprintf("at %.3f kb on %s strand", $row->{pos}/1000, $row->{strand});
    if ($locusId ne "") {
      my $displayName = $genesh{$locusId}{gene} || $genesh{$locusId}{sysName} || $locusId;
      $label .= ", within $displayName";
    }
    $svg .= $plot->point($row->{pos}/1000, $row->{avg},
                         'size' => $locusId eq "" ? 2 : 3,
                         'symbol' => $row->{strand},
                         'fill' => $fill, 'color' => $fill,
                         'title' => $label);
  }
}

if ($format eq "svg") {
  die "No experiment specified\n" if scalar(@expNames) == 0;
  # Rewrite onclick to include the full URL
  my $base = url();
  $base =~ s/strainTable.cgi.*//;
  $svg =~ s!href=&apos;!href=&apos;$base!g;
  print
    "Content-Type:image/svg+xml\n",
    "Content-Disposition: attachment; filename=strainTable_${orgId}.svg\n\n",
    join("\n",
         qq{<?xml version="1.0"?>},
         qq{<svg xmlns="http://www.w3.org/2000/svg" xmlns:xlink="http://www.w3.org/1999/xlink" width="890" height="$plotBottom">},
         "<defs>",
         $defs,
         "</defs>",
         qq{<g transform="scale(1)">},
         $svg,
         "</g>",
         "</svg>"),
    "\n";
  exit(0);
}

# else make web page
print header;
my @expShow = ();
foreach my $expName (@expNames) {
  push @expShow, a({href => "exp.cgi?orgId=$orgId&expName=$expName", title => $expName}, $expinfo->{$expName}{expDesc});
}
my $begComma = Utils::commify($begin);
my $endComma = Utils::commify($end);
print
  Utils::start_page("Strain Fitness in $genome"),
  q{<div id="ntcontent">},
  h2("Strain Fitness in ",
     a({-href => "org.cgi?orgId=$orgId"}, "$genome"),
     defined $locusSpecShow ? "around " . a({-href => "singleFit.cgi?orgId=$orgId&locusId=$locusSpec"}, $locusSpecShow)
     : " at $scaffoldId: $begComma to $endComma");
if (@expShow == 0) {
  print p("No experiments selected");
} elsif (@expShow == 1) {
  print p("Experiment:", @expShow);
} else {
  print p("Averaging across", scalar(@expShow), "experiments:", @expShow);
}
if ($help) {
  print qq[<div class="helpbox">
        <b><u>About this page:</u></b><BR><ul>
        <li>View the fitness of genes in various strains under a condition.</li>
        <li> Data on + strands are colored green, data on - strands are colored red, and data not affiliated with a gene are colored gray. Hover on points or blue links to see more information. Use the buttons to navigate the genome.</li>
        <li>To get to this page, click on the colored fitness boxes across the site (mostly on any pages relating to genes).</li>
        <li>To see the average strain fitness across replicates, type an condition (like "cisplatin") into the "Add experiments" box and hit "Add."
        <ul><li>You can remove unwanted experiments by scrolling down to the per-strain table and clicking on "remove."</ul>
        </ul></div>];
}

print
  start_form(-name => 'input', -method => 'GET', -action => 'strainTable.cgi'),
  hidden( -name => 'orgId', -value => $orgId, -override => 1),
  hidden( -name => 'scaffoldId', -value => $scaffoldId, -override => 1),
  hidden( -name => 'begin', -value => $begin, -override => 1),
  hidden( -name => 'end', -value => $end, -override => 1),
  hidden( -name => 'object', -value => $objspec, -override => 1),
  join("\n", map { hidden( -name => 'expName', -value => $_, -override => 1) } @expNames),
  p({-class => "buttons", style=>"align:left; white-space:nowrap; line-height:40px;"}, "Add experiment(s): ",
    textfield(-name => 'addexp', -default => "", -override => 1, -size => 20, -maxLength => 100),
    submit('Add','Add') ),
  p({-class => "buttons", style=>"max-width:500px; line-height:40px; white-space:nowrap;"},
    "Zoom:", submit('zoom','in'), submit('zoom','out'), "\tPan:", submit('pan','left'), submit('pan','right')),
  end_form;
print 
  p(small(qq{Only strains with sufficient reads to estimate fitness are shown,
                   but the strain fitness values are still rather noisy and may be biased towards zero.
                   Strains near the edge of a gene are not shown as being associated with that
                   gene (they are in grey).
                   Strains in the central 10-90% of a gene are color coded by the insertion's strand.
                   Usually, "+" means that the selectable marker is encoded on the forward strand,
                   i.e., transcribed rightward.}
         ))
  if scalar(@expNames) > 0;

print p("No genes in range.") if @$genes == 0 && ! $objspec;

if (@$strains == 0) {
  print "No fitness data for strains within " . Utils::commify($begin) . " to " . Utils::commify($end) . "\n";
} else {
}

my $baseURL = "strainTable.cgi?orgId=$orgId&scaffoldId=$scaffoldId&begin=$begin&end=$end";
$baseURL .= "&object=$objspec" if $objspec;
my $svgURL = join("&", $baseURL, "format=svg", map { "expName=$_" } @expNames);
print div({-style => "float:right;"}, small(a({-href => $svgURL}, "download SVG"))), "\n";
print join("\n",
           "<center>",
           # Do not use class arrows, as that would limit the height
           # The gene region is at 80 to 880; use a slightly higher width so that
           # points at the far right do not get clipped
           qq{<svg width="890" height="$plotBottom">},
           "<defs>",
           $defs,
           "</defs>",
           qq{<g transform="scale(1.0)">},
           $svg,
           "</g>",
           "</svg>", "</center>")
  . "\n";

my @trows = (); # the HTML table
# header row
my @headings = qw{Position Strand Gene LocusTag};
push @headings, a({-title => "Fractional position within gene"}, "Fraction");
my @base_headings = @headings;
foreach my $expName (@expNames) {
  push @headings, a({-href => "exp.cgi?orgId=$orgId&expName=$expName", -title => $expName},
                    $expinfo->{$expName}{expDesc});
}
push @trows, Tr({-align => 'CENTER', -valign=>'TOP'}, th(\@headings));

my %locusIds = map { $_->{locusId} => 1 } @$strains;
my %genes = ();                 # locusId => row
foreach my $locusId (keys %locusIds) {
  next if $locusId eq "";
  my $gene = $dbh->selectrow_hashref("SELECT * FROM Gene WHERE orgId = ? AND locusId = ?",
                                     {
                                     }, $orgId, $locusId);
  die "Unknown locusId $locusId" unless exists $gene->{locusId};
  $genes{$locusId} = $gene;
}

Utils::endHtml($cgi) if scalar(@expNames) == 0 || scalar(@$strains) == 0;

# add row of links for removing items
my @remove_row = map { td("") } @base_headings; # 
foreach my $expName (@expNames) {
  my @otherExps = grep { $_ ne $expName } @expNames;
  my @otherExpSpec = map { "expName=$_" } @otherExps;
  push @remove_row, td( a({ -title => "remove $expName : $expinfo->{$expName}{expDesc}",
                            -href => join("&", $baseURL, @otherExpSpec) },
                          "remove") );
}
push @trows, Tr({-align => 'CENTER', -valign=>'TOP'}, @remove_row);

my @avgFits = ();
foreach my $row (@$strains) {
  my $locusId = $row->{locusId};
  my $locusShow = "";
  my $gene = undef;
  if ($locusId ne "") {
    $gene = $genes{$locusId};
    $locusShow = $gene->{sysName} || $gene->{locusId};
  }
  my @row = ( a({-title => "barcode $row->{barcode}"}, Utils::commify($row->{pos})),
              $row->{strand},
              $genesh{$row->{locusId}}{gene},
              $locusId eq "" ? "" : a({-title => $gene->{desc}, -href => "singleFit.cgi?orgId=$orgId&locusId=$locusId"},
                                      $locusShow),
              $locusId eq "" ? "" : sprintf("%.2f",
                                            ($row->{pos} - $gene->{begin}) / ($gene->{end} - $gene->{begin} + 1))
            );
  @row = map { td($_) } @row;
  my $totalFit = 0;             #gather the total for averaging
  my $ind = 0;                  #gather number of entries 
  foreach my $expName (@expNames) {
    my $fit = $row->{ $expName };
    $totalFit += $fit;
    $ind += 1;
    push @row, td( { -bgcolor => Utils::fitcolor($fit) }, sprintf("%+.1f",$fit));
  }
  # above we give up if no experiments are specified, so this no longer causes divide by 0
  push @avgFits, $totalFit/$ind;

  push @trows, Tr({-align => 'CENTER', -valign=>'TOP'}, @row);
}
my $tableURL = join("&", $baseURL, "format=tsv", map { "expName=$_" } @expNames);
print div({-style => "float:right;"},
          p({-style => "text-align:right; font-size: small;"},
            a({-href => $tableURL }, "download strain data"))),
  "\n";
print h3("Per-strain Table"), "\n";
print small(table({ cellspacing => 0, cellpadding => 3, }, @trows)), "\n";

print p("Or see this region's",
        a({ -href => "getNtSeq.cgi?orgId=$orgId&scaffoldId=$scaffoldId&begin=$begin&end=$end" },
          "nucleotide sequence"));

$dbh->disconnect();
Utils::endHtml($cgi);
