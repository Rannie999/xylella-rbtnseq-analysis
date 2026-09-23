#!/usr/bin/perl -w

#######################################################
## domains.cgi -- sequence analysis of a protein
##
## Copyright (c) 2015 University of California
##
## Authors:
## Victoria Lo and Morgan Price
#######################################################
#
# Required CGI parameters:
# orgId -- which organism
# locusId (or gene) -- which locus

use strict;
use CGI qw(:standard Vars start_ul);
use CGI::Carp qw(warningsToBrowser fatalsToBrowser);
use DBI;
use IO::Handle; # for autoflush
use Bio::SeqIO;
use LWP::Simple qw{get};

use lib "../lib";
use Utils;

my $cgi=CGI->new;

my $orgId = $cgi->param('orgId') || die "No orgId parameter";
my $locusId = $cgi->param('locusId') || $cgi->param('gene') || die "No locusId or gene parameter";

my $dbh = Utils::get_dbh();
my $orginfo = Utils::orginfo($dbh);
die "Unknown organism" unless $orgId eq "" || exists $orginfo->{$orgId};
my $gene = $dbh->selectrow_hashref("SELECT * FROM Gene WHERE orgId=? AND locusId=?",
				   {}, $orgId, $locusId);
if (!defined $gene) {
  my $xref = $dbh->selectall_arrayref("SELECT * from LocusXref WHERE orgId=? AND xrefId = ?",
                                      { Slice => {} }, $orgId, $locusId);
  if (@$xref == 1) {
    $locusId = $xref->[0]{locusId};
    $gene = $dbh->selectrow_hashref("SELECT * FROM Gene WHERE orgId=? AND locusId=?",
                                    undef, $orgId, $locusId);
  }
}
die "Unknown gene" unless defined $gene->{locusId};
my $sys = $gene->{sysName} || $gene->{locusId};

if ($gene->{type} != 1) {
    # this can be reached by some links from gene descriptions
    # redirect to the gene overview page if not a protein
    print redirect(-url => "geneOverview.cgi?orgId=$orgId&gene=$locusId");
    exit(0);
}

# write the title
my $title = "Protein Info for $sys in $orginfo->{$orgId}{genome}";
my $start = Utils::start_page("$title");

my $tabs = Utils::tabsGene($dbh,$cgi,$orgId,$locusId,0,1,"protein");

# domains and families (PFam and TIGRFam)
my $cond = $dbh->selectall_arrayref(qq{SELECT * FROM GeneDomain WHERE orgId=? AND locusId=? ORDER BY begin;},
    { Slice => {} },
    $orgId, $locusId);

# protein features
my $features = $dbh->selectall_arrayref(qq{SELECT * FROM GeneFeature WHERE orgId=? AND locusId=? ORDER BY begin},
                                        { Slice => {} },
                                        $orgId, $locusId);
# only show transmambrane and signal peptide
my @featuresShow = grep $_->{featureType} eq "signal peptide" || $_->{featureType} eq "transmembrane", @$features;

# Fetch amino acid sequence
my $tmpDir = Utils::tmp_dir();
my $seqFile = "$tmpDir/$orgId+$locusId.fasta";
my $myDB = Utils::blast_db();
my $id = join(":",$orgId,$locusId);
my $fastacmd = '../bin/blast/fastacmd';
system($fastacmd,'-d',$myDB,'-s',$id,'-o',$seqFile)==0 || die "Error running $fastacmd -d $myDB -s $id -o $seqFile -- $!";
my $in = Bio::SeqIO->new(-file => $seqFile,-format => 'fasta');
my $seq = $in->next_seq()->seq;
my $seqLen = length($seq);
unlink($seqFile) || die "Error deleting $seqFile: $!";


# Reannotation information, if any
my $reanno = $dbh->selectrow_hashref("SELECT * from Reannotation WHERE orgId = ? AND locusId = ?",
                                     {}, $orgId, $locusId);

my @toplines = ();
push @toplines, "Name: $gene->{gene}" if $gene->{gene} ne "";
my $desc_title = "Annotation";
my $descUpdated = $gene->{desc} || "no description";
if ($reanno->{new_annotation}) {
  my $comment = $reanno->{comment};
  $comment =~ s/^# *//;
    push @toplines,
      "Updated annotation (from data): $reanno->{new_annotation}",
      "Rationale:  $comment";
    $desc_title = "Original annotation";
    $descUpdated = $reanno->{new_annotation};
}
push @toplines, "${desc_title}: " . ($gene->{desc} || "no description");

print
  header, $start, $tabs, '<div id="tabcontent">',
  h2("Protein Info for $sys in " . $cgi->a({href => "org.cgi?orgId=$orgId"},$orginfo->{$orgId}{genome})),
  p(join("<BR>", @toplines)),
  p(small(qq{These analyses and tools can help you predict a protein's function, but be skeptical.
             For enzymes, over 10% of annotations from KEGG or SEED are probably incorrect.
             For other types of proteins, the error rates may be much higher.
             MetaCyc and Swiss-Prot have low error rates,
             but the best hits in these databases are often quite distant,
             so this protein's function may not be the same.
             TIGRFam has low error rates.
             Finally, many experimentally-characterized proteins are not in any of these databases.
             To find relevant papers, use PaperBLAST.})),
  h3("Protein Families and Features");

my %ecall = (); # ec number => source => 1
# (source is one of TIGRFam, KEGG, SEED, reanno, MetaCyc)

if (@$cond == 0 && @featuresShow == 0) {
    print p("No protein families (PFam or TIGRFam), signal peptides, or transmembrane helices were found in this protein.");
} else {
  # Build a combined data structure, with one row per type of feature or domain
  # First signal peptide, then transmembrane, then domains sorted by 1st begin
  # Since signal peptide should be before transmembrane, just sort featuresShow that way too
  my %nameToRows = ();
  foreach my $row (@$cond) {
    push @{ $nameToRows{ $row->{domainId} } }, $row;
  }
  foreach my $row (@featuresShow) {
    $row->{domainDb} = 'phobius';
    push @{ $nameToRows{ $row->{featureType} } }, $row;
  }
  my @rows = sort { ($a->[0]{domainDb} ne 'phobius') <=> ($b->[0]{domainDb} ne 'phobius')
                      || $a->[0]{begin} <=> $b->[0]{begin} } (values %nameToRows);

  # Build an svg
  # Layout: 1 row per domain or type of feature
  # x axis: $padleft to $padLeft + seqWidth has the boxes
  #    there's also extra right space for labels
  # y axis: each row is of height $rowHeight (0 at top); and extra space at bottom for a labeled axis
  my $rowHeight = 45;
  my $totHeight = (@rows+0.65) * $rowHeight;
  my $seqWidthLow = 300;
  my $seqWidthHigh = 900;
  my $lenLow = 200;
  my $lenHigh = 800;
  my $seqWidth;
  if ($seqLen <= $lenLow) {
    $seqWidth = $seqWidthLow;
  } elsif ($seqLen > $lenHigh) {
    $seqWidth = $seqWidthHigh;
  } else {
    $seqWidth = $seqWidthLow + ($seqLen - $lenLow) * ($seqWidthHigh-$seqWidthLow) / ($lenHigh - $lenLow);
  }
  my $padLeft = 10;
  my $aaWidth = $seqWidth / $seqLen;
  my $totWidth = $padLeft + $seqWidthHigh + 100;
  my @svg = (); # the components within the svg

  # Add axis at bottom
  my $axisY = (scalar(@rows) + 0.15) * $rowHeight;
  my $tickBottom = (scalar(@rows) + 0.25) * $rowHeight;
  my $axisLabelTop = (scalar(@rows) + 0.45) * $rowHeight;
  my $left1 = $padLeft + $aaWidth;
  my $right = $padLeft + $seqWidth;
  push @svg, qq{<line x1="$left1" y1="$axisY" x2="$right" y2="$axisY" style="stroke:black; stroke-width:1"/>};
  my @tickAt = ();
  my $tickSpacing = 50;
  if ($seqLen < 50) {
    @tickAt = (1, $seqLen);
    $tickSpacing = $seqLen - 1;
  } else {
    $tickSpacing = 100 if $seqLen > 1000;
    $tickSpacing = 200 if $seqLen > 2000;
    $tickSpacing = 500 if $seqLen > 3000;
    $tickSpacing = 25 if $seqLen < 100;
    @tickAt = (1);
    for (my $i = $tickSpacing; $i < $seqLen; $i+= $tickSpacing) {
      push @tickAt, $i;
    }
    push @tickAt, $seqLen if $seqLen != $tickAt[-1];
  }
  foreach my $aa (@tickAt) {
    my $x = $padLeft + $aa * $aaWidth;
    push @svg, qq{<line x1="$x" y1="$axisY" x2="$x" y2="$tickBottom" style="stroke:black; stroke-width:1"/>};
    push @svg, qq{<text x="$x" y="$axisLabelTop" text-anchor="middle" alignment-baseline="central" style="font-size: 90%">$aa</text>}
      unless $aa > $seqLen - $tickSpacing/2 && $aa < $seqLen;
  }

  foreach my $i (0..(scalar(@rows)-1)) {
    my $list = $rows[$i];
    my $first = $list->[0];
    my $rowTop = $i * $rowHeight;
    #my $rowBottom = ($i+1) * $rowHeight;
    my $textBottom = ($i + 0.3) * $rowHeight;
    my $objTop = ($i + 0.4) * $rowHeight;
    my $objBottom = ($i + 0.75) * $rowHeight;
    my $objHeight = $objBottom - $objTop;
    my ($title, $titleURL, $color);

    if ($first->{domainDb} eq "phobius") {
      $color = $first->{featureType} eq "transmembrane" ? "red" : "magenta";
      $title = $first->{featureType};
    } elsif ($first->{domainDb} eq "TIGRFam") {
      $title = "$first->{domainId}: $first->{definition}";
      $titleURL = "https://www.ncbi.nlm.nih.gov/Structure/cdd/" . $first->{domainId};
      $color = "chocolate";
    } elsif ($first->{domainDb} eq "PFam") {
      $title = "$first->{domainId}: $first->{domainName}";
      $titleURL = "https://www.ebi.ac.uk/interpro/entry/pfam/$first->{domainId}/";
      $color = "lightblue";
    }
    if ($titleURL) {
      push @svg, qq{<a xlink:href="$titleURL">};
    }
    push @svg, qq{<text x="0.5" y="$textBottom">$title</text>"};
    if ($titleURL) {
      push @svg, qq{</a>};
    }

    # add rounded rectangles for each object and a small label with e-values below
    foreach my $obj (@$list) {
      my $beginScaled = $padLeft + $obj->{begin} * $aaWidth;
      my $widthScaled = ($obj->{end} - $obj->{begin}) * $aaWidth;
      my ($tooltip, $alignURL);
      my $nAA = $obj->{end} - $obj->{begin} + 1;
      if ($obj->{domainDb} eq "phobius") {
        $tooltip = "amino acids $obj->{begin} to $obj->{end} ($nAA residues), see Phobius details";
        $alignURL = "myPhobius.cgi?name=${sys}&seq=${seq}";
      } else {
        $tooltip = "amino acids $obj->{begin} to $obj->{end} ($nAA residues), $obj->{score} bits, see alignment";
        $alignURL = "showHmmAlign.cgi?orgId=$orgId&locusId=$locusId&domainDb=$obj->{domainDb}&domainName=$obj->{domainName}";
      }
      my $tooSmallForE = $widthScaled <= 40;
      $tooltip .= " (E = $obj->{evalue})" if exists $obj->{evalue} && $tooSmallForE;
      push @svg, qq{<a xlink:href="$alignURL">};
      push @svg, qq{<title>$tooltip</title>};
      push @svg, qq{<rect x="$beginScaled" y="$objTop" width="$widthScaled" height="$objHeight"
                      rx="5" ry="5"
                      style="fill:$color; stroke-width: 1; stroke: black;" /> };
      my $midScaled = $beginScaled + $widthScaled/2;
      my $objMiddle = ($objTop + $objBottom)/2;
      if (exists $obj->{evalue} && ! $tooSmallForE) {
        my $showE = $widthScaled > 70 ? "E=$obj->{evalue}" : $obj->{evalue};
        push @svg, qq{<text dominant-baseline="middle" text-anchor="middle"
                          x="$midScaled" y="$objMiddle"
                          style="font-size: 75%;">$showE</text>};
      }
      push @svg, qq{</a>};
    }
  }
  my $svg = join("\n",
                 qq{<svg width="${totWidth}" height="${totHeight}" style="position: relative; left: 1em;">},
                 qq{<g transform="scale(1)">},
                 @svg,
                 "</g>",
                 "</svg>");
  print "\n", div($svg);
}
autoflush STDOUT 1; # so preliminary results appear
print "\n";

print h3("Best Hits");

# SwissProt hit information, if any
my $bhSprot = $dbh->selectrow_hashref("SELECT * from BestHitSwissProt
                                         JOIN SwissProtDesc USING (sprotAccession,sprotId)
                                         WHERE orgId = ? AND locusId = ?",
                                      {}, $orgId, $locusId);
if (defined $bhSprot->{sprotAccession}) {
    my $acc = $bhSprot->{sprotAccession};
    print
        p("Swiss-Prot:",
            sprintf("%.0f%% identical to", $bhSprot->{identity}),
            a({ href => "http://www.uniprot.org/uniprot/$acc",
                title => $bhSprot->{sprotAccession} },
              $bhSprot->{sprotId} ) . ":",
            $bhSprot->{desc},
            $bhSprot->{geneName} ? "($bhSprot->{geneName})" : "",
            "from",
            $bhSprot->{organism} );
}
print "\n";

# KEGG information, if any
my $kegg = Utils::kegg_info($dbh, $orgId, $locusId);
if (defined $kegg) {
    my $ko = $kegg->{ko};
    my @kopieces = ();
    if (scalar(@{ $kegg->{ko} }) > 0) {
        foreach my $row (@$ko) {
            my $ecs = $row->{ec};
            my @ecdesc = ();
            foreach my $ec (@$ecs) {
                if ($ec =~ m/-/) {
                    push @ecdesc, $ec;
                } else {
                    push @ecdesc, a({-href => "http://www.kegg.jp/dbget-bin/www_bget?ec:$ec"}, $ec);
                }
                $ecall{$ec}{"KEGG"} = 1;
            }
            my $ecdesc = join(" ", @ecdesc);
            $ecdesc = " [EC: $ecdesc]" if $ecdesc ne "";

            push @kopieces,
              join("",
                   a({-href => "http://www.kegg.jp/dbget-bin/www_bget?ko:" . $row->{kgroup}},
                     $row->{kgroup}),
                   ", ",
                   $row->{desc} || "(no description)",
                   $ecdesc);
        }
    } else {
      push @kopieces, "None";
    }
    my $keggOrg = $kegg->{keggOrg};
    my $keggId = $kegg->{keggId};
    print join(" ", "KEGG orthology group:", @kopieces,
               sprintf("(inferred from %.0f%% identity to ", $kegg->{identity}),
               a({-href => "http://www.kegg.jp/dbget-bin/www_bget?$keggOrg:$keggId"},
                 "$keggOrg:$keggId").")");
}
print "\n";

# MetaCyc information, if any
my $bhMetacyc = $dbh->selectall_arrayref(qq{ SELECT * FROM BestHitMetacyc
                                             LEFT JOIN MetacycReaction USING (rxnId)
                                             WHERE orgId = ? AND locusId = ?
                                             ORDER BY rxnName DESC },
                                         { Slice => {} }, $orgId, $locusId);
my %metacycrxn = (); # which reactions are linked to
if (scalar(@$bhMetacyc) > 0) {
  my @showrxns = ();
  my @showecs = ();
  my %ecSeen = ();
  # Multiple hits if gene is linked to >1 reaction
  foreach my $row (@$bhMetacyc) {
    next if $row->{rxnId} eq "";
    $metacycrxn{ $row->{rxnId} } = 1;
    my $rxnName = $row->{rxnName};
    my $ecnums = $dbh->selectcol_arrayref("SELECT ecnum from MetacycReactionEC WHERE rxnId = ?",
                                          {}, $row->{rxnId});
    foreach my $ec (@$ecnums) {
      next if exists $ecSeen{$ec} || $ec eq "";
      $ecSeen{$ec} = 1;
      $ecall{$ec}{"MetaCyc"} = 1;
      push @showecs, $ec;
      # and maybe update the name
      # If no reaction name, or it is just some EC number(s) --
      # replace with KEGG description of first ec
      if (!defined $rxnName || $rxnName eq "" || $rxnName =~ m/^[0-9][.][0-9.,]+$/) {
        ($rxnName) = $dbh->selectrow_array("SELECT ecdesc FROM ECInfo WHERE ecnum = ?",
                                           {}, $ec);
      }
    }
    my $showrxn = a({ -href => "http://metacyc.org/META/NEW-IMAGE?type=REACTION&object=" . $row->{rxnId},
                      -title => "see $row->{rxnId} in MetaCyc" }, $rxnName || $row->{rxnId});
    if (@showecs > 0) {
      @showecs = map a({ -href => "https://metacyc.org/META/NEW-IMAGE?type=EC-NUMBER&object=EC-$_" }, $_), @showecs;
      $showrxn .= " [EC: " . join(", ", @showecs) . "]";
    }
    push @showrxns, $showrxn;
  }
  my $row = $bhMetacyc->[0];
  my $URL = "http://metacyc.org/gene?orgId=META&id=" . $row->{protId};
  my $desc = a({ href => $URL }, $row->{desc});
  $desc .= br() . small(join("; ", @showrxns)) if @showrxns > 0;
  print p( "MetaCyc:",
           sprintf("%.0f%% identical to", $row->{identity}),
           $desc);
}
print "\n";

# SEED information, if any
my @show_classes = ();
my ($seed_desc,$seed_classes) = Utils::seed_desc($dbh,$orgId,$locusId);
foreach my $row (@$seed_classes) {
    my ($type,$num) = @$row;
    my $url_pre = $type == 1 ? "http://www.kegg.jp/dbget-bin/www_bget?ec:"
        : "http://www.tcdb.org/search/result.php?tc=";
    my $text = ($type == 1 ? "EC " : "TC ").$num;
    push @show_classes, ($num =~ m/-/ ? $text
                         : a({-href => $url_pre . $num}, $text));
    $ecall{$num}{"SEED"} = 1 if $type == 1;
}
my $subsysShow = "";
if ($seed_desc) {
  my $subsysList = $dbh->selectcol_arrayref(qq{ SELECT DISTINCT subsystem FROM SeedAnnotationToRoles
                                                JOIN SEEDRoles USING (seedrole)
                                                WHERE seed_desc = ? },
                                            {}, $seed_desc);
  my @subsysShow = ();
  foreach my $subsys (@$subsysList) {
    my $nice = $subsys; $nice =~ s/_/ /g;
    push @subsysShow, a({ -href => "seedsubsystem.cgi?orgId=$orgId&subsystem=$subsys" },
                        $nice);
  }
  $subsysShow = " in subsystem " . join(" or ", @subsysShow) if @subsysShow > 0;
}

print
    h3("Predicted SEED Role"),
    p(defined $seed_desc ?  '"' . $seed_desc . '"' . $subsysShow : "No annotation",
      @show_classes > 0 ? "(" . join(", ", @show_classes) . ")" : "");
print "\n";

# And add reannotated EC numbers to %ecall
# (Above we already did TIGRFam, KEGG, SEED)
my $reannoEc = $dbh->selectcol_arrayref(qq{ SELECT ecnum FROM ReannotationEC
                                            WHERE orgId = ? AND locusId = ? },
                                        {}, $orgId, $locusId);
foreach my $ec (@$reannoEc) {
  $ecall{$ec}{"reanno"} = 1;
}

# Also link to metacyc reactions via SEED roles
my $rxnsFromSEED = $dbh->selectcol_arrayref(qq{ SELECT rxnId FROM SEEDAnnotation
                                           JOIN SEEDAnnotationToRoles USING (seed_desc)
                                           JOIN SeedRoleReaction USING (seedrole)
                                           JOIN SEEDReaction USING (seedrxnId)
                                           JOIN MetacycReaction USING (keggrxnId)
                                           WHERE orgId = ? AND locusId = ? AND keggrxnId <> "" },
                                            {}, $orgId, $locusId);
foreach my $rxnId (@$rxnsFromSEED) {
  $metacycrxn{$rxnId} = 1;
}

# And expand the list of potential MetaCyc reactions by using EC numbers
foreach my $ec (keys %ecall) {
  my $ecrxns = $dbh->selectcol_arrayref("SELECT rxnId FROM MetacycReactionEC WHERE ecnum = ?",
                                        {}, $ec);
  foreach my $rxnId (@$ecrxns) {
    $metacycrxn{$rxnId} = 1;
  }
}

# Now that we have all the EC numbers and the potential MetaCyc reactions --
# link to MetaCyc pathways and to KEGG maps

if (keys(%metacycrxn) > 0) {
  # Find all relevant MetaCyc pathways, using MetacycPathwayReaction
  my %pathways = (); # pathwayId => pathway object
  foreach my $rxnId (keys %metacycrxn) {
    my $paths = $dbh->selectall_arrayref(qq{ SELECT pathwayId, pathwayName, nFound, nSteps
                                             FROM MetacycPathwayReaction
                                             JOIN MetacycPathway USING (pathwayId)
                                             JOIN MetacycPathwayCoverage USING (pathwayId)
                                             WHERE rxnId = ? AND orgId = ? },
                                         { Slice => {} }, $rxnId, $orgId);
    foreach my $path (@$paths) {
      $pathways{ $path->{pathwayId} } = $path;
    }
  }
  foreach my $path (values %pathways) {
    $path->{score} = Utils::MetacycPathwayScore($path->{nSteps}, $path->{nFound});
  }
  my @path = sort { $b->{score} <=> $a->{score}
                      || $a->{nSteps} <=> $b->{nSteps}
                        || $a->{pathwayName} cmp $b->{pathwayName} } values %pathways;

  my @pathList = ();
  foreach my $path (@path) {
    my $pathId = $path->{pathwayId};
    my $nFound = $path->{nFound};
    my $nSteps = $path->{nSteps};
    push @pathList, li(a( {-href => "pathway.cgi?orgId=$orgId&pathwayId=$pathId"},
                          $path->{pathwayName} ),
                       "($nFound/$nSteps steps found)");
  }
  print h3("MetaCyc Pathways"), start_ul(), join("\n", @pathList), end_ul(), "\n"
    if @pathList > 0;
}
print "\n";

# KEGG maps and isozymes
if (keys(%ecall) > 0) {
    my @ec = sort keys %ecall;

    my @ecspec = map { "'" . $_ . "'" } @ec;
    my $ecspec = join(",", @ecspec);
    my $maps = $dbh->selectall_arrayref(qq{SELECT DISTINCT mapId,title
                                           FROM KEGGConf JOIN KEGGMap USING (mapId)
                                           WHERE type=1 AND objectId IN ( $ecspec )
                                           ORDER BY title});
    if (scalar(@$maps) > 0) {
        my @rendered = ();
        foreach my $row (@$maps) {
            my ($mapId,$mapdesc) = @$row;
            push @rendered, li(a({href => "keggmap.cgi?orgId=$orgId&mapId=$mapId&ec="
                                      . join(",", @ec)},
                                 $mapdesc));
        }
        print
            h3("KEGG Metabolic Maps"),
            "<ul>",
            join("", @rendered),
            "</ul>";
    }

    my @links = ();
    my $ecGenes = Utils::EcToGenes($dbh, $orgId, \@ec);

    foreach my $ec (@ec) {
        my @iso = keys %{ $ecGenes->{$ec} };
        my @iso2 = grep { $_ ne $locusId } @iso;
        if (@iso2 > 0) {
            my @locusSpecs = map "locusId=$_", @iso;
            push @links, a( { -href => "genesFit.cgi?orgId=$orgId&" . join("&", @locusSpecs) },
                            $ec );
        }
    }
    my @curatedlinks = map a({-href => "http://papers.genomics.lbl.gov/cgi-bin/genomeSearch.cgi?gdb=FitnessBrowser&gid=$orgId&orgId=$orgId&query=$_&word=1" },
                             $_), @ec;
    print
        h3("Isozymes"),
        scalar(@links) > 0 ? p("Compare fitness of predicted isozymes for:", join(", ", @links))
          : "No predicted isozymes",
        p("Use Curated BLAST to search for", join(" or ", @curatedlinks));
}
print "\n";

# print sequence
$seq =~ s/(.{60})/$1\n/gs;
my %positiveDivisions = map { $_ => 1 } qw{Actinobacteria Firmicutes};
my %archaeaDivisions = map { $_ => 1 } qw{Euryarchaeota Crenarchaeota};
my ($psortType, $psortShow) = ("negative", "Gram-negative bacteria");
if (exists $positiveDivisions{ $orginfo->{$orgId}{division} }) {
  ($psortType, $psortShow) = ("positive", "Gram-positive bacteria");
} elsif (exists $archaeaDivisions{ $orginfo->{$orgId}{division} }) {
  ($psortType, $psortShow) = ("archaea", "archaea");
}

my $newline = "%0A";
print
  h3("Sequence Analysis Tools"),
  p(a({-href => "http://papers.genomics.lbl.gov/cgi-bin/litSearch.cgi?query=>${sys}$newline$seq"},
      "PaperBLAST"),
    "(search for papers about homologs of this protein)"),
  p(a({-href => "http://www.ncbi.nlm.nih.gov/Structure/cdd/wrpsb.cgi?seqinput=>${sys}$newline$seq"},
      "Search CDD"),
    "(the Conserved Domains Database, which includes COG and superfam)"),

  # See documentation of HMMer web server (version 1.0) at
  # https://hmmer-web-docs.readthedocs.io/en/latest/searches.html
  # Set E= and domE= to cause weak hits to appear (they are still labeled insignificant)
  # As of October 2022, it does not seem to be functional, so removed
#  p(start_form(-name => "PfamForm", -id => "PfamForm",
#               -method => "POST", -action => "https://www.ebi.ac.uk/Tools/hmmer/search/hmmscan"),
#    hidden('seq', ">$sys\n$seq"),
#    hidden('hmmdb', 'pfam'),
#    hidden('E', '1'),
#    hidden('domE', '1'),
#    submit(-style => "display: none;", -id => 'PfamButton'),
#    end_form,
#    a({-href => "javascript:document.getElementById('PfamForm').submit()"},
#       "Search PFam"),
#    "(including for weak hits, up to E = 1)"),

  p("Compare to",
    a({ -href => "http://papers.genomics.lbl.gov/cgi-bin/pdbBlast.cgi?query=>${sys}%0A$seq",
        -title => "Find similar proteins with known structures" },
      "protein structures")),
  p("Predict protein localization: ",
    a({-href => "https://papers.genomics.lbl.gov/cgi-bin/psortb.cgi?name=${sys}&type=${psortType}&seq=${seq}",
       -title => "PSORTb v3.0 for $psortShow"},
      "PSORTb"),
    "($psortShow)"),
  p("Predict transmembrane helices and signal peptides:",
    a({-href => "myPhobius.cgi?name=${sys}&seq=${seq}"},
      "Phobius")),
  p("Check the current SEED with",
    # %0A encodes "\n" so that it looks like fasta input.
    a({-href => "http://pubseed.theseed.org/FIG/seedviewer.cgi?page=FigFamViewer&fasta_seq=>${sys}%0A$seq"},
      "FIGfam search")),
  p("Find homologs in",
    a({-href => "https://fast.genomics.lbl.gov/cgi/findHomologs.cgi?seqDesc=$sys&seq=$seq",
       -title => "Quickly find homologs in representative genomes"},
      i("fast.genomics")),
    "or the",
    a({-href => "https://iseq.lbl.gov/genomes/seqsearch?sequence=>${sys}%0A$seq"},
      "ENIGMA genome browser")),
  "\n";
  # Formerly used
  # <SCRIPT>document.getElementById("InterPro").innerHTML = "$html"</SCRIPT>
  # to set this with a delay
  #qq{<DIV id="InterPro"><p>&nbsp;</p></DIV>},

# Show link to InterPro, if possible
my $xrefU = $dbh->selectrow_hashref("SELECT * FROM LocusXref WHERE xrefDb = 'uniprot' AND orgId = ? AND locusId = ?",
                                    {}, $orgId, $locusId);
if (defined $xrefU && exists $xrefU->{xrefId}) {
  my $uniprotId = $xrefU->{xrefId};
  # (Earlier version of the tables had the long id; db_setup.pl no longer does this)
  # my $url = "https://rest.uniprot.org/uniprotkb/$uniprotId.fasta";
  # my $fasta = get($url); # a bit slow
  # $shortId = $1 if $fasta && $fasta =~ m/^>[a-z]+[|]([A-Z][A-Z0-9]+)[|]/;
  print p("See $uniprotId at",
          a({-href => "https://www.uniprot.org/uniprotkb/$uniprotId/entry"}, "UniProt"),
          "or",
          a({-href => "https://www.ebi.ac.uk/interpro/protein/$uniprotId"}, "InterPro"));
} else {
  print p("Find the",
          a({-href => "https://fast.genomics.lbl.gov/cgi/bestHitUniprot.cgi?query=>${sys}%0A$seq",
             -title => "See UniProt's annotation, the predicted structure, and protein families from InterPro"},
            "best match"),
          "in UniProt");
}

print
  h3("Protein Sequence ($seqLen amino acids)"),
  pre(">$sys $descUpdated ($orginfo->{$orgId}{genome})\n$seq"),
  "\n";

print '</div>';
$dbh->disconnect();
Utils::endHtml($cgi);
