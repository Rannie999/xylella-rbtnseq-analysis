#!/usr/bin/perl -w
#######################################################
## geneOverview.cgi
##
## Copyright (c) 2015 University of California
##
## Authors:
## Victoria Lo, Wenjun Shao (wjshao@berkeley.edu) and 
## Morgan Price
#######################################################
#
# Required CGI parameters:
# gene -- a locusId, sysName, or gene name to match on
#	(may show multiple hits)
# Optional CGI parameters:
# orgId -- which organism to search in
# showAll -- 1 if showing all fitness values instead of just the most extreme ones

use strict;

use CGI qw(:standard Vars -nosticky);
use CGI::Carp qw(warningsToBrowser fatalsToBrowser);
use URI::Escape qw{uri_escape};
use DBI;
use List::Util qw{min max};

use lib "../lib";
use Utils;

my $cgi=CGI->new;
print $cgi->header;
my $style = Utils::get_style();

my $orgSpec = $cgi->param('orgId') || "";
my $geneSpec = $cgi->param('gene') || "";
my $showAll = $cgi->param('showAll') ? 1 : 0;
my $start = Utils::start_page(uri_escape("Overview for $geneSpec ($orgSpec)"));

$geneSpec =~ s/ *$//;
$geneSpec =~ s/^ *//;

if (!defined $geneSpec || $geneSpec eq "") {
    Utils::fail($cgi, "you must enter the gene name or locus tag");
}

# check user input
Utils::fail($cgi, "$geneSpec is invalid. Please enter correct gene name!") unless ($geneSpec =~ m/^[A-Za-z0-9_.-]*$/);

# connect to database
my $dbh = Utils::get_dbh();
my $orginfo = Utils::orginfo($dbh);
die "Invalid orgSpec $orgSpec" unless $orgSpec eq "" || exists $orginfo->{$orgSpec};

my $query = qq{SELECT * FROM Gene
		WHERE ( locusId = ? OR sysName = ? OR upper(gene) = upper(?) )};
my $hits;
if ($orgSpec) {
    $query .= " AND orgId = ?";
    $hits = $dbh->selectall_arrayref($query, { Slice => {} }, $geneSpec, $geneSpec, $geneSpec, $orgSpec);
} else {
    $hits = $dbh->selectall_arrayref($query, { Slice => {} }, $geneSpec, $geneSpec, $geneSpec);
}

# and add another column for whether there is fitness data
foreach my $gene (@$hits) {
    $gene->{has_fitness} = Utils::gene_has_fitness($dbh,$gene->{orgId},$gene->{locusId});
}

if (@$hits == 0) {
	print $start,
    print $cgi->h3("No gene found for $geneSpec",
		   (exists $orginfo->{$orgSpec}{genome} ? " in " . $orginfo->{$orgSpec}{genome} : ""));
} elsif (@$hits > 1) {
	print $start,

	h3("Genes found for $geneSpec:");
    my @trows = ();
    push @trows, $cgi->Tr({-align=>'CENTER',-valign=>'TOP'},
			  $cgi->th( [ 'geneId','sysName','geneName','description','genome','fitness' ] ) );
    foreach my $gene (@$hits) {
	my ($fitstring, $fittitle) = Utils::gene_fit_string($dbh, $gene->{orgId}, $gene->{locusId});
	my @trow = map $cgi->td($_), ($gene->{locusId}, $gene->{sysName}, $gene->{gene}, $gene->{desc},
				      $cgi->a({href => "org.cgi?orgId=". $orginfo->{$gene->{orgId}}->{orgId}}, "$orginfo->{$gene->{orgId}}->{genome}"),
				      a( {href => "myFitShow.cgi?orgId=$gene->{orgId}&gene=$gene->{locusId}", title => $fittitle, },
					 $fitstring));
	push @trows, $cgi->Tr(@trow);
    }
    
    print $cgi->table( { cellspacing=>0, cellpadding=>3 }, @trows);
    # PRINT TABLE

} else {
    # just 1 hit
    my $gene = $hits->[0];
    my $orgId = $gene->{orgId};
    my $locusId = $gene->{locusId};
    my $scaffold = $gene->{scaffoldId};
    my $begin = $gene->{begin};
    my $end = $gene->{end};
    my $strand = $gene->{strand};
    my $type = $gene->{type};
    my $typeName = "";
    $typeName = "Protein-coding gene" if $type == 1;
    $typeName = "23S (large subunit) ribosomal RNA" if $type == 2;
    $typeName = "16S (small subunit) ribosomal RNA" if $type == 3;
    $typeName = "5S ribosomal RNA" if $type == 4;
    $typeName = "Transfer RNA" if $type == 5;
    $typeName = "Other non-coding RNA" if $type == 6;
    $typeName = "Pseudogene derived from a protein-coding gene" if $type == 7;
    $typeName = "Pseudogene derived from an RNA gene" if $type == 8;
    $typeName = "CRISPR" if $type == 9;
    $typeName = "CRISPR spacer" if $type == 10;
    $typeName = "Antisense RNA" if $type == 11;
    $typeName = "Unclassified feature (possibly a pseudogene)" if $type == 99;

    #nearby  genes, 5 on each side
    my $idShow = $gene->{sysName} || $gene->{locusId};
    my %spacingDesc = (); # locusId => spacing description
    my $centralId = $locusId;
    my $scgenes = $dbh->selectall_arrayref("SELECT * from Gene where orgId = ? AND scaffoldId = ? ORDER BY begin",
					   { Slice => {} }, $orgId, $gene->{scaffoldId});
    die "Cannot find genes for $gene->{scaffoldId}" unless scalar(@$scgenes) > 0;
    my ($iCentral) = grep { $scgenes->[$_]{locusId} eq $centralId } (0..(scalar(@$scgenes)-1));
    die if !defined $iCentral;
    my $i1 = max(0, $iCentral - 5);
    my $i2 = min(scalar(@$scgenes)-1, $iCentral + 5);
    my @locusIds = map { $scgenes->[$_]{locusId} } ($i1..$i2);

    my @genes = ();
    foreach my $locusId2 (@locusIds) {
      my $gene = $dbh->selectrow_hashref("SELECT * FROM Gene WHERE orgId = ? AND locusId = ?",
                                         {}, $orgSpec, $locusId2);
      die "No such locus $locusId2 in org $orgId" if !defined $gene->{locusId};
      # expName => "fit" => fitness value
      $gene->{fit} = $dbh->selectall_hashref(qq{SELECT expName,fit,t FROM GeneFitness
	                                               WHERE orgId = ? AND locusId = ?},
                                             "expName", {}, $orgId, $locusId2);
      $gene->{nExps} = scalar(keys %{ $gene->{fit} });
      push @genes, $gene;
    }

    my $start = Utils::start_page("Overview for $idShow ($orginfo->{$orgId}{genome})");
    my $title = "Gene Info for $idShow in $orginfo->{$orgId}{genome}";
    my $tabs = Utils::tabsGene($dbh,$cgi,$orgSpec,$geneSpec,$showAll,$gene->{type},"gene");
    my $beginC = Utils::commify($begin);
    my $endC = Utils::commify($end);

    # Reannotation information, if any
    my $reanno = $dbh->selectrow_hashref("SELECT * from Reannotation WHERE orgId = ? AND locusId = ?",
                                         {}, $orgId, $locusId);

    my @lines = ();
    push @lines, "Name: $gene->{gene}" if $gene->{gene} ne "";
    my $desc_title = "Annotation";
    if ($reanno->{new_annotation}) {
      my $comment = $reanno->{comment};
      $comment =~ s/^# *//;
      push @lines,
        "Updated annotation (from data): <B>$reanno->{new_annotation}</B>",
          "Rationale:  $comment";
      $desc_title = "Original annotation";
    }
    my $desc = $gene->{desc};
    $desc = b($desc) if $desc ne "" && ! $reanno->{new_annotation};
    $desc = "no description" if $desc eq "";
    push @lines, "${desc_title}: $desc";

    push @lines, "Type $type: $typeName";
    push @lines, "Located on scaffold $scaffold, $strand strand, nucleotides $beginC to $endC";
    push @lines, join(" ",
                      "See nucleotide sequence around",
                      a({-href => "getNtSeq.cgi?orgId=$orgId&locusId=$locusId" }, "$idShow"),
                      "or for all of",
                      a({-href => "getNtSeq.cgi?orgId=$orgId&scaffoldId=$scaffold"}, "scaffold $scaffold"));

    print
      $start, $tabs,
      h2("Gene info for $idShow in " . $cgi->a({href => "org.cgi?orgId=$orgId"}, "$orginfo->{$orgId}{genome}")),
      p(join(br(), @lines));

    my @links = ();
    if ($gene->{locusId} =~ m/^\d+$/) {
	push @links, $cgi->a({href => "http://www.microbesonline.org/cgi-bin/fetchLocus.cgi?locus=$gene->{locusId}"},
			     "MicrobesOnline");
    }
    if ($orgId eq "Keio" && $gene->{sysName} =~ m/^b\d+$/) {
	push @links, $cgi->a({href => "http://ecocyc.org/ECOLI/search-query?type=GENE&gname=$gene->{sysName}"}, "EcoCyc");
    }

    # Linking to NCBI is tricky -- for example, Cupriavidus basilensis 4G11 is in NCBI but the locus tags are not indexed yet.
    # Similarly PS417_05815 from P. fluorescens WCS417 -- hits P. simaie recA instead if search Protein, wierd.
    # So assume that it will work if the genome has an NCBI taxonomyId. This is a hack.

    push @links, $cgi->a({href => "http://www.ncbi.nlm.nih.gov/protein/?term=$gene->{sysName}#see-all"}, "NCBI")
        if $gene->{type} == 1 && $gene->{sysName} ne "" && $orginfo->{$orgId}{taxonomyId} ne "";

    # Add external links
    my $xrefs = $dbh->selectall_arrayref(qq{SELECT * from LocusXref WHERE orgId = ? AND locusId = ? ORDER by xrefDb,xrefId},
                                         { Slice => {} }, $orgId, $locusId);
    foreach my $xref (@$xrefs) {
      my $xrefDb = $xref->{xrefDb};
      my $xrefId = $xref->{xrefId};
      my $URL = undef;
      if ($xrefDb eq "uniprot") {
        $URL = "https://www.uniprot.org/uniprotkb/$xrefId/entry";
      } elsif ($xrefDb eq "IMG") {
        $URL = "https://img.jgi.doe.gov/cgi-bin/m/main.cgi?section=GeneDetail&page=geneDetail&gene_oid="
          . $xrefId;
      }
      my $show = "$xrefDb: $xrefId";
      $show = a({-href => $URL}, $show) if defined $URL;
      push @links, $show;
    }

    print $cgi->p("Other databases: " . join(", ", @links)) if (@links > 0);

    print $cgi->h3({style=>'text-align:center'},"Nearby Genes");
    print Utils::geneArrows($dbh, \@genes, $locusId, undef, undef);

    my @trows = ( Tr({-valign => "top", -align => 'center'},
                     th([ 'Locus', 'Name', 'Description', 'Strand',
                          'Phenotype', ])));

    foreach my $row (@genes) {
      my $label = $row->{gene} || $row->{sysName} || ""; #|| $row->{locusId};
      my $label2 =  $row->{gene} || $row->{sysName} || $row->{locusId};
      my $label3 = $label2;
      $label3 =~ s/^.*_/_/ if $row->{sysName} || $row->{locusId};
      $label2 = $row->{sysName}. ": " . $label2 if $row->{sysName};
      my $bgcolor = undef;
      $bgcolor = "#FFFFFF" if $row->{locusId} eq $locusId;
      my ($phen, $tip) = Utils::gene_fit_string($dbh,$orgSpec,$row->{locusId});
      push @trows, Tr({ -valign => 'top', -align => 'left', -bgcolor => $bgcolor},
                      td([ Utils::gene_link($dbh, $row, "name", "geneOverview.cgi"),
                           a({href => "geneOverview.cgi?orgId=$orgId&gene=$row->{locusId}"},$row->{gene} || $row->{sysName}), 
                           Utils::gene_link($dbh, $row, "desc", "domains.cgi"),
                           $row->{strand},
                           a({href => "myFitShow.cgi?orgId=$orgId&gene=$row->{locusId}", title=>$tip},$phen),
                         ]));
    }
    print table({cellspacing => 0, cellpadding => 3}, @trows);

    my $browsebegin = $genes[0]{begin};
    my $browseend = $genes[-1]{end};
    print p("Or",
            a({-href => "genomeBrowse.cgi?orgId=$orgId&scaffoldId=$scaffold&begin=$browsebegin&end=$browseend&zoom=out"},
              "browse"),
            "nearby genes");
    print p("Or see the",
            a({-href => "cmpbrowser.cgi?anchorOrg=$orgId&anchorLoci=$locusId" },
              "comparative fitness browser"))
      if $gene->{type} == 1;
}

$dbh->disconnect();
Utils::endHtml($cgi);
