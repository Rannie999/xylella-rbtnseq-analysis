#!/usr/bin/perl -w
#######################################################
## blastOrg.cgi
#
# Required CGI parameters: orgId
# Optional CGI paramaters:
#  query (protein fasta or sequence)
#  sixframe (1 if search the 6-frame translation of the genome)

use strict;

use CGI qw(:standard Vars);
use CGI::Carp qw(warningsToBrowser fatalsToBrowser);
use Time::HiRes qw(gettimeofday);
use DBI;
use IO::Handle; # for autoflush
use HTML::Entities; # for encode_entities()
use URI::Escape; # for uri_escape()
use HTML::Entities; # for encode_entities()
use List::Util qw{min max};

use lib "../lib";
use FEBA_Utils qw{ReadFastaEntry};
use Utils;

my $cgi=CGI->new;
my $dbh = Utils::get_dbh();

my $fastacmd = '../bin/blast/fastacmd';
my $formatdb = '../bin/blast/formatdb';
my $blastall = '../bin/blast/blastall';
my $usearch = '../bin/usearch';
foreach my $x ($fastacmd, $formatdb, $blastall, $usearch) {
  die "No such executable: $x\n" unless -x $x;
}

my $orgId = $cgi->param('orgId') || die "No orgId parameter";
my $orginfo = Utils::orginfo($dbh);
Utils::fail($cgi, "Unknown organism: $orgId") unless exists $orginfo->{$orgId};

my $query = $cgi->param('query') || "";

if ($query =~ m/[A-Za-z]/) {
  # parse the input sequence
  my $seq = "";
  my $def = "";
  my @lines = split /[\r\n]+/, $query;
  $def = shift @lines if @lines > 0 && $lines[0] =~ m/^>/;
  $def =~ s/^>//;
  foreach (@lines) {
    s/[ \t]//g;
    s/^[0-9]+//; # leading digit/whitespace occurs in UniProt format
    next if $_ eq "//";
    Utils::fail($cgi,"Error: more than one sequence was entered.") if m/^>/;
    Utils::fail($cgi,"Unrecognized characters in $_") unless m/^[a-zA-Z*]*$/;
    s/[*]/X/g;
    $seq .= uc($_);
  }
  $def = substr($seq,0,20) . "..." if $def eq "";
  my $defEncoded = encode_entities($def); # the defline may have wierd characters
  my $defURI = uri_escape($def);

  my $procId = $$;
  my $timestamp = int (gettimeofday * 1000);
  my $tmpDir = Utils::tmp_dir();
  my $tmpPre = "$tmpDir/blastOrg.$procId.$timestamp";

  my $sixframe = param('sixframe') || 0;
  my $text = $sixframe ? "BLAST against the six-frame translation of" : "BLAST hits in";
  print
     header,
     Utils::start_page($text . " " . $orginfo->{$orgId}{genome}),
     '<div id="ntcontent">',
     h3($text, a({-href => "org.cgi?orgId=$orgId"}, $orginfo->{$orgId}{genome})),
     "\n";

  my $pGenes = $dbh->selectall_arrayref("SELECT * from Gene WHERE orgId = ? AND type = 1 ORDER BY scaffoldId, begin",
                                        { Slice => {} }, $orgId);
  my %subjectLen;
  my %subjectToGene;
  my %scGenes; # protein-coding genes by scaffold
  my %scSeq;
  if ($sixframe) {
    my $scSeqs = $dbh->selectall_arrayref("SELECT * FROM ScaffoldSeq WHERE orgId = ?",
                                          { Slice => {} }, $orgId);
    %scSeq = map { $_->{scaffoldId} => $_->{sequence } } @$scSeqs;
    die "No scaffolds for $orgId" unless @$scSeqs > 0;
    open(my $fhFna, ">", "$tmpPre.fna") || die "Cannot write to $tmpPre.fna";
    foreach my $row (@$scSeqs) {
      my $sc = $row->{scaffoldId};
      my $seq = $row->{sequence};
      print $fhFna Utils::formatFASTA($sc,$seq)."\n";
    }
    my $usearchCmd = "$usearch -fastx_findorfs $tmpPre.fna -mincodons 25 -orfstyle 7 -aaout $tmpPre.subject -quiet";
    system($usearchCmd) == 0 || die "Error running $usearchCmd\n-- $!";
    foreach my $gene (@$pGenes) {
      push @{ $scGenes{ $gene->{scaffoldId} } }, $gene;
    }
  } else {
    %subjectToGene = map { $_->{locusId} => $_ } @$pGenes;
    # AASeq is not in the sqlite3 database, so need to use fastacmd to get them all
    die "No proteins in $orgId" unless @$pGenes > 0;
    open(my $fhList, ">", "$tmpPre.list") || die "Error writing to $tmpPre.list";
    foreach my $gene (@$pGenes) {
      print $fhList "$orgId:" . $gene->{locusId}. "\n";
    }
    close($fhList) || die "Error writing to $tmpPre.list";
    my @fastaCmd = ($fastacmd, '-d', Utils::blast_db(), '-i', "$tmpPre.list", '-o', "$tmpPre.subject");
    system(@fastaCmd) == 0
      || die "Error running " . join(" ", @fastaCmd) . "\n -- $!";
  }

  # Read lengths of subjects
  open (my $fhS, "<", "$tmpPre.subject") || die "Cannot read $tmpPre.subject";
  my $state = {};
  while (my ($id, $seq) = ReadFastaEntry($fhS, $state)) {
    $id =~ s/ .*//;
    if ($sixframe) {
      $subjectLen{$id} = length($seq);
    } else {
      my (undef, $locusId) = split /:/, $id;
      die unless defined $locusId && exists $subjectToGene{$locusId};
      $subjectLen{$locusId} = length($seq);
    }
  }
  close($fhS) || die "Error reading $tmpPre.subject";

  if ($sixframe) {
    print p("Searching", scalar(keys %subjectLen), "potential reading frames for homologs of", $defEncoded);
  } else {
    print p("Searching the", scalar(keys %subjectLen), "annotated proteins for homologs of", $defEncoded);
  }
  print "\n";

  open(my $fhQ, ">", "$tmpPre.query") || die "Cannot write to $tmpPre.query";
  print $fhQ ">query\n$seq\n";
  close($fhQ) || die "Error writing to $tmpPre.query";

  my @formatCmd = ($formatdb, '-p', 'T', '-i', "$tmpPre.subject");
  system(@formatCmd); # ignore spurious errors

  my @blastCmd = ($blastall, '-p', 'blastp',
                  '-i', "$tmpPre.query",
                  '-d', "$tmpPre.subject",
                  '-o', "$tmpPre.hits",
                  '-E', 1,
                  '-m', 8,
                  '-F', 'm S');
  system(@blastCmd) == 0 || die "Error running " . join(" ", @blastCmd) . "\n -- $!";

  my @hits;
  open (my $fhHits, "<", "$tmpPre.hits") || die "Cannot read $tmpPre.hits";
  while(my $line = <$fhHits>) {
    chomp $line;
    my @F = split /\t/, $line;
    push @hits, \@F;
  }
  foreach my $suffix (qw{fna list subject subject.pin subject.psi subject.phr subject.psd subject.psq query hits}) {
    unlink("$tmpPre.$suffix");
  }
  if (@hits == 0) {
    print p("No hits found with E < 1"), "\n";
  } else {
    print p("Found", scalar(@hits), "hits with E < 1"), "\n";
    my $seqlen = length($seq);
    my @header_IdCovE = ( a({-title => "Percent identity", -style=>"color: black;"}, '%Id'),
                          a({-title => "Percent coverage of query", -style=>"color: black;"}, 'Cov'),
                          a({-title => "expected number of hits this strong", -style => "color: black;"}, 'E') );
    if ($sixframe) {
      my @header = ('Scaffold', 'Frame', 'Range', 'Gene', @header_IdCovE);
      print qq{<TABLE style="border-bottom: none;" cellspacing=0 cellpadding=3 >},
        Tr(map th($_), @header),
        "\n";
      foreach my $hit (@hits) {
        my ($queryId,$subjectId,$percIdentity,$alnLength,$mmCnt,$gapCnt,
            $queryStart,$queryEnd,$subjectStart,$subjectEnd,$eVal,$bitScore) = @$hit;
        my ($scaffoldId, $code) = split /[|]/, $subjectId;
        die "Unknown scaffold $scaffoldId in $subjectId" unless exists $scSeq{$scaffoldId};
        $code =~ m/^([+-]\d):(\d+)-(\d+)[(]/ || die "Cannot parse code $code from orf id $subjectId";
        my ($frame, $orfBegin, $orfEnd) = ($1,$2,$3);
        my $subjectLen = $subjectLen{$subjectId};
        my ($s, $hitBegin, $hitEnd);
        if ($frame > 0) {
          $s = 1;
          # If the hit is 1 a.a., then subjectStart = subjectEnd and end = begin + 2
          $hitBegin = $orfBegin + ($subjectStart-1) * 3;
          $hitEnd = $orfBegin + ($subjectEnd-1) * 3 + 2;
        } else {
          $s = -1;
          $hitBegin = $orfEnd - ($subjectEnd-1) * 3 - 2;
          $hitEnd = $orfEnd - ($subjectStart-1) * 3;
        }
        
        my $orfObject = join(":", "b", $orfBegin, "e", $orfEnd, "s", $s, "n", "ORF");
        my $hitObject = join(":", "b", $hitBegin, "e", $hitEnd, "s", $s, "n", "hit");
        # Find the annotated protein, if any, that likely corresponds to this. Unfortunately findx_orfs has wierd coordinates so do
        # not try to find an exact match; instead look for ORFs in the same frame that contain most of the alignment
        my $overlapGene = undef;
        foreach my $gene (@{ $scGenes{$scaffoldId} }) {
          next unless $hitBegin % 3 == $gene->{begin} % 3; # same frame
          if ($gene->{strand} eq "+") {
            next unless $frame > 0;
          } else {
            next unless $frame < 0;
          }
          next unless $gene->{begin} <= $hitEnd && $gene->{end} >= $hitBegin;
          my $overlapBeg = max($hitBegin, $gene->{begin});
          my $overlapEnd = min($hitEnd, $gene->{end});
          if ($overlapEnd - $overlapBeg + 1 >= 0.5 * ($hitEnd-$hitBegin+1)) {
            $overlapGene = $gene;
            last;
          }
        }
        my @row = ($scaffoldId, $frame,
                   $cgi->a({ href => "genomeBrowse.cgi?orgId=$orgId&scaffoldId=$scaffoldId&object=$orfObject&object=$hitObject" },
                           $frame > 0 ? "${orfBegin}:${orfEnd}" : "${orfEnd}:${orfBegin}"),
                   $overlapGene ? Utils::gene_link($dbh, $overlapGene, "name", "geneOverview.cgi") : "",
                   $cgi->a({title=>"$bitScore bits", -style=>"color:black;"}, int(0.5 + $percIdentity)),
                   $cgi->a({title => "$queryStart:$queryEnd / $seqlen of query and $subjectStart:$subjectEnd / $subjectLen of subject",
                            style => "color:black;"},
                           int(0.5 + 100*abs($queryEnd - $queryStart + 1)/length($seq))),
                   $cgi->a({title=>"bits: $bitScore", style=>"color:black;"}, $eVal));
        my @td = map td($_), @row;
        print $cgi->Tr({ -align => 'left', '-valign' => 'top', bgcolor => 'white'}, @td);
        print "\n";
      }
      print "</TABLE>\n";
    } else { # not 6-frame
      my @header = ('Gene', 'Name', 'Description', 'Fitness', @header_IdCovE);
      my @widths = ("14%", "5%", "37%", "6%", "5%", "5%", "8%");
      my @th = map th({width => $widths[$_]}, small($header[$_])), (0..$#header);
      print qq{<TABLE style="border-bottom: none;" cellspacing=0 cellpadding=3 width=99% >},
        Tr(@th),
        "\n";
      foreach my $hit (@hits) {
        my ($queryId,$subjectId,$percIdentity,$alnLength,$mmCnt,$gapCnt,
            $queryStart,$queryEnd,$subjectStart,$subjectEnd,$eVal,$bitScore) = @$hit;
        my (undef,$locusId) = split /:/, $subjectId;
        my $gene = $subjectToGene{$locusId} || die "Unknown locusId $locusId in $orgId";
        $bitScore =~ s/ +//;
        my $cov = int(0.5 + 100*abs($queryEnd - $queryStart + 1)/length($seq));
        $percIdentity = int(0.5 + $percIdentity);
        my ($fitstring, $fittitle) = Utils::gene_fit_string($dbh, $orgId, $locusId);
        my $showId = $gene->{sysName} || $gene->{locusId};
        my $aln_URL = "showAlign.cgi?query=$defURI&querySequence=$seq&subject=${orgId}:${locusId}";
        my $slen = $subjectLen{$locusId} || die $locusId;
        my $cov_title = "amino acids $queryStart:$queryEnd / $seqlen of query"
          . " are similar to a.a. $subjectStart:$subjectEnd / $slen of $showId";
        my $covShow = defined $aln_URL ? a({ title => $cov_title, href => $aln_URL }, $cov)
          : a({ title => $cov_title }, $cov);
        my @row = (Utils::gene_link($dbh, $gene, "name", "geneOverview.cgi"),
                   $gene->{gene},
                   Utils::gene_link($dbh, $gene, "desc", "domains.cgi"),
                   $cgi->a({href => "myFitShow.cgi?orgId=$orgId&gene=$locusId", title => $fittitle }, $fitstring ),
                   $cgi->a({title=>"$bitScore bits", style=>"color:black;"},$percIdentity),
                   $covShow,
                   $cgi->a({title=>"$bitScore bits", -style=>"color:black;"}, $eVal));
        my @td = map td($_), @row;
        print $cgi->Tr({ -align => 'left', -valign => 'top', bgcolor=>'white' }, @td);
        print "\n";
      }
      print "</TABLE>\n";
    } # else not 6-frame
  } # else has hits
  print p("Or",
          a({-href => "mySeqSearch.cgi?query=>$defURI" . uri_escape("\n") . $seq},
            "BLAST against annotated proteins in", scalar(keys %$orginfo), "genomes")),
        "\n";
} else { # no query
  print
     header,
     Utils::start_page("BLAST against " . $orginfo->{$orgId}{genome}),
     div({-id=>"ntcontent"},
         h3("BLAST against", a({-href => "org.cgi?orgId=$orgId"}, $orginfo->{$orgId}{genome})),
         start_form( -class => 'search', -name    => 'input', -method  => 'GET', -action  => 'blastOrg.cgi' ),
         hidden('orgId'),
         p("Enter a protein sequence in FASTA or Uniprot format: ",
           br(),
           textarea( -name  => 'query', -value => '', -cols  => 68, -rows  => 10 )),
         p("Search against",
           popup_menu( -name => 'sixframe', -values => { 0 => 'annotated proteins', 1 => 'six-frame translation'},
                       '-default' => 0)),
         p(submit("Search")),
         end_form);
}
$dbh->disconnect();
Utils::endHtml($cgi);
