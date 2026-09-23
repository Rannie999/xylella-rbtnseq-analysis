#!/usr/bin/perl -w
# Given two annotationss for the same or close-related genomes, uses uesarch to line them up
use strict;
use Getopt::Long;
use FindBin qw{$RealBin};
use List::Util qw{min};
use lib "$RealBin/../lib";
use FEBA_Utils qw{ReadFastaEntry};

my $minIdentity = 0.98;
my $minCoverage = 0.8;
my $usearch = "$RealBin/usearch";

my $usage = <<END
matchProt.pl -in proteins1.faa proteins2.faa > matches.tsv
  Finds matching proteins from two different annotations of
  the same genome or closely-related genomes.

  The input files should be protein sequences in fasta format.

  The output table has four columns: id1, id2, identity, and coverage
  (the lower of the two values), where identity and coverage are
  fractions.

Optional arguments:
-identity $minIdentity (minimum fraction, 0.5 or higher)
-coverage $minCoverage (minimum fraction)
-usearch usearch_executable
   default: $usearch
END
;

my @inFiles;
die $usage
  unless GetOptions('in=s{2,2}' => \@inFiles,
                    'identity=f' => \$minIdentity,
                    'coverage=f' => \$minCoverage,
                    'usearch=s' => \$usearch)
  && @ARGV == 0;

die $usage unless @inFiles == 2;
die "No such executable: $usearch\n" unless -x $usearch;
die "identity threshold must be between 0.5 and 1\n"
  unless $minIdentity >= 0.5 && $minIdentity <= 1;
die "coverage threshould must be between 0 and 1\n"
  unless $minCoverage >= 0 && $minCoverage <= 1;

my $tmpDir = $ENV{TMPDIR} || "/tmp";
die "No such directory, change TMPDIR environmental variable: $tmpDir\n"
  unless -d $tmpDir;
my $tmpPre = "$tmpDir/matchProt.$$";

my @lenHash = ();
my @tmpFaaFiles = ();
foreach my $i (0,1) {
  my %len = ();
  my $inFile = $inFiles[$i];
  my $faaFile = "$tmpPre.$i.faa";
  push @tmpFaaFiles, $faaFile;
  open(my $fh, "<", $inFile) || die "Cannot read $inFile\n";
  open(my $fhOut, ">", "$faaFile") || die "Cannot write to $faaFile";
  my $state = {};
  while (my ($id, $seq) = ReadFastaEntry($fh, $state)) {
    $id =~ s/ .*//; # ignore any description
    die "Duplicate sequence for $id in $inFile\n" if exists $len{$id};
    $seq =~ s/[*]//g;
    die "Invalid sequence for $id\n" unless $seq =~ m/^[A-Z]+$/;
    $len{$id} = length($seq);
    print $fhOut ">$id\n$seq\n";
  }
  close($fh) || die "Error reading $inFile";
  close($fhOut) || die "Error writing to $faaFile";
  push @lenHash, \%len;
}

my ($lenH1, $lenH2) = @lenHash;
my $matchFile = "$tmpPre.blast6out";
my @cmd = ("$usearch", "-usearch_global",
           $tmpFaaFiles[0], "-db", $tmpFaaFiles[1],
           "-id", $minIdentity,
           "-maxaccepts", 20,
           "-maxrejects", 20,
           "-query_cov", $minCoverage,
           "-target_cov", $minCoverage,
           "-blast6out", $matchFile,
           "-quiet");
print STDERR "Running " . join(" ", @cmd)."\n";
system(@cmd) == 0 || die "Error running usearch: $!\n";

unlink($tmpFaaFiles[0]);
unlink($tmpFaaFiles[1]);

my %hits = (); # query to list of [score, subject, identity fraction, coverage]
my %revHits = (); # subject to list of [score, query, identity fraction, coverage]
open(my $fh, "<", $matchFile) || die "Cannot read $matchFile\n";
while (my $line = <$fh>) {
  chomp $line;
  my ($query, $subject, $identity, $alen, $mm, $gap, $qbeg, $qend, $sbeg, $send, $eval, $score)
     = split /\t/, $line;
  die "Invalid blast6out line from usearch\n$line\n" unless defined $score && $score ne "";
  die "Unknown query $query" unless exists $lenH1->{$query};
  die "Unknown subject $subject" unless exists $lenH2->{$subject};
  my $identityFraction = $identity / 100;
  my $cov = min(($qend-$qbeg+1) / $lenH1->{$query},
                ($send-$sbeg+1) / $lenH2->{$subject});
  my $scoreUse = $identityFraction * (($qend-$qbeg+1) / $lenH1->{$query});
  push @{ $hits{$query} }, [ $scoreUse, $subject, $identityFraction, $cov ];
  push @{ $revHits{$subject} }, [ $scoreUse, $query, $identityFraction, $cov ];
}
close($fh) || die "Error reading $matchFile";

unlink($matchFile);

print join("\t", "query", "subject", "identity", "coverage")."\n";
my $nOut  = 0;
foreach my $query (sort keys %hits) {
  my @hits = sort { $b->[0] <=> $a->[0] } @{ $hits{$query} };
  my ($score, $subject, $identity, $cov) = @{ $hits[0] };
  my @revHits = sort { $b->[0] <=> $a->[0] } @{ $revHits{$subject} };
  if (@revHits > 0 && $revHits[0][1] eq $query) {
    print join("\t", $query, $subject, $identity, $cov)."\n";
    $nOut++;
  }
}
print STDERR sprintf("Input: %d and %d sequences; %d pairs found\n",
                     scalar(keys %$lenH1), scalar(keys %$lenH2), $nOut);
