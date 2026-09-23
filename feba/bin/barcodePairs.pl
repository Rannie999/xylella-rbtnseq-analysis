#!/usr/bin/perl -w
# Identify pairs of barcodes from a fastq file
use strict;
use Getopt::Long;
use FindBin qw{$RealBin};
use lib "$RealBin/../lib";
use FEBA_Utils qw{reverseComplement};
sub findPreSeq($$$);

my $nPreExpected = 9;
my $preSeq1 = "GAGGTCTCT";
my $postSeq1 = "TACCGTTC";
my $nBetweenExpected = 17;
my $preSeq2 = "CGAACGGTA";
my $postSeq2 = "CGAGCTCG";
my $wobble = 2;
my $bcLen = 20;

my $usage = <<END
Usage: barcodePairs.pl -in fastq.gz > pairs.tsv
  Each barcodes is expected to be 20 nucleotides.
  If the input file suffix ends with .gz, it is
  decompressed with zcat.

Optional arguments:
-pre1 $preSeq1 -- constant sequence before the first barcode
-post1 $postSeq1 -- constant sequence after the first barcode
-pre2 $preSeq2 -- constant sequence before the second barcode
-post2 $postSeq2 -- constant sequence after the second barcode
-nPre $nPreExpected -- number of nucleotides before the 1st sequence
-between $nBetweenExpected -- number of nucleotides between
  post1 and pre2
-wobble $wobble -- amount of wobble to allow in position of either
                   barcode
-bcLen $bcLen -- length of each barcode
-rc -- reverse complement the reads before processing
-limit N -- process just the first N reads
END
;

my ($inFile, $nLimit, $rc, $debug);
die $usage
  unless GetOptions('inFile=s' => \$inFile,
                    'nPre=s' => \$nPreExpected,
                    'pre1=s' => \$preSeq1,
                    'post1=s' => \$postSeq1,
                    'between=i' => \$nBetweenExpected,
                    'pre2=s' => \$preSeq2,
                    'post2=s' => \$postSeq2,
                    'wobble=i' => \$wobble,
                    'bcLen=i' => \$bcLen,
                    'limit=i' => \$nLimit,
                    'rc' => \$rc,
                    'debug' => \$debug)
  && @ARGV == 0
  && defined $inFile;
die "Invalid -wobble: must be nonnegative\n" if $wobble < 0;
die "Invalid -nPre: must be nonnegative\n" if $nPreExpected < 0;
die "Invalid pre1\n" unless $preSeq1 =~ m/^[ACGT]+$/;
die "Invalid post1\n" unless $postSeq1 =~ m/^[ACGT]+$/;
die "Invalid pre2\n" unless $preSeq2 =~ m/^[ACGT]+$/;
die "Invalid post2\n" unless $postSeq2 =~ m/^[ACGT]+$/;
die "Invalid bcLen\n" unless $bcLen > 0;
$nLimit = undef if defined $nLimit && $nLimit < 1;

die "No such file: $inFile\n" unless -e $inFile;
my $fh;
my $pipe;
if ($inFile =~ m/[.gz]$/) {
  open($fh, "zcat $inFile |") || die "Cannot run zcat on $inFile\n";
  $pipe = 1;
} else {
  open ($fh, "<", $inFile) || die "Cannot read $inFile\n";
}

my $nRead = 0;
my $nUsed = 0;
my %n = (); # barcode1 => barcode2 => count
while (my $read = <$fh>) {
  chomp $read;
  $read =~ m/^@/ || die "Not fastq input: $read\n";
  my $seq = <$fh> || die "Truncated fastq\n";
  chomp $seq;
  $seq =~ m/^[ACGTN]+$/ || die "Invalid characters in fastq sequence (must be ACGTN)\n";
  my $line = <$fh> || die "Truncated fastq\n";
  chomp $line;
  $line eq "+" || die "Not fastq: line after sequence is $line\n";
  my $quality = <$fh> || die "Truncated fastq\n";
  chomp $quality;
  die "quality string is wrong length\n" unless length($quality) == length($seq);
  $nRead++;

  $seq = reverseComplement($seq) if defined $rc;
  my $pos1 = findPreSeq($seq, $preSeq1, $nPreExpected);
  if (defined $pos1) {
    my $start1 = $pos1 + length($preSeq1);
    if (length($seq) >= $start1 + $bcLen + length($postSeq1)) {
      my $postSeen1 = substr($seq, $start1 + $bcLen, length($postSeq1));
      if ($postSeen1 eq $postSeq1) {
        my $barcode1 = substr($seq, $start1, $bcLen);
        my $pos2 = findPreSeq($seq, $preSeq2, $start1 + $bcLen + length($postSeq1) + $nBetweenExpected);
        if (defined $pos2) {
          my $start2 = $pos2 + length($preSeq2);
          if (length($seq) >= $start2 + $bcLen + length($preSeq2)) {
            my $postSeen2 = substr($seq, $start2 + $bcLen, length($postSeq2));
            if ($postSeen2 eq $postSeq2) {
              my $barcode2 = substr($seq, $start2, $bcLen);
              print STDERR join("\t", $barcode1, $barcode2, $read)."\n"
                if defined $debug;
              $n{$barcode1}{$barcode2}++;
              $nUsed++;
            }
          }
        }
      }
    }
  }
  last if defined $nLimit && $nRead >= $nLimit;
}
if (defined $inFile) {
  close($fh) || $pipe || die "Error reading $inFile\n";
}
print STDERR "Read $nRead reads, used $nUsed\n";
print join("\t", qw{barcode1 barcode2 n})."\n";
foreach my $barcode1 (sort keys %n) {
  my $hash = $n{$barcode1};
  foreach my $barcode2 (sort keys %$hash) {
    print join("\t", $barcode1, $barcode2, $hash->{$barcode2})."\n";
  }
}

# Looks for a sbsequence matching $preSeq starting around $at.
# Returns the 0-based index of the match, or undef
sub findPreSeq($$$) {
  my ($seq, $preSeq, $at) = @_;
  my @off = 0;
  for (my $w = 1; $w <= $wobble; $w++) {
    push @off, (-$w, $w);
  }
  foreach my $off (@off) {
    my $i = $at + $off;
    return $i
      if $i >= 0
        && length($seq) >= length($preSeq) + $i
        && substr($seq, $i, length($preSeq)) eq $preSeq;
  }
  return undef;
}
