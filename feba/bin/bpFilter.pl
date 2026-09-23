#!/usr/bin/perl -w
use strict;
use lib "/usr2/people/mprice/src/PaperBLAST/lib";
use pbutils qw{ReadTable};
sub commify($);

die <<END
Usage: bpFilter.pl pool1.withgenes pool2.withgenes bcpairs > subset
   Assumes that bcpairs is tab-delimited with the fields
   rcbarcode2 barcode1 n. If name ends with .gz it is gunzipped.
END
  unless @ARGV==3;

my ($pool1File, $pool2File, $pairFile) = @ARGV;
my @pool1 = ReadTable($pool1File, [qw{barcode rcbarcode scaffold strand pos locusId f}]);
my @pool2 = ReadTable($pool2File, [qw{barcode rcbarcode scaffold strand pos locusId f}]);
my %bc1 = map { $_->{barcode} => $_ } @pool1;
my %rbc2 = map { $_->{rcbarcode} => $_ } @pool2;

my $nPairsTot = 0;
my $nReadsTot = 0;
my $nPairsMatch = 0;
my $nReadsMatch = 0;

my $gz = $pairFile =~ m/[.]gz$/;
my $fh;
if ($gz) {
  open($fh, "zcat $pairFile |") || die "Cannot zcat $pairFile\n";
  $gz = 1;
} else {
  open($fh, "<", $pairFile) || die "Cannot read $pairFile\n";
}
print join("\t", qw{barcode1 scaffold1 strand1 pos1 locusId1 f1 barcode2 scaffold2 strand2 pos2 locusId2 f2 n})."\n";
while (my $line = <$fh>) {
  chomp $line;
  my ($rbc2, $bc1, $nReads) = split /\t/, $line;
  die "Invalid input: $line\n" unless defined $nReads && $nReads ne "";
  next if $nReads eq "n";
  die "Invalid reads $nReads" unless $nReads =~ m/^\d+$/;
  $nPairsTot++;
  $nReadsTot += $nReads;
  if (exists $bc1{$bc1} && exists $rbc2{$rbc2}) {
    $nPairsMatch++;
    $nReadsMatch += $nReads;
    my $p1 = $bc1{$bc1};
    my $p2 = $rbc2{$rbc2};
    print join("\t",
               $p1->{barcode}, $p1->{scaffold}, $p1->{strand}, $p1->{pos}, $p1->{locusId}, $p1->{f},
               $p2->{barcode}, $p2->{scaffold}, $p2->{strand}, $p2->{pos}, $p2->{locusId}, $p2->{f},
               $nReads)."\n";
  }
}
close($fh) || die "Error reading $pairFile\n";
print STDERR join(" ",
                  "Read", commify($nPairsTot), "pairs with", commify($nReadsTot), "reads.",
                  "Kept", commify($nPairsMatch), "pairs with", commify($nReadsMatch), "reads.")."\n";

sub commify($) {
    local $_  = shift;
    1 while s/^(-?\d+)(\d{3})/$1,$2/;
    return $_;
}
