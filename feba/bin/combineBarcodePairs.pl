#!/usr/bin/perl -w
use strict;

my $usage = <<END
Usage: combinePairs.pl pairs1.tsv ... pairsN.tsv > combined.tsv

combineBarcodePairs.pl combines multiple outputs of barcodePairs.pl
into one table. Both input and output have the fields barcode1,
barcode2, and n, with a header row.
END
;

die $usage if @ARGV == 0 || $ARGV[0] =~ m/^-/;

my @inFiles = @ARGV;
foreach my $file (@inFiles) {
  die "No such file: $file\nRun with no arguments for usage\n" unless -e $file;
}

my %n = (); # bc1 to bc2 to n
foreach my $file (@inFiles) {
  open(my $fh, "<", $file) || die "Cannot read $file\n";
  while(my $line = <$fh>) {
    chomp $line;
    my @F = split /\t/, $line;
    die "Should have 3 columns" unless @F == 3;
    my ($bc1, $bc2, $n) = @F;
    next if $bc1 eq "barcode1";
    $n{$bc1}{$bc2} += $n;
  }
  close($fh) || die "Error reading $file\n";
}

my $nPairs = 0;
foreach my $barcode1 (sort keys %n) {
  my $hash = $n{$barcode1};
  foreach my $barcode2 (sort keys %$hash) {
    print join("\t", $barcode1, $barcode2, $hash->{$barcode2})."\n";
    $nPairs++;
  }
}
print STDERR "Wrote $nPairs pairs of barcodes\n";
