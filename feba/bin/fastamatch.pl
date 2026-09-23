#!/usr/bin/perl -w
use strict;
use Getopt::Long;
use FindBin qw{$RealBin};
use lib "$RealBin/../lib";
use FEBA_Utils qw{ReadFastaEntry};
sub reverseComplement($);

my $usage = "Usage: fastamatch.pl [-strands] [-v] regexp < input.faa > subset.faa\n";
my $v = undef;
my $strands = undef;
die $usage
  unless GetOptions('v' => \$v, 'strands' => \$strands)
  && @ARGV == 1;
my ($regexp) = @ARGV;

my $state = {};
my $fh = *STDIN;
while (my ($header,$seq) = ReadFastaEntry($fh, $state)) {
  my $match = $seq =~ $regexp;
  if (defined $strands && ! $match) {
    $seq = reverseComplement($seq);
    $match = $seq =~ $regexp;
  }
  print ">$header\n$seq\n"
    if (defined $v ? ! $match : $match);
}

sub reverseComplement($)
{
    my $seq = shift;
    chomp $seq;
	my $origSeq=$seq;

    die "Invalid sequence \"$origSeq\" in reverseComplement" if ( ($seq =~ 
tr/RYKMSWBDHVNATCGXrykmswbdhvnatcg-/YRMKWSVHDBNTAGCXyrmkwsvhdbntagc-/) != 
length($seq) );
    $seq = reverse $seq;
    return $seq;
}
