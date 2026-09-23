#!/usr/bin/perl -w
# Test routines for linking barseq read files to the correct index
use strict;
use Getopt::Long;
use FindBin '$Bin';
use lib "$Bin/../lib";
use FEBA_Utils qw{FastqByIndex FastqToIndex};

my $usage = <<END
Usage: fastqIndex.pl -dir directory_with_fastqs
       fastqIndex.pl -fq filename
       fastqIndex.pl -fqs file_listing_fqs

Outputs lines with the index and the file name

Optional arguments:
 -bywell -- expect well names like A01 instead of IT001
 -bynumber -- index S11 is coded by number (_11_)
 -offset 23 -- the sample number is offset by this much from the
   index number
END
;

my ($dir, $fq, $fqsFile);
my ($bywell, $bynumber, $offset);

die $usage unless GetOptions('dir=s' => \$dir,
                             'fq=s' => \$fq,
                             'fqs=s' => \$fqsFile,
                             'bywell' => \$bywell,
                             'bynumber' => \$bynumber,
                             'offset=i' => \$offset)
  && @ARGV == 0
  && (defined $dir xor defined $fq xor defined $fqsFile);

if (defined $dir) {
  die "No such directory: $dir\n" unless -d $dir;
  my $fastqByIndex = FastqByIndex('in' => $dir,
                                  'bywell' => defined $bywell,
                                  'bynumber' => defined $bynumber,
                                  'offset' => $offset || 0);
  foreach my $index (sort keys %$fastqByIndex) {
    foreach my $fq (@{ $fastqByIndex->{$index} }) {
      print join("\t", $index, $fq)."\n";
    }
  }
} else {
  my @fqs;
  if (defined $fq) {
    @fqs = ($fq);
  } else {
    if ($fqsFile eq "-") {
      @fqs = <STDIN>;
    } else {
      open(IN, "<", $fqsFile) || die "Cannot read $fqsFile\n";
      @fqs = <IN>;
      close(IN) || die "Error reading $fqsFile\n";
    }
    map chomp, @fqs;
  }
  foreach my $fq (@fqs) {
    my $fastq = FastqToIndex($fq, defined $bynumber, defined $bywell, $offset || 0);
    print join("\t", $fq, $fastq || "none")."\n";
  }
}
