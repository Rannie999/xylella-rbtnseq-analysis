#!/usr/bin/perl -w
use strict;
use Getopt::Long;
use FindBin qw{$RealBin};
use lib "$RealBin/../lib";
use FEBA_Utils qw{ReadFasta};

my $phobiusDir = "$RealBin/phobius";
my $usage = <<END
Usage: myPhobius.pl -in in.faa -out prefix
       myPhobius.pl -batch -in in.faa

The input fasta should contain just one sequence
unles run with -batch.

In batch mode, writes to prefix.raw and prefix.tsv
Otherwise, writes to prefix.plp, prefix.R, prefix.tsv, and prefix.svg

Options:
  -dir $phobiusDir -- the phobius directory
   It must contain the decodeanhmm executable,
   phobius.options, and phobius.model
END
;

my ($inFaa, $out, $batch);
die $usage unless GetOptions('in=s' => \$inFaa,
                             'out=s' => \$out,
                             'batch' => \$batch,
                             'dir=s' => \$phobiusDir)
  && @ARGV == 0
  && defined $inFaa && defined $out;

die "No such directory: $phobiusDir\n" unless -d $phobiusDir;
my $decode = "$RealBin/phobius/decodeanhmm";
die "No such executable: $decode\n" unless -x $decode;
my $optFile = "$RealBin/phobius/phobius.options";
my $modelFile = "$RealBin/phobius/phobius.model";
foreach my $file ($optFile, $modelFile) {
  die "No such file: $file\n" unless -e $file;
}

my $seqs = ReadFasta($inFaa);
die "$inFaa should have exactly one sequence\n" unless keys(%$seqs) == 1 || defined $batch;

my $rawFile = defined $batch ? "$out.raw" : "$out.plp.$$.tmp";

system("($decode -f $optFile $modelFile -plp < $inFaa > $rawFile) 2> /dev/null") == 0
  || die "phobius decodeanhmm failed: $!";

if (defined $batch) {
  my @stateNames = ("cytoplasmic", "non-cytoplasmic", "transmembrane", "signal peptide");
  my %state = (); # sequence => list of indexes of preferred states by position
  open (my $fhRaw, "<", $rawFile) || die "Cannot read $rawFile";
  my $name; # which protein we are reading results for
  while(my $line = <$fhRaw>) {
    chomp $line;
    if ($line =~ m/^# (\S+)$/) {
      $name = $1;
    } elsif ($line =~ m/^# +i +o +O/) {
      # skip
    } elsif ($line =~ m/^[A-Z] +/) {
      die "No sequence specifier before line" unless defined $name;
      my @F = split / +/, $line;
      my @p = ($F[1], $F[2]+$F[3], $F[4], $F[5]+$F[6]+$F[7]+$F[8]); # i, o, M, S
      my $i = 0;
      my $maxP = $p[0];
      my $maxI = 0;
      foreach my $i (1..3) {
        if ($p[$i] > $maxP) {
          $maxP = $p[$i];
          $maxI = $i;
        }
      }
      push @{ $state{$name} }, $maxI;
    }
  }
  close($fhRaw) || die "Error reading $rawFile\n";
  open(my $fhT, ">", "$out.tsv") || die "Cannot write to $out.tsv";
  foreach my $name (sort keys %state) {
    my $states = $state{$name};
    my $len = scalar(@$states);
    my @regions = (); # stored 0-based
    foreach my $i (0..($len-1)) {
      if ($i > 0 && $states->[$i] == $states->[$i-1]) {
        die unless $states->[$i] == $states->[ $regions[-1][0] ];
        $regions[-1][1] = $i;
      } else {
        push @regions, [$i, $i];
      }
    }
    foreach my $region (@regions) {
      my ($beg,$end) = @$region;
      print $fhT join("\t", $name, $beg+1, $end+1, $stateNames[ $states->[$beg] ])."\n";
    }
  }
  close($fhT) || die "Error writing $out.tsv\n";
  print STDERR "Wrote $out.raw\n";
  print STDERR "Wrote $out.tsv\n";
  exit(0);
}

# else
my ($seq) = values(%$seqs);
my $seqLen = length($seq);
my ($id) = keys(%$seqs);
$id =~ s/ .*//;
$id =~ s/["']//g;


# Convert the decodeanhmm output, which will look like
# # BT2157
#         i         o         O         M         n         h         c         C   
# M      0.00002   0.00001   0.00009   0.00000   0.99988   0.00000   0.00000   0.00000
# K      0.00002   0.00001   0.00009   0.00000   0.99988   0.00000   0.00000   0.00000
# K      0.00002   0.00001   0.00009   0.00000   0.99945   0.00043   0.00000   0.00000
# into gnuplot format
# # BT2157
# #pos	aa	i	o	M	S
# 1	M	0.00002	0.0001	0.00000	0.99988
# 2	K	0.00002	0.0001	0.00000	0.99988
# 3	K	0.00002	0.0001	0.00000	0.99988

open(my $fhRaw, "<", $rawFile) || die "Cannot read $rawFile";
open(my $fhP, ">", "$out.plp") || die "Cannot write to $out.plp";

my $header1 = <$fhRaw>;
print $fhP $header1;
my $header2 = <$fhRaw>;
die "Invalid line 2 from decodeanhmm: $header2" unless $header2 =~ m/# +i +/;
print $fhP join("\t", qw{pos aa i o M S})."\n";
my $iLine = 0;
while(my $line = <$fhRaw>) {
  $iLine++;
  chomp $line;
  my @F = split / +/, $line;
  die "Invalid data line from decodeanhmm: $line" unless @F == 9;
  my $S = 
  print $fhP join("\t", $iLine, $F[0], $F[1], $F[2]+$F[3], $F[4], $F[5]+$F[6]+$F[7]+$F[8])."\n";
}
close($fhRaw) || die "Error reading $rawFile";
unlink($rawFile);
close($fhP) || die "Error writing $out.plp";

open(my $fhR, ">", "$out.R")
  || die "Cannot write to $out.gnuplot";
my $Rcmd = qq{
id = "$id";
out = "$out";
};
$Rcmd .= q{

data = read.delim(paste0(out,".plp"), skip=1);
colors = c("red", "magenta", "blue", "green");
labels = c("transmembrane", "signal peptide", "cytoplasmic", "non-cytoplasmic");
columns = c("M", "S", "i", "o");

# compute sequential regions
class = apply(data[,columns], 1, which.max);
rle = rle(class);
d = data.frame(constEnd = cumsum(rle$lengths),
               class = rle$values);
d$constBegin = d$constEnd - rle$lengths + 1;

write.table(data.frame(begin=d$constBegin, end=d$constEnd, state=labels[d$class]),
            paste0(out,".tsv"),
            sep="\t", row.names=F, quote=F);

svg(paste0(out,".svg"), width=7, height=5);
par(mar=c(5,3,2.5,1), mgp=c(2,0.8,0), bty="l", xaxs="i", yaxs="i");
plot(c(1,max(data$pos)), c(0,1), col=NA,
  xlab="",
  ylab="Posterior probability", main=paste("Phobius analysis of\n", id));

tm = subset(d, class==1); # transmembrane regions
if (nrow(tm) > 0) rect(tm$constBegin, 0, tm$constEnd, 1, col="pink", border=NA);

for(i in 1:length(labels))
  lines(data$pos, data[[ columns[i] ]], col=colors[i], lwd=2, xpd=T);
legend(0, -4 * strheight("A"), labels, lty=1, col=colors, lwd=2, xpd=T, ncol=2, bty="n",
  fill=c("pink","white","white","white"), border=NA);

# arrows for sequential regions
y = -3.5 * strheight("A");
arrows(d$constBegin, y, d$constEnd, y,
  length=0.05, code=3, col=colors[d$class], xpd=T, lwd=2);
invisible( dev.off() );

};

print $fhR $Rcmd;
close($fhR) || die "Error writing to $out.R";

system("Rscript", "$out.R") == 0
  || die "R failed on $out.R -- $!";

