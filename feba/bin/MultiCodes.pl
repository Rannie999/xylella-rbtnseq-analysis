#!/usr/bin/perl -w
use strict;
use Getopt::Long;
use FindBin qw($Bin);
use lib "$Bin/../lib";
use FEBA_Utils qw{ReadTable};

my $minQuality = 10;
my $preseqOff1 = 0;

my $usage = <<END
Usage: MultiCodes.pl -out out_prefix -index name < fastq
Or, if not yet multiplexed,
       MultiCodes.pl -out out_prefix -primers PrimerIndexTable < fastq
Optional arguments:
    [ -debug ] [ -dntag ] [ -limit maxReads ] [ -minQuality $minQuality ]
    [ -n25 | -bs3 | -bs4 ]
    [ -preseq CAGCGTACG -postseq AGAGACCTC -nPreExpected 9 ]
    [ -index2 GCACTA -at2 2 ] [ -preseqOff1 $preseqOff1 ]
    [ -barcodeLength 20 ]

    -n25 or -bs3 indicate newer multiplexing designs:
    -n25 means 11:14 nt before the pre-sequence, corresponding to a read with
	2:5 Ns, GTCGACCTGCAGCGTACG, N20, AGAGACC
    -bs3 means 1:4 + 6 + 11 = 18:21 nt before the pre-sequence, corresponding to
	1:4 Ns, index2, GTCGACCTGCAGCGTACG, N20, AGAGACC
	where nN and index2 is specified in ../primers/barseq3.index2
    -bs4 means 1:4 + 8 + 11 = 20:23 before the pre-sequence, corresponding to
	1:4 Ns, index2, GTCGACCTGCAGCGTACG, N20, AGAGACC
	where nN and index2 is specified in ../primers/barseq4.index2
    -nPreExpected can also be a range, i.e. 11:14
    -index2 and -at2 can check for a specific sequence somewhere in the read
    -preseqOff1 $preseqOff1 -- if non-zero, allow a 1-nt error in the pre-sequence

    PrimerIndexTable should be tab-delimited with multiplex name and a primer like nACGACG
    The fastq file should be fastq with phred+33 ("sanger") encoding of quality scores
    (as in MiSeq or 2012+ HiSeq)

This script analyzes multiplexed random barcode sequences, such as
nATCACGAG GTCGACCTGCAGCGTACG 20N AGAGACCTCGTGGACATCAGATC
for primer nATCACGAG (nLeading=1, indexseq=ATCACGAG), where the 20 Ns are the barcode

or, if -dntag is specified,
nATCACGAG CGGTGTCGGTCTCGTAG 20N CGATGAATTCGAGCTCGTT

or, if the index is specified, there is no multiplexing, and the read is
nnnnn GTCGACCTGCAGCGTACG 20N AGAGACC (where the leading ns are random and are ignored)

or, if -preseq -postseq -nPreExpected are all used, it expects a read of the form
   nPreExpected characters (any, although this includes the multiplexing tag unless -index is used)
   preseq
   barcode
   postseq

For a barcode to be counted, the preseq (9 nt upstream of the barcode)
much match exactly; the postseq (9 nt downstream of the barcode) much
also be correct unless minQuality = 0 (in which case the post-sequence
is not checked and there is no guarantee that the barcode is the
correct length).

If -minQuality > 0, then it also tries to estimate the diversity of the sample,
and the out_prefix.good file lists the "good" barcodes (not off-by-1 from
a more abundant barcode and with at least two reads). The "match10"
field shows if the first 10 or last 10 nucleotides of the barcode match
to a much-more abundant barcode; these are likely chimeras.
END
    ;

sub Variants($); # return all 1-nt variants for a sequence
sub FindPrefix($$); # sequence, list of index prefixes => index of prefix or -1
sub FindBarcode($$$); # sequence, quality, first post-prefix index => (spacingOk, barcode) or undef
sub sum(@); # sum of a list of values

my $preseq = undef; # before the barcode
my %preseqAccept = (); # accepted sequencsE (may include 1-offs)
my $postseq = undef; # after (for 50 nt reads this needs to be shortened to AGA)
my $nPreExpected; # nt between prefix and preseq
my $nPreExpectedMin;
my $nPreExpectedMax;
my $debug = 0;
my $dntag = 0;
my $iname = undef;
my $doOff1 = undef; # check for off-by-1 barcodes
my $barcodeLength = 20;
{
    my ($indexfile,$out,$nLimit,$n25,$bs3,$bs4);
    my ($index2,$index2At); # check that sequence has $index2 starting at (1-based) $index2At
    (GetOptions('primers=s' => \$indexfile,
		'out=s' => \$out,
		'index=s' => \$iname,
                'minQuality=i' => \$minQuality,
                'limit=i' => \$nLimit,
                'nPreExpected=s' => \$nPreExpected,
		'preseq=s' => \$preseq,
		'postseq=s' => \$postseq,
                'dntag' => \$dntag,
                'debug' => \$debug,
                'n25' => \$n25,
                'bs3' => \$bs3,
                'bs4' => \$bs4,
                'index2=s' => \$index2,
                'at2=i' => \$index2At,
                'off1=i' => \$doOff1,
                'preseqOff1=i' => \$preseqOff1,
                'barcodeLength=i' => \$barcodeLength)
     && defined $out)
        || die $usage;
    die $usage unless (defined $indexfile xor defined $iname);
    die "Invalid barcode length\n" unless $barcodeLength > 1;
    $doOff1 = $minQuality > 0 unless defined $doOff1;

    if (defined $preseq) {
	die "Missing -postseq: $usage" unless defined $postseq;
	die "Missing -nPreExpected: $usage" unless defined $nPreExpected;
    } else {
	if ($dntag) {
	    die "-index with -dntag not supported" if defined $iname;
	    $preseq = "GTCTCGTAG";
	    $nPreExpected = 8 unless defined $nPreExpected;
	    $postseq = "CGATGAATT";
	} else {
	    $preseq = "CAGCGTACG";
	    $nPreExpected = (defined $iname ? 14 : 9)
              unless defined $nPreExpected;
	    $postseq = "AGAGACCTC";
	}
    }

    die "Cannot specify multiple protocols\n"
      if (defined $n25) + (defined $bs3) > 1;
    if (defined $n25) {
      die "Cannot specify -n25 unless -index is set\n" unless defined $iname;
      $nPreExpected = '11:14';
    }

    if (defined $bs3 || defined $bs4) {
      die "Cannot specify -bs3 unless -index is set\n" unless defined $iname;
      my $ifile = defined $bs3 ? "$Bin/../primers/barseq3.index2"
        : "$Bin/../primers/barseq4.index2";
      die "No such file: $ifile\n" unless -e $ifile;
      my @tab = FEBA_Utils::ReadTable($ifile, ["index_name","index2","nN"]);
      my @tabMatch = grep { $_->{index_name} eq $iname } @tab;
      if (@tabMatch == 0) {
        print STDERR "Warning! Ignoring the second index -- index_name $iname does not appear in $ifile\n";
        $nPreExpected = '16:19';
      } else {
        die "Multiple matches for index $iname in $ifile" if @tabMatch > 1;
        my $row = $tabMatch[0];
        $index2 = $row->{index2};
        $index2At = $row->{nN}+1;
        $nPreExpected = $index2At + length($index2) + 9;
      }
    }

    if ($nPreExpected =~ m/^(\d+):(\d+)$/) {
        $nPreExpectedMin = $1;
        $nPreExpectedMax = $2;
    } else {
        die $usage unless $nPreExpected =~ m/^\d+$/;
        $nPreExpectedMin = $nPreExpected - 2;
        $nPreExpectedMax = $nPreExpected + 2;
    }

    my $index2len;
    if (defined $index2) {
      die "Must specify -at2 with -index2\n" unless defined $index2At;
      die "Invalid index2 sequence $index2" unless $index2 =~ m/^[ACGT]+$/;
      die "Invalid $index2At" unless $index2At =~ m/^\d+$/ && $index2At >= 1;
      $index2len = length($index2);
      print STDERR "Checking for index2 $index2 at position $index2At\n";
    }

    if ($minQuality > 0 && defined $bs4 && length($postseq) > 3) {
      $postseq = substr($postseq, 0, 3);
      print STDERR "Only checking for the 1st 3 nt of the postseq: $postseq\n";
    }

    my $nReads = 0;
    my $nMulti = 0; # number with prefix identified
    my $nWithBarcode = 0; # number with a (quality) barocde identified
    my %codes = (); # barcode => number of times barcode is seen, with an array of one entry per multiplex tag
    my @prefix = (); # list of [ nLeading, indexseq (no leading Ns), name ]

    my @prefixNames = ();
    if (defined $indexfile) {
	open(INDEX, "<", $indexfile) || die "Error reading $indexfile";
	while (my $line = <INDEX>) {
	    chomp $line;
	    my ($name,$index) = split /\t/, $line;
	    next if $name =~ m/name/i;
	    my $nLeading;
	    my $indexseq;
	    if ($name eq "none" && $index eq "") {
		$nLeading = 0;
		$indexseq = "";
	    } else {
		die "Invalid index sequence $index" unless $index =~ m/^([nN]*)([ACGT]+)$/;
		$nLeading = length($1);
		$indexseq = $2;
	    }
	    die $line if !defined $nLeading || !defined $indexseq || !defined $name;
	    push @prefix, [ $nLeading, $indexseq, $name ];
	}
	close(INDEX) || die "Error reading $indexfile";
	print STDERR "Read " . scalar(@prefix) . " indices from $indexfile\n";
	@prefixNames = map { $_->[2] } @prefix;
    }
    if (defined $iname) {
	@prefix = ( [0, "", $iname] );
	@prefixNames = ( $iname );
    }

    # All 1-nt errors of the pre-sequence
    %preseqAccept = ( $preseq => 1 );
    if ($preseqOff1 != 0) {
      foreach my $variant (Variants($preseq)) {
        $preseqAccept{$variant} = 1;
      }
    }

    my $nWrongIndex2 = 0;
    while(my $header = <STDIN>) {
      chomp $header;
      my $seq = <STDIN>;
      chomp $seq;
      $nReads++;
      <STDIN>;
      my $quality = <STDIN>;
      last if defined $nLimit && $nReads >= $nLimit;

      if (defined $index2) {
        my $part = substr($seq, $index2At-1, $index2len);
        if ($part ne $index2) {
          $nWrongIndex2++;
          next;
        }
      }

      my $iPrefix = defined $iname ? 0 : FindPrefix($seq, \@prefix);
      next if $iPrefix < 0;
      $nMulti++;
      my ($nLeading, $indexseq, $prefixName) = @{ $prefix[$iPrefix] };
      my ($barcode,$okSpacing) = FindBarcode($seq, $quality, $nLeading+length($indexseq));
      next unless defined $barcode;
      $nWithBarcode++;
      $codes{$barcode} = [] unless exists $codes{$barcode};
      $codes{$barcode}[$iPrefix]++;
    }

    my $nUniq = scalar(keys %codes);
    print STDERR "Reads $nReads Multiplexed $nMulti Usable($barcodeLength) $nWithBarcode "
      . sprintf("(%.1f%%)", 100 * $nWithBarcode/$nReads) . " unique codes $nUniq \n" if $nReads > 0;
    print STDERR sprintf("Failed to match index2: %d reads (%.3f%%)\n", $nWrongIndex2, 100*$nWrongIndex2/$nReads)
      if $nReads > 0 && defined $index2;

    my %nPerCount = (); # number of codes with that count
    open(CODES, ">", "$out.codes.tmp") || die "Cannot write to $out.codes.tmp";
    print CODES join("\t","barcode", @prefixNames)."\n";
    while (my ($code,$count) = each %codes) {
        my @counts = ();
        my $nPrefix = scalar(@prefix);
        for (my $i=0; $i < $nPrefix; $i++) {
            push @counts, $count->[$i] || 0;
        }
        print CODES join("\t",$code,@counts)."\n";
        my $nAll = sum(@counts);
        $nPerCount{$nAll}++;
    }
    close(CODES) || die "Error writing to $out.codes.tmp";
    rename("$out.codes.tmp","$out.codes") || die "Cannot rename to $out.codes";
    print STDERR "Wrote " . scalar(keys %codes) . " unique barcodes to $out.codes\n";

    open(COUNTS, ">", "$out.counts") || die "Cannot write to $out.counts";
    print COUNTS join("\t", "Count", "nCodes", "Frac")."\n";
    foreach my $count (sort {$a<=>$b} keys %nPerCount) {
        print COUNTS join("\t",$count,$nPerCount{$count},$nPerCount{$count}/$nUniq)."\n";
    }
    close(COUNTS) || die "Error writing to $out.counts";
    print STDERR "Wrote number of barcodes seen a certain number of times to $out.counts";
    print STDERR "; nOnce = $nPerCount{1}" if (exists $nPerCount{1});
    print STDERR "\n";

    # and look for off-by-1 cases
    my %offby1 = (); # barcode => 1 for likely off-by-1 errors
    if ($doOff1) {
      open(CLOSE, ">", "$out.close") || die "Cannot write to $out.close";
      print CLOSE join("\t",qw{code1 count1 code2 count2})."\n";
      my $nCases = 0;
      my $nOff1Reads = 0;

      while (my ($code,$count) = each %codes) {
        my @variants = Variants($code);
        foreach my $variant (@variants) {
          if (($code cmp $variant) > 0 && exists $codes{$variant}) {
            my $n1 = sum(@$count);
            my $n2 = sum(@{$codes{$variant}});
            print CLOSE join("\t",$code,$n1,$variant,$n2)."\n";
            $offby1{ $n1 < $n2 ? $code : $variant } = 1;
            $nCases++;
            $nOff1Reads += $n1 < $n2 ? $n1 : $n2;
          }
        }
      }
      close(CLOSE) || die "Error writing to $out.close";
      my $fOff1 = sprintf("%.3f", $nOff1Reads / ($nWithBarcode || 1));
      print STDERR "Wrote $nCases off-by-1 pairs ($nOff1Reads reads, fraction $fOff1) to $out.close\n";

      # And write out good barcodes with >1 count, in order of abundance
      my %goodCount = (); # barcode to total
      while (my ($code, $count) = each %codes) {
        if (!exists $offby1{$code}) {
          my $tot = sum(@$count);
          $goodCount{$code} = $tot if $tot > 1;
        }
      }

      # If the 1st 10 or last 10 nt match, is potentially a chimera or other artefact
      my %half1 = ();
      my %half2 = ();
      foreach my $code (keys %goodCount) {
        my $half1 = substr($code, 0, 10);
        push @{ $half1{$half1} }, $code;
        my $half2 = substr($code, 10, 10);
        push @{ $half2{$half2} }, $code;
      }

      open(my $fhGood, ">", "$out.good") || die "Cannot write to $out.good";
      print $fhGood join("\t", "barcode", "n", "match10")."\n";
      foreach my $code (sort { $goodCount{$b} <=> $goodCount{$a} || $a cmp $b } (keys %goodCount)) {
        my $match1 = $half1{ substr($code, 0, 10) };
        my $match2 = $half2{ substr($code, 10, 10) };
        my $matchTo = "";
        foreach my $sim (@$match1, @$match2) {
          $matchTo = $sim
            if $goodCount{$sim} >= 100 * $goodCount{$code};
        }
        print $fhGood join("\t", $code, $goodCount{$code}, $matchTo)."\n";
      }
      close($fhGood);
      print STDERR "Wrote $out.good\n";
    }

    # and estimate diversity
    if ($minQuality > 0 && $nWithBarcode >= 1000 && $nPerCount{2} >= 10) {
	foreach my $fNoise (0, 0.005, 0.01, 0.02) {
	    my $nNoise = int(0.5 + $nWithBarcode * $fNoise);
	    next if $nNoise >= $nPerCount{1};
	    print STDERR sprintf("If %.1f%% of reads are noise: diversity %.1f K from total barcodes %.1f K seen once %.1f K seen twice %.1f K\n",
				 $fNoise*100,
				 ( $nUniq - $nNoise + ($nPerCount{1} - $nNoise)**2 / (2 * $nPerCount{2}) )/1000.0,
				 ($nUniq-$nNoise)/1000.0, ($nPerCount{1} - $nNoise)/1000.0, $nPerCount{2}/1000.0);
	}
    }
    if ($minQuality > 0 && $nWithBarcode >= 1000 && $doOff1) {
        my $nGoodCodes = 0;
        my $nGoodReads = 0;
        while (my ($code,$count) = each %codes) {
            my $tot = sum(@$count);
            if (!exists $offby1{$code} && $tot > 1) {
                $nGoodCodes++;
                $nGoodReads += $tot;
            }
        }
        print STDERR sprintf("Aside from singletons and off-by-1s, see %.1f K barcodes (%.1f%% of reads)\n",
                             $nGoodCodes/1000.0, $nGoodReads * 100.0 / $nWithBarcode);
    }
    # and estimate bias
    if ($minQuality > 0 && $nUniq >= 5000) {
	# What fraction of reads are accounted for by the top 1% of strains?
	my $f = 0.01;
	my $nReadsSofar = 0;
	my $nCodesSofar = 0;
	my $countSofar = -1;
	foreach my $count (sort {$b <=> $a} keys(%nPerCount)) {
	    $nReadsSofar += $count * $nPerCount{$count};
	    $nCodesSofar += $nPerCount{$count};
	    $countSofar = $count;
	    last if $nCodesSofar >= $f * $nUniq;
	}
	print STDERR sprintf("Barcodes with >= %d reads each: %.2f%% of codes (%.2f K), %.2f%% of reads (%.1f K)\n",
                             $countSofar,
                             (100 * $nCodesSofar)/$nUniq, $nCodesSofar/1000,
                             (100.0 * $nReadsSofar)/$nWithBarcode, $nReadsSofar/1000);
    }
}

sub FindPrefix($$) {
    my ($seq, $indexes) = @_;
    my @matches = ();
    for (my $i = 0; $i < scalar(@$indexes); $i++) {
        my ($nLeading,$indexseq,$name) = @{ $indexes->[$i] };
        if (substr($seq, $nLeading, length($indexseq)) eq $indexseq) {
            push @matches, $i;
        }
    }
    return $matches[0] if @matches == 1;
    print STDERR "No prefix for $seq\n" if $debug;
    return -1;
}

# returns barcode and if spacing to pre-sequence is correct, or undef
sub FindBarcode($$$) {
  my ($seq,$quality,$offset) = @_;
  my $seq2 = substr($seq, $offset);
  my $quality2 = substr($quality, $offset);

  my $prepos;
  for ($prepos = $nPreExpectedMin; $prepos <= $nPreExpectedMax; $prepos++) {
    my $preseqCandidate = substr($seq2, $prepos, length($preseq));
    if (exists $preseqAccept{$preseqCandidate}) {
      last;
    }
  }
  unless($prepos <= $nPreExpectedMax) {
    print STDERR "seq2 $seq2 has invalid index-of-preseq $prepos\n" if $debug;
    return undef; # preseq not found
    # Older code returned undef,0 to indicate that preseq was at the wrong location,
    # but this is no longer supported
  }
  my $barcode = substr($seq2, $prepos+length($preseq), $barcodeLength);
  unless(length($barcode) == $barcodeLength && $barcode =~ m/^[A-Z]+$/) {
    # Incorrect length should be rare but is possible if read is too short.
    # Barcodes with ambiguous sequences are ignored
    print STDERR "seq2 $seq2 has invalid barcode sequence $barcode\n" if $debug;
    return undef;
  }

  if ($minQuality > 0) {
    # check postSeq
    my $seqAfter = substr($seq2, $prepos + length($preseq) + $barcodeLength, length($postseq));
    my $postSeqShort = substr($postseq, 0, length($seqAfter));
    if ($seqAfter ne $postSeqShort) {
      print "seq $seq2 barcode $barcode postseq not found\n" if $debug;
      return undef;
    }

    # check barcode quality
    # the sanger encoding is !"#$%&'()*+,-./0123456789:;<=>?@ABCDEFGHIJ
    # so the quality score is encoded as 33+
    my $barqual = substr($quality2, $prepos+length($preseq), $barcodeLength);
    my @scores = map { $_ - 33 } unpack("%C"x$barcodeLength, $barqual);
    foreach my $score (@scores) {
      die "Invalid score $score from barcode $barcode quality $barqual" if $score < 0 || $score > 100;
      if ($score < $minQuality) {
          print STDERR "Low quality $score for barcode $barcode in $seq2\n" if $debug;
          return undef;
        }
    }
  }
  return($barcode, 1);
}

sub Variants($) {
    my ($baseseq) = @_;
    my @out = ();
    $baseseq = uc($baseseq);
    my $len = length($baseseq);
    foreach my $i (0..($len-1)) {
        my $pre = substr($baseseq,0,$i);
        my $char = substr($baseseq,$i,1);
        my $post = substr($baseseq,$i+1);
        next unless $char eq "A" || $char eq "C" || $char eq "G" || $char eq "T";
        foreach my $newchar (qw{A C G T}) {
            push @out, $pre . $newchar . $post unless $newchar eq $char;
        }
    }
    return(@out);
}

sub sum( @ ) {
    my @x = @_;
    my $sum = 0;
    foreach (@x) { $sum += $_ if defined $_; }
    return $sum;
}
