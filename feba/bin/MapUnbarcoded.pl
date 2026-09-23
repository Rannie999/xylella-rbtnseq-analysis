#!/usr/bin/perl -w
# Given TnSeq data from a library of transposon insertions without random barcodes, for each usable read,
# identify the barcode and the location in the genome.

use strict;
use Getopt::Long;
use FindBin qw($Bin);
use lib "$Bin/../lib";

my $wobbleAllowed = 4; # uncertainty in location of end of transposon
my $nExtra = 0;
my $tmpdir = defined $ENV{TMPDIR} ? $ENV{TMPDIR} : "/tmp";
my $minIdentity = 90; # minimum %identity for mapping to genome or past-end
my $minScore = 15; # minimum score for mapping to genome or past-end
my $delta = 5; # minimum difference in score for considering a hit unique
my $debug = undef;
# Parameters for BLAT. Thanks to Judy Savitskaya for these additions
my $tileSize = 11; # size of an alignment tile
my $stepSize = 11; # distance between the starting bases of alignment tiles (will overlap if stepSize<tileSize)

# Given BLAT rows (as list of lists), output the mapping for the read
sub HandleGenomeBLAT($$);

# Global so that HandleGenomeBLAT() can use them
my $nMapped = 0;
my $nMapUnique = 0;
my $nPastEndIgnored = 0; # weak hit to past-end ignored
my $nPastEndTrumps = 0; # hit to past-end (close to) as good as hit to genome

my $usage = <<END
Usage: MapUnbarcoded.pl [ -debug ] [ -limit maxReads ]
            [ -minIdentity $minIdentity ] [ -minScore $minScore ] [ -delta $delta ]
            [ -wobble $wobbleAllowed ] [ -extra 0 ]
            [ -tileSize $tileSize ] [ -stepSize $stepSize ]
            [-tmpdir $tmpdir ]
            -genome fasta_file -model model_file -first fastq_file > output_file

    The fastq file should have phred+33 ("sanger") encoding of quality scores
    (as in MiSeq or 2012+ HiSeq). If it is named .gz, it will be gunzipped before reading.
    If it contains paired-end reads, the second read (name matching " 2:") will be ignored.

    The model file contains 2 lines -- the first shows what a typical
    read should look like up until the junction with the genome, e.g.
    CGCCCTGCAGGGATGTCCACGAGGTCTCTAGGTCGACGGCCGGCCAGACCGGGGACTTATCAGCCAACCTGT

    The second line in the model file contains the sequence "past the
    end" of the transposon that might arise from residual intact
    plasmid. All characters in the model read must be ACGT.

    This script does not handle sample demultiplexing.

    The output file is tab-delimited and contains, for each usable
    read, the read name, which scaffold the insertion
    lies in, the position of the insertion, the strand that the read
    matched, a boolean flag for if this mapping location is unique or
    not, the beginning and end of the hit to the genome in the read
    after trimming the transposon sequence, the bit score, and the
    %identity.

    extra is the number of extra characters before the model (0 default).

    tileSize and stepSize are parameters for BLAT.  If the portion of
    the read after the junction is short, try a lower stepSize.  Also
    see BLAT's documentation.
END
    ;

sub FindModelEnd($$$);
sub FindSubstr($$$$);
sub BLAT8($$$$$$$$); # BLAT to a blast8 format file

{
    my $limit = undef;
    my $fastqFile = undef;
    my $genomeFile = undef;
    my $modelFile = undef;
    my $blatcmd = -e "$Bin/blat" ? "$Bin/blat" : "blat";
    my $minGenomeId = 90;
    my $nExtra = 0;

    (GetOptions('debug' => \$debug,
                'limit=i' => \$limit,
                'minIdentity=i' => \$minIdentity,
                'minScore=i' => \$minScore,
                'delta=i' => \$delta,
                'tileSize=i' => \$tileSize,
                'stepSize=i' => \$stepSize,
                'tmpdir=s' => \$tmpdir,
                'blat=s' => \$blatcmd,
                'wobble=i' => \$wobbleAllowed,
                'extra=i' => \$nExtra,
                'genome=s' => \$genomeFile,
                'model=s' => \$modelFile,
                'first=s' => \$fastqFile)
     && @ARGV==0) || die $usage;
    die $usage unless defined $modelFile && defined $fastqFile;

    die "Cannot read $genomeFile" unless -r $genomeFile;
    die "Cannot read $fastqFile" unless -r $fastqFile;
    die "Not a directory: $tmpdir" unless -d $tmpdir;
    die "minScore must be at least 10" if $minScore < 10;
    die "tileSize must be at least 7" if $tileSize < 7;
    die "stepSize must be at least 1" if $stepSize < 1;
    die "Invalid minimum identity" if $minIdentity < 50 || $minIdentity > 100;
    die "delta cannot be negative\n" if $delta < 0;

    open(MODEL, "<", $modelFile) || die "Cannot read $modelFile";
    my $model = <MODEL>;
    $model =~ s/[\r\n]+$//;
    die "Invalid model: $model" unless $model =~ m/^[ACGT]+$/;
    my $pastEnd = <MODEL>;
    die "No pastEnd\n" unless defined $pastEnd;
    chomp $pastEnd;
    die "Invalid past-end sequence: $pastEnd" unless $pastEnd =~ m/^[ACGT]+$/;
    close(MODEL) || die "Error reading $modelFile";

    print STDERR "Parsed model $modelFile\n";
    my $nCheck = 8;
    my $checkAt = $nExtra + length($model) - $nCheck;
    my $checkString = substr($model, length($model) - $nCheck);

    print STDERR "Expecting pre-junction of $checkAt +/- $wobbleAllowed\n";

    my $pipe = 0;
    if ($fastqFile =~ m/[.]gz$/) {
        $pipe = 1;
        open(FASTQ, '-|', 'zcat', $fastqFile) || die "Cannot run zcat on $fastqFile";
        print STDERR "Reads (gzipped) from $fastqFile\n";
    } elsif ($fastqFile =~ m/[.zip]$/) {
        $pipe = 1;
        open(FASTQ, "7za -so e $fastqFile | zcat |") || die "Cannot run 7za and zcat on $fastqFile";
        print STDERR "Reads (zipped) from $fastqFile\n";
    } else {
        open(FASTQ, "<", $fastqFile) || die "Cannot read $fastqFile";
        print STDERR "Reads (raw) from $fastqFile\n";
    }

    my $rand = rand();
    my $tmpPre = $tmpdir . "/MapUnbarcoded_" . $$ . "_$rand";
    my $tmpFna = "$tmpPre.fna";
    open(TMPFNA, ">", $tmpFna) || die "Cannot write to $tmpFna";

    my $nReads = 0;
    my $nTryToMap = 0;

    # Find barcodes and end of transposon and write remaining sequence to TMPFNA
    while(!defined $limit || $nReads < $limit) {
      my $name = <FASTQ>;
      (last, next) unless $name;
      chomp $name;
      die "Sequence name line does not start with @" unless $name =~ m/^@/;

      my $seq = <FASTQ>;
      chomp $seq;
      die "Sequence line is empty"
        unless defined $seq && length($seq) > 0;
      die "Sequence line contains invalid characters: $seq"
        unless $seq =~ m/^[A-Z]+$/;
      my $line = <FASTQ>;
      die "Third line does not start with +"
        unless defined $line && $line =~ m/^[+]/;
      my $quality = <FASTQ>;
      chomp $quality;
      die "Quality line is of wrong length"
        unless defined $quality && length($quality) == length($seq);

      next if $name =~ m/^\S+ 2:/; # ignore second side of paired-end reads
      $nReads++;

      my $postJunction;
      for (my $i = $checkAt - $wobbleAllowed; $i <= $checkAt + $wobbleAllowed; $i++) {
        if (length($seq) >= $i + $nCheck && substr($seq, $i, $nCheck) eq $checkString) {
          $postJunction = substr($seq, $i + $nCheck);
          last;
        }
      }
      next unless defined $postJunction;

      next unless length($postJunction) >= $minScore;

      my $shortName = $name; $shortName =~ s/ .*$//;
      print STDERR "Try to map $shortName with $postJunction\n" if $debug;
      print TMPFNA ">$shortName\n$postJunction\n";
      $nTryToMap++;
    }

    # Prematurely closing the pipe gives an error
    close(FASTQ) || ($pipe && defined $limit) || die "Error reading from $fastqFile: $!";
    close(TMPFNA) || die "Error writing to $tmpFna";
    print STDERR "Read $nReads reads\n";

    if ($nTryToMap == 0) {
      print STDERR "None of the reads are candidates for mapping (none match the model and are long enough)\n";
      unlink($tmpFna);
      exit(0);
    }

    my %hitsPastEnd = (); # short name to score

    # Map to past-end-of transposon
    my $endFna = "$tmpPre.endfna";
    open(END, ">", $endFna) || die "Cannot write to $endFna";
    print END ">pastend\n$pastEnd\n";
    close(END) || die "Error writing to $endFna";
    my $blat8 = BLAT8($tmpFna, $endFna, $tmpdir, $blatcmd, $minScore, $minIdentity,$tileSize,$stepSize);
    print STDERR "Parsing past-end hits to $blat8\n" if defined $debug;
    open(BLAT, "<", $blat8) || die "Cannot read $blat8";
    while(<BLAT>) {
      chomp;
      my @F = split /\t/, $_;
      my ($query, $subject, $identity, $len, $mm, $gaps, $qBeg, $qEnd, $sBeg, $sEnd, $eval, $score) = @F;

      $hitsPastEnd{$query} = $score unless exists $hitsPastEnd{$query} && $hitsPastEnd{$query} > $score;
    }
    close(BLAT) || die "Error reading $blat8";
    unlink($blat8);
    unlink($endFna);

    # Map to the genome
    $blat8 = BLAT8($tmpFna, $genomeFile, $tmpdir, $blatcmd, $minScore, $minIdentity, $tileSize, $stepSize);
    print STDERR "Parsing $blat8\n" if defined $debug;
    open(BLAT, "<", $blat8) || die "Cannot read $blat8";
    my @lines = ();
    while(<BLAT>) {
      chomp;
      my @F = split /\t/, $_;
      my ($query, $subject, $identity, $len, $mm, $gaps, $qBeg, $qEnd, $sBeg, $sEnd, $eval, $score) = @F;
      if (exists $hitsPastEnd{$query}) {
        if ($hitsPastEnd{$query} < $score) {
          $hitsPastEnd{$query} = 0; # mark as ignored
          $nPastEndIgnored++;
        } else {
          next; # skip this hit
        }
      }
      if (@lines == 0 || $query eq $lines[0][0]) {
        push @lines, \@F;
      } else {
        HandleGenomeBLAT(\@lines, \%hitsPastEnd);
        @lines = \@F;
      }
    }
    HandleGenomeBLAT(\@lines, \%hitsPastEnd);

    close(BLAT) || die "Error reading $blat8";
    unlink($blat8);

    unlink($tmpFna);

    # print out hits-past-end
    while (my ($read,$score) = each %hitsPastEnd) {
        next unless $score > 0; # score=0 means it was ignored above
        print join("\t", $read, "pastEnd", "", "", "", "", "", "", "")."\n";
    }

    print STDERR "Reads processed $nReads\n";
    print STDERR "Mapping attempted for $nTryToMap Mapped $nMapped Uniquely $nMapUnique\n";
    print STDERR "Hits past end of transposon: " . (scalar(keys %hitsPastEnd) - $nPastEndIgnored) . "\n";
    print STDERR sprintf("Proportions: Attempted %.3f Mapped %.3f Past-end %.3f\n",
			 $nTryToMap/$nReads,
			 $nMapped/$nReads,
			 (scalar(keys %hitsPastEnd) - $nPastEndIgnored)/$nReads)
	if $nReads > 0;
}

# returns the name of the file
sub BLAT8($$$$$$$$) {
    my ($queriesFile,$dbFile,$tmpdir,$blatcmd,$minScore,$minIdentity,$tileSize,$stepSize) = @_;

    my $blat8File = "$tmpdir/MapUnbarcoded.$$.".rand() . ".psl";
    # prevent writing to stdout by BLAT
    open(OLD_STDOUT, ">&STDOUT");
    print OLD_STDOUT ""; # prevent warning
    open(STDOUT, ">", "/dev/null");
    my @cmd = ($blatcmd, "-out=blast8", "-t=dna", "-q=dna", "-minScore=$minScore", "-minIdentity=$minIdentity", "-maxIntron=0", "-noTrimA", "-tileSize=$tileSize", "-stepSize=$stepSize", $dbFile, $queriesFile, $blat8File);
    print STDERR join(" ", "Running:", @cmd)."\n" if $debug;
    system(@cmd) == 0 || die "Cannot run $blatcmd: $!";
    close(STDOUT);
    open(STDOUT, ">&OLD_STDOUT") || die "Cannot restore stdout: $!";
    return($blat8File);
}

sub HandleGenomeBLAT($$) {
    my ($rows, $hitsPastEnd) = @_;
    return if scalar(@$rows) == 0;
    my $read = $rows->[0][0];
    die if !defined $read;

    my @hits = ();

    # indexes within besthits entries
    my ($SCAFFOLD,$POSITION,$STRAND,$SCORE,$QBEG,$QEND) = (0,1,2,3,4,5);
    my @besthits = ();

    foreach my $row (@$rows) {
        my ($read2, $subject, $identity, $len, $mm, $gaps, $qBeg, $qEnd, $sBeg, $sEnd, $eval, $score) = @$row;
        die unless $read2 eq $read;
        if (scalar(@besthits) == 0 || $score >= $besthits[0][$SCORE] - $delta) {
            # convert from 0-based to 1-based position, and note that sBeg always < sEnd so flip if stranded
            push @besthits,  [ $subject, $sBeg, ($sBeg < $sEnd ? "+" : "-"), $score,
                               $identity, $qBeg, $qEnd ];
            print STDERR "Hit for $read:$qBeg:$qEnd to $subject:$sBeg:$sEnd score $score identity $identity\n"
                if $debug;
        }
    }

    die if @besthits == 0;
    my ($scaffold, $position, $strand, $score, $identity, $qBeg, $qEnd) = @{ $besthits[0] };

    # and output a mapping row (or none)
    if (exists $hitsPastEnd->{$read}) {
        if ($hitsPastEnd->{$read} >= $score - $delta) {
            $nPastEndTrumps++;
            print STDERR "Past end trumps for $read\n" if $debug;
            return;
        } else {
            $nPastEndIgnored++; # ignore weak hit to past-end sequence
            print STDERR "Ignoring hit past end of transposon for $read\n" if $debug;
            $hitsPastEnd->{$read} = 0; # so we do not print out a line for it later on
        }
    }
    # else if no trumping
    print join("\t", $read, $scaffold, $position, $strand, @besthits == 1 ? 1 : 0,
               $qBeg, $qEnd, $score, $identity)."\n";
    $nMapUnique++ if @besthits == 1;
    $nMapped++;
    return;
}
