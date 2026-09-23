#!/usr/bin/perl -w
# Process a small set of barseq test samples, including MultiCodes.pl, combineBarSeq.pl, and BarSeqR.pl
use strict;
use warnings;
use Getopt::Long;
use FindBin qw($Bin);
use lib "$Bin/../lib";
use FEBA_Utils qw{CodesByIndex FastqByIndex};

my $usage =<<END
Usage: BarSeqTest.pl -org organism [ -n25 | -bs3 | -bs4 ] -index S1:S2:S3 -desc Time0:Lactate:Glucose
	    -fastqdir directory

Assumes that g/organism includes all the typical files, including a pool.n10
or pool file, and writes to g/organism/barseqtest by default.

The fastq directory should include file(s) for each index. (See
RunBarSeqLocal.pl for details)

Optional arguments:
  -pool g/organism/pool.n10
  -outdir g/organism/barseqtest
  -hascodes -- use the *.codes files that already exist in the fastq directory
  -test -- do not actually run any commands
END
    ;

sub maybeRun($); # run command unless $test is defined
my $test = undef; # check for files, but do no work

{
    my ($org, $indexSpec, $descSpec, $fastqdir, $pool, $outdir);
    my ($n25, $bs3, $bs4, $hascodes);

    GetOptions('org=s' => \$org,
	       'index=s' => \$indexSpec,
	       'desc=s' => \$descSpec,
	       'fastqdir=s' => \$fastqdir,
	       'pool=s' => \$pool,
	       'test' => \$test,
               'n25' => \$n25,
               'bs3' => \$bs3,
               'bs4' => \$bs4,
               'hascodes' => \$hascodes,
	       'outdir=s' => \$outdir) || die $usage;
    @ARGV == 0 || die $usage;
    die $usage unless defined $org && defined $indexSpec && defined $descSpec && defined $fastqdir;
    die "No such directory: g/$org" unless -d "g/$org";
    die "No such directory: $fastqdir" unless -d $fastqdir;
    die "Specify just one of -n25 -bs4 -bs4" if (defined $n25) + (defined $bs3) + (defined $bs4) > 1;

    if (!defined $pool) {
	$pool = "g/$org/pool.n10";
	unless (-e $pool) {
	    $pool = "g/$org/pool";
	    die "No such file: g/$org/pool.n10 or g/$org/pool" unless -e $pool;
	}
    } else {
	die "No such file: $pool" unless -e $pool;
    }

    # parse the index and description fields
    my @indexes = split /:/, $indexSpec, -1;
    die "Invalid -index $indexSpec -- must have at least two, colon-separated " if @indexes < 2;
    my @desc = split /:/, $descSpec, -1;
    die "Invalid -desc $descSpec -- must match -index" if scalar(@desc) != scalar(@indexes);
    my @time0s = grep { $_ eq "Time0" } @desc;
    die "Must have both Time0 and non Time0 samples" if scalar(@time0s) < 1 || scalar(@time0s) == scalar(@desc);

    my %index;
    my %numToIndex;
    foreach my $index (@indexes) {
      $index{$index} = 1;
      if ($index =~ m/^[A-Z]+0*(\d+)$/) {
        $numToIndex{$1} = $index;
      }
    }
    die "Non-unique indexes\n" unless scalar(keys %index) == scalar(@indexes);
    if (scalar(@indexes) != scalar(keys %numToIndex)) {
      print STDERR "Warning: cannot identify unique sample numbers for each index; not matching by number\n";
      %numToIndex = ();
    }

    my $bynumber = 0;
    if (defined $hascodes) {
      # bynumber if index is of style S1
      $bynumber = 1 if $indexes[0] =~ m/^S\d+$/;
    } else {
      $bynumber = 1 if defined $bs4;
    }

    my $indexToCodes = {};
    my $indexToFq = {};
    if (defined $hascodes) {
      $indexToCodes = CodesByIndex('in' => $fastqdir,
                                   'bynumber' => $bynumber,
                                   'bywell' => 0,
                                   'offset' => undef);
    } else {
      $indexToFq = FastqByIndex('in' => $fastqdir,
                                'bynumber' => $bynumber,
                                'bywell' => 0,
                                'offset' => undef);
    }
    $outdir = "g/$org/barseqtest" if !defined $outdir;
    mkdir($outdir) if ! -d $outdir && !defined $test;
    my @codesFiles = ();
    foreach my $index (@indexes) {
      if (defined $hascodes) {
        die "No codes files for $index in $fastqdir\n" unless exists $indexToCodes->{$index};
        push @codesFiles, @{ $indexToCodes->{$index} };
      } else {
        if (!exists $indexToFq->{$index}) {
          print STDERR "Warning! No fastq files for $index in $fastqdir\n";
          next;
        }
        my @filenames = @{ $indexToFq->{$index} };
        my $pathSpec = join(" ", @filenames);
        push @codesFiles, "$outdir/$index.codes";
        my $extraopt = "";
        $extraopt = "-n25" if defined $n25;
        $extraopt = "-bs3" if defined $bs3;
        $extraopt = "-bs4" if defined $bs4;
        &maybeRun("zcat $pathSpec | $Bin/MultiCodes.pl $extraopt -minQuality 0 -index $index -out $outdir/$index");
      }
    }
    &maybeRun("$Bin/combineBarSeq.pl -all $outdir/test $pool " . join(" ", @codesFiles));
    unless (defined $test) {
	# Write the metadata file
	open(EXPS, ">", "$outdir/exps_table") || die "Cannot write to $outdir/exps_table";
	print EXPS join("\t", qw{SetName Index Person Date_pool_expt_started Description})."\n";
	foreach my $i (0..(scalar(@indexes)-1)) {
	    print EXPS join("\t", "test", $indexes[$i], "someone", "somedate", $desc[$i])."\n";
	}
	close(EXPS) || die "Error writing to $outdir/exps_table";
    }
    # Remove any strain usage files from the output directory, otherwise they may confuse
    # BarSeqR.pl
    if (!defined $test) {
      unlink("$outdir/strainusage.barcodes");
      unlink("$outdir/strainusage.genes");
      unlink("$outdir/strainusage.genes12");
    }
    &maybeRun("$Bin/BarSeqR.pl -org $org -exps $outdir/exps_table  -pool $pool -indir $outdir -outdir $outdir -genes g/$org/genes.GC test");
}

sub maybeRun($) {
    my ($cmd) = @_;
    if (defined $test) {
        print STDERR "Would run: $cmd\n";
    } else {
	print STDERR "Running: $cmd\n";
        system($cmd) == 0 || die "script failed: $cmd";
    }
}
