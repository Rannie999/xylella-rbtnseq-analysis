#!/usr/bin/perl -w
# RunBarSeqLocal.pl -- process a lane of barseq into a table of counts for each strain in each sample
# issuing jobs locally using submitter.pl
use strict;
use warnings;
use Getopt::Long;
use FindBin '$Bin';
use File::Basename;
use lib "$Bin/../lib";
use FEBA_Utils qw{FastqByIndex};

my $minQuality = 0;

my $usage = <<END
Usage:
RunBarSeqLocal.pl [ -n25 | -bs3 | -bs4 | -hascodes ]
    -in directory_with_fastq.gz_files
    [ -sets organism1_libraryname1_setName1,...,organismN_libraryN_setNameN | -nopool ]

Examples:
feba/bin/RunBarSeqLocal.pl -n25 -in HiSeq_barcodes/FEBA_BS_186 -sets MR1_ML3_set12,Koxy_ML2_set9
feba/bin/RunBarSeqLocal.pl -bs3 -in HiSeq_barcodes/FEBA_BS_265 -nopool

RunBarSeq.pl has two phases -- the first phase counts the barcodes in
the *.fastq.gz or *.fq.gz files in the specified directory, using
MultiCodes.pl. Usually these files are in the directory, with the
primer name within the file name. If there are no fastq.gz files
(other than "Undetermined" files) in the directory, it will look
across all subdirectories. The -n25, -bs3, or -bs4 options indicate
which primer design was used and are passed to MultiCodes.pl. Use
-hascodes to skip the 1st phase.

The output files from the first phase are the *.codes files in that
same directory.

The second phase aggregates the counts of these
barcodes. RunBarSeqLocal.pl ignores barcodes that do not match the
pool definition, which should be in g/organism/pool or
g/organism/pool.n10

The second phase is skipped if -nopool is used.

Output files from the 2nd phase are in g/organism/ --
setName.colsum -- total parsed reads per index
setName.poolcount -- counts for strains in the pool
    based on joining to the pool or pool.n10 file.
setName.codes.ignored -- all the counts for the barcodes
	that did not match the pool.

Options for identifying the index from the file name:
  -offset -- If IT001 is named S97, use -offset 96
  -bywell -- fastq.gz file names begin with A01_ for IT001, etc.

Other options:
  -minQuality -- Set the minimum quality for each nucleotide in the barcode
	(default is $minQuality)
  -test -- do not do any work, just show what commands would be run
  -limit -- limit the #reads analyzed per sample (mostly for debugging)
  -preseq, -postseq, -nPreExpected -- see MultiCodes.pl
  -cmdsFile -- store the commands in this file
  -force -- overwrite existing poolcount files
END
  ;

sub maybeRun($);                # run command unless $test is defined

my $test = undef;
{
  my ($setspec, $inDir, $nopool);
  my ($n25, $bs3, $bs4, $hasCodes);
  my ($preseq, $postseq, $nPreExpected, $limitReads);
  my ($bywell, $offset);
  my $cmdsFile;
  my $force;

  die $usage
    unless GetOptions('test' => \$test,
                      'n25' => \$n25,
                      'bs3' => \$bs3,
                      'bs4' => \$bs4,
                      'hascodes' => \$hasCodes,
                      'minQuality=i' => \$minQuality,
                      'limit=i' => \$limitReads,
                      'preseq=s' => \$preseq,
                      'postseq=s' => \$postseq,
                      'nPreExpected=s' => \$nPreExpected,
                      'in=s' => \$inDir,
                      'sets=s' => \$setspec,
                      'offset=i' => \$offset,
                      'nopool' => \$nopool,
                      'bywell' => \$bywell,
                      'cmdsFile=s' => \$cmdsFile,
                      'force' => \$force)
    && defined $inDir
    && (defined $nopool | defined $setspec);
  die "No such directory: $inDir\n" unless -d $inDir;

  # Parse the sets specifier
  my @setNames = ();
  my %setToGdir = ();
  my %gdirs = ();               # unique gdir => first set
  # If there are multiple sets for the same organism in the same sequencing run,
  # run the sripts for just the first one, and the make a symbolic link to create the
  # other polcount files.
  my %dupSets = ();             # duplicateSet name => baseSet
  die "Cannot use both -sets and -nopool\n" if defined $setspec && defined $nopool;
  if (defined $setspec) {
    die $usage unless defined $inDir && @ARGV == 0;
    @setNames = split /[,;]\s*/, $setspec;
    my @misc = grep m/^misc_set|htcomp|bartn7/i, @setNames;
    print STDERR "Not running combineBarSeq for miscellaneous sets: @misc\n" if @misc > 0;
    @setNames = grep !m/^misc_set|htcomp|bartn7/i, @setNames;
    foreach my $setName (@setNames) {
      my @setparts = split /_/, $setName;
      die "Invalid set name: $setName" unless @setparts > 1;
      my $gdir = undef;
      if ($setName =~ m/^(.*)_ML/i) {
        my $org = $1;
        $gdir = "g/$org" if -d "g/$org";
      }
      my $beg = join("_", $setparts[0], $setparts[1]);
      if (!defined $gdir) {
        $gdir = "g/$beg" if -d "g/$beg";
      }
      if (!defined $gdir) {
        # Try to use glob to find the mapping from the first two parts of the library name to an organism
        my @hits = glob("g/*/$beg*");
        die "Could not convert $setName to an organism, and gdir does not exist"
          unless @hits > 0;
        my $hit1 = $hits[0];
        my @parts = split "/", $hit1;
        $gdir = "g/$parts[1]";
        print STDERR "Set $setName matches directory $gdir\n";
      }
      die unless -d $gdir;
      $setToGdir{$setName} = $gdir;
      if (exists $gdirs{$gdir}) {
        die "Duplicate entries for $setName which is not allowed\n"
          if $setName eq $gdirs{$gdir};
        print STDERR "Multiple sets for $gdir: will create symbolic link from $setName to $gdirs{$gdir}\n";
        $dupSets{$setName} = $gdirs{$gdir};
      } else {
        $gdirs{$gdir} = $setName;
      }
    }
    # and check that all of these are new
    foreach my $setName (@setNames) {
      my $gdir = $setToGdir{$setName};
      my $out = "$gdir/$setName";
      foreach my $file ("$out.poolcount") {
        if (-s $file) {         # if not empty
          if (defined $force) {
            print STDERR "Warning: overwriting $file\n";
          } else {
            die join("\n",
                     "Error: output for $gdir $setName already exists!",
                     "See file $file",
                     "Please remove the colsum and poolcount files or correct the set name",
                     "");
          }
        }
      }
    }
  }

  my %gdirToPool = ();
  foreach my $gdir (sort keys %gdirs) {
    my $poolfile = "$gdir/pool";
    $poolfile = "$poolfile.n10" if !-e $poolfile;
    die "Cannot find $poolfile (or withut .n10)" unless -e $poolfile;
    $gdirToPool{$gdir} = $poolfile;
  }

  my $indexToFq;
  unless (defined $hasCodes) {
    $indexToFq = FastqByIndex('in' => $inDir,
                              'bynumber' => defined $bs4,
                              'bywell' => defined $bywell,
                              'offset' => $offset);
    die "No input reads found in $inDir\n"
      if %$indexToFq == 0;
  }
  my @indexes = sort keys %$indexToFq;

  # Check that we can read all the input files
  foreach my $index (@indexes) {
    foreach my $inFile (@{ $indexToFq->{$index} }) {
      open(my $fh, "<", $inFile) || die "Cannot read from $inFile\n";
      close($fh) || die "Error reading from $inFile";
    }
  }

  my @outPre = ();        # the output files without the .codes suffix
  # Are the fastqs in the directory (inBase = true), or an extra level down?
  my ($fqDir, $inBase);
  if (defined $hasCodes) {
    my @codes = glob("$inDir/*.codes");
    if (@codes > 0) {
      $inBase = 1;
    } else {
      my @codes = glob("$inDir/*/*.codes");
      die "No codes files in $inDir/*.codes or $inDir/*/*.codes\n"
        unless @codes > 0;
      $inBase = 0;
    }
  } else {
    $fqDir = dirname($indexToFq->{ $indexes[0] }[0]);
    my $inDir2 = $inDir; $inDir2 =~ s!/*$!!;
    $inBase = $fqDir eq $inDir2;
  }

  unless (defined $hasCodes) {
    # build the list of commands, and the list of codes files we'll make
    # (Also remove any codes files that would be overwritten)
    $cmdsFile = "$inDir/BarSeq.codecmds"
      unless defined $cmdsFile;
    maybeRun("rm -f $cmdsFile*") if -e $cmdsFile;
    open(CMDS, ">", $cmdsFile) || die "Cannot write to $cmdsFile";
    foreach my $index (@indexes) {
      foreach my $inFile (@{ $indexToFq->{$index} }) {
        print STDERR "Index $index fq $inFile\n" if $test;
        my $outPre = $inFile;
        $outPre =~ s/[.]f[astq]+[.]gz$//;
        push @outPre, $outPre;
        my $corecmd = "zcat $inFile | $Bin/MultiCodes.pl -minQuality $minQuality -index $index -out $outPre";
        $corecmd .= " -n25" if defined $n25;
        $corecmd .= " -bs3" if defined $bs3;
        $corecmd .= " -bs4" if defined $bs4;
        $corecmd .= " -limit $limitReads" if defined $limitReads;
        $corecmd .= " -preseq $preseq" if defined $preseq;
        $corecmd .= " -postseq $postseq" if defined $postseq;
        $corecmd .= " -nPreExpected $nPreExpected" if defined $nPreExpected;
        print CMDS "$corecmd >& $outPre.log", "\n";
      }
    }
    close(CMDS) || die "Error writing to $cmdsFile";

    if ($inBase) {
      maybeRun("rm -f $inDir/*.codes");
      maybeRun("rm -f $inDir/*.codes.gz");
    } else {
      maybeRun("rm -f $inDir/*/*.codes");
      maybeRun("rm -f $inDir/*/*.codes.gz");
    }

    # run them all
    maybeRun("$Bin/submitter.pl $cmdsFile");

    if (!defined $test) {
      # Parse all the log files
      my @logs = map "$_.log", @outPre;
      my $nReads = 0;
      my $nMulti = 0;
      my $nUsable = 0;
      foreach my $log (@logs) {
        if (-e $log) {
          open(LOG, "<", $log) || die "Cannot read $log";
          my $nLines = 0; # numbers of lines parsed with counts of reads
          while (<LOG>) {
            if (m/^Reads\s+(\d+)\s+Multiplexed\s+(\d+)\s+Usable\S*\s+(\d+)\b/) {
              $nReads += $1;
              $nMulti += $2;
              $nUsable += $3;
              $nLines++;
            }
          }
          close(LOG) || die "Error reading $log";
          print STDERR "Warning: no Reads/Multiplexed/Usable line in $log\n"
            if $nLines == 0;
        } else {
          print STDERR "Warning: no log file $log\n";
        }
      }
      print STDERR sprintf("Total reads (in millions) %.1f Multi %.1f (%.1f%%) Usable %.1f (%.1f%%) from %d files\n",
                           $nReads/1e6,
                           $nMulti/1e6, 100.0*$nMulti/($nReads+0.1),
                           $nUsable/1e6, 100.0*$nUsable/($nReads+0.1),
                           scalar(@outPre));
    }
  }                             # end if !defined $hasCodes

  if (defined $nopool) {
    print STDERR "codes files were computed; skipping comparison to pool(s)\n"
      unless defined $test || defined $hasCodes;
  } else {
    foreach my $setName (@setNames) {
      next if exists $dupSets{$setName};
      my $gdir = $setToGdir{$setName};
      my $path = $gdir . "/" . $setName;
      print STDERR "Combining codes files into $path\n" unless defined $test;
      my @codesFiles = map "$_.codes", @outPre;
      my $cmd = "$Bin/combineBarSeq.pl $path $gdirToPool{$gdir}";
      if ($inBase) {
        $cmd .= " $inDir/*.codes";
      } else {
        $cmd .= " $inDir/*/*.codes";
      }
      maybeRun($cmd);
    }
    foreach my $setName (sort keys %dupSets) {
      my $gdir = $setToGdir{$setName};
      my $baseSet = $dupSets{$setName};
      unless (defined $test) {
        print STDERR "Linking $setName to $baseSet\n";
        unlink("$gdir/$setName.poolcount");
      }
      maybeRun("ln -s $baseSet.poolcount $gdir/$setName.poolcount");
    }
    print STDERR "Finished " . scalar(@setNames) . " sets\n";
  }
}

sub maybeRun($) {
  my ($cmd) = @_;
  if (defined $test) {
    print STDERR "Would run: $cmd\n";
  } else {
    system($cmd) == 0 || die "script failed: $cmd";
  }
}
