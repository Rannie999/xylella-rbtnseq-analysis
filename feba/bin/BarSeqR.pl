#!/usr/bin/perl -w
# BarSeqR.pl -- combine the individual sets from combineBarSeq.pl along with the genes table
# and use FEBA.R to produce the results
#
# Limitations:
# Genes that wrap around the origin (i.e., begin > end) are ignored (no strain will map within them)
use strict;
use Getopt::Long;
use FileHandle;
use FindBin qw($Bin);
use lib "$Bin/../lib";
use FindGene;                   # for LocationToGene()
use FEBA_Utils;                 # for ReadTable(), ReadColumnNames()
use Compounds;

my $metadir = "$Bin/../metadata";
my $usage = <<END
Usage: BarSeqR.pl -org organism [ -indir g/organism ]
    [ -exps indir/FEBA_BarSeq.tsv ] [ -genes indir/genes.GC ]
    [ -pool indir/pool.n10 or indir/pool ]
    [  -outdir html/organism ]
    [ -metadir $metadir ]
    [ -test | -noR ]
    [ setnames ]

    By default, the input directory includes FEBA_BarSeq.tsv,
    genes.GC, and setname.poolcount, and the "all.poolcount" file is
    also written to this directory.  This script requires the genes
    table to include the fields locusId, sysName, scaffoldId, begin,
    end, and strand, and the experiments table to include the fields
    SetName, Index, Description, and Date_pool_expt_started:

      SetName -- which poolcount file the data for this sample is in
        for example, "Keio_ML9_set1" means look in
        g/Keio/Keio_ML9_set1.poolcount
      Index -- what the column name is in the poolcount file
      Group -- set to "Time0" for control samples
      Date_pool_expt_started -- each experiment is compared to the Time0
        sample(s) with a matching value. If possible, only Time0
        samples within the same set are used.

    By default, all sets in the experiments table are processed except
    for test sets, which are ignored.  If no set.poolcount is
    available, a warning is issued. Optionally, the specific set(s) to
    process can be specified.

    Sometimes a set is sequenced on multiple lanes. The extra files
    should have suffixes. E.g. for set1, the files would be "set1b" or
    "set1_repN" or "set1_seqNN". A warning is issued whenever this
    occurs. Otherwise, the SetName in the experiments table should
    exactly match to indir/SetName.poolcount and indir/SetName.colsum.

    Creates files in the output directory that can be used by
    RunFEBA.R to compute fitness values and create a mini web site in
    the output directory. These files are genes, pool, exps, and
    all.poolcount, which contains the barseq data for each strain in
    the pool, along with whether the strain is in a gene.

    Note -- the R step is parallel, use the MC_CORES environment
    variable to control it.

    Note -- by default, strainusage files in the input directory are
    used if present. Use FEBA_NO_STRAIN_USAGE=1 to turn this off.

END
  ;

{
  my ($org, $indir, $expsfile, $genesfile, $poolfile, $outdir);
  my $noR;                      # set if skip running FEBA.R
  my $test;                     # check for poolcount files, etc.

  die $usage unless
    GetOptions('org=s' => \$org,
               'indir=s' => \$indir,
               'metadir=s' => \$metadir,
               'exps=s' => \$expsfile,
               'genesfile=s' => \$genesfile,
               'poolfile=s' => \$poolfile,
               'outdir=s' => \$outdir,
               'noR' => \$noR,
               'test' => \$test )
      && defined $org;
  my @sets = @ARGV;

  $indir = "g/$org" unless defined $indir;
  die "No such directory: $indir" unless -d $indir;

  die "No such directory : $metadir" unless -d $metadir;

  $expsfile = "$indir/FEBA_BarSeq.tsv" unless defined $expsfile;
  die "No such experiments file: $expsfile" unless -e $expsfile;

  $genesfile = "$indir/genes.GC" unless defined $genesfile;
  die "No such genes file: $genesfile" unless -e $genesfile;

  if (!defined $poolfile) {
    $poolfile = "$indir/pool.n10";
    $poolfile = "$indir/pool" if ! -e $poolfile;
    die "No pool file: $poolfile.n10 or $poolfile" if ! -e $poolfile;
  } else {
    die "No such pool file: $poolfile" if ! -e $poolfile;
  }
  my @poolcols = ReadColumnNames($poolfile);
  my %poolcols = map { $_ => 1 } @poolcols;
  foreach my $col (qw{barcode rcbarcode scaffold strand pos n nTot}) {
    print STDERR "Warning: no column named $col in $poolfile\n" unless exists $poolcols{$col};
  }

  if (!defined $outdir) {
    die "No such directory: html" if ! -d "html";
    $outdir = "html/$org";
    mkdir($outdir) if ! -d $outdir;
  } else {
    die "No such directory: $outdir" if ! -d $outdir;
  }

  my $Rscript = "$Bin/RunFEBA.R";
  die "No such script file: $Rscript" if ! -e $Rscript && !defined $noR;

  LoadCompounds($metadir);
  LoadMedia($metadir);

  my @exps = &ReadTable($expsfile, qw{SetName Index Description Date_pool_expt_started});
  my $alpha = chr(206).chr(177);
  my $Delta = chr(206).chr(148);
  # Wierd dash that sometimes gets introduced into metadata
  my $dash = chr(0xe2).chr(0x88).chr(0x92);
  foreach my $exp (@exps) {
    # extra spaces at end are common data entry errors
    $exp->{Description} =~ s/ *$//;
    $exp->{SetName} =~ s/ *$//;
    $exp->{Index} =~ s/ *$//;
    # Fix non-ascii characters
    foreach my $field (qw{Description Condition_1 Condition_2 Condition_3 Condition_4}) {
      next if !exists $exp->{$field};
      $exp->{$field} =~ s/^ +//;
      $exp->{$field} =~ s/ +$//;
      $exp->{$field} =~ s/$alpha/a/g;
      if ($field eq "Description") {
        $exp->{$field} =~ s/$Delta/&Delta;/g;
      } else {
        $exp->{$field} =~ s/$Delta/Delta/g;
      }
      $exp->{$field} =~ s/$dash/-/g;
      # replace other non-ascii characers with .
      $exp->{$field} =~ s/[^[:ascii:]]+/./g;
    }
  }

  @exps = grep { $_->{Description} ne "" } @exps;
  my $prespec_sets = scalar(@sets) > 0;
  if ($prespec_sets) {
    # ignore experiments not in pre-specified sets
    my %sets = map { $_ => 1 } @sets;
    @exps = grep { exists $sets{ $_->{SetName} } } @exps;
    die "No experiments in specified sets (having Description filled out)\n" if scalar(@exps) == 0;
  } else {
    # ignore tests
    my %sets = map { $_->{SetName} => 1 } @exps;
    foreach my $set (keys %sets) {
      if ($set =~ m/test/) {
        print STDERR "Ignoring SetName = $set\n";
      }
    }
    @exps = grep { $_->{SetName} !~ /test/ } @exps;
    die "No experiments after filtering out empty Description and test sets\n" if scalar(@exps) == 0;
    #  and set up @sets
    my %setsSeen = ();
    @sets = ();
    foreach my $exp (@exps) {
      my $set = $exp->{SetName};
      push @sets, $set unless exists $setsSeen{$set};
      $setsSeen{$set} = 1;
    }
  }

  # check the exps that are under consideration
  my %unknownMedia = ();
  my %unknownCompound = ();
  my @noMedia = ();
  foreach my $exp (@exps) {
    # ignore lines not filled out or dropped
    next if ($exp->{Index} eq "" && $exp->{SetName} eq "");
    if ($exp->{Drop}) {
      my $drop = $exp->{Drop};
      $drop =~ s/ //g;
      if (uc($drop) eq "TRUE" || uc($drop) eq "DROP") {
        print STDERR "Dropping $exp->{SetName} $exp->{Index} $exp->{Description}\n";
        $exp->{Drop} = "TRUE";
      } elsif (uc($drop) eq "NA") {
        $exp->{Drop} = "";
      } elsif (uc($drop) ne "FALSE") {
        print STDERR "Unknown drop code $drop for $exp->{SetName} $exp->{Index}\n";
        $exp->{Drop} = "";
      }
    }

    # Clean up media and warn if not a known media
    if (defined $exp->{Media} && $exp->{Media} ne "") {
      $exp->{Media} =~ s/^ +//;
      $exp->{Media} =~ s/ +$//;
      $unknownMedia{ $exp->{Media} } = 1 if !defined GetMediaComponents( $exp->{Media} );
    } else {
      push @noMedia, $exp->{SetName}.".".$exp->{Index};
    }

    # Clean up Group
    if (defined $exp->{Group}) {
      $exp->{Group} =~ s/^ +//;
      $exp->{Group} =~ s/ +$//;
      if (lc($exp->{Group}) eq "time0") {
        $exp->{Group} = "Time0";
      } else {
        $exp->{Group} = lc($exp->{Group}) unless $exp->{Group} eq "pH";
      }
    }

    # Clean up condition_1,2,3 and warn if not a known compound -- but skip this for Time0 samples
    # as sometimes they have condition set to Time0
    foreach my $field (qw{Condition_1 Condition_2 Condition_3 Condition_4}) {
      if (exists $exp->{$field} && $exp->{Group} ne "Time0") {
        $exp->{$field} = "" if lc($exp->{$field}) eq "none";
        $exp->{$field} = "" if $exp->{$field} eq "NA";
        if ($exp->{$field} ne "") {
          my $compound = FindCompound($exp->{$field});
          $compound = $exp->{$field} if !defined $compound
            && ( defined GetMediaComponents($exp->{$field})
                 || defined GetMixComponents($exp->{$field}) );
          if (!defined $compound) {
            $unknownCompound{$exp->{$field}} = 1
              unless $exp->{$field} =~ m/ exudates?$/ || $exp->{$field} =~ m/^supernatant; /i;
          } else {
            $exp->{$field} = $compound;
          }
        }
      }
    }
  }
  print STDERR join("\t","No media for",@noMedia)."\n" if scalar(@noMedia) > 0;
  print STDERR join("\t","Unknown media", sort keys %unknownMedia)."\n" if scalar(keys %unknownMedia) > 0;
  print STDERR join("\t","Unknown compound", sort keys %unknownCompound)."\n" if scalar(keys %unknownCompound) > 0;

  my @genes = &ReadTable($genesfile, qw{locusId scaffoldId sysName begin end strand});
  my %geneScaffolds = map { $_->{scaffoldId} => 1 } @genes;
  my %genesSorted = ();    # scaffold to list of genes sorted by begin
  foreach my $gene (@genes) {
    push @{ $genesSorted{$gene->{scaffoldId}} }, $gene;
  }
  foreach my $scaffold (keys %genesSorted) {
    my @sorted = sort { $a->{begin} <=> $b->{begin} } @{ $genesSorted{$scaffold} };
    $genesSorted{$scaffold} = \@sorted;
  }
  CheckGeneLocations(\%genesSorted); # writes to STDERR

  # Next, can we find all set files?
  my %setpre = (); # set name to list of prefixes having a poolcount and ideally also a colsum file
  my @pcfiles = ();             # list of candidate poolcount files
  opendir(DIR, $indir) || die $!;
  while (my $file = readdir(DIR)) {
    if ($file =~ m/^(.*)[.]poolcount$/) {
      push @pcfiles, $1;
    }
  }
  closedir(DIR);
  die "No *.poolcount files in $indir\n" if @pcfiles == 0;
  my %sets = map { $_ => 1 } @sets;
  my %setFiles = ();     # set to list of prefixes for poolcount files
  my %pcToSet = ();      # pcFile to set
  foreach my $set (@sets) {
    my @pcfileThis = ();
    foreach my $pcfile (@pcfiles) {
      if ($pcfile eq $set) {
        push @pcfileThis, $pcfile;
        $pcToSet{$pcfile} = $set;
      } elsif (!exists $sets{$pcfile} && lc(substr($pcfile, 0, length($set))) eq lc($set)) {
        my $postfix = substr($pcfile, length($set));
        if ($postfix eq "" || $postfix =~ m/^_?[a-zA-Z]$/ || $postfix =~ m/^_rep[a-zA-Z0-9]+$/ || $postfix =~ m/seq[a-zA-Z0-9]+$/ || $postfix =~ m/^_re$/) {
          print STDERR "Found extra file $pcfile for $set\n";
          push @pcfileThis, $pcfile;
          $pcToSet{$pcfile} = $set;
        }
      }
    }
    if (scalar(@pcfileThis) == 0) {
      print STDERR "No poolcount file for $set, skipping it\n";
      @exps = grep { $_->{SetName} ne $set } @exps;
      die "No experiments remaining, giving up\n" if scalar(@exps) == 0;
    } else {
      $setFiles{$set} = \@pcfileThis;
    }
  }
  @sets = grep { exists $setFiles{$_} } @sets;

  if (! $prespec_sets) {
    # look for set files that are not in the metadata
    foreach my $pcfile (@pcfiles) {
      if (!exists $pcToSet{$pcfile} && $pcfile !~ m/test/) {
        print STDERR "Warning: poolcount file with no metadata: $pcfile\n";
      }
    }
  }

  # Build list of experiments for each set
  my %setExps = ();
  foreach my $exp (@exps) {
    push @{ $setExps{$exp->{SetName}} }, $exp;
  }
  # And check that each index is unique for each set
  while (my ($set,$exps) = each %setExps) {
    my %indexSeen = ();
    foreach my $exp (@$exps) {
      my $index = $exp->{Index};
      if (exists $indexSeen{$index}) {
        die "Duplicate experiment entries for index $index in set $set";
      }
      #else
      $indexSeen{$index} = 1;
    }
  }

  print STDERR sprintf("%s %d experiments, %d genes, %d sets for $org\n",
                       defined $test ? "Test found" : "Processing",
                       scalar(@exps), scalar(@genes), scalar(@sets));

  # now, we need to map strains to genes, compute "f", combine counts, etc.
  # read each file in parallel
  my %setFh = (); # set to list of file handles reading counts for that set
  my %setIndex = (); # set to list of count fields, may be more than is in the metadata
  # each hash as index name => field number (0-based)
  my @metaFields = qw{barcode rcbarcode scaffold strand pos};
  my $nmeta = scalar(@metaFields);
  my ($BARCODE,$RCBARCODE,$SCAFFOLD,$STRAND,$POS) = (0..4);
  foreach my $set (sort keys %setFiles) {
    my $filelist = $setFiles{$set};
    my @noField; # experiments with no data
    foreach my $file (@$filelist) {
      my @fields = ReadColumnNames("$indir/$file.poolcount");
      die "Too few fields in $file" unless scalar(@fields) >= $nmeta;
      foreach my $i (0..($nmeta-1)) {
        die "Expected $metaFields[$i] but field is $fields[$i]" unless $fields[$i] eq $metaFields[$i];
      }
      my @indexNew = @fields[($nmeta..(scalar(@fields)-1))];
      if (!exists $setIndex{$set}) {
        # first file for this set
        $setIndex{$set} = \@indexNew;
        # In the experiment metadata, replace S0001 with S1 if necessary.
        # Warn about any Indexes in @exps for this set that are not present.
        my %index = map { $_ => 1 } @indexNew;
        foreach my $exp (@{ $setExps{$set} }) {
          die unless $exp->{SetName} eq $set;
          my $index = $exp->{Index};
          if (!exists $index{$index}) {
            my $index2 = $index;
            $index2 =~ s/([a-zA-Z])0+([0-9]+)$/$1$2/;
            if (exists $index{$index2}) {
              $index = $index2;
            } elsif ($index =~ m/^\d+$/) {
              # If the experiment metadata has Index = "23" instead of "S23", allow that
              my $index3 = "S".$index;
              $index = $index3 if exists $index{$index3};
            }
          }
          if (exists $index{$index}) {
            $exp->{Index} = $index;
          } else {
            push @noField, $exp->{Index};
          }
        }
      } else {
        # additional file for this set
        my @indexOld = @{ $setIndex{$set} };
        my %indexOld = map { $_ => 1 } @indexOld;
        my @indexNew = @fields[$nmeta..(scalar(@fields)-1)];
        my %indexNew = map { $_ => 1 } @indexNew;
        my %drop = (); # fields from setIndex to be ignored
        foreach my $index (@indexOld) {
          if (!exists $indexNew{$index}) {
            print STDERR "Warning: index $index is missing from in $file and will be ignored\n";
            $drop{$index} = 1;
          }
        }
        foreach my $index (@indexNew) {
          if (!exists $indexOld{$index}) {
            print STDERR "Warning: index $index is new in $file and will be ignored\n";
          }
        }
        my @indexKept = grep !exists $drop{$_}, @indexOld;
        $setIndex{$set} = \@indexKept;
      }
    }
    if (@noField > 0) {
      print STDERR "No data for " . scalar(@noField) . " fields for $set -- " . join(" ", @noField) . "\n";
    }
  }

  # Final list of experiments
  my %setIndexHash = (); # set => index => 1 if it has counts
  while (my ($set, $indexList) = each %setIndex) {
    foreach my $index (@$indexList) {
      $setIndexHash{$set}{$index} = 1;
    }
  }
  @exps = grep exists $setIndexHash{$_->{SetName}}{$_->{Index}}, @exps;

  if (defined $test) {
    print STDERR "All poolcount file headers verified with " . scalar(@exps) . " experiments\n";
    exit(0);
  }

  # barcode index (0-based) to list of metadata fields (the fields listed in @metaFields)
  my @meta = ();
  # Read the barcode counts files
  my %counts = (); # set => index => list of values, one per barcode
  while (my ($set, $filelist) = each %setFiles) {
    my $indexList = $setIndex{$set};
    foreach my $file (@$filelist) {
      open(my $fh, "<", "$indir/$file.poolcount") || die "Cannot read $indir/$file.poolcount";
      my $firstFile = scalar(@meta) == 0; # set if this is the first file for *any* set
      my $headerLine = <$fh>;
      chomp $headerLine;
      my @fields = split /\t/, $headerLine;
      die "Not enough fields in $file\n" unless scalar(@fields) >= $nmeta;
      my %fields = map { $fields[$_] => $_ } (0..(scalar(@fields)-1));
      my $iLine = 0;
      while(my $line = <$fh>) {
        chomp $line;
        my @values = split /\t/, $line;
        die "Wrong number of fields in $file\n"
          unless scalar(@values) == scalar(@fields);
        my @lineMeta = @values[0..($nmeta-1)];
        if ($firstFile) {
          $meta[$iLine] = \@lineMeta;
        } else {
          foreach my $i (0..($nmeta-1)) {
            die "Non-matching $meta[$i] = $lineMeta[$i] in $set line $iLine, expected $meta[$iLine][$i]\n"
              . "You may need to rerun combineBarSeq.pl with the new pool\n"
                . "or use bin/resortPoolCount.pl to make it match\n"
                  unless $meta[$iLine][$i] eq $lineMeta[$i];
          }
        }
        foreach my $index (@$indexList) {
          my $value = $values[$fields{$index}];
          die "No value for $index in $file at line $iLine" if !defined $value || $value eq "";
          if ($firstFile) {
            $counts{$set}{$index}[$iLine] = $value;
          } else {
            $counts{$set}{$index}[$iLine] += $value;
          }
        }
        $iLine++;
      } # end loop over lines
      if (!$firstFile && $iLine != scalar(@meta)) {
        die "Wrong number of lines in $file -- expected " . scalar(@meta) . " data lines\n";
      }
      close($fh) || die "Error reading $indir/$file.poolcount";
    } # end loop over files
  } # end loop over sets

  # Write out the all.poolcount file
  my @allfields = qw{barcode rcbarcode scaffold strand pos locusId f};
  foreach my $exp (@exps) {
    push @allfields, $exp->{SetName} . "." . $exp->{Index};
  }

  my $allfile = "$outdir/all.poolcount";
  print STDERR "Processed all the poolcount files; writing $allfile\n";
  open(ALL, ">", $allfile) || die "Cannot write to $allfile";
  print ALL join("\t",@allfields)."\n";
  my $nInGene = 0;
  foreach my $i (0..(scalar(@meta)-1)) {
    my $metaValues = $meta[$i];
    my ($locusId, $f) = &LocationToGene($metaValues->[$SCAFFOLD], $metaValues->[$POS], \%genesSorted);
    $nInGene++ if $locusId ne "";
    my @countsOut = ();
    foreach my $exp (@exps) {
      my $count = $counts{ $exp->{SetName} }{ $exp->{Index} }[$i];
      die unless defined $count;
      push @countsOut, $count;
    }
    print ALL join("\t", @$metaValues, $locusId, $f, @countsOut)."\n";
  }
  close(ALL) || die "Error writing to $allfile";
  print STDERR "Wrote data for " . scalar(@meta) . " barcodes ($nInGene in genes) to $allfile\n";
  die join("\n",
           "No insertions in genes!",
           "Please check that $genesfile contains genes,",
           "$poolfile contains strains,",
           "and that the scaffold identifiers match",
           "") if $nInGene == 0;

  # write outdir/exps, with only the exps for these sets and with non-empty Description
  my @expCols = ReadColumnNames($expsfile); # write the columns in a reasonable order
  my $expsout = "$outdir/exps";
  open(EXPSOUT, ">", $expsout) || die "Cannot write to $expsout";
  print EXPSOUT join("\t",@expCols)."\n";
  foreach my $exp (@exps) {
    my @out = map { $exp->{$_} } @expCols;
    print EXPSOUT join("\t", @out)."\n";
  }
  close(EXPSOUT) || die "Error writing to $expsout";
  print STDERR "Wrote " . scalar(@exps) . " experiments to $expsout\n";

  # write outdir/pool, outdir/genes
  system("cp", $poolfile, "$outdir/pool") == 0 || die $!;
  system("cp", $genesfile, "$outdir/genes") == 0 || die $!;

  # copy over the strain usage files if they exist unless FEBA_NO_STRAIN_USAGE is set
  if (-e "$indir/strainusage.barcodes" && ! $ENV{FEBA_NO_STRAIN_USAGE}) {
    die if ! -e "$indir/strainusage.genes";
    die if ! -e "$indir/strainusage.genes12";
    if ($indir ne $outdir) {
      print STDERR "Copying over strain usage files\n";
      system("cp",
             "$indir/strainusage.barcodes", "$indir/strainusage.genes", "$indir/strainusage.genes12",
             $outdir) == 0 || die $!;
    }
  } else {
    unlink("$outdir/strainusage.barcodes");
    unlink("$outdir/strainusage.genes");
    unlink("$outdir/strainusage.genes12");
  }

  my $Rcmd = "Rscript $Rscript $org $outdir $Bin/.. > $outdir/log";
  if (defined $noR) {
    print STDERR "Skipping the R step: $Rcmd\n";
  } else {
    print STDERR "Running R: $Rcmd\n";
    unlink("$outdir/.FEBA.success");
    system($Rcmd);
    die "R failed, see $outdir/log\n" unless -e "$outdir/.FEBA.success";
    print STDERR "Success for $org, please visit $outdir/index.html\n";
  }
}
