
# FEBA_Utils.pm -- Utilities for the command-line scripts.
# Copyright (c) 2017 All Rights Reserved by UC Berkeley
# Authors: Morgan Price at LBL

package FEBA_Utils;
require Exporter;
use strict;
use File::stat;
use File::Basename;
our (@ISA,@EXPORT);
@ISA = qw(Exporter);
@EXPORT = qw( ReadTable ReadColumnNames ReadFasta NewerThan ReadFastaDesc ReadFastaEntry reverseComplement
              sqlite_quote FastqByIndex CodesByIndex FastqToIndex );

# filename and list of required fields => list of hashes, each with field->value
# The list can be a single name or a reference to a list
# (It used to be an actual list; not sure how that went wrong with perl prototypes?)
sub ReadTable($*) {
    my ($filename,@required) = @_;
    if (scalar(@required) == 1 && ref $required[0]) {
        @required = @{ $required[0] };
    }
    open(IN, "<", $filename) || die "Cannot read $filename";
    my $headerLine = <IN>;
    $headerLine =~ s/[\r\n]+$//; # for DOS
    # Check for Mac style files -- these are not supported, but give a useful error
    die "Tab-delimited input file $filename is a Mac-style text file, which is not supported\n"
        . "Use\ndos2unix -c mac $filename\nto convert it to a Unix text file.\n"
        if $headerLine =~ m/\t/ && $headerLine =~ m/\r/;
    my @cols = split /\t/, $headerLine, -1;
    my %cols = map { $cols[$_] => $_ } (0..(scalar(@cols)-1));
    foreach my $field (@required) {
	die "No field $field in $filename" unless exists $cols{$field};
    }
    my @rows = ();
    while(my $line = <IN>) {
	$line =~ s/[\r\n]+$//;
	my @F = split /\t/, $line, -1;
	die "Wrong number of columns in:\n$line\nin $filename"
	    unless scalar(@F) == scalar(@cols);
	my %row = map { $cols[$_] => $F[$_] } (0..(scalar(@cols)-1));
	push @rows, \%row;
    }
    close(IN) || die "Error reading $filename";
    return @rows;
}

# filename to list of column names
sub ReadColumnNames($) {
    my ($filename) = @_;
    open(IN, "<", $filename) || die "Cannot read $filename";
    my $line = <IN>;
    close(IN) || die "Error reading $filename";

    $line =~ s/[\r\n]+$//; # for DOS
    my @cols = split /\t/, $line;
    return @cols;
}

# returns a reference to a hash of name => sequence
sub ReadFasta ($) {
    my ($filename) = @_;
    open(IN, "<", $filename) || die "Cannot read $filename";
    my %seqs = ();
    my $name = undef;
    while(<IN>) {
	chomp;
	if (m/>(\S+)/) {
	    die "Duplicate entry for $1 in filename" if exists $seqs{$1};
	    $name = $1;
	} else {
	    die "Unexpected line $_" unless defined $name;
	    $seqs{$name} .= $_;
	}
    }
    close(IN) || die "Error reading $filename";
    return(\%seqs);
}

# Returns a hash containing either "error"
# or hashes of "seq", "desc", and "len"
sub ReadFastaDesc($) {
    my ($file) = @_;
    my %seq = ();
    my %desc = ();
    my $name = undef;
    open(FAA, "<", $file) || return('error' => "Cannot read $file" );
    while(<FAA>) {
        s/[\r\n]+$//;
        if (m/^>(.*)$/) {
            my $header = $1;
            if ($header =~ m/^(\S+)\s+(\S.*)$/) {
                $name = $1;
                $desc{$name} = $2;
            } else {
                return('error' => "bad header for sequence:\n$header\n") unless $header =~ m/^\S+$/;
                $name = $header;
                $desc{$name} = $header;
            }
            return('error' => "Duplicate sequence id:\n$name\n") if exists $seq{$name};
            $seq{$name} = "";
        } else {
            return('error' => "sequence before header:\n$_\n") unless defined $name;
            s/\s//g;
            $seq{$name} .= $_;
        }
    }
    close(FAA) || return('error' => "Error reading $file");
    my %len = ();
    while (my ($name,$seq) = each %seq) {
        $len{$name} = length($seq);
        return('error' => "No sequence for id:\n$name\n") if ($len{$name} == 0);
    }
    return("desc" => \%desc, "seq" => \%seq, "len" => \%len);
}

# Read one entry at a time from a fasta file
# The first argument is a hash to keep track of saved state, i.e.:
#   my $state = {};
#   while(my ($header,$sequence) = ReadFastaEntry($fh,$state)) { ... }
# (header will have the ">" removed)
sub ReadFastaEntry {
  my ($fh, $state) = @_;
  die unless ref $state;
  return () if exists $state->{DONE}; # end state
  # initialization
  if (!defined $state->{header}) {
    $state->{header} = "";
    $state->{sequence} = "";
  }
  while (my $line = <$fh>) {
    chomp $line;
    if ($line =~ m/^>(.*)/) {
      my $old_header = $state->{"header"};
      my $old_sequence = $state->{"sequence"};
      $state->{"header"} = $1;
      die "Empty header in $line" if $state->{header} eq "";
      $state->{"sequence"} = "";
      return ($old_header, $old_sequence) if $old_header ne "";
    } else {
      die "Unexpected sequence with no header" if $state->{"header"} eq "";
      $line =~ s/ //g;
      $line = uc($line);
      # allow - or . as used in alignments and * as used for stop codons
      die "Invalid sequence line $line" unless $line =~ m/^[A-Z*.-]*$/;
      $state->{sequence} .= $line;
    }
  }
  # reached end of file
  $state->{DONE} = 1;
  return () if $state->{header} eq ""; # empty file
  return ($state->{header}, $state->{sequence});
}

sub NewerThan($$) {
    my ($file1, $file2) = @_;
    die "Invalid arguments to NewerThan" unless defined $file1 && defined $file2;
    die "No such file: $file2" unless -e $file2;
    return 0 unless -e $file1;
    return stat($file1)->mtime >= stat($file2)->mtime ? 1 : 0;
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

# This is for uploading to sqlite3. 
# sqlite3 expects CVS format, not exactly tab delimited format
# So, need to replace any " with "" and surround the field with quotes.
sub sqlite_quote($) {
  my ($in) = @_;
  return $in unless $in =~ m/"/;
  $in =~ s/"/""/g;
  return '"' . $in . '"';
}

# Despite the name, this should also work on a .codes filename
sub FastqToIndex($$$$) {
  my ($fq, $bynumber, $bywell, $offset) = @_;
  my $base = basename($fq);
  my $sampleNum;
  if ($bynumber) {
    # example: FEBA_BS_425_IT014_S494_L006_R1_001.fastq.gz
    # example: FEBA_BS_448_289_S1_R1.fastq.gz
    # example: YH_BS_run-001_100_R1.fastq.gz
    die "No index number in file $fq in bynumber mode\n"
      unless $base =~ m/_(\d+)_S\d+_[LR]\d+[_.]/
      || $base =~ m/_(\d+)_R1[.]/;
    $sampleNum = $1;
  } elsif ($bywell) {
    # example: A01_16_L001_R1_001.fastq.gz
    # example: BbD1_1_C10_C10_S34_R1_001_merge.fastq.gz
    $base =~ m/^([A-Z])(\d+)_/
      || $base =~ m/_([A-Z])(\d+)_S\d+_/
      || $base =~ m/-([A-Z])(\d+)_S\d+_/
      || die "Cannot find well in file $fq in bywell mode\n";
    my ($row,$column) = ($1,$2);
    $row = ord($row) - ord("A");
    $sampleNum = $row * 12 + $column;
    die "Invalid sample number $sampleNum from file name $fq for well $row$column in bywell mode\n"
      if $sampleNum < 1 || $sampleNum > 96;
  }
  my $index;
  if (defined $sampleNum) {
    die "Invalid sampleNum $sampleNum -- must be greater than offset $offset"
      unless $sampleNum > $offset;
    $sampleNum -= $offset;
    if ($bynumber) {
      $index = "S".$sampleNum;
    } else {
      $index = sprintf("IT%03d", $sampleNum);
    }
  }
  if (!defined $index) {
    # Look for the ITnnn index number within the file name
    my @pieces = split /[._-]/, $base;
    my @indexes = grep m/^IT\d+$/, @pieces;
    # sometimes ITO (capital O) not IT0 (numeral zero)
    my @indexes2 = grep m/^ITO\d+$/, @pieces;

    if (@indexes == 1) {
      $index = $indexes[0];
    } elsif (@indexes2 == 1) {
      # e.g. HiSeq_barcodes/FEBA_BS_81/Sample_FEBA_BS_81_ITO88/FEBA_BS_81_ITO88_GGCCTG_L001_R1_001.fastq.gz
      $index = $indexes2[0];
      $index =~ s/ITO/IT0/;
    } elsif ($base =~ m/_(\d+)_[ACGT][ACGT][ACGT][ACGT][ACGT][ACGT]_/
             || $base =~ m/_Index(\d+)_[ACGT][ACGT][ACGT][ACGT][ACGT][ACGT]_/i
             || $base =~ m/_Index(\d+)_S\d+_/
             || $fq =~ m!_IT(\d\d\d)[/_.]!) {
      # e.g. FEBA_BS_60_10_TAGCTT_L001_R1_001.fastq.gz
      # e.g. FEBA_BS_125_Index10_S10_L001_R1_001.fastq.gz
      # e.g. FEBA_BS_195_IT001/FEBA_BS_195_S97_L002_R1_001.fastq.gz
      $sampleNum = $1;
    } elsif ($base =~ m/_(\d+)_S\d+_L\d+_/) {
      # e.g. FEBA_BS_117_24_S120_L002_R1_001.fastq.gz is IT024
      $sampleNum = $1;
    } elsif ($base =~ m/^(\d+)_S\d+_L\d+_R\d+_/) {
      # e.g. 18_S1_L001_R1_001.fastq.gz is IT018
      $sampleNum = $1;
    } elsif ($base =~ m/^[A-Z]*\d*_S(\d+)_*L*\d*_R\d+/) {
      # e.g. A10_S106_L003_R1_001.fastq.gz
      # e.g. A10_S10_R1_001.fastq.gz
      $sampleNum = $1;
    } elsif ($base =~ m/_\d+_S(\d+)_R\d+/) {
      # e.g. FEBA_BS_434_100_S100_R1.fastq.gz
      $sampleNum = $1;
    }
  }
  if (defined $sampleNum && !defined $index) {
    die "Invalid sampleNum $sampleNum -- must be greater than offset $offset"
      unless $sampleNum > $offset;
    $sampleNum -= $offset;
    $index = sprintf("IT%03d", $sampleNum);
  }
  $index = "IT0" . $1 if defined $index && $index =~ m/^IT(\d\d)$/;
  return $index;
}

# Like FastqByIndex but returns codes files
# Takes as an argumnte a hash, with the input directory in 'in'
# Optional arguments:
#   bynumber -- 1 if look for sammples by index number (usually for BarSeq4)
#   bywell -- 1 if look for samples by well, i.e. A01 not IT001
#   offset -- i.e., set to 96 if IT001 is named S97
# Returns a reference to a hash of index to list of codes files
sub CodesByIndex {
  my (%params) = @_;
  die "No input directory" unless defined $params{in};
  my $inDir = $params{in};
  die "Not a directory: $inDir" unless -d $inDir;
  my @codes = glob("$inDir/*.codes");
  die "No codes file found for $inDir\n" if @codes == 0;
  my %indexToFile = (); # index to list of fastq files
  foreach my $file (@codes) {
    my $index = FastqToIndex($file, $params{bynumber}, $params{bywell}, $params{offset} || 0);
    die "Cannot parse index from filename\n$file\n" unless defined $index;
    push @{ $indexToFile{$index} }, $file;
  }
  return \%indexToFile;
}

# Takes as an argumnte a hash, with the input directory in 'in'
# Optional arguments:
#   bynumber -- 1 if look for sammples by index number (usually for BarSeq4)
#   bywell -- 1 if look for samples by well, i.e. A01 not IT001
#   offset -- i.e., set to 96 if IT001 is named S97
# Returns a reference to a hash of index to list of fastqs
sub FastqByIndex {
  my (%params) = @_;
  die "No input directory" unless defined $params{in};
  my $inDir = $params{in};
  die "Not a directory: $inDir" unless -d $inDir;
  my @globs = ("$inDir/*.fastq.gz",
               "$inDir/*.fq.gz",
               "$inDir/*/*.fastq.gz",
               "$inDir/*/*.fq.gz");
  my @fqs;
  foreach my $glob (@globs) {
    @fqs = glob($glob);
    # Remove Undetermined reads before deciding if subdirectories need to be searched
    @fqs = grep !m/Undetermined_/i, @fqs;
    last if @fqs > 0;
  }
  die "No fastq files found for $inDir\n" if @fqs == 0;
  # Remove 2nd reads, if any
  @fqs = grep !m/_2[.]f[astq]+[.]gz$/, @fqs;
  @fqs = grep !m/_R2_\d+[.]f[astq]+[.]gz$/, @fqs;

  my %indexToFq = (); # index to list of fastq files
  foreach my $fq (@fqs) {
    next if $fq =~ m/_R2[.]/; # ignore any second reads
    my $index = FastqToIndex($fq, $params{bynumber}, $params{bywell}, $params{offset} || 0);
    die "Cannot parse index from filename\n$fq\n" unless defined $index;
    push @{ $indexToFq{$index} }, $fq;
  }
  return \%indexToFq;
}
1;
