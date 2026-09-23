#!/usr/bin/perl -w

# Given the results of all-vs-all BLASTp on a database with sequence identifiers of the form
# organism:locus, make a table of "orthologs"

# Throughout, "locusId" is of the form organism:locus, not the second half only

use strict;
use Getopt::Long;
if ($ENV{MEMUSAGE}) {
  use Devel::Size qw{total_size};
}
use FindBin qw{$RealBin};
use lib "$RealBin/../lib";
use FEBA_Utils qw{ReadFasta};

sub ReadBlastp($$$);

my $minCoverage = 0.8; # default requirement for coverage both ways
my $minRatio = 0.25; # default minimum bit score ratio and %identity
my $usage = <<END
Usage: bbh.pl [ -coverage $minCoverage ] [ -ratio $minRatio ]
              -files blastp_files -faa file1 ... fileN -out out

blastp_files is a list of output files from blastp, in tab-delimited format,
  and locus identifiers must have the form orgId:locusId

Each fasta file must have identifiers of the form orgId:locusId

Writes to out.scores, with one line per pair of BBHs and fields tax1,
locus1, tax2, locus2, and ratio, where ratio is the bit score divided
by locus1's self-score.
END
  ;

my ($SUBJECT, $BITS) = (0,1); # indexes into the best hit data structure

{
    my $compare = 0;
    my $listFile;
    my @faaFiles;
    my $prefix = undef;

    (GetOptions('coverage=f' => \$minCoverage,
                'ratio=f' => \$minRatio,
		'files=s' => \$listFile,
                'faa=s{1,}' => \@faaFiles,
		'out=s' => \$prefix)
     && @ARGV==0
     && defined $listFile && defined $prefix) || die $usage;
    die "Must specify faa files\n$usage" unless @faaFiles > 0;

    my @blastpFiles = ();
    open (my $fhList, "<", $listFile) || die "Cannot read $listFile\n";
    while (my $line = <$fhList>) {
      chomp $line;
      die "No such file: $line\n" unless -e $line;
      push @blastpFiles, $line;
    }
    close($fhList) || die "Error reading $listFile";

    my %seqs = ();
    foreach my $faaFile (@faaFiles) {
      my $hash = ReadFasta($faaFile); # hash of name => sequence
      while (my ($name, $seq) = each %$hash) {
        die "$name in $faaFile is duplicated from a previous file"
          if exists $seqs{$name};
        die "Invalid sequence name $name in $faaFile\n"
          unless $name =~ m/^[a-zA-Z0-9._-]+:[a-zA-Z90-9._-]+$/;
        $seqs{$name} = $seq;
      }
    }
    print STDERR "Read " . scalar(keys %seqs) . " sequences\n";

    print STDERR "Reading " . scalar(@blastpFiles) . " blastp files\n";
    # best hits: locusId => subject orgId => [subject,score,qbegin,qend,sbegin,send]
    # (Ensures that there is no tie for the best hit, and that the best hit meets the coverage
    #  requirement.)
    my %bestHit = ();
    foreach my $file (@blastpFiles) {
      ReadBlastp($file, \%seqs, \%bestHit);
    }
    print STDERR "Read hits for " . scalar(keys %bestHit) . " proteins\n";

    open(my $fhOut, ">", "$prefix.scores") || die "Cannot write to $prefix.scores\n";
    print $fhOut join("\t", qw{tax1 locus1 tax2 locus2 ratio})."\n";
    foreach my $query (sort keys %bestHit) {
      my ($qOrgId, $qLocus) = split /:/, $query;
      my $hash = $bestHit{$query};
      # Usually, the lack of a best hit in its own genome means that the protein has an identical paralog,
      # in which case there would not be an BBHs anyway.
      next unless exists $hash->{$qOrgId};
      my $selfHit = $hash->{$qOrgId};
      my $selfBits = $selfHit->[$BITS];
      foreach my $sOrgId (sort keys %$hash) {
        next if $sOrgId eq $qOrgId;
        my ($subject,$bits,$qBeg,$qEnd,$sBeg,$sEnd) = @{ $hash->{$sOrgId} };
        if (exists $bestHit{$subject}{$qOrgId} && $bestHit{$subject}{$qOrgId}[$SUBJECT] eq $query) {
          # bidirectional best hit
          my $ratio = $bits / $selfBits;
          print $fhOut join("\t", $qOrgId, $query, $sOrgId, $subject, $bits / $selfBits) . "\n"
            if $ratio >= $minRatio;
        }
      }
    }
    close($fhOut) || die "Error writing to $prefix.scores";
    if (exists $ENV{MEMUSAGE}) {
      print STDERR sprintf("Memory usage for best hit data structure: %.1f MB\n",
                           total_size(\%bestHit) / 1024**2);
    }
    print STDERR "Wrote $prefix.scores\n";
}

# Reads the file, and uses the sequence information to ignore low-coverage hits and weaker hits
# (within a taxon).
# Saves high-coverage top hits to $hits->{$query}{$subjectOrgId}
sub ReadBlastp($$$) {
  my ($file, $seqs, $bestHit) = @_;
  my %bestBits = (); # subject => hit taxId => best score
  open(my $fh, "<", $file) || die "Cannot read $file\n";
  while (my $line = <$fh>) {
    chomp $line;
    my ($query, $subject,
        $identity, $aln_length, $mismatch, $gaps,
        $qBeg, $qEnd, $sBeg, $sEnd, $evalue, $bits) = split /\t/, $line;
    die "Wrong number of columns in\n$line\n..." unless defined $bits && $bits =~ m/^ *[0-9.e+]+$/;
    $bits =~ s/^ +//;
    my ($qOrgId,$query2) = split /:/, $query; # to get taxId and locusId
    my ($sOrgId,$subject2) = split /:/, $subject;
    die "Invalid query $query" unless defined $query2;
    die "Invalid subject $subject" unless defined $subject2;
    die "No sequence for query $query in $file\n" unless exists $seqs->{$query};
    die "No sequence for subject $subject in $file\n" unless exists $seqs->{$subject};
    if (exists $bestBits{$query}{$sOrgId}) {
      # If there is a tie for best hit, unset the best hit
      if ($bits == $bestBits{$query}{$sOrgId}
          && exists $bestHit->{$query}{$sOrgId}
          && $bestHit->{$query}{$sOrgId}[$SUBJECT] ne $subject) {
        delete $bestHit->{$query}{$sOrgId};
        next;
      }
      # Ignore unless bit score is better
      next unless $bits > $bestBits{$query}{$sOrgId};
    }
    $bestBits{$query}{$sOrgId} = $bits;
    my $qLen = length($seqs->{$query});
    my $cov = ($qEnd - $qBeg + 1) / $qLen;
    $bestHit->{$query}{$sOrgId} = [$subject,$bits,$qBeg,$qEnd,$sBeg,$sEnd]
      if $cov >= $minCoverage;
  }
  close($fh) || die "Error reading $file";
}
