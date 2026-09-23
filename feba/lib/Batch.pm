#
#######################################################
## Batch.pm -- utilities for batch Fitness BLAST
##
## Copyright (c) 2016 All Right Reserved by UC Berkeley
##
## Authors:
## Morgan Price
#######################################################

package Batch;
require Exporter;

use strict;
use warnings;

use DBI;
use Carp;

# Checks the validity of the argument as well
# Returns undef if it cannot load the database
sub get_batch_dbh {
    my ($jobId) = @_;
    die "Invalid jobId" unless $jobId =~ m/^[a-zA-Z0-9_:-]+$/;
    die "No directory for jobId $jobId" unless -d "../job_data/$jobId";
    my $dbfile = "../job_data/$jobId/job.db";
    return undef unless -s $dbfile;
    my $bdb = DBI->connect("dbi:SQLite:dbname=$dbfile","","");
    return undef unless $bdb;
    my $query = $bdb->selectrow_arrayref("SELECT * FROM Query limit 1;");
    return undef unless @$query;
    return $bdb;
}

sub get_job_name {
    my ($jobId,$bdb) = @_;
    my ($id1) = $bdb->selectrow_array("SELECT queryId FROM Query ORDER BY queryId limit 1;");
    my ($idL) = $bdb->selectrow_array("SELECT queryId FROM Query ORDER BY queryId DESC limit 1;");
    my ($nProt) = $bdb->selectrow_array("SELECT COUNT(*) FROM Query;");
    return "$id1..$idL ($nProt)";
}

1;
