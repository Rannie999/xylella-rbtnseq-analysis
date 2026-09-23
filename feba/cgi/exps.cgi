#!/usr/bin/perl -w
#######################################################
## exps.cgi -- search for experiments, or list all in one organism
##
## Copyright (c) 2015 University of California
##
## Authors:
## Morgan Price and Victoria Lo
#######################################################
#
# Key parameters: orgId and query, for which organism and to look up experiments
#	At least one must be meaningful (present and not empty)
# OR, orgId and expGroup
# OR, orgId and expGroup AND condition1 (condition1 may be empty, but is still used to restrict the results)
# OR, orgId and date
# OR, orgId and media
# If orgId only is specified, forwards to the organism page

use strict;
use CGI qw(:standard Vars);
use CGI::Carp qw(warningsToBrowser fatalsToBrowser);
use DBI;

use lib "../lib";
use Utils;
use URI::Escape;
use HTML::Entities;

my $cgi=CGI->new;

my $orgId = $cgi->param('orgId');
$orgId = "" if !defined $orgId;

my $expSpec = $cgi->param('query');
$expSpec = "" if !defined $expSpec;
$expSpec =~ s/ +$//;
$expSpec =~ s/^ +$//;
$expSpec =~s/[\"\n\r\;\\]//g; # ' is allowed, handled in matching_exps

my $expGroup = $cgi->param('expGroup');
my $condition1 = $cgi->param('condition1');
my $dateSpec = $cgi->param('date');
my $mediaSpec = $cgi->param('media');

my $dbh = Utils::get_dbh();

# Redirect to org.cgi if displaying all experiments from one organism,
# or to orgAll.cgi if displaying all experiments
if ($orgId ne "" && !defined $expGroup && !defined $dateSpec && !defined $mediaSpec
    && ($cgi->param("All experiments") || $expSpec eq "")) {
    print redirect(-url=>"org.cgi?orgId=$orgId");
    exit(0);
} elsif ($orgId eq "" && !$expSpec && !$expGroup && !$condition1) {
  print redirect(-url=>"orgAll.cgi");
  exit(0);
}

my $specShow = encode_entities($expSpec);
if ($specShow eq "") {
    $specShow = defined $condition1 ? join(" ", encode_entities($expGroup), encode_entities($condition1))
      : encode_entities($expGroup);
    $specShow = "" if !defined $specShow;
    $specShow = "date " . encode_entities($dateSpec) if defined $dateSpec && $dateSpec ne "";
    $specShow = "media " . encode_entities($mediaSpec) if defined $mediaSpec && $mediaSpec ne "";
}
my $start = Utils::start_page("Experiments for $specShow");
$expSpec = "" if $cgi->param("All experiments");

print $cgi->header, $start, '<div id="ntcontent">';


# Make sure orgId parameter is safe
my $orginfo = Utils::orginfo($dbh);
Utils::fail($cgi, "Unknown organism: $orgId") unless $orgId eq "" || exists $orginfo->{$orgId};

my $exps;
if ($orgId ne "" && defined $expGroup && !defined $condition1) {
  $exps = $dbh->selectall_arrayref(qq{SELECT * from Experiment WHERE orgId = ? AND expGroup = ?
                                      ORDER BY expGroup, condition_1, expName},
                                   { Slice => {} },
                                   $orgId, $expGroup);
} elsif ($orgId ne "" && defined $expGroup && defined $condition1) {
  $exps = $dbh->selectall_arrayref(qq{SELECT * from Experiment WHERE orgId = ? AND expGroup = ? AND condition_1 = ?
                                      ORDER BY expGroup, condition_1, expName},
				    { Slice => {} },
				    $orgId, $expGroup, $condition1);
    Utils::fail($cgi, "No experiments for specified organism, group, and condition") if scalar(@$exps) == 0;
} elsif ($orgId eq "" && defined $expGroup && defined $condition1) {
  $exps = $dbh->selectall_arrayref(qq{SELECT * from Experiment JOIN Organism Using (orgId) 
                                      WHERE expGroup = ? AND condition_1 = ?
                                      ORDER BY genus, species, strain, expGroup, condition_1, concentration_1, expDesc;},
				    { Slice => {} },
				    $expGroup, $condition1);
    Utils::fail($cgi, "No experiments for specified organism, group, and condition") if scalar(@$exps) == 0;
} elsif ($orgId eq "" && defined $expSpec) {
  $exps = Utils::matching_exps($dbh, "", $expSpec);

} elsif ($orgId ne "" && defined $dateSpec && $dateSpec ne "") {
  $exps = $dbh->selectall_arrayref("SELECT * from Experiment WHERE orgId = ? AND dateStarted = ?"
                                   . " ORDER BY expGroup,condition_1,expDesc",
                                   { Slice => {} }, $orgId, $dateSpec);
} elsif ($orgId ne "" && defined $mediaSpec && $mediaSpec ne "") {
  $exps = $dbh->selectall_arrayref("SELECT * from Experiment WHERE orgId = ? AND media = ?"
                                   . " ORDER BY expGroup,condition_1,expDesc",
                                   { Slice => {} }, $orgId, $mediaSpec);
} else {
    $exps = Utils::matching_exps($dbh, $orgId, $expSpec);
}

# sort experiments by organism, or by lc(condition)
if ($orgId eq "") {
    my @exps = sort { $orginfo->{$a->{orgId}}{genome} cmp $orginfo->{$b->{orgId}}{genome} } @$exps;
    $exps = \@exps;
}
if ($orgId ne "" && defined $expGroup) {
    my @exps = sort { lc($a->{condition_1}) cmp lc($b->{condition_1}) } @$exps;
    $exps = \@exps;
}

if (@$exps == 0) {
  my $msg = qq{No experiment found matching "} . encode_entities($expSpec) . qq{"};
  $msg .= " in $orginfo->{$orgId}{genome}" if $orgId;
  print $cgi->h3($msg);
} else {
  my $heading = "Experiments";
  $heading .= " in ". $cgi->a({href => "org.cgi?orgId=$orgId"}, "$orginfo->{$orgId}{genome}") if $orgId ne "";
  if ($expSpec ne "") {
      $heading .= qq{ matching "$expSpec"};
  } elsif ($specShow ne "") {
      $heading .= " for $specShow";
  }
  print $cgi->h2($heading);

  my $exp1 = $exps->[0];
  if ($exp1->{expGroup} ne ""
      && ($expSpec ne "" || (defined $expGroup && defined $condition1))) {
      print p(a( { -href => "orthCond.cgi?expGroup=" . uri_escape($exp1->{expGroup})
                       . "&condition1=" . uri_escape($exp1->{condition_1}) },
                 "Or see specific phenotypes for $exp1->{expGroup} $exp1->{condition_1} across organisms"));
  }
}

# experiment object to link, including hover text
sub ExpLink($) {
  my ($exp) = @_;
  return a({ -href => "exp.cgi?orgId=$exp->{orgId}&expName=$exp->{expName}",
             -title => "$exp->{expName} from $exp->{dateStarted} by $exp->{person}" },
           $exp->{expName});
}

my @trows = ();
# remove organism column if already specified
if (defined $expGroup && $orgId ne "" && !defined $condition1){
   push @trows, $cgi->Tr({-valign => "top"}, $cgi->th(['Name', 'Group', 'Condition', 'Description' ]));
  foreach my $row (@$exps) {
    my @cond = map $row->{"condition_".$_}, (1..4);
    @cond = grep $_ ne "", @cond;
    push @trows, $cgi->Tr({-valign => "top"},
                          $cgi->td([ ExpLink($row),
                                     $row->{expGroup},
                                     join("; ", @cond),
                                     $row->{expDesc} ]));
    }
} else {
  push @trows, $cgi->Tr({-valign => "top"}, $cgi->th([ 'Organism', 'Name', 'Group', 'Condition', 'Description' ]));
  foreach my $row (@$exps) {
      push @trows, $cgi->Tr({-valign => "top"},
             $cgi->td([ $cgi->a({ -href => "org.cgi?orgId=". $row->{orgId} },
                                $orginfo->{$row->{orgId}}{genome}),
                        ExpLink($row),
                        $row->{expGroup}, $row->{condition_1}, $row->{expDesc} ]));
  }
}

print table({cellspacing => 0, cellpadding => 3}, @trows),
    "<br><br>";
print p("Or see all specific phenotypes in ",
        a({href => "spec.cgi?orgId=$orgId&expGroup=" . uri_escape($expGroup)}, $expGroup),
        " experiments for $orginfo->{$orgId}{genome}")
    if defined $expGroup && @$exps > 0 && $orgId ne "";

$dbh->disconnect();
Utils::endHtml($cgi);
