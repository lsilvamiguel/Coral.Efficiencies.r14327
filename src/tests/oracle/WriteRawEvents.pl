#!/usr/local/bin/perl -w

require 5.004;

use strict;
use diagnostics;
use DBI;

$ENV{ORACLE_CERN}        = "/afs/cern.ch/project/oracle";
$ENV{ORACLE_MOUNT}       = "$ENV{ORACLE_CERN}/\@sys";
$ENV{TNS_ADMIN}          = "$ENV{ORACLE_CERN}/admin";
$ENV{ORACLE_HOME}        = "$ENV{ORACLE_MOUNT}/920";
$ENV{ORA_NLS33}          = "$ENV{ORACLE_HOME}/ocommon/nls/admin/data";
$ENV{DIR_ORACLE}         = "$ENV{ORACLE_HOME}";
$ENV{NLS_LANG}           = "american_america.WE8ISO8859P9";
$ENV{NLS_DATE_FORMAT}    = "DD-MON-FXYYYY";
$ENV{STAGE_HOST}         = "stagecompass";
$ENV{PATH}               = "$ENV{DIR_ORACLE}/bin:$ENV{PATH}";
$ENV{STAGE_POOL_EXCLUDE} = "/afs/cern.ch/user/o/objsrvvy/castor/.StageHelperPoolsExclude.inf";

unless( defined($ENV{LD_LIBRARY_PATH}) ) {
    $ENV{LD_LIBRARY_PATH} = "$ENV{DIR_ORACLE}/lib";
}
else {
    $ENV{LD_LIBRARY_PATH} = "$ENV{DIR_ORACLE}/lib:$ENV{LD_LIBRARY_PATH}";
}

#if($#ARGV != 0) {
#    die("Usage: WriteRawEvents.pl <events list file name>\n");
#}


#my $event_list = $ARGV[0];

my $castor_dir         = "/afs/cern.ch/compass/scratch/d17/objsrvvy/J_Psi/";
my $work_dir           = "/afs/cern.ch/compass/scratch/d17/objsrvvy/J_Psi";
my $exec_program       = "$work_dir/Linux/WriteRawEvents";
my $event_list         = "$work_dir/dumpjpsi.dat";
my $exclude_event_list = "$work_dir/J_Psi.list";
my $complete_file_name = "$work_dir/events_to_read.list";


unless(open(EXCLUDE_LIST, $exclude_event_list)) {
    die "Cannot open file: $exclude_event_list\n";
}

my %run_list;
my @event_list_to_read;
my @event_list_to_exclude;
my @complete_event_list;

unless(open(LIST_FILE, $event_list)) {
    die "Cannot open file: $event_list\n";
}

while(<LIST_FILE>) {
    if($_ =~ /(\d+)\s+(\d+)\s+(\d+)/){
	push(@event_list_to_read,"$1 $2 $3");
    }
}

close LIST_FILE;
unless(open(EXCLUDE_LIST, $exclude_event_list)) {
    die "Cannot open file: $exclude_event_list\n";
}

while(<EXCLUDE_LIST>) {
    if($_ =~ /(\d+)\s+(\d+)\s+(\d+)/){
	push(@event_list_to_exclude,"$1 $2 $3");
    }
}

close EXCLUDE_LIST;

unless(open(COMPLETE_LIST, "> $complete_file_name")) {
    die "Cannot create file: $complete_file_name\n";
}

foreach my $etr (@event_list_to_read) {
    my $to_be_exclude = 0;
    foreach my $ete (@event_list_to_exclude) {
	if($etr eq $ete) {
	    $to_be_exclude = 1;
	    last;
	}
    }
    unless($to_be_exclude) {
	$etr =~ /^(\d+)/;
	my $run_number = $1;
	print COMPLETE_LIST "$etr\n";
	unless(defined $run_list{$run_number}) {
	    $run_list{$run_number} = 1;
	}
	else {
	    $run_list{$run_number} += 1;
	}
    }
}

close COMPLETE_LIST;

foreach my $run_number (keys %run_list) {
    print("Write $run_list{$run_number} event(s) from run number $run_number\n");
    my $ret_val = system($exec_program, $run_number, $complete_file_name, $castor_dir);
    if($ret_val == 0) {
	print "OK!\n";
    }
    else {
	print "ERROR: error code = $ret_val\n";
    }
}

