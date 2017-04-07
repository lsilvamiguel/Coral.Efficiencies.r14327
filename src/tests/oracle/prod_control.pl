#!/usr/local/bin/perl -w

require 5.004;

use strict;
use diagnostics;
use Time::Local;
use IO::Handle;


#if($#ARGV != 0) {
#    die("Usage: prod_control.pl <chunk name>\n");
#}

my $slot_number   = $ENV{DST_SLOT_NUMBER};
my $chunk_name    = "$ENV{DST_RUN_NUMBER}-$ENV{DST_CHUNK_NUMBER}";
my $log_file_dir  = "/afs/cern.ch/compass/scratch/d17/objsrvvy/production/logs";
my $log_file_name = "$log_file_dir/$chunk_name-$slot_number.log";

open(LOG_OUTPUT,"> $log_file_name") || die "Cannot create LOG file: $log_file_name\n";
#LOG_OUTPUT->autoflush(1);
select((select(LOG_OUTPUT), $| = 1)[0]);
$| = 1;

my $job_start_time    = 0;
my $event_start_time  = 0;
my $first_event_start = 0;
my $last_event_number = 0;
my $return_code       = 0;
my $is_summary        = 0;

while(<STDIN>) {
    print $_;
    if($_ =~ /CORAL RETURN CODE\: (\d+)/) {
	my $time = timelocal(localtime);
	my $run_time = $time-$job_start_time;
	print LOG_OUTPUT "CORAL END:     $time\t\tERROR CODE: $1\n\n";
	print LOG_OUTPUT "RUN TIME:      $run_time sec.\n";
	$return_code = $1;
	$is_summary = 0;
    }
    elsif($is_summary != 0) {
	print LOG_OUTPUT $_;
    }
    elsif($_ =~ /Total initialization time/) {
	my $time = timelocal(localtime);
	print LOG_OUTPUT "CORAL START:   $time\n\n";
	my $job_start_time = $time;
    }
    elsif($_ =~ /DATABASE INIT started/) {
	my $time = timelocal(localtime);
	print LOG_OUTPUT "DB_INIT START: $time\n\n";
    }
    elsif($_ =~ /DATABASE INIT ended/) {
	my $time = timelocal(localtime);
	print LOG_OUTPUT "DB_INIT END:   $time\n\n";
    }
    elsif($_ =~ /DATABASE SCAN started/) {
	my $time = timelocal(localtime);
	print LOG_OUTPUT "DB_SCAN START: $time\n\n";
    }
    elsif($_ =~ /DATABASE SCAN ended/) {
	my $time = timelocal(localtime);
	print LOG_OUTPUT "DB_SCAN END:   $time\n\n";
    }
    elsif($_ =~ /1ST EVENT UPLOAD started/) {
	my $time = timelocal(localtime);
	print LOG_OUTPUT "FIRST_EVENT:   $time\n";
	$event_start_time  = $time;
	$first_event_start = $time;
	$last_event_number = 1;
    }
    elsif($_ =~ /\[\d+\:\d+\:\d+\]\s+[\d\.]+\s+\-\-\s+Run\s\#\s\d+\s+Event\s+\#\s+\d+\s+\((\d+)\)/) {
	my $event_number = $1;
	my $time = timelocal(localtime);
	my $time_per_event = $event_start_time-$time;
	$time_per_event /= ($event_number-$last_event_number);
	my $total_tpe = ($time-$first_event_start)/$event_number;
	print LOG_OUTPUT "EVENT $event_number: $time\t\t$time_per_event\t\t$total_tpe\n";
	$event_start_time = $time;
	$last_event_number = $event_number;
    }
    elsif($_ =~ /Total time in CORAL/) {
	$is_summary = 1;
    }
    select((select(LOG_OUTPUT), $| = 1)[0]);
    $| = 1;
}

exit($return_code);
