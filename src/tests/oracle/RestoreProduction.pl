#!/usr/bin/perl -w

require 5.004;

use strict;
use diagnostics;
use DBI;
use DBD::Oracle qw(:ora_types);
use File::stat;
use Time::Local;

$ENV{ORACLE_CERN}     = "/afs/cern.ch/project/oracle";
$ENV{ORACLE_MOUNT}    = "$ENV{ORACLE_CERN}/\@sys";
$ENV{TNS_ADMIN}       = "$ENV{ORACLE_CERN}/admin";
$ENV{ORACLE_HOME}     = "$ENV{ORACLE_MOUNT}/8174";
$ENV{ORA_NLS33}       = "$ENV{ORACLE_HOME}/ocommon/nls/admin/data";
$ENV{DIR_ORACLE}      = "$ENV{ORACLE_HOME}";
$ENV{NLS_LANG}        = "american_america.WE8ISO8859P9";
$ENV{NLS_DATE_FORMAT} = "DD-MON-FXYYYY";
$ENV{STAGE_HOST}      = "stagecompass";
$ENV{PATH}            = "$ENV{DIR_ORACLE}/bin:$ENV{PATH}";

unless( defined($ENV{LD_LIBRARY_PATH}) ) {
    $ENV{LD_LIBRARY_PATH} = "$ENV{DIR_ORACLE}/lib";
}
else {
    $ENV{LD_LIBRARY_PATH} = "$ENV{DIR_ORACLE}/lib:$ENV{LD_LIBRARY_PATH}";
}

if($#ARGV < 1) {
    die "USAGE: RestoreProduction <run number> <slot number> [ALL]";
}

my $run_number  = $ARGV[0];
my $slot_number = $ARGV[1];
my $answer      = $ARGV[2];
unless(defined $answer) {
    $answer = "no";
}
my $user_name   = `whoami`;
chomp $user_name;
my $period_name = "";

my @chunk_list      = ();
my @time_list       = ();
my @nevent_list     = ();
my @cur_speed_list  = ();
my @prod_speed_list = ();
my @status_list     = ();
my @logfile_list    = ();

my $debug_flag  = 0;

my $run = $run_number;
my $slot = $slot_number;

sub db_open;
sub db_close;
sub getInfo;

&main;

exit 0;

#################################################################################################

sub main {
    my $dbh = &db_open("compass_all","compass","compdb6");
    my $sql_query = qq { select file_name, period from all_files where run_number=$run and file_type='RAW' order by file_name };
    my $sth = $dbh->prepare($sql_query);
    $sth->execute();
    my $file_name = "";
    $sth->bind_columns( undef, \$file_name, \$period_name );
    while( $sth->fetch() ) {
	if($file_name =~ /cdr(\d+)\-(\d+)\.raw/) {
	    push(@chunk_list,"$2-$1");
	}
    }
    $sth->finish();
    
    &db_close($dbh);

    foreach my $chunk (@chunk_list) {
	my $log_file = "/afs/cern.ch/compass/scratch/d17/objsrvvy/production/logs/$chunk-$slot.log";
	my @ls = `ls -l $log_file 2> /dev/null`;
	if(defined $ls[0] && $ls[0] =~ /$log_file$/ ) {
	    my @log_info = getInfo($log_file);
	    push(@nevent_list,$log_info[0]);
	    push(@cur_speed_list,$log_info[1]);
	    push(@prod_speed_list,$log_info[2]);
	    push(@status_list,$log_info[3]);
	    push(@time_list,$log_info[4]);
	    push(@logfile_list,$log_file);
	}
	else {
	    push(@time_list,"0");
	    push(@nevent_list,"\&nbsp;");
	    push(@cur_speed_list,"\&nbsp;");
	    push(@prod_speed_list,"\&nbsp;");
	    push(@status_list,"not available");
	    push(@logfile_list,"\&nbsp;");
	}
    }
    my $cur_time = time();
    my @stuck_chunks = ();
    my $is_killed = 0;
    foreach my $chunk (@chunk_list) {
	my $ctime      = shift(@time_list);
	my $nevent     = shift(@nevent_list);
	my $cur_speed  = shift(@cur_speed_list);
	my $prod_speed = shift(@prod_speed_list);
	my $status     = shift(@status_list);
	my $logref     = shift(@logfile_list);
	$chunk =~ /(\d+)\-(\d+)/;
	my $chunk_name = "$2-$1";
	if($cur_time-$ctime > 3600 && !($status =~ /CORAL stopped/ || $status =~ /Exited with error code/ || $status =~ /undefined/)) {
	    print "Production for the chunk $chunk seems to be stuck, do you wish kill the job? (y/n)\n";
	    unless($answer =~ /a/i) { 
		$answer = <STDIN>;
	    }
	    else {
		print "ALL\n";
	    }
	    if($answer =~ /^y/i || $answer =~ /^a/i) {
		my $job_info = `/usr/local/lsf/bin/bjobs -u $user_name -w -J $chunk_name | grep $user_name`;
		if(defined $job_info) {
		    print "$job_info";
		    if($job_info =~ /(\d+)\s+$user_name\s+(\w+)\s+\S+\s+\S+\s+\S+\s+$chunk_name/) {
			my $job_id = $1;
			my $stat   = $2;
			if($stat =~ /PEND/) {
			    print "The job just pended.\n";
			}
			else {
			    my $ret_value = system("/usr/local/lsf/bin/bkill $1");
			    if($ret_value == 0) {
				$is_killed = 1;
				print "Job $chunk is killed\n";
			    }
			    else {
				print "Cannot kill job $chunk_name\n";
			    }
			}
		    }
		    else {
			print "Cannot get information about job.\n";
		    }
		}
		else {
		    print "Sorry, cannot get information about the job.\n";
		}
	    }
	}
    }
    print "Do you wish restart production for the run $run_number in slot $slot_number? (y/n)\n";
    unless($answer =~ /^a/i) {
	my $answer = <STDIN>;
    }
    else {
	print "ALL\n";
    }
    if($answer =~ /^y/i || $answer =~ /^a/i) {
	if($is_killed) {
	    print "Wait 5 minutes, please, until all job will terminated.\n";
	    sleep(300);
	}
	`./prodPeriod $run_number 1 $slot_number rich=no`;
    }
}

#################################################################################################

###########========================== Subroutines: ===============================###############

sub db_open {
    my($user,$passwd,$database)=@_;
    my($tmp_dbh);
# this automatically loads DBD::Oracle
    $tmp_dbh=DBI->connect("dbi:Oracle:".$database,
			  $user,$passwd,{
			      PrintError => 1,
			      AutoCommit => 0,
			      RaiseError => 1
			      }
			  )
	or die "Can't connect to $database";
return $tmp_dbh;
}

sub db_close {
    my($dbh)=@_;
# as Oracle automatically automatically commits any outstanding changes
# at disconnect time, we rollback any uncommited transaction
    $dbh->rollback;
    $dbh->disconnect;
}

sub getInfo {
    my $log_file = $_[0];
    my @ret_values = ();
    my $nevent = 0;
    my $cur_speed = 0;
    my $prod_speed = 0;
    my $status = "undefined";
    my $time = "\&nbsp;";
    if(open(LOG_FILE,$log_file)) {
	while(<LOG_FILE>) {
	    if($_ =~ /CORAL START\:\s+(\d+)/) {
		my $event_time = $1;
		$status = "CORAL started";
		$time = $event_time;
	    }
	    elsif($_ =~ /DB_SCAN START\:\s+(\d+)/) {
		my $event_time = $1;
		$status = "DB scan started";
		$time = $event_time;
	    }
	    elsif($_ =~ /DB_SCAN END\:\s+(\d+)/) {
		my $event_time = $1;
		$status = "DB scan finished";
		$time = $event_time;
	    }
	    elsif($_ =~ /FIRST_EVENT\:\s+(\d+)/) {
		my $event_time = $1;
		$status = "First event read";
		$time = $event_time;
	    }
	    elsif($_ =~ /EVENT\s+(\d+)\:\s+(\d+)\s+([\d\.]+)\s+([\d\.]+)/) {
		my $event_time = $2;
		$nevent = $1;
		$status = "In production";
		$time = $event_time;
		$cur_speed = $3;
		$prod_speed = $4;
	    }
	    elsif($_ =~ /Number of Processed Events\s+\:\s+(\d+)/) {
		my $event_time = stat($log_file)->ctime;
		$status = "Production finished";
		$nevent = $1;
		$time = $event_time;
	    }
	    elsif($_ =~ /CORAL END\:\s+(\d+)\s+ERROR CODE\:\s+(\d+)/) {
		my $event_time = $1;
		if($2 == 0) {
		    $status = "CORAL stopped: OK";
		}
		else {
		    $status = "Exited with error code: $2";
		}
		$time = $event_time;
	    }
	}
	close LOG_FILE;
    }
    push(@ret_values,$nevent);
    push(@ret_values,$cur_speed);
    push(@ret_values,$prod_speed);
    push(@ret_values,$status);
    push(@ret_values,$time);
    return @ret_values;
}
