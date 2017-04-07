#!/usr/local/bin/perl -w

require 5.004;

use strict;
use diagnostics;
use DBI;
use DBD::Oracle qw(:ora_types);

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

if($#ARGV != 1) {
    die("Usage: removeTBnames.pl <run number> <slot number>\n");
}

my $run_number  = $ARGV[0];
my $slot_number = $ARGV[1];

my $debug_flag = 0;
#my $debug_flag = 1;

sub GetPeriod;
sub db_open;
sub db_close;
sub RemoveTBNamesString;

&main;
exit 0;

#################################################################################################

sub main {
    my@my_period = GetPeriod($run_number);
    if($#my_period == 7) {
	my $period_year = $my_period[1];
	my $database    = $my_period[2];
	my $period_name = $my_period[5];
	$debug_flag && print("Period $period_name : $period_year, $database\n");
	my $user = "compdst_$period_name";
	my $password = "dstwriter";
# database handler used for all DBI operations
	my $dbh=&db_open($user,$password,$database);
	
	RemoveTBNamesString($run_number,$slot_number,$dbh);
	&db_close($dbh);
    }
    else {
	die("Run number $run_number is out of database\n");
    }
}

#################################################################################################

###########========================== Subroutines: ===============================###############

sub GetPeriod {
    my $dbh_all = &db_open("compass_all","compass","compdb6");
    my $run_number = $_[0];

    my @ret_period = ();
    
    my $sql_query = qq { select DB_NAME, PERIOD, T_MIN, T_MAX from runs_all_mv where RUN_NUMBER=$run_number };
    my $sth = $dbh_all->prepare($sql_query);
    $sth->execute();

    my $db_name;
    my $run_year = "20";
    my $db_user = "companl_";
    my $db_password = "compass";
    my $period;
    my $t_min;
    my $t_max;
    $sth->bind_columns( undef, \$db_name, \$period, \$t_min, \$t_max );
    if( $sth->fetch() ) {
	$debug_flag && print("Run:      $run_number\nDatabase: $db_name\nPeriod:   $period\nt_min:    $t_min\nt_max:    $t_max\n");
	$db_user = "$db_user$period";
	$period =~ /^(\d\d)/;
	my $run_year = "$run_year$1";
	push(@ret_period,$run_number);
	push(@ret_period,$run_year);
	push(@ret_period,$db_name);
	push(@ret_period,$db_user);
	push(@ret_period,$db_password);
	push(@ret_period,$period);
	push(@ret_period,$t_min);
	push(@ret_period,$t_max);
    }
    $sth->finish();
    &db_close($dbh_all);
    return @ret_period;
}

sub db_open {
    my($user,$passwd,$database)=@_;
    my($tmp_dbh);
# this automatically loads DBD::Oracle
    $tmp_dbh=DBI->connect("dbi:Oracle:".$database,
			  $user,$passwd,{
			      PrintError => 1,
			      AutoCommit => 0,
			      LongTruncOk => 1,
			      ora_auto_lob => 1,
			      LongReadLen => 20000,			      
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

sub Remove;

sub RemoveTBNamesString {
    my $run_num  = $_[0];
    my $slot_num = $_[1];
    my $dbh      = $_[2];

    my $sqlx = qq { select LOG_INFO from RUNS where RUN_NUMBER=$run_num };
    my $sthx = $dbh->prepare( $sqlx,, );
    $sthx->execute || die ("Select failed for $sqlx: . $DBI::errstr\n . $!");
    my $sql = qq { update RUNS set LOG_INFO = ? where  RUN_NUMBER=$run_num };
    my $sth = $dbh->prepare($sql,,);
    if(my $pointer = $sthx->fetchrow_hashref) {
	if(defined $pointer->{'LOG_INFO'}) {
	    $debug_flag && print("OLD: $pointer->{'LOG_INFO'}\n");
	    my $log_info = $pointer->{'LOG_INFO'};
	    $pointer->{'LOG_INFO'} = Remove($log_info,$slot_num);
	    $debug_flag && print("NEW: $pointer->{'LOG_INFO'}\n");
	    $sth->execute($pointer->{'LOG_INFO'}) || die ("Select failed for $sql: . $DBI::errstr\n . $!");
	}
    }
    $sth->finish();
    $sthx->finish();
    $dbh->commit;
}

sub Remove {
    my $log_info = $_[0];
    my $slot_num = $_[1];

    my $start_string = ":Start of TBNames string:";
    my $stop_string  = ":End of TBNames string:";

    my $tbn_slot1 = "";
    my $tbn_slot2 = "";
    my $tbn_slot3 = "";
    $debug_flag && print("log_info: $log_info\n\n\nslot_number: $slot_num\n\n\n");
    if( defined $log_info ) {
	if($log_info =~ /^$start_string(\w*)$stop_string/) {
	    $tbn_slot1 = $1;
	    $debug_flag && print("SLOT1: $tbn_slot1\n");
	}
	if($log_info =~ /^$start_string(\w*)$stop_string$start_string(\w*)$stop_string/) {
	    $tbn_slot2 = $2;
	    $debug_flag && print("SLOT2: $tbn_slot2\n");
	}
	if($log_info =~ /^$start_string(\w*)$stop_string$start_string(\w*)$stop_string$start_string(\w*)$stop_string/) {
	    $tbn_slot3 = $3;
	    $debug_flag && print("SLOT3: $tbn_slot3\n");
	}
    }
    $log_info = "";
    if($slot_num == 1) {
	$tbn_slot1 = "";
    }
    elsif($slot_num == 2) {
	$tbn_slot2 = "";
    }
    elsif($slot_num == 3) {
	$tbn_slot3 = "";
    }
    my $is_tbnames = 0;
    if($tbn_slot3 ne "") {
	$log_info = "$start_string$tbn_slot1$stop_string$start_string$tbn_slot2$stop_string$start_string$tbn_slot3$stop_string";
	$is_tbnames = 1;
    }
    elsif($tbn_slot2 ne "") {
	$log_info = "$start_string$tbn_slot1$stop_string$start_string$tbn_slot2$stop_string";
	$is_tbnames = 1;
    }
    elsif($tbn_slot1 ne "") {
	$log_info = "$start_string$tbn_slot1$stop_string";
	$is_tbnames = 1;
    }
    if($is_tbnames == 0) {
	$log_info = "";
    }
    $debug_flag && print("New TBNames string: $log_info\n");
    return $log_info;
 }
