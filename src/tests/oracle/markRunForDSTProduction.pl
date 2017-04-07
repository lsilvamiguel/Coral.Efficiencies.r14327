#!/usr/bin/perl -w

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

if($#ARGV != 0) {
    die("Usage: markRunForDSTProduction.pl <run's list file name>\n");
}

my $file_name = $ARGV[0];
my $debug_flag   = 0; 

if(open(RUNS_LIST,$file_name)) {
    my $dbh;
    my $period_name = "";
    while(<RUNS_LIST>) {
	unless($_ =~ /\w/) {
	    next;
	}
	if($_ =~ /^\#/) { # it's commentary
	    next;
	}
	if($_ =~ /PERIOD\s+(\w+)/i) {
	    if($period_name ne "") {
		&db_close($dbh);
	    }
	    $period_name = $1;
	    my @my_period = GetPeriodRuns($period_name);
	    if($#my_period == 6) {
		my $dst_period  = $my_period[0];
		my $first_run   = $my_period[1];
		my $last_run    = $my_period[2];
		my $database    = $my_period[3];
		my $user_name   = $my_period[4];
		my $password    = $my_period[5];
		my $year        = $my_period[6];
		$dbh = &db_open($user_name,$password,$database);
	    }
	    else {
		die "Incorrect period's name: $period_name\n";
	    }
	    next;
	}
	if($_ =~ /^\s*(\d+)\s*/) { # Run number
	    if($period_name eq "") {
		die "Period is not defined.\n";
	    }
	    my $run_number = $1;
	    print("Set status=1 for the run $run_number.\n");
	    my $sql_query = qq { update runs set status=1 where run_number=$run_number };
	    my $sth = $dbh->prepare($sql_query);
	    $sth->execute();
	    $sth->finish();
	    $dbh->commit;
	}
    }
    if(defined $dbh) {
	&db_close($dbh);
    }
    close(RUNS_LIST);
}
else {
    die "Cannot open file: $file_name\n";
}

#################################################################################################

###########========================== Subroutines: ===============================###############

sub GetPeriodRuns {

    my $dbh_all = &db_open("compass_all","compass","compdb6");
    my $period = $_[0];

    my @ret_period = ();
    
    my $sql_query = qq { select min(RUN_NUMBER), max(RUN_NUMBER) from runs_all_mv where PERIOD='$period' };
    my $sth = $dbh_all->prepare($sql_query);
    $sth->execute();
    
    my $min_run_num;
    my $max_run_num;
    $sth->bind_columns( undef, \$min_run_num, \$max_run_num );
    if( $sth->fetch() ) {
	$debug_flag && print("Period: $period, Min run number: $min_run_num, Max run number: $max_run_num\n");
	$sth->finish();
	$sql_query = qq {select unique DB_NAME, USER_NAME from runs_all_mv where PERIOD='$period' };
	$sth = $dbh_all->prepare($sql_query);
	$sth->execute();
	my $db_name;
	my $user_name;
	$sth->bind_columns( undef, \$db_name, \$user_name );
	if( $sth->fetch() ) {
	    $debug_flag && print("Database: $db_name, Username: $user_name\n");
	    $period =~ /^(\d\d)(\S+)/;
	    my $year = "20$1";
	    $period  = $2;
	    $user_name =~ s/ANL/DST/;
	    push(@ret_period,$period);
	    push(@ret_period,$min_run_num);
	    push(@ret_period,$max_run_num);
	    push(@ret_period,$db_name);
	    push(@ret_period,$user_name);
	    print "Password: ";
	    my $password = <STDIN>;
	    chomp $password;
	    push(@ret_period,$password);
	    push(@ret_period,$year);
	}
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
#			      LongTruncOk => 1,
#			      ora_auto_lob => 1,
#			      LongReadLen => 20000,			      
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
