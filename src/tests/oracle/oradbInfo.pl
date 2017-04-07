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
    die("Usage: oradbInfo.pl <period name> <slot number>     - example: oradbInfo.pl 03P2J 1\n");
}

my $period_name = $ARGV[0];
my $slot_number = $ARGV[1];

my $debug_flag = 0;
sub GetPeriodRuns;
sub db_open;
sub db_close;

&main;

exit 0;

#################################################################################################
sub main {
    my @my_period = GetPeriodRuns($period_name);
    if($#my_period == 6) {
	my $period      = $my_period[0];
	my $first_run   = $my_period[1];
	my $last_run    = $my_period[2];
	my $database    = $my_period[3];
	my $user_name   = $my_period[4];
	my $password    = $my_period[5];
	my $year        = $my_period[6];
	print "*** Period $period of $year statistic:\n\n";
	print "Year:                       $year\n";
	print "Period:                     $period\n";
	print "Slot:                       $slot_number\n";
	print "DB username:                $user_name\n";
	print "DB password:                $password\n";
	print "Database:                   $database\n\n\n";
	print "Run \#\tStatus\tN.Ch\tN.DST\ttodo\tN.mDST\tmerged\tPhast ver.\n\n";
	
	my @merged_mdst_list = `/usr/local/bin/nsls /castor/cern.ch/compass/data/$year/oracle_dst/$period/mDST/`;
	my @mdst_chunk_list = `/usr/local/bin/nsls /castor/cern.ch/compass/data/$year/oracle_dst/$period/mDST.chunks/`;
	my $dbh = &db_open($user_name,$password,$database);
	
	my $sql_query = qq { select RUN_NUMBER, STATUS from RUNS order by run_number };
	my @run_numbers = ();
	my @run_status = ();
	my $run_number = 0;
	my $num_runs = 0;
	my $status = 0;
	my $sth = $dbh->prepare($sql_query);
	$sth->execute();
	$sth->bind_columns( undef, \$run_number, \$status);
	while( $sth->fetch() ) {
	    $debug_flag && print("$run_number\n");
	    push(@run_numbers, $run_number);
	    push(@run_status, $status);
	    $num_runs++;
	}
	$sth->finish();
	if($num_runs != 0) {
	    my $i = 0;
	    foreach $run_number (@run_numbers) {
		my $status = shift(@run_status);
		$sql_query = qq {select count(*) from FILE_MAPS where RUN_NUMBER=$run_number and FILE_TYPE='RAW' };
		my $sth = $dbh->prepare($sql_query);
		$sth->execute();
		my $n_chunks = 0;
		$sth->bind_columns( undef, \$n_chunks);
		$sth->fetch();
		$sth->finish();
		if($n_chunks == 0) {
		    next;
		}
		$sql_query = qq {select count(*) from DST_FILES where RUN_NUMBER=$run_number 
				     and DST_VERSION=$slot_number };
		$sth = $dbh->prepare($sql_query);
		$sth->execute();
		my $n_dst_chunks = 0;
		$sth->bind_columns( undef, \$n_dst_chunks);
		$sth->fetch();
		$sth->finish();
		my $todo = $n_chunks-$n_dst_chunks;
		my $n_mdst_chunks = 0;
		foreach my $chunk (@mdst_chunk_list) {
		    if($chunk =~ /mDST\-$run_number\-\d+\-$slot_number\-\d\.root/) {
			$n_mdst_chunks++;
		    }
		}
		my $phast_ver = "";
		my $is_merged = "NO";
		foreach my $m_file (@merged_mdst_list) {
		    if($m_file =~ /mDST\-$run_number\-$slot_number\-(\d)\.root/) {
			$phast_ver = "$1";
			$is_merged = "YES";
		    }
		}
		my $message = "";
		if($todo == 0 && $is_merged eq "YES") {
		    $message = "OK";
		}
		elsif($todo == 0 && $is_merged eq "NO") {
		    $message = " <---- Not merged yet";
		}
		elsif($todo != 0 && $is_merged eq "YES") {
		    $message = " <---- DST production is not finished yet";
		}
		elsif($todo < ($n_chunks/20)) {
		    $message = " <---- Less than 5 percent of the chunks left to be done";
		}
		elsif($n_dst_chunks != $n_mdst_chunks) {
		    $message = " <---- N DST is not equal N mDST";
		}
		$i++;
		if(($i%25) == 0) {
		    print "\n\nRun \#\tStatus\tN.Ch\tN.DST\ttodo\tN.mDST\tmerged\tPhast ver.\n\n";
		}
		printf("%5d\t%5d\t%3d\t%3d\t%4d\t%   3d\t$is_merged\t\t$phast_ver\t$message\n",
		       $run_number,$status,$n_chunks,$n_dst_chunks,$todo,$n_mdst_chunks);
	    }
	}
	else {
	    print "No runs for period $period\n";
	}
	&db_close($dbh);
    }
    else {
	die "Cannot retrieve period information from DB\n";
    }
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
	    push(@ret_period,$period);
	    push(@ret_period,$min_run_num);
	    push(@ret_period,$max_run_num);
	    push(@ret_period,$db_name);
	    push(@ret_period,$user_name);
	    push(@ret_period,"compass");
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

