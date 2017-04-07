#!/usr/local/bin/perl -w

require 5.004;

use strict;
use diagnostics;
use DBI;
use DBD::Oracle qw(:ora_types);

my $dst_user = `whoami`;
chomp $dst_user;

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

unless( defined ($ENV{ROOTSYS}) ) {
    $ENV{ROOTSYS} = "/afs/cern.ch/sw/root/v3.05.07/rh73_gcc2952/root";
}

unless( defined($ENV{LD_LIBRARY_PATH}) ) {
    $ENV{LD_LIBRARY_PATH} = "$ENV{DIR_ORACLE}/lib:$ENV{ROOTSYS}/lib";
}

unless ($ENV{LD_LIBRARY_PATH} =~ /$ENV{DIR_ORACLE}\/lib/) {
    $ENV{LD_LIBRARY_PATH} = "$ENV{DIR_ORACLE}/lib:$ENV{LD_LIBRARY_PATH}";
}
unless ($ENV{LD_LIBRARY_PATH} =~ /$ENV{ROOTSYS}\/root/) {
    $ENV{LD_LIBRARY_PATH} = "$ENV{ROOTSYS}/lib:$ENV{LD_LIBRARY_PATH}";
}
unless ($ENV{PATH} =~ /$ENV{ROOTSYS}\/root/ ) {
    $ENV{PATH} = "$ENV{PATH}:$ENV{ROOTSYS}/root/bin";
}

unless (defined $ENV{CUR_DIR}) {
    my $cur_dir = `pwd`;
    chomp $cur_dir;
    $ENV{CUR_DIR} = $cur_dir;
}

if($#ARGV != 2) {
    die("Usage: clear_compass009d.pl <run number> <slot number> <PHAST version>\n");
}

my $run_number    = $ARGV[0];
my $slot_number   = $ARGV[1];
my $phast_version = $ARGV[2];

my $debug_flag = 0;
sub GetPeriodRuns;
sub db_open;
sub db_close;
sub GetMdstFileName;
sub GetFileSize;
sub GetTrafdicFile;

&main;

exit 0;

#################################################################################################
sub main {
    my @my_period = GetPeriodRuns($run_number);
    if($#my_period == 4) {
	my $period      = $my_period[0];
	my $database    = $my_period[1];
	my $user_name   = $my_period[2];
	my $password    = $my_period[3];
	my $year        = $my_period[4];
	print "*** Period $period of $year statistic:\n\n";
	print "Year:                       $year\n";
	print "Period:                     $period\n";
	print "Slot:                       $slot_number\n";
	print "DB username:                $user_name\n";
	print "DB password:                $password\n";
	print "Database:                   $database\n";
	my $castor_dir = "/castor/cern.ch/compass/data/$year/oracle_dst/$period/mDST.chunks";
	my $merged_dir = "/castor/cern.ch/compass/data/$year/oracle_dst/$period/mDST";
	my $histos_dir = "/castor/cern.ch/compass/data/$year/oracle_dst/$period/histos/slot$slot_number";
	my $dbh = &db_open($user_name,$password,$database);
	my $sql_query = qq {select FILE_NAME from FILE_MAPS where RUN_NUMBER=$run_number and FILE_TYPE='RAW' };
	my $file_name = "";
	my @chunk_numbers = ();
	my $sth = $dbh->prepare($sql_query);
	$sth->execute();
	$sth->bind_columns( undef, \$file_name);
	while( $sth->fetch() ) {
	    if($file_name =~ /^cdr(\d\d\d\d\d)\-(\d\d\d\d\d)\.raw/) {
		my $chunk = $1;
		push(@chunk_numbers,$chunk);
	    }
	}
	$sth->finish();
	my @dst_chunk_size = ();
	my @mdst_chunk_size_in_compass009d = ();
	my @mdst_chunk_size_in_castor = ();
	my @mdst_file_name_in_compass = ();
	my @trafdic_file_name = ();
	my @trafdic_file_size = ();
	foreach my $chunk (@chunk_numbers) {
	    $sql_query = qq { select FILE_SIZE from DST_FILES where FILE_NAME like '\%$chunk-$run_number\%' and DST_VERSION = $slot_number };
	    my $file_size = 0;
	    my $sth = $dbh->prepare($sql_query);
	    $sth->execute();
	    $sth->bind_columns( undef, \$file_size);
	    unless($sth->fetch()) { $file_size = 0; }
	    push(@dst_chunk_size,$file_size);
	    $sth->finish();
	    my $mdst_file_name_compass = GetMdstFileName($run_number,$chunk,$slot_number,$phast_version);
	    my $mdst_file_size_compass = 0;
	    if($mdst_file_name_compass ne "") {
		$mdst_file_size_compass = GetFileSize($mdst_file_name_compass);
	    }
	    push(@mdst_chunk_size_in_compass009d,$mdst_file_size_compass);
	    push(@mdst_file_name_in_compass,$mdst_file_name_compass);

	    my $mdst_file_name_castor = "$castor_dir/mDST-$run_number-$chunk-$slot_number-$phast_version.root";
	    my $mdst_file_size_castor = GetFileSize($mdst_file_name_castor);


	    if($mdst_file_size_compass != $mdst_file_size_castor && $mdst_file_size_compass != 0) {
		print "Copy  to castor\n";
		if(system("/usr/local/bin/rfcp $mdst_file_name_compass $mdst_file_name_castor") == 0) {
		    print "Done!\n";
		    $mdst_file_size_castor = GetFileSize($mdst_file_name_castor);
		}
	    }
	    
	    push(@mdst_chunk_size_in_castor,$mdst_file_size_castor);

	    my $trafdic_file_name = GetTrafdicFile($run_number,$chunk,$slot_number);
	    my $trafdic_file_size = 0;
	    if($trafdic_file_name ne "") {
		$trafdic_file_size = GetFileSize($trafdic_file_name);
	    }
	    push(@trafdic_file_name,$trafdic_file_name);
	    push(@trafdic_file_size,$trafdic_file_size);
	}
	&db_close($dbh);
	if($debug_flag) {
	    for(my $i = 0; $i <= $#chunk_numbers; $i++) {
		print "$chunk_numbers[$i]  $dst_chunk_size[$i]\t\t$mdst_chunk_size_in_compass009d[$i]\t\t$mdst_chunk_size_in_castor[$i]\t\t$trafdic_file_name[$i]\n";
	    }
	}

	my $is_good_run = 1;
	my $is_root_files_are_ready_for_merging = 1;
	for(my $i = 0; $i <= $#chunk_numbers; $i++) {
	    if($mdst_chunk_size_in_compass009d[$i] != $mdst_chunk_size_in_castor[$i] && $mdst_chunk_size_in_compass009d[$i] != 0) {
		$is_good_run = 0;
		$is_root_files_are_ready_for_merging = 0;
	    }
	    elsif($trafdic_file_size[$i] == 0) {
		$is_root_files_are_ready_for_merging = 0;
	    }
	}
	if($is_good_run) {
	    my $mdst_file_size = `/usr/local/bin/rfdir $merged_dir/mDST-$run_number-$slot_number-$phast_version.root |  awk \'{ print \$5 }\'`;
	    unless(defined $mdst_file_size) { $mdst_file_size = 0; }
	    else { chomp $mdst_file_size; }
	    if($mdst_file_size != 0) {
		foreach my $ch (@mdst_file_name_in_compass) {
		    my $command = "/usr/local/bin/rfrm $ch";
		    print "$command\n";
		    `$command`;
		}
	    }
	    else {
		print "Mini-DSTs are not merged yet.\n";
	    }
	}
	if($is_root_files_are_ready_for_merging) {
	    print "TRAFDIC files are ready for merging\n";
	    my @cp_fn = ();
	    my $is_OK = 1;
	    foreach my $f (@trafdic_file_name) {
		if(system("rfcp $f ./") == 0) {
		    $f =~ /($run_number\S+\.root)$/;
		    push(@cp_fn, "$1");
		}
		else {
		    print("Cannot copy file $f to local directory.\n");
		    $is_OK = 0;
		}
	    }
	    if($is_OK) {
		$debug_flag && print("$ENV{CUR_DIR}/HistAdder histsum-$run_number.root @cp_fn\n");
		if(system("$ENV{CUR_DIR}/HistAdder histsum-$run_number-$slot_number.root @cp_fn") == 0) {
		    print("TRAFDIC files are merged successfully.\n");
		    `rm @cp_fn`;
		    if(system("/usr/local/bin/rfcp histsum-$run_number-$slot_number.root $histos_dir") == 0) {
			print "Merged TRAFDIC file is copied on the castor successfuly.\n";
			`rm histsum-$run_number-$slot_number.root`;
			print "Remove TRAFDIC files from compass009d.\n";
			$debug_flag && print "rfrm @trafdic_file_name\n";
			if(system("rfrm @trafdic_file_name") == 0) {
			    print "Done!\n";
			}
			else {
			    print "Cannot remove TRAFDIC files from compass009d\n";
			}
		    }
		}
	    }
	}
    }
    else {
	die "Run number $run_number is not found in database.\n"
    }
}

#################################################################################################

###########========================== Subroutines: ===============================###############

sub GetPeriodRuns {

    my $dbh_all = &db_open("compass_all","compass","compdb6");
    my $run = $_[0];

    my @ret_period = ();
    
    my $sql_query = qq { select PERIOD from runs_all_mv where run_number=$run };
    my $sth = $dbh_all->prepare($sql_query);
    $sth->execute();
    my $period = "";

    $sth->bind_columns( undef, \$period );
    if( $sth->fetch() ) {
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

sub GetMdstFileName {
    my ($run_number,$chunk,$slot_number,$phast_version) = @_;
    for(my $i = 1; $i <= 5; $i++) {
	my $file_name = "/shift/compass009d/data0$i/production/mDST-$run_number-$chunk-$slot_number-$phast_version.root";
	my $ret_val = system("/usr/local/bin/rfstat $file_name 2> /dev/null");
	if($ret_val == 0) {
	    return $file_name;
	}
    }
    return "";
}

sub GetFileSize {
    my ($file_name) = @_;
    my $mdst_file_size = `/usr/local/bin/rfdir $file_name | awk \'{ print \$5 }\'`;
    chomp $mdst_file_size;
    unless(defined $mdst_file_size) { $mdst_file_size = 0; }
    return $mdst_file_size;
}

sub GetTrafdicFile {
    my ($run_number,$chunk,$slot_number) = @_;
    my $file_name = "/shift/compass009d/data01/production/TRAFDIC/$run_number-$chunk-$slot_number.root";
    my $ret_val = system("/usr/local/bin/rfstat $file_name 2> /dev/null");
    if($ret_val == 0) {
	return $file_name;
    }
    return "";
}
