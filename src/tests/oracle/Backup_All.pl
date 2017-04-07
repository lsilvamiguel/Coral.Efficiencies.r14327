#!/usr/local/bin/perl -w
#BSUB -L /bin/sh
#BSUB -q prod400
#BSUB -R "type=LINUX7 select[swp>400&&mem>200&&maxmem>500&&cpuf>1.5] order[swp:mem] rusage[swp=400:mem=200]"
#BSUB -J BACKUP_ALL
#BSUB -o /afs/cern.ch/compass/scratch/d17/objsrvvy/2003/P1H/slot1/BACKUP_ALL.log

require 5.004;

use strict;
use diagnostics;
use DBI;

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

##########################################
#Following variables have to be redefined for each new production:
$ENV{CUR_DIR}         = "/afs/cern.ch/compass/scratch/d17/objsrvvy/2003/P1J/slot1";
my $production_year   = 2003;
my $production_period = "P1J";
my $slot_number       = 1;
my $phast_version     = 6;

unless( defined ($ENV{ROOTSYS}) ) {
    $ENV{ROOTSYS} = "/afs/cern.ch/sw/root/v3.05.07/rh73_gcc2952/root";
}
##########################################

#my @compass009d_subdirs = ( "data01", "data02", "data03", "data04", "data05" ); 
#srand;
#my $subdir = rand(@compass009d_subdirs);
#my $mdst_chunks_dir   = "/shift/compass009d/$subdir/production";

my $scratch_dir       = "/afs/cern.ch/compass/scratch/mdst1";
my $username          = `whoami`;
chomp $username;

unless( defined ($ENV{CORAL}) ) {
    $ENV{CORAL} = "$ENV{CUR_DIR}/coral";
}

unless( defined($ENV{LD_LIBRARY_PATH}) ) {
    $ENV{LD_LIBRARY_PATH} = "$ENV{DIR_ORACLE}/lib:$ENV{CORAL}/lib/Linux:$ENV{ROOTSYS}/lib";
}

unless ($ENV{LD_LIBRARY_PATH} =~ /$ENV{DIR_ORACLE}\/lib/) {
    $ENV{LD_LIBRARY_PATH} = "$ENV{DIR_ORACLE}/lib:$ENV{LD_LIBRARY_PATH}";
}
unless ($ENV{LD_LIBRARY_PATH} =~ /$ENV{CORAL}\/lib\/Linux/) {
    $ENV{LD_LIBRARY_PATH} = "$ENV{CORAL}/lib/Linux:$ENV{LD_LIBRARY_PATH}";
}
unless ($ENV{LD_LIBRARY_PATH} =~ /$ENV{ROOTSYS}\/root/) {
    $ENV{LD_LIBRARY_PATH} = "$ENV{ROOTSYS}/lib:$ENV{LD_LIBRARY_PATH}";
}
unless ($ENV{PATH} =~ /$ENV{ROOTSYS}\/root/ ) {
    $ENV{PATH} = "$ENV{PATH}:$ENV{ROOTSYS}/root/bin";
}

sub IsDstProductionFinished;
sub IsMdstMerged;
sub IsHistosMerged;
sub GetMdstFileName;
sub FixBug;
sub BackupLogs;
sub GetPeriod;
sub db_open;
sub db_close;

my @chunk_list      = ();
my @dst_file_list   = ();
my @is_dst_produced = ();

my @run_dirs        = `ls -1 -d $ENV{CUR_DIR}/Run_*`;

foreach my $dir (@run_dirs) {
    chomp $dir;
#    print "$dir\n";
    if($dir =~ /Run_(\d+)/) {
	my $run_number = $1;
	my $ret_value = system("ls $dir/Run_$run_number.backup.tgz 2> /dev/null");
	if($ret_value == 0) {
	    next;
	}
	if(IsDstProductionFinished($run_number,$slot_number) == 1) {
	    if(IsMdstMerged($run_number,$slot_number,$phast_version) == 1) {
		if(IsHistosMerged($run_number,$slot_number) == 1) {
		    BackupLogs($run_number, $slot_number, $production_period, $production_year);
		}
		else {
		    system("$ENV{CUR_DIR}/clear_compass009d.pl $run_number $slot_number $phast_version");
		    if(IsHistosMerged($run_number,$slot_number) == 1) {
			BackupLogs($run_number, $slot_number, $production_period, $production_year);
		    }
		}
	    }
	    else {
		FixBug($run_number, $slot_number, $phast_version);
		if(IsMdstMerged($run_number,$slot_number,$phast_version) == 1) {
		    if(IsHistosMerged($run_number,$slot_number) == 1) {
			BackupLogs($run_number, $slot_number, $production_period, $production_year);
		    }
		    else {
			system("$ENV{CUR_DIR}/clear_compass009d.pl $run_number $slot_number $phast_version");
			if(IsHistosMerged($run_number,$slot_number) == 1) {
			    BackupLogs($run_number, $slot_number, $production_period, $production_year);
			}
		    }
		}
	    }
	}
    }
}


##################################################################################

sub IsDstProductionFinished {
    my $is_production_finished = 1;
    @chunk_list      = ();
    @dst_file_list   = ();
    @is_dst_produced = ();
    my ($run_number, $slot_number) = @_;
    my@my_period = GetPeriod($run_number);
    if($#my_period == 7) {
	my $period_year = $my_period[1];
	my $database    = $my_period[2];
	my $period_name = $my_period[5];
	$period_name =~ /\d\d([\d\w]+)/;
	my $dst_period = $1;
	my $user = "compdst_$period_name";
	my $password = "dstwriter";
# database handler used for all DBI operations
	my $dbh=&db_open($user,$password,$database);

	my $sql_query = qq { select FILE_NAME from FILE_MAPS where RUN_NUMBER=$run_number and FILE_TYPE='RAW' 
			      order by FILE_NAME };
	my $sth = $dbh->prepare( $sql_query );
	$sth->execute();
	my $file_name = "";
	$sth->bind_columns( undef, \$file_name );
	while( $sth->fetch() ) {
	    if($file_name =~ /cdr(\d+)\-(\d+)\.raw/) {
		push(@chunk_list,"$1-$2");
	    }
	}
	$sth->finish();
	foreach my $chunk (@chunk_list) {
	    $sql_query = qq { select FILE_DIR, FILE_NAME from DST_FILES where FILE_NAME like '\%$chunk\%' 
				   and DST_VERSION = $slot_number };
	    $sth = $dbh->prepare( $sql_query );
	    $sth->execute();
	    my $file_dir  = "";
	    my $file_name = "";
	    $sth->bind_columns( undef, \$file_dir, \$file_name );
	    if($sth->fetch()) {
		push(@dst_file_list,"$file_dir$file_name");
		push(@is_dst_produced,1);
	    }
	    else {
		push(@dst_file_list,"");
		push(@is_dst_produced,0);
	    }
	    $sth->finish();
	}
	&db_close($dbh);
	for(my $i = 0; $i <= $#dst_file_list; $i++) {
	    unless($is_dst_produced[$i]) {
		$is_production_finished = 0;
		next;
	    }
	    my $dst_file = $dst_file_list[$i];
	    my $info = `/usr/local/bin/nsls -l $dst_file`;
	    if($info =~ /^\-/ ) { # chunk is not on the tape for the moment (prbably due to CASTOR error)
		print "Chunk $dst_file is not on the tape yet.\n";
		$is_dst_produced[$i] = 0;
		$is_production_finished = 0;
	    }
	}
    }
    if($is_production_finished) {
	print "Production for the run $run_number SLOT$slot_number is finished.\n";
	return 1;
    }
    print "Production for the run $run_number SLOT$slot_number is not finished yet.\n";
    return 0;
}

sub IsMdstMerged {
    my ($run_number,$slot_number,$phast_version) = @_;
    my $castor_dir = "";
    my@my_period = GetPeriod($run_number);
    if($#my_period == 7) {
	my $period_year = $my_period[1];
	my $database    = $my_period[2];
	my $period_name = $my_period[5];
	$period_name =~ /\d\d([\d\w]+)/;
	my $dst_period = $1;
	$castor_dir = "/castor/cern.ch/compass/data/$period_year/oracle_dst/$dst_period";
    }
    my $info = `/usr/local/bin/nsls -l $castor_dir/mDST/mDST-$run_number-$slot_number-$phast_version.root`;
    if($info =~ /^mrw/) {
	print "mDST for the run $run_number SLOT$slot_number has been merged successfully.\n";
	return 1;
    }
    elsif($info =~ /^\-rw/) {
	print "merged mDST file $castor_dir/mDST/mDST-$run_number-$slot_number-$phast_version.root is not on tape yet.\n";
	return 2;
    }
    print "mDST for the run $run_number SLOT$slot_number is not merged yet.\n";
    return 0;
}

sub IsHistosMerged {
    my ($run_number,$slot_number) = @_;
    my $castor_dir = "";
    my@my_period = GetPeriod($run_number);
    if($#my_period == 7) {
	my $period_year = $my_period[1];
	my $database    = $my_period[2];
	my $period_name = $my_period[5];
	$period_name =~ /\d\d([\d\w]+)/;
	my $dst_period = $1;
	$castor_dir = "/castor/cern.ch/compass/data/$period_year/oracle_dst/$dst_period";
    }
    my $info = `/usr/local/bin/nsls -l $castor_dir/histos/slot$slot_number/histsum-$run_number-$slot_number.root`;
    if($info =~ /^mrw/) {
	print "TRAFDIC histos for the run $run_number SLOT$slot_number have been merged successfully.\n";
	return 1;
    }
    elsif($info =~ /^\-rw/) {
	print "merged histos file $castor_dir/histos/slot$slot_number/histsum-$run_number-$slot_number.root is not on tape yet.\n";
	return 2;
    }

    print "TRAFDIC histos for the run $run_number SLOT$slot_number are not merged yet.\n";
    return 0;
}

sub FixBug {
    my ($run_number,$slot_number, $phast_version) = @_;
    my @test_jobs = `/usr/local/lsf/bin/bjobs -w | grep $run_number`;
    if($#test_jobs >= 0) {
	return 0;
    }
    my $castor_dir = "";
    my@my_period = GetPeriod($run_number);
    if($#my_period == 7) {
	my $period_year = $my_period[1];
	my $database    = $my_period[2];
	my $period_name = $my_period[5];
	$period_name =~ /\d\d([\d\w]+)/;
	my $dst_period = $1;
	$castor_dir = "/castor/cern.ch/compass/data/$period_year/oracle_dst/$dst_period";
    }
    $ENV{DST_RUN_NUMBER}  = $run_number;
    $ENV{DST_SLOT_NUMBER} = $slot_number;
    $ENV{PHAST_VERSION}   = $phast_version;
    $ENV{STAGE_POOL}      = "compassmdst";
    $ENV{STAGE_HOST}      = "stagecompass";
    $ENV{MERGING_DIR}     = $scratch_dir;

    my $work_dir  = $ENV{CUR_DIR};
    my $ret_value = system("$work_dir/mDstInfo $work_dir/Run_$run_number/.info");
    print "Return value: $ret_value\n";
    if($ret_value == 256) {
        $ENV{STAGE_POOL} = "compassmdst";
	$ret_value = system("/usr/local/bin/rfcp $scratch_dir/mDST-$run_number-$slot_number-$phast_version.root $castor_dir/mDST/");
	$ENV{STAGE_POOL} = "compass004d";
	if($ret_value == 0) {
	    unlink "$scratch_dir/mDST-$run_number-$slot_number-$phast_version.root";
	    system("cd $scratch_dir; $work_dir/clear_compass009d.pl $run_number $slot_number $phast_version");
	}
    }
    else {
	open INFO, "> $work_dir/Run_$run_number/.info";
	for(my $i = 0; $i <= $#chunk_list; $i++) {
	    my $chunk = $chunk_list[$i];
	    if($chunk =~ /(\d+)-$run_number/) {
		my $chunk_number = $1;
		my $file_name = GetMdstFileName($run_number,$chunk_number,$slot_number,$phast_version);
		my $info = `/usr/local/bin/rfdir $file_name`;
		if(($info =~ /^\S+\s+\S+\s+$username\s+\S+\s+(\d+)/) && $1 != 0) {
		    print INFO "$file_name 1\n";
		}
		else {
		    $info = `/usr/local/bin/nsls -l $castor_dir/mDST.chunks/mDST-$run_number-$chunk_number-$slot_number-$phast_version.root`;
		    if(($info =~ /^\S+\s+\S+\s+$username\s+\S+\s+(\d+)/) && $1 != 0) {
			print INFO "$castor_dir/mDST.chunks/mDST-$run_number-$chunk_number-$slot_number-$phast_version.root 1\n";
		    }
		}
	    }
	}
	close INFO;
	my $ret_value = system("$work_dir/mDstInfo $work_dir/Run_$run_number/.info");
	if($ret_value == 256) {
	    $ENV{STAGE_POOL} = "compassmdst";
	    $ret_value = system("/usr/local/bin/rfcp $scratch_dir/mDST-$run_number-$slot_number-$phast_version.root $castor_dir/mDST/");
	    $ENV{STAGE_POOL} = "compass004d";
	    if($ret_value == 0) {
		unlink "$scratch_dir/mDST-$run_number-$slot_number-$phast_version.root";
		system("cd $scratch_dir; $work_dir/clear_compass009d.pl $run_number $slot_number $phast_version");
	    }
	}
	else {
	    print "Cannot fix bug for the run $run_number!!!\n";
	}
    }
    $ENV{STAGE_POOL}      = "compass004d";
}

sub BackupLogs {
    my ($run_number, $slot_number, $period, $year ) = @_;
    my $short_logs_dir = "/afs/cern.ch/compass/scratch/mdst2/production/logs";

    if(-e "$ENV{CUR_DIR}/Run_$run_number.backup.tgz") {
	print "Backup has been done already.\n";
	return 1;
    }
    `mv $ENV{CUR_DIR}/logFiles/optdst*-$run_number.log $ENV{CUR_DIR}/Run_$run_number/`;
    `mv $short_logs_dir/$year/$period/slot$slot_number/$run_number/*.log $ENV{CUR_DIR}/Run_$run_number/short_logs/`;
    my $ret_value = system("cd $ENV{CUR_DIR}; tar cf - ./Run_$run_number/* | gzip -c > $ENV{CUR_DIR}/Run_$run_number.backup.tgz");
    if($ret_value == 0) {
	`rm $ENV{CUR_DIR}/Run_$run_number/*`;
	`rm $ENV{CUR_DIR}/Run_$run_number/.info`;
	`mv $ENV{CUR_DIR}/Run_$run_number.backup.tgz $ENV{CUR_DIR}/Run_$run_number/`;
	print "Backup of log files has been done successfully.\n";
    }
    else {
	return 0;
    }

    return 1;
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
	#print("Run:      $run_number\nDatabase: $db_name\nPeriod:   $period\nt_min:    $t_min\nt_max:    $t_max\n");
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
