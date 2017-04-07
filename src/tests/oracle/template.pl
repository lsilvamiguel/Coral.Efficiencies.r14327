#!/usr/local/bin/perl -w
#BSUB -L /bin/sh
#BSUB -q prod400
#BSUB -R "type=LINUX7 select[swp>400&&mem>200&&maxmem>500&&cpuf>1.5] order[swp:mem] rusage[swp=400:mem=200]"

# The queue could be modified by prodPeriod program,
# use its parameter queue=queue_name

require 5.004;

use strict;
use diagnostics;

# Following variables should be defined by prodPeriod program:
### $current_directory
### $run_number
### $chunk_number
### $slot_number
### $dst_period
### $dst_year

unless(defined $current_directory) {
    die "Current directory is not defined.\n";
}
unless(defined $run_number) {
    die "Run number is not defined.\n";
}
unless(defined $chunk_number) {
    die "Chunk number is not defined.\n";
}
unless(defined $slot_number) {
    die "Slot number is not defined.\n";
}
unless(defined $dst_period) {
    die "Period name is not defined.\n";
}
unless(defined $dst_year) {
    die "Year of run is not defined.\n";
}

$ENV{ORACLE_CERN}        = "/afs/cern.ch/project/oracle";
$ENV{ORACLE_MOUNT}       = "$ENV{ORACLE_CERN}/\@sys";
$ENV{TNS_ADMIN}          = "$ENV{ORACLE_CERN}/admin";
$ENV{ORACLE_HOME}        = "$ENV{ORACLE_MOUNT}/920";
$ENV{ORA_NLS33}          = "$ENV{ORACLE_HOME}/ocommon/nls/admin/data";
$ENV{DIR_ORACLE}         = "$ENV{ORACLE_HOME}";
$ENV{NLS_LANG}           = "american_america.WE8ISO8859P9";
$ENV{NLS_DATE_FORMAT}    = "DD-MON-FXYYYY";
$ENV{STAGE_HOST}         = "stagecompass";
$ENV{STAGE_POOL_EXCLUDE} = "/afs/cern.ch/user/o/objsrvvy/castor/.StageHelperPoolsExclude.inf";
$ENV{CORAL}              = "$current_directory/coral";
$ENV{ROOTSYS}            = "/afs/cern.ch/sw/root/v3.05.07/rh73_gcc2952/root";
$ENV{LD_LIBRARY_PATH}    = "$ENV{DIR_ORACLE}/lib:$ENV{CORAL}/lib/Linux:$ENV{ROOTSYS}/lib";
$ENV{PATH}               = "$ENV{DIR_ORACLE}/bin:$ENV{PATH}:$ENV{ROOTSYS}/root/bin";
$ENV{DST_SLOT_NUMBER}    = $slot_number;
$ENV{DST_RUN_NUMBER}     = $run_number;
$ENV{DST_CHUNK_NUMBER}   = $chunk_number;
$ENV{MERGING_DIR}        = "/afs/cern.ch/compass/scratch/mdst1";
my $phast_version        = 6;
$ENV{PHAST_VERSION}      = $phast_version;
$ENV{CUR_DIR}            = $current_directory;

my $mdst_merging_dir     = "/afs/cern.ch/compass/scratch/mdst1";
my $online_log_dir       = "/afs/cern.ch/compass/scratch/mdst2/production/logs/$dst_year/$dst_period/slot$slot_number";
my $control_prog         = "$current_directory/prodControl";
my $opt_file             = "$current_directory/Run_$run_number/$chunk_number-$run_number.opt";
my $mdst_file_name       = "mDST-$run_number-$chunk_number-$slot_number-$phast_version.root";
my $gfile_file_name      = "histodst$slot_number-$chunk_number-$run_number.gfile";
my $gfile_dir            = "/castor/cern.ch/compass/data/$dst_year/oracle_dst/$dst_period/RICH";
my $mdst_info_file       = "$current_directory/Run_$run_number/.info";
my $mdst_eor_prog        = "$current_directory/mDstInfo";
my $castor_dir           = "/castor/cern.ch/compass/data/$dst_year/oracle_dst/$dst_period";
my $production_prog      = "$current_directory/coral.exe";
### setenv PRODUCTION_PROGRAM ${CUR_DIR}/coral.exe   <-- needs for prodPeriod program

my $username             = `whoami`;
chomp $username;

my $mdst_chunks_dir = "";
srand;
my $subdir = int(rand 5) + 1;;
$mdst_chunks_dir   = "/shift/compass009d/data0$subdir/production";

my $xxx = $run_number%100;
$xxx = int $xxx/10;
$xxx *= 10;
if($xxx == 0) {
    $xxx = "00";
}

my $dst_file_name = "$castor_dir/$xxx/slot$slot_number/cdr$chunk_number-$run_number.dst$slot_number";
my $raw_file_name = "/castor/cern.ch/compass/data/$dst_year";
if($dst_year == 2002 || $dst_year == 2001) {
    $raw_file_name = "$raw_file_name/raw_migration/$dst_period/$xxx/cdr$chunk_number-$run_number.raw";
}
else {
    $raw_file_name = "$raw_file_name/raw/$dst_period/cdr$chunk_number-$run_number.raw";
}

my $job_test = `/usr/local/lsf/bin/bjobs -w -u u_vy | grep \"$run_number\\-$chunk_number\" | wc -l`;
chomp $job_test;
if(defined $job_test && $job_test != 0) {
    die "The job $run_number-$chunk_number is running already.\n";
}

print("Production user $username on host $ENV{HOST}\n");
my $test_raw = `/usr/local/bin/stageqry -M $raw_file_name | grep STAGED`;

if(defined $test_raw && $test_raw =~ /\S+\s+STAGED\s+\S+\s+\S+\s+(\S+)/) {
    $ENV{STAGE_POOL} = $1;
}
else {
    $ENV{STAGE_POOL} = `/afs/cern.ch/user/o/objsrvvy/public/stagepool.pl $ENV{STAGE_HOST}`;
    if(defined $ENV{STAGE_POOL}) {
	chomp $ENV{STAGE_POOL};
    }
}

unlink "$current_directory/logFiles/optdst$slot_number-$chunk_number-$run_number.log";

my $execute_command = "(time $production_prog $opt_file ; echo \"CORAL RETURN CODE: \$?\") | $control_prog $online_log_dir";

my $ret_value = system($execute_command);

if($ret_value != 0) {
    unlink <core*>;
    unlink $mdst_file_name;
    unlink "myfile-data.gfile";
    unlink "trafdic.root";
    `/usr/local/bin/rfrm $dst_file_name`;
    `/usr/local/bin/stageclr -h $ENV{STAGE_HOST} -M $dst_file_name`;
    die "Production error: $ret_value\n";
}

$ENV{STAGE_POOL} = `/afs/cern.ch/user/o/objsrvvy/public/stagepool.pl $ENV{STAGE_HOST}`;
chomp $ENV{STAGE_POOL};
print "Current castor staging pool is $ENV{STAGE_POOL}\n";

print "Copy trafdic.root  to $mdst_chunks_dir/TRAFDIC/$run_number-$chunk_number-$slot_number.root\n";
$ret_value = system("/usr/local/bin/rfcp trafdic.root $mdst_chunks_dir/TRAFDIC/$run_number-$chunk_number-$slot_number.root");

if($ret_value != 0) { #in case of error try again after 10 min
    print "Waiting 10 minutes...\n";
    sleep(600);
    $ENV{STAGE_POOL} = `/afs/cern.ch/user/o/objsrvvy/public/stagepool.pl $ENV{STAGE_HOST}`;
    chomp $ENV{STAGE_POOL};
    print "Current castor staging pool is $ENV{STAGE_POOL}\n";
    print "Try again...\n";
    $ret_value = system("/usr/local/bin/rfcp trafdic.root $mdst_chunks_dir/TRAFDIC/$run_number-$chunk_number-$slot_number.root");
}
if($ret_value == 0) {
    print "done!\n";
}
else {
    print "Cannot copy trafdic.root file to $mdst_chunks_dir/TRAFDIC/\n";
}

unlink "trafdic.root";

my $mdst_chunk_location = "";

print "Copy $mdst_file_name to $castor_dir/mDST.chunks/\n";

$ret_value = system("/usr/local/bin/rfcp $mdst_file_name $castor_dir/mDST.chunks/");
if($ret_value != 0) {
    print "Waiting 10 minutes...\n";
    sleep(600);
    $ENV{STAGE_POOL} = `/afs/cern.ch/user/o/objsrvvy/public/stagepool.pl $ENV{STAGE_HOST}`;
    chomp $ENV{STAGE_POOL};
    print "Current castor staging pool is $ENV{STAGE_POOL}\n";
    print "Try again...\n";
    $ret_value = system("/usr/local/bin/rfcp $mdst_file_name $castor_dir/mDST.chunks/");
}
if($ret_value == 0) {
    print "done!\n";
    $mdst_chunk_location = "$castor_dir/mDST.chunks/$mdst_file_name";
}
else {
    print "Cannot copy $mdst_file_name to $castor_dir/mDST.chunks/\n";
    `/usr/local/bin/stageclr -h $ENV{STAGE_HOST} -M $castor_dir/mDST.chunks/$mdst_file_name`;
}

print "Copy $mdst_file_name to $mdst_chunks_dir/\n";
$ret_value = system("/usr/local/bin/rfcp $mdst_file_name $mdst_chunks_dir/");
if($ret_value != 0) { 
    $subdir = int(rand 5) + 1;
    $mdst_chunks_dir   = "/shift/compass009d/data0$subdir/production";
    print "Try again...\n";
    $ret_value = system("/usr/local/bin/rfcp $mdst_file_name $mdst_chunks_dir/");
}
if($ret_value != 0) {
    $subdir = int(rand 5)+1;
    $mdst_chunks_dir   = "/shift/compass009d/data0$subdir/production";
    print "Try again...\n";
    $ret_value = system("/usr/local/bin/rfcp $mdst_file_name $mdst_chunks_dir/");
}
if($ret_value == 0) {
    print "done!\n";
    $mdst_chunk_location = "$mdst_chunks_dir/$mdst_file_name";
}
else {
    print "Cannot copy $mdst_file_name to $mdst_chunks_dir/\n";
}

if($mdst_chunk_location eq "") {
    if(system("cp $mdst_file_name $ENV{MERGING_DIR}/") == 0) {
	$mdst_chunk_location = "$ENV{MERGING_DIR}/$mdst_file_name";
    }
}

print "Mini-DST chunk location: $mdst_chunk_location\n";

if($mdst_chunk_location ne "") {
    my $ret_val = system("echo $mdst_chunk_location 1 >> $mdst_info_file");
    if($ret_val != 0) {
	for(my $i = 0; $i < 10; $i++) {
	    sleep(10);
	    $ret_val = system("echo $mdst_chunk_location 1 >> $mdst_info_file");
	    if($ret_val == 0) {
		last;
	    }
	}
    }
}
else {
    print "ERROR:   rfcp error mini-DST is missing.\n";
}

unlink $mdst_file_name;

$ret_value = system("$mdst_eor_prog $mdst_info_file");

print "Return value = $ret_value\n"; 

if($ret_value == 256) { #return value of mdstInfo program == 1 (0x100 for PERL)
    $ENV{STAGE_POOL} = "compassmdst";
    $ret_value = system("/usr/local/bin/rfcp $mdst_merging_dir/mDST-$run_number-$slot_number-$phast_version.root $castor_dir/mDST/");
    if($ret_value == 0) {
	unlink "$mdst_merging_dir/mDST-$run_number-$slot_number-$phast_version.root";
	`$current_directory/clear_compass009d.pl $run_number $slot_number $phast_version`;
    }
    else {
	print("ERORR: rfcp error. mDST merging file will be left on $mdst_merging_dir/\n");
	`stageclr -h $ENV{STAGE_HOST} -M $castor_dir/mDST/mDST-$run_number-$slot_number-$phast_version.root`;
    }
}

$ENV{STAGE_POOL} = `/afs/cern.ch/user/o/objsrvvy/public/stagepool.pl $ENV{STAGE_HOST}`;
chomp $ENV{STAGE_POOL};

if(-x "myfile-data.gfile") {
    $ret_value = system("/usr/local/bin/rfcp myfile-data.gfile $gfile_dir/$gfile_file_name");
    if($ret_value != 0) {
	print "ERROR: rfcp error. Cannot copy gfile to castor\n";
	`stageclr -h $ENV{STAGE_HOST} -M $gfile_dir/$gfile_file_name`;
    }
    unlink "myfile-data.gfile";
}


