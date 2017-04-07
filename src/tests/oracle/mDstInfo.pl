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

unless( defined ($ENV{ROOTSYS}) ) {
    $ENV{ROOTSYS} = "/afs/cern.ch/sw/root/v3.05.07/rh73_gcc2952/root";
}

unless( defined($ENV{LD_LIBRARY_PATH}) ) {
    $ENV{LD_LIBRARY_PATH} = "$ENV{DIR_ORACLE}/lib";
}
else {
    $ENV{LD_LIBRARY_PATH} = "$ENV{DIR_ORACLE}/lib:$ENV{LD_LIBRARY_PATH}";
}

unless ($ENV{LD_LIBRARY_PATH} =~ /$ENV{ROOTSYS}\/root/) {
    $ENV{LD_LIBRARY_PATH} = "$ENV{ROOTSYS}/lib:$ENV{LD_LIBRARY_PATH}";
}

if($#ARGV != 0) {
    die("Usage: mDstInfo.pl <run info file>\n");
}

if(!(defined $ENV{MERGING_DIR}) || $ENV{MERGING_DIR} eq "") {
    $ENV{MERGING_DIR} = "/afs/cern.ch/compass/scratch/mdst1";
}

if(!(defined $ENV{PHAST_VERSION}) || $ENV{PHAST_VERSION} eq "") {
    $ENV{PHAST_VERSION} = 6;
}

if(!(defined $ENV{CUR_DIR}) || $ENV{CUR_DIR} eq "") {
    $ENV{CUR_DIR} = "/afs/cern.ch/compass/scratch/d17/objsrvvy/2003/P1J/slot1";
}

my $mdst_info_file = $ARGV[0];
my $run_number     = $ENV{DST_RUN_NUMBER};
my $slot_number    = $ENV{DST_SLOT_NUMBER};
my $work_dir       = $ENV{CUR_DIR};
my $phast_version  = $ENV{PHAST_VERSION};
my $phast          = "$work_dir/phast.exe";
my $scratch_dir    = $ENV{MERGING_DIR};
my $merging_file   = "$scratch_dir/mDST-$run_number-$slot_number-$phast_version.root";
my $return_code    = 0;

my $dbh;

sub GetPeriod;
sub db_open;
sub db_close;
sub mdst_alarm;

my $debug_flag = 0;

my@my_period = GetPeriod($run_number);
if($#my_period == 7) {
    my $database = $my_period[2];
    my $user     = $my_period[3];
    my $password = $my_period[4];
# database handler used for all DBI operations
    $dbh=&db_open($user,$password,$database);
}
else {
    mdst_alarm("Run $run_number is out of DB.\n");
    exit (-1);
}

&main;

&db_close($dbh);

if($return_code == 0) {
    mdst_alarm("Run $run_number DST production is finished.\n");
    unless(open(INFO_FILE, "$mdst_info_file")) {
	mdst_alarm("Cannot find info file: $mdst_info_file\n");
	exit( -2);
    }
    my @file_list = ();
    while(<INFO_FILE>) {
	if($_ =~ /(\S+)\s+(\d)/) {
	    my $file_name = $1;
	    my $is_OK = $2;
	    if($is_OK && ($file_name =~ /mDST-(\d+)-(\d+)-(\d)-\d\.root/)) {
		if($1 == $run_number && $3 == $slot_number) {
		    push(@file_list,$file_name);
		}
	    }
	}
    }
    close(INFO_FILE);
    if($#file_list >= 0) {
	my $list_file = "$scratch_dir/__Run_$run_number\_$slot_number.list";
	if(open(LIST_FILE, ">$list_file")) {
	    foreach my $fn (@file_list) {
		print LIST_FILE "$fn\n";
	    }
	    close(LIST_FILE);
	}
	else {
	    mdst_alarm("Cannot open file: $list_file\n");
	    exit(-3);
	}
	if(system("$phast -m -o $merging_file -l $list_file") != 0) {
	    mdst_alarm("Run: $run_number. Merging error\n");
	    $return_code = 7;
	}
	else {
	    $return_code = 1; # Merging is finished successfully
	}
    }
}

exit $return_code;

#################################################################################################
sub main {
    my $sql_query = qq { select FILE_NAME from FILE_MAPS where FILE_TYPE='RAW' and RUN_NUMBER=$run_number };
    my $sth = $dbh->prepare($sql_query);
    $sth->execute();
    my $raw_file_name = "";
    my %chunk_list;
    $sth->bind_columns( undef, \$raw_file_name);
    while( $sth->fetch() ) {
	$debug_flag && print("$raw_file_name\n");
	$raw_file_name =~ /cdr(\d\d\d\d\d)\-(\d\d\d\d\d)\.raw/;
	if($2 != $run_number) {
	    mdst_alarm "Incorrect RAW file name: $raw_file_name\n";
	    $return_code = 8;
	    return $return_code;
	}
	$chunk_list{$1} = 1;
    }
    if(keys(%chunk_list) == 0) {
	mdst_alarm "No RAW files for the run $run_number.\n";
	$return_code = 6;
	return $return_code;
    }

    $sth->finish();
    $sql_query = qq { select file_name from dst_files where run_number= $run_number
			  and DST_VERSION=$slot_number };
    $sth = $dbh->prepare($sql_query);
    $sth->execute();
    my $dst_file_name;
    $sth->bind_columns( undef, \$dst_file_name);
    while( $sth->fetch() ) {
	$debug_flag && print("$dst_file_name\n");
	$dst_file_name =~ /cdr(\d\d\d\d\d)\-(\d\d\d\d\d)\.dst(\d)/;
	if($2 != $run_number || $3 != $slot_number) {
	    mdst_alarm "Incorrect run or|and slot number in DST file name: $dst_file_name\n";
	    $return_code = 2;
	    return $return_code;
	}
	unless(defined $chunk_list{$1}) {
	    mdst_alarm "DST file $dst_file_name without RAW file.\n";
	    $return_code = 3;
	    return $return_code;
	}
	$chunk_list{$1} += 10;
    }
    unless(open(MDST_INFO, $mdst_info_file)) {
	mdst_alarm "Cannot open mini-DST info file: $mdst_info_file\n";
	$return_code = 4;
	return $return_code;
    }
    my $i = 0;
    my @mdst_list = ();
    while(<MDST_INFO>) {
	$i++;
	if($i == 1 && $_ =~ /Total number of chunks:/) {
	    next; # for compatibility with previous version of prodPeriod
	}
	unless($_ =~ /mDST/) {
	    next; # skip commentary etc.
	}
	if($_ =~ /mDST\-(\d\d\d\d\d)\-(\d\d\d\d\d)\-(\d)\-[\S]+\s+(\d)/) {
	    if($1 != $run_number || $3 != $slot_number) {
#		print ("Error in line $i of $mdst_info_file, line skipped.\n");
		next;
	    }
	    if($4 == 0) { # the mini-DST file is not produced due to an error
		next;
	    }
	    if($chunk_list{$2} != 11) {
		mdst_alarm "mini-DST: $_  without DST file\n";
		next;
	    }
	    if($chunk_list{$2} == 111) {
		mdst_alarm "2 mini-DST files for the chunk $1-$2\n";
		next;
	    }
	    $chunk_list{$2} += 100;
	    push(@mdst_list,$_);
	}
    }
    close MDST_INFO;
    foreach my $chunk_num (keys (%chunk_list)) {
	if($chunk_list{$chunk_num} != 111) {
	    $return_code = 5;
	    last;
	}
    }
    return $return_code;
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

sub mdst_alarm {
    my @arg = @_;
    my %arg=qw();

    foreach (@arg){ $arg{$_}=1};
    @arg=keys(%arg);

    exit if ($#arg<0); # nothing to send
    if(-e "$work_dir/mail_address.enable") {
	my $to_address = `cat $work_dir/mail_address.enable`;
	chomp($to_address);
	my $subject = "DST production: @arg";
	open MAIL,"| /usr/sbin/sendmail -t -oi";
	print MAIL "From: DST_production <na58mdst1\@mail.cern.ch>\n";
	print MAIL "To: $to_address \n";
	print MAIL "Subject: $subject \n";
	print MAIL "\n";
	print MAIL "No message body to be sent along with GSM\n";
	close MAIL;
	print "mdst_alarm: mail sent: $subject\n";
    }
    else {
	print "mdst_alarm: mail disabled\n";
    }
}

