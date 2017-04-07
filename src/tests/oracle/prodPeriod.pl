#!/usr/local/bin/perl -w

require 5.004;

use strict;
use diagnostics;
use DBI;

sub Usage;

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
$ENV{CUR_DIR}         = $ENV{PWD};

unless( defined ($ENV{CORAL}) ) {
    $ENV{CORAL} = "$ENV{PWD}/coral";
}

unless( defined ($ENV{ROOTSYS}) ) {
    $ENV{ROOTSYS} = "/afs/cern.ch/sw/root/v3.05.07/rh73_gcc2952/root";
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

if($#ARGV < 2) {
    Usage();
    exit;
}

my $first_run_number = $ARGV[0];
my $n_runs           = $ARGV[1];
my $slot_number      = $ARGV[2];
my $read_parameter   = "read=raw";   # by default read RAW file
my $test_parameter   = "test=no";    # by default no test of program (submit job)
my $mdst_parameter   = "mdst=yes";   # produce mDST by default
my $rich_parameter   = "rich=rand";  # random RICH gfile production by default
my $dst_parameter    = "dst=yes";    # produce DST by default
my $send_to_top      = 0;            # move job to top of queue (default value is NO)
my $prod_queue       = "";           # production queue
my $producer_mail    = "Jiawei.Zhao\@cern.ch";

if($#ARGV > 2) {
    my $i = $#ARGV;
    while($i > 2) {
	my $parameter = $ARGV[$i];
	if($parameter =~ /read\=/) {
	    $read_parameter = $parameter;
	}
	elsif($parameter =~ /queue\=(\S+)/) {
	    $prod_queue = $1;
	}
	elsif($parameter =~ /top/i) {
	    $send_to_top = 1;   
	}
	elsif($parameter =~ /test\=/) {
	    $test_parameter = $parameter;
	}
	elsif($parameter =~ /mdst\=/) {
	    $mdst_parameter = $parameter;
	}
	elsif($parameter =~ /rich\=/) {
	    $rich_parameter = $parameter;
	}
	elsif($parameter =~ /dst\=/) {
	    $dst_parameter = $parameter;
	}
	elsif($parameter =~ /mail\=(\S+)/) {
	    $producer_mail = $1;
	}
	$i--;
    }
}
else {
    Usage();
    exit;
}

sub GetPeriod;
sub db_open;
sub db_close;
sub TestTBNamesString;
sub SetTBNamesString;

my $debug_flag = 0;

if($debug_flag == 1) {
    $test_parameter = "test=yes";
}

&main;
exit 0;

#################################################################################################

sub main {
    my@my_period = GetPeriod($first_run_number);
    if($#my_period == 7) {
	my $period_year = $my_period[1];
	my $database    = $my_period[2];
	my $period_name = $my_period[5];
	$period_name =~ /\d\d([\d\w]+)/;
	my $dst_period = $1;
	$debug_flag && print("Period $period_name : $period_year, $database\n");
	my $user = "compdst_$period_name";
	my $password = "dstwriter";
	
# database handler used for all DBI operations
	my $dbh=&db_open($user,$password,$database);
	
# Get list of runs from DB:
############ Last version: (12/11/2003)
	my $sql_query = qq { select RUN_NUMBER from RUNS where STATUS=1 and RUN_NUMBER>=$first_run_number order by RUN_NUMBER }; #select only runs marked as GOOD for DST production
	my $sth = $dbh->prepare( $sql_query );
	$sth->execute();
	
	my( $run_number);
	$sth->bind_columns( undef, \$run_number);
	my @run_list = ();
	while( $sth->fetch() && $#run_list < ($n_runs-1) ) {
	    $debug_flag && print("$run_number\n");
	    push(@run_list,$run_number);
	}
	$sth->finish();
	$debug_flag && print("@run_list\n");
	if($#run_list < 0) {
	    die "No any runs for that period.\n"
	    }

	my @prod_file_list = ();

	foreach my $rn (@run_list) {
	    if(TestTBNamesStrings($rn,$dbh,$period_year,$slot_number) == 0) {
		$debug_flag && print("No TBName string for the slot $slot_number\n");
		if(SetTBNamesString($rn,$slot_number,$period_year,$period_name) != 0) { 
		    &db_close($dbh); die "Cannot set TBName string.\n";
		}
	    }
	    $sql_query = qq { select FILE_DIR, FILE_NAME from FILE_MAPS where RUN_NUMBER=$rn and FILE_TYPE='RAW' 
			      order by FILE_NAME };
	    $sth = $dbh->prepare( $sql_query );
	    $sth->execute();
	    my ($file_dir,$file_name);
	    $sth->bind_columns( undef, \$file_dir, \$file_name );
	    my @raw_file_list = ();
	    while( $sth->fetch() ) {
		my $full_file_name = "$file_dir$file_name";
		$debug_flag && print("$full_file_name\n");
		push(@raw_file_list,$full_file_name);
	    }
	    $sth->finish();
	    my @dst_file_list = ();
	    $sql_query = qq { select FILE_NAME from DST_FILES where RUN_NUMBER=$rn and DST_VERSION=$slot_number 
				  order by FILE_NAME };
	    $sth = $dbh->prepare( $sql_query );
	    $sth->execute();
	    $sth->bind_columns( undef, \$file_name );
	    while( $sth->fetch() ) {
		$debug_flag && print("$file_name\n");
		push(@dst_file_list,$file_name);
	    }
	    $sth->finish();
	    foreach my $fn (@raw_file_list) {
		$fn =~ /\/cdr(\d\d\d\d\d\-\d\d\d\d\d)\.raw/;
		my $chunk_name = $1;
		my $is_produced = 0;
		foreach my $dst_fn (@dst_file_list) {
		    if($dst_fn =~ /cdr$chunk_name\.dst/) {
			$is_produced = 1;
			last;
		    }
		}
		if($is_produced == 0) { # the DST is not produced for the file yet
		    # to check is it in production already:
		    print "Chunk: $chunk_name:\n";
		    my @prod_list = `/usr/local/lsf/bin/bjobs -w -u all -J $chunk_name`;
		    $debug_flag && print("PROD_LIST: @prod_list\n#PROD_LIST: $#prod_list\n");
		    unless($#prod_list < 0) {
			$is_produced = 1;
		    }
		}
		if($is_produced == 0) {
		    push(@prod_file_list,$fn);
		}
	    }
	}
	if($debug_flag) {
	    print "List of RAW files for production:\n";
	    foreach my $fn (@prod_file_list) {
		print("$fn\n");
	    }
	}
	&db_close($dbh);
	unless($#prod_file_list < 0) {
	    my $cur_dir = $ENV{PWD};
	    my $log_dir_name = "$cur_dir/logFiles";
	    unless(-d $log_dir_name) {
		if(system("mkdir $log_dir_name") != 0) {
		    die "Cannot create directory: $log_dir_name\n";
		}
	    }
	    foreach my $fn (@prod_file_list) {
		$fn =~ /cdr(\d\d\d\d\d)\-(\d\d\d\d\d)\.raw$/;
		my $chunk_num = $1;
		my $run_num   = $2;
		$run_num =~/(\d)(\d)$/;
		my $XXX = "$1"."0";
		my $odd_even = 1;
		if($2 =~ /[02468]/) {
		    $odd_even = 2;
		}
		my $rich_gfile_production = 0;
		unless($rich_parameter =~ /rich\=no/) {
		    if($rich_parameter =~ /rich\=yes/) {
			$rich_gfile_production = 1;
		    }
		    elsif(rand(20) < 1.0) { # RICH gfile will be produced for 5% of RAW files
			$rich_gfile_production = 1;
		    }
		}

		my $log_file_name = "$log_dir_name/optdst$slot_number-$chunk_num-$run_num.log";
		if(defined $chunk_num && defined $run_num) {
		    my $opt_dir_name = "$cur_dir/Run_$run_num";
		    unless(-d $opt_dir_name) {
			if(system("mkdir $opt_dir_name") != 0) {
			    die "Cannot create directory: $opt_dir_name\n";
			}
		    }
		    unless(open(TEMPL_OPT, "$ENV{PWD}/template.opt")) {
			die("Cannot open file: $ENV{PWD}/template.opt\n");
		    }
		    unless(open(TEMPL_LSF, "$ENV{PWD}/template.lsf")) {
			die("Cannot open file: $ENV{PWD}/template.lsf\n");
		    }
		    my $opt_file_name = "$opt_dir_name/$chunk_num-$run_num.opt";
		    my $lsf_file_name = "$opt_dir_name/$chunk_num-$run_num.lsf";
		    unless(open(JOB_OPT, "> $opt_file_name")) {
			die("Cannot open file: $opt_file_name\n");
		    }
		    unless(open(JOB_LSF, "> $lsf_file_name")) {
			die("Cannot open file: $lsf_file_name\n");
		    }
		    
		    my $phast_version = 5;
		    while(<TEMPL_LSF>){
			if($_ =~ /\#BSUB\ \-q/ && $prod_queue ne "") {
			    print JOB_LSF "#BSUB -q $prod_queue\n";
			    next;
			}
			if($_ =~ /\#BSUB\ \-R/) {
			    print JOB_LSF $_;
			    print JOB_LSF "#BSUB -o $log_file_name\n";
			    print JOB_LSF "#BSUB -J $chunk_num-$run_num\n\n";
			}
			elsif($_ =~ /\#\#\#\s+\$current_directory/) {
			    print JOB_LSF "my \$current_directory = \"$ENV{PWD}\";\n"
			}
			elsif($_ =~ /\#\#\#\s+\$run_number/) {
			    print JOB_LSF "my \$run_number = $run_num;\n";
			}
			elsif($_ =~ /\#\#\#\s+\$chunk_number/) {
			    print JOB_LSF "my \$chunk_number = \"$chunk_num\";\n";
			}
			elsif($_ =~ /\#\#\#\s+\$slot_number/) {
			    print JOB_LSF "my \$slot_number = $slot_number;\n";
			}
			elsif($_ =~ /\#\#\#\s+\$dst_period/) {
			    print JOB_LSF "my \$dst_period = \"$dst_period\";\n";
			}
			elsif($_ =~ /\#\#\#\s+\$dst_year/) {
			    print JOB_LSF "my \$dst_year = $period_year;\n";
			}
			else {
			    print JOB_LSF $_;
			}
			if($_ =~ /my \$phast_version\s+\=\s+([\S]+)\;/) {
			    $phast_version = $1;
			}
		    }
		    close TEMPL_LSF;
		    close JOB_LSF;

		    my $mdst_file_name = "mDST-$run_num-$chunk_num-$slot_number-$phast_version.root";
		    print JOB_OPT "//DST production with Oracle DB:\n";
		    print JOB_OPT "Data CsStore oracle\n";
		    print JOB_OPT "Data container cdr$chunk_num-$run_num\n";
		    print JOB_OPT "Data year $period_year\n";
		    print JOB_OPT "Data period $dst_period\n";
		    print JOB_OPT "Data run select $run_num\n";
		    print JOB_OPT "Dst select slot $slot_number\n";
		    if($read_parameter =~ /read\=dst/) {
			print JOB_OPT "Data type dst\n";
		    }
		    else {
			print JOB_OPT "Data type raw\n";
		    }
		    unless($mdst_parameter =~ /mdst\=no/) {
			print JOB_OPT "mDST file $mdst_file_name\n";
		    }
		    print JOB_OPT "// ********************\n\n";
		    while(<TEMPL_OPT>) {
			if($rich_gfile_production == 0 && $_ =~ /RICHONE/) {
			    next;
			}
			if($dst_parameter =~ /dst\=no/ && $_ =~ /make\s+DST/) {
			    next;
			}
			print JOB_OPT $_;
		    }
		    close TEMPL_OPT;
		    close JOB_OPT;
		    `touch $opt_dir_name/.info`;
		    unless($test_parameter =~ /test\=yes/) {
			my $ret_string = `/usr/local/lsf/bin/bsub -u $producer_mail < $lsf_file_name`;
			print "$ret_string";
			if($send_to_top) {
			    if($ret_string =~ /\<(\d+)\>/) {
				system("/usr/local/lsf/bin/btop $1");
			    }
			}
		    }
		    else {
			print "/usr/local/lsf/bin/bsub -u $producer_mail < $lsf_file_name\n";
		    }
		}
		else {
		    print("Incorrect RAW file name: $fn . Skipped.\n");
		}
	    }
	}
	else {
	    print "All DST chunks for slot $slot_number have been produced already.\n";
	}
    }
    else {
	die("Run number $first_run_number is out of database\n");
    }
}

#################################################################################################

###########========================== Subroutines: ===============================###############

sub GetPeriod {
    my $dbh_all = &db_open("compass_all","compass","compdb6");
    my $run_number = $_[0];

    my @ret_period = ();
    
    my $sql_query = qq { select DB_NAME, PERIOD, T_MIN, T_MAX from runs_all_mv where RUN_NUMBER=$run_number };
    $debug_flag && print("$sql_query\n");
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

sub TestTBNamesStrings {
    my $run_number   = $_[0];
    my $dbh          = $_[1];
    my $period       = $_[2];
    my $test_slot    = $_[3];
    my $start_string = ":Start of TBNames string:";
    my $stop_string  = ":End of TBNames string:";
    my $slot         = 0;

    my $sql_query = qq { select LOG_INFO from RUNS where RUN_NUMBER=$run_number };
    $dbh->{LongReadLen} = 20 * 1024 ;
    $dbh->{LongTruncOk} = 1;
    my $sth = $dbh->prepare( $sql_query );
    $sth->execute();

    my( $log_info);
    $sth->bind_columns( undef, \$log_info);
    my $is_found = 0;
    if( $sth->fetch() && defined $log_info) {
	$debug_flag && print("$log_info\n");
	if($log_info =~ /^$start_string([\S]+)$stop_string/) {
	    $debug_flag && print("SLOT1: $1\n");
	    $slot = 1;
	    if($slot == $test_slot) {
		$is_found = 1;
	    }
	}
	if($log_info =~ /^$start_string([\S]*)$stop_string$start_string([\S]+)$stop_string/) {
	    $debug_flag && print("SLOT2: $2\n");
	    $slot = 2;
	    if($slot == $test_slot) {
		$is_found = 1;
	    }
	}
	if($log_info =~ /^$start_string([\S]*)$stop_string$start_string([\S]*)$stop_string$start_string([\S]+)$stop_string/) {
	    $debug_flag && print("SLOT3: $3\n");
	    $slot = 3;
	    if($slot == $test_slot) {
		$is_found = 1;
	    }
	}

    }
    $sth->finish();
    $debug_flag && print("Last slot number: $slot\n");
    return $is_found;
}

sub SetTBNamesString {

    my $run_number  = $_[0];
    my $slot_number = $_[1];
    my $year        = $_[2];
    my $period      = $_[3];

    unless(open(TEMPL_OPT, "$ENV{PWD}/template.opt")) {
	print("Cannot open file: $ENV{PWD}/template.opt\n");
	return 1;
    }
    unless(-d "$ENV{PWD}/Run_$run_number") {
	$debug_flag && print("Directory $ENV{PWD}/Run_$run_number does not exist. Create it.\n");
	if(system("mkdir $ENV{PWD}/Run_$run_number") != 0) {
	    print("Cannot create directory: $ENV{PWD}/Run_$run_number\n");
	    close TEMPL_OPT;
	    return 2;
	}
    }
    unless(open(TBNAMES_OPT, "> $ENV{PWD}/Run_$run_number/TBNamesUpload.opt")) {
	print("Cannot create file: $ENV{PWD}/Run_$run_number/TBNamesUpload.opt\n");
	close TEMPL_OPT;
	return 3;
    }

    print TBNAMES_OPT "//DST production with Oracle DB:\n";
    print TBNAMES_OPT "Data CsStore oracle\n";
    print TBNAMES_OPT "Data run select $run_number\n";
    print TBNAMES_OPT "Dst select slot $slot_number\n";
    print TBNAMES_OPT "Data year $year\n";
    print TBNAMES_OPT "Data period $period\n";
    print TBNAMES_OPT "Data type raw\n";
    print TBNAMES_OPT "store TBNames on DST only\n";
    print TBNAMES_OPT "// ********************\n\n";
    while(<TEMPL_OPT>) {
	if($_ =~ /histograms package ROOT/ || $_ =~ /histograms home  trafdic.root/) {
	    next;
	}
	print TBNAMES_OPT $_;
    }
    close TBNAMES_OPT;
    close TEMPL_OPT;

    my $prod_program;

    unless(open(TEMPL_LSF, "$ENV{PWD}/template.lsf")) {
	print("Cannot find file: $ENV{PWD}/template.lsf\n");
	return 4;
    }
    my $search_str = "setenv PRODUCTION_PROGRAM";
    while(<TEMPL_LSF>) {
	if($_ =~ /$search_str\s+(\S+)/) {
	    $prod_program = "$1";
	    $debug_flag && print("Production program: $prod_program\n");
	}
    }
    close TEMPL_LSF;
    unless(defined $prod_program) {
	print("template.lsf: cannot find variable \$prod_program\n");
	return 5;
    }
    if(system(". /afs/cern.ch/project/oracle/script/setoraenv.sh -s 920; $prod_program $ENV{PWD}/Run_$run_number/TBNamesUpload.opt") != 0) {
	return 6;
    }
    return 0;
}

sub Usage {
    print "USAGE: prodPeriod <run number> <number of runs> <slot number> [options]\n";
    print "       where options:  read={raw|dst}\n";
    print "                       test={no|yes}\n";
    print "                       mdst={yes|no}\n";
    print "                        dst={yes|no}\n";
    print "                       rich={rand|yes|no}\n";
    print "                      queue={prod100|prod200|prod400|...}\n";
    print "                        top      - move pending jobs to the top of the queue\n";
    print "                       mail=<your e-mail address>\n";
}
