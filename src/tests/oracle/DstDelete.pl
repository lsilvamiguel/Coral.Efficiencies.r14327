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

my $min_phast_version = 5;
my $max_phast_version = 7;

my $debug_flag   = 0; 

my $run_number   = 0;
my $slot_number  = 0;
my $chunk_name   = " ";
my $chunk_number = 0;
my $clear_all    = 0;
my $password     = "";

foreach my $param (@ARGV) {
    if($param =~ /RUN\=(\d+)/i ) {
	$run_number = $1;
    }
    elsif ($param =~ /CHUNK\=\D*(\d+)\-(\d+)/i ) {
	$chunk_number = $1;
	$run_number   = $2;
	$chunk_name   = "$1-$2";
    }
    elsif ($param =~ /SLOT\=(\d)/i ) {
	$slot_number = $1;
    }
    elsif ($param =~ /PASSWORD\=(\S+)/i ) {
	$password = $1;
    }
    elsif ($param =~ /CLEAR\=ALL/i ) {
	$clear_all = 1;
    }
}

if($debug_flag == 1) {
    print("run_number = $run_number\nchunk_name=$chunk_name\nslot_number=$slot_number\nclear_all=$clear_all\n");
}

if($run_number == 0 || $slot_number == 0) {
    die("Usage: DstDelete { RUN=<run number> | CHUNK=<chunk name>} SLOT=<slot number> [CLEAR=ALL] [password=*******]\n");    
}

my $username = `whoami`;
chomp($username);
my $dbh;

sub GetPeriod;
sub db_open;
sub db_close;
sub DeleteChunk;

my $year;
my $period;

my@my_period = GetPeriod($run_number);
if($#my_period == 7) {
    $year        = $my_period[1];
    my $database = $my_period[2];
    my $user     = $my_period[3];
    my $pwd      = $my_period[4];
    $period      = $my_period[5];
    if($password eq "") {
	print("Password: ");
	$password = <STDIN>;
	chomp $password;
	unless($pwd eq $password) {
	    die "Incorrect password.\n";
	}
    }
# database handler used for all DBI operations
    $dbh=&db_open($user,$password,$database);
}
else {
    die("Run $run_number is out of DB.\n");
}

&main;

&db_close($dbh);
exit(0);


########################### MAIN ##############################################################

sub main {
    unless($chunk_name eq " ") {
	DeleteChunk($chunk_name);
    }
    else { #delete all chunk of run
	my $sql_query = qq { select FILE_NAME from FILE_MAPS where RUN_NUMBER=$run_number };
	$debug_flag && print("SQL query: $sql_query\n");
	my $sth = $dbh->prepare($sql_query);
	$sth->execute();
	my @file_names = ();
	my $file_name = "";
	$sth->bind_columns( undef, \$file_name );
	while( $sth->fetch() ) {
	    push(@file_names,$file_name);
	}
	$sth->finish();
	foreach my $f (@file_names) {
	    $f =~ /cdr(\d\d\d\d\d)-(\d\d\d\d\d)/;
	    my $chunk_n = $1;
	    my $run_n   = $2;
	    if($run_n == $run_number) {
		DeleteChunk("$chunk_n-$run_n");
	    }
	}
	if($clear_all) {
	    # Remove merged files:
	    my $castor_dir = "/castor/cern.ch/compass/data/$year/oracle_dst/$period";
	    
	    # Remove merged mini-DST:
	    for(my $phast_version = $min_phast_version; $phast_version <= $max_phast_version; $phast_version++) { 
		my $ret_value = system("/usr/local/bin/rfrm $castor_dir/mDST/mDST-$run_number-$slot_number-$phast_version.root 2> /dev/null");
		if($ret_value == 0) { # file is found and removed successfully
		    print "File $castor_dir/mDST/mDST-$run_number-$slot_number-$phast_version.root has been removed\n";
		    last;
		}
	    }
	    #Remove merged histogram file:
	    my $ret_value = system("/usr/local/bin/rfrm $castor_dir/histos/slot$slot_number/histsum-$run_number-$slot_number.root");
	    if($ret_value == 0) {
		print "File $castor_dir/histos/slot$slot_number/histsum-$run_number-$slot_number.root has been removed\n";
	    }
	    #Remove Run_XXXXX directory:
	    $ret_value = system("rm -rf ./Run_$run_number");
	    if($ret_value == 0) {
		print "Directory $ENV{PWD}/Run_$run_number has been removed.\n";
	    }

	    # Remove TBnames string for the slot:

	    $ret_value = system("./removeTBnames $run_number $slot_number");
	    if($ret_value == 0) {
		print "TBname string is removed.\n";
	    }
	}
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
    my $db_user = "compdst_";
    my $db_password = "dstwriter";
    my $period;
    my $t_min;
    my $t_max;
    $sth->bind_columns( undef, \$db_name, \$period, \$t_min, \$t_max );
    if( $sth->fetch() ) {
	$debug_flag && print("Run:      $run_number\nDatabase: $db_name\nPeriod:   $period\nt_min:    $t_min\nt_max:    $t_max\n");
	$db_user = "$db_user$period";
	$period =~ /^(\d\d)(\w+)/;
	$period = $2;
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

sub DeleteChunk {
    my $chunk_name = $_[0];
    $chunk_name =~ /(\d\d\d\d\d)-(\d\d\d\d\d)/;
    my $chunk_number = $1;
    my $run_number   = $2;
    my $sql_query = qq { select FILE_ID, FILE_DIR, FILE_NAME from DST_FILES where  FILE_NAME like '\%$chunk_name\%'
			     and RUN_NUMBER=$run_number and DST_VERSION=$slot_number };
    $debug_flag && print("SQL query: $sql_query\n");
    my $sth = $dbh->prepare($sql_query);
    $sth->execute();
    my $file_id = 0;
    my $file_name = "";
    my $file_dir = "";
    $sth->bind_columns( undef, \$file_id, \$file_dir, \$file_name );
    if( $sth->fetch() ) {	
	$sth->finish();
	print("Delete chunk $file_name from database.\n");
	$sql_query = qq { delete from DST_FILES where FILE_ID=$file_id };
	$debug_flag && print("SQL query: $sql_query\n");
	$sth = $dbh->prepare($sql_query);
	$sth->execute();
	$sth->finish();
	$dbh->commit;
	print "Chunk $chunk_name is removed from Oracle DB.\n";
    }
    else {
	print "Cannot find chunk $chunk_name in DB\n";
    }
    if($clear_all) {
	my $mdst_chunk_name = "";
	if($file_name ne "") {
	    my $ret_value = system("/usr/local/bin/rfrm $file_dir/$file_name");
	    if($ret_value == 0) {
		print "$file_dir/$file_name removed.\n";
	    }
	}
	my $ret_value = system("rm ./Run_$run_number/$chunk_number-$run_number.*");
	if($ret_value == 0) {
	    print "The files $ENV{PWD}/Run_$run_number/$chunk_number-$run_number.* have been deleted.\n";
	}
	$ret_value = system("rm /afs/cern.ch/compass/scratch/mdst2/production/logs/$year/$period/slot$slot_number/$run_number/$run_number-$chunk_number-$slot_number.log");
	if($ret_value == 0) {
	    print "The file /afs/cern.ch/compass/scratch/mdst2/production/logs/$year/$period/slot$slot_number/$run_number/$run_number-$chunk_number-$slot_number.log has been deleted.\n";
	}
#
# update ./Run_NNNNN/.info file:
#
	if(open(INFO_FILE,"./Run_$run_number/.info")) {
	    if(open(NEW_FILE,"> ./Run_$run_number/.info.new")) {
		while(<INFO_FILE>) {
		    if($_ =~ /$run_number-$chunk_number/ ) {
			if($_ =~ /^(\S+)\s+1/) {
			    $mdst_chunk_name = $1;
			}
			next;
		    }
		    print NEW_FILE $_;
		}
		close NEW_FILE;
		close INFO_FILE;
		$ret_value = system("rm ./Run_$run_number/.info; mv ./Run_$run_number/.info.new ./Run_$run_number/.info");
		if($ret_value == 0) {
		    print "Chunk $run_number-$chunk_number is removed from .info file.\n";
		}
	    }
	    else {
		print("Cannot open temporary .info file. Skip update.\n");
		close INFO_FILE;
	    }
	}
	else {
	    print("Cannot open .info file. Update is skipped.\n");
	}

	$ret_value = system("rm ./logFiles/optdst$slot_number-$chunk_number-$run_number.log");
	if($ret_value == 0) {
	    print "File ./logFiles/optdst$slot_number-$chunk_number-$run_number.log has been removed.\n";
	}
		
	$ret_value = system("/usr/local/bin/rfrm /castor/cern.ch/compass/data/$year/oracle_dst/$period/RICH/histodst$slot_number-$chunk_name.gfile");
	if($ret_value == 0) {
	    print "File /castor/cern.ch/compass/data/$year/oracle_dst/$period/RICH/histodst$slot_number-$chunk_name.gfile has been removed.\n";
	}

	$ret_value = system("/usr/local/bin/rfrm /castor/cern.ch/compass/data/$year/oracle_dst/$period/TRAFDIC/$run_number-$chunk_number-$slot_number.root");
	if($ret_value == 0) {
	    print "File /castor/cern.ch/compass/data/$year/oracle_dst/$period/TRAFDIC/$run_number-$chunk_number-$slot_number.root has been removed.\n";
	}

	for(my $phast_version = $min_phast_version; $phast_version <= $max_phast_version; $phast_version++) {
	    $ret_value = system("/usr/local/bin/rfrm /castor/cern.ch/compass/data/$year/oracle_dst/$period/mDST.chunks/mDST-$run_number-$chunk_number-$slot_number-$phast_version.root 2> /dev/null");
	    if($ret_value == 0) {
		print "/castor/cern.ch/compass/data/$year/oracle_dst/$period/mDST.chunks/mDST-$run_number-$chunk_number-$slot_number-$phast_version.root removed.\n";
		last;
	    }
	}

#Clean up compass009d: 
	if($mdst_chunk_name ne "") {
	    $ret_value = system("/usr/local/bin/rfrm $mdst_chunk_name 2> /dev/null");
	    if($ret_value == 0) {
		print "$mdst_chunk_name removed.\n";
	    }
	    my $trafdic_file = "/shift/compass009d/data05/TRAFDIC/$run_number-$chunk_number-$slot_number.root";
	    $ret_value = system("/usr/local/bin/rfrm $trafdic_file");
	    if($ret_value == 0) {
		print "$trafdic_file removed\n";
	    }
	}
    }
}
