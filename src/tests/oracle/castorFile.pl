#!/usr/local/bin/perl -w

require 5.004;

use strict;
use diagnostics;
use DBI;
use DBD::Oracle qw(:ora_types);
use Getopt::Long;
use Time::Local;

my $noStageQuery=0;
my $help;
my $useOracleVersion = '10102gcc323';

GetOptions(
			  "noStageQuery" => \$noStageQuery,
			  "help" => \$help,
			  "useOracleVersion=s" => \$useOracleVersion,
			 ) || die;



#&initOracleEnv('9206');
&initOracleEnv($useOracleVersion); # call with oracle client version

my $run_number;

$#ARGV == 0 && ($run_number=$ARGV[0]);

die "Unknown option $run_number" unless (!defined($run_number) || $run_number eq "-help" || ($run_number=~/\d+/) );

unless ( defined $run_number )  {
  print("USAGE: \ncastorFile <run_number> [-noStageQuery] look for the corresponding DATE run\ncastorFile -help                        print this help\n");
  exit;
}

sub GetPeriod;
sub db_open;
sub db_close;
sub Cut2Human;

my $debug_flag = 0;

my $dbh;

my@my_period = GetPeriod($run_number);
if($#my_period == 7) {
    my $run_num     = $my_period[0];
    my $run_year    = $my_period[1];
    my $db_name     = $my_period[2];
    my $db_user     = $my_period[3];
    my $db_password = $my_period[4];
    my $db_period   = $my_period[5];
    my $run_t_min   = $my_period[6];
    my $run_t_max   = $my_period[7];
    my $t_min = Cut2Human($run_t_min);
    my $t_max = Cut2Human($run_t_max);
    print("Run number:                 $run_num\n");
    print("Year:                       $run_year\n");
    print("Period:                     $db_period\n");
    print("DB username:                $db_user\n");
    print("DB password:                $db_password\n");
    print("Database:                   $db_name\n");
    print("Start time:                 $t_min\n");
    print("Stop  time:                 $t_max\n");
    print("----------------------------------------------------------------------------------------------------\n");
    print("   CASTOR file name                                                                 Size  Location  \n");
    print("----------------------------------------------------------------------------------------------------\n");
# database handler used for all DBI operations
    $dbh=&db_open($db_user,$db_password,$db_name);
}
else {
    print("Run $run_number is out of DB.\n");
    exit (-1);
}

&main;

&db_close($dbh);

print("\n\n---\nYou can always copy locally the file with the rfcp command:\n");
print("rfcp castor_file_name local_path\n\n");
print("Please do not copy data file on the AFS but use a scratch directory of yours\nand do not forget to clean up the copy of data files when you are finished.\n");
exit(0);

#################################################################################################
sub main {
  my $sql_query = qq { select FILE_DIR, FILE_NAME, FILE_SIZE from FILE_MAPS where FILE_TYPE='RAW' and run_number=$run_number order by FILE_NAME };
  my $sth = $dbh->prepare($sql_query);
  $sth->execute();
  my $file_dir;
  my $file_name;
  my $file_size;
  
  $sth->bind_columns( undef, \$file_dir, \$file_name, \$file_size);
  
  my $total_size = 0;
  my $n_files = 0;
  
  while( $sth->fetch() ) {
	 my $ret = "onDISK";
	 unless ( $noStageQuery ) {
		my $status = `stageqry -M $file_dir$file_name`;
		$ret = "onTAPE" unless $status =~ / STAGED /;
	 } else {
		$ret = "UNCHECKED";
	 }
	 
	 printf("%-74s   %11d  %s\n", "$file_dir$file_name",  $file_size, $ret);
	 $total_size += $file_size;
	 $n_files++;
  }
  print ("----------------------------------------------------------------------------------------------------\n");
  print ("Total number of files: $n_files\n");
  printf("Total size:            %.3f GB (1GB = 1073741824 bytes)\n",$total_size/1024./1024./1024.);
  print ("----------------------------------------------------------------------------------------------------\n");
  
  $sth->finish();
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

sub Cut2Human {
    my ($sec,$min,$hour,$mday,$mon,$year,$wday,$yday,$isdst) =
        localtime($_[0]);
    return sprintf("%02d/%02d/%4d %02d:%02d:%02d",$mday, $mon+1, $year+1900, $hour, $min, $sec);
}

sub initOracleEnv {
  my $oraVersion = shift @_ || '10102gcc323';
  my $CoW = 'Coral Weekly';

  $ENV{ORACLE_CERN} = "/afs/cern.ch/project/oracle";
  unless ( -e $ENV{ORACLE_CERN} ) {
	 die "ERROR. Cannot determine where oracle is at CERN\n";
  }
  $ENV{ORACLE_MOUNT} = "$ENV{ORACLE_CERN}/\@sys";
  $ENV{TNS_ADMIN} = "$ENV{ORACLE_CERN}/admin";
  unless ( -e $ENV{TNS_ADMIN} ) {
	 die "ERROR. Cannot set TNS_ADMIN\n. Please report to $CoW\n";
  }
  unless (defined $ENV{ORACLE_HOME} ) {
	 $ENV{ORACLE_HOME}     = "$ENV{ORACLE_MOUNT}/$oraVersion";
  } # otherwise use the one defined in the environment
  unless ( -e $ENV{ORACLE_HOME} ) {
	 die "ERROR. Cannot use ORACLE_HOME=$ENV{ORACLE_HOME}\. Please report to $CoW\n";
  }
  unless ( defined $ENV{ORA_NLS33} ) {
	 if ( -e "$ENV{ORACLE_HOME}/nls/data" ) {
		$ENV{ORA_NLS33} = "$ENV{ORACLE_HOME}/nls/data";
	 } elsif ( -e "$ENV{ORACLE_HOME}/nls/admin/data" ) {
		$ENV{ORA_NLS33} = "$ENV{ORACLE_HOME}/nls/admin/data";
	 } elsif ( -e "$ENV{ORACLE_HOME}/ocommon/nls/admin/data/" ) {
		$ENV{ORA_NLS33} = "$ENV{ORACLE_HOME}/ocommon/nls/admin/data/";
	 } else {
		die "ERROR. Oracle has changed again ORA_NLS33 path convention. Please report to $CoW\n";
	 }
  } elsif ( ! -e $ENV{ORA_NLS33} ) {
	 die "ERROR. Currently defined ORA_NLS33 does not exist!\n";
  }
  $ENV{DIR_ORACLE}      = "$ENV{ORACLE_HOME}";
  if ( $oraVersion =~ m%^(:?8174|9206)$% ) {
	 $ENV{NLS_LANG} = "american_america.WE8ISO8859P9";
  }
  $ENV{NLS_DATE_FORMAT} = "DD-MON-FXYYYY";
  $ENV{STAGE_HOST}      = "stagecompass";
  $ENV{PATH}            = "$ENV{DIR_ORACLE}/bin:$ENV{PATH}";
  
  $ENV{LD_LIBRARY_PATH} = "$ENV{DIR_ORACLE}/lib" . ($ENV{LD_LIBRARY_PATH} || "");
}
