#!/usr/bin/perl -w

require 5.004;

use strict;
use diagnostics;
use DBI;
use DBD::Oracle qw(:ora_types);
use CGI qw(param);

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

my $first_phast_version = 5;
my $last_phast_version  = 6;

my $period_name = param("period");
my $slot        = param("slot");
my $run         = param("run");
my $year        = 0;
my $dst_period  = "";


my $debug_flag  = 0;

if($debug_flag) {
    $period_name = "02P2D";
    $slot        = 1;
    $run         = 22018;
}

my @html_data   = ();

sub GetPeriodRuns;
sub db_open;
sub db_close;

print  <<END_of_Multiline_Text;
Content-type: text/html

<html>

    <head>
    <title>Run $run SLOT: $slot PERIOD: $period_name</title>
    </head>

    <body>

    <p><i><b><font size="4">DST production for the run $run of <a href="periods?period=$period_name&slot=$slot">$period_name</a> (SLOT$slot)</font></b></i></p>
<table border="1">
  <tr>
    <td width="10%" align="center"><b><font size="4">Chunk Name</font></b></td>
    <td width="10%" align="center"><b><font size="4">RAW file size</font></b></td>
    <td width="10%" align="center"><b><font size="4">N RAW events</font></b></td>
    <td width="10%" align="center"><b><font size="4">DST file size</font></b></td>
    <td width="10%" align="center"><b><font size="4">N DST events</font></b></td>
    <td width="10%" align="center"><b><font size="4">Mini-DST file size</font></b></td>
  </tr>

END_of_Multiline_Text


&main;

my $n_chunks           = 0;
my $tot_raw_file_size  = 0;
my $n_raw_events       = 0;
my $tot_dst_file_size  = 0;
my $n_dst_events       = 0;
my $n_dst_files        = 0;
my $n_mdst_files       = 0;
my $mdst_file_size     = 0;
my $merged_file_size   = 0;

foreach my $h (@html_data) {
    if($h =~ /^\d+ \S+ (\S+) (\d+) (\d+)\s*\d*\s*\d*\s*\S*\s*\S*\s*(\d*)\s*(\d*)/ ) {
	my $raw_file_name = $1;
	my $raw_file_size = $2;
	my $n_events      = $3;
	my $dst_file_size = $4;
	my $dst_n_events  = $5;
	$raw_file_name =~ /cdr(\d+)\-(\d+)\.raw/;
	my $chunk_number = $1;
	my $run_number   = $2;
	$n_chunks++;
	$tot_raw_file_size += $raw_file_size;
	$n_raw_events += $n_events;
	my $mini_dst_size = 0;
	for(my $phast_version = $first_phast_version; $phast_version <= $last_phast_version; $phast_version++) {
	    my $mini_dst_info = `/usr/local/bin/nsls -l /castor/cern.ch/compass/data/$year/oracle_dst/$dst_period/mDST.chunks/mDST-$run_number-$chunk_number-$slot-$phast_version\.root`; 
	    if($mini_dst_info =~ /root/) {
		$mini_dst_info =~ /^\S+\s+\S+\s+\S+\s+\S+\s+(\d+)/;
		$mini_dst_size = $1;
		$mdst_file_size += $mini_dst_size;
		$n_mdst_files++;
		last;
	    }
	}
	print "<tr>";
	print "<td width=\"10\%\" align=\"center\">\&nbsp; $chunk_number\-$run_number </td>";
	print "<td width=\"10\%\" align=\"center\">\&nbsp; $raw_file_size </td>";
	print "<td width=\"10\%\" align=\"center\">\&nbsp; $n_events </td>";

	if($dst_file_size != 0) {
	    print "<td width=\"10\%\" align=\"center\">\&nbsp; $dst_file_size </td>";
	    print "<td width=\"10\%\" align=\"center\">\&nbsp; $dst_n_events </td>";
	    $tot_dst_file_size += $dst_file_size;
	    $n_dst_files++;
	    $n_dst_events += $dst_n_events;
	}
	else {
	    print "<td width=\"10\%\" align=\"center\">\&nbsp;</td>";
	    print "<td width=\"10\%\" align=\"center\">\&nbsp;</td>";
	}
	if($mini_dst_size != 0) {
	    print "<td width=\"10\%\" align=\"center\">\&nbsp; $mini_dst_size</td>";
	}
	else {
	    print "<td width=\"10\%\" align=\"center\">\&nbsp;</td>";
	}
	print "</tr>";
    }
}

my $merged_mdst_file_name = "";

for(my $phast_version = $first_phast_version; $phast_version <= $last_phast_version; $phast_version++) {
    my $merged_info = `/usr/local/bin/nsls -l /castor/cern.ch/compass/data/$year/oracle_dst/$dst_period/mDST/mDST-$run-$slot-$phast_version\.root`;
    if($merged_info =~ /root/) {
	$merged_info =~ /^\S+\s+\S+\s+\S+\s+\S+\s+(\d+)/;
	$merged_file_size = $1;
	$merged_mdst_file_name = "/castor/cern.ch/compass/data/$year/oracle_dst/$dst_period/mDST/mDST-$run-$slot-$phast_version\.root";
	last;
    }
}

print  "</table>";

print  "<table border=0>";
print  "<tr><th align=left>Total number of chunks for the run:  \&nbsp;</th><td align=right> $n_chunks</td></tr>";
printf("<tr><th align=left>Total size of RAW files (GB):        \&nbsp;</th><td align=right> %.3f</td></tr>",$tot_raw_file_size/1024/1024/1024);
print  "<tr><th align=left>Total number of RAW events:          \&nbsp;</th><td align=right> $n_raw_events</td></tr>";
print  "<tr><th align=left>Total number of DST chunks:          \&nbsp;</th><td align=right> $n_dst_files</td></tr>";
printf("<tr><th align=left>Total size of DST files (MB):        \&nbsp;</th><td align=right> %.3f</td></tr>",$tot_dst_file_size/1024/1024);
print  "<tr><th align=left>Total number of DST events:          \&nbsp;</th><td align=right> $n_dst_events</td></tr>";
print  "<tr><th align=left>Total number of mini-DST files:      \&nbsp;</th><td align=right> $n_mdst_files</td></tr>";
printf("<tr><th align=left>Total size of mini-DST files (MB):   \&nbsp;</th><td align=right> %.3f</td></tr>",$mdst_file_size/1024/1024);
printf("<tr><th align=left>Size of merged mini-DST file (MB):   \&nbsp;</th><td align=right> %.3f</td></tr>",$merged_file_size/1024/1024);

print "</table>";

if($merged_mdst_file_name ne "") {
    print "<p>Merged mini-DST file name: $merged_mdst_file_name</p>";
}

print <<All_Done;

    </body>

</html>


All_Done
    exit 0;
#################################################################################################
sub main {
    my @my_period = GetPeriodRuns($period_name);
    if($#my_period == 6) {
	$dst_period  = $my_period[0];
	my $first_run   = $my_period[1];
	my $last_run    = $my_period[2];
	my $database    = $my_period[3];
	my $user_name   = $my_period[4];
	my $password    = $my_period[5];
	$year        = $my_period[6];

	my $dbh = &db_open($user_name,$password,$database);

	my @raw_data = ();
	my @rawev_data = ();

	my $sql_query = qq { select FILE_ID, FILE_DIR, FILE_NAME, FILE_SIZE from FILE_MAPS where RUN_NUMBER=$run and FILE_TYPE='RAW' order by FILE_NAME };
	my $file_id    = 0;
	my $file_dir   = "";
	my $file_name  = "";
	my $file_size  = 0;
	my $sth = $dbh->prepare($sql_query);
	$sth->execute();
	$sth->bind_columns( undef, \$file_id, \$file_dir, \$file_name, \$file_size);
	while( $sth->fetch() ) {
	    my $raw = "$file_id $file_dir $file_name $file_size";
	    push(@raw_data,$raw);
	}
	$sth->finish();
	foreach my $r (@raw_data) {
	    $r =~ /^(\d+)/;
	    $file_id = $1;
	    $sql_query = qq { select count(*) from EVENT_HEADERS where RAWEV_FILE_ID=$file_id };
	    my $n_events = 0;
	    $sth = $dbh->prepare($sql_query);
	    $sth->execute();
	    $sth->bind_columns( undef, \$n_events);
	    $sth->fetch();
	    push(@rawev_data,"$r $n_events");
	    $sth->finish();
	}
	my @dst_data = ();

	my $rawev_file_id = 0;

	$sql_query = qq { select RAWEV_FILE_ID, FILE_ID, FILE_DIR, FILE_NAME, FILE_SIZE from DST_FILES where RUN_NUMBER=$run and DST_VERSION=$slot };
	$sth = $dbh->prepare($sql_query);
	$sth->execute();
	$sth->bind_columns( undef, \$rawev_file_id, \$file_id, \$file_dir, \$file_name, \$file_size);
	while( $sth->fetch() ) {
	    my $dst = "$rawev_file_id $file_id $file_dir $file_name $file_size";
	    push(@dst_data,$dst);
	}
	$sth->finish();
	

	my @dstev_data = ();

	foreach my $d (@dst_data) {
	    $d =~ /^\d+ (\d+)/;
	    $file_id = $1;
	    $sql_query = qq { select count(*) from DST where FILE_ID=$file_id };
	    my $n_events = 0;
	    $sth = $dbh->prepare($sql_query);
	    $sth->execute();
	    $sth->bind_columns( undef, \$n_events);
	    $sth->fetch();
	    push(@dstev_data,"$d $n_events");
	    $sth->finish();
	}

	&db_close($dbh);

	foreach my $raw_file (@rawev_data) {
	    $raw_file =~ /^(\d+)/;
	    $rawev_file_id = $1;
	    my $has_dst = 0;
	    foreach my $dst_file (@dstev_data) {
		$dst_file =~ /^(\d+)/;
		if($1 == $rawev_file_id) {
		    my $html = "$raw_file $dst_file";
		    push(@html_data,$html);
		    $has_dst = 1;
		    last;
		}
	    }
	    if($has_dst == 0) {
		push(@html_data,$raw_file);
	    }
	}
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
