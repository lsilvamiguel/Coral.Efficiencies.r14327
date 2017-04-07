#!/usr/bin/perl -w

require 5.004;

use strict;
use diagnostics;
use DBI;
use DBD::Oracle qw(:ora_types);
use CGI qw(param);
use File::stat;
use Time::Local;

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

my $slot        = param("slot");
my $run         = param("run");
my $period_name = "";
my $year        = "";

my @chunk_list      = ();
my @time_list       = ();
my @nevent_list     = ();
my @cur_speed_list  = ();
my @prod_speed_list = ();
my @status_list     = ();
my @logfile_list    = ();

my $debug_flag  = 0;

if($debug_flag) {
    $period_name = "03P1H";
    $slot        = 1;
    $run         = 31216;
}

sub db_open;
sub db_close;
sub getInfo;
sub convertTime;

&main;

print  <<END_of_Multiline_Text;
Content-type: text/html

<html>

    <head>
    <title>Production of run $run (SLOT$slot) check list</title>
    </head>

    <body>

    <p><i><b><font size="4">DST production for the run <a href="runs?period=$period_name&slot=$slot&run=$run">$run</a> of <a href="periods?period=$period_name&slot=$slot">$period_name</a> (SLOT$slot) check-list</font></b></i></p>

<table border="1">
  <tr>
    <td width="10%" align="center"><b><font size="4">Chunk Name</font></b></td>
    <td width="10%" align="center"><b><font size="4">Time of last changing</font></b></td>
    <td width="10%" align="center"><b><font size="4">Number of events have been reconstructed</font></b></td>
    <td width="10%" align="center"><b><font size="4">Current production speed (sec/event)</font></b></td>
    <td width="10%" align="center"><b><font size="4">Production speed (mean value sec/events)</font></b></td>
    <td width="10%" align="center"><b><font size="4">Status</font></b></td>
    <td width="10%" align="center"><b><font size="4">Log file</font></b></td>
  </tr>

END_of_Multiline_Text
    my $cur_time = time();

my $n_chunks   = 0;
my $n_finished = 0;
my $n_events   = 0;

foreach my $chunk (@chunk_list) {
    my $time       = shift(@time_list);
    my $nevent     = shift(@nevent_list);
    my $cur_speed  = shift(@cur_speed_list);
    my $prod_speed = shift(@prod_speed_list);
    my $status     = shift(@status_list);
    my $logref     = shift(@logfile_list);
    $n_chunks++;
    if($status =~ /CORAL stopped: OK/) {
	$n_finished++;
    }
    unless($nevent =~ /nbsp/) {
	$n_events += $nevent;
    }
    $logref =~ /$run\/(\S+)$/;
    my $logfile    = $1;
    my $ctime = $time;
    unless($ctime =~ /nbsp/) {
	if(defined $ctime && $ctime != 0) {
	    $time = convertTime($ctime);
	}
	else {
	    $time = "\&nbsp;";
	    $ctime = 0;
	}
    }
    print "<tr>";
    print "<td width=\"10\%\" align=\"center\">\&nbsp; $chunk </td>";
    if($cur_time-$ctime > 3600 && !($status =~ /CORAL stopped/ || $status =~ /Exited with error code/) && $ctime != 0) {
	print "<td width=\"10\%\" align=\"center\" bgcolor=\"#FF4000\">$time</td>";
	$status = "STUCK";
    }
    else {
	print "<td width=\"10\%\" align=\"center\">\&nbsp; $time</td>";
    }

    print "<td width=\"10\%\" align=\"center\">\&nbsp; $nevent</td>";
    print "<td width=\"10\%\" align=\"center\">\&nbsp; $cur_speed</td>";
    print "<td width=\"10\%\" align=\"center\">\&nbsp; $prod_speed</td>";
    if($status eq  "In production") {
	print "<td width=\"10\%\" align=\"center\" bgcolor=\"#10FF00\">\&nbsp; $status</td>";
    }
    elsif ($status eq "STUCK") {
	print "<td width=\"10\%\" align=\"center\" bgcolor=\"#FF4000\">\&nbsp; $status</td>";
    }
    elsif ($status =~ /Exited with error code/) {
	print "<td width=\"10\%\" align=\"center\" bgcolor=\"#FF0000\">\&nbsp; $status</td>";
    }
    else {
	print "<td width=\"10\%\" align=\"center\">\&nbsp; $status</td>";
    }
    unless(defined $logfile) {
	$logfile = "\&nbsp;";
    }
    unless($logfile =~ /nbsp/ ) {
	print "<td width=\"10\%\" align=\"center\">\&nbsp;<a href=\"show?file=$logref\">$logfile</a></td>";
    }
    else {
	print "<td width=\"10\%\" align=\"center\">\&nbsp; $logfile</td>";
    }
    print "</tr>";    
}

print "</table>";

print "<table border=0>";
print "<tr><th align=left>Total number of chunks for the run:  \&nbsp;</th><td align=right> $n_chunks</td></tr>";
print "<tr><th align=left>Number of fineshed chunks:   \&nbsp;</th><td align=right> $n_finished</td></tr>";
print "<tr><th align=left>Number of produced events:   \&nbsp;</th><td align=right> $n_events</td></tr>";
print "</table>";
print <<All_Done;

    </body>

</html>

All_Done

exit 0;

#################################################################################################

sub main {
    my $dbh = &db_open("compass_all","compass","compdb6");
    my $sql_query = qq { select file_name, period from all_files where run_number=$run and file_type='RAW' order by file_name };
    my $sth = $dbh->prepare($sql_query);
    $sth->execute();
    my $file_name = "";
    $sth->bind_columns( undef, \$file_name, \$period_name );
    while( $sth->fetch() ) {
	if($file_name =~ /cdr(\d+)\-(\d+)\.raw/) {
	    push(@chunk_list,"$2-$1");
	}
    }
    $sth->finish();
    
    $period_name =~ /^(\d\d)(\S+)/;
    my $period = $2;
    $year = "20$1";

    &db_close($dbh);

    foreach my $chunk (@chunk_list) {
	my $log_file = "logs/$year/$period/slot$slot/$run/$chunk-$slot.log";
	my @ls = `ls -l $log_file`;
	if(defined $ls[0] && $ls[0] =~ /$log_file$/ ) {
	    my @log_info = getInfo($log_file);
	    push(@nevent_list,$log_info[0]);
	    push(@cur_speed_list,$log_info[1]);
	    push(@prod_speed_list,$log_info[2]);
	    push(@status_list,$log_info[3]);
	    push(@time_list,$log_info[4]);
	    push(@logfile_list,$log_file);
	}
	else {
	    push(@time_list,"\&nbsp;");
	    push(@nevent_list,"\&nbsp;");
	    push(@cur_speed_list,"\&nbsp;");
	    push(@prod_speed_list,"\&nbsp;");
	    push(@status_list,"not available");
	    push(@logfile_list,"\&nbsp;");
	}
    }
}

#################################################################################################

###########========================== Subroutines: ===============================###############

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

sub getInfo {
    my $log_file = $_[0];
    my @ret_values = ();
    my $nevent = 0;
    my $cur_speed = 0;
    my $prod_speed = 0;
    my $status = "undefined";
    my $time = "\&nbsp;";
    if(open(LOG_FILE,$log_file)) {
	while(<LOG_FILE>) {
	    if($_ =~ /CORAL START\:\s+(\d+)/) {
		my $event_time = $1;
		$status = "CORAL started";
		$time = $event_time;
	    }
	    elsif($_ =~ /DB_SCAN START\:\s+(\d+)/) {
		my $event_time = $1;
		$status = "DB scan started";
		$time = $event_time;
	    }
	    elsif($_ =~ /DB_SCAN END\:\s+(\d+)/) {
		my $event_time = $1;
		$status = "DB scan finished";
		$time = $event_time;
	    }
	    elsif($_ =~ /FIRST_EVENT\:\s+(\d+)/) {
		my $event_time = $1;
		$status = "First event read";
		$time = $event_time;
	    }
	    elsif($_ =~ /EVENT\s+(\d+)\:\s+(\d+)\s+([\d\.]+)\s+([\d\.]+)/) {
		my $event_time = $2;
		$nevent = $1;
		$status = "In production";
		$time = $event_time;
		$cur_speed = $3;
		$prod_speed = $4;
	    }
	    elsif($_ =~ /Number of Processed Events\s+\:\s+(\d+)/) {
		my $event_time = stat($log_file)->ctime;
		$status = "Production finished";
		$nevent = $1;
		$time = $event_time;
	    }
	    elsif($_ =~ /CORAL END\:\s+(\d+)\s+ERROR CODE\:\s+(\d+)/) {
		my $event_time = $1;
		if($2 == 0) {
		    $status = "CORAL stopped: OK";
		}
		else {
		    $status = "Exited with error code: $2";
		}
		$time = $event_time;
	    }
	}
	close LOG_FILE;
    }
    push(@ret_values,$nevent);
    push(@ret_values,$cur_speed);
    push(@ret_values,$prod_speed);
    push(@ret_values,$status);
    push(@ret_values,$time);
    return @ret_values;
}

sub convertTime {
    my $ret_val = 0;
    my ($sec,$min,$hour,$mday,$mon,$year,$wday,$yday,$isdst) =
        localtime($_[0]);
    $year += 1900;
    $mon++;
    $ret_val = sprintf("%02d.%02d.%4d      %02d:%02d:%02d",$mday,$mon,$year,$hour,$min,$sec);
    return $ret_val;
}
