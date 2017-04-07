#!/usr/bin/perl -w

require 5.004;

use strict;
use diagnostics;
use CGI qw(param);
use Time::Local;

sub convertTime;

my $file_name    = param("file");
my $is_exist     = 1;
my $run_number   = 0;
my $chunk_number = 0;
my $slot_number  = 0;

if($file_name =~ /(\d+)\-(\d+)\-(\d+)\.log$/) {
    $run_number = $1;
    $chunk_number = $2;
    $slot_number = $3;
}
else {
    $is_exist = 0;
}

my @content = ();
if(open(LOG_FILE,$file_name)) {
    while(<LOG_FILE>) {
	push(@content,$_);
    }
    close LOG_FILE;
}
else {
    $is_exist = 0;
}

print  <<END_of_Multiline_Text;
Content-type: text/html

<html>

    <head>
    <title>Production of chunk $run_number-$chunk_number log file</title>
    </head>

    <body>

END_of_Multiline_Text

    if($is_exist == 0) {
	print "<p><font size=\"5\">File $file_name cannot be shown</font></p>";
    }
else {
    print "<table border=\"0\">";
    my $time = "";
    foreach my $line (@content) {
	if($line =~ /CORAL START:\s+(\d+)/) {
	    $time = convertTime($1);
	    print "<tr>";
	    print "<td width=\"10\%\" align=\"left\"><b>CORAL started</b></td><td><i> $time </i></td>";
	    print "</tr>";
	}
	elsif($line =~ /DB_INIT START:\s+(\d+)/) {
	    $time = convertTime($1);
	    print "<tr>";
	    print "<td width=\"10\%\" align=\"left\"><b>Oracle DB initialization started</b></td><td><i> $time </i></td>";
	    print "</tr>";
	}
	elsif($line =~ /DB_INIT END:\s+(\d+)/) {
	    $time = convertTime($1);
	    print "<tr>";
	    print "<td width=\"10\%\" align=\"left\"><b>Oracle DB initialization finished</b></td><td><i> $time </i></td>";
	    print "</tr>";
	}
	elsif($line =~ /DB_SCAN START:\s+(\d+)/) {
	    $time = convertTime($1);
	    print "<tr>";
	    print "<td width=\"10\%\" align=\"left\"><b>Oracle DB scanning started</b></td><td><i> $time </i></td>";
	    print "</tr>";
	}
	elsif($line =~ /DB_SCAN END:\s+(\d+)/) {
	    $time = convertTime($1);
	    print "<tr>";
	    print "<td width=\"10\%\" align=\"left\"><b>Oracle DB scanning finished</b></td><td><i> $time </i></td>";
	    print "</tr>";
	}
	elsif($line =~ /FIRST_EVENT:\s+(\d+)/) {
	    $time = convertTime($1);
	    print "<tr>";
	    print "<td width=\"10\%\" align=\"left\"><b>First event reconstruction started</b></td><td><i> $time </i></td>";
	    print "</tr>";
	}
	elsif($line =~ /EVENT\s+(\d+):\s+(\d+)\s+([\d\.]+)\s+([\d\.]+)/) {
	    $time = convertTime($2);
	    print "<tr>";
	    print "<td width=\"25\%\" align=\"left\"><b>Reconstruction of $1 events finished</b></td><td><i> $time </i></td><td>\&nbsp;$3</td><td>\&nbsp;$4</td>";
	    print "</tr>";
	}
	elsif($line =~ /CORAL END:\s+(\d+)\s+ERROR CODE:\s+(\d+)/) {
	    $time = convertTime($1);
	    print "<tr>";
	    print "<td width=\"25\%\" align=\"left\"><b>CORAL is finished with error code $2</b></td><td><i> $time </i></td>";
	    print "</tr>";
	}
	elsif($line =~ /RUN TIME:\s+(\d+)/) {
	    print "<tr>";
	    print "<td width=\"25\%\" align=\"left\"><b>TOTAL RUN TIME (real time)</b></td><td><i> $1 sec </i></td>";
	    print "</tr>";
	}
	elsif($line =~ /\w+/) {
	    chomp $line;
	    print "<tr>";
	    print "<td width=\"25\%\" align=\"left\"></td><td><pre><font size=\"-1\"><b>$line</b></font></pre></td>";
	    print "</tr>";
	}
    }
}

print <<All_Done;

    </body>

</html>

All_Done

exit 0;

######################################################################################
sub convertTime {
    my $ret_val = 0;
    my ($sec,$min,$hour,$mday,$mon,$year,$wday,$yday,$isdst) =
        localtime($_[0]);
    $year += 1900;
    $mon++;
    $ret_val = sprintf("\&nbsp;  %02d.%02d.%4d \&nbsp; %02d:%02d:%02d",$mday,$mon,$year,$hour,$min,$sec);
    return $ret_val;
}
