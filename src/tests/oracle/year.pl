#!/usr/bin/perl -w

use CGI qw(param);

my $YEAR = param("year");

my @cat_file = `cat ./periods_statistic.dat`;
my @config = ();

foreach my $str (@cat_file) {
    chomp $str;
    if($str =~ /^$YEAR/ ) {
	push(@config,$str);
    }
}


print  <<END_of_Multiline_Text;
Content-type: text/html

<html>

    <head>
    <title>COMPASS DST Production Page</title>
    </head>

    <body>

    <p><i><b><font size="4">$YEAR</font></b></i></p>
    <P><i><font size="4">Select slot number of $YEAR year period</i></P>

<table border="1">
  <tr>
    <td width="10%" align="center"><b><font size="4">Period</font></b></td>
    <td width="10%" align="center"><b><font size="4">Database</font></b></td>
    <td width="10%" align="center"><b><font size="4">DST Slots</font></b></td>
    <td width="10%" align="center"><b><font size="4">Date</font></b></td>
    <td width="10%" align="center"><b><font size="4">Runs</font></b></td>
    <td width="10%" align="center"><b><font size="4">N of runs</font></b></td>
    <td width="10%" align="center"><b><font size="4">N of chunks</font></b></td>
    <td width="10%" align="center"><b><font size="4">N RAW events</font></b></td>
    <td width="10%" align="center"><b><font size="4">Size (GB)</font></b></td>
  </tr>

END_of_Multiline_Text

my $n_periods    = 0;
my $n_runs_tot   = 0;
my $n_raw_chunks = 0;
my $n_raw_events = 0;
my $tot_size     = 0;

foreach my $line (@config) {
    if($line =~ /(\d+)\s+(\w+)\s([\d\/\-]+)\s+(\w+)\s+(\d+)\-(\d+)\s+(\d+)\s+(\d+)\s+(\d+)\s+([\d+\.]+)/) {
	my $year      = $1;
	my $period    = $2;
	my $date      = $3;
	my $database  = $4;
	my $first_run = $5;
	my $last_run  = $6;
	my $n_runs    = $7;
	my $n_chunks  = $8;
	my $n_events  = $9;
	my $raw_size  = $10;

	$n_periods++;
	$n_runs_tot   += $n_runs;
	$n_raw_chunks += $n_chunks;
	$n_raw_events += $n_events;
	$tot_size     += $raw_size;

	print "<tr>";
	print "<td width=\"10\%\" align=\"center\">$period</td>";
	print "<td width=\"10\%\" align=\"center\">$database</td>";
	print "<td width=\"10\%\" align=\"center\">";
	if($period ne "03T01" && $period ne "03T02") {
	    print "<a href=\"periods?period=$period\&slot=1\">1</a>\&nbsp;\&nbsp;\&nbsp;\&nbsp;";
	    print "<a href=\"periods?period=$period\&slot=2\">2</a>\&nbsp;\&nbsp;\&nbsp;\&nbsp;";
	    print "<a href=\"periods?period=$period\&slot=3\">3</a></td>";
	}
	else {
	    print "<p>Not available</p>";
	}
	print "<td width=\"10\%\" align=\"center\">\&nbsp; $date </td>";
	print "<td width=\"10\%\" align=\"center\">\&nbsp; $first_run-$last_run </td>";
	print "<td width=\"10\%\" align=\"center\">\&nbsp; $n_runs </td>";
	print "<td width=\"10\%\" align=\"center\">\&nbsp; $n_chunks </td>";
	print "<td width=\"10\%\" align=\"center\">\&nbsp; $n_events </td>";
	print "<td width=\"10\%\" align=\"center\">\&nbsp; $raw_size </td>";
	print "</tr>";
    }
}
print("</table>");
print  "<table border=0>";

print  "<tr><th align=left>Total number of periods for the year: \&nbsp;</th><td align=right> $n_periods</td></tr>";
print  "<tr><th align=left>Total number of runs in DB:           \&nbsp;</th><td align=right> $n_runs_tot</td></tr>";
print  "<tr><th align=left>Total number of RAW chunks:           \&nbsp;</th><td align=right> $n_raw_chunks</td></tr>";
print  "<tr><th align=left>Total number of events in DB:         \&nbsp;</th><td align=right> $n_raw_events</td></tr>";
printf("<tr><th align=left>Total size of RAW files (TB):         \&nbsp;</th><td align=right> %.3f</td></tr>",$tot_size/1024);
print "</table>";

print <<All_Done;

</body>

</html>
All_Done

