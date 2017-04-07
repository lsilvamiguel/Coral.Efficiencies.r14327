#!/usr/bin/perl -w
use CGI qw(param);

sub IsHistosMerged;

my $period = param("period");
$period =~ /^(\d\d)(\S+)/;
my $year = "20$1";
my $db_period = $2;
my $slot_number = param("slot");
my @ls = `./oradbInfo.pl $period $slot_number`;

print  <<END_of_Multiline_Text;
Content-type: text/html

<HTML>
	<HEAD>
	<TITLE>Period $db_period of $year</TITLE>
	</HEAD>
	<BODY>
	<H1>$db_period period of $year (SLOT$slot_number)</H1>


END_of_Multiline_Text

print "<table border=1>";
print "<font face=\"Courier New\" size=\"1\">";
print "<tr>";
print "<th>Run #</th><th>Status</th><th>N.chunks</th><th>N.DST</th><th>todo</th><th>N.mDST</th><th>Is merged</th><th>histo file</th><th>Phast version</th><th>Commentaries</th>";
print "<th>Check production</th>";
print "</tr>";

my $n_runs_tot      = 0;
my $n_chunks_tot    = 0;
my $n_dst_tot       = 0;
my $n_mdst_tot      = 0;
my $n_mdst_merged   = 0;
my $n_runs_finished = 0;
my $n_warnings      = 0;
my $n_good_runs     = 0;
my $n_good_chunks   = 0;

foreach my $f (@ls) {
    chomp $f;
    if($f =~ /^(\d+)\s+(\d+)\s+(\d+)\s+(\d+)\s+(\d+)\s+(\d+)\s+(\w+)\s*(\d*)\s*([\<\-\ \w\%\d]*)/) {
	my $run_number    = $1;
	my $run_status    = $2;
	my $n_chunks      = $3;
	my $n_dsts        = $4;
	my $todo          = $5;
	my $n_mdst        = $6;
	my $is_merged     = $7;
	my $phast_version = $8;
	my $commentary    = $9;

	$n_runs_tot++;
	$n_chunks_tot += $n_chunks;
	$n_dst_tot    += $n_dsts;
	$n_mdst_tot   += $n_mdst;
	my $is_histos = IsHistosMerged($run_number);
	if($is_merged eq "YES") {
	    $n_mdst_merged++;
	}
	if($todo == 0) {
	    $n_runs_finished++;
	}

	unless($phast_version =~ /\S+/) {
	    $phast_version = "<td>\&nbsp;</td>";
	}
	else {
	    $phast_version = "<td align=center>$phast_version</td>";
	}
	unless($commentary =~ /\S+/) {
	    $commentary = "<td>\&nbsp;</td>";
	}
	elsif($commentary =~ /DST production is not finished yet/) {
	    $commentary = "<td align=center bgcolor=\"#FF8000\">$commentary</td>";
	    $n_warnings++;
	}
	elsif($commentary =~ /N DST is not equal N mDST/) {
	    $commentary = "<td align=center bgcolor=\"#FF1000\">$commentary</td>";
	    $n_warnings++;
	}
	elsif($commentary =~ /Not merged yet/) {
	    $commentary = "<td align=center bgcolor=\"#10FF00\">$commentary</td>";
	    $n_warnings++;
	}
	elsif($commentary =~ /Less than 5 percent of the chunks left to be done/) {
	    $commentary = "<td align=center bgcolor=\"#FFFF00\">$commentary</td>";
	    $n_warnings++;
	}
	else {
	    $commentary = "<td align=center>$commentary</td>";
	}
	print "<tr>";
	my $is_good_run = $run_status%10;
	if($is_good_run) {
	    $n_good_runs++;
	    $n_good_chunks += $n_chunks;
	    print "<td align=center bgcolor=\"#A0FFA0\">";
	}
	else {
	    print "<td align=center>";
	}
	print "<a href=\"http://compass-dst-production.web.cern.ch/compass-dst-production/cgi-bin/runs?period=$period&slot=$slot_number&run=$run_number\">";
	print "$run_number</a></td>";
	print "<td align=center>$run_status</td>";
	print "<td align=center>$n_chunks</td>";
	if($n_dsts != ($n_chunks-$todo)) {
	    print "<td align=center bgcolor=\"#A08080\">$n_dsts</td>";
	}
	else {
	    print "<td align=center>$n_dsts</td>";
	}
	print "<td align=center>$todo</td>";
	if($n_mdst != $n_dsts) {
	    print "<td align=center bgcolor=\"#FFA040\">$n_mdst</td>";
	}
	else {
	    print "<td align=center>$n_mdst</td>";
	}
	if($is_good_run && $is_merged eq "YES") {
	    print "<td align=center  bgcolor=\"#A0FFA0\">";
	}
	elsif($is_good_run && $is_merged eq "NO") {
	    print "<td align=center  bgcolor=\"#FFA0A0\">";
	}
	else {
	    print "<td align=center>";
	}
	print "$is_merged</td>";
	if($is_merged eq "YES" && $is_histos eq "NO") {
	    print("<td align=center  bgcolor=\"#FFFF80\">$is_histos</td>");
	}
	else {
	    print("<td align=center>$is_histos</td>");
	}
	print $phast_version;
	print $commentary;
	print "<td align=center><a href=production?run=$run_number&slot=$slot_number> \&nbsp;check\&nbsp;</td>";
	print "</tr>";
    }
}
print "</font>";
print "</table>";
print "<table>";
print "<table border=0>";
print "<tr><th align=left>Total number of runs for the period: \&nbsp;</th><td align=right> $n_runs_tot</td></tr>";
print "<tr><th align=left>Total number of chunks:              \&nbsp;</th><td align=right> $n_chunks_tot</td></tr>";
print "<tr><th align=left>Number of good runs:                 \&nbsp;</th><td align=right> $n_good_runs</td></tr>";
print "<tr><th align=left>Number of good chunks:               \&nbsp;</th><td align=right> $n_good_chunks</td></tr>";
print "<tr><th align=left>Total number of produced DST chunks: \&nbsp;</th><td align=right> $n_dst_tot</td></tr>";
print "<tr><th align=left>Total number of mini-DST chunks:     \&nbsp;</th><td align=right> $n_mdst_tot</td></tr>";
print "<tr><th align=left>Total number of merged mini-DST's:   \&nbsp;</th><td align=right> $n_mdst_merged</td></tr>";
print "<tr><th align=left>Total number of completed runs:      \&nbsp;</th><td align=right> $n_runs_finished</td></tr>";
print "<tr><th align=left>Number of warnings:                  \&nbsp;</th><td align=right> $n_warnings</td></tr>";
print "</table>";
print <<All_Done;

	</BODY>
</HTML>
All_Done


    exit 0;
#################################################################################################

sub IsHistosMerged {
    my ($run_number) = @_;
    my $result = `/usr/local/bin/nsls -l /castor/cern.ch/compass/data/$year/oracle_dst/$db_period/histos/slot$slot_number/histsum-$run_number-$slot_number.root`;
    if($result =~ /^mrw\S+\s+\S+\s+\S+\s+\S+\s+(\d+)/) {
	if($1 > 0) {
	    return "YES";
	}
    }
    return "NO";
}
