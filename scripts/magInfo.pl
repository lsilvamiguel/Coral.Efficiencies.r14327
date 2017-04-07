#!/usr/bin/perl -w

# $Id: magInfo.pl 13608 2012-11-12 22:59:44Z ybedfer $

# GET MAGNET AND TARGET RELATED INFO from MySQLDB, viz.:
#  SOLENOID CURRENT, DIPOLE CURRENT, POLAR UP and DOWN +/- ERROR and SM2 NMR
# Usage: magInfo.pl <min run>..<max run>

#  All runs are output, whether they do have relevant info or not.
#  Except for magnet currents, if info is missing, the most recent available
# info is printed. This can be misleading. For clarity then, the run# to which
# the printed info pertains is also printed, which may turn out not informative
# enough: better would be to print both the run# AND a time elapsed since the
# printed value was measured (can vary greatly depending upon the #spills of
# runs involved).
# => In order to benefit from this behaviour right from the start, include some
#   safety margin in your specification of the run# range, ahead of the 1st run
#   you want to target. (Note: A margin of 20 runs is systematically applied.)

#  Output may be filed, possibly edited, and then used as a "CsMagInfo File"
# file in coral's options file.
#  This may be interesting:
# - When working on computer system having no access to a COMPASS server.
# - When the info returned by MySQLDB does not fit, for whatever reason, it
#  provides a work around (e.g. because allowing to update the magnet info by
#  editing the file).

#  Derived from Colin's "CSoft/Perl/polar.pl"
#  Modified to:
#  - Include NMR.
#  - List all runs, assigning infoless runs the info of previous one.
#  - Retrieve solenoid/dipole currents from "tb_tgtcurrent" table.
#  - Retrieve SM1/2 currents from "tb_beamvalues".

# $Log: magInfo.pl,v $
# Revision 1.4  2010/11/23 16:50:11  ybedfer
#  Cosmetics...
#
# Revision 1.3  2010/11/21 00:29:39  ybedfer
#  - Bug fix: Check the dimension of "row4" before testing it [2] component.
#  - Compute and print averages.
#
# Revision 1.2  2010/08/02 10:38:58  ybedfer
#  Get the sign of the solenoid current from "SOLN-DCCT_ai", when available.
#
# Revision 1.1  2010/06/16 23:30:24  ybedfer
# Initial revision
#

use strict;
use DBI;
##use lib qw(/afs/cern.ch/compass/detector/monitor/Perl/lib/perl5
##           /afs/cern.ch/compass/detector/monitor/Perl/lib/perl5/site_perl/5.00503/i386-linux);

if ($#ARGV<0 || $ARGV[0] =~ /--help/ || $ARGV[0] =~ /-h/) {
    print "magInfo:\a GET MAGNET AND TARGET RELATED INFO from MySQLDB, viz.:\n";
    print "  SOLENOID/DIPOLE CURRENTS, POLAR UP/DOWN+/-ERROR, SM2 NMR, SM1/2 CURRENTS\n";
    print "Usage : magInfo.pl <minRun..maxRun>\n";
    print "\n";
print " All runs are output, whether they do have relevant info or not.\n";
print " Except for magnet currents, if info is missing, the most recent available info\n";
print "is printed (the run# to which the printed info pertains being specified.)\n";
print "=> In order to benefit from this behaviour right from the start, include some\n";
print "  safety margin in your specification of the run# range, ahead of the 1st run\n";
print "  you want to target.\n";
    exit(1);
}

my $dbname = "runlb";
my $host = "wwwcompass.cern.ch";
my $db = "dbi:mysql:database=runlb;host=$host";
my $dbh=0;

$dbh = DBI->connect($db,'anonymous','')
    || die "Can't connect to the DB : $DBI::errstr\n";

my $runs = shift;
my $minrun = 0;
my $maxrun = 0;

########## PARSING ARGUMENT ##########
if($runs) {
    my $now = localtime(time);
    if($runs =~ /^(\d+)?\.{2}(\d+)?$/) {
	
	$1 && ($minrun = $1);
	$2 && ($maxrun = $2);
	if ($maxrun<$minrun) {
	    print "** magInfo: Wrong arg list: maxRun < minRun!\n\n";
	    print "Usage : magInfo.pl <minRun>..<maxRun>\n";
	    exit(1);
	}
	print "/* magInfo.pl: runs $minrun..$maxrun - $now */\n";
    }
    elsif($runs =~ /^(\d+)$/) {
	$minrun = $1;
	$maxrun = $1;
	print "/* magInfo.pl:run $minrun... - $now */\n";  # N.B.: Leaving a blank space ahead of "*/"
    }
    else {
	print "magInfo:\a GET MAGNET AND TARGET RELATED INFO from MySQLDB, viz.:\n";
	print "  SOLENOID/DIPOLE CURRENTS, POLAR UP/DOWN+/-ERROR, SM2 NMR, SM1/2 CURRENTS\n";
	print "Usage : magInfo.pl <minRun>..<maxRun>\n";
	print "\n";
print " All runs are output, whether they do have relevant info or not.\n";
print " Except for magnet currents, if info is missing, the most recent available info\n";
print "is printed (the run# to which the printed info pertains being specified.)\n";
print "=> In order to benefit from this behaviour right from the start, include some\n";
print "  safety margin in your specification of the run# range, ahead of the 1st run\n";
print "  you want to target.\n";
	exit(1);
    }
}

########## START SOMEWHAT AHEAD OF $minrun,
# FOR IN CASE NO DB INFO FOR run==$minrun
my $strtrun = $minrun - 20;

########## LIST of VARIABLES RETRIEVED from DB
my $run     = 0; my $type     = 0;
my $soleno  = 0; my $dipole   = 0;  # "tb_tgtcurrent"
my $uppol   = 0; my $updpol   = 0;  # "tb_offlinepolar"
my $dwnpol  = 0; my $dwndpol  = 0;
my $polInfo = 0; my $polRef = 0; # run# from which info comes if absent for current one
my $nmr     = 0;                    # "tb_NMRSM2"
my $nmrInfo = 0; my $nmrRef = 0; # run# from which info comes if absent for current one
my $sm1     = 0; my $sm2      = 0;  # "tb_beamvalues"
##### Average values
# Average over all runs
my $nRuns = 0;
my $types = -2; # -2: not yet set, -1: ambiguous, else cf. runlb/tb_runtype
my $sSoleno = 0; my $sDipole = 0 ; my $nTargets = 0;
my $sUppol = 0; my $sDwnpol = 0; my $sUpdpol2 = 0; my $sDwndpol2 = 0; my $nPols = 0;
my $sNMR = 0; my $nNMRs = 0; my $sSM1 = 0; my $sSM2 = 0; my $nSMs = 0;
# Average over physics runs (Note: same specifiers as supra, but capitalized)
my $physicsTypes = 0xce61; # All physics types, up to bit 2^15 = DVCS.
my $NRuns = 0; my $Types = -2;
my $SSoleno = 0; my $SDipole = 0 ; my $NTargets = 0;
my $SUppol = 0; my $SDwnpol = 0; my $SUpdpol2 = 0; my $SDwndpol2 = 0; my $NPols = 0;
my $SNMR = 0; my $NNMRs = 0; my $SSM1 = 0; my $SSM2 = 0; my $NSMs = 0;

########## HEADER ##########
print "/* run      solen   dipole   up    +/-    down   +/-  (<run#)      NMR  (<run#)    SM1     SM2 */\n";
########## FORMATTING (NOTA BENE the terminating "." dot.) ##########
format STDOUT =
 @<<<<<<@<  @<<<<<< @<<<<<<  @<<<<< @<<<<  @<<<<< @<<<<  @<<<<<<  @<<<<<< @<<<<<<  @<<<<<< @<<<<<<
$run,  $type, $soleno,$dipole, $uppol, $updpol,$dwnpol,$dwndpol,$polInfo,$nmr,  $nmrInfo,$sm1,  $sm2
.

########## FETCHING RUN NUMBERS BETWEEN minrun AND maxrun ##########
########## RETRIEVE SM1/SM2 CURRENTS ##########
my $query = "SELECT runnb,runtypeid,sm1,sm2 FROM tb_run LEFT JOIN tb_beamvalues USING ( bvalueid ) WHERE runnb>=$strtrun and runnb<=$maxrun;";
my $sth = $dbh->prepare($query); $sth->execute();

while ( my @row = $sth->fetchrow() ) {

    ############### LOOPING on ALL RUNS PREVIOUSLY FOUND ###############

    $run = $row[0]; $type = $row[1]; $sm1 = $row[2]; $sm2 = $row[3];

    my $isPhysics = (1<<$type-1)&$physicsTypes;
    if ($run>=$minrun) {
	if    ($types==-2)                 { $types = $type; }
	elsif ($types>=0 && $type!=$types) { $types = -1; }
	$nRuns++; $sSM1 += $sm1; $sSM2 += $sm2; $nSMs++;
	if ($isPhysics) {
	    $NRuns++; $SSM1 += $sm1; $SSM2 += $sm2; $NSMs++;
	    if    ($Types==-2)                 { $Types = $type; }
	    elsif ($Types>=0 && $type!=$Types) { $Types = -1; }
	}
    }

    ##### RETRIEVE POLAR+dPOLAR FROM "tb_offlinepolar" WITHIN...#####
    # ...    [start,stop]
    #     or [start,start+30'] if no stop
    $query = "SELECT Ups_D_polar,Error_Ups_D_polar,Dwn_D_polar,Error_Dwn_D_polar,Solenoid_current  FROM tb_offlinepolar,tb_run WHERE Timing>starttime and Timing<IFNULL(stoptime, starttime + INTERVAL 30 MINUTE) && runnb=$run;";
    my $sth2 = $dbh->prepare($query); $sth2->execute();
    #####                         AVERAGING ON ALL POLAR'S VALID FOR CURRENT RUN
    my $Suppol  = 0; my $Supdpol2  = 0;
    my $Sdwnpol = 0; my $Sdwndpol2 = 0;
    my $Ssoleno = 0;
    my $nmeas   = 0;
    while ( my @row2 = $sth2->fetchrow() ) {
	# print "$row[0]\t @row2\n";
	$Suppol  += $row2[0]; $Supdpol2  += $row2[1]*$row2[1];
	$Sdwnpol += $row2[2]; $Sdwndpol2 += $row2[3]*$row2[3];
	$Ssoleno += $row2[4];
	$nmeas++;
    }
    if ($nmeas) {
	$uppol  = $Suppol/$nmeas;  $updpol  = sqrt($Supdpol2)/$nmeas;
	$dwnpol = $Sdwnpol/$nmeas; $dwndpol = sqrt($Sdwndpol2)/$nmeas;
	$soleno = $Ssoleno/$nmeas;
	$polRef = $run;
	if ($run>=$minrun) { # Averaging over all runs
	    $sUppol  += $Suppol;  $sUpdpol2  += $Supdpol2; $nPols += $nmeas;
	    $sDwnpol += $Sdwnpol; $sDwndpol2 += $Sdwndpol2;
	    if ($isPhysics) { # Averaging over physics runs
		$SUppol  += $Suppol;  $SUpdpol2  += $Supdpol2; $NPols += $nmeas;
		$SDwnpol += $Sdwnpol; $SDwndpol2 += $Sdwndpol2;
	    }
	}
    }
    $dipole = 1.e9; # No info about target dipole yet
    if ($polRef==$run) {
	$polInfo = 0;
    }
    else {
	$polInfo = $polRef;
    }

    ##### RETRIEVE NMR FROM "tb_NMRSM2" WITHIN...#####
    # ...    [start,stop]
    #     or [start,start+30'] if no stop
    $query = "SELECT nmr FROM tb_NMRSM2,tb_run WHERE date>starttime and date<IFNULL(stoptime, starttime + INTERVAL 30 MINUTE) && runnb=$row[0];";
    $sth2 = $dbh->prepare($query); $sth2->execute();
    my $Snmr = 0; $nmeas = 0;
    while ( my @row3 = $sth2->fetchrow() ) {
	# print "$row[0]\t @row3\n";
	$Snmr += $row3[0]; $nmeas++;
    }
    #####                           AVERAGING ON ALL NMR'S VALID FOR CURRENT RUN
    if ($nmeas) {
	$nmr = $Snmr/$nmeas; $nmrRef = $run;
	if ($run>=$minrun) { # Averaging over all runs
	    $sNMR += $Snmr; $nNMRs += $nmeas;
	    if ($isPhysics) { # Averaging over physics runs
		$SNMR += $Snmr; $NNMRs += $nmeas;
	    }
	}
    }
    if ($nmrRef==$run) {
	$nmrInfo = 0;
    }
    else {
	$nmrInfo = $nmrRef;
    }

    if ($run>=33000) { ##### IF LATER OR EQUAL 2004
	##### ...RETRIEVE TARGET CURRENTS FROM "tb_tgtcurrent" WITHIN...#####
	# ...    [start,stop]
	#     or [start,start+30'] if no stop
	#Excerpt from "$CORAL/src/condb/ysqldb/MyInterface/MySQLDBInterface.cc":
	#//from 2006 the solencur column is signless, the good sign is to be taken from the SOLN-DCCT_ai column
	$query = "SELECT solencur,dipolcur,`SOLN-DCCT_ai` FROM tb_tgtcurrent,tb_run WHERE mestime>starttime and mestime<IFNULL(stoptime, starttime + INTERVAL 30 MINUTE) && runnb=$run;";
	$sth2 = $dbh->prepare($query); $sth2->execute();
	$soleno = 1.e9; $dipole = 1.e9;
	my $Ssoleno = 0; my $Sdipole = 0; $nmeas = 0;
	while (my @row4 = $sth2->fetchrow() ) {
	    my $solSign = 1;
	    if ( $#row4>2 ) { # Meaning there was indeed a "SOLN-DCCT_ai" info.
		if ($row4[2]<0) { $solSign = -1; }
	    }
	    $Ssoleno += $row4[0]*$solSign;
	    $Sdipole += $row4[1]; $nmeas++;
	}
	#####             AVERAGING ON ALL TARGET CURRENTS VALID FOR CURRENT RUN
	if ($nmeas) {
	    $soleno = $Ssoleno/$nmeas; $dipole = $Sdipole/$nmeas;
	    if ($run>=$minrun) { # Averaging over all runs
		$sSoleno += $Ssoleno; $sDipole += $Sdipole; $nTargets += $nmeas;
		if ($isPhysics) { # Averaging over physics runs
		    $SSoleno += $Ssoleno; $SDipole += $Sdipole; $NTargets += $nmeas;
		}
	    }
	}
    }
    elsif ($nmeas!=0) {
	if ($run>=$minrun) { # Averaging over all runs
	    $sSoleno += $Ssoleno; $nTargets += $nmeas;
	    if ($isPhysics) { # Averaging over physics runs
		$SSoleno += $Ssoleno; $NTargets += $nmeas;
	    }
	}
    }
    

    if ($run>=$minrun) { # Write only those runs w/in the requested range
	write;
    }
}

printf "-------------------------------------------------\n";

#### PRINT AVERAGES on ALL RUNS
$run = $nRuns; $type = $types;
if ($nSMs!=0)  { $sm1 = $sSM1/$nSMs; $sm2 = $sSM2/$nSMs; }
else           { $sm1 = 0;           $sm2 = 0; }
if ($nNMRs!=0) { $nmr = $sNMR/$nNMRs; }
else           { $nmr = 0; }
if ($nTargets!=0) {
    $soleno = $sSoleno/$nTargets; $dipole = $sDipole/$nTargets;
}
else {
    $soleno = 0; $dipole = 0;
}
if ($nPols!=0) {
    $uppol  = $sUppol/$nPols;  $updpol  = sqrt($sUpdpol2)/$nPols;
    $dwnpol = $sDwnpol/$nPols; $dwndpol = sqrt($sDwndpol2)/$nPols;
}
else {
    $uppol = 0; $dwnpol = 0; $updpol = 0; $dwndpol = 0;
}
write;

##### PRINT AVERAGES on PHYSICS RUNS
if ($NRuns!=$nRuns) {
    $run = $NRuns; $type = $Types;
    if ($NSMs!=0)  { $sm1 = $SSM1/$NSMs; $sm2 = $SSM2/$NSMs; }
    else           { $sm1 = 0;           $sm2 = 0; }
    if ($NNMRs!=0) { $nmr = $SNMR/$NNMRs; }
    else           { $nmr = 0; }
    if ($NTargets!=0) {
	$soleno = $SSoleno/$NTargets; $dipole = $SDipole/$NTargets;
    }
    else {
	$soleno = 0; $dipole = 0;
    }
    if ($NPols!=0) {
	$uppol  = $SUppol/$NPols;  $updpol  = sqrt($SUpdpol2)/$NPols;
	$dwnpol = $SDwnpol/$NPols; $dwndpol = sqrt($SDwndpol2)/$NPols;
    }
    else {
	$uppol = 0; $dwnpol = 0; $updpol = 0; $dwndpol = 0;
    }
    write;
}
