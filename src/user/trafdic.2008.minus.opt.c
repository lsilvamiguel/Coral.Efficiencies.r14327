C2OPT_HEADER
#if 0

// $Id: trafdic.2008.minus.opt.c 13326 2012-03-23 17:57:53Z ybedfer $

// coral options c file for TraFDic on 2008 data
//
// NOTA.BENE.: c file is to be processed with:
//
//        "gcc -E -P -C -D<defines>"
//
// before being supplied to coral. 

// Originally copied from "./trafdic.2008.opt.c,v1.33"

#endif
#if !defined C2OPT_COMPUTERFARM
#  define C2OPT_COMPUTERFARM 1
#endif
#if   C2OPT_COMPUTERFARM == 1
// TraFDic on 2008 <0 hadron data (i.e. run# <= 71018) @ CERN.
#elif C2OPT_COMPUTERFARM == 2
// TraFDic on 2008 <0 hadron data (i.e. run# <= 71018) @ gridKa.
#elif C2OPT_COMPUTERFARM == 3
// TraFDic on 2008 <0 hadron data (i.e. run# <= 71018) @ Lyon.
#elif C2OPT_COMPUTERFARM == 4
// TraFDic on 2008 <0 hadron data (i.e. run# <= 71018) @ Compass ONLINE.
#else
// TraFDic on 2008 <0 hadron data (i.e. run# <= 71018).
#endif

#ifdef C2OPT_BARE_FILE
//  ``Bare'' options file to be included in "$CORAL/src/alignment/traf.*.opt"
// files for processing 2008 FIELD ON data for alignment purposes:
//  - It lacks input data specification (this must be entered in the main
//   "traf.*.opt" file)...
//  - ...as well as mDST output, all high level (vertex and the like)
//   reconstruction options, etc..
//  - More generally, this file should not be edited: any modification to the
//   option entries it contains must instead be implemented in the main
//   "traf.*.opt" by overwriting.
#else
//  Restricted to magnets ON data.
// (Cf. also dedicated file for proton data)

//  This file sums up the main options. Other options are specified elsewhere
// and referenced here by "include <path_to_2ndary_options_file>" statements.
// All options can be entered several times, whether in main or secondary files.
// The last one, in the stream of entries and "include" statements, always wins.

//  This file can serve as the basis for the template of mass production (
// i.e. so-called "template.opt" file). It has to be adapted though, e.g. in
// order to specify placeholders for the IO. Cf. in fine the block entitled:
//         "TEMPLATE FOR MASS PRODUCTION".

//  For any other application, it may have to be edited or, better, ``included
// and overwritten''. Cf. the TWiki page:
// wwwcompass.cern.ch/twiki/bin/view/DataReconstruction/CoralSoftware#Options_files

//  !!!NOTA BENE!!!: Working aroud RUN-TO-RUN MISALIGNMENT, cf. dedicated block
// of options infra, is disabled by default, since several special run-by-run
// alignment subdirectories are available for a genuine fix, in the official
// 2008 geometry directory. But neither of these covers the all of 2008.
//  => Have to check whether data about to be precessed are indeed covered. If
//    not, the work-around block has to be enabled.
#endif

#ifndef C2OPT_BARE_FILE
//                                         ***** PHAST OUTPUT
// - Active only if used on a PHAST-compatible executable, cf. "./README"
// - Cf. "$PHAST/coral/README" for the details of the options
mDST file		${HOME}/phast.root
// Information reduction bit pattern: "mDST select".
mDST select		9	// At least one vertex or particle | First and last helices only.
//mDST hits		HL	// MegaDST...
//mDST Digits		EC	// GigaDST... (raw information after calibrations are applied, according to documentation only useful for non-tracking detectors, use mDST hits for tracking detectors)
//mDST DAQdigits	SC	// GigaDST... (raw information as read from raw data stream, no calibrations applied)

//					// ***** WRITE-BACK of RAW DATA
//				<File name>	<Type>
//CsKalmanFitting WriteBack	evtDump.raw	Physics		// As of 08/02: excl. 5/7 pions, K0, Lambda, J/psi(?)
//CsEvent         WriteBack	evtDump.raw	Randoms		// <File name> may differ from CsKalmanFitting's
#endif

//					   ***** HISTOGRAMMING
histograms package	ROOT
histograms home		${HOME}/trafdic.root

#ifdef C2OPT_BARE_FILE
//					   ***** THIS IS A DATA RUN...
//					   ***** ...W/O ANY INPUT SPECIFIED!
#else
//					   ***** THIS IS A DATA RUN
Data job
#endif
hadron run					// ...of the hadron type
#ifndef C2OPT_BARE_FILE
//					   ***** INPUT DATA...
#  if   C2OPT_COMPUTERFARM == 1
Data file	/castor/cern.ch/compass/data/2008/raw/W37/cdr24011-70460.raw
#  elif C2OPT_COMPUTERFARM == 2
Data file 	/grid/fzk.de/compass/compass2/data/raw/W24/cdr13009-66806.raw
#  elif C2OPT_COMPUTERFARM == 3
Data file	cchpsscompass:/hpss/in2p3.fr/group/compass/data/2008/raw/W37/cdr24011-70460.raw
#  elif C2OPT_COMPUTERFARM >= 4
Data file	/castor/cern.ch/compass/data/2008/raw/W37/cdr24011-70460.raw
#  endif
#endif


// 					   ***** DETECTOR TABLE...
//  A full directory (path name terminated by a slash "/") or the full path to
// a particular file can be specified.
//  Here, one of the RUN-BY-RUN subdirectories is specified. An alternative
// being suggested, commented out by default.
//  If neither of these RUN-BY-RUN sub-dir's are retained, one has to enable the
// RUN-TO-RUN MISALIGNMENT work-around block infra.
detector table  	$COMPASS_FILES/geometry/2008/W37/  	// RUN-BY-RUN sub-directory
//detector table  	$COMPASS_FILES/geometry/2008/W37b/  	// Alternative RUN-BY-RUN sub-dir


// 					   ***** ROOTGeometry
CsROOTGeometry file	$COMPASS_FILES/geometry/2008/ROOTGeometry/detectors.lh2.r354.C	// LH2: this the file used in 08W37/slot4 mass prod.
CsROOTGeometry massDefault	.1396	// for muons use .1057 instead
TraF	ELossStraggling 	1


//					   ***** CALIBRATION
#if C2OPT_COMPUTERFARM > 4
use calibration
CDB use  	FileDB
CDB location  	$COMPASS_FILES/calibrations/2008
#else
use calibration
CDB use 	MySQLDB
#  if   C2OPT_COMPUTERFARM == 1
CDB server	wwwcompass	// This is, as of 07/05, pointing to "compass02" located in computer center (reachable outside Cern)
#  elif C2OPT_COMPUTERFARM == 2
CDB server compass.gridka.de	// Located at gridKa (internal network)
CDB specialplace GRIDKA 	// Location of calibration files at gridKa
#  elif C2OPT_COMPUTERFARM == 3
CDB server	cccompassdb	// Located at Lyon   (internal network)
CDB portnumber  23306		// Port to be used in place of the std one
CDB specialplace LYON	 	// Location of calibration files at Lyon
#  elif C2OPT_COMPUTERFARM == 4
CDB server	pccodb00
CDB specialplace DAQ	 	// Location of calibration files in ONLINE
#  endif
// "CDB entrytime Y-M-D-H:M:S" can be entered to reject all calibrations issued
// after specified date. CDB being constantly evolving, leaving "CDB entrytime"
// disabled may give rise to unpredictable behaviour.
//CDB	entrytime	2007-12-01-00:00:00
// "entrytime" entry can be made to affect only a (set of) specific TBName(s):
// "TBname(or TB) CDBentrytime(sic!)"
//DW	CDBentrytime	2007-12-01-00:00:00	// Disregarding the calibs entered on Dec. 14th, 2007
// CDB can be replaced by a File DB, for any given (set of) TBname(s):
#  if   C2OPT_COMPUTERFARM == 1
//DW FileDB	/afs/cern.ch/user/l/leberig/public/marcin/
#  else
//DW FileDB	$COMPASS_FILES/fileDBs/DW/marcin
#  endif
#endif

//					   ***** MAPPING
decoding map	$COMPASS_FILES/maps/2008.xml


//					   ***** EVENT SELECTION
events to read 30000			// # of events to read
//events to skip 130//23014-70499#32526749 //2328//25003-69853#7367980 //576//25003-69853#7346956 //4 //143//25003-69853#7341760 //119//25003-69853#7341472
// BoS Veto. 
//events BOS_skip 1.151	// Status of hardware veto unknown as of 08/06 => Keep software veto commented out.

//selection trigger mask f		// Trigger selection (hexadecimal)
//selection trigger strict


// seed for random number generations:
random number engine JamesEngine
random number seed 19990102
//reset random seed every new event

DC 	make always two clusters	// Even if they coincide, so that their re-evaluation w/ an event time differing from the trigger might yield 2 distinct clusters
DR 	make always two clusters
DW 	make always two clusters
MC 	make always two clusters
ST 	make always two clusters

pattern method 1        // not used yet
tracking method 1       // not used yet

// 					   ***** DECODING
make decoding 	 	// <nothing>, MCExact

//  	 	 	 	 	   ***** CLUSTERING
make clustering 	// <nothing>, MCExact, MCSmeared, MCQuantized

//					   ***** RECONSTRUCTION SCHEMA
reconstruction schema 1

//					   ***** TRACKING
make tracking
track prepattern method traffic
track bridging   method traffic
track fitting    method traffic

#ifndef C2OPT_BARE_FILE
//					   ***** Beam RECONSTRUCTION
//make beam reconstruction
//include ${CORAL}/src/pkopt/beam_2004.opt
//                   Overwrite what's set in "beam_2004.opt"
//BeamRecons	useTRAFFIC	1	// >0: Traffic is used, =0: standalone beam package is used
//BeamRecons	AllowMultCand	2
//                   Beam PID
PID_doBeamID	selAllBeams	1


//					   ***** Vertex RECONSTRUCTION
make	vertex	reconstruction
vertex pattern method averaging
vertex fitting method kalman
include ${CORAL}/src/pkopt/vertex.2002.opt
//                   Overwrite what's set in "vertex.opt"
CsAverPattern  Hist	1	// 0 - Off, 1 - ON
CsAverPattern  findSec	1	// 0 - Off, 1 - ON
CsAverPattern  Print [ 0 ]	0	// Prefilter info.
CsAverPattern  Print [ 1 ]	0	// Error info.
// Track is considered as accepted if its momentum is smaller
// than "AcceptTill" percents of beam track
CsAverPattern	AcceptTill	105	// in %
CsAverPattern	TimePrimCut	10 	// in sigmas
// 	 	 	Z Dist 	DCA(number of sigmas)
CsAverPattern CUTS 	2550 	10 	// Z cut so that vertex search domain includes GP01
CsKalmanFitting Specials	0	// 0 - Off, 1 -  mu' takes precedence
CsAverPattern	Refit		1	// 0 - Off, 1 - ON
CsAverPattern	Retrack		2.5	// 0 - Off, !0 = cut on event time
CsKalmanFitting	Hist		1	// 0 - Off, 1 - ON
CsKalmanFitting	Print [ 0 ]	0	// Kalman filter info.
CsKalmanFitting	Print [ 1 ]	0	// Global fit info.
CsKalmanFitting	Print [ 2 ]	0	// IKF info.
CsKalmanFitting	RefPlane	500	// in mm
CsKalmanFitting	RefBeam	-1500	// in mm
CsKalmanFitting	RefMargin	500	// in mm
// Chi2 increment cuts:		DirChi2Trk	InvChi2Trk
CsKalmanFitting	CUTS		7		17
// Vertex Chi2 cut:
CsKalmanFitting	Chi2VCut	10.5
CsKalmanFitting BeamP0		190	// Beam reference momentum (for histo purposes)

//					   ***** RICH1
make	rich1	reconstruction
include ${CORAL}/src/pkopt/rich1.2008.opt
CsRICH1UpGrade  CONFIG		// RICH upgrade: 12 APV 4 MAPMT
RICHONE 	PrintDetMirr	0
RICHONE 	AcknoMethods	0	// prints methods in use
RICHONE 	PrintKeys	0	// 1 = print keys read,  0 = No print
RICHONE 	PrintConsts	0	// 1 = print rec. constants read,  0 = No print
RICHONE 	DoThetaLikeMax 	NO 	// Faster but "the output buffer of rich won't have the angle that maximizes the likelihood anymore" (according to Giulia, cf. mail from 2009/06/29)

//					   ***** CALORIMETERS
make	calorimeters	reconstruction
// To select the reconstruction algorithm ot be used change the file included,
// or set appropriate MoreRecoOptions lines below
include ${CORAL}/src/pkopt/calorim_2008.opt
// if for the period you are looking at no LED/Laser corrections are available,
// uncomment the next two lines
//EC01P1__ MoreRecoOptions USE_LED_REF_CALIB_CORRECTION=NO USE_LED_REF_CALIB_CORRECTION_IN_SPILLS=NO
//EC02P1__ MoreRecoOptions USE_LED_REF_CALIB_CORRECTION=NO USE_LED_REF_CALIB_CORRECTION_IN_SPILLS=NO

//CsBuildPart	Hist		5	// Calo histos. Disabled by default, because too big (several MB). But programmed to be enabled in mass production: cf. in fine.

//					   ***** RICH WALL + CALORIMETER
make	RW	charge reconstruction	// RW+calo charge reco
include	${CORAL}/src/pkopt/rw_ch_2008.opt

//					   ***** FLUX SCALERS
include ${CORAL}/src/pkopt/trigger.2008.opt
#else
//					   ***** NO Beam RECONSTRUCTION

//					   ***** NO Vertex RECONSTRUCTION

//					   ***** NO RICH1

//					   ***** NO CALORIMETERS

//					   ***** NO MU' PID
#endif

//  	 	 	 	 	   ***** GEOMETRICAL ZONES
define zone      -350    3500   before M1
define zone      3500   21000   between M1 and M2
define zone     21000   32800   between M2 and Muon Wall
define zone     32800   99999   after Muon Wall

define zone     -8000    -350   before the target


//  	  	  	  	  	   ***** MAGNETIC FIELD MAPS
CsField SM1m_field_measured $COMPASS_FILES/maps/mag_fields/SM1m/SM1M.map.172.data
CsField SM2_field	 $COMPASS_FILES/maps/mag_fields/SM2/FSM.map.5000.data

//					   ***** MAGNETs INFO
#if C2OPT_COMPUTERFARM <= 4
CsMagInfo	MySQLDB
#else
CsMagInfo	File	/afs/cern.ch/user/y/ybedfer/public/magInfo.2008.txt	// ...if MySQL not avail
#endif
CsMagInfo	SM2	1 	// Get NMR from MySQLDB and re-scale SM2 w/ NMR / NMR_map * correction factor. The correction factor is here set =1. This yields, for an average run of 2008/W37, a re-scaling of 1.0042, i.e. slightly higher than the 1.0041 that's typically specified in detectors.dat's (and was derived I don't know how).



//					   ***** LOGGER
error logger log level		error	// debugging, verbose, info, anomaly, warning, error, fatal
error logger store level	none
error logger verbosity		normal	// low, normal, high. 

//	  	  	  	  	   ***** TraFDic (+DET OPTIONS+...)

include ${CORAL}/src/pkopt/trafdic.2008.opt

//	********** OVERWRITE WHAT'S SET IN "../pkopt/trafdic.2008.opt" *********

// ***** WORK-AROUND for RUN-TO-RUN MISALIGNMENT
// (To be enabled when RUN-BY-RUN alignment is not available, or for some
// reason, not retained.)
// - Detector table:
//detector table  	$COMPASS_FILES/geometry/2008/	// General purpose alignment directory
// ==> ==> Implications:
// - RUN-TO-RUN misalignment can be worked around, somewhat, by artificiallly
//  blowing up SI routes and uncertainties:
//TraF	dCut	[84]	.0050	// SI route enlargement w.r.t. resolution
//TraF	dCut	[85]	.0040	// SI position uncertainties correction term

// ***** HADRON SPECIFICS in the track reconstruction scheme:
// - Bridging over target:
Traf	ReMode	[26]	1	// Always try bridging over target
// - Enable double bridging:
TraF	ReMode	[18]	15	// Backtrack. &0x1: Zone 1 (&0x2: even w/out P). &0x4: Zone 2, &0x8: Double Bridging
// - Fast fit in Bridging conditioned by RICH pipe:
TraF	ReMode	[46]	1	// Track fit in Bridging: fast fit (except for tracks through RICH pipe) terminated one-pass full KF. 
TraF	ReMode	[14]	15	// Dico fit: 0xf = in PR and in Bridging (for init'ing P)
TraF	dCut	[97]	25.	// Preliminary chi2/Nhits for track candidate bridged w/ FullKF
TraF	dCut	[ 9]	25.	// Final chi2/Nhits cut for bridged track candidate
// - Arbitration between 0x3 and 0x6 bridges
TraF	dCut	[94]	0	// Disabled
// - SI space point cleaning
TraF	iCut	[36]	4	// Check consistency of SI space points w/ multiplicity >= 4.
//TraF	ReMode	[47]	3	// CPU saving in Bridging: 0x1: Recycle FullKF fit, 0x2: No systematic refit after cleaning.

// Redefine the TBNAMES of DETECTORS to be EXCLUDED FROM PATTERN RECOGNITION:
// (Caveat: one has to supply the full list of excluded detectors each time one re-specifies this "TraF DetNameOff" option.)
// - In any case, some detectors ("VI/VO/DW02X2") have to be left turned off , cf. explanations in "../pkopt/trafdic.2008.opt".
// - But one may want/need to turn off some other detector...
// - E.g. "DW02Y1" broke one of its wire end of August 2008, probably between 08W35#70113 and 08W35#70116, cf. "wwwcompass.cern.ch//runLogbook/dirphp/show_comment.php?commnb=24375" => It has to be turned off depending upon the period (>=08W35) and the mass production of the real data (check the options file used therein) one wants to emulate, as is done in the following, commented out, entry.
//TraF	DetNameOff	VO	VI	DW02X2	DW02Y1	HI05
// ==> ==> Implications:
// - DETECTORS being EXCLUDED may imply updating some PR options, e.g.:
//TraF	ReMode	[17]	0	// Disable GEM amplitude correlations, if any GEM is among the turned-off detectors.
LR associate  DW02Y1__   DW02Y2__ 45.0	  0.5 	0	// To be disabled if "DW02Y1" is turned off

//TraF	Graph	[0]	0	// coral's graphics switches off automatically under BATCH @ CERN, gridKa and Lyon. This option forces it to do so in any case.
//TraF	iCut	[0]	7	// Trigger selection in TraFDic
//TraF	Print	[0]	3	// Event-by-event numbering
//TraF	Print	[8]	3	// EoE

#ifndef C2OPT_BARE_FILE
//			SMOOTHING POINTS
// The points @ RICH (resp. Calos) are to used by RICH (resp. Calo) softwares.
TraF	SmoothPos [0]	745.		// RICH
// From: Vladimir Kolosov
// Sent: Thu 2/28/2008 9:14 AM
// Let's decide to calculate track position 10 cm upstream front surface of a calorimeter.
// calo  EC01P1__  1399.500 - 45./2       - 10. = 1367.0
// calo  EC02P1__  3325.200 - 45./2       - 10. = 3292.7
// calo  HC01P1__  1573.700 -100./2 -22.5 - 10. = 1501.2
// calo  HC02P1__  3580.000 -120./2       - 10. = 3510.0
// 2011-04-28: updated to match current detectors.dat
TraF	SmoothPos [1]	1367.0	// Smooth track at ECAL1
TraF	SmoothPos [2]	1501.2	// Smooth track at HCAL1
TraF	SmoothPos [3]	3292.7	// Smooth track at ECAL2
TraF	SmoothPos [4]	3510.0	// Smooth track at HCAL2
#endif

//                 STRAW options: X-ray and signal propagation correction
//  X-ray can only be enabled if alignment (cf. "detector table" entry supra)
// did also take X-ray correction into account.
//  And in any case, is incompatible w/ "VarPitch", cf. infra.
// T0 per channel calib are available for 2009, cf. mail from J-F.:
// ======================================================================
// On Mon, 2009-06-22 at 16:45 +0200, Jean-Francois Rajotte wrote:
// I have calibrated the T0s for 2009.  I have put the calibration to the
// database and can be switched on with this line in coral opt file :
// 
// STRAW settings spacers=YES signal_propagation=YES T0=YES
// 
// The calibration changes the T0 by 3-4 ns on average, [...]
// ======================================================================
// I checked T0=YES improves average track chi2 by ~2.5% on 08W30:25003-69853
// NOTA BENE: A unique entry for all "STRAW settings".
STRAW settings spacers=YES signal_propagation=YES T0=YES

//                 Enabling variable pitch
//  This can be enabled, for any given detector, only if alignment (cf.
// "detector table" entry supra) takes VarPitch correction into account.
//  And in any case, as far as straws are concerned, is incompatible w/ X-ray
// correction, cf. supra.
//VarPitch    DetOffVarP	DC	ST	DW	FI	GM
//VarPitch	table	~corallib/public/PitchTable/PitchTable_FIGMDCSTDW_2003.dat

#if   C2OPT_FLAG == 1
// ********** SPECIAL SETTING TO HISTOGRAM GEMs Amplitude Correlation **********
// (in "TEv::RDMonitor", which need be recompiled w/ "RDMon_DIGIT_DATA" and
// "RDMon_DIGIT_GM", and, in a first iteration, "RDMon_GM_DEFAULT" defined)
TraF	Hist	[ 1]	17	// Require high level histo in TraFDic
TraF	Hist	[16]	0	// Turn off histo'ing active detectors 
TraF	Hist	[17]	14	// OR of zones histo'ed = 0xe
TraF	Hist	[18]	2	// Tracks used: zones &= 0x2
// (in "TEv::Quadruples", which need be recompiled w/ "Quadruples_HISTO" defined)
TraF	Hist	[ 6]	1	// In PrePattern et al. (FindSpace, Quadruples)
#elif C2OPT_FLAG == 2
// ********** SPECIAL SETTING TO HISTOGRAM RESIDUALS IN HIs AND MBs **********
TraF	Hist	[ 1]	17	// Monitor, &0x10:Residuals
TraF	Hist	[17]	12	// Residuals: |Groups of inactive det's histo'ed
TraF	Hist	[18]	14	// Residuals: &Groups of tracks histo'ed
TraF	DetNameOff	VO	VI	BM	SI03	H	MB
TraF	SmoothPos [1]	2600.	// Smoothing point, to extrapolate from, when in zone 0x4
#elif C2OPT_FLAG == 3 || C2OPT_FLAG == 4
// ********** SPECIAL SETTING TO TRACK OFF-TIME MUONS **********
// ********** AND HISTOGRAM their residuals in FI01:2s and SIs **********
// *****  I) Tracks reconstructed in spectrometer
// ***** II) Tracks reconstructed in both spectro and scifi/Si telescope,...
//          ...bridged over target, refit with the sole FI02s and SI01XYs.
FI01X1__ HitTime [0-1] 50 200
FI01Y1__ HitTime [0-1] 50 200
FI02X1__ HitTime [0-1] 50 200
FI02Y1__ HitTime [0-1] 50 200
FI03X1__ HitTime [0-1] 50 200
FI03Y1__ HitTime [0-1] 50 200
FI03U1__ HitTime [0-1] 50 200
FI04X1__ HitTime [0-1] 50 200
FI04Y1__ HitTime [0-1] 50 200
FI04U1__ HitTime [0-1] 50 200
FI05X1__ HitTime [0-1] 100 400
FI05Y1__ HitTime [0-1] 100 400
FI06X1__ HitTime [0-1] 100 400
FI06Y1__ HitTime [0-1] 100 400
FI06V1__ HitTime [0-1] 100 400
FI07X1__ HitTime [0-1] 100 400
FI07Y1__ HitTime [0-1] 100 400
FI08X1__ HitTime [0-1] 100 400
FI08Y1__ HitTime [0-1] 100 400
#  if C2OPT_FLAG == 3
// *****  I) Tracks reconstructed in spectrometer
TraF	Hist	[ 1]	17	// Monitor, &0x10:Residuals
TraF	Hist	[16]	16	// Residuals: |Groups of det's  histo'ed
TraF	Hist	[17]	0	// Residuals: |Groups of inactive det's histo'ed
TraF	Hist	[18]	7	// Residuals: &Groups of tracks histo'ed
TraF	iPRpar [  0 -  5 ]	2   3   6	2  3  6
TraF	iPRpar [ 10 - 15 ]	2   3   5	2  3  5
TraF	iPRpar [ 20 - 25 ]	2   2   4	2  2  4
TraF	iPRpar [ 30 - 35 ]	3   1   9	3  1  7
TraF	iPRpar [ 40 - 45 ]	2   4  11	3  2 10		// W/out FI15
TraF	DetNameOff	V	M	D	ST	P	H	GM01	GM02	GM03	GM04	GM05	GM06
#  elif C2OPT_FLAG == 4
// ***** II) Tracks reconstructed in both spectro and scifi/Si telescope,...
TraF	ReMode	[26]	1	// != 0 - Bridging over target
TraF	Hist	[ 1]	17	// Monitor, &0x10:Residuals
TraF	Hist	[16]	16	// Residuals: |Groups of det's  histo'ed
TraF	Hist	[17]	0	// Residuals: |Groups of inactive det's histo'ed
TraF	Hist	[18]	23	// Residuals: &Groups of tracks histo'ed
TraF	DetNameOff	V	M	D	ST	P	H	GM01	GM02	GM03	GM04	GM05	GM06	SI03
TraF	Det2Go2Fit	FI02	SI02X	SI02Y	FI03	FI04	FI05	FI06	FI07	FI08	GM
#  endif
#endif

//  ===== TEMPLATE FOR MASS PRODUCTION (uncomment, and modify if needed) =====

//					***** PLACEHOLDERS for I/O
// (Note: One has to uncomment next 6 entries.)
//events to read 200000
//Data		file	TMPRAWFILE
//mDST		file	TMPMDSTFILE
//histograms	home	TMPHISTFILE
//CsKalmanFitting WriteBack	testevtdump.raw	Physics
//TMPRICHFILE


//					***** mDST OPTIONS
// (Note: One has to systematically uncomment next 3 entries. For special MegaDST production, one has to uncomment also the 4th one. I don't know what to recommend for the 5th one.)
//mDST	select	0	// 0: no selection
//mDST	DAQdigits	RP CE
//mDST	selectTrigger	3072	// Random trigger: unconditional output (Note: the option is not useful when "mDST select 0", i.e. all-out unconditional output.)

//mDST	hits	ALL		// SPECIAL OPTION! for the production of MegaDST
//mDST	Digits	EC


//					***** CDB OPTIONS
// (Note: One has to freeze the CDB info accessible to the production, so that the latter be independent of CDB later developments, and hence stable. A reasonable choice for the "entrytime" is the date of the start of production. One may still want to be more restrictive on all of the CDB. Or on the part of it concerning a particular detector or set of detectors: in that case, use the option "<generic name> CDBentrytime" (cf. examples in the body of this options file).)
//CDB	entrytime	2009-2-17-00:00:00	// TO BE MODIFIED!


//					***** DETECTOR MONITORING
// (Note: Have to uncomment next 4 entries.)
//Monitoring	Directory	ON	// Paolo's data monitor graphs.
//Monitoring	Residuals	ON	// Albert's monitoring histos.
//Monitoring	Efficiency	ON
//CsBuildPart	Hist		5	// Calorimeter


//					***** SETUP FILES
// (Note: One has to freeze the setup info accessible to the production, so that the latter be independent of the modifications later brought to the setup files, and hence stable. This is achieved by copying the setup files to sub-directories of the production directory. E.g. infra are listed what need be done for the production of 2008W37/slot2.) 
//decoding	map	/afs/cern.ch/compass/scratch/d17/objsrvvy/2008/W37/slot2/maps
//detector	table	/afs/cern.ch/compass/scratch/d17/objsrvvy/2008/W37/slot2/alignement/
//TraF		Dicofit	/afs/cern.ch/compass/scratch/d17/objsrvvy/2008/W37/slot2/dico/
//CsROOTGeometry file	/afs/cern.ch/compass/scratch/d17/objsrvvy/2008/W37/slot2/ROOTGeometry/detectors.lh2.r???.C

//					***** RECONSTRUCTION OPTION
// (Note: Depending upon the choices made supra, some alternative reconstuction option may be activated.) 
// - RUN-TO-RUN MISALIGNMENT has to be worked around by artificiallly blowing up SI routes and uncertainties, if there's no RUN-BY-RUN alignment available:
//TraF	dCut	[84]	.0050	// SI route enlargement w.r.t. resolution
//TraF	dCut	[85]	.0040	// SI position uncertainties correction term

// ===========================================================================
//end
