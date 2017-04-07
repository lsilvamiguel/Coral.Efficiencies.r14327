C2OPT_HEADER
#if 0

// $Id: trafdic.2015.opt.c 14094 2015-11-06 15:28:48Z lsilva $

// coral options c file for TraFDic on 2015 (DY) data
//
// NOTA.BENE.:c file is to be processed with:
//
//        "gcc -E -P -C -D<defines>"
//
// before being supplied to coral. 

// Derived from "./trafdic.2014.opt.c,v13981"

#endif
#if !defined C2OPT_COMPUTERFARM
#  define C2OPT_COMPUTERFARM 1
#endif
#if   C2OPT_COMPUTERFARM == 1
// TraFDic on 2015 DY data @ CERN.
#elif C2OPT_COMPUTERFARM == 2
// TraFDic on 2015 DY data @ gridKa.
#elif C2OPT_COMPUTERFARM == 3
// TraFDic on 2015 DY data @ Lyon.
#elif C2OPT_COMPUTERFARM == 4
// TraFDic on 2015 DY data @ Compass ONLINE.
#else
// TraFDic on 2015 DY data.
#endif
// (Default=PION, w/ traditional beamT and vertex options.
//  A set of alternative options is available, see "BEAM&VERTEX BLOCK" infra. 
//  MUON case is provided for in comments, see "MUON BLOCK" infra.
//  Special settings for "EARLIER DATA TAKING" are also available.)

#ifdef C2OPT_BARE_FILE
//  ``Bare'' options file to be included in "$CORAL/src/alignment/traf.*.opt"
// files for processing 2015 FIELD ON data for alignment purposes:
//  - It lacks input data specification (this must be entered in the main
//   "traf.*.opt" file)...
//  - ...as well as mDST output, all high level (vertex and the like)
//   reconstruction options, etc..
//  - More generally, this file should not be edited: any modification to the
//   option entries it contains must instead be implemented in the main
//   "traf.*.opt" by overwriting.
#else
//  This file sums up the main options. Other options are specified elsewhere
// and referenced here by "include <path_to_2ndary_options_file>" statements.
// All options can be entered several times, whether in main or secondary files.
// The last one, in the stream of entries and "include" statements, always wins.

//  For the preparation of a mass production, cf. the dedicated block in fine
// entitled "TEMPLATE FOR MASS PRODUCTION".
//  For any other particular application, this file has to be edited or, better,
// ``included and overwritten''. Cf. the TWiki page:
// "wwwcompass.cern.ch/twiki/bin/view/DataReconstruction/CoralSoftware#Options_files"
#endif

#ifndef C2OPT_BARE_FILE
//                                         ***** PHAST OUTPUT
// - Active only if used on a PHAST-compatible executable, cf. "./README"
// - Cf. "$PHAST/coral/README" for the details of the options
mDST file		phast.root
// - Information reduction bit pattern: "mDST select".
mDST select		0	// No selection
// - One may want to output more info, as exemplified infra:
//mDST hits		ALL	// MegaDST...
//mDST Digits		EC	// Raw information after calibrations are applied: only useful for non-tracking detectors, else "mDST hits" is good enough.
//mDST DAQdigits  	CA FI15	HO04 	// Raw information as read from raw data stream, no calibrations applied)
//mDST selectTrigger	1024	// Random trigger: unconditional output (Note: the option is not useful when "mDST select 0", i.e. all-out unconditional output.)
//mDST AllScalersTriggerMask	1024	// Write out all beam scalers for random trigger

//					// ***** WRITE-BACK of RAW DATA
//				<File name>	<Type>
//CsKalmanFitting WriteBack	evtDump.raw	Physics		// D*, excl. phi, K0, J/psi
//CsEvent         WriteBack	evtDump.raw	Randoms		// All 0xc00 trigger. <File name> may differ from CsKalmanFitting's
#endif

//					   ***** HISTOGRAMMING
histograms package	ROOT
histograms home		trafdic.root

#ifdef C2OPT_BARE_FILE
//					   ***** THIS IS A DATA RUN...
//					   ***** ...W/O ANY INPUT SPECIFIED!
#else
//					   ***** THIS IS A DATA RUN
Data job
#endif
#ifndef C2OPT_BARE_FILE
//					   ***** INPUT DATA...
#  if   C2OPT_COMPUTERFARM == 1
Data file 	/castor/cern.ch/compass/data/2015/raw/W10/cdr12022-261974.raw
#  elif C2OPT_COMPUTERFARM == 2
Data file 	/grid/fzk.de/compass/compass2/data/raw/tmp/cdr12102-257467.raw
#  elif C2OPT_COMPUTERFARM == 3
Data file	/hpss/in2p3.fr/group/compass/data/2015/raw/cdr12002-257467.raw
#  elif C2OPT_COMPUTERFARM >= 4
Data file	/castor/cern.ch/compass/data/2015/raw/tmp/cdr13022-258130.raw
#  endif
#endif


// 					   ***** DETECTOR TABLE...
//  A full directory (path name terminated by a slash "/") or the full path to
// a particular file can be specified.
detector table   	$COMPASS_FILES/geometry/2015/

// 					   ***** ROOTGeometry
CsROOTGeometry file	$COMPASS_FILES/geometry/2015/ROOTGeometry/detectors.C
CsROOTGeometry massDefault	.105658367
CsROOTGeometry simpleELoss	0
CsROOTGeometry ELossStraggling	0

//					   ***** CALIBRATION
#if C2OPT_COMPUTERFARM > 4
use calibration
CDB use  	FileDB
CDB location  	$COMPASS_FILES/calibrations/2012
#else
use calibration
CDB use 	MySQLDB
#  if   C2OPT_COMPUTERFARM == 1
CDB server	wwwcompass	// This is an alias, pointing, as of 10/11, to "compass02" located in computer center (reachable outside Cern)
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
//CDB	entrytime	2010-12-01-00:00:00
// "entrytime" entry can be made to affect only a (set of) specific TBName(s):
// "TBname(or TB) CDBentrytime(sic!)"
//FI	CDBentrytime 	2010-10-31-23:59:59	// Disregarding the FI calibs entered in 201/11: not yet finalised.
// CDB can be replaced by a File DB, for any given (set of) TBname(s):
#  if   C2OPT_COMPUTERFARM == 1
//DW FileDB	/afs/cern.ch/user/l/leberig/public/marcin/
#  else
//DW FileDB	$COMPASS_FILES/fileDBs/DW/marcin
#  endif
#endif

//					   ***** MAPPING
decoding map	$COMPASS_FILES/maps/2015.xml


//					   ***** EVENT SELECTION
events to read 500000			// # of events to read
//events to skip 7
// BoS Veto:
// - It has to be enforced software-wise to ensure data reco be uniform, if forgotten in hardware. (Although could be argued that this could be done in PHAST.)
// - It's commented out here.
//events BOS_skip	1.151
//selection trigger mask 20f		// Trigger selection (hexadecimal)
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
// (Muon data: have to enable BMS reconstruction, cf. MUON BLOCK infra)
PID_doBeamID	selAllBeams	1	// Allow beams w/o BMS in the vertexing

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
CsAverPattern	AcceptTill	105	// In %
CsAverPattern	TimePrimCut	10 	// Time diff between track and beam, in sigmas
// 	 	 	Z Dist 	DCA(number of sigmas)
CsAverPattern	CUTS	1600	10
CsKalmanFitting Specials	0	// 0 - Off, 1 -  mu' takes precedence
CsAverPattern	Refit		1	// 0 - Off, 1 - ON
CsAverPattern	Retrack		3	// 0 - Off, !0 = cut on event time
CsKalmanFitting	Hist		1	// 0 - Off, 1 - ON
CsKalmanFitting	Print [ 0 ]	0	// Kalman filter info.
CsKalmanFitting	Print [ 1 ]	0	// Global fit info.
CsKalmanFitting	Print [ 2 ]	0	// IKF info.
CsKalmanFitting	RefPlane	1350	// in mm
CsKalmanFitting	RefBeam 	-8000	// in mm
CsKalmanFitting	RefMargin	500	// in mm
// Chi2 increment cuts:		DirChi2Trk	InvChi2Trk
CsKalmanFitting	CUTS		7		17
// Vertex Chi2 cut:
CsKalmanFitting	Chi2VCut	15
CsKalmanFitting BeamP0		190	// Beam reference momentum (for histo purposes)


//					   ***** RICH1
make rich1 reconstruction
include ${CORAL}/src/pkopt/rich1.2014.opt
// RICH upgrade: 12 APV 4 MAPMT
// (Note: not sure this is need for RD, as opposed to MC)
CsRICH1UpGrade	CONFIG
// Option for faster processing: but "the output buffer of rich won't have the angle that maximizes the likelihood anymore" (according to Giulia, cf. mail from 2009/06/29)
RICHONE 	DoThetaLikeMax	NO

//					   ***** CALORIMETERS
//make	calorimeters	reconstruction
//include ${CORAL}/src/pkopt/calorim_2012.opt


//					   ***** MU' PID +  FLUX SCALERS
include ${CORAL}/src/pkopt/trigger.2014.opt
#else
//					   ***** NO Beam RECONSTRUCTION

//					   ***** NO Vertex RECONSTRUCTION

//					   ***** NO RICH1

//					   ***** NO CALORIMETERS

//					   ***** NO MU' PID
#endif

//  	 	 	 	 	   ***** GEOMETRICAL ZONES
define	 zone	 -700	  3500	 before M1
define	 zone	 3500	 17000	 between M1 and M2
define	 zone	17000	 32500	 between M2 and Muon Wall
define	 zone	32500	 99999	 after Muon Wall

define	 zone	-8000	 -3500	 before the target
define	 zone	-3500	  -700	 zone of the vertex detector

//  	  	  	  	  	   ***** MAGNETIC FIELD MAPS
CsField SOL_field	 $COMPASS_FILES/maps/mag_fields/SOL/OD_dipole.fieldmap
CsField SM1m_field_measured $COMPASS_FILES/maps/mag_fields/SM1m/SM1M.map.172.data
CsField SM2_field	 $COMPASS_FILES/maps/mag_fields/SM2/FSM.map.5000.data

//					   ***** MAGNETs INFO
#if C2OPT_COMPUTERFARM <= 4
CsMagInfo	MySQLDB
#else
CsMagInfo	File	~ybedfer/public/magInfo.2012.txt	// ...if MySQL not avail
#endif
CsMagInfo	SM2	1



//					   ***** LOGGER
error logger log level		error	// debugging, verbose, info, anomaly, warning, error, fatal
error logger store level	none
error logger verbosity		normal	// low, normal, high. 

//	  	  	  	  	   ***** TraFDic (+DET OPTIONS+...)

include $CORAL/src/pkopt/trafdic.2015.opt

//		***** OVERWRITE WHAT'S SET IN "../pkopt/trafdic.2015.opt" *****

// ***** TEMPORARY SETTING FOR NEW DICO STYLE
// (The inclusion of DC05 necessitates a modification of the structure of the
// dico files (those of 2015, and only those). The modified files are stored in
// an alternative directory (viz.: "dico2"). This, temporarily, until most of
// the users can reasonably be assumed to have updated their private copy of
// coral. Then "dico2" is going to be renamed "dico", w/ a symlink to the
// directory name "dico2", so as to ensure backward compatibility)
TraF	Dicofit		$COMPASS_FILES/geometry/2015/dico2/

// ***** MUON BLOCK (to be enabled for muon data)
// I) BMS RECONSTRUCTION
//make beam reconstruction
//include ${CORAL}/src/pkopt/beam_2014.opt
// II) Beam ID and momentum spread
//TraF	iCut	[30]	0	// Beam particle: 0=muon
//TraF	dCut	[ 5]	1.80e-4 // Beam cop spread (c/GeV). Assigned to d(1/p). Here corresponds to the 6.5GeV/sqrt(12) of the 190 Gev beam.

// ***** SPECIAL SETTINGS FOR DY 
// - Enlarge road width in zone 0x2
//  There are 3 alternatives: LIP (bold), Saclay (modest), Marcia (intermediate)
//  (Note that anyway dPRpar[11] is redefined infra)
//TraF	dPRpar	[ 10 - 13 ]	2.5  0.500 	0.110	0.050	// Saclay
//TraF	dPRpar	[ 10 - 13 ]	10.  0.500      0.100   0.060	// LIP
TraF	dPRpar	[ 10 - 13 ]	 3.  0.500 	0.100	0.060	//from Marcia
// - Special ForeTracking to MA/HG (in "ForeTrack2MA"):
//   - Candidate (for extension into MAs) tracks not required to be reco'd in (
//    or beyond) PS01s, even if they fall into their acceptance (and other
//    relaxed cceptance cuts).
//TraF	ReMode	[49]	1
//   - If in addition ReMode[49]&0x2: looser reco criteria + enlarged road widths
//TraF	ReMode	[49]	3
// - Options for Vertex Detector: as of 2015/06, let's disable the special reco
TraF	ReMode	[49]	3	// !&=0x4 disables BackTracking to VD
TraF	DY_VD_time	5.0
TraF	DY_VD_Chi2	0.3
// Loosen hits search in MAs (see "TEv::ForeTrack2MA")
TraF	ReMode	[50]	1

// ***** BEAM&VERTEX BLOCK
// (It's a set of options meant to improved the reconstruction of the
// beamTelescope and of the primary vertex.)
// * DETECTOR OPTIONS
// - Time gates on HMs: temporarily opened, to make around possible bad timing,
//  as reported by Catarina in 2015/09.
//HM04X1_d	HitTime [0-1] -10 10
//HM05X1_d	HitTime [0-1] -10 10
//HM04Y1_d	HitTime [0-1] -10 10
//HM05Y1_d	HitTime [0-1] -10 10
//HM04X1_u	HitTime [0-1] -10 10
//HM05X1_u	HitTime [0-1] -10 10
//HM04Y1_u	HitTime [0-1] -10 10
//HM05Y1_u	HitTime [0-1] -10 10
// - Time gates on HG: temporarily relaxed, to explore the impact.
//HG01Y1__   HitTime [0-1] -16.  16.
//HG02Y1__   HitTime [0-1] -16.  16.
//HG02Y2__   HitTime [0-1] -16.  16.
// - Time gates on the FIs of the BeamT: opened (temporarily, very widely) to
//  allow for the reconstruction of accidentally coincident beams, in particular
//  decay muons, which can later be associated to their continuation in the
//  spectrometer and hence removed from BpV in a more documented way.
//FI01X1__ HitTime [0-1]  -10.0  10.0	// SciFI: Time Cut (in ns)
//FI01Y1__ HitTime [0-1]  -10.0  10.0
//FI15X1__ HitTime [0-1]  -10.0  10.0
//FI15Y1__ HitTime [0-1]  -10.0  10.0
//FI15U1__ HitTime [0-1]  -10.0  10.0
//FI03X1__ HitTime [0-1]  -10.0  10.0
//FI03Y1__ HitTime [0-1]  -10.0  10.0
//FI03U1__ HitTime [0-1]  -10.0  10.0
// - Scifi clusterisation
//FI	ClusterDTime	6	// Max. Time Diff in cluster in sigma (assuming sigma = 0.6s)
// * TRACKING in BEAMT
// - Pattern Recognition requirements:
//   - dPRpar:  angular range (now taking tilt due to chicane into account)
//    narrowed (and so are the extra parameters enlarging the road width in
//    later iterations).
//   - iPRpar: more demanding min. #hits.
//TraF	dPRpar [ 40 - 43 ]	3.   0.004	0.005	0.005
//TraF	iPRpar [ 40 - 47 ]	2   3   7	2  2  6  	2 2     // 1 FI(XY) + 2FIs(XYU)
// - Weight of timing
//TraF	iCut [37] 2		// >0: Scifi specifics in PR of BeamT (except enhancement), >1: Allow for one out of time scifi hit.	(TAlgo::FindSpace)
//TraF	iCut [38] 5		// Weight (in %) of timing in CleanTrackList for BeamT (TEv::CleanTrackList)
//TraF	dCut	[63]	4.	// Coeff. of timing in Space PR Q function (applies only to beamT)
// * VERTEXING
//CsKalmanFitting	BpVDY		1	// Special BestPrimaryVertex for DY	
//CsKalmanFitting	BpVTimeW	1	// Weight of timing in BestPrimaryVertex for DY
//CsAverPattern	TimePrimCut	6 	// Time diff between track and beam, in sigmas
// Let's not refit/retrack at the pattern recognition stage, where no good
// determination of the BestprimaryVertex is available. Let's postpone the
// decision to retrack to the fitting stage.
//CsAverPattern	Refit		0	// 0 - Off, 1 - ON
//CsAverPattern	Retrack		0	// 0 - Off, !0 = cut on event time
//CsKalmanFitting	Retrack		2	// 0 - Off, !0 = cut on event time

// ***** SPECIAL DY2014, temporarily all ported to DY2015
// - TENTATIVE SETTINGS TO MAKE FOR LOW REDUNDANCY in...
//  ...BEAM TELESCOPE
TraF	iPRpar [ 10 - 17 ]      3   2   9	3  2  7		2 2	// between SM1 and SM2
// - LARGER BEAM SPREAD (due to open collimators?)
TraF	dCut	[ 5]	1.0e-4	// Beam cop spread (c/GeV). If >0: assigned to d(1/p). Here corresponds to ~4 GeV spread on p.	(TEv::TracksFit2)

// ***** SPECIAL TEMPORARY SETTING for MPs
TraF	dCut	[98]	.1	// Enlarge search road along v (MPs not yet v-aligned)
TraF	ReMode	[48]	2	// Strict angle cut in proj. search (MPs still noisy or not yet time-cut)

// ***** SPECIAL DY2015 EARLIER DATA TAKING
// - If DC5 not yet calibrated
//DC05Y1__ IS_MWPC
//DC05Y2__ IS_MWPC
//DC05U1__ IS_MWPC
//DC05U2__ IS_MWPC
//DC05V1__ IS_MWPC
//DC05V2__ IS_MWPC
// - In earliest data: ST3 in bad shape, while DC5 not working =>
//   - ZONE 0x2: loose requirements
//TraF    iPRpar [ 10 - 17 ]      3   2   7	3  2  7		2 2	// between SM1 and SM2
// - So loose requirements means susceptibility to ghost combinations.
//  => In order to limit the appearance of those, let's limit the number of
//    combinations considered by tightening the max. angle
//TraF	dPRpar	[ 11 ]	0.100
//   - The following, meant to save CPU, turned out to have a beneficial effect
//    on tracking efficiency. 
//TraF	ReMode	[47]	3	// 0x1: Recycle FullKF fit, 0x2: No systematic refit after cleaning.


// Redefine the TBNAMES of DETECTORS to be EXCLUDED FROM PATTERN RECOGNITION:
// - Some detectors may have to be left turned off in any case, cf. explanations in "../pkopt/trafdic.2015.opt".
// - But one may want, or may have, to turn off some more, on a period by period basis.
// - Also, in coral jobs dedicated to detectors studies, one usually turns off the particular detector under exam. Note that then your redefinition will overwrite all "DetNameOff" lists previously entered: take care to include in it all the detectors that have to be turned off in any case.
//TraF	DetNameOff	DC05X	// X planes not working for all of 2015

//		==> ==> IMPLICATIONS:
// DETECTORS being EXCLUDED may imply updating some PATTERN RECOGNITION options:
//TraF	ReMode	[17]	0	// Disable GEM amplitude correlations, if any GEM is among the turned-off detectors.

// Activate Dead Zones of DCs and GEMs: to process magnets-ON alignment data
//TraF	DZisActive	DC	GM	// Activate DZ of DCs and GEMs

// Graphics and printing options
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
// calo  EC01P1__  1111.500 - 45./2 - 10. = 1111.500 - 32.5 = 1079.0
// calo  EC02P1__  3325.200 - 45./2 - 10. = 3325.200 - 32.5 = 3292.7
// calo  HC01P1__  1267.500 -100./2 - 10. = 1267.500 - 60.0 = 1207.5
// calo  HC02P1__  3576.000 -112./2 - 10. = 3576.000 - 66.0 = 3510.0
// (NOTA BENE: The positions of E/HCAL1 may in fact differ from what is quoted supra and corresponds to the 2006/7 setup.)
TraF	SmoothPos [1]	1079.0	// Smooth track at ECAL1
TraF	SmoothPos [2]	1207.5	// Smooth track at HCAL1
TraF	SmoothPos [3]	3292.7	// Smooth track at ECAL2
TraF	SmoothPos [4]	3510.0	// Smooth track at HCAL2
#endif


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
TraF	DetNameOff	SI03	H	MB
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
TraF	DetNameOff	M	D	ST	P	H	GM01	GM02	GM03	GM04	GM05	GM06
#  elif C2OPT_FLAG == 4
// ***** II) Tracks reconstructed in both spectro and scifi/Si telescope,...
TraF	ReMode	[26]	1	// != 0 - Bridging over target
TraF	Hist	[ 1]	17	// Monitor, &0x10:Residuals
TraF	Hist	[16]	16	// Residuals: |Groups of det's  histo'ed
TraF	Hist	[17]	0	// Residuals: |Groups of inactive det's histo'ed
TraF	Hist	[18]	23	// Residuals: &Groups of tracks histo'ed
TraF	DetNameOff	M	D	ST	P	H	GM01	GM02	GM03	GM04	GM05	GM06	SI03
TraF	Det2Go2Fit	FI02	SI02X	SI02Y	FI03	FI04	FI05	FI06	FI07	FI08	GM
#  endif
#endif

//  ===== TEMPLATE FOR MASS PRODUCTION (uncomment, and modify if needed) =====

//					***** PLACEHOLDERS for I/O
// (Note: One has to uncomment next 7 entries.)
//events to read 200000
//Data		file	TMPRAWFILE
//mDST		file	TMPMDSTFILE
//histograms	home	TMPHISTFILE
//CsEvent         WriteBack	testevtdump.raw	Randoms		// All 0x1ffefc00 triggers
//CsKalmanFitting WriteBack	testevtdump.raw	Physics		// D*, excl. phi, K0, J/psi
//TMPRICHFILE

//					***** mDST OPTIONS
// (Note: One has to systematically uncomment next 2 entries. For special MegaDST production, one has to uncomment also the 3rd one.)
//mDST	select	0	// 0: no selection
//mDST	selectTrigger	1152	// Random trigger: unconditional output (Note: the option is not useful when "mDST select 0", i.e. all-out unconditional output.)

//mDST	hits	ALL		// SPECIAL OPTION! for the production of MegaDST


//					***** CDB OPTIONS
// (Note: One has to freeze the CDB info accessible to the production, so that the latter be independent of CDB later developments, and hence stable. A reasonable choice for the "entrytime" is the date of the start of production. One may still want to be more restrictive on all of the CDB. Or on the part of it concerning a particular detector or set of detectors: in that case, use the option "<generic name> CDBentrytime" (cf. examples in the body of this options file).)
//CDB	entrytime	2007-12-01-00:00:00	// TO BE MODIFIED!


//					***** DETECTOR MONITORING
// (Note: Have to uncomment next 4 entries.)
//Monitoring	Directory	ON	// Paolo's data monitor graphs.
//Monitoring	Residuals	ON	// Albert's monitoring histos.
//Monitoring	Efficiency	ON
//CsBuildPart	Hist		5	// Calorimeter


//					***** SETUP FILES
// (Note: One has to freeze the setup info accessible to the production, so that the latter be independent of the modifications later brought to the setup files, and hence stable. This is achieved by copying the setup files to sub-directories of the production directory. E.g. infra are listed what need be done for the test production long2011test9.) 
//decoding	map	/afs/cern.ch/compass/scratch/d17/objsrvvy/generalprod/testcoral/long2011test9/maps
//detector	table	/afs/cern.ch/compass/scratch/d17/objsrvvy/generalprod/testcoral/long2011test9/alignement/
//TraF		Dicofit	/afs/cern.ch/compass/scratch/d17/objsrvvy/generalprod/testcoral/long2011test9/dico/
//CsROOTGeometry file	/afs/cern.ch/compass/scratch/d17/objsrvvy/generalprod/testcoral/long2011test9/ROOTGeometry/det_2010.C

//					***** DETECTORS EXCLUDED FROM TRACKING
// (Note: One has to uncomment the following entry depending upon the data taking period one is about to produce. It allows to exclude some detectors from the track reconstruction, had they been failing during a, large, part (then production gets more stable) or all (tracking performs better if it is told not to expect hits from a particular detector) of the period. (Caveat: one has to supply the full list of excluded detectors each time one re-specifies this "TraF DetNameOff" option.))

//					***** UPDATE PATTERN RECOGNITION options
// (Note: Excluding detectors may imply updating some pattern recognition options. In the present case, there's nothing we can do to make for the loss of DC04Y.)
// ===========================================================================
//end
