C2OPT_HEADER
#if 0

// $Id: trafdic.2004.opt.c,v 1.27 2011/02/15 23:59:07 ybedfer Exp $

// coral options c file for TraFDic on 2004 data
//
// NOTA.BENE.: c file is to be processed with:
//
//        "gcc -E -P -C -D<defines>"
//
// before being supplied to coral. 

// From "./trafdic.2004.opt,v 1.3"
// Modified to be made a c file, able to produce opt files for various
// COMPUTER SYSTEMs: (1) lxplus, (2) gridKa, (3) Lyon
// $Log: trafdic.2004.opt.c,v $
// Revision 1.27  2011/02/15 23:59:07  ybedfer
//  - Add re-tracking in "CsAverPattern", in line w/ what's been done for
//   later data taking years.
//  - Cosmetics...
//
// Revision 1.26  2010/12/31 14:37:18  suhl
// * remove obsolete option 'mDST digits', replace it with 'mDST Digits' and 'mDST DAQdigits' in the template options file
//
// Revision 1.25  2010/11/04 23:55:41  ybedfer
//  - Bug fix: "BeamRecons BackTrack" entry supplied. The default for 2004 is
//   selected, for want of re-calibration (as was done for 2003).
//  - Several settings aligned w/ what is done in other years:
//    - Drift "make always two clusters".
//    - "DW CDBentrytime".
//    - Allow beams w/o BMS in the vertexing.
//    - "CsAverPattern AcceptTill": 90 -> 105.
//    - "CsAverPattern TimePrimCut 10".
//    - STRAW settings spacers: NO -> YES. (Caveat: not tested.)
//
// Revision 1.24  2010/05/10 10:51:35  ybedfer
//  Add COMPUTERFARM == 4 for the only system, w/ server pccodb00.
//
// Revision 1.23  2009/11/05 02:47:46  ybedfer
//  - Less verbose RICH.
//  - Smooth points 10 cm upstream front surface of a calorimeter.
//
// Revision 1.20  2007/05/02 01:01:37  ybedfer
//  - A unique entry for all "STRAW settings".
//
// Revision 1.19  2007/03/13 08:53:59  ybedfer
//  - "STRAW" options:
//    - Enable entry for "signal_propagation". This is now required by the
//     newly committed CsStraw software.
//    - X-ray correction is kept disabled: have to be commented back in if
//     required, i.e. if relevant detectors.dat file takes it into account.
//
// Revision 1.18  2007/03/13 07:46:05  ybedfer
//  By default, get input data from castor (instead of from oracle).
//
// Revision 1.17  2006/12/06 08:30:54  ybedfer
//  "rich1.2003.opt" -> "rich1.2004.opt".
//
// Revision 1.16  2006/11/30 01:39:22  ybedfer
//  "rich1.2003.opt" -> "rich1.2004.opt"
//
// Revision 1.15  2006/08/23 00:05:53  ybedfer
//  Do not book mDSTs ifndef "C2OPT_BARE_FILE".
//
// Revision 1.14  2006/08/22 10:41:35  ybedfer
//  - Material map: 2003.01 -> 2004.00.
//  - "C2OPT_BARE_FILE" macro.
//
// Revision 1.13  2006/08/05 03:47:13  ybedfer
//  Use the 2004 version of the beam options file.
//
// Revision 1.12  2006/08/01 04:41:55  ybedfer
//  - Add an entry for enabling X-ray correction (commented out by default).
//  - Add a cpp conditonal ("C2OPT_VOID_INPUT") for creating a version of the
//   options file w/o input data specification.
//  - Remove smooth point @ ECAL2.
//
// Revision 1.11  2006/03/14 17:29:28  ybedfer
//  Added smoothing points @ the abscissa of the calorimeters.
//
// Revision 1.10  2005/10/27 02:03:55  ybedfer
//  Add post-vertexing tracks refit ("CsAverPattern Refit").
//
// Revision 1.9  2005/07/11 12:59:14  ybedfer
//  "cchpsscompass" tape server @ Lyon.
//
// Revision 1.8  2005/05/22 23:54:12  ybedfer
//  Option to update solenoid map from DB as a function of time. Useful for
// processing target rotation data. Disabled (commented out) by default.
//
// Revision 1.7  2005/01/19 14:10:34  ybedfer
//  Calib files moved to ``specialplace'' @ Lyon.
//
// Revision 1.6  2004/09/21 23:56:22  ybedfer
//  - No more CDB entry time limit for DWs and STs.
//  - Specification of vertexing ref. planes, overwriting
//   "../pkopt/vertex.2002.opt".
//
// Revision 1.5  2004/08/31 10:17:41  ybedfer
//  New mySQL server.
//
// Revision 1.4  2004/06/17 23:19:08  ybedfer
//  - Undo previous change: do not output all tracks' helices by default.
//  - Correct wrong comment concerning mDST select(ion).
//
// Revision 1.3  2004/06/17 11:58:56  ybedfer
//  mDST select: have "3" as the default, i.e. all booked helices go to output.
//
// Revision 1.2  2004/06/17 10:44:46  ybedfer
//  PHAST selection: new syntax, valid for version >= 7.
//
// Revision 1.1  2004/06/08 02:48:17  ybedfer
//  Template options file.
//

#endif
#if !defined C2OPT_COMPUTERFARM
#  define C2OPT_COMPUTERFARM 1
#endif
#if   C2OPT_COMPUTERFARM == 1
// TraFDic on 2004 data @ CERN.
#elif C2OPT_COMPUTERFARM == 2
// TraFDic on 2004 data @ gridKa.
#elif C2OPT_COMPUTERFARM == 3
// TraFDic on 2004 data @ Lyon.
#elif C2OPT_COMPUTERFARM == 4
// TraFDic on 2004 data @ Compass ONLINE.
#else
// TraFDic on 2004 data.
#endif

#ifdef C2OPT_BARE_FILE
//  ``Bare'' options file to be included in "$CORAL/src/alignment/traf.*.opt"
// files for processing 2004 FIELD ON data for alignment purposes:
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

//  The file is to be edited for:
//        Enabling mass-production specific templates and options, cf. the
//       block entitled "TEMPLATE FOR MASS PRODUCTION" in fine.
// Or otherwise for:
//    I)  Specifying input data one wants to process instead of the default
//   II)  Setting accordingly:
//        - Target magnetic field, cf. the (inappropriately named)
//         "CsField SOL_field" option infra.
//        - In case of transversity, full path to detector table and Dico file.
//  III)  Redefining the # of events to process, cf. "events to read" infra.
//   IV)  Adding whatever options are needed to enable a given output:
//      - mDST (enabled by default, but active only if executable is
//       PHAST-compatible, cf. "./README").
//        Default mDST selection is:
//        - Event with at least one vertex or one particle.
//        - Only first and last helix of every track is stored.
//       Cf. $PHAST/coral/README for how to change this behaviour.
//      - Write-back of raw data.
//    V)  Selection criteria: Trigger selection, BoS rejection, etc...
//   VI)  CDB entrytime.
//  VII)  Enabling X-ray correction: cf. "STRAW" entry infra. It's disabled by
//      default. Should only be enabled if alignment is also X-ray corrected.
//       => Check "detector table".
// VIII)  Enabling VarPitch correction: cf. "VarPitch" entry infra.
//   IX)  Enabling the updating the solenoid field, if processing a target
//      rotation run.


// Tip: Set your "emacs" in "c++-mode" to view this ".opt" file. And have this
//     done automatically by setting the appropriate hook in your "~/.emacs.el"
//     init file, cf. "~ybedfer/.emacs.el".
#endif

#ifndef C2OPT_BARE_FILE
//                                         ***** PHAST OUTPUT
// - Active only if used on a PHAST-compatible executable, cf. "./README"
// - Cf. "$PHAST/coral/README" for the details of the options
mDST file ${HOME}/phast.root
// - Information reduction bit pattern: "mDST select".
mDST select		9	// At least one vertex or particle | First and last helices only.
// - One may want to output more info, as exemplified infra:
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
//					   ***** INPUT DATA...
#  if   C2OPT_COMPUTERFARM == 1

//						...either raw file
Data file 	/castor/cern.ch/compass/data/2004/raw/W39/cdr16001-41063.raw

//						...or retrieved from oracle DB
//Database select 	oracle
//Data type		raw

//Data year 		2004
//Data run select 	41063
//Data container  	cdr16001-41063
#  elif C2OPT_COMPUTERFARM == 2
Data file 	/grid/fzk.de/compass/compass2/data/raw/2004/xxx.raw
#  elif C2OPT_COMPUTERFARM == 3
Data file 	cchpsscompass:/hpss/in2p3.fr/group/compass/data/2004/W28/raw/cdr16001-37275.raw
#  elif C2OPT_COMPUTERFARM >= 4
Data file 	/castor/cern.ch/compass/data/2006/raw/W40/cdr22001-52613.raw
#  endif
#endif
// 					   ***** DETECTOR TABLE...
//  A full directory (path name terminated by a slash "/") or the full path to
// a particular file can be specified.
detector table   	$COMPASS_FILES/geometry/2004/
//detector table $COMPASS_FILES/geometry/2004/detectors.39486.transv.dat


// 					   ***** MATERIAL MAPS
CsMaterialMap Zone_1		$COMPASS_FILES/maps/material/2004.00/material.target.m2004.00.map
CsMaterialMap Zone_2		$COMPASS_FILES/maps/material/2004.00/material.saclay.m2004.00.map
CsMaterialMap Zone_3		$COMPASS_FILES/maps/material/2004.00/material.straw.m2004.00.map
CsMaterialMap Zone_4		$COMPASS_FILES/maps/material/2004.00/material.muf2.m2004.00.map
CsMaterialMap Zone_5		$COMPASS_FILES/maps/material/2004.00/material.muf3.m2004.00.map
CsMaterialMap Zone_6		$COMPASS_FILES/maps/material/2004.00/material.hcal1.m2004.00.map
CsMaterialMap Zone_7		$COMPASS_FILES/maps/material/2004.00/material.muf1.m2004.00.map
CsMaterialMap Zone_8		$COMPASS_FILES/maps/material/2004.00/material.rich.m2004.00.map
CsMaterialMap ELossTarget	$COMPASS_FILES/maps/material/2004.00/material.target_de.m2004.00.map

//					   ***** CALIBRATION
#if C2OPT_COMPUTERFARM > 4
use calibration
CDB use  	FileDB
CDB location  	$COMPASS_FILES/calibrations/2004
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
//CDB	entrytime	2007-12-01-00:00:00
// "entrytime" entry can be made to affect only a (set of) specific TBName(s):
// "TBname(or TB) CDBentrytime(sic!)"
DW CDBentrytime 2007-12-01-00:00:00 	// Disregarding the calibs entered by error on Dec. 14th, 2007
// CDB can be replaced by a File DB, for any given (set of) TBname(s):
#  if   C2OPT_COMPUTERFARM == 1
//DW FileDB	/afs/cern.ch/user/l/leberig/public/marcin/
#  else
//DW FileDB	$COMPASS_FILES/fileDBs/DW/marcin
#  endif
#endif

//					   ***** MAPPING
decoding map	$COMPASS_FILES/maps/2004.xml


//					   ***** EVENT SELECTION
events to read 30000 	// # of events to read
//events to skip 843 //139 //705 //229 //843
//events BOS_skip 1.22			// Skip first ~400(? cf. J�rgen) ms in spill (disabled in 2004)
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
make beam reconstruction
include ${CORAL}/src/pkopt/beam_2004.opt
//                   Overwrite what's set in "beam.opt"
BeamRecons	useTRAFFIC	1	// >0: Traffic is used, =0: standalone beam package is used
BeamRecons	AllowMultCand	2
PID_doBeamID	selAllBeams	1	// Allow beams w/o BMS in the vertexing
//			Bactracking BMS <- beamTelescope
// Specify the parameters w/ help of "BackTrack ResBMS/Corr" infra. They tell the resolution of the backtracking at the various BMS planes (only 4 of the planes, used so far, as of 2010/04) and the offset corrections to be applied.
BeamRecons BackTrack Default		// Defaults are available for 2004.

//					   ***** Vertex RECONSTRUCTION
make vertex reconstruction
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
CsAverPattern	CUTS	1600	10
CsAverPattern	Refit	1	// 0 - Off, 1 - ON
CsAverPattern	Retrack	3	// 0 - Off, !0 = cut on event time
CsKalmanFitting  Hist	1	// 0 - Off, 1 - ON
CsKalmanFitting  Print [ 0 ]	0	// Kalman filter info.
CsKalmanFitting  Print [ 1 ]	0	// Global fit info.
CsKalmanFitting  Print [ 2 ]	0	// IKF info.
CsKalmanFitting  RefPlane	500	// in mm
CsKalmanFitting  RefBeam	-1500	// in mm
CsKalmanFitting  RefMargin	500	// in mm
// Chi2 increment cuts:		DirChi2Trk	InvChi2Trk
CsKalmanFitting  CUTS		7		10.75
// Vertex Chi2 cut:
CsKalmanFitting  Chi2VCut	10.5
CsKalmanFitting BeamP0 160 	// Beam reference momentum (for histo purposes)

//					   ***** RICH1
make rich1 reconstruction
include ${CORAL}/src/pkopt/rich1.2004.opt
// Have RICH somewhat less vebose than what's set in "rich1.2004.opt"
RICHONE PrintDetMirr 0
RICHONE AcknoMethods 0 	// prints methods in use
RICHONE PrintKeys 0 	// 1 = print keys read,  0 = No print
RICHONE PrintConsts 0 	// 1 = print rec. constants read,  0 = No print

//					   ***** CALORIMETERS
make calorimeters reconstruction
include ${CORAL}/src/pkopt/calorim.opt

//					   ***** MU' PID
include ${CORAL}/src/pkopt/trigger.2003.opt
#else
//					   ***** NO Beam RECONSTRUCTION 

//					   ***** NO Vertex RECONSTRUCTION

//					   ***** NO RICH1

//					   ***** NO CALORIMETERS

//					   ***** NO MU' PID
#endif

//  	 	 	 	 	   ***** GEOMETRICAL ZONES
define zone 0 3500 before M1
define zone 3500 17000 between M1 and M2
define zone 17000 32500 between M2 and Muon Wall
define zone 32500 99999 after Muon Wall

define zone -8000 0 before the target


//  	  	  	  	  	   ***** MAGNETIC FIELD MAPS
// Select entry relevant for longitudinal or transverse, among following 2:
CsField SOL_field	 $COMPASS_FILES/maps/mag_fields/SOL/smctgt_solenoide.fieldmap
//CsField SOL_field	 $COMPASS_FILES/maps/mag_fields/SOL/smctgt_dipole.fieldmap

CsField SM1m_field_measured $COMPASS_FILES/maps/mag_fields/SM1m/SM1M.map.172.data
CsField SM2_field	 $COMPASS_FILES/maps/mag_fields/SM2/FSM.map.4000.data

//					   ***** MAGNETs INFO
#if C2OPT_COMPUTERFARM <= 3
CsMagInfo	MySQLDB
#else
CsMagInfo	File	~ybedfer/public/magInfo.2004.txt	// ...if MySQL not avail
#endif
CsMagInfo	SM2	1		// ==0: Do NOT use NMR, !=0: Rescaling factor (applied to NMR corrected map)
//CsMagInfo Update solenoid		// To be enabled during target rotation



//					   ***** LOGGER
error logger log level		error	// debugging, verbose, info, anomaly, warning, error, fatal
error logger store level	none
error logger verbosity		normal	// low, normal, high. 

//	  	  	  	  	   ***** TraFDic (+DET OPTIONS+...)

include ${CORAL}/src/pkopt/trafdic.2004.opt

//		***** OVERWRITE WHAT'S SET IN "../pkopt/trafdic.2004.opt" *****

//  Specify full path to some alternative dico file
//TraF	Dicofit		$COMPASS_FILES/geometry/2004/dico/dico.39486.transv

//TraF	Graph	[0]	0	// coral's graphics switches off automatically under BATCH @ CERN, gridKa and Lyon. This option forces it to do so in any case. 
//TraF	iCut	[0]	7	// Trigger selection in TraFDic
//TraF	Print	[0]	3	// Event-by-event numbering
//TraF	Print	[8]	3	// EoE
#ifndef C2OPT_BARE_FILE
//		   Smoothing points
// The points @ RICH (resp. Calos) are to used by RICH (resp. Calo) software.
TraF	SmoothPos [0]	745.		// RICH
// From: Vladimir Kolosov
// Sent: Thu 2/28/2008 9:14 AM
// Let's decide to calculate track position 10 cm upstream front surface of a calorimeter.
// calo  HC01P1__  1267.500 -100./2 - 10. = 1267.500 - 60.0 = 1207.5
// calo  HC02P1__  3576.000 -112./2 - 10. = 3576.000 - 66.0 = 3510.0
TraF	SmoothPos [1]	1207.5		// Smooth track at HCAL1
TraF	SmoothPos [2]	3510.0		// Smooth track at HCAL2
#endif

//                 STRAW options: X-ray and signal propagation correction
//  X-ray can only be enabled if alignment (cf. "detector table" entry supra)
// did also take X-ray correction into account.
//  And in any case, is incompatible w/ "VarPitch", cf. infra.
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
// "RDMon_DIGIT_GM" defined)
TraF	Hist	[ 1]	17	// MCMonitor (&2:D0Kpi,&4:Scat'd mu,&8:Lppi)
TraF	Hist	[18]	2	// Residuals: &Groups of tracks histo'ed
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
// (Note: One has to systematically uncomment next 2 entries. For special MegaDST production, one has to uncomment also the 3rd one.)
//mDST	select	0	// 0: no selection
//mDST	selectTrigger	3072	// Random trigger: unconditional output (Note: the option is not useful when "mDST select 0", i.e. all-out unconditional output.)

//mDST	hits	ALL		// SPECIAL OPTION! for the production of MegaDST


//					***** CDB OPTIONS
// (Note: One has to freeze the CDB info accessible to the production, so that the latter be independent of CDB later developments, and hence stable. A reasonable choice for the "entrytime" is the date of the start of production. One may still want to be more restrictive on all of the CDB. Or on the part of it concerning a particular detector or set of detectors: in that case, use the option "<generic name> CDBentrytime" (cf. examples in the body of this options file).)
//CDB	entrytime	2009-2-17-00:00:00	// TO BE MODIFIED!


//					***** DETECTOR MONITORING
// (Note: Have to uncomment next 3 entries. 4th one is optional.)
//Monitoring	Directory	ON	// Paolo's data monitor graphs.
//Monitoring	Residuals	ON	// Albert's monitoring histos.
//Monitoring	Efficiency	ON
//CsBuildPart	Hist		5	// Calorimeter


//					***** SETUP FILES
// (Note: One has to freeze the setup info accessible to the production, so that the latter be independent of the modifications later brought to the setup files, and hence stable. This is achieved by copying the setup files to sub-directories of the production directory. E.g. infra are listed what need be done for the production of 2003P1D/slot3.) 
//decoding	map	/afs/cern.ch/compass/scratch/d17/objsrvvy/2003/P1D/slot3/maps
//detector	table	/afs/cern.ch/compass/scratch/d17/objsrvvy/2003/P1D/slot3/alignement/
//TraF		Dicofit	/afs/cern.ch/compass/scratch/d17/objsrvvy/2003/P1D/slot3/dico/

//					***** BEAM BACKTRACKING PARAMETERS
// (Note: One has to uncomment the following two entries depending upon the data taking period one is about to produce. Their purpose is to specify the parameters for the backtracking BMS <- beam Telescope. Defaults have been tuned on P1E#29975. It was checked that they are approximately valid for P1D#28893 => therefore no large drifting observed, although the impact of the small yet finite drifting on the backtracking chi2 has not been quantified.)
//BeamRecons BackTrack ResBMS 4.66 3.98  2.55  2.20	// Tuned on 03P1E#29975
//BeamRecons BackTrack Corr   6.53 3.52 -1.94 -5.66	// Tuned on 03P1E#29975

// ===========================================================================
//end
