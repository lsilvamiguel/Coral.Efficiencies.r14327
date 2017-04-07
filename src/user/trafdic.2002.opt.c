C2OPT_HEADER
#if 0

// $Id: trafdic.2002.opt.c,v 1.13 2011/02/15 22:53:42 ybedfer Exp $

// coral options c file for TraFDic on 2002 data
//
// NOTA.BENE.:c file is to be processed with:
//
//        "gcc -E -P -C -D<defines>"
//
// before being supplied to coral. 

// Able to produce opt files for various
// COMPUTER SYSTEMs: (1) lxplus, (2) gridKa, (3) Lyon
// From "./trafdic.2003.opt,v 1.5" and "./trafdic.2002.opt,v1.9".
// $Log: trafdic.2002.opt.c,v $
// Revision 1.13  2011/02/15 22:53:42  ybedfer
//  Cosmetics...
//
// Revision 1.12  2010/12/31 14:37:18  suhl
// * remove obsolete option 'mDST digits', replace it with 'mDST Digits' and 'mDST DAQdigits' in the template options file
//
// Revision 1.11  2010/05/10 10:51:35  ybedfer
//  Add COMPUTERFARM == 4 for the online system, w/ server pccodb00.
//
// Revision 1.9  2007/05/02 01:01:37  ybedfer
//  - A unique entry for all "STRAW settings".
//  - New mySQL server.
//
// Revision 1.8  2007/03/13 08:50:29  ybedfer
//  "STRAW" options:
//   - Enable entry for "signal_propagation". This is now required by the
//    newly committed CsStraw software.
//   - X-ray correction is kept disabled: have to be commented back in if
//    required, i.e. if relevant detectors.dat file takes it into account.
//
// Revision 1.6  2006/08/22 10:43:17  ybedfer
//  - Material map: 2002.02 -> 2004.04.
//  - More smoothing points.
//  - "C2OPT_BARE_FILE" macro.
//
// Revision 1.3  2005/04/22 20:35:48  ybedfer
//  - Outer straws and MBs handled as drift detectors (all "IS_MWPC" entries
//   canceled).
//  - "ST CDBentrytime" set for straws.
//
// Revision 1.1  2004/09/21 23:46:38  ybedfer
//  Initial version.
//

#endif
#if !defined C2OPT_COMPUTERFARM
#  define C2OPT_COMPUTERFARM 1
#endif
#if   C2OPT_COMPUTERFARM == 1
// TraFDic on 2002 data @ CERN.
#elif C2OPT_COMPUTERFARM == 2
// TraFDic on 2002 data @ gridKa.
#elif C2OPT_COMPUTERFARM == 3
// TraFDic on 2002 data @ Lyon.
#elif C2OPT_COMPUTERFARM == 4
// TraFDic on 2002 @ Compass ONLINE.
#else
// TraFDic on 2002 data.
#endif

#ifdef C2OPT_BARE_FILE
//  ``Bare'' options file to be included in "$CORAL/src/alignment/traf.*.opt"
// files for processing 2002 FIELD ON data for alignment purposes:
//  - It lacks input data specification (this must be entered in the main
//   "traf.*.opt" file)...
//  - ...as well as mDST output, all high level (vertex and the like)
//   reconstruction options, etc..
//  - More generally, this file should not be edited: any modification to the
//   option entries it contains must instead be implemented in the main
//   "traf.*.opt" by overwriting.
#else
//  Restricted to magnets ON data.
//  Both longitudinal and transverse are catered for. Take care still to set
// the specification of the target magnetic field as need be: cf.
// "CsField SOL_field" infra.

// This options file is to be edited for:
//   I)  Specifying input data one wants to process instead of the default
//     (and setting target magnetic field accordingly).
//  II)  Redefining the # of events to process, cf. "events to read" infra.
// III)  Adding whatever options are needed to enable a given output, e.g.
//     mDST (default, but active only if executable is phast-compatible, cf.
//     "./README").
//       Default mDST "information reduction" bit pattern is:
//        - Event with at least one vertex or one particle.
//        - Only first and last helix of every track is stored.
//     Cf. $PHAST/coral/README for how to change this behaviour.
//  IV)  Selection criteria: Trigger selection, BoS rejection, etc...
//   V)  CDB entrytime.
//  VI)  Enabling X-ray correction: cf. "STRAW" entry infra. It's disabled by
//     default. Should only be enabled if alignment is also X-ray corrected.
//      => Check "detector table".


// Tip: Set your "emacs" in "c++-mode" to view this ".opt" file. And have this
//     done automatically by setting the appropriate hook in your "~/.emacs.el"
//     init file, cf. "~ybedfer/.emacs.el".
#endif

#ifndef C2OPT_BARE_FILE
//                                         ***** PHAST OUTPUT
// - Active only if used on a PHAST-compatible executable, cf. "./README"
// - Cf. "$PHAST/coral/README" for the details of the options
mDST file		${HOME}/phast.root
// - Information reduction bit pattern: "mDST select".
mDST select		9	// 0x9 = At least one vertex or particle | First and last helices only.
// - One may want to output more info, as exemplified infra:
//mDST hits		HL	// MegaDST...
//mDST Digits		EC	// GigaDST... (raw information after calibrations are applied, according to documentation only useful for non-tracking detectors, use mDST hits for tracking detectors)
//mDST DAQdigits	SC	// GigaDST... (raw information as read from raw data stream, no calibrations applied)
#endif

//					   ***** HISTOGRAMMING
histograms package ROOT

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
Data file 	/castor/cern.ch/compass/data/2002/raw_migration/P2D/10/cdr19006-22019.raw

//						...or retrieved from oracle DB
//Data type   raw
//Database select oracle

//Data year   2002
//Data run select 22019
//Data container cdr19006-22019
#  elif C2OPT_COMPUTERFARM == 2
Data file 	/grid/fzk.de/compass/compass2/data/raw/2002/P2D/cdr19006-22019.raw
#  elif C2OPT_COMPUTERFARM == 3
Data file	cchpsscompass:/hpss/in2p3.fr/group/compass/data/2002/P2D/raw/cdr19006-22019.raw
#  elif C2OPT_COMPUTERFARM >= 4
Data file 	/castor/cern.ch/compass/data/2002/raw_migration/P2D/10/cdr19006-22019.raw
#  endif
#endif
// 					   ***** DETECTOR TABLE...
detector table   	$COMPASS_FILES/geometry/2002/


// 					   ***** MATERIAL MAPS
CsMaterialMap Zone_1		$COMPASS_FILES/maps/material/2002.04/material.target.m2002.04.map
CsMaterialMap Zone_2		$COMPASS_FILES/maps/material/2002.04/material.saclay.m2002.04.map
CsMaterialMap Zone_3		$COMPASS_FILES/maps/material/2002.04/material.straw.m2002.04.map
CsMaterialMap Zone_4		$COMPASS_FILES/maps/material/2002.04/material.muf2.m2002.04.map
CsMaterialMap Zone_5		$COMPASS_FILES/maps/material/2002.04/material.muf3.m2002.04.map
CsMaterialMap Zone_6		$COMPASS_FILES/maps/material/2002.04/material.hcal1.m2002.04.map
CsMaterialMap Zone_7		$COMPASS_FILES/maps/material/2002.04/material.muf1.m2002.04.map
CsMaterialMap Zone_8		$COMPASS_FILES/maps/material/2002.04/material.rich.m2002.04.map
CsMaterialMap ELossTarget	$COMPASS_FILES/maps/material/2002.04/material.target_de.m2002.04.map


//					   ***** CALIBRATION
#if C2OPT_COMPUTERFARM > 4
use calibration
CDB use  	FileDB
CDB location  	$COMPASS_FILES/calibrations/2002
#else
use calibration
CDB use 	MySQLDB
#  if   C2OPT_COMPUTERFARM == 1
CDB server	wwwcompass	// This is, as of 07/05, pointing to "lxfsrb6103" located in computer center (reachable outside Cern)
#  elif C2OPT_COMPUTERFARM == 2
CDB server compass.gridka.de	// Located at gridKa (internal network)
CDB specialplace GRIDKA 	// Location of calibration files at GridKA
#  elif C2OPT_COMPUTERFARM == 3
CDB server	cccompassdb	// Located at Lyon   (internal network)
CDB portnumber  23306		// Port to be used in place of the std one
CDB specialplace LYON	 	// Location of calibration files at Lyon
#  elif C2OPT_COMPUTERFARM == 4
CDB server	pccodb00
CDB specialplace DAQ	 	// Location of calibration files in ONLINE
#  endif
// "CDB entrytime Y-M-D-H:M:S" can be entered to reject all calibrations issued
// after specified date. E.g., in order to reject all trigger matrix files:
// N.B.: CDB is constantly evolving => Leaving "CDB entrytime" disabled may
// give rise to unprecdictable behaviour.

// "entrytime" entry can be made to affect only a (set of) specific TBName(s):
// "TBname(or TB) CDBentrytime(sic!)"
// This is enabled here so that calibs from 2003/09 are selected, for they are,
// marginally, better, as measured by coral summary of perfs on run 22019.
ST CDBentrytime 2003-10-01-00:00:00	// Files registered before specif'd date

// CDB can be replaced by a File DB, for any given (set of) TBName(s):
#  if   C2OPT_COMPUTERFARM == 1
//DW FileDB	/afs/cern.ch/user/l/leberig/public/marcin/
#  else
//DW FileDB	$COMPASS_FILES/fileDBs/DW/marcin
#  endif
#endif

//					   ***** MAPPING
decoding map	$COMPASS_FILES/maps/2002.xml


//					   ***** EVENT SELECTION
events to read 26000			// # of evts to read
//events to skip 843 //139 //705 //229 //843
events BOS_skip 1.22			// skip first ~400(? cf. Jürgen) ms in spill
//selection trigger mask f		// Trigger selection (hexadecimal)
//selection trigger strict


// seed for random number generations:
random number engine JamesEngine
random number seed 19990102
//reset random seed every new event

//ST IS_MWPC
//MB IS_MWPC
//DW IS_MWPC

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
include ${CORAL}/src/pkopt/beam_2002.opt
//                   Overwrite what's set in "beam.opt"
BeamRecons    useTRAFFIC  1    // >0: Traffic is used, =0: standalone beampackage is used
BeamRecons AllowMultCand  2


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
// track is considered as accepted if its momentum is smaller
// than "AcceptTill" percents of beam track
CsAverPattern	AcceptTill	90	// in %
//                     Z Dist(from centre)  DCA(number of sigmas)
CsAverPattern	CUTS		1600			10
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


//					   ***** RICH1
make rich1 reconstruction
include ${CORAL}/src/pkopt/rich1.2002.opt

//					   ***** CALORIMETERS
make calorimeters reconstruction
include ${CORAL}/src/pkopt/calorim.opt

//					   ***** MU' PID
include ${CORAL}/src/pkopt/trigger.2002.opt
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
CsMagInfo	File	~ybedfer/public/magInfo.2002.txt	// ...if MySQL not avail
#endif
CsMagInfo	SM2	1		// ==0: Do NOT use NMR, !=0: Rescaling factor (applied to NMR corrected map)




//					   ***** LOGGER
error logger log level		error	// debugging, verbose, info, anomaly, warning, error, fatal
error logger store level	none
error logger verbosity		normal	// low, normal, high. 

//	  	  	  	  	   ***** TraFDic (+DET OPTIONS+...)

include ${CORAL}/src/pkopt/trafdic.2002.opt

//                   Overwrite what's set in "../pkopt/trafdic.2002.opt"
//TraF	Graph	[0]	0		// For INTERACTIVE coral w/o graphics
//TraF	iCut	[0]	7		// Trigger selection in TraFDic
//TraF	Print	[0]	3		// Event-by-event numbering
//TraF	Print	[8]	3		// EoE
#ifndef C2OPT_BARE_FILE
//		   Smoothing points
// The points @ RICH (resp. Calos) are to used by RICH (resp. Calo) software.
TraF	SmoothPos [0]	745.		// RICH
TraF	SmoothPos [1]	1266.5		// HCAL1
TraF	SmoothPos [2]	3575.		// HCAL2
#endif

//                 STRAW options: X-ray and signal propagation correction
//  X-ray can only be enabled if alignment (cf. "detector table" entry supra)
// did also take X-ray correction into account.
//  And in any case, is incompatible w/ "VarPitch", cf. infra.
// NOTA BENE: A unique entry for all "STRAW settings".
STRAW settings spacers=NO signal_propagation=YES

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
//end
