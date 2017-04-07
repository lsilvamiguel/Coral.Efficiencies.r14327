C2OPT_HEADER
#if 0

// $Id: trafdic.2007.opt.c,v 1.52 2010/12/31 14:37:18 suhl Exp $

// coral options c file for TraFDic on 2007 data
//
// NOTA.BENE.:c file is to be processed with:
//
//        "gcc -E -P -C -D<defines>"
//
// before being supplied to coral. 

// From "./trafdic.2006.opt.c,v1.28"
// Modified...
// $Log: trafdic.2007.opt.c,v $
// Revision 1.52  2010/12/31 14:37:18  suhl
// * remove obsolete option 'mDST digits', replace it with 'mDST Digits' and 'mDST DAQdigits' in the template options file
//
// Revision 1.51  2010/12/08 13:58:05  ybedfer
//  Update of the block en entries for histo'ing GEM Amplitude Correlations.
//
// Revision 1.50  2010/11/24 14:21:47  tnagel
// * remove CsTimingInfo from trafdic.*.opt.c (it's not needed there)
// * updated CsTimingInfo options in pkopt/calorim_2009*.opt and added some documentation
//
// Revision 1.49  2010/11/22 20:17:44  ybedfer
//  Useless mDST options removed from the "TEMPLATE for MASS PROD" block.
//
// Revision 1.48  2010/11/22 16:33:15  ybedfer
//  Add "DW02X2" to the alternative list of "DetNameOff" (commented out by
// default).
//
// Revision 1.47  2010/11/21 13:10:05  ybedfer
//  Enable re-scaling SM2 w/ NMR. I checked that this yields a re-scaling
// factor close to that stored in detectors.dat for a typical run of the
// period 2007/W39.
//
// Revision 1.46  2010/08/12 18:05:42  ybedfer
//  "CsTimingInfo":
//   - Disabled by default.
//   - Yet, in the block "TEMPLATE FOR MASS PRODUCTION", it's reminded that
//    onr should try and enable it for nay mass production.
//
// Revision 1.45  2010/06/10 01:03:15  ybedfer
//  No longer specify explicit detectors.dat and dico files in the transverse
// case: coral is able to get along w/o them even in that case.
//
// Revision 1.44  2010/06/01 00:51:49  ybedfer
//  Drift detectors: enable "make always two clusters".
//
// Revision 1.43  2010/05/27 15:06:07  ybedfer
//  - New values for Backtracking to BMS: they were tuned on 07/W30/#59891.
//  - Enable T0 in STRAW option.
//
// Revision 1.42  2010/05/26 14:16:20  ybedfer
//  Add temporary values for the Bactracking BMS <- beamTelescope.
//
// Revision 1.41  2010/05/10 10:44:52  ybedfer
//  pccodb01 -> pccodb00.
//
// Revision 1.40  2010/02/07 03:45:14  ybedfer
//  Comment out the "always two clusters" entries, which were introduced by
// error in v1.39: not yet tested.
//
// Revision 1.39  2010/02/07 03:27:17  ybedfer
//  More detailed comments...
//
// Revision 1.38  2009/02/24 14:12:46  ybedfer
//  Template for mass production: add a, commented out, entry for setting the
// # of events to preocess to a large number.
//
// Revision 1.35  2009/02/19 01:46:54  ybedfer
//  Upgrade RICH options file.
//
// Revision 1.34  2009/02/19 01:31:27  ybedfer
//  - Add a block "template for production", to be uncommented/edited in
//   order to produce a file (so-called "template.opt") valid for mass
//   production.
//  - "mDST select" uncommented.
//
// Revision 1.30  2008/06/17 23:25:26  ybedfer
//  "CsAverPattern AcceptTill": Increase it from 95 to 105%. The rationale is
// the same as for allowing vertices w/ BMS, viz. the definition of the best
// vertex has to be done w/o any a priori.
//
// Revision 1.28  2008/03/09 20:09:19  ybedfer
//  "DetNameOff": Recall that, as yet, "HQ01Y1_m" has to be turned off (in
// the commented-out entry suggested for overriding what's set in
// "../pkopt/trafdic.2007.opt").
//
// Revision 1.27  2008/03/01 23:51:41  ybedfer
//  Smoothing point @ calorimeters as in Vladimir's recommendation.
//
// Revision 1.26  2008/03/01 21:13:52  ybedfer
//  Use version "2007.00" of material maps.
//
// Revision 1.23  2008/01/20 18:25:10  ybedfer
//  Bug fix: "CDBentrytime" in ONE word.
//
// Revision 1.22  2008/01/20 18:07:47  ybedfer
//  Max. entry time for DWs calibs.
//
// Revision 1.21  2008/01/08 20:37:58  ybedfer
//  - New default input data in the longitudinal case: 58926.
//  - "CDB entrytime" entry blocking the new calib of DWs until it's
//   evaluated.
//
// Revision 1.19  2007/11/11 18:32:40  ybedfer
//  "BOS_skip" enabled by default for transversity.
//
// Revision 1.18  2007/11/11 01:02:55  ybedfer
//  Reference to 2007 transverse raw data.
//
// Revision 1.17  2007/11/09 22:29:35  ybedfer
//  "CaAverPattern": Time cut now expressed in # of sigmas.
//
// Revision 1.15  2007/10/02 17:07:57  ybedfer
//  Cancel the enabling the re-scaling SM2 by NMR, for it's useless (no NMR
// data available for 2007).
//
// Revision 1.14  2007/10/02 10:04:00  ybedfer
//  - Detector table and dico: defaults paths are set = official dirs.
//  - Cancel the printing RICH reco constants.
//
// Revision 1.13  2007/08/09 09:18:52  ybedfer
//  Comment out entry for BoS rejection: it's implemented hardware-wise since
// 07/09, i.e. middle of W27.
//
// Revision 1.10  2007/07/23 01:22:25  ybedfer
//  - BoS rejection enabled (after BoS source code has been upgraded, cf.
//   "../event/CsEvent.cc").
//  - New default input data.
//  - Geometry files set accordingly.
//  - Less verbose RICH.
//
// Revision 1.9  2007/07/20 16:51:14  ybedfer
//  2007.00 versions of target material maps.
//
// Revision 1.8  2007/07/10 20:22:14  korzenev
// trigger.2007.opt was added
//
// Revision 1.7  2007/06/12 04:13:13  ybedfer
//  - Transversity:
//   - Add a dedicated "C2OPT_TRANSVERSE" macro. W/
//     - Specific detectors.dat and dico files.
//     - Dipole field map.
//  - Mention MM03Y (it's included in an alternative entry, commented out by
//   default) and cancel GM01.
//

#endif
#if !defined C2OPT_COMPUTERFARM
#  define C2OPT_COMPUTERFARM 1
#endif
#if   C2OPT_COMPUTERFARM == 1
// TraFDic on 2007 data @ CERN.
#elif C2OPT_COMPUTERFARM == 2
// TraFDic on 2007 data @ gridKa.
#elif C2OPT_COMPUTERFARM == 3
// TraFDic on 2007 data @ Lyon.
#elif C2OPT_COMPUTERFARM == 4
// TraFDic on 2007 data @ Compass ONLINE.
#else
// TraFDic on 2007 data.
#endif

#ifdef C2OPT_BARE_FILE
//  ``Bare'' options file to be included in "$CORAL/src/alignment/traf.*.opt"
// files for processing 2007 FIELD ON data for alignment purposes:
//  - It lacks input data specification (this must be entered in the main
//   "traf.*.opt" file)...
//  - ...as well as mDST output, all high level (vertex and the like)
//   reconstruction options, etc..
//  - More generally, this file should not be edited: any modification to the
//   option entries it contains must instead be implemented in the main
//   "traf.*.opt" by overwriting.
#else
//  Restricted to magnets ON data.
#  ifdef C2OPT_TRANSVERSE
//  Special transversity.
#  else
//  Both longitudinal and transverse setup are covered. Longitudinal being
// default, and transverse specifics being provided as commented out
// alternatives (cf. "./trafdic.2007.transv.opt" where these are commented in).
#  endif

//  This file sums up the main options. Other options are specified elsewhere
// and referenced here by "include <path_to_2ndary_options_file>" statements.
// All options can be entered several times, whether in main or secondary files.
// The last one, in the stream of entries and "include" statements, always wins.

//  The file is to be edited for:
//        Enabling mass-production specific templates and options, cf. the
//       block entitled "TEMPLATE FOR MASS PRODUCTION" in fine.
// Or otherwise for:
//    I)  Specifying input data one wants to process instead of the default.
//   II)  Setting accordingly:
//       - The full path to the relevant detectors.dat and dico files if they
//        differ from the ones guessed from the default geometry and dico
//        directories by coral's automatic procedure, or if the latter fails.
#  ifndef C2OPT_TRANSVERSE
//       - Target magnetic field, cf. the (inappropriately named)
//        "CsField SOL_field" option infra.
#  endif
//  III)  Trimming the list of turned off detectors (cf. "DetNameOff" entry
//       infra): have to be turned off software-wise all detectors dead (or
//       gravely ill):
//       - for a significant fraction of the run being analysed,
//       - or, when reconstructing in view of an asymmetry analysis, for a
//        significant fraction of a spin configuration and then they have to be
//        turned off for the reconstruction of that period and of that period's
//        counterpart in the asymmetry combined analysis.
//         E.g. SI03X/U.
//   IV)  Redefining the # of events to process, cf. "events to read" infra.
//    V)  Adding whatever options are needed to enable a given output:
//      - mDST (enabled by default, but active only if executable is
//       PHAST-compatible, cf. "./README").
//        Default mDST selection is:
//        - Event with at least one vertex or one particle.
//        - Only first and last helix of every track is stored.
//       Cf. $PHAST/coral/README for how to change this behaviour.
//      - Write-back of raw data.
//   VI)  Selection criteria: Trigger selection, BoS rejection, etc...
//  VII)  CDB entrytime.
// VIII)  Disabling X-ray correction: cf. "STRAW" entry infra. It's enabled by
//      default. Should only be disabled if alignment is NOT X-ray corrected.
//       => Check "detector table".
//   IX)  Enabling VarPitch correction: cf. "VarPitch" entry infra.
//    X)  Enabling the updating the solenoid field, if processing a target
//      rotation run.
//   XI)  Using alternative MM resolution (cf. MM description in detector table)
//      and clusterisation cuts.

// Tip: Set your "emacs" in "c++-mode" to view this ".opt" file. And have this
//     done automatically by setting the appropriate hook in your "~/.emacs.el"
//     init file, cf. "~ybedfer/.emacs.el".
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
//CsKalmanFitting WriteBack	evtDump.raw	Physics		// As of 07/01: D*, excl. phi, K0, J/psi
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
#  if   C2OPT_COMPUTERFARM == 1 && defined C2OPT_TRANSVERSE
Data file	/castor/cern.ch/compass/data/2007/raw/W27/cdr32001-58926.raw
#  elif C2OPT_COMPUTERFARM == 1
Data file 	/castor/cern.ch/compass/data/2007/raw/W33/cdr27002-61061.raw
#  elif C2OPT_COMPUTERFARM == 2
Data file 	/grid/fzk.de/compass/compass2/data/raw/2007/cdr22001-56435.raw
#  elif C2OPT_COMPUTERFARM == 3 && defined C2OPT_TRANSVERSE
Data file	cchpsscompass:/hpss/in2p3.fr/group/compass/data/2007/raw/W27/cdr32001-58926.raw
#  elif C2OPT_COMPUTERFARM == 3
Data file 	cchpsscompass:/hpss/in2p3.fr/group/compass/data/2007/tmp/raw/cdr22001-56435.raw
#  elif C2OPT_COMPUTERFARM >= 4
Data file	/castor/cern.ch/compass/data/2007/raw/W27/cdr32001-58926.raw
#  endif
#endif


// 					   ***** DETECTOR TABLE...
//  A full directory (path name terminated by a slash "/") or the full path to
// a particular file can be specified.
detector table   	$COMPASS_FILES/geometry/2007/
//detector table   	$COMPASS_FILES/geometry/2007/detectors.56285.2.5T.plus.dat

// 					   ***** MATERIAL MAPS
CsMaterialMap Zone_1 $COMPASS_FILES/maps/material/2007.00/material.target.m2007.00.dat
CsMaterialMap Zone_2 $COMPASS_FILES/maps/material/2007.00/material.saclay.m2007.00.map
CsMaterialMap Zone_3 $COMPASS_FILES/maps/material/2007.00/material.straw.m2007.00.map
CsMaterialMap Zone_4 $COMPASS_FILES/maps/material/2007.00/material.muf2.m2007.00.map
CsMaterialMap Zone_5 $COMPASS_FILES/maps/material/2007.00/material.muf3.m2007.00.map
CsMaterialMap Zone_6 $COMPASS_FILES/maps/material/2007.00/material.hcal1.m2007.00.map
CsMaterialMap Zone_7 $COMPASS_FILES/maps/material/2007.00/material.muf1.m2007.00.map
CsMaterialMap Zone_8 $COMPASS_FILES/maps/material/2007.00/material.rich.m2007.00.map
CsMaterialMap ELossTarget $COMPASS_FILES/maps/material/2007.00/material.target_de.m2007.00.dat

//					   ***** CALIBRATION
#if C2OPT_COMPUTERFARM > 4
use calibration
CDB use  	FileDB
CDB location  	$COMPASS_FILES/calibrations/2007
#else
use calibration
CDB use 	MySQLDB
#  if   C2OPT_COMPUTERFARM == 1
CDB server	wwwcompass	// This is, as of 07/05, pointing to "lxfsrb6103" located in computer center (reachable outside Cern)
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
DW	CDBentrytime	2007-12-01-00:00:00	// Disregarding the calibs entered on Dec. 14th, 2007
// CDB can be replaced by a File DB, for any given (set of) TBname(s):
#  if   C2OPT_COMPUTERFARM == 1
//DW FileDB	/afs/cern.ch/user/l/leberig/public/marcin/
#  else
//DW FileDB	$COMPASS_FILES/fileDBs/DW/marcin
#  endif
#endif

//					   ***** MAPPING
decoding map	$COMPASS_FILES/maps/2007.xml


//					   ***** EVENT SELECTION
events to read 30000			// # of events to read
//events to skip 7//1048661 //38//1049002 //75//1049607 //39//1049013 //111//1049805 //272//161520181 //199//161519378 //843 //139
// BoS Veto. Cf. Heiner mail from 07/12 ms:
//  - The hardware value is determined to be = 1.151s w.r.t. EWE.
//  - It had not been initially (before 07/09, i.e. middle of W27) applied to
//   all trigger inputs.
//  - Therefore it has to be enforced software-wise to ensure data taking is
//   uniform in the W27+W28 batch of transversity.
#ifdef C2OPT_TRANSVERSE
//  - It's done here for all of 2007 transversity.
events BOS_skip	1.151
#else
//  - It's commented out here.
//events BOS_skip	1.151
#endif
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
//                   Overwrite what's set in "beam_2004.opt"
BeamRecons	useTRAFFIC	1	// >0: Traffic is used, =0: standalone beam package is used
BeamRecons	AllowMultCand	2
PID_doBeamID	selAllBeams	1	// Allow beams w/o BMS in the vertexing
//			Bactracking BMS <- beamTelescope
// Specify the parameters w/ help of "BackTrack ResBMS/Corr" infra. They tell the resolution of the backtracking at the various BMS planes (only 4 of the planes, used so far, as of 2010/04) and the offset corrections to be applied.
//BeamRecons BackTrack Default		// No built-in defaults available for 2007.
BeamRecons BackTrack ResBMS  5.31 3.80  2.29   1.92	// Tuned on 07W30#59891
BeamRecons BackTrack Corr   14.39 6.92 -5.25 -13.01	// Tuned on 07W30#59891

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
CsAverPattern	CUTS	1600	10
CsAverPattern	Refit	1	// 0 - Off, 1 - ON
CsAverPattern	Retrack	3	// 0 - Off, !0 = cut on event time
CsKalmanFitting	Hist	1	// 0 - Off, 1 - ON
CsKalmanFitting	Print [ 0 ]	0	// Kalman filter info.
CsKalmanFitting	Print [ 1 ]	0	// Global fit info.
CsKalmanFitting	Print [ 2 ]	0	// IKF info.
CsKalmanFitting	RefPlane	500	// in mm
CsKalmanFitting	RefBeam	-1500	// in mm
CsKalmanFitting	RefMargin	500	// in mm
// Chi2 increment cuts:		DirChi2Trk	InvChi2Trk
CsKalmanFitting	CUTS		7		10.75
// Vertex Chi2 cut:
CsKalmanFitting	Chi2VCut	10.5
CsKalmanFitting BeamP0		160	// Beam reference momentum (for histo purposes)


//					   ***** RICH1
make	rich1	reconstruction
include ${CORAL}/src/pkopt/rich1.2007.opt
CsRICH1UpGrade  CONFIG		// RICH upgrade: 12 APV 4 MAPMT
RICHONE 	PrintDetMirr	0
RICHONE 	AcknoMethods	0	// prints methods in use
RICHONE 	PrintKeys	0	// 1 = print keys read,  0 = No print
RICHONE 	PrintConsts	0	// 1 = print rec. constants read,  0 = No print
RICHONE 	DoThetaLikeMax 	NO 	// Faster but "the output buffer of rich won't have the angle that maximizes the likelihood anymore" (according to Giulia, cf. mail from 2009/06/29)

//					   ***** CALORIMETERS
make	calorimeters	reconstruction
include ${CORAL}/src/pkopt/calorim.opt

//                                         ** Rich Wall + ECAL1 reconstruction
make	RW	reconstruction                 
include ${CORAL}/src/pkopt/rw_2007.opt

//					   ***** MU' PID
include ${CORAL}/src/pkopt/trigger.2007.opt
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
#ifdef C2OPT_TRANSVERSE
CsField SOL_field	 $COMPASS_FILES/maps/mag_fields/SOL/OD_dipole.fieldmap
#else
// Select entry relevant for longitudinal or transverse, among following 2:
CsField SOL_field	 $COMPASS_FILES/maps/mag_fields/SOL/SOL_map_fabrice.dat
//CsField SOL_field	 $COMPASS_FILES/maps/mag_fields/SOL/OD_dipole.fieldmap
#endif
CsField SM1m_field_measured $COMPASS_FILES/maps/mag_fields/SM1m/SM1M.map.172.data
CsField SM2_field	 $COMPASS_FILES/maps/mag_fields/SM2/FSM.map.4000.data

//					   ***** MAGNETs INFO
#if C2OPT_COMPUTERFARM <= 4
CsMagInfo	MySQLDB
#else
//CsMagInfo	MySQLDB
CsMagInfo	File	/afs/cern.ch/user/y/ybedfer/public/magInfo.2007.txt	// ...if MySQL not avail
#endif
CsMagInfo	SM2	1	// Get NMR from MySQLDB and re-scale SM2 w/ NMR / NMR_map * correction factor. The correction factor is here set =1. This yields, for an average run of 2007/W39, a re-scaling of 1.0098, i.e. slightly higher than the 1.0085 that's typically specified in detectors.dat's (and was derived for a 2006 period).
//CsMagInfo	Update solenoid	// To be enabled during target rotation



//					   ***** LOGGER
error logger log level		error	// debugging, verbose, info, anomaly, warning, error, fatal
error logger store level	none
error logger verbosity		normal	// low, normal, high. 

//	  	  	  	  	   ***** TraFDic (+DET OPTIONS+...)

include ${CORAL}/src/pkopt/trafdic.2007.opt

//		***** OVERWRITE WHAT'S SET IN "../pkopt/trafdic.2007.opt" *****
//  Specify full path to some alternative dico file (if needed).
//TraF	Dicofit		$COMPASS_FILES/geometry/2007/dico/dico.56285.2.5T.plus

// Redefine the TBNAMES of DETECTORS to be EXCLUDED FROM PATTERN RECOGNITION:
// (Caveat: one has to supply the full list of excluded detectors each time one re-specifies this "TraF DetNameOff" option.)
// - In any case, some detectors ("VI/VO/DC04Y/HQ01Y1") have to be left turned off , cf. explanations in "../pkopt/trafdic.2007.opt".
// - MM03Y was dead from 2007/06/06 until 2007/06/12. It has to be turned off software-wise for the whole of that part of W23 that was in transverse mode. (How early does this start is not clear.) 
// - DW02X2 is operated at a lower HV (for whatever reason), since at least W35#61459, cf. "wwwcompass.cern.ch/runLogbook/dirphp/show_histrun.php?runnb=61459&cfgnb=0".
// - Note that some detectors are ignored altogether, cf. "Det2Ignore" in "../pkopt/trafdic.2007.opt".
//TraF	DetNameOff	VO	VI	DC04Y	HQ01Y1	DW02X2	MM03Y

//TraF	ReMode	[17]	0	// Disable GEM amplitude correlations, if any GEM is among the turned-off detectors.

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
TraF	SmoothPos [1]	1079.0	// Smooth track at ECAL1
TraF	SmoothPos [2]	1207.5	// Smooth track at HCAL1
TraF	SmoothPos [3]	3292.7	// Smooth track at ECAL2
TraF	SmoothPos [4]	3510.0	// Smooth track at HCAL2
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
// (Note: One has to systematically uncomment the next entry. For special MegaDST production, one has to uncomment also the 2nd one.)
//mDST	selectTrigger	3072	// Random trigger: unconditional output (Note: the option is not useful when "mDST select 0", i.e. all-out unconditional output.)

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
// (Note: One has to freeze the setup info accessible to the production, so that the latter be independent of the modifications later brought to the setup files, and hence stable. This is achieved by copying the setup files to sub-directories of the production directory. E.g. infra are listed what need be done for the production of 2008W37/slot2.) 
//decoding	map	/afs/cern.ch/compass/scratch/d17/objsrvvy/2008/W37/slot2/maps
//detector	table	/afs/cern.ch/compass/scratch/d17/objsrvvy/2008/W37/slot2/alignement/
//TraF		Dicofit	/afs/cern.ch/compass/scratch/d17/objsrvvy/2008/W37/slot2/dico/


//					***** DETECTORS EXCLUDED FROM TRACKING
// (Note: One has to uncomment the following entry depending upon the data taking period (or transversity super-period) one is about to produce. It allows to exclude some detectors from the track reconstruction, had they been failing during a, large, part (production is then more stable) or all (tracking performs better if it is told which are the bad guys and not to expect hits from them) of the period.
// In 2007, MM03Y has to be so excluded in period W23, DW02X2 in periods W35 et sqq. In any case, care to treat equally the 2 legs of a transversity super-period.)
// (Caveat: one has to supply the full list of excluded detectors each time one re-specifies this "TraF DetNameOff" option.)
//TraF	DetNameOff	VO	VI	DC04Y	HQ01Y1	DW02X2	MM03Y

// ===========================================================================
//end
