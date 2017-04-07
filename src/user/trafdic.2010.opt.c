C2OPT_HEADER
#if 0

// $Id: trafdic.2010.opt.c,v 1.21 2011/02/15 23:01:36 ybedfer Exp $

// coral options c file for TraFDic on 2010 data
//
// NOTA.BENE.:c file is to be processed with:
//
//        "gcc -E -P -C -D<defines>"
//
// before being supplied to coral. 

// From "./trafdic.2007.opt.c,v1.40"
// Modified...
// - Transversity is the default.
// - Setup and mapping 2007 -> 2010. 
// - 2007 entries otherwise retained for most year dependent options, but DetNameOff, RICH and calo options taking into account later developments.
// $Log: trafdic.2010.opt.c,v $
// Revision 1.21  2011/02/15 23:01:36  ybedfer
//  Cancel CDBentrytime restriction on FIs, following the updating the CDB
// done by Christopher.
//
// Revision 1.20  2010/12/31 14:37:18  suhl
// * remove obsolete option 'mDST digits', replace it with 'mDST Digits' and 'mDST DAQdigits' in the template options file
//
// Revision 1.19  2010/12/08 13:58:05  ybedfer
//  Update of the block en entries for histo'ing GEM Amplitude Correlations.
//
// Revision 1.18  2010/12/06 03:15:10  ybedfer
//  Put a temporary CDBentrytime on FI: latest calibs, entered in 2001/11,
// are bugged.
//
// Revision 1.16  2010/12/02 10:48:41  fsozzi
// Reference to new RICH option file src/pkopt/rich1.2010.opt
//
// Revision 1.15  2010/11/30 22:44:23  ybedfer
//  New BMS<->beamT tuning, adjusted on 10W27 data.
//
// Revision 1.14  2010/11/24 14:21:48  tnagel
// * remove CsTimingInfo from trafdic.*.opt.c (it's not needed there)
// * updated CsTimingInfo options in pkopt/calorim_2009*.opt and added some documentation
//
// Revision 1.11  2010/10/29 14:16:55  ybedfer
//  - DC04Y.
//
// Revision 1.9  2010/10/20 16:21:52  fsozzi
// update rich option
//
// Revision 1.8  2010/08/12 18:05:42  ybedfer
//  "CsTimingInfo":
//   - Disabled by default.
//   - Yet, in the block "TEMPLATE FOR MASS PRODUCTION", it's reminded that
//    one should try and enable it for any mass production.
//
// Revision 1.7  2010/07/23 13:47:48  ybedfer
//  Reference to the newly created "../pkopt/trigger.2010.opt".
//
// Revision 1.6  2010/06/10 01:03:15  ybedfer
//  No longer specify explicit detectors.dat and dico files in the transverse
// case: coral is able to get along w/o them even in that case.
//
// Revision 1.5  2010/06/01 00:50:54  ybedfer
//  - Default input data = 10T19/cdr22012-84473.raw.
//  - Backtrack'ing coeffs added, copied from 2007.
//  - NMR enabled.
//
// Revision 1.4  2010/05/11 21:19:48  ybedfer
//  Reference to newly created "../pkopt/trafdic.2010.opt".
//
// Revision 1.1  2010/05/04 16:21:29  ybedfer
//  Initial version. Before dataking even started. Not finalised.
//

#endif
#if !defined C2OPT_COMPUTERFARM
#  define C2OPT_COMPUTERFARM 1
#endif
#if   C2OPT_COMPUTERFARM == 1
// TraFDic on 2010 data @ CERN.
#elif C2OPT_COMPUTERFARM == 2
// TraFDic on 2010 data @ gridKa.
#elif C2OPT_COMPUTERFARM == 3
// TraFDic on 2010 data @ Lyon.
#elif C2OPT_COMPUTERFARM == 4
// TraFDic on 2010 data @ Compass ONLINE.
#else
// TraFDic on 2010 data.
#endif

#ifdef C2OPT_BARE_FILE
//  ``Bare'' options file to be included in "$CORAL/src/alignment/traf.*.opt"
// files for processing 2010 FIELD ON data for alignment purposes:
//  - It lacks input data specification (this must be entered in the main
//   "traf.*.opt" file)...
//  - ...as well as mDST output, all high level (vertex and the like)
//   reconstruction options, etc..
//  - More generally, this file should not be edited: any modification to the
//   option entries it contains must instead be implemented in the main
//   "traf.*.opt" by overwriting.
#else
//  Restricted to magnets ON data.
#  ifdef C2OPT_LONGITUDINAL
//  Special longitudinal.
#  else
//  Both longitudinal and transverse setup are covered. Transverse being
// default, and longitudinal specifics being provided as commented out
// alternatives.
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
#  ifndef C2OPT_LONGITUDINAL
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
mDST select		9	// At least one vertex or particle | First and last helices only.
// - One may want to output more info, as exemplified infra:
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
#  if   C2OPT_COMPUTERFARM == 1 && defined C2OPT_LONGITUDINAL
Data file	/castor/cern.ch/compass/data/2010/raw/T19/cdr22012-84473.raw
#  elif C2OPT_COMPUTERFARM == 1
Data file 	/castor/cern.ch/compass/data/2010/raw/W27/cdr22004-85821.raw
#  elif C2OPT_COMPUTERFARM == 2
Data file 	/grid/fzk.de/compass/compass2/data/raw/2010/cdr22012-84473.raw
#  elif C2OPT_COMPUTERFARM == 3 && defined C2OPT_LONGITUDINAL
Data file	cchpsscompass:/hpss/in2p3.fr/group/compass/data/2010/raw/T19/cdr22012-84473.raw
#  elif C2OPT_COMPUTERFARM == 3
Data file 	cchpsscompass:/hpss/in2p3.fr/group/compass/data/2010/T19/raw/cdr22012-84473.raw
#  elif C2OPT_COMPUTERFARM >= 4
Data file 	/castor/cern.ch/compass/data/2010/raw/W27/cdr22004-85821.raw
#  endif
#endif


// 					   ***** DETECTOR TABLE...
//  A full directory (path name terminated by a slash "/") or the full path to
// a particular file can be specified.
detector table   	$COMPASS_FILES/geometry/2010/
//detector table   	$COMPASS_FILES/geometry/2010/detectors.84847.ECal2CaLib.dat

// 					   ***** MATERIAL MAPS
// (Temporarily take 2007 files)
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
CDB location  	$COMPASS_FILES/calibrations/2010
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
decoding map	$COMPASS_FILES/maps/2010.xml


//					   ***** EVENT SELECTION
events to read 30000			// # of events to read
//events to skip 7
// BoS Veto:
// - It has to be enforced software-wise to ensure data reco be uniform, if forgotten in hardware. (Although could be argued that this could be done in PHAST.)
// - It's commented out here.
//events BOS_skip	1.151
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
// Specify the parameters w/ help of "BackTrack ResBMS/Corr" infra. They tell the resolution of the backtracking at the various BMS planes (only 4 of the planes, used so far, as of 2010/11) and the offset corrections to be applied.
// The values infra were tuned on 10W27#85821, cf. "http://wwwcompass.cern.ch/twiki/bin/view/DataReconstruction/CoralSoftware#Back_propagation_beamTelescope_B".
BeamRecons BackTrack ResBMS 4.74 4.20  2.02   2.11
BeamRecons BackTrack Corr  12.23 5.31 -5.69 -12.81

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
make rich1 reconstruction
include ${CORAL}/src/pkopt/rich1.2010.opt	

CsRICH1UpGrade  CONFIG		// RICH upgrade: 12 APV 4 MAPMT
RICHONE 	PrintDetMirr	0
RICHONE 	AcknoMethods	0	// prints methods in use
RICHONE 	PrintKeys	0	// 1 = print keys read,  0 = No print
RICHONE 	PrintConsts	0	// 1 = print rec. constants read,  0 = No print
RICHONE 	DoThetaLikeMax 	NO 	// Faster but "the output buffer of rich won't have the angle that maximizes the likelihood anymore" (according to Giulia, cf. mail from 2009/06/29)

//					   ***** CALORIMETERS
make	calorimeters	reconstruction
include ${CORAL}/src/pkopt/calorim_2009.opt
// no LED/Laser corrections available for 2010??!
EC01P1__ MoreRecoOptions USE_LED_REF_CALIB_CORRECTION=NO USE_LED_REF_CALIB_CORRECTION_IN_SPILLS=NO
EC02P1__ MoreRecoOptions USE_LED_REF_CALIB_CORRECTION=NO USE_LED_REF_CALIB_CORRECTION_IN_SPILLS=NO

//                                         ** Rich Wall + ECAL1 reconstruction
// (Note: Don't know what to do here, whether to use 2007 or 2008 version. The latter are specified as commented out alternatives.)
make RW reconstruction                 
include ${CORAL}/src/pkopt/rw_2007.opt
//make	RW	charge reconstruction	// RW+calo charge reco
//include	${CORAL}/src/pkopt/rw_ch_2008.opt

//					   ***** MU' PID +  FLUX SCALERS
include ${CORAL}/src/pkopt/trigger.2010.opt
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
#ifdef C2OPT_LONGITUDINAL
CsField SOL_field	 $COMPASS_FILES/maps/mag_fields/SOL/SOL_map_fabrice.dat
#else
// Select entry relevant for longitudinal or transverse, among following 2:
//CsField SOL_field	 $COMPASS_FILES/maps/mag_fields/SOL/SOL_map_fabrice.dat
CsField SOL_field	 $COMPASS_FILES/maps/mag_fields/SOL/OD_dipole.fieldmap
#endif
CsField SM1m_field_measured $COMPASS_FILES/maps/mag_fields/SM1m/SM1M.map.172.data
CsField SM2_field	 $COMPASS_FILES/maps/mag_fields/SM2/FSM.map.4000.data

//					   ***** MAGNETs INFO
#if C2OPT_COMPUTERFARM <= 4
CsMagInfo	MySQLDB
#else
CsMagInfo	File	~ybedfer/public/magInfo.2010.txt	// ...if MySQL not avail
#endif
CsMagInfo	SM2	1		// Do rescale SM2 w/ NMR (retrieved from mySQLDB) \times correcting factor specified as argument to this "CsMagInfo SM2" option.



//					   ***** LOGGER
error logger log level		error	// debugging, verbose, info, anomaly, warning, error, fatal
error logger store level	none
error logger verbosity		normal	// low, normal, high. 

//	  	  	  	  	   ***** TraFDic (+DET OPTIONS+...)

include ${CORAL}/src/pkopt/trafdic.2010.opt

//		***** OVERWRITE WHAT'S SET IN "../pkopt/trafdic.2010.opt" *****
//  Specify full path to some alternative dico file (if needed).
//TraF	Dicofit		$COMPASS_FILES/geometry/2010/dico/dico.84819.ECAL1

// Redefine the TBNAMES of DETECTORS to be EXCLUDED FROM PATTERN RECOGNITION:
// - In any case, some detectors ("VI/VO/DW02X2/DW02Y1") have to be left turned off , cf. explanations in "../pkopt/trafdic.2010.opt".
// - But one may want, or may have, to turn off some more.
//   - "DC04Y1/2" are not working during the last part of 2010 data taking. => Have to be turned off in reconstruction. For mass production, cf. the "DetNameOff" rntry in the block "TEMPLATE FOR MASS PRODUCTION" in fine. Otherwise, uncomment the "DetNameOff" option infra.
// - Also, in coral jobs dedicated to detectors studies, one usually turns off the particular detector under exam. Note that then your redefinition will overwrite all "DetNameOff" lists previously entered: take care to include in it all the detectors that have to be turned off in any case.
//TraF	DetNameOff	VO	VI	BM	DW02Y1	DW02X2	MP	DC04Y

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
// (NOTA BENE: The positions of E/HCAL1 may in fact differ from what is quoted supra and corresponds to the 2006/7 setup.)
TraF	SmoothPos [1]	1079.0	// Smooth track at ECAL1
TraF	SmoothPos [2]	1207.5	// Smooth track at HCAL1
TraF	SmoothPos [3]	3292.7	// Smooth track at ECAL2
TraF	SmoothPos [4]	3510.0	// Smooth track at HCAL2
#endif

//                 STRAW options: X-ray and signal propagation correction
// - X-ray can only be enabled:
//   - if alignment (cf. "detector table" entry supra) did also take X-ray correction into account.
//   - if it's is compatible w/ the "VarPitch" option chosen infra.
// - Assuming T0 per channel calib are available for 2010.
// - NOTA BENE: A unique entry for all "STRAW settings".
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
// (Note: One has to systematically uncomment next 4 entries. For special MegaDST production, one has to uncomment also the 5th one.)
//mDST	select	0	// 0: no selection
//mDST	Digits	EC
//mDST	DAQdigits	RP CE
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
// (Note: One has to uncomment the following entry depending upon the data taking period one is about to produce. It allows to exclude some detectors from the track reconstruction, had they been failing during a, large, part (then production gets more stable) or all (tracking performs better if it is told not to expect hits from a particular detector) of the period. In 2010, DC04Y1/2 have to be so excluded in periods W42 et sqq. (last periods of the year). (Caveat: one has to supply the full list of excluded detectors each time one re-specifies this "TraF DetNameOff" option.))
//TraF	DetNameOff	VO	VI	BM	DW02Y1	DW02X2	MP	DC04Y	// SPECIAL SETTING for PERIOD W42 et sqq. => TO BE UNCOMMENTED ONLY THEN.

//					***** UPDATE PATTERN RECOGNITION options
// (Note: Excluding detectors may imply updating some pattern recognition options. In the present case, there's nothing we can do to make for the loss of DC04Y.)
// ===========================================================================
//end
