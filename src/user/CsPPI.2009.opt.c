C2OPT_HEADER
#if 0

#endif
#if !defined C2OPT_COMPUTERFARM
#  define C2OPT_COMPUTERFARM 1
#endif
#if   C2OPT_COMPUTERFARM == 1
//		CsPPI on 2009 hadron data @ CERN.
#elif C2OPT_COMPUTERFARM == 2
//		CsPPI on 2009 hadron data @ gridKa.
#elif C2OPT_COMPUTERFARM == 3
//		CsPPI on 2009 hadron data @ Lyon.
#elif C2OPT_COMPUTERFARM == 4
//		CsPPI on 2009 hadron data @ Compass ONLINE.
#else
//		CsPPI on 2009 hadron data
#endif


#ifndef C2OPT_BARE_FILE
//                                         ***** PHAST OUTPUT
// - Active only if used on a PHAST-compatible executable, cf. "./README"
// - Cf. "$PHAST/coral/README" for the details of the options
mDST file		/dev/null
#endif

//					   ***** HISTOGRAMMING
histograms package	ROOT
histograms home		/dev/null

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
Data file  /castor/cern.ch/compass/data/2009/raw/W45/cdr22002-81912.raw
#  elif C2OPT_COMPUTERFARM == 2
Data file 	/grid/fzk.de/compass/compass2/data/XXX.raw
#  elif C2OPT_COMPUTERFARM == 3
Data file	cchpsscompass:/hpss/in2p3.fr/group/compass/data/XXX.raw
#  elif C2OPT_COMPUTERFARM == 4
Data file	/castor/cern.ch/compass/data/2009/raw/W24/cdr26014-74892.raw
#  endif
#endif


// 					   ***** DETECTOR TABLE...
//  A full directory (path name terminated by a slash "/") or the full path to
// a particular file can be specified.
detector table	$COMPASS_FILES/geometry/2009/


// 					   ***** MATERIAL MAPS
CsMaterialMap	Zone_1		$COMPASS_FILES/maps/material/2008.04/material.target.h2008.04.dat
CsMaterialMap	Zone_2		$COMPASS_FILES/maps/material/2008.04/material.saclay.h2008.04.map
CsMaterialMap	Zone_3		$COMPASS_FILES/maps/material/2008.04/material.straw.h2008.04.map
CsMaterialMap	Zone_4		$COMPASS_FILES/maps/material/2008.04/material.rich.h2008.04.map
CsMaterialMap	Zone_5		$COMPASS_FILES/maps/material/2008.04/material.hcal1.h2008.04.dat
CsMaterialMap	Zone_6		$COMPASS_FILES/maps/material/2008.04/material.muf1.h2008.04.dat
CsMaterialMap	Zone_7		$COMPASS_FILES/maps/material/2008.04/material.muf2.h2008.04.map
CsMaterialMap	Zone_8		$COMPASS_FILES/maps/material/2008.04/material.muf3.h2008.04.map
CsMaterialMap	ELossTarget	$COMPASS_FILES/maps/material/2008.04/material.target_de.h2008.04.dat

//					   ***** CALIBRATION
#if C2OPT_COMPUTERFARM > 4
use calibration
CDB use  	FileDB
CDB location  	$COMPASS_FILES/calibrations/2009
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
CDB server	pccodb01
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
decoding map	$COMPASS_FILES/maps/2009.xml


//					   ***** EVENT SELECTION
//events to read 30000			// # of events to read
//events to skip 1945//41049-78172#116438559 //714//41049-78172#116413939 473//41049-78172#116409119 443//41049-78172#116408519 473//116409119 //2744//41049-78172#116454539
// BoS Veto. 
//events BOS_skip 1.151	// Status of hardware veto unknown as of 09/12 => Keep software veto commented out.

//selection trigger mask f		// Trigger selection (hexadecimal)
selection trigger mask 7ff		// Trigger selection (hexadecimal)
//selection trigger strict


// seed for random number generations:
random number engine JamesEngine
random number seed 19990102
//reset random seed every new event


// 					   ***** DECODING
make decoding 	 	// <nothing>, MCExact



CsPPI_EC02time_DB       host  wwwcompass
CsPPI_EC02time_DB       db    runlb
CsPPI_EC02time_DB       user  anonymous
CsPPI_EC02time_DB       pass  
CsPPI_EC02time_DB       table_tcalib  tb_ECALref
CsPPI_EC02time_TIMEpar  output_folder /afs/cern.ch/compass/scratch/weekly/d03/ppi
//
CsPPI_EC02time_TIMEpar time_hist_width 5
CsPPI_EC02time_TIMEpar cfd_threshold 20
CsPPI_EC02time_TIMEpar cfd_delay 2
CsPPI_EC02time_TIMEpar cfd_ampl 2
CsPPI_EC02time_TIMEpar baseline 50
//




//					   ***** RECONSTRUCTION SCHEMA
reconstruction schema 0




//  	  	  	  	  	   ***** MAGNETIC FIELD MAPS
CsField SM1m_field_measured $COMPASS_FILES/maps/mag_fields/SM1m/SM1M.map.172.data
CsField SM2_field	 $COMPASS_FILES/maps/mag_fields/SM2/FSM.map.5000.data
#ifdef C2OPT_BARE_FILE
// keep this line for the bare option files used for example in alignment
// 4000A is for example used during DVCS
// CsField SM2_field	 $COMPASS_FILES/maps/mag_fields/SM2/FSM.map.4000.data
#endif

//					   ***** MAGNETs INFO
#if C2OPT_COMPUTERFARM <= 4
#else
CsMagInfo	File	~ybedfer/public/magInfo.2007.txt	// ...if MySQL not avail
#endif
CsMagInfo	SM2	0 	// =0: No NMR. !=0: Use NMR w/ numerical parameter as rescaling factor. (Note: From: Catarina Marques Quintans. Sent: Thu 2/5/2009 7:02 PM.  The DCS group submitted to the mySQL database the SM2 NMR readings for the 2008 data taking. But the evaluation of their impact is not clear: in order to enable them, we would probably need to specify a non-unity scaling factor. Yet to be defined...)
CsMagInfo	MySQLDB
//CsMagInfo	Update solenoid	// To be enabled during target rotation



//					   ***** LOGGER
error logger log level		error	// debugging, verbose, info, anomaly, warning, error, fatal
error logger store level	none
error logger verbosity		normal	// low, normal, high. 

//	  	  	  	  	   ***** TraFDic (+DET OPTIONS+...)


//		***** OVERWRITE WHAT'S SET IN "../pkopt/trafdic.2009.opt" *****

// Specify the beam charge, ID and momentum (Will be X-checked against magnet polarities from the detector table)
TraF	iCut	[15]	-1	// Beam charge.
TraF	iCut	[30]	1	// Beam particle ID: 1=hadron.
TraF	dCut	[ 4]	190	// Beam momentum (GeV).



// ===========================================================================
//end
