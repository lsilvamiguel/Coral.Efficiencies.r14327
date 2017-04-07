// $Id: CsMicroMegaDetector.cc,v 1.96 2010/02/11 15:31:19 suhl Exp $

/*!
   \file    CsMicroMegaDetector.cc
   \brief   Compass MicroMega like detector Class.
   \author  Benigno Gobbo
   \version $Revision: 1.96 $
   \date    $Date: 2010/02/11 15:31:19 $
*/

// $Log: CsMicroMegaDetector.cc,v $
// Revision 1.96  2010/02/11 15:31:19  suhl
// removed unnecessary include of CLHEP random number classes
//
// Revision 1.95  2010/01/28 12:51:25  tnagel
// dropped gcc 2.x support (removing ~1700 lines of code)
//
// Revision 1.94  2010/01/24 16:10:41  suhl
// * replacing "unsigned int" or "int" used to store timestamps with "time_t"
//   to fix the strict-aliasing compiler warning
//
// Revision 1.93  2010/01/22 15:15:14  tnagel
// create symbolic link in include dir for CsComgNtCommons.h
//
// Revision 1.91  2008/11/12 20:57:39  ybedfer
//  - MC Trigger offset subtracted from hit time and "extraTimeWidth" alloted
//   to the time windows.
//  - Relative jitter between CsHit's originating from a same CsMCDigit set
//   equal to zero, for simplicity's sake (was randomly distributed according
//   to absolute detector resolution before).
//
// Revision 1.90  2008/01/08 18:10:45  ybedfer
//  - Major change: Time cut:
//    - The hit time gate is now defined in ns (previously was F1 unit).
//    - It is applied in 2 steps: firstly in decoding, secondly in
//     clustering, after the master trigger has been made available and used
//     to redefine the extra-time accommodating the trigger jitter (and to
//     define the, possible, trigger offset).
//    - The offset of the alternative master trigger is taken into account if
//     defined.
//    - Default hit time gate is 1.5 times cluster time gate (Same factor
//     when deriving hit gate from off-time cluster gate).
//
// Revision 1.89  2007/12/11 13:54:44  ybedfer
//  In block conditioned by "uM_KEEP_CLOSEST_HIT", i.e. retaining the CsDigit
// closest to T0: a "fabs" was forgotten in the comparison of the time of
// new candidate hits w/ that of best candidate, so that <0 times tended to
// be retained in priority.
//
// Revision 1.88  2007/03/12 10:24:00  ybedfer
//  Time resolution: input via option (can be made to vary from one plane to
// the next), accessor. (Note: In the MC simulation: a built-in of 9 ns is
// still used.)
//
// Revision 1.87  2006/12/06 05:55:52  ybedfer
//  "Pos2Wire" now accepts as argument the raw distance to 1st wire, not
// corrected for "wirD_". This change occured w/ the committing version v1.54
// of "CsDetector.cc", but was not documented in the log. (It's meant to
// restore the symmetry between the "Wire2Pos" and "Pos2Wire" functions.)
//
// Revision 1.86  2006/07/11 15:21:21  ybedfer
//  Working around a, putative, bug in "CsEvent::getNextEvent" that fails
// to handle properly decoding troubleshootings. The bug leads to a multiple
// reentering "CsMM::DecodeChipDigit" w/o executing "CsMM::clusterize" in
// the mean time. Which, in turn, leads to "CsMM::oDigits_" keeping track of
// a reference to a no longer existing "CsDigit". The fix, until CsEvent is
// itself fixed, is to clear "oDigits" in "DecodeChipDigit", when it turns
// out that "(CsDet)myDigits_" is empty, instead of (in addition to, rather)
// "clusterize".
//  NOTA BENE: I (Y.B.) haven't investigated the full consequences of the
// fix, so that it may be wiser to disable it when the CsEvent bug is itself
// fixed: en(dis)abling by (un)defining "CsMM_MULTIDECODE_WO_CLUSTERING".
//  The abovementioned bug can be observed in event #82837507 of raw data
// file "03P1H/cdr14008-31520" (one reaches it directly by "events to skip
// 29294"): enable "CsMM_DEBUG_MULTIDECODE" in order to trace the origin
// of the bug.
//
// Revision 1.85  2005/04/16 16:04:32  ybedfer
//  Bug fix: init "cutTriangle_" to 0, so that "TriCut" is called only if
// specified by option!
//
// Revision 1.84  2004/11/20 21:12:59  ybedfer
//  Make use of Colin's "MMLibrary".
//
// Revision 1.83  2003/11/16 19:11:59  ybedfer
//  Extra, trigger dependent, width to hit time window to account for trigger
// jitter.
//
// Revision 1.82  2003/06/08 22:39:46  zvyagin
// Use updated DaqDataDecoding TriggerTime class.
//
// Revision 1.81  2003/04/17 09:20:21  benigno
// --> gcc3
//
// Revision 1.78  2003/03/13 15:41:43  cbernet
// added 	hit crates histograms in MM class.
// 	hit ctoff
//
// hit crates is a profile histogram filled only with hits being inside an off time window
// specified with the following line in the option file
//
// MM01V* HitTimeOff [0-1] 300 800
//
// Revision 1.77  2003/03/07 12:24:36  ybedfer
//  Reintroduce the upper limit on the # of hits to be considered
// for the determination of the centroid (or any other attribute) of
// the cluster. The upper limit is specified by a cpp macro:
// "uM_RETAIN_MAX_HITS". I checked impact on reco perfromances
// on #22019-19006, w/ and w/o limit:
//               Tracks                  Vertex
//      toto   LAS   SAS    mu'     %    #Tracks  mu/mu'
// w/o  3.15  2.54  1.80  0.906   38.7%   3.259    31.0%
// <7   3.16  2.54  1.81  0.906   38.6%   3.261    31.1%
//
// Revision 1.76  2003/02/16 11:45:26  cbernet
// new option :
//
// MM02V* FlagChannels	false/true
//
// when false, channels flagged as bad will be kept. If true or missing, these
// channels are removed at decoding stage
//
// Revision 1.75  2003/02/12 07:59:28  cbernet
// triangle cut in the tot VS t plane can now be defined for MM.
// add this line to the option file :
//
// MM* TriCut  0       0       4       4       400
//
// this cut defines 2 straight lines with slopes 4 and 4.
// 1st 2 numbers are the (t,tot) coordinates of the intersection point
// 400 is tot_max.
//
// this cut is not yet in use.
//
// Revision 1.74  2003/01/30 22:08:21  valexakh
// Storage of hits from charged products added
//
// Revision 1.73  2003/01/23 13:19:42  cbernet
// possibility do define an off trigger time range in the option file :
// MM02* ClusterTime [0-1] -10000 10000
//
// if this cut is passed, new histograms are filled :
// [TBName]_crates		cluster rates
// [TBName]_ctoff         	off trigger cluster time distribution
//
// #ifdef uM_DETECTOR_STUDIES
// is now only used for some cuts. histograms don't depend on defines anymore,
// but on the histogram level
//
// Revision 1.72  2002/12/19 14:27:26  neyret
// Huge clean-ups in calibration code for preparation of the future calib
// MySQLDB
// + some bugs fixed
//
// Revision 1.70  2002/10/22 15:15:39  cbernet
// calibration flags are now supported. Format of calibration files :
// chan#	T0	flag
//
// flag=0 : ok
// flag=1 : noisy (will be removed)
// flag=0 : missing (treated as ok)
//
// - Old calibration files are still accepted.
// - flag 0k (ie 0) is not mandatory.
//
// Revision 1.69  2002/10/18 16:41:15  cbernet
// - only 3 histograms left in "Normal " histogramming level, for quality
// check purposes :
// 	hit chan
// 	cluster t
// 	cluster tot
//
// cluster tot will be replaced by cluster amplitude, for clusters with a
// size > 2 (cf coool).
//
// Revision 1.68  2002/10/18 16:05:56  cbernet
// - created a new 2d histo hit time vs chan for calibration purpose
// (histogramminfg level = high)
//
// - range of histos related to hit time depends on calibration usage.
//
// Revision 1.67  2002/07/01 17:16:46  ybedfer
//  "COMPASS_SETUP == 2001' conditions "uM_PARTIALLY_EQUIPPED".
//
// Revision 1.66  2002/06/26 04:30:29  hpereira
// CsMicroMega: commented some debug printout colin forgot.
//
// Revision 1.65  2002/06/26 00:31:19  cbernet
// - bad channels for 2001, and xlsat parity lead are hidden behind a define
//
// Revision 1.64  2002/06/25 05:02:13  ttoeda
// Move the entry of '#include "CDB.h"' to the .cc file
//
// Revision 1.63  2002/06/17 21:58:09  ybedfer
//  - Replace "CsEvent::getTriggerTime()" by "ChipF1::GetTT().GetTimeNorm()",
//   assuming mM F1 chips have ``normal'' precision. For the conversion to
//   ns: nothing changed.
//  - cpp macro "COMPASS_SETUP" instead of "SETUP_2001".
//
// Revision 1.62  2002/06/13 16:41:31  ybedfer
//  "uM_PARTIALLY_EQUIPPED" conditionned by "SETUP_2001" (and bug fix when its not set).
//
// Revision 1.61  2002/05/24 02:05:37  ybedfer
//  "uM_PARTIALLY_EQUIPPED": properly indented, defined only once...
//
// Revision 1.60  2002/04/22 16:27:54  ybedfer
//  Modifications in the handling of time:
//   - Hit time cut:
//      - from a HitTime range no longer centered on 0,
//      - modifiable by option.
//   - Hit time cut in MC:
//      - according to HitTime range,
//      - on TDC time (instead of MC time).
//   - MC TDC time:
//      - retain TDC closest to trigger (instead of earliest).
//
// Revision 1.59  2002/03/28 14:08:11  ybedfer
//  "calibration are found" -> "CsErrLog(elInfo)".
//
// Revision 1.58  2002/03/15 11:07:19  ttoeda
// Modify readCalibration and DecodeChipDigit
//
// Revision 1.57  2002/03/08 16:08:24  hpereira
// /src/geom/: almost all detector modified to unify histogram switches. Option format is to be found in /pkopt/hist.opt
// CsMicroMegaDetector.h and .cc: more modifications to allow call to bookHistograms after call to calibrationDB. DeadChannels and partiallyEquiped flags are now set in ReadCalib.
// CsDetector.h/.cc: added method ReadHistLevel(), GetHistLevel() and member hLevel.
// GetHistLevel() Returns either CsDetector::None (ie 0) | CsDetector::Normal (ie 1) | CsDetector::High (ie 2) of type CsDetector::histogramLevel (enum) or unsigned int through cast
//
// Revision 1.56  2002/02/08 11:01:58  benigno
// Splitted SciFi dets in J & D types
//
// Revision 1.55  2002/02/01 08:06:33  benigno
// Fixed time set on cluster (on MCs).
//
// Revision 1.54  2002/01/31 16:33:26  benigno
// Added time to MC digit (with sigma of 9 ns)
//
// Revision 1.53  2002/01/22 16:00:12  benigno
// Not cluster setting in CsRecoEvent is centralized in CsEvent...
//
// Revision 1.52  2002/01/21 03:49:41  ybedfer
//  - Access the "CsInit" info about CDB in use instead of scanning the options file
//   every time.
//  - Check for the consistency between the size of the calibration file and the
//   number of wires once and for all, when calibration data are read in.
//  - Apply Hit Time cut consistentl, whether calibration iis from CDB or not.
//  - Recasting of some error messages, using "CsErrLog".
//
// Revision 1.50  2001/12/05 21:34:38  ybedfer
//  - "uM_PARTIALLY_EQUIPPED" enabled.
//
// Revision 1.49  2001/12/05 21:27:36  ybedfer
//  - Major changes:
//    - Code for ``associated'' clusterisation of MC data.
//    - Macro "uM_PARTIALLY_EQUIPPED" (enabled).
//    - Macro "uM_DETECTOR_STUDIES" introducing some cuts on clusters to
//     select golden tracks (disabled).
//  - Also:
//    - Cuts on time (hit, cluster) made non static. Can be modified by
//     options.
//    - Histo's: "htot_" cancelled, replaced by "ctot_", "ct_" added,
//     Edge times, cluster profile conditioned by "uM_DETECTOR_STUDIES".
//    - "xsiz_", "ysiz_" now set by "CsDetector::AddSubDetector".
//    - Options retrieving moved out of "clusterize".
//    - Indent.
//
// Revision 1.46  2001/11/30 08:47:07  benigno
// Now all CsEvent methods start with small letter
//
// Revision 1.45  2001/11/29 18:28:01  bernet
// MC implemented for vaiable pitch MM and hodoscopes.
//
// Revision 1.44  2001/11/22 16:32:36  bernet
//
// New function CsDetector::Wire2Pos(int wire) implemented for clustering
// in variable pitch detectors.
//
// this function is now used also in CsTriggerHodoDetector::clusterize()
//
// new histogram added to CsTriggerHodoDetector : cluster position
//
// Revision 1.43  2001/10/26 12:08:23  ybedfer
//  Several bugs fixed:
//    - Normalisation of ToT when exponentiated was expressed in F1 units:
//     needs ns.
//    - "Wire2Pos" uncorrectly accounting for strips spacing at slices
//     boundaries: introduce alternative method "W2P".
//    - Cut on cluster time moved out of "uM_WEIGHTED_ClTIME" block.
//
// Revision 1.42  2001/10/23 13:41:21  cbernet
// a couple of assert's removed in CsMicroMegaDetector.cc
//
// Revision 1.41  2001/10/18 22:50:05  cbernet
// - cut on cluster time in CsMicroMegaDetector
// - cluster time histo in         "
// - markers removed, replaced by CsHist pointers
//
// Revision 1.40  2001/10/17 18:05:39  ybedfer
//   Optimisation of clusterisation:
//     - Validation of the taking into account up to 5 hits.
//     - Weigthing w/ exp(ToT/1000)
//     - Weigthed average for cluster time.
//

#include "CsMicroMegaDetector.h"

#include <math.h>
#include <string.h>
#include <strings.h>
#include <queue>
#include <functional>
#include "CsInit.h"
#include "CsZone.h"
#include "CsDigit.h"
#include "CsOpt.h"
#include "CsEvent.h"
#include "CsDigit.h"
#include "CsMCDigit.h"
#include "CsComgNtCommons.h"
#include "CsGeom.h"
#include <cstdlib>
#include "CsMCTrkHit.h"
#include "CsRandom.h"
#include "CsGeant3.h"
#include "CsMCTrack.h"
#include "DaqDataDecoding/ChipF1.h"
#include "CDB.h"
//LS Eff
#include <stdexcept>
//#define DEBUG_EFF
#include "CsRandom.h"

using namespace std;
using namespace CLHEP;

extern QhitType Qhit;

const float CsMicroMegaDetector::leadtWght_ = 0.6;      
const float CsMicroMegaDetector::f1Tick_    = 0.12892;  // ns

///////////////////////////////////////////////////////////////////////////////////
CsMicroMegaDetector::CsMicroMegaDetector( const int    row,
					  const int    id,    const char* name,  const char *TBname,
					  const int    unit,  const int    type,
					  const double rdLen, const double xsiz,
					  const double ysiz,  const double zsiz,
					  const double xcm,   const double ycm,
					  const double zcm,   const HepMatrix rotDRS,
					  const HepMatrix rotWRS,
					  const double wirD,  const double ang,
					  const int    nWir,  const double wirP,
					  const double eff,   const double bkg,
					  const double tGate )
  :CsDetector( row, id, name, TBname,unit, type, rdLen, xsiz, ysiz, zsiz,
	       xcm, ycm, zcm, rotDRS, rotWRS, wirD, ang,
	       nWir, wirP, eff, bkg, tGate ),
  parityLead_(0),decodeCard_(false), useCalib_(false), lastChan_(-1),
  digitData_(new double[2]), removeBadChannels_(true) {

  // ********** READ OPTIONS FILE **********
  string tag = "";
  list<string> options; list<string>::iterator Is;
  {
    // Decoding options
    string key = "make decoding";
    if (CsOpt::Instance()->getOpt( tag, key )) {
      if (CsOpt::Instance()->getOpt( tag, key, options )) {
	for( Is=options.begin(); Is!=options.end(); Is++ ) {
	  if( *Is == "MM" || *Is == "MicroMega" || *Is == "all" )
	    decodeCard_ = true;
	}
      }
      else  decodeCard_ = true;
    }
  }

  {
    string MM = "MM";  // tag is "MM" or TB name.

    //   ********** DETECTOR CHARACTERISTICS **********

    cresolution_ = -1;
    string key = "cluster_resolution";
    if( CsOpt::Instance()->getOpt(MM,key,cresolution_) ) {
      CsErrLog::msg(elInfo,__FILE__,__LINE__,"\"%s\" Cluster Resolution = %f",
		    TBname,cresolution_);      
    }

    //     ***** Time resolution *****
    tRes_ = 9;  // Default is 9 ns.
    // This option entry can be made dependent upon detector's TB name, e.g:
    // MM       TimeResolution 9
    // MM01V1__ TimeResolution 12
    key = "TimeResolution"; if (CsOpt::Instance()->getOpt(TBname,key,tRes_) ||
				CsOpt::Instance()->getOpt(MM,key,tRes_))
      CsErrLog::msg(elInfo,__FILE__,__LINE__,"\"%s\" Time Resolution = %f",
		    TBname,tRes_);

    //   ********** CLUSTERISATION OPTIONS **********

    associate_ = true;
    key = "no associated digits";
    if (CsOpt::Instance()->getOpt( tag, key )) {
      if (CsOpt::Instance()->getOpt( tag, key, options )) {
	for( Is=options.begin(); Is!=options.end(); Is++ ) {
	  if( *Is == "MM" || *Is == "MicroMega" || *Is == "all" )
	    associate_ = false;
	}
      }
      else  associate_ = false;
    }

    //       ***** Cluster Splitting 1/0 ****
    key = "Split";
    splitClustersMin_ = splitClustersMax_ = -1;
    vector<double> v;
    if (CsOpt::Instance()->getOpt(MM,key,v) ) {
      if(v.size() != 2) 
	CsErrLog::mes(elFatal,"Error syntax in specifying SplitClusters");
      
      splitClustersMin_ = static_cast<int>(v[0]);
      splitClustersMax_ = static_cast<int>(v[1]);
      CsErrLog::msg(elInfo,__FILE__,__LINE__,"\"%s\" Cluster Splitting ON from channel %d to $d",
		    TBname,splitClustersMin_, splitClustersMax_);
      
    }

    //         ***** CUT ON Cluster Time *****
    // This cut can be made dependent upon detector's TB name, e.g:
    // MM ClusterTime [0-1] -20 20
    // MM01V1__ ClusterTime [0-1] -75 85
    hitTCut_ = 500;
    clTMin_ = -40; clTMax_ = 40;         // Default: -40<T<40
    key = "ClusterTime";
    v.clear();
    if (CsOpt::Instance()->getOpt(TBname,key,v) ||
	CsOpt::Instance()->getOpt(MM,key,v)) {
      if (v.size()!=2)
	CsErrLog::mes(elFatal,"Error syntax in specifying Cluster Time cut");
      clTMin_ = v[0]; clTMax_ = v[1];
      CsErrLog::msg(elInfo,__FILE__,__LINE__,"\"%s\" %f < ClusterTime < %f",
		    TBname,clTMin_,clTMax_);
      hitTCut_ = fabs(clTMin_)<fabs(clTMax_) ?
	fabs(clTMax_)*1.5 : fabs(clTMin_)*1.5;
      if (hitTCut_<500) hitTCut_ = 500;
    }

    //     ***** CUT ON Cluster Time to get off time clusters *****
    cltoffMin_ = clTMax_; cltoffMax_ = 2*clTMax_; 
    key = "ClusterTimeOff";
    v.clear();
    if (CsOpt::Instance()->getOpt(TBname,key,v) || 
	CsOpt::Instance()->getOpt(MM,key,v)) {
      if (v.size()!=2)
	CsErrLog::mes(elFatal,"Error syntax in specifying Cluster Time cut");
      cltoffMin_ = v[0]; cltoffMax_ = v[1];
      CsErrLog::msg(elInfo,__FILE__,__LINE__,"\"%s\" %f < ClusterTimeOff < %f",
		    TBname,cltoffMin_,cltoffMax_);
      
      // option has been specified. User definitely wants off-time clusters
      // for rate measurements -> cut on hit times must be opened. 
      hitTCut_ = fabs(cltoffMin_)<fabs(cltoffMax_) ?
	fabs(cltoffMax_)*1.5 : fabs(cltoffMin_)*1.5;
      hitTCut_ *= 2; // we want something OUTSIDE the cuts
    }

    if (CsInit::Instance()->IsAMonteCarloJob()) hitTCut_ = tGate_/2;
    hitTMin_ = -hitTCut_; hitTMax_ = hitTCut_;


    //             ***** CUT ON Hit Time *****
    // This cut can be made dependent upon detector's TB name, cf. supra:
    key = "HitTime";
    v.clear();
    if (CsOpt::Instance()->getOpt(TBname,key,v) ||
	CsOpt::Instance()->getOpt(MM,key,v)) {
      if (v.size()!=2)
	CsErrLog::mes(elFatal,
		      "Error syntax in specifying Hit Time cut");
      hitTMin_ = v[0]; hitTMax_ = v[1];
      if (hitTMin_<-tGate_/2 || hitTMax_>tGate_/2)
	CsErrLog::msg(elWarning,__FILE__,__LINE__,
  "\"%s\" HitTime window [%f,%f] > tGate %f",TBname,hitTMin_,hitTMax_,tGate_);
      else
	CsErrLog::msg(elInfo,__FILE__,__LINE__,"\"%s\" %f < HitTime < %f",
		      TBname,hitTMin_,hitTMax_);
    }


    //       ***** CUT ON Hit Time to get off time clusters *****
    hitTOffMin_ = hitTMax_; hitTOffMax_ = 2*hitTMax_; 
    key = "HitTimeOff";
    v.clear();
    if (CsOpt::Instance()->getOpt(TBname,key,v) || 
	CsOpt::Instance()->getOpt(MM,key,v)) {
      if (v.size()!=2)
	CsErrLog::mes(elFatal,"Error syntax in specifying HitTimeOff cut");
      hitTOffMin_ = v[0]; hitTOffMax_ = v[1];
      CsErrLog::msg(elInfo,__FILE__,__LINE__,"\"%s\" %f < HitTimeOff < %f",
		    TBname,hitTOffMin_,hitTOffMax_);
    }


    //          ***** 2D CUT ON ToT vs Time *****
   
    key = "TriCut";
    v.clear();
    cutTriangle_ = NULL;
    if (CsOpt::Instance()->getOpt(TBname,key,v) ) {
      if (v.size()!=5)
	CsErrLog::mes(elFatal,
		      "Syntax error in specifying triangle cut");
      cutTriangle_ = new MM::ECutTriangle(v);
      CsErrLog::msg(elInfo,__FILE__,__LINE__,"\"%s\" Triangle cut created : %f %f %f %f %f",
		    TBname, v[0], v[1], v[2], v[3], v[4]);
    }
   
    //#define uM_DETECTOR_STUDIES
#ifdef uM_DETECTOR_STUDIES
    // This macro allows to set strict conditions on any given uM from
    // entries in the options file: e.g.
    // MM ToT 90
    // MM NClusters 2
    // MM ClusterSize 6
    // MM01V1__ ToT 0
    // MM01V1__ NClusters 32
    // MM01V1__ ClusterSize 32
    // To be enabled for detector studies dedicated analyses.

    // ***** CUT ON ToT
    ToTMin_ = 0;                         // Default is 0 <-> No cut
    key = "ToT";
    if (CsOpt::Instance()->getOpt(TBname,key,ToTMin_) ||
	CsOpt::Instance()->getOpt(MM,key,ToTMin_)) {
      CsErrLog::msg(elWarning,__FILE__,__LINE__,"\"%s\" Cut ToT > %f",
		    TBname,ToTMin_);
    }
    //#  define uM_SINGLE_TRACK
#  ifdef uM_SINGLE_TRACK
    //         ***** CUT ON MULTIPLICITY *****
    nclMax_ = 32;
    key = "NClusters";
    if (CsOpt::Instance()->getOpt(TBname,key,nclMax_) ||
	CsOpt::Instance()->getOpt(MM,key,nclMax_)) {
      CsErrLog::msg(elWarning,__FILE__,__LINE__,"\"%s\" NClusters < %d",
		    TBname,nclMax_);
    }
#  endif
    //        ***** CUT ON Cluster Size *****
    clsMax_ = 32;
    key = "ClusterSize";
    if (CsOpt::Instance()->getOpt(TBname,key,clsMax_) ||
	CsOpt::Instance()->getOpt(MM,key,clsMax_)) {
      CsErrLog::msg(elWarning,__FILE__,__LINE__,"\"%s\" ClusterSize < %d",
		    TBname,clsMax_);
    }
#endif
  }

  // calibrations
  useCalib_ = CsInit::Instance()->useCalibration();
  if (CsInit::Instance()->IsAMonteCarloJob())  useCalib_ = false;

  //LS Eff
  //                                          ***** eff. MAP ENABLED ? *****
  mcEffMapsEnabled_ = false;
  if (CsOpt::Instance()->getOpt("ALL"      ,            "mcEffMapsEnabled") ||
      CsOpt::Instance()->getOpt(GetTBName(),            "mcEffMapsEnabled") ||
      CsOpt::Instance()->getOpt(GetTBName().substr(0,2),"mcEffMapsEnabled")) {
    mcEffMapsEnabled_ = true;
#ifdef DEBUG_EFF
    cout<<"**** LS *** : CsMicroMegaDetector:: "<< mcEffMapsEnabled_<<endl;  
#endif
  }
  else CsErrLog::msg(elWarning,__FILE__,__LINE__,
                     "%s: eff. Map will NOT be used",TBname);

  // disable bad channels flagging
  string flag;
  if (CsOpt::Instance()->getOpt(TBname,"FlagChannels",flag)) {
    if(flag == "true") {
      removeBadChannels_ = true;
    } else if(flag == "false") {
      removeBadChannels_ = false;
      string msg =  GetTBName(); 
      msg += " : bad channels will not be ignored.";
      CsErrLog::Instance()->mes( elWarning, msg);
    }
    else {
      removeBadChannels_ = true;
      CsErrLog::Instance()->mes( elError, "FlagChannels bad option : must be true or false.");
    }
  }        
}


///////////////////////////////////////////////////////////////////////////////////
void CsMicroMegaDetector::BookHistograms() {


  //=== check if histograms are to be booked ===
  CsDetector::ReadHistLevel();

  float tmin = -800;
  float tmax = 800;

  // position of the time peak w/r trigger, before calibration
  float tpeakpos = -9730;

  if (CsInit::Instance()->IsAMonteCarloJob()){
    hLevel_ = None;
    return;
  }
  else {
    if (! CsInit::Instance()->useCalibration()) {
      // real data, no calibrations -> time peak is not at 0
      tmin += tpeakpos;
      tmax += tpeakpos; 	
    }
  }

  string tbn = GetTBName();

  if(hLevel_ >= Normal ) {
    CsHistograms::SetCurrentPath("/MM");

    // hit profile
    hch_ = new CsHist1D(tbn+"_hch",tbn+", hit profile",
			nWir_,0,nWir_);
    hists1D_.push_back(hch_);


    // hit rates 
    hrates_ = new CsHist1D(tbn+"_hrates",tbn+", Hit Rates",
			   nWir_,0,nWir_);
    hists1D_.push_back(hrates_);
      

    // hit times off trigger
    htoff_ = new CsHist1D(tbn+"_htoff",tbn+" , Hit Time - off trigger (ns)",
			  200,tmin,tmax);
    hists1D_.push_back(htoff_);
      

    // Cluster time
    ct_ = new CsHist1D(tbn+"_ct",tbn+" , Cluster Time (ns)",
		       200,-100,100);
    hists1D_.push_back(ct_);

    // Cluster ToT
    ctot_ = new CsHist1D(tbn+"_ctot",tbn+", Cluster ToT",
			 100,0,400);
    hists1D_.push_back(ctot_);

    // Cluster size
    cs_  = new CsHist1D(tbn+"_cs",tbn+", Cluster Size",32,.5,32.5);
    hists1D_.push_back(cs_);

    // Cluster profile

    int cposnbins = 400;
    double cposmin = Wire2Pos( -1 ) + wirD_;
    double cposmax = Wire2Pos( getNWir() ) + wirD_;

    cch_ = new CsHist1D(tbn+"_cch",tbn+", Cluster position (LWRS)",
			cposnbins,cposmin,cposmax);
    hists1D_.push_back(cch_);


    if(hLevel_ >= High) {
      
      // cluster rates 
      crates_ = new CsHist1D(tbn+"_crates",tbn+", Cluster Rates",
			     cposnbins,cposmin,cposmax);
      hists1D_.push_back(crates_);
      
      // cluster times off trigger
      ctoff_ = new CsHist1D(tbn+"_ctoff",tbn+" , Cluster Time - off trigger (ns)",
			    800,-400,400);
      hists1D_.push_back(ctoff_);

      // cluster tot vs time
      ctotvst_ = new CsHist2D(tbn+"_ctotvst",tbn+", Cluster ToT vs Time",
			      50,-100,100, 50, 0, 400);
      hists2D_.push_back(ctotvst_);
      
      // Leading time
      lt_ = new CsHist1D(tbn+"_lt",tbn+", leading time", 200,-12500,-7000);
      hists1D_.push_back(lt_);
      
      // Trailing time
      tt_ = new CsHist1D(tbn+"_tt",tbn+", trailing time", 200,-12500,-7000);
      hists1D_.push_back(tt_);

      // hit time
      ht_ = new CsHist1D(tbn+"_ht",tbn+", hit time", 200,tmin,tmax);
      hists1D_.push_back(ht_);

      // hit time vs chan
      htvsch_ = new CsHist2D(tbn+"_tvsch",tbn+", hit time vs channel",
			     getNWir(),0,getNWir(),100,tmin,tmax);
      hists2D_.push_back(htvsch_);
    }

    CsHistograms::SetCurrentPath("/");
  }


  {
    // ***** W2P PARAMETERISATION
    map<int,double>::iterator ipitch; int npitches;
    for (ipitch=wirPs_.begin(), npitches = 0;
	 ipitch!=wirPs_.end(); ipitch++, npitches++ ) {
      switch (npitches) {
      case 0:
	w2pLow = nWir_-.5;
	w2pP1 = ipitch->second; w2pO1 = 0;
	w2pUp=w2pO2=w2pP2=w2pO3=w2pP3 = 0;  // Defaults
	break;
      case 1:
	w2pUp = nWir_-.5;
	w2pLow = ipitch->first-.5;
	w2pP2 = ipitch->second; w2pO2 = w2pLow*(w2pP1-w2pP2);
	break;
      case 2:
	w2pUp = ipitch->first-.5;
	w2pP3 = ipitch->second;
	w2pO3 = w2pLow*(w2pP1-w2pP3)+(w2pUp-w2pLow)*(w2pP2-w2pP3);
	break;
      }
    }
  }
}


///////////////////////////////////////////////////////////////////////////////////
CsMicroMegaDetector::~CsMicroMegaDetector() {

  for(unsigned i=0; i<hists1D_.size(); i++)
    delete hists1D_[i];
  hists1D_.clear();


  for(unsigned i=0; i<hists2D_.size(); i++)
    delete hists2D_[i];
  hists2D_.clear();

  delete[] digitData_;
}


///////////////////////////////////////////////////////////////////////////////////
bool CsMicroMegaDetector::operator==( const CsMicroMegaDetector& det ) const {
  return( CsDetector::operator==(det) );
}

///////////////////////////////////////////////////////////////////////////////////
bool CsMicroMegaDetector::operator<( const CsMicroMegaDetector& det ) const {
  return( CsDetector::operator<(det) );
}

///////////////////////////////////////////////////////////////////////////////////
void CsMicroMegaDetector::makeMCDecoding() {

  // should I proceed?
  if( !decode_ && !decodeCard_ ) return;

  // Already done?
  if( decodingDone_ ) return;

  // clear
  myDigits_.clear();

  // get a link to relevant pieces...
  CsEvent* event = CsEvent::Instance();

  // this can be done only with zebra binary files...

  double extra = // Extra time width added to compensate for the trigger jitter
    // of the current event being larger than reference trigger jitter.
    CsEvent::Instance()->getExtraTimeWidth();
  double nsmin =  hitTMin_-extra, nsmax = hitTMax_+extra; // Time window
  if( !CsGeant3::Instance()->isAnNtFile() ) {
    double triggerOffset = // Time offset of the trigger w.r.t. the event
      CsEvent::Instance()->getTriggerMCOffset();
    list<CsMCHit*>::iterator Ih;
    for( Ih=myMCHits_.begin(); Ih!=myMCHits_.end(); Ih++ ) { // loop on hits

      CsMCTrkHit* thit = dynamic_cast<CsMCTrkHit*>(*Ih);
      if( thit == 0 ) return;

      // Only charged particles or charged products
      if (((thit->getMCTrack())->getParticle())->getCharge() ||
	  thit->getOrigin()) {

	double t  = thit->getDTime();  // Delay time (ns)
	double ui = thit->getUin();    // Hit in point (DRS)
	double vi = thit->getVin();
	double wi = thit->getWin();
	double uo = thit->getUout();   // Hit out point (DRS)
	double vo = thit->getVout();
	double wo = thit->getWout();
	int id    = thit->GetID();

	t -= // Trigger time (and hence offset) is subtracted from hit times
	  triggerOffset;
	//         ********** TIME -> TDC and TIME CUT **********
	// - This time cut is to be supplemented by a cut performed in the
	//  "clusterize" method. (Which, in the MC case, is done against the
	//  very same time gate as the present one: no equivalent here of the
	//  ambiguities affecting the definition of the master time in RD).
	// - Therefore a loose cut, performed on "t" (as opposed to "tdc"),
	//  could have been done here, saving the CPU it takes to randomize all
	//  hit times, including those far outside time gate. We prefer to do it
	//  already at the MC decoding stage, so that the random generation be
	//  not affected by the actual setting of time gates: the call to
	//  "CsRandom" is done in any case.
	// - This implicitly sets = 0 the relative jitter betweeen CsHit's
	//  originating from a same CsMCDigit. Would the randomisation have been
	//  applied at CsHit instantiation time, the other extreme option, i.e.
	//  relative jitter = (in average) absolute jitter, would have been
	//  taken. Truth lies in between...
	//double tdc = t + tRes_ * CsRandom::gauss();
	double tdc = t + 9 * CsRandom::gauss();
	if (tdc<nsmin || nsmax<tdc) continue;

	// Find center of this hit's subdetector
	double xcm=xcm_, ycm=ycm_, zcm=zcm_;
	map<int, double>::iterator i;
	if((i=xcms_.find(id))!=xcms_.end())
	  xcm = i->second;
	else
	  cout<<"Warning : CsMicroMegaDetector::makeMCDecoding : No xcm for subdetector with id "<<id<<endl;
	if((i=ycms_.find(id))!=ycms_.end())
	  ycm = i->second;
	else
	  cout<<"Warning : CsMicroMegaDetector::makeMCDecoding : No ycm for subdetector with id "<<id<<endl;
	if((i=zcms_.find(id))!=zcms_.end())
	  zcm = i->second;
	else
	  cout<<"Warning : CsMicroMegaDetector::makeMCDecoding : No zcm for subdetector with id "<<id<<endl;

	int err;
	HepMatrix iRotM(3,3); iRotM = rotWRS_.inverse( err );
	double xi = rotDRS_(1,1)*ui+rotDRS_(1,2)*vi+rotDRS_(1,3)*wi+
	  xcm - _xshift;
	double yi = rotDRS_(2,1)*ui+rotDRS_(2,2)*vi+rotDRS_(2,3)*wi+
	  ycm - _yshift;
	double zi = rotDRS_(3,1)*ui+rotDRS_(3,2)*vi+rotDRS_(3,3)*wi+
	  zcm;
	double xo = rotDRS_(1,1)*uo+rotDRS_(1,2)*vo+rotDRS_(1,3)*wo+
	  xcm - _xshift;
	double yo = rotDRS_(2,1)*uo+rotDRS_(2,2)*vo+rotDRS_(2,3)*wo+
	  ycm - _yshift;
	double zo = rotDRS_(3,1)*uo+rotDRS_(3,2)*vo+rotDRS_(3,3)*wo+
	  zcm;
	double Ui = iRotM(1,1)*xi+iRotM(1,2)*yi+iRotM(1,3)*zi; // WRS
	//double Vi = iRotM(2,1)*xi+iRotM(2,2)*yi+iRotM(2,3)*zi;
	double Wi = iRotM(3,1)*xi+iRotM(3,2)*yi+iRotM(3,3)*zi;
	double Uo = iRotM(1,1)*xo+iRotM(1,2)*yo+iRotM(1,3)*zo; // WRS
	//double Vo = iRotM(2,1)*xo+iRotM(2,2)*yo+iRotM(2,3)*zo;
	double Wo = iRotM(3,1)*xo+iRotM(3,2)*yo+iRotM(3,3)*zo;

	int wireF = Pos2Wire(Ui-wirD_);  // First wire for this hit
	int wireL = Pos2Wire(Uo-wirD_);  // Last  wire for this hit
	if( wireL < wireF ) {      // Have first first!
	  int tmp = wireL; wireL = wireF; wireF = tmp;
	}
	if (wireL<0 || wireF>=nWir_) continue;
	wireF = wireF<0      ? 0       : wireF;
	wireL = wireL>=nWir_ ? nWir_-1 : wireL;

	for (int i = wireF; i<=wireL; i++) {
	  // Look if a digit on this detector with this wire already exists. 
	  // Boring, but no idea how to do in other vay... :(
	  list<CsDigit*>::iterator Id; bool found;
	  for (Id=myDigits_.begin(), found = false;
	       Id!=myDigits_.end() && !found; Id++) {
	    if (i==(*Id)->getAddress()) {
	      found = true; 	      // Here it is, add this hit to it...
	      dynamic_cast<CsMCDigit*>(*Id)->addHit(*(*Ih));
	      if (fabs(tdc)<fabs((*Id)->getDatum()))
		// Here we select the hit which is closest to event time. Guess
		// this is what is being done for real (as opposed to MC) data (
		// to be checked). To be checked also that this is the best
		// strategy, particularly in the presence of a trigger offset.
		(*Id)->replaceDatum(tdc);
	    }
	  }
	  if (!found) { // No digits found with these wire and detector
	    CsDigit* digit = new CsMCDigit( *this, i, &tdc );
	    dynamic_cast<CsMCDigit*>(digit)->addHit( *(*Ih) );
	    // add this digit to my list
	    myDigits_.push_back( digit );
	  }
          //LS Eff
          if( mcEffMapsEnabled_ && !MCEff_.empty() ){
            float mceff = MCEff_[i];
	    float random = (float) CsRandom::flat();
#ifdef DEBUG_EFF
            cout<<"**** LS *** : CsMicroMegaDetector::makeMCDecoding : "<< GetTBName() <<" - "<<
              "MCEff_[" <<i<<"]= "<< mceff <<", random "<<random<<endl;
#endif
            if ( random > mceff ) break;
          }
	}
      }
    }
  }
  else {
    if( Qhit.ndig == 0 && Qhit.nhit != 0 ) {
      // No way to make Digits from ntuple files...
      string str = "Decoding not possible on MC Ntuple and no Digits available on this file";
      CsErrLog::Instance()->mes( elFatal, str );
    }

    // Well, digits are in the n-tuple file...
    int CGvers = atoi( (CsGeom::Instance()->getGeomVers()).c_str() + 1 );
    for( int i=0; i<Qhit.ndig; i++ ) { // loop on digits...

      int detn = (Qhit.ip1dig[i]   & 0x0000ffff );
      int ip1  = (Qhit.ip1dig[i-1] & 0xffff0000 ) >> 16;
      int ip2  = (Qhit.ip1dig[i]   & 0xffff0000 ) >> 16;
      if( ip1 == 0 ) ip1 = ip2;
      int dig1 = (Qhit.ip2dig[i]   & 0x0000ffff );
      int dig2 = (Qhit.ip2dig[i]   & 0xffff0000 ) >> 16;

      int mynum;
      if( CGvers < 5 ) {
	mynum = row_;
      }
      else {
	mynum = GetID();
      }
      if( detn != mynum ) continue;

      CsDigit* digit = new CsMCDigit( *this, dig1-1 );
      list<CsMCHit*>::iterator Ih;
      for( int k=ip1-1; k<ip2-1; k ++ ) {
	int hitn =  Qhit.jpdig[k];
	if( hitn != 0 ) {
	  for( Ih=myMCHits_.begin(); Ih!=myMCHits_.end(); Ih++ ) {
	    if( ( (*Ih)->getX() == (Qhit.hit[hitn-1][0]*10.) ) &&
		( (*Ih)->getY() == (Qhit.hit[hitn-1][1]*10.) ) ) {
	      dynamic_cast<CsMCDigit*>(digit)->addHit( *(*Ih) );
	    }
	  }
	}
      }
      // add this digit to its detector list
      myDigits_.push_back( digit );
    }
  }
  decodingDone_ = true;
}


namespace {
  struct CsMicroMegaDigitInfo {  // can't be local.... then just anonymized
    double wtot;
    int    wire;
    CsDigit*    digit;
    CsMicroMegaDigitInfo() : wtot(-1), wire(-1), digit(0) {};
    CsMicroMegaDigitInfo(double wtot, int wire, CsDigit* rfdigit) : wtot(wtot), wire(wire), digit(rfdigit) {};
    bool operator<(const CsMicroMegaDigitInfo& di) const { return this->wtot < di.wtot; }
  };
}

///////////////////////////////////////////////////////////////////////////////////
void CsMicroMegaDetector::clusterize() {

  typedef map<int, list<CsDigit*>::iterator>::iterator IOD;

  clearClusterList();

  int err; HepMatrix iRotM(3,3); iRotM = rotWRS_.inverse( err );
  double wireDCorr = wirD_ + iRotM(1,1) * _deltax + iRotM(1,2) *_deltay;

  double extra = // Extra time width added to compensate for the trigger jitter
    // or the current event being larger than reference trigger jitter.
    CsEvent::Instance()->getExtraTimeWidth();
  double t0    =   // Wrong master trigger was selected by Decoding library
    CsEvent::Instance()->getTriggerTimeCorr();
  double wtmin = t0+hitTMin_-extra, wtmax = t0+hitTMax_+extra;
  static int debug = 1; extra *= debug;
  double ctmin = t0+clTMin_-extra,  ctmax = t0+clTMax_+extra;

  if (associate_) {
    if (!CsEvent::Instance()->isAMonteCarloEvent()) {

      // ********** ``ASSOCIATED'' CLUSTERISATION: I) REAL DATA **********

      if( oDigits_.empty() ) return;
      for(IOD iod=oDigits_.begin();iod!=oDigits_.end();) {

	// priority_queue<CsMicroMegaDigitInfo> digitqueue;
	CsDigit *digit = *(iod->second);

	if (digit->isClusterized())         // ***** SKIP IF ALREADY CLUSTERIZED
	  { iod++; continue; }
	assert (digit->getDataSize() >= 2);
	double wtime = digit->getData()[0];
	if (wtime<wtmin || wtmax<wtime)                 // ***** CUT ON HIT TIME
	  { iod++; continue; }

	digit->setClusterized();
	double wtot = digit->getData()[1];
	int lastwire = digit->getAddress();

	int delta = 0;

	// digitqueue.push(CsMicroMegaDigitInfo(wtot, lastwire, digit));

	list<CsDigit*> lstdigit;
	lstdigit.push_back(digit);
	
	vector<MM::Digit> mmdigits;
	mmdigits.push_back(MM::Digit(lastwire,wtime,wtot,Pitch(lastwire)>.4));
	

	do {
	  iod++;
	  if( iod != oDigits_.end() ) {
	    CsDigit *dtmp = *(iod->second);
	    delta = dtmp->getAddress() - lastwire;
	    if( delta <= 2 ) {
	      assert (dtmp->getDataSize() >= 2);
	      wtime = dtmp->getData()[0];
	      if (wtime<wtmin || wtmax<wtime) continue; // ***** CUT ON HIT TIME
	      wtot = dtmp->getData()[1];
	      lastwire = dtmp->getAddress();
	      // digitqueue.push(CsMicroMegaDigitInfo(wtot, lastwire, dtmp));
	      lstdigit.push_back(dtmp);
	      mmdigits.push_back
		(MM::Digit(lastwire,wtime,wtot,Pitch(lastwire)>.4));
	    }
	  }
	} while( iod != oDigits_.end() && delta <= 2 );

	MM::Cluster *mmcluster = new MM::Cluster(mmdigits);

	vector<MM::Cluster*> mmclusters ;
	
	if(splitClustersMin_ < mmcluster->wire && 
	   mmcluster->wire  <  splitClustersMax_)  
	  mmclusters = mmcluster->Split();
	else 
	  mmclusters.push_back(mmcluster);	


	for(unsigned immc = 0; immc < mmclusters.size(); immc++) {
	  
	  double wire = mmclusters[immc]->wire;
	  double cltime = mmclusters[immc]->t;
	  double clToT = mmclusters[immc]->tot;
	  int nwires = mmclusters[immc]->s;
	  
	  //#define ALTERNATIVE_W2P
#ifdef ALTERNATIVE_W2P
	  // `cause I suspect "Wire2Pos" to be wrong
	  double u = W2P(wire);
#else
	  double u = Wire2Pos(wire);
#endif
	  
	  u += wireDCorr;
	  double v = 0;
	  double w = zcm_;
	  
	  if (!(ctmin<cltime && cltime<ctmax)) continue;
	  if (cutTriangle_ && ! cutTriangle_->IsInside(cltime,clToT) ) continue;
	
#ifdef uM_DETECTOR_STUDIES
	  if (clToT<ToTMin_) continue;
	  if (nwires>=clsMax_) continue;
#endif

	  if( hLevel_ >= Normal ) {
	    cch_ -> Fill(u);  // cluster pos in *WRS*
	    ct_->Fill(cltime); 
	    ctot_->Fill(clToT); 
	    cs_->Fill((double)nwires);
	  }
	  if( hLevel_ >= High ) {
	    if(cltime > cltoffMin_ && cltime < cltoffMax_ ) {
	      ctoff_ -> Fill(cltime);
	      crates_ -> Fill(u);
	    }
	    ctotvst_ -> Fill(cltime,clToT);
	  }
	  
	  
	  // Set errors:
	  HepMatrix cov(3,3,0);  // Zero matrix
	  
	  if(cresolution_ < 0)
	    cov(1,1) = pow( Pitch(wire) / sqrt(12.0), 2 ); // wire pitch / sqrt(12)
	  else {
	    cov(1,1) = pow( cresolution_ , 2 );
	  }
	  
	  if (fabs(ang_-90)<1) {           // wire length/2
	    cov(2,2) = pow( (getXsiz())/2., 2 );
	  }
	  else {
	    cov(2,2) = pow( (getYsiz())/2./cos(ang_/180.*(M_PI)), 2 );
	  }
	  cov(3,3) = 1.;               // assume 1 mm resolution in Z

	  // Save the cluster:

	  const vector<MM::Digit>& mmdigvec = mmclusters[immc]->GetDigits();
	  assert( !mmdigvec.empty() );
	  float minchan = mmdigvec[0].ch;
	  float maxchan = mmdigvec.back().ch;
	  
	  CsCluster* cluster = new CsCluster( u, v, w, cov );
	  for (list<CsDigit*>::iterator it = lstdigit.begin();
	       it!=lstdigit.end(); it++) {
	    int ch = (*it)->getAddress();
	    if(ch<minchan) continue;
	    if(ch>maxchan) break;
	    cluster->addDigit(**it);
	  }
	  
	  cluster->addDet( *this );
	  cluster->setTime(cltime);
	  cluster->setAnalog(clToT);
	  addCluster( *cluster );
	  
	  delete mmclusters[immc];
	}
      }

#ifdef uM_SINGLE_TRACK
      if (getMyClusters().size()>=(unsigned int)nclMax_) {
	for (unsigned int i = 1; i<getMyClusters().size(); i++) {
	  CsEvent::Instance()->subtractCluster();
	}
	clearClusterList();
      }
#endif
      oDigits_.clear();
    }
    else {

      // ********** ``ASSOCIATED'' CLUSTERISATION: II) MC DATA **********

      list<CsDigit*>::iterator Id;
      list<CsDigit*> digits = getMyDigits();

      vector<list<CsDigit*>::iterator> iterators;
      iterators.clear();

      // protection
      if( digits.empty() ) return;

      for( Id=digits.begin(); Id!=digits.end(); ) {

	// skip if already clusterized
	if( (*Id)->isClusterized() ) continue;
	(*Id)->setClusterized();

	iterators.clear();
	int firstwire    = int((*Id)->getAddress());
	int lastwire     = firstwire;
	int previouswire = firstwire;
	int delta = 0;
	double cltime = (*Id)->getDatum();
	iterators.push_back( Id );

	do {
	  Id++;
	  if( Id != digits.end() ) {
	    delta = abs(int((*Id)->getAddress()) - previouswire);
	    if( delta == 1 ) {
	      lastwire = int((*Id)->getAddress());
	      previouswire = lastwire;
	      // Time: preliminary. The following does not reproduce
	      // faithfully what's happening with RD.
	      double htime = (*Id)->getDatum();
	      if (fabs(htime)<fabs(cltime)) cltime = htime;
	      iterators.push_back( Id );
	    }
	  }
	} while( Id != digits.end() && delta == 1 );

	if (!(clTMin_<cltime && cltime<clTMax_)) continue;

	double wire = double( lastwire + firstwire ) / 2.;
	int nwires = abs( lastwire - firstwire ) + 1;

	// Set the values in the WRS
	// u: perpendicular to the detector wires direction
	// v: parallel to the detector wires direction
	// w: = Z coordinate in MRS
	// origin: the origin of the MRS

	double u = wireDCorr + Wire2Pos(wire);
	double v = 0;
	double w = zcm_;

	// Set errors:
	HepMatrix cov(3,3,0);  // Zero matrix

	if(cresolution_ < 0)
	  cov(1,1) = pow( Pitch(wire) / sqrt(12.0), 2 ); // wire pitch / sqrt(12)
	else {
	  cov(1,1) = pow( cresolution_ , 2 );
	}

	if (fabs(ang_-90)<1) {           // wire length/2
	  cov(2,2) = pow( (getXsiz())/2., 2 );
	}
	else {
	  cov(2,2) = pow( (getYsiz())/2./cos(ang_/180.*(M_PI)), 2 );
	}
	cov(3,3) = 1.;               // assume 1 mm resolution in Z

	// Save the cluster:
	CsCluster* cluster = new CsCluster( u, v, w, cov );
	for( unsigned int i=0; i<iterators.size(); i++ ) {
	  cluster->addDigit( *(*(iterators[i])) );
	}
	cluster->setTime(cltime);
	cluster->addDet( *this );
	addCluster( *cluster );
      }
    }
  }
  else {

    // ********** NON ``ASSOCIATED'' CLUSTERISATION **********

    list<CsDigit*>::iterator Id;
    list<CsDigit*> digits = getMyDigits();

    vector<list<CsDigit*>::iterator> iterators;
    iterators.clear();

    // protection
    if( digits.empty() ) return;

    for( Id=digits.begin(); Id!=digits.end(); Id++ ) {

      CsDigit *digit = *Id;

      // skip if already clusterized
      if( digit->isClusterized() ) continue;
      digit->setClusterized();

      int wire = digit->getAddress();

      // Set the values in the WRS
      // u: perpendicular to the detector wires direction
      // v: parallel to the detector wires direction
      // w: = Z coordinate in MRS
      // origin: the origin of the MRS

#ifdef ALTERNATIVE_W2P
	// `cause I suspect "Wire2Pos" to be wrong
	double u = W2P(wire);
#else
	double u = wireDCorr + Wire2Pos(wire);
#endif
      double v = 0;
      double w = zcm_;

      // Set errors:
      HepMatrix cov(3,3,0);  // Zero matrix

      if(cresolution_ < 0)
	cov(1,1) = pow( Pitch(wire) / sqrt(12.0), 2 ); // wire pitch / sqrt(12)
      else {
	cov(1,1) = pow( cresolution_ , 2 );
      }

      if (fabs(ang_-90)<1) {           // wire length/2
	cov(2,2) = pow( (getXsiz())/2., 2 );
      }
      else {
	cov(2,2) = pow( (getYsiz())/2./cos(ang_/180.*(M_PI)), 2 );
      }
      cov(3,3) = 1.;               // assume 1 mm resolution in Z

      // Save the cluster:
      CsCluster* cluster = new CsCluster( u, v, w, cov );
      cluster->addDigit( *digit );
      cluster->setTime( digit->getDatum() );
      cluster->addDet( *this );
      addCluster( *cluster );
    }
  } // !!!!!!!!!!!!!
  sortClusters();
  setClusteringDone();
}


///////////////////////////////////////////////////////////////////////////////////
void CsMicroMegaDetector::DecodeChipDigit(const CS::Chip::Digit &digit)
{

  // ------------------------- DecodeChipDigit -------------------------

  //   I) Instantiates CsDigits from chip digits
  //     - Associating ONE trailing edge chip digit to ONE leading edge such.
  //     - Applying calibration and trigger time.
  //  II) Fills "(CsDet)this->myDigits_" list w/ them.
  // III) Fills also the CsMM's "oDigits_" map of references to the CsDigits
  //     indexed by their channel #. This is used to single out cases where
  //     2 hits for a same channel, and retain the most appropriate (which
  //     is most appropriate is decided depending upon "uM_KEEP_CLOSEST_HIT").

  const CS::ChipF1::Digit *d = dynamic_cast<const CS::ChipF1::Digit *>(&digit);
  if( d==NULL )
    throw CS::Exception("CsMicroMegaDetector::DecodeRawData(): Wrong digit type");

#define CsMM_MULTIDECODE_WO_CLUSTERING
#ifdef CsMM_MULTIDECODE_WO_CLUSTERING
  //  The "oDigits" map used to be only cleared by "CsMM::clusterize" method.
  // Which is OK as long as coral (more precisely it "CsEvent::getNextEvent"
  // engine) executes systematically "clusterize" once per "DecodeChipDigit". It
  // turns out not to be the case (for CsEvent.cc,v1.243), cf. event #82837507
  // of raw data file "03P1H/cdr14008-31520" (one reaches it directly by
  // "events to skip 29294") mentioned by Albert. Which gets coral crash upon
  // dereferencing a oDigits_ reference to a longer existing CsDigits.
  //  This multiple "DecodeChipDigit" w/o "clusterize" must originate in
  // a bad handling of decoding troubleshootings, for the above mentioned crash
  // occurs after a series of the following type of error messages:
  // ERROR, on Tue, 11/Jul/2006 14:42:24.340052 (GMT)
  //  from: CsEvent.cc   912
  // `+----------------------------------------------------------+ 
  //     Event skipped due to decoding troubles. 
  //     Run: 31520, Event: 81826343 
  //  +----------------------------------------------------------+'
  // A long series by ant standard: the troubleshootings start w/ #81809337
  // and last until the abovementioned #82837507.
  // NOTA BENE: I (Y.B.) haven't investigated the impact of the present patch.
  //           Therefore it may be wiser to disable it once the bug in CsEvent
  //           at handling decoding troubleshooting is fixed.
  if (!myDigits_.size()) oDigits_.clear();
#endif

  int chan = d->GetChannel();

  //#define CsMM_DEBUG_MULTIDECODE
#ifdef CsMM_DEBUG_MULTIDECODE
  //  This piece of code conditioned by "CsMM_DEBUG_MULTIDECODE" traces the
  // occurence of the crash mentioned supra: it's in MM01U1 channel #555.
  static int debug_multidecode;
  if (GetTBName()=="MM01U1__" && debug_multidecode==2) {
    map<int,list<CsDigit*>::iterator>::iterator d555 = 
      oDigits_.find(555);
    printf("\nEntering \"DecodeChipDigit\" for \"MM01U1__\": %d myDigits_, %d oDigits_ including ref to CsDigit @ channel 555\n",
	   myDigits_.size(),oDigits_.size());
  }
#endif

  if (useCalib_) {                                // ***** REMOVE NOISY CHANNELS
    if (CsInit::Instance()->useCDB())
      if(removeBadChannels_ && calib_data_ext[chan].flag == 1) return;
  }

  if (chan<0 || chan>=getNWir())                   // ***** CHANNEL OUT OF RANGE
    throw CS::Exception
      ("CsMicroMegaDetector::DecodeRawData(): Unexpected wire number %d for %s",
       d->GetChannel(),GetTBName().c_str());

  double time = d->GetTimeDecoded()/d->GetTimeUnit();// ***** TRIGGER TIME CORR.

  if (IsLeading(d->GetAmplitude())) {                   // ***** LEADING EDGE...
    digitData_[0] = time; lastChan_ = chan;    // Store
  }
  else {                                            // ***** ...TRAILING EDGE...
    if (chan==lastChan_) {                 // ***** ... => INSTANTIATE A CsDigit
      lastChan_ = -1;                          // Close search for trailing
      digitData_[1] = time - digitData_[0];    // ToT
      if( hLevel_ >= High ) { lt_->Fill(digitData_[0]); tt_->Fill(time); }
      digitData_[0]=digitData_[0]*leadtWght_ + time*(1-leadtWght_); // Hit time

      if (useCalib_) {                                // ***** APPLY CALIBRATION
	if (CsInit::Instance()->useCDB()) {
	  // CASE 1: Any of the "./src/condb" DB
	  digitData_[0] -= calib_data_ext[chan].t0;
	}
      }

      if (hLevel_>=High) {                // ***** CsDigit HIGH LEVEL HISTOGRAMS
	htvsch_-> Fill(chan,digitData_[0]); ht_-> Fill(digitData_[0]);
      }


      // ATTENTION FROM HERE ON NANOSECONDS !!!!
      //                                           ***** CONVERSION F1 UNIT-> ns
      // (Simpler not to use "CS::ChipF1::GetTimeUnit()"...anyway it does
      // not matter so much)
      digitData_[0] = digitData_[0] * f1Tick_;
      digitData_[1] = digitData_[1] * f1Tick_;

      double extra =// ***** EVENT MAY HAVE LARGE TRIGGER JITTER? ENLARGE t GATE
	// (The expectation for the trigger jitter may be refined later on,
	// after all data including the trigger pattern TDC have been decoded:
	// the present "getExtraTimeWidth" is an upper bound and is used to
	// cut away hits that can't in any case be retained, in order to speed
	// up the processing.)
	CsEvent::Instance()->getExtraTimeWidth();
      if (digitData_[0]<hitTMin_-extra ||
	  hitTMax_+extra<digitData_[0]) return;         // ***** CUT ON HIT TIME

      if (hLevel_>= Normal) {
	hch_-> Fill(chan);
	if (hitTOffMin_<digitData_[0] && digitData_[0]<hitTOffMax_) {
	  hrates_ -> Fill(chan); htoff_ -> Fill(digitData_[0]);
	} 
      }

      myDigits_.push_back(new CsDigit(*this,chan,digitData_,2));

      //                                           ***** SEVERAL DIGITS PER WIRE
      list<CsDigit*>::iterator digiti = --(myDigits_.end());

#ifdef CsMM_DEBUG_MULTIDECODE    // Tracing bug in MM01U1 channel #555...
      if (GetTBName()=="MM01U1__") {
	printf("TBname %s:  New CsDigit @ channel #%d w/ data = %.4f\n",
	       GetTBName().c_str(),chan,myDigits_.back()->getData()[0]);
	if (oDigits_.size()==0) debug_multidecode = 0;
	if (chan==555) debug_multidecode = 1;
	else if (debug_multidecode) {
	  debug_multidecode = 2;
	  map<int,list<CsDigit*>::iterator>::iterator d555 = 
	    oDigits_.find(555);
	  printf("while hit @ #555 already encountered, w/ data = %.4f\n",
		 (*(d555->second))->getData()[0]);
	}
      }
#endif

      //      ********** RETAIN A SINGLE HIT per CHANNEL **********
      map<int,list<CsDigit*>::iterator>::iterator previous=oDigits_.find(chan);
      if (previous!=oDigits_.end()) {	
#define uM_KEEP_CLOSEST_HIT
#ifdef uM_KEEP_CLOSEST_HIT               // ***** KEEP THE CsDigit CLOSEST to T0
  	if (fabs((*(previous->second))->getData()[0])>fabs(digitData_[0])) {
  	  previous->second = digiti;
  	}	
#else                                // ***** RETAIN THE CsDigit WITH HIGHER ToT
  	if (((*(previous->second))->getData())[1]<digitData_[1]) {
  	  previous->second = digiti;
  	}	
#endif
      }
      else oDigits_[chan]= digiti;        // ***** FILL MAP channel# <-> CsDigit
    }
  }
}


/////////////////////////////////////////////////////////////////////////////////////////

double CsMicroMegaDetector::W2P(double wire) {
  double x;
  if      (wire<w2pLow) return (wire*w2pP1+w2pO1);
  else if (wire<w2pUp ) return (wire*w2pP2+w2pO2);
  else                  return (wire*w2pP3+w2pO3);
}


//----------------------------------------------------------------------------

void CsMicroMegaDetector::readCalibration(time_t timePoint){

  // prior to calibration reading, some build in checks are to be done
  // for partialy equipped MM detectors...
  // In previous versions this was done in BookHistograms...

  string tbn = GetTBName();
  channel0_ = 0; nWtot_ = nWir_;

  // !!!!!!!!!!!!! temporary: This may not be the right place
  // for the following blocks but is the best for the time being: one would
  // need a method that is systematically called upon starting
  // processing a new run.

  //=== Read the calibration DB
  CDB::Time tp(timePoint,0);

  string strdata("");
  cdb_->read(GetTBName(),strdata,tp);
  istringstream istrdata(strdata);
  istrdata >> calib_data_ext;

  if (getenv("COND_DB_DEBUG")!=0 || getenv("COND_DB_DEBUG_MM")!=0) {
//     cout <<"string read: "<<strdata<<endl;
    vector<MM::ChannelCalibExt>::iterator it;
    for (it = calib_data_ext.begin(); it != calib_data_ext.end(); it++) {
      cout <<GetTBName()<<" data_ext: ch "<<(*it).ch<<" t0 "<<(*it).t0<<" flag "<<(*it).flag<<endl; 
    }
  }

  if (calib_data_ext.size() == 0) {
    tm *t = localtime(&tp.first);
    cout << GetTBName() << ", no calibration for local time "
	 << t <<" in CDB"<< endl;
    useCalib_ = false;
  } else {
    if (calib_data_ext.size() != (unsigned) nWtot_) {
      CsErrLog::msg(elError,__FILE__,__LINE__,
		    "Size of Calibration not correct! Should be %d, Is %d",
		    nWtot_,calib_data_ext.size());
      useCalib_ = false;
    }
  }
}


//LS Eff                                                                                  
void CsMicroMegaDetector::readMCEffMaps(time_t timePoint)
{    
  if( mccdb_ == NULL ) return;
  CDB::Time tp(timePoint,0);
  tm *t = localtime(&tp.first);
  string data("");
  mccdb_->read(GetTBName(),data,tp);
  istringstream is(data);
  string line;
  while( getline(is,line) ){
    if( line.length()==0 || line[0]=='#' ) 
      continue;
    int c;
    float mcEf, mcEf_err;
    if( 3!=sscanf(line.c_str(),"%d %g %g",&c,&mcEf,&mcEf_err) )
      throw std::logic_error("CsMicroMegaDetector_MCEff::read(): bad line structure.");
    MCEff_.push_back(mcEf);
    MCEff_err_.push_back(mcEf_err); }
}
