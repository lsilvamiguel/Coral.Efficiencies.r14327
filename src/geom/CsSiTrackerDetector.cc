// $Id: CsSiTrackerDetector.cc 14069 2015-09-17 20:44:46Z lsilva $

/*!
   \file    CsSiTrackerDetector.cc
   \brief   Compass SiTracker like detector Class.
   \author  Benigno Gobbo
   \version $Revision: 14069 $
   \date    $Date: 2015-09-17 22:44:46 +0200 (Thu, 17 Sep 2015) $
*/

#include "CsSiTrackerDetector.h"

#include <math.h>
#include <string.h>
#include <strings.h>
#include <functional>
#include <TMath.h>
#include "CsZone.h"
#include "CsDigit.h"
#include "CsOpt.h"
#include "CsInit.h"
#include "CsEvent.h"
#include "CsDigit.h"
#include "CsMCDigit.h"
#include "CsComgNtCommons.h"
#include "CsGeom.h"
#include <cstdlib>
#include "CsRandom.h"
#include "CsMCTrkHit.h"
#include "CsGeant3.h"
#include "CsMCTrack.h"


#include "DaqDataDecoding/ChipAPV.h"
#include "CDB.h"
#include "silicon_timing.c"
#include "CsSiTrackerResCor.h"

using namespace std;
using namespace CLHEP;

extern QhitType Qhit;
// ------------------------------------------------------------
CsSiTrackerDetector::CsSiTrackerDetector( const int    row,
		        const int    id,    const char* name,    const char *TBname,
			const int    unit,  const int    type,
			const double rdLen, const double xsiz,
			const double ysiz,  const double zsiz,
			const double xcm,   const double ycm,
			const double zcm,   const HepMatrix rotDRS,
			const HepMatrix rotWRS,
			const double wirD,  const double ang,
			const int    nWir,  const double wirP,
			const double eff,   const double bkg,
			const double tGate ) :

  // initialisation of member variables

  CsDetector( row, id, name, TBname, unit, type, rdLen, xsiz, ysiz, zsiz,
              xcm, ycm, zcm, rotDRS, rotWRS, wirD, ang,
              nWir, wirP, eff, bkg, tGate ),

  //added for MC simulation
  spSig_(0.),
  eGain_(0.),
  eGSig_(0.),
  sWidth_(0.),
  tRes_(0.),
  //
  cind_clus_setup(false),
  has_calib_data(false),
  isX(TString(TBname).Contains("X")),
  isU(TString(TBname).Contains("U")),
  isY(TString(TBname).Contains("Y")),
  isV(TString(TBname).Contains("V")),
  rescor(0.)


{
// -------------------------------------------------------------

  mateDet=NULL;

  // calculate number of chips (8 or 10)
  ChipsPerPlane=getNWir()/ChipChannels+!!(getNWir()%ChipChannels);

  if (getenv("COND_DB_DEBUG")!=0 || getenv("COND_DB_DEBUG_SI")!=0) {
    cout << "CsSiTrackerDetector[" << GetTBName() << "]::CsSiTrackerDetector()" << endl;
    cout << "ChipsPerPlane: " << ChipsPerPlane << endl;



  }

  const string SI = "SI";  // tag is "SI" or TB name.

  // initiate Common Mode Noise counter
  for(int i = 0;i<10;i++){
    fEventNr[i]=0;

  }


  // Options
  // =======

  // 1. decoding card
  // ----------------
  decodeCard_ = false;
  string tag = "";
  string key = "make decoding";
  bool status = CsOpt::Instance()->getOpt( tag, key );
  if( status ) {
    list<string> options;
    list<string>::iterator Is;
    status = CsOpt::Instance()->getOpt( tag, key, options );
    if( status ) {
      for( Is=options.begin(); Is!=options.end(); Is++ ) {
	if( *Is == "SI" || *Is == "SiTracker" || *Is == "all" ) {
	  decodeCard_ = true;
	}
      }
    }
    else {
      decodeCard_ = true;
    }
  }


  // 2. timing cuts
  // --------------

  // Two methodes can be applied:
  // - ratio cut (default)
  // - time cut, needs two extra variables, namely the limits
  //
  // modes can be switched on/off independently:
  // SI      ClusterTime [0] -1:Off, 0,1: ratio, 2: time, 3 ratio && time
  //
  // cut can be applied to projections sepeartely:
  // SI       ClusterTime [1-2] -100 100
  // SI01V1__ ClusterTime [1-2] -15 15


  // default: ratio cut only
  fRatioCut  = 1; fTimeCut   = 0; fLCut = -1e6; fRCut= 1e6;

  key = "ClusterTime";
  vector<double> v;
  if (CsOpt::Instance()->getOpt(TBname,key,v) ||
      CsOpt::Instance()->getOpt(SI,key,v)) {

    if (v.size()>0){
      switch ((int)v[0]) {
      case -1:
	fTimeCut=0;
	fRatioCut=0;
	break;
      case 0:
	// "0" has to be forced to be default,
	// in case it is not specially initialized
	// * default values *
	break;
      case 1:
      case 2:
      case 3:
        fTimeCut  = ((int)v[0])&2;
	fRatioCut = ((int)v[0])&1;
	break;
      default:
	CsErrLog::mes(elFatal,"Wrong time cut in option file: SI ClusterTime [0]: {-1,0,1,2,3}!") ;
	break;
      }

      if (v.size()==3){
	fLCut = v[1];
	fRCut = v[2];
      	CsErrLog::
	  msg(elInfo,__FILE__,__LINE__,"\"%s\" %.1f < ClusterTime < %.1f",
	      TBname,fLCut,fRCut);
      }
      if (v.size()==2)
	CsErrLog::mes(elFatal,"Wrong number of variables for time cut in option file SI ClusterTime!") ;

    }
  }
  CsErrLog::
    msg(elInfo,__FILE__,__LINE__,
	"\"%s\" Using %d RatioCut and %d TimeCut",
	TBname,fRatioCut,fTimeCut);


  // 3. multiplicity cut
  // -------------------

  // special cut for checking amplitude spectra and
  // timing behaviour in cleaner events
  key = "MultiplicityCut";
  fMultCut=1000;

  double m;

  if (CsOpt::Instance()->getOpt(TBname,key,m) ||
      CsOpt::Instance()->getOpt(SI,key,m)) {
    if (m<1)
      CsErrLog::mes(elFatal,
		    "Error syntax in specifying Multiplicity cut");
    fMultCut = m;
    CsErrLog::msg(elInfo,__FILE__,__LINE__,"\"%s\" Multiplicity cut: %.1f",
		  TBname,fMultCut);
  }


  // 4. clustering selection
  // -------------------

  // possibility to select between cinderella_timing and clustering with correct
  // error propagation of silicon time, or use of old ratio clustering

  key = "ClusteringSelect";

  cluster_cind_opt=true; //switch in option file
  cluster_mate_opt=true;

  list<string> options_cl;
  list<string>::iterator Idl;

  if (CsOpt::Instance()->getOpt(TBname,key,options_cl) ||
      CsOpt::Instance()->getOpt(SI,key,options_cl)){

    for( Idl=options_cl.begin(); Idl!=options_cl.end(); Idl++ ) {
      if( *Idl == "ratio"  ) {
	cluster_cind_opt=false;
	cluster_mate_opt=false;
	CsErrLog::msg(elInfo,__FILE__,__LINE__,"\"%s\" Ratio clustering is used ",
		  TBname);
      }
      else if( *Idl == "cind"){
	cluster_cind_opt=true;
	cluster_mate_opt=false;
	CsErrLog::msg(elInfo,__FILE__,__LINE__,"\"%s\" Cinderella clustering is used ",
		  TBname);
      }
      else if( *Idl == "mate"){
	cluster_cind_opt=true;
	cluster_mate_opt=true;
	CsErrLog::msg(elInfo,__FILE__,__LINE__,"\"%s\" Mate clustering is used, if not possible fall back to Cinderella ",
		  TBname);
      }


      else
	CsErrLog::mes(elFatal,
		      "Error syntax in specifying silicon clustering selection in option file\n, choose 'SI ClusteringSelect ratio' or 'SI ClusteringSelect cind'!! ");

    }
  }


  // 5. E18 timing calibration files for testing
  // possibility to select between official time calibrations checked in in the DB or
  // calibration files located in Munich /nfs/hicran/project/compass/silicon/time_calib/timingfiles/
  // this is necessary to check the new calibrations before they are comitted to the database
  key = "TimeCalibMunich";
  time_calib_db=true; //should be the default
  time_calib_munich =false; // testing of calibfiles has to be swiched on explicitly;
  list<string> options_calib;
  list<string>::iterator Id_calib;

  if (CsOpt::Instance()->getOpt(TBname,key,options_calib) ||
      CsOpt::Instance()->getOpt(SI,key,options_calib)){

    for( Id_calib=options_calib.begin(); Id_calib!=options_calib.end(); Id_calib++ ) {
      if( *Id_calib == "on"  ) {
	time_calib_db=false;
	time_calib_munich =true;
	CsErrLog::msg(elError,__FILE__,__LINE__,"\"%s\" TIME CALIBRATIONS READ FROM FILE !! ",
		      TBname);
      }
      else if( *Id_calib == "off"){
	time_calib_db=true;
	time_calib_munich =false;
	CsErrLog::msg(elInfo,__FILE__,__LINE__,"\"%s\" Cinderella clustering is used ",
		      TBname);
      }
      else{
	calib_dir = *Id_calib;
	//string calib_dir2 = string(*Id_calib);	
	time_calib_munich =true;
	time_calib_db=false;
	string fileinin = "/nfs/hicran/project/compass/silicon/time_calib/timingfiles/"+calib_dir+"/"+GetTBName();
	CsErrLog::msg(elError,__FILE__,__LINE__,"\"%s\" CALIBRATIONS from FILE are used %s",
		      TBname,fileinin.c_str());
      }
    }
  }


  // 6. setup Sergei's alignment patch
  // -------------------
  // Alignment patch for hadron runs 42320 - 43350,
  // values in CsSiTrackerResCor.h can be switchen on and off
  // per default on

  key = "AlignmentPatch";

  align_patch_run = 0;
  align_patch = true;

  list<string> options_align;
  list<string>::iterator Id_align;

  if (CsOpt::Instance()->getOpt(TBname,key,options_align) ||
      CsOpt::Instance()->getOpt(SI,key,options_align)){

    for( Id_align=options_align.begin(); Id_align!=options_align.end(); Id_align++ ) {
      if( *Id_align == "off"  ) {
	align_patch=false;
	CsErrLog::msg(elInfo, __FILE__, __LINE__, "\"%s\" Silicon alignment patch "
                      "for runs 42320 - 43350 is NOT used!", TBname);
      }
    }
  }



  // 7. Correction of timing by a half clockcycles
  // --------------------------------------------------------
  // switch on/off the shift (used e.g. for 2009 data)
  // of the timing by a half clock cycle
  // default setting is: off

  key = "Time12nsJumps";

  time_12nsjumps = false;

  list<string> options_tjumps;
  list<string>::iterator Id_tjumps;

  if (CsOpt::Instance()->getOpt(TBname,key,options_tjumps) ||
      CsOpt::Instance()->getOpt(SI,key,options_tjumps)){

    for( Id_tjumps=options_tjumps.begin(); Id_tjumps!=options_tjumps.end(); Id_tjumps++ ) {

      if ( *Id_tjumps == "on" ){
	time_12nsjumps = true;
	CsErrLog::msg(elInfo,__FILE__,__LINE__,"\"%s\" 12ns TimeJump correction active!",
		      TBname);
      }
      else if ( *Id_tjumps == "off" )      time_12nsjumps = false;

      else if (time_calib_munich){                  // if TimeCalibMunich was set before Time12nsJumps in option file AND another argument than "on" or "off" is given, jumps will be read from same folder as calibrations
	string filein_tjumps = "/nfs/hicran/project/compass/silicon/time_calib/timingfiles/"
	                    +calib_dir+"/"+GetTBName()+"12nsjumps";
	CsErrLog::msg(elError,__FILE__,__LINE__,"\"%s\" TimeCalibMunich switched on and TimeCorrection active:"
		      " reading correction from %s",
		      TBname, filein_tjumps.c_str());
      }

      else   CsErrLog::msg(elError,__FILE__,__LINE__,"\"%s\" unrecognised option for timing jumps\nNo jumps will be read!",
			    TBname);

    }
  }


  // 8. 2-strip cluster position correction
  // This provides the option to switch off the position correction of 2-strip clusters
  // that is obtained from the residual vs. eta histogram fit.

  // Default is on, switch off to create new eta correction (via "SI EtaCorrection off" in the option file)

  key = "EtaCorrection";
  eta_correction = true;
  eta_calib = true;

  list<string> options_etacorr;
  list<string>::iterator Id_etacorr;

  if (CsOpt::Instance()->getOpt(TBname,key,options_etacorr) ||
      CsOpt::Instance()->getOpt(SI,key,options_etacorr)){

    for( Id_etacorr=options_etacorr.begin(); Id_etacorr!=options_etacorr.end(); Id_etacorr++ ) {
      if ( *Id_etacorr == "on" ){
	eta_correction = true;
	CsErrLog::msg(elInfo,__FILE__,__LINE__,"\"%s\" Charge sharing correction on",
		      TBname);
      }
      else if ( *Id_etacorr == "off" ){
	eta_correction = false;
	CsErrLog::msg(elError,__FILE__,__LINE__,"\"%s\" Charge sharing correction OFF! Use ONLY for creating new corrections!",
		      TBname);
      }
      else if ( *Id_etacorr == "hard" ){
	eta_correction = true;
	eta_calib = false;
	CsErrLog::msg(elError,__FILE__,__LINE__,"\"%s\" Hard-coded eta corrections are used! (Instead of calib)\nAre you sure you know what you're doing?",
		      TBname);
      }
      else{
	eta_correction = true;
	CsErrLog::msg(elError,__FILE__,__LINE__,"\"%s\" Charge sharing correction: unrecognized option! (\"on\" assumed)",
		      TBname);
      }

    }
  }
  //9. Stuff for MC

  //key = "amplitudeMCdecoding";

  amplitudeMCDecoding_ = false;
  if (CsOpt::Instance()->getOpt(TBname,"amplitudeMCDecoding") ||
      CsOpt::Instance()->getOpt(SI,"amplitudeMCDecoding")) {
    amplitudeMCDecoding_ = true;
    CsErrLog::msg(elWarning,__FILE__,__LINE__,
  		  "%s: amplitudes simulation in MC decoding will be used",TBname);
  }

  // ***** MONTE CARLO: Amlitudes simulation parameters *****
  //
  // spSig_ detector space resolution (mm)
  // eGain_ effective gain - should be tuned to reproduce cluster amlitudes
  // eGSig_ gain sigma (a.u.) for amplitude correlation, for example.
  // sWidth_ signal width (mm) (should be tuned to have correct number of strips/cluster)
  // tRes_ time resolution (ns)
  // These pars can be made dependent upon detector's TB name, cf. supra:
  // GM ampParsMC [0-4] 0.05  2000. 20.  0.30  12.
  // GM01X1__ ampParsMC [0-4] 0.05  2000. 20. 0.30  12.

  key = "ampParsMC";
  if (CsOpt::Instance()->getOpt(TBname,key,v) ||
      CsOpt::Instance()->getOpt(SI,key,v)) {
    if (v.size()!=5) // 5??
      CsErrLog::msg(elFatal,__FILE__,__LINE__,
		    "%s: Syntax error in specifying MC amplitude parameters!",TBname);
    spSig_ = v[0]; eGain_ = v[1]; eGSig_ = v[2]; sWidth_ = v[3]; tRes_ = v[4];
    CsErrLog::msg(elInfo,__FILE__,__LINE__,"%s: space res. %f mm, gain %f,sigma gain %f , signal width %f mm, time res. %f ns",
                  TBname,spSig_,eGain_,eGSig_,sWidth_,tRes_);
  }

  if(amplitudeMCDecoding_) {
    if( sWidth_ <= 0. || eGain_ <=0. ){
     CsErrLog::msg(elFatal,__FILE__,__LINE__,
                    "%s: amplitude simulation impossible, wrong "
                    "width=%f or eGain=%f  ",
                     TBname,sWidth_,eGain_);
    }
  }

  key ="cls3ParsMC";
  if (CsOpt::Instance()->getOpt(TBname,key,v) ||
      CsOpt::Instance()->getOpt(SI,key,v)) {
    if (v.size()!=6)
      CsErrLog::msg(elFatal,__FILE__,__LINE__,
		    "%s: Syntax error in specifying MC clustersize parameters!",TBname);
    cls3_[0] = v[0]; cls3_[1] = v[1]; cls3_[2] = v[2]; cls3_[3] = v[3]; cls3_[4] = v[4]; cls3_[5] = v[5];
  }

  key ="clsParsMC";
  if (CsOpt::Instance()->getOpt(TBname,key,v) ||
      CsOpt::Instance()->getOpt(SI,key,v)) {
    if (v.size()!=7)
      CsErrLog::msg(elFatal,__FILE__,__LINE__,
		    "%s: Syntax error in specifying MC clustersize parameters!",TBname);
    cls_[0] = v[0]; cls_[1] = v[1]; cls_[2] = v[2]; cls_[3] = v[3]; cls_[4] = v[4]; cls_[5] = v[5];   cls_[6] = v[6];
  }


  key ="etaParsMC";
  if (CsOpt::Instance()->getOpt(TBname,key,v) ||
      CsOpt::Instance()->getOpt(SI,key,v)) {
    if (v.size()!=8)
      CsErrLog::msg(elFatal,__FILE__,__LINE__,
		    "%s: Syntax error in specifying MC eta parameters!",TBname);
    eta_[0] = v[0]; eta_[1] = v[1]; eta_[2] = v[2]; eta_[3] = v[3];
    eta_[4] = v[4]; eta_[5] = v[5]; eta_[6] = v[6]; eta_[7] = v[7];
  }
}


// ------------------------------------------------------------
void CsSiTrackerDetector::BookHistograms() {
// ------------------------------------------------------------

  // check if histograms are to be booked
  CsDetector::ReadHistLevel();
  if (hLevel_==None) return;

  const string &name=GetTBName();

  CsHistograms::SetCurrentPath("/Silicon");

  TString sn(name);

  switch (hLevel_) {
  case High:

    T_SIdeb_Ccl = new TTree(sn+"SI_CINDARELLA_cluster","example");

    T_SIdeb_Ccl->Branch("tcs_cor_c", &tcs_cor_c, "tcs_cor_c");
    T_SIdeb_Ccl->Branch("time_deb",&time_deb,"time_deb/F");
    T_SIdeb_Ccl->Branch("time_deb_t0",&time_deb_t0,"time_deb_t0/F");
    T_SIdeb_Ccl->Branch("time_deb_t1",&time_deb_t1,"time_deb_t1/F");
    T_SIdeb_Ccl->Branch("time_clus",&time_clus,"time_clus/F");
    T_SIdeb_Ccl->Branch("time_clus_t0",&time_clus_t0,"time_clus_t0/F");
    T_SIdeb_Ccl->Branch("time_clus_t1",&time_clus_t1,"time_clus_t1/F");
    T_SIdeb_Ccl->Branch("time_out",&time_out,"time_out/F");
    T_SIdeb_Ccl->Branch("pos_out",&pos_out,"pos_out/F");
    T_SIdeb_Ccl->Branch("cluster_size_out",&cluster_size_out,"cluster_size_out/I");
    T_SIdeb_Ccl->Branch("multCut_thiemo",&i_cl,"i_cl/I");
    T_SIdeb_Ccl->Branch("amp2_out", &amp2_out, "amp2_out");
    T_SIdeb_Ccl->Branch("amp1_out", &amp1_out, "amp1_out");
    T_SIdeb_Ccl->Branch("amp0_out", &amp0_out, "amp0_out");
    T_SIdeb_Ccl->Branch("r1_out", &r1_out, "r1_out");
    T_SIdeb_Ccl->Branch("r0_out", &r0_out, "r0_out");

    T_SIdeb_Ccl->Branch("clustertime_err_out", &clustertime_err_out, "clustertime_err_out");

    h(name+"_ch_vs_amp",getNWir(),0,getNWir(),300,0,600); // Channel vs Amp
    h(name+"_r0",35,-5,30,150,0,1.5);                     // amplitude ratio A0/A2
    h(name+"_r1",35,-5,30,150,0,1.5);                     // amplitude ratio A1/A2
    h(name+"_time0",130,-50,80);                          // time from ratio r0
    h(name+"_time1",130,-50,80);                          // time from ratio r1
    // common noise (all chips), only even numbers
    h(name+"_CoNo",ChipsPerPlane*200,150+200,ChipsPerPlane*400+200+150);
    h(name+"_tcsPhase",80,-10.3,30.9);                    // trigger wrt TCS clock
    h(name+"_time_deb",80,-50,80);
    h(name+"_time_deb_t0",80,-50,80);
    h(name+"_time_deb_t1",80,-50,80);

    h(name+"_time_clus",80,-50,80);
    h(name+"_time_clus_t0",80,-50,80);
    h(name+"_time_clus_t1",80,-50,80);

    h(name+"_time_out",130,-50,80);
    h(name+"_pos_out",400,0,1200);
    h(name+"_cluster_size_out",20,0,20);


    // to avoid problems with binning it is recommended to use
    // multiple of the time_unit of the time measuremnt:
    //
    // F1-unit = 0.1289231  ns in standard mode
    // F1-unit = 0.06446156 ns in high resolution mode
    // F1-unit = 4.899079   ns in latch mode (MWPC)
    //

    // no break here!

  case Normal:
    h(name+"_ch",getNWir(),0,getNWir());          // hits per strip
    h(name+"_strips",200,0,200);                  // number of strips above threshold
    h(name+"_cAmp",300,0,600);                    // cluster amplitude
    h(name+"_cAmpCut",300,0,600);                 // cluster amplitude with time cut
    h(name+"_cPos",getNWir()/4,0,getNWir());      // cluster position
    h(name+"_cSize",10,0,10);                     // cluster size
    h(name+"_cMult",100,0,100);                   // cluster multiplicity old clust.
    h(name+"_cMultCut",100,0,100);                // cluster multiplicity after time cut old clust.
    h(name+"_cMultCut_thiemo",100,0,100);         // cluster multiplicity after thiemos algorithm
    h(name+"_cTime",100,-50,50);                  // cluster time

    // for mate clustering
    h(name+"_csplitTime ",100,-50,50);			  // cluster time of splitted clusters
    h(name+"Nsplit",10,0,10);				  // number of new added clusters
    h(name+"csplitPos",500,-40,40);			  // position of new added clusters
    h(name+"cAmpsplit",250,0,500);				  // cluster amp of new added clusters
    h(name+"cSizesplit",10,0,10);				  // cluster size of new added clusters
    h(name+"cA2_1",250,0,500);
    h(name+"cA2_2",250,0,500);
    h(name+"cA2_wide",250,0,500);
    h(name+"cA2_diff",500,-500,500);
    h(name+"findbestcalled",3,0,3);
    h(name+"findbestsuccess",3,0,3);
    h(name+"a21:a22", 250,0,500,250,0,500);
    h(name+"diff_vs_a2sum", 250,-250,250,250,0,500);
    h(name+"large_vs_a2sum", 250,0,500,250,0,500);


    h(name+"3strip",10,0,10 );

    h(name+"_finalizePos",500,-40,40);

    h(name+"_delta_called",2,0,2);
    h(name+"_delta_A2_1",250,0,500);
	h(name+"_delta_A2_2",250,0,500);
	h(name+"_delta_A2_diff",500,-500,500);
    h(name+"_delta_cls_1",10,0,10);
	h(name+"_delta_cls_2",10,0,10);
    h(name+"_delta_cls_in",30,0,30);
    h(name+"_delta_pos_in",1000,-40,40);
	h(name+"_delta_pos_1",400,-40,40);
	h(name+"_delta_pos_2",400,-40,40);
	h(name+"_ndigits",40,0,40);
	h(name+"after_mate",40,0,40);

   h(name+"mc_clusteramp", 400,0,400);
   h(name+"mc_eta", 100,0,1);
   h(name+"mc_eta_pos", 200, 0, 60, 100, 0, 1);
   h(name+"mc_dist", 200,0,60);
   h(name+"mc_dist1", 200,0,60);
   h(name+"mc_dist2", 200,0,60);
   h(name+"mc_cls", 5,0,5);
   h(name+"mc_wire", 1400,0,1400);



    break;
  default:
    return;
  }
}


// ============================================================
// operators

bool CsSiTrackerDetector::operator==( const CsSiTrackerDetector& det ) const {

  return( CsDetector::operator==(det) );
}

bool CsSiTrackerDetector::operator<( const CsSiTrackerDetector& det ) const {
  return( CsDetector::operator<(det) );
}
// ============================================================

// ============================================================
// histo accessors
CsHist1D* CsSiTrackerDetector::Histograms::operator()(string name,int N,double from,double to)
{
  CsHist1D* hist=new CsHist1D(name.c_str(),name.c_str(),N,from,to);
  hist1[name]=hist;
  return hist;
}

CsHist2D* CsSiTrackerDetector::Histograms::operator()(string name,int Nx,double fromx,double tox,int Ny,double fromy,double toy)
{
  CsHist2D* hist=new CsHist2D(name.c_str(),name.c_str(),Nx,fromx,tox,Ny,fromy,toy);
  hist2[name]=hist;
  return hist;
}

void CsSiTrackerDetector::Histograms::operator()(string name,double value)
{
  if(hist1.find(name)!=hist1.end()) {
    hist1[name]->Fill(value);
  }
}

void CsSiTrackerDetector::Histograms::operator()(string name,double x,double y)
{
  if(hist2.find(name)!=hist2.end()) {
    hist2[name]->Fill(x,y);
  }
}
// ============================================================


// --------------------------------------------------
istream & operator >> (istream &in,vector<CsSiTrackerDetector::Calib> &c) {
  // --------------------------------------------------

  while(1)
    {
      char s[111];
      in.getline(s,sizeof(s)); // read line with information about chips
      if( !in )
        break; // end of data stream

      if( s[0]=='#' )
        {
          char s2[sizeof(s)];
          int a,src_id,adc_id;
          unsigned chip=9999999;

          if( 5!=sscanf(s,"%s %d %d %d %d",s2,&a,&src_id,&adc_id,&chip) )
            throw CS::Exception("CsSiTrackerDetector::Calib::operator>>: bad string \"%s\".",s);

          // mapping: (adc_id, chip id) ---> consecutive chip #
          if(adc_id % 2 == 0){
            chip= chip-3;
          }else{
            chip= 10-chip;
          }

          if( chip >= c.size() )
            throw CS::Exception("CsSiTrackerDetector::Calib::operator>>: "
                                "Bad chip number %d, min=0, max=%d", chip, c.size());

          c[chip].src_id      = src_id;
          c[chip].adc_id      = adc_id;
          c[chip].initialised = true;

          for( int i=0; i<CsSiTrackerDetector::ChipChannels; i++ )
            {
              in.getline(s,sizeof(s)); // read line with information about chips

              if( !in )
                throw CS::Exception("CsSiTrackerDetector::Calib::operator>>: Error in reading data.");

              c[chip].channels.push_back(CsSiTrackerDetector::Calib::Channel(s));
            }

          // Debugging output
          if (getenv("COND_DB_DEBUG")!=0 || getenv("COND_DB_DEBUG_SI")!=0) {
            vector<CsSiTrackerDetector::Calib>::iterator it;
            vector<CsSiTrackerDetector::Calib::Channel>::iterator itc;
            for (it = c.begin(); it != c.end(); it++) {
              cout<<"next Calib structure:\n";
              vector<CsSiTrackerDetector::Calib::Channel>& vch = (*it).channels;
              cout << " vch.size() " << vch.size()  << endl;
              for (itc = vch.begin(); itc != vch.end(); itc++) {
                (*itc).Print();
              }
            }
          }
        }
      else
        throw CS::Exception("CsSiTrackerDetector::Calib::operator>>: format error.");
    }//end while loop, sorting the chips and channels per plane

  in.clear(); // Why should I do this? I don't know...
  return in;
} //end istream operator >>



// ------------------------------------------------------------
void CsSiTrackerDetector::readCalibration(time_t timePoint){
// ------------------------------------------------------------

  if (getenv("COND_DB_DEBUG")!=0 || getenv("COND_DB_DEBUG_SI")!=0)
    cout << "CsSiTrackerDetector[" << GetTBName() << "]::readCalibration()" << endl;

  // create space for per-chip calibration data
  for (size_t i=0; i<ChipsPerPlane; i++) {
    calib_data.push_back(CsSiTrackerDetector::Calib());
  }

  CDB::Time tp(timePoint,0);


  // Read calibration data for each strip
  string strdata("");
  cdb_->read(GetTBName(), strdata, tp);

  istringstream istrdata(strdata);

  // push calibration string into calibration vector
  istrdata >> calib_data;

  has_calib_data = true;

  CsErrLog::msg(elInfo,__FILE__,__LINE__,
                "%s yes! calibration data found in CDB for local time ????",
                GetTBName().c_str());


  // Debugging output
  if (getenv("COND_DB_DEBUG")!=0 || getenv("COND_DB_DEBUG_SI")!=0) {
    cout << GetTBName() <<endl;
    cout << "calib_data.size() " << calib_data.size() << endl;
    vector<CsSiTrackerDetector::Calib>::iterator itcal;
    int chip;
    for (itcal = calib_data.begin(), chip=0; itcal != calib_data.end(); itcal++, chip++) {
      cout << "calib_data[" << chip << "]:" << endl;
      if (!itcal->IsInitialised()) {
        cout << "NOT INITIALISED" << endl;
        continue;
      }
      const vector<CsSiTrackerDetector::Calib::Channel>& vchn = itcal->GetChannels();
      cout << " vchn.size() " << vchn.size() << endl;
      vector<CsSiTrackerDetector::Calib::Channel>::const_iterator itchn;
      for (itchn = vchn.begin(); itchn != vchn.end(); itchn++) {
        itchn->Print();
      }
    }
    cout<<endl;
  }


  //---------------------------------------
  // read timing calibrations
  if (time_12nsjumps){
    CsErrLog::msg(elInfo,__FILE__,__LINE__,
		  "%s timing: Correcting for 12ns jumps by spill!",
		  GetTBName().c_str());

    string filein_tjumps ="/nfs/hicran/project/compass/silicon/time_calib/timingfiles/"+calib_dir+"/"
      +GetTBName()+"12nsjumps";
    std::ifstream t_jumps(filein_tjumps.c_str());  //in case of TimeCalibMunich on, this will be used

    if (!t_jumps||time_calib_db||!time_calib_munich){    //otherwise (or in case file does not exist), read from DB
      CsErrLog::msg(elInfo,__FILE__,__LINE__,
		    "%s: reading time jumps from Database",GetTBName().c_str());
      string strdata3("");
      cdb_->read(GetTBName(), strdata3, tp, "12nsjumps");
      istringstream istrdata3(strdata3);
      istrdata3 >> timing_spill_jumps;
    }
    else{
      t_jumps >>  timing_spill_jumps;
      CsErrLog::msg(elInfo,__FILE__,__LINE__,
		    "%s reading timing jumps from file %s, \n%u integers read",
		    GetTBName().c_str(),filein_tjumps.c_str(), timing_spill_jumps.size());


		    }
    if (timing_spill_jumps.size()==0){
      CsErrLog::msg(elInfo,__FILE__,__LINE__,
		    "%s: while reading timing jumps: No jumps were found (assuming zero for all spills)",
		    GetTBName().c_str(),filein_tjumps.c_str(), timing_spill_jumps.size());

    }

  }


  string filein = "/nfs/hicran/project/compass/silicon/time_calib/timingfiles/"+calib_dir+"/"+GetTBName();
  std::ifstream t_calib(filein.c_str());

  if (!t_calib&&time_calib_munich){
    CsErrLog::msg(elFatal,__FILE__,__LINE__,
		  "%s NO time calibrations are found in %s \n YOU TRY TO READ FROM FILE!! Set option SI TimeCalibMunich off;",GetTBName().c_str(),filein.c_str());
  }

  if (!t_calib||time_calib_db||!time_calib_munich){
    CsErrLog::msg(elInfo,__FILE__,__LINE__,
		  "%s reading timing files from Database",GetTBName().c_str());
    //cout << "in der if_Bedingung, reading from DB " << GetTBName() << endl;
    string strdata2("");
    cdb_->read(GetTBName(), strdata2, tp, "timing");
    istringstream istrdata2(strdata2);
    istrdata2 >> timing_calib;
  }
   else{
     t_calib >>  timing_calib;
     //cout << "in der else_Bedingung, reading from file " << GetTBName() << endl;
     CsErrLog::msg(elError,__FILE__,__LINE__,
		  "%s reading time calibrations from file %s, \n ONLY USE FOR TESTING NEW TIME CALIBRATIONS!!;",
		   GetTBName().c_str(),filein.c_str());
   }

  bool cindCalibration(false);
  if (timing_calib.size()>10){
    CsErrLog::msg(elInfo,__FILE__,__LINE__,
		  "%s time calibrations for cinderella clustering are found;",
		  GetTBName().c_str());
   cindCalibration=true;
   }

  if (timing_calib.size()<=10){
    CsErrLog::msg(elInfo,__FILE__,__LINE__,
		  "%s only time calibrations for ratio clustering are found;",
		  GetTBName().c_str());
   }

  if (cluster_cind_opt && !cindCalibration) 
    CsErrLog::msg(elFatal, __FILE__, __LINE__,
                  "CsSiTrackerDetecor::readCalibration: %s: cinderella clustering requested, but not suitable time calibrations found.",
                  GetTBName().c_str());




  if (timing_calib.size()!=0) {
    CsErrLog::msg(elInfo,__FILE__,__LINE__,
		  "%s calibration OK.",GetTBName().c_str());

    // Just debug printouts
    if (getenv("COND_DB_DEBUG")!=0 || getenv("COND_DB_DEBUG_SI")!=0) {
      const double& r1_min = this->timing_calib[0];
      const double& r1_max = this->timing_calib[1];
      const double& r2_min = this->timing_calib[2];
      const double& r2_max = this->timing_calib[3];
      const double& t0_1   = this->timing_calib[4];
      const double& t0_2   = this->timing_calib[5];
      const double& sl_1   = this->timing_calib[6];
      const double& sl_2   = this->timing_calib[7];
      const double& tc_1   = this->timing_calib[8];
      const double& tc_2   = this->timing_calib[9];

      p_ratio1.r0 = this->timing_calib[12];
      p_ratio1.t0 = this->timing_calib[13];
      p_ratio1.a  = this->timing_calib[14];
      p_ratio1.b  = this->timing_calib[15];
      p_ratio1.c  = this->timing_calib[16];
      p_ratio1.d  = this->timing_calib[17];

      p_ratio2.r0 = this->timing_calib[20];
      p_ratio2.t0 = this->timing_calib[21];
      p_ratio2.a  = this->timing_calib[22];
      p_ratio2.b  = this->timing_calib[23];
      p_ratio2.c  = this->timing_calib[24];
      p_ratio2.d  = this->timing_calib[25];

      cout << endl <<"#######Timing calibrations read in for "<< GetTBName()<< endl;
      cout << GetTBName() << endl;
      cout << "r1_min = " << r1_min << endl;
      cout << "r1_max = " << r1_max << endl;
      cout << "r2_min = " << r2_min << endl;
      cout << "r2_max = " << r2_max << endl;
      cout << "t0_1   = " << t0_1   << endl;
      cout << "t0_2   = " << t0_2   << endl;
      cout << "sl_1   = " << sl_1   << endl;
      cout << "sl_2   = " << sl_2   << endl;
      cout << "tc_1   = " << tc_1   << endl;
      cout << "tc_2   = " << tc_2   << endl;

      cout << "p_ratio1.r0 " << p_ratio1.r0 << endl;
      cout << "p_ratio1.t0 " << p_ratio1.t0 << endl;
      cout << "p_ratio1.a  " << p_ratio1.a  << endl;
      cout << "p_ratio1.b  " << p_ratio1.b  << endl;
      cout << "p_ratio1.c  " << p_ratio1.c  << endl;
      cout << "p_ratio1.d  " << p_ratio1.d  << endl;
	      				
      cout << "p_ratio2.r0 " << p_ratio2.r0  << endl;
      cout << "p_ratio2.t0 " << p_ratio2.t0  << endl;
      cout << "p_ratio2.a  " << p_ratio2.a   << endl;
      cout << "p_ratio2.b  " << p_ratio2.b   << endl;
      cout << "p_ratio2.c  " << p_ratio2.c   << endl;
      cout << "p_ratio2.d  " << p_ratio2.d   << endl;

      cout << "timing_calib.size: " <<timing_calib.size()<<endl;
      cout << endl;

    }
  } else {
    tm *t = localtime(&tp.first);
    cout << GetTBName()+"timing"<< ", no calibration for local time "
	 << t <<" in CDB"<< endl;
    CsErrLog::msg(elFatal,__FILE__,__LINE__,
		  "%s No time calibration.",GetTBName().c_str());
  }
  /* Charge fraction position corrections for 2-strip-clusters */
  eta_corr.clear();
  if (!eta_correction){ // if switched off in option file (used to create new calibrations)
    eta_corr.push_back(0.);
    eta_corr.push_back(0.);
    eta_corr.push_back(0.);
    eta_corr.push_back(0.);
    eta_corr.push_back(0.);
    eta_corr.push_back(0.);
  }
  else if (eta_calib){ //read from calib, default
    string strdata4;
    cdb_->read(GetTBName(), strdata4, tp, "eta_correction");
    istringstream istrdata4(strdata4);
    istrdata4 >> eta_corr;
    if (eta_corr.size()!=0) while (eta_corr.size()<6) eta_corr.push_back(0.); // for shorter polynomials, just put zeroes for the last coefficients
  }

  if (eta_corr.size()==0 && eta_correction){  //if no calibration was found, or "SI EtaCorrection hard" given in the option file
    if (TString(GetTBName()).Contains("1U")){
      eta_corr.push_back(18.1812);
      eta_corr.push_back(-342.077);
      eta_corr.push_back(1762.25);
      eta_corr.push_back(-3997.68);
      eta_corr.push_back(4248.65);
      eta_corr.push_back(-1707.42);
    }
    else if (TString(GetTBName()).Contains("1V")){
      eta_corr.push_back(5.52724);
      eta_corr.push_back(-141.362);
      eta_corr.push_back(914.437);
      eta_corr.push_back(-2376.42);
      eta_corr.push_back(2700.69);
      eta_corr.push_back(-1110.55);
    }
    else if (TString(GetTBName()).Contains("1X")){
      eta_corr.push_back(-0.598299);
      eta_corr.push_back(-144.379);
      eta_corr.push_back(844.126);
      eta_corr.push_back(-1945.36);
      eta_corr.push_back(2078.34);
      eta_corr.push_back(-836.669);
    }
    else if (TString(GetTBName()).Contains("1Y")){
      eta_corr.push_back(11.5317);
      eta_corr.push_back(-229.534);
      eta_corr.push_back(1331.18);
      eta_corr.push_back(-3233.54);
      eta_corr.push_back(3505.67);
      eta_corr.push_back(-1394.97);
    }
    else if (TString(GetTBName()).Contains("2U")){
      eta_corr.push_back(17.1205);
      eta_corr.push_back(-298.721);
      eta_corr.push_back(1791.42);
      eta_corr.push_back(-4478.7);
      eta_corr.push_back(4902.9);
      eta_corr.push_back(-1946.1);
    }
    else if (TString(GetTBName()).Contains("2V")){
      eta_corr.push_back(9.60402);
      eta_corr.push_back(-222.222);
      eta_corr.push_back(1372.38);
      eta_corr.push_back(-3542.4);
      eta_corr.push_back(4003.88);
      eta_corr.push_back(-1638.24);
    }
    else if (TString(GetTBName()).Contains("2X")){
      eta_corr.push_back(0);
      eta_corr.push_back(0);
      eta_corr.push_back(0);
      eta_corr.push_back(0);
      eta_corr.push_back(0);
      eta_corr.push_back(0);
    }
    else if (TString(GetTBName()).Contains("2Y")){
      eta_corr.push_back(12.217);
      eta_corr.push_back(-279.923);
      eta_corr.push_back(1712.22);
      eta_corr.push_back(-4329.05);
      eta_corr.push_back(4780.66);
      eta_corr.push_back(-1912.03);
    }
    else if (TString(GetTBName()).Contains("3U")){
      eta_corr.push_back(10.5748);
      eta_corr.push_back(-280.52);
      eta_corr.push_back(1452.93);
      eta_corr.push_back(-3217.4);
      eta_corr.push_back(3361.43);
      eta_corr.push_back(-1336.36);
    }
    else if (TString(GetTBName()).Contains("3V")){
      eta_corr.push_back(7.60036);
      eta_corr.push_back(-180.71);
      eta_corr.push_back(1073.38);
      eta_corr.push_back(-2671.91);
      eta_corr.push_back(2970.55);
      eta_corr.push_back(-1209.47);
    }
    else if (TString(GetTBName()).Contains("3X")){
      eta_corr.push_back(15.4977);
      eta_corr.push_back(-310.73);
      eta_corr.push_back(1627.39);
      eta_corr.push_back(-3687.75);
      eta_corr.push_back(3915.24);
      eta_corr.push_back(-1574.9);
    }
    else if (TString(GetTBName()).Contains("3Y")){
      eta_corr.push_back(13.7328);
      eta_corr.push_back(-309.451);
      eta_corr.push_back(1776.31);
      eta_corr.push_back(-4307.52);
      eta_corr.push_back(4661.97);
      eta_corr.push_back(-1851.6);
    }
    else if (TString(GetTBName()).Contains("4U")){
      eta_corr.push_back(18.3886);
      eta_corr.push_back(-306.2);
      eta_corr.push_back(1667.81);
      eta_corr.push_back(-4009.04);
      eta_corr.push_back(4323.76);
      eta_corr.push_back(-1714.73);
    }
    else if (TString(GetTBName()).Contains("4V")){
      eta_corr.push_back(11.1792);
      eta_corr.push_back(-253.678);
      eta_corr.push_back(1542.16);
      eta_corr.push_back(-3956.47);
      eta_corr.push_back(4457.71);
      eta_corr.push_back(-1820.99);
    }
    else if (TString(GetTBName()).Contains("4X")){
      eta_corr.push_back(15.4242);
      eta_corr.push_back(-256.929);
      eta_corr.push_back(1568.61);
      eta_corr.push_back(-3974.14);
      eta_corr.push_back(4399.6);
      eta_corr.push_back(-1764.97);
    }
    else if (TString(GetTBName()).Contains("4Y")){
      eta_corr.push_back(14.1581);
      eta_corr.push_back(-215.118);
      eta_corr.push_back(1229.91);
      eta_corr.push_back(-2964.85);
      eta_corr.push_back(3165.16);
      eta_corr.push_back(-1236.36);
    }
    else if (TString(GetTBName()).Contains("5U")){
      eta_corr.push_back(19.3904);
      eta_corr.push_back(-307.173);
      eta_corr.push_back(1676.29);
      eta_corr.push_back(-4031.59);
      eta_corr.push_back(4355.18);
      eta_corr.push_back(-1732.27);
    }
    else if (TString(GetTBName()).Contains("5V")){
      eta_corr.push_back(10.7923);
      eta_corr.push_back(-184.352);
      eta_corr.push_back(1131.05);
      eta_corr.push_back(-2922.87);
      eta_corr.push_back(3305.16);
      eta_corr.push_back(-1352.32);
    }
    else if (TString(GetTBName()).Contains("5X")){
      eta_corr.push_back(17.0162);
      eta_corr.push_back(-241.83);
      eta_corr.push_back(1306.21);
      eta_corr.push_back(-3130.99);
      eta_corr.push_back(3381.23);
      eta_corr.push_back(-1347.57);
    }
    else if (TString(GetTBName()).Contains("5Y")){
      eta_corr.push_back(13.3911);
      eta_corr.push_back(-256.382);
      eta_corr.push_back(1461.16);
      eta_corr.push_back(-3533.71);
      eta_corr.push_back(3781.12);
      eta_corr.push_back(-1478.4);
    }
  }

}



// Setup Sergei's alignment patch for the current run number, if not already
// done so.  According to Yann, it is possible that the run number changes
// within the same invocation of CORAL.
// --------------------------------------------------
void CsSiTrackerDetector::SetupAlignmentPatch() {
// --------------------------------------------------
  if (!align_patch) {
    rescor = 0.;
    return;
  }

  unsigned int runnb;
  runnb = CsEvent::Instance()->getRunNumber();

  // check whether we have already setup the patch for this run number
  if (runnb == align_patch_run) {
    return;
  }

  if (runnb < 42320 || 43350 < runnb) {
    rescor = 0.;
    CsErrLog::msg(elInfo, __FILE__, __LINE__,
                  "%s: No alignment patch foreseen for run %u", GetTBName().c_str(), runnb);
    return;
  }

  // find reference run
  int refrun = 0;
  int diffrun = -1000000;
  for( int rn=0; rn < ResN; rn++){
    const int thisrun = (unsigned int) SIResCor[rn][0]; //runnumber of correction
    const int thisdiffrun = runnb-thisrun;
    if( abs( thisdiffrun)< abs(diffrun)){
      diffrun = runnb - thisrun ;
      refrun=rn;
    }
  }

  if (refrun == 0)
    throw CS::Exception("%s %i: Alignment patch: No reference run for %i!",
                        __FILE__, __LINE__, runnb);

  if      (TString(GetTBName()).Contains("SI01X")) {rescor = SIResCor[refrun][2];}
  else if (TString(GetTBName()).Contains("SI01Y")) {rescor = SIResCor[refrun][3];}
  else if (TString(GetTBName()).Contains("SI01U")) {rescor = SIResCor[refrun][4];}
  else if (TString(GetTBName()).Contains("SI01V")) {rescor = SIResCor[refrun][5];}
  else if (TString(GetTBName()).Contains("SI02X")) {rescor = SIResCor[refrun][6];}
  else if (TString(GetTBName()).Contains("SI02Y")) {rescor = SIResCor[refrun][7];}
  else if (TString(GetTBName()).Contains("SI02U")) {rescor = SIResCor[refrun][8];}
  else if (TString(GetTBName()).Contains("SI02V")) {rescor = SIResCor[refrun][9];}
  else if (TString(GetTBName()).Contains("SI03X")) {rescor = SIResCor[refrun][10];}
  else if (TString(GetTBName()).Contains("SI03Y")) {rescor = SIResCor[refrun][11];}
  else if (TString(GetTBName()).Contains("SI03U")) {rescor = SIResCor[refrun][12];}
  else if (TString(GetTBName()).Contains("SI03V")) {rescor = SIResCor[refrun][13];}
  else if (TString(GetTBName()).Contains("SI04X")) {rescor = SIResCor[refrun][14];}
  else if (TString(GetTBName()).Contains("SI04Y")) {rescor = SIResCor[refrun][15];}
  else if (TString(GetTBName()).Contains("SI04U")) {rescor = SIResCor[refrun][16];}
  else if (TString(GetTBName()).Contains("SI04V")) {rescor = SIResCor[refrun][17];}
  else if (TString(GetTBName()).Contains("SI05X")) {rescor = SIResCor[refrun][18];}
  else if (TString(GetTBName()).Contains("SI05Y")) {rescor = SIResCor[refrun][19];}
  else if (TString(GetTBName()).Contains("SI05U")) {rescor = SIResCor[refrun][20];}
  else if (TString(GetTBName()).Contains("SI05V")) {rescor = SIResCor[refrun][21];}
  else throw CS::Exception("%s %i: %s: No alignment patch for that plane!",
                           __FILE__, __LINE__, GetTBName().c_str());

    /*
      Extra corrections which takes into account changes in pivot planes (03X, 03Y, 05X, 05Y) positions
      with respect to whole spectrometer after global alignment which had been done for runs 42789 and 43035

      Pivot planes positions extracted from corresponding detectors.dat [cm]:

      /afs/cern.ch/compass/detector/geometry/2004/detectors.42655.hadron.dat
      (this is global alignment which had been used to perform internal alignment
      according to Anna-Maria mail of 20.08.07)

      SI03X  -202.2550    0.63710    0.08650
      SI03Y  -202.2570    0.63710    0.08650

      SI05X  -100.5500    0.84598+   0.10855
      SI05Y  -100.5520    0.84598    0.10855+

      /afs/cern.ch/compass/detector/geometry/2004/detectors.42789.hadron.dat

      SI03X  -202.2550    0.63710    0.08650
      SI03Y  -202.2570    0.63710    0.08650

      SI05X  -100.5500    0.84563*   0.10880
      SI05Y  -100.5520    0.84563    0.10880*

      /afs/cern.ch/compass/detector/geometry/2004/detectors.43035.hadron.dat

      SI03X  -202.2550    0.63710    0.08650
      SI03Y  -202.2570    0.63710    0.08650

      SI05X  -100.5500    0.84553*   0.10946
      SI05Y  -100.5520    0.84553    0.10946*

      (+) - reference offsets of pivot planes
      (*) - changed position w.r.t. spectrometer


     as position of pivot planes 03X and 03Y are the same, the only transformation
     which has to be applied is rotation of all Si planed around Si03 position on the angle
     defined by global shift of other pivot planes 05X and 05Y

                                                              ges@mail.cern.ch
    */

  double alphaX = 0;
  double alphaY = 0;

  if(runnb >= 42789 && runnb < 43035) {
    alphaX = (0.84563 - 0.84598) / (-100.5500 + 202.2550);
    alphaY = (0.10880 - 0.10855) / (-100.5520 + 202.2570);
  }
  if(runnb >= 43035 && runnb <= 43348) {
    alphaX = (0.84553 - 0.84598) / (-100.5500 + 202.2550);
    alphaY = (0.10946 - 0.10855) / (-100.5520 + 202.2570);
  }

  double z0 = -2022.560;
  double z  = this->getZcm();
  double dx = (z-z0)*alphaX;
  double dy = (z-z0)*alphaY;
  //printf("%s  Z = %12.5f  dx = %12.9f  dy = %12.9f \n", this->GetTBName().c_str(), z, dx, dy);
  // transform correction to WRS
  HepMatrix mW;
  mW=this->getRotWRS();
  double sina = -mW[0][1];
  double cosa =  mW[1][1];
  double du   =  cosa*dx + sina*dy;
  //printf("rotation = %12.9f  du = %12.9f \n", atan(sina/cosa)*180./3.14159265358979312, du);

  rescor -= du; // apply rotation correction

  align_patch_run = runnb;
}



// Setup the variables params_noise[1300], params_pedestal[1300],
// pedestal_mean_chip[10] which are used in Cinderella clustering, if not
// already done so.
// --------------------------------------------------
void CsSiTrackerDetector::SetupCinderellaClustering() {
// --------------------------------------------------


  if (!(CsInit::Instance()->IsAMonteCarloJob())){
  const unsigned nwires = ChipsPerPlane * ChipChannels;

  double *initialised = new double[nwires];     // true if wire has been initialised
  for (unsigned int i = 0; i < nwires; i++)
    initialised[i] = false;

  // check decoding maps
  const CS::Chip::Maps &daq_maps = CsInit::Instance()->getDaqMaps();
  if (daq_maps.size() == 0) {
    throw CS::Exception("%s %i: %s: Daq maps not found!",
                        __FILE__, __LINE__, GetTBName().c_str());
  }

  // Check if calibrations were loaded
  if (!has_calib_data) {
    CsErrLog::msg(elError, __FILE__, __LINE__,
                  "%s: No calibrations found!", GetTBName().c_str());
  } else {

    // Loop over all 8-10 chips and get calibration data for each chip
    for (unsigned chip_count=0; chip_count<calib_data.size(); chip_count++) {
      Calib &c = calib_data[chip_count];

      // check whether calibrations for that chip were loaded ok
      if (!c.IsInitialised()) {
        CsErrLog::msg(elWarning, __FILE__, __LINE__,
                      "%s: No calibration for chip %i!",
                      GetTBName().c_str(), chip_count);
        continue;
      }
	
      // Get Source ID and ADC ID of chip
      const int src_id=c.GetSrcId();
      const int adc_id=c.GetAdcId();

      if (getenv("COND_DB_DEBUG")!=0 || getenv("COND_DB_DEBUG_SI")!=0) {
        cout << "Looking for source Id " << src_id
             << ", ADC Id " << adc_id
             << ", chip Id " << chip_count << endl;
      }

      // reverse mapping: (adc_id, consecutive chip #) ---> chip_id
      unsigned chip = chip_count;
      if(adc_id % 2 == 0){
        chip= chip+3; // adc id is even, chip number starts with 3,4,5,6
      }else{
        chip= 10-chip; // adc_id odd, chip number starts with 6,5,4,3,2,1
      }

      assert(chip < ChipsPerPlane);

      if (getenv("COND_DB_DEBUG")!=0 || getenv("COND_DB_DEBUG_SI")!=0) {
        cout << "Looking for source Id " << src_id
             << ", ADC Id " << adc_id
             << ", chip Id " << chip <<" after mod" << endl;
      }	


      // Create a short name for iterator through map
      typedef CS::Chip::Maps::const_iterator m_it;

      // Loop over all channels of chip and find the associated wire
      for (unsigned chip_chan=0; chip_chan<c.GetChannels().size(); chip_chan++) {

        // All maps with given data ID
        CS::ChipAPV::DataID data_id(src_id, adc_id, chip, chip_chan);

        //  cout << "Data Id " << data_id << endl;
        const pair<m_it,m_it> m_range = daq_maps.equal_range(data_id);

        // Read map
        bool mapped=false;
        unsigned int wire;
        for( m_it cc=m_range.first; cc!=m_range.second; cc++ ) {
          const CS::ChipAPV::Digit *digit1 =
            dynamic_cast<CS::ChipAPV::Digit*>(cc->second);

          if( digit1==NULL ){
            CsErrLog::msg(elWarning,__FILE__,__LINE__,
                          "%s: ChipAPV wrong map!",GetTBName().c_str());
            continue;
          }

          if( GetTBName()!=digit1->GetDetID().GetName() ){
            CsErrLog::msg(elError,__FILE__,__LINE__,
                          "%s: Inconsistency between mapping file and calibration file!",
                          GetTBName().c_str());
            CsErrLog::msg(elError,__FILE__,__LINE__,
                          "%s from calibration file does not match TBName",
                          digit1->GetDetID().GetName().c_str());
            continue;
          }

          // Apply mapping
          try {
            wire=digit1->GetChannel();
          }
          catch(...) {
            CsErrLog::msg(elError,__FILE__,__LINE__,
                          "%s: Could not apply mapping!",
                          GetTBName().c_str());
            continue;
          }
          mapped=true;
          break;
        } // Loop over maps

        // Assert correct mapping
        if (!mapped)
          throw CS::Exception("%s %i: %s: No mapping for "
                              "srcid %i, adcid %i, chip %i, chipchan %i",
                              __FILE__, __LINE__, GetTBName().c_str(),
                              src_id, adc_id, chip, chip_chan);

        if (nwires <= wire)
          throw CS::Exception("%s %i: %s: Wire %i is outside valid range: "
                              "srcid %i, adcid %i, chip %i, chipchan %i",
                              __FILE__, __LINE__, GetTBName().c_str(), wire,
                              src_id, adc_id, chip, chip_chan);

        if (initialised[wire])
          throw CS::Exception("%s %i: %s: Duplicate mapping for wire %i: "
                              "srcid %i, adcid %i, chip %i, chipchan %i",
                              __FILE__, __LINE__, GetTBName().c_str(), wire,
                              src_id, adc_id, chip, chip_chan);

        if (getenv("COND_DB_DEBUG")!=0 || getenv("COND_DB_DEBUG_SI")!=0) {
          cout << GetTBName().c_str() << ": "
               << "src_id: " << src_id
               << "  adc_id: " << adc_id
               << "  chip: " << chip
               << "  chipchan: " << chip_chan
               << "  => wire: " << wire
               << " ..... ";

          cout << setprecision(6) << c.GetChannels()[chip_chan].flag << " ";
          cout << setprecision(6) << c.GetChannels()[chip_chan].pedestal_mean << " ";
          cout << setprecision(6) << c.GetChannels()[chip_chan].pedestal_sigma << " ";
          cout << setprecision(6) << c.GetChannels()[chip_chan].calibration_mean << " ";
          cout << setprecision(6) << c.GetChannels()[chip_chan].calibration_sigma <<endl;
        }

        /* Factor 0.5 is to scale the error which is in raw ADC units to
           match the digits which are lowest-bit truncated ADC units. */
        params_noise[wire]      = 0.5 * c.GetChannels()[chip_chan].pedestal_sigma;
        params_pedestal[wire]   = c.GetChannels()[chip_chan].pedestal_mean;
        initialised[wire]       = true;

      } // Loop over chip channels
    } // Loop over chips
  }

  // Apply dummy calibrations
  unsigned dummy_count = 0;
  for (unsigned int wire=0; wire < nwires; wire++) {
    if (!initialised[wire]) {
      params_noise[wire]      = 2.;
      params_pedestal[wire]   = 750.;
      dummy_count++;
    }
    if (getenv("COND_DB_DEBUG")!=0 || getenv("COND_DB_DEBUG_SI")!=0) {
      cout << GetTBName()
           << " wire " << wire
           << " params_noise " << setprecision(3) << params_noise[wire] << endl;
    }
  }

  // Display warning/error concerning dummy calibrations
  if (0 < dummy_count && dummy_count <= 128)
    CsErrLog::msg(elWarning, __FILE__, __LINE__,
                  "%s: Using dummy calibrations for %i wires!",
                  GetTBName().c_str(), dummy_count);

  if (dummy_count > 128)
    CsErrLog::msg(elError, __FILE__, __LINE__,
                  "%s: Using dummy calibrations for %i wires!",
                  GetTBName().c_str(), dummy_count);

  // calculate for each chip the mean value of the pedestal
  for (unsigned int nbc = 0; nbc < ChipsPerPlane; nbc++) {
    double this_ped_mean = 0;
    for (int imc = 0; imc < ChipChannels; imc++) {
      this_ped_mean += params_pedestal[nbc*ChipChannels + imc];
    }
    pedestal_mean_chip[nbc] = this_ped_mean / 128.;
  }

  delete [] initialised;

  // following code moved from CsSiTrackerDetector::clusterize()

      // variables for struct time_reconstruction_options_t for Cinderellas timing function
      const double TCS_corr = 40;// Calibration constant (Maximum of TCS distr hist.)
      CsEvent* event = CsEvent::Instance();

     /* time_reconstruction_option_t params_time; 		*/					
      params_time.tcsphase = (event->getTCSPhaseTime()-TCS_corr);	

	
      tcs_cor_c =  params_time.tcsphase;


      // see comments to constants in silicon_timing_types.h
      params_time.tcscycle        = 25.77;
      // kind of ratio cut, very little means nearly no cut, neccessary for technical reasons
      params_time.delta           = 0.01;
      params_time.decay_time      = 300.0; // decay_time of a signal in the silicon disable with 100000ns
      params_time.min_early_error = 30; // disable with 0,
      params_time.param[0]        = &p_ratio1;
      params_time.param[1]        = &p_ratio2;
      params_time.region_size_shift[0] = shft0;
      params_time.region_size_shift[1] = shft1;

      if (getenv("COND_DB_DEBUG")!=0 || getenv("COND_DB_DEBUG_SI2")!=0) {
	cout << GetTBName()       << "  "  <<  GetTBName() << "  " << GetTBName() << endl;
	cout << "GetEventnumber in Run : " << event->getEventNumberInRun() << endl;
	cout << setprecision(7)
	     << " tcsphase_raw "  << event->getTCSPhaseTime()
	     << " tcsphase_corr " <<  params_time.tcsphase  << endl;
	cout << "P0: "
	     << "r0 " << params_time.param[0]->r0
	     << "  t0 " << params_time.param[0]->t0
	     << "  a " << params_time.param[0]->a
	     << "  b " << params_time.param[0]->b
	     << "  c " << params_time.param[0]->c
	     << "  d " << params_time.param[0]->d
	     << endl;
	cout << "P1: "
	     << "  r0 " << params_time.param[1]->r0
	     << "  t0 " << params_time.param[1]->t0
	     << "  a " << params_time.param[1]->a
	     << "  b " << params_time.param[1]->b
	     << "  c " << params_time.param[1]->c
	     << "  d " << params_time.param[1]->d
	     << endl;
      }

  }
  else{
    unsigned int nwires = (isU || isX) ? 1280 : 1024;
    double *initialised = new double[nwires];     // true if wire has been initialised

    // in case of Monte Carlo, just dummy pedestal/sigmas
    unsigned dummy_count = 0;
    for (unsigned int wire=0; wire < nwires; wire++) {
      //if (!initialised[wire]) {
	params_noise[wire]      = 1.5;
	params_pedestal[wire]   = 650.;
	initialised[wire]=true;
	dummy_count++;
	//}
    }

    // since we have only dummy calibrations, the mean pedestal value is rather simple
    for (unsigned int nbc = 0; nbc < ChipsPerPlane; nbc++) {
      pedestal_mean_chip[nbc] = 650;
    }

    delete [] initialised;

      // following code moved from CsSiTrackerDetector::clusterize()

     /* time_reconstruction_option_t params_time; 		*/					
      params_time.tcsphase = 0.;	
      tcs_cor_c =  params_time.tcsphase;

      // see comments to constants in silicon_timing_types.h
      params_time.tcscycle        = 25.77;
      // kind of ratio cut, very little means nearly no cut, neccessary for technical reasons
      params_time.delta           = 0.01;
      params_time.decay_time      = 300.0; // decay_time of a signal in the silicon disable with 100000ns
      params_time.min_early_error = 30; // disable with 0,
      params_time.param[0]        = &p_ratio1;
      params_time.param[1]        = &p_ratio2;
      params_time.region_size_shift[0] = 0.;
      params_time.region_size_shift[1] = 0.;

    }


  cind_clus_setup = true;
}


// --------------------------------------------------
void CsSiTrackerDetector::DecodeChipDigit(const CS::Chip::Digit &digit) {
// --------------------------------------------------
  /*
    gets digits from CsDet,
    decodes digits for Silicon
    some checks are performed

    currently only sparse mode is supported
  */

  // Quality check
  // =============


  // 1. check for APV digit
  // ----------------------
  const CS::ChipAPV::Digit *d = dynamic_cast<const CS::ChipAPV::Digit *>(&digit);
  if( d==NULL )
    throw CS::Exception("CsSiTrackerDetector::DecodeChipDigit(): Wrong digit type");

  // 2. reasonable address
  // ---------------------

  if( d->GetChip()==0 ||
      d->GetChip()>ChipsPerPlane ||
      d->GetChipChannel()>=ChipChannels ||
      d->GetChannel() >= int(getNWir()) ) {
    d->Print(cerr,"BAD DIGIT:  ");
    throw
      CS::Exception("CsSiTrackerDetector::DecodeChipDigit()::238: wrong address of digit! \n"
                    "\t\t%s: chip=%d(max=%d) chan=%d(max=%d) wire=%d(max=%d)",GetTBName().c_str(),
		    d->GetChip(),ChipsPerPlane,d->GetChipChannel(),ChipChannels,d->GetChannel(),
		    getNWir());
  }

  if( d->IsSingleFrame() )
    throw CS::Exception("CsSiTrackerDetector::DecodeChipDigit()::244: single frame is not yet supported!");

  if( d->IsSparsifed() )
    sparse=true;
  else
    throw CS::Exception("CsSiTrackerDetector::DecodeChipDigit()::471: found digits without sparsification. Latch all events not yet supported!");


  // Exctract data
  // =============

  int iWire=d->GetChannel();
  //amplitude and commonmode noise correction
  double ampl[12]={d->GetAmplitude()[0],
		   d->GetAmplitude()[1],
		   d->GetAmplitude()[2],
		   -1,1000,3400,
		   double (d->GetCoNo()[0] ),
		   double (d->GetCoNo()[1] ),
		   double (d->GetCoNo()[2] ),
		   0,0,0};//place holder for pedestal and sigma, pedestal_mean

  /* cout <<  iWire  << "common mode correction: "   << " "
       << d-> GetCoNo()[0] << " ampl0: " << d->GetAmplitude()[0] << "  "
       << d-> GetCoNo()[1] << " ampl1: " << d->GetAmplitude()[1] << "  "
       << d-> GetCoNo()[2] << " ampl2: " << d->GetAmplitude()[2] << "  "
       << endl;*/



  // Sort in histograms

  if (hLevel_>=Normal) {              // Fill histograms
    const string &name=GetTBName();
    h(name+"_ch",iWire);
    if (hLevel_==High) {
      h(name+"_ch_vs_amp",iWire,ampl[2]);
      CsEvent* event = CsEvent::Instance();
      int chip = d->GetChannel() / 128;
      if (fEventNr[chip] != event->getEventNumberInRun()) {
	fEventNr[chip] = event->getEventNumberInRun();
	uint16 cono= chip*400 + (d->GetCoNo())[2];
	h(name+"_CoNo",cono);
      }
    }
  }


  // Create CORAL digit
  // ==================
  myDigits_.push_back( new CsDigit(*this, iWire, ampl, 12) );
}



//
// Local class for detector responce simulation
// (Used in CsSiTrackerDetector::makeMCDecoding())
//

// tempopary digits
class SIdig
{
public:

  SIdig(CsDigit* dig) {
    wire=dig->getAddress();
      amp1=dig->getData()[0];
      amp2=dig->getData()[1];
      amp3=dig->getData()[2];
      d = dig;
  }
  SIdig() {
    wire=0;
      amp1=0;
      amp2=0;
      amp3=0;
  }

  int    wire;  // wire #
  double amp1;  // amp1
  double amp2;  // amp2
  double amp3;  // amp3
  CsMCHit* ref; // reference to MCHit
  CsDigit* d;   // reference to the CsDigit
  bool operator < (const SIdig& sd) const
  {
    return (wire < sd.wire);
  };
};


void getMCAmp(const double& t, double amps[] )
{
  // third sampling time
  const float tsamp3= 75.;

  // t1,t2(ns) signal time dependence is A(t)~(1-exp(-t/t1))*exp(-t/t2)
  // for t>0 and A(t)=0. if t<=0. Integral(A(t))[0,inf[=1.
  const double t1= 38.3;
  const double t2= 129.5;
  const double atnorm=(t1+t2)/t2/t2;    // time function Integral^-1

  // timing calibrations , default ones.
  const double t0_1  =  -40.;
  const double sl_1  =  22.;
  const double tc_1  = 1.65;

  const double t0_2  = 0.0;
  const double sl_2  = 27.;
  const double tc_2  = 1.23;

// normalisation of the 3-d amplitude
  double tt=tsamp3-t;
  if(tt<0.) {
    amps[2]=0.;
  }
  else {
    amps[2]=atnorm*(1.-exp(-tt/t1))*exp(-tt/t2);
  }

  // inverted ratios from getGEMtime()
  double r13=tc_1/(exp((t-t0_1)/sl_1)+1.);
  double r23=tc_2/(exp((t-t0_2)/sl_2)+1.);

  amps[0]=amps[2]*r13;
  amps[1]=amps[2]*r23;

  return;
}






// --------------------------------------------------
void CsSiTrackerDetector::makeMCDecoding() {
// --------------------------------------------------
  if (!decode_ && !decodeCard_) return;   // Should I proceed?
  if (decodingDone_) return;              // Already done?
  myDigits_.clear();                      // Clear

const string & name = GetTBName();


  // The following can be done only with zebra (as opposed to NTuple) files...
  if (!CsGeant3::Instance()->isAnNtFile()) {
    double triggerOffset = // Time offset of the trigger w.r.t. the event
      CsEvent::Instance()->getTriggerMCOffset();

    list<CsMCHit*>::iterator Ih;

      //simplistic old MC decoding:
      if (!amplitudeMCDecoding()) {

	for (Ih = myMCHits_.begin(); Ih!=myMCHits_.end(); Ih++) { // loop on hits
	  CsMCTrkHit* thit = dynamic_cast<CsMCTrkHit*>(*Ih); if (!thit) return;
	  if (!thit->getMCTrack()->getParticle()->getCharge() &&
	      !thit->getOrigin()) continue;
	
	  // ***** LOOP ON HITS from charged particles or charged products *****
	
	  double t  = thit->getDTime();  // Delay time (ns)
	  double ui = thit->getUin();    // Hit in point (DRS)
	  double vi = thit->getVin();
	  double wi = thit->getWin();
	  double uo = thit->getUout();   // Hit out point (DRS)
	  double vo = thit->getVout();
	  double wo = thit->getWout();
	  double eloss = thit->getELos()*1.E6; // dE in detector, keV
	


	  t -= // Trigger time (and hence offset) is subtracted from hit times
	    triggerOffset;
	


	  //    ***** CHECK TIME GATE *****
	  // The time gate is taken from the detector table (i.e. so called
	  // "detectors.dat"), contrary to what's done in most other CsDetector's
	  // (where the time gate is specified in the options file). And one
	  // expects it to be quite loose. This situation mimicks the RD case, where
	  // the CsSi class does not apply any time gate.
	  if (t<-tGate_/2 ||  tGate_/2<t) continue;
	  //double tdc = t + tRes_ * CsRandom::gauss();
	  double tdc = t;

      int err;
      HepMatrix iRotM(3,3); iRotM = rotWRS_.inverse( err );
      double xi = rotDRS_(1,1)*ui+rotDRS_(1,2)*vi+rotDRS_(1,3)*wi + xcm_ - _xshift;
      double yi = rotDRS_(2,1)*ui+rotDRS_(2,2)*vi+rotDRS_(2,3)*wi + ycm_ - _yshift;
      double zi = rotDRS_(3,1)*ui+rotDRS_(3,2)*vi+rotDRS_(3,3)*wi + zcm_;
      double xo = rotDRS_(1,1)*uo+rotDRS_(1,2)*vo+rotDRS_(1,3)*wo + xcm_ - _xshift;
      double yo = rotDRS_(2,1)*uo+rotDRS_(2,2)*vo+rotDRS_(2,3)*wo + ycm_ - _yshift;
      double zo = rotDRS_(3,1)*uo+rotDRS_(3,2)*vo+rotDRS_(3,3)*wo + zcm_;

      double Ui = iRotM(1,1)*xi+iRotM(1,2)*yi+iRotM(1,3)*zi; // WRS
      double Vi = iRotM(2,1)*xi+iRotM(2,2)*yi+iRotM(2,3)*zi;
      double Wi = iRotM(3,1)*xi+iRotM(3,2)*yi+iRotM(3,3)*zi;
      double Uo = iRotM(1,1)*xo+iRotM(1,2)*yo+iRotM(1,3)*zo; // WRS
      double Vo = iRotM(2,1)*xo+iRotM(2,2)*yo+iRotM(2,3)*zo;
      double Wo = iRotM(3,1)*xo+iRotM(3,2)*yo+iRotM(3,3)*zo;

      int wireF, wireL;  // First and last wire for this hit
      if ((Ui-wirD_)/wirP_<0) wireF = int( (Ui-wirD_)/wirP_-0.5 );
      else                    wireF = int( (Ui-wirD_)/wirP_+0.5 );
      if ((Uo-wirD_)/wirP_<0) wireL = int( (Uo-wirD_)/wirP_-0.5 );
      else                    wireL = int( (Uo-wirD_)/wirP_+0.5 );
      if (wireL<wireF) {
	  int tmp = wireL; wireL = wireF; wireF = tmp;
      }

	  for (int i = wireF; i<=wireL; i++) {
	    // Look if a digit on this detector with this wire already exists.
	    // Boring, but no idea how to do in other vay... :(
	    list<CsDigit*>::iterator Id; bool found;
	    for (Id = myDigits_.begin(), found = false;
		 Id!=myDigits_.end() && !found; Id++) {
	      if (i==(*Id)->getAddress()) {
		found = true; 	      // Here it is, add this hit to it...
		dynamic_cast<CsMCDigit*>(*Id)->addHit(*(*Ih));
		if (tdc<(*Id)->getDatum()) (*Id)->replaceDatum(tdc);
	      }
	    }
	    if( !found ) { // No digits found with these wire and detector
	      if (i<0 || nWir_<=i) {
		CsErrLog::msg(elAnomaly,__FILE__,__LINE__,
			      "%s: Wire# =%d outside range [0,%d]",GetTBName().c_str(),i,nWir_);
	      }
	      else {
		CsDigit* digit = new CsMCDigit( *this, i , &tdc);
		dynamic_cast<CsMCDigit*>(digit)->addHit( *(*Ih) );
		myDigits_.push_back( digit );   // Add this digit to my list
	      }
	    }
	  }
	}

    } else { // new amplitude MC decoding
      // amplitude MC decoding delivers digits as in the real data: with 3 amplitudes and a wire number
      // also, the implemented charge spread produces multi-strip-clusters which were, though not excluded, very rare in the old MC decoding

      /*
        status:
        - charge sharing for planes with and without intermediate strips implemented
        - 3 strip clusters are created depending on the total a2sum
        - the number of 1 and 2 strip clusters resembles the real data case,
          this depends on the impact point, namely the distance from the next strip

        todo:
        - eta distribution of planes with intermediate strips shows some small features
          not seen in the data. This should have low priority since the it is really small
        - For the clustering only one calibration (hard coded) is used, a connection to
          MySql would be nice. This way also the exact eta correction could by applied
        - The parameters of the functions for clustersize and charge sharing need to
          be determined for all detectors seperately
        - a TCS phase would be nice
        - the 3 strip case needs to be revisited. For once fraction against the clusteramplitude
          drops for very high clusteramps, which is not right mimiced right now. Further the 3
          amplitudeds are just randomly distributed. A pattern in the real data needs to be found
          and implemented
        - The dependence of the clustersize distribution on the main cinderella cut ( tear point sigmas)
          seems to be too small
        - Also clusters of real deltas, meaning size 4 or lager should be implemented
        - Noisy strips would be nice, the problem was the distribution in space
        - Also dead strips would be nice.

        M.L.

      */

      SetupCinderellaClustering();

      list<SIdig> ld; //list of temporary digits (amplitude sets, arising in 1 strip from 1 hit)

      for (Ih = myMCHits_.begin(); Ih!=myMCHits_.end(); Ih++) { // loop on hits
        CsMCTrkHit* thit = dynamic_cast<CsMCTrkHit*>(*Ih);
        if (!thit)
          return;
        if (!thit->getMCTrack()->getParticle()->getCharge() && !thit->getOrigin())
          continue;


        // ***** LOOP ON HITS from charged particles or charged products *****

        double t  = thit->getDTime();  // Delay time (ns)
        double ui = thit->getUin();    // Hit in point (DRS)
        double vi = thit->getVin();
        double wi = thit->getWin();
        double uo = thit->getUout();   // Hit out point (DRS)
        double vo = thit->getVout();
        double wo = thit->getWout();


        // used to "force" amplitude correlation
        // hit of track in same wafer is found and then the elosses are averaged

        list< CsMCHit * > MChits = thit->getMCTrack()->getMCHits();
        list< CsMCHit * >::iterator i_mchit;

        bool amplitude_correlation =true;
        bool found=false;

        if(amplitude_correlation) {
          string name_mate=name;
          if (TString(name).Contains("V")) name_mate.replace(4,1,"U");
          if (TString(name).Contains("U")) name_mate.replace(4,1,"V");
          if (TString(name).Contains("X")) name_mate.replace(4,1,"Y");
          if (TString(name).Contains("Y")) name_mate.replace(4,1,"X");

          for ( i_mchit = MChits.begin(); i_mchit!= MChits.end(); i_mchit++) {
            if ( (*i_mchit)->getDet()->GetTBName() == name_mate ) { found =true; break; }
          }
        }

        double eloss; // dE in detector, keV

        if(amplitude_correlation && found) {
          eloss= ( (thit->getELos()+ dynamic_cast<CsMCTrkHit*>(*i_mchit)->getELos()) )  /2.  *1.E6;
        } else
          eloss = thit->getELos()*1.E6;




        h(name+"eloss", eloss);
        h(name+"mc_time", t);
        t -= // Trigger time (and hence offset) is subtracted from hit times
            triggerOffset;



        //    ***** CHECK TIME GATE *****
        // The time gate is taken from the detector table (i.e. so called
        // "detectors.dat"), contrary to what's done in most other CsDetector's
        // (where the time gate is specified in the options file). And one
        // expects it to be quite loose. This situation mimicks the RD case, where
        // the CsSi class does not apply any time gate.
        if (t<-tGate_/2 ||  tGate_/2<t) continue;
        double tdc = t + tRes_ * CsRandom::gauss();
        //double tdc = t;

        int err;
        HepMatrix iRotM(3,3); iRotM = rotWRS_.inverse( err );
        double xi = rotDRS_(1,1)*ui+rotDRS_(1,2)*vi+rotDRS_(1,3)*wi + xcm_ - _xshift;
        double yi = rotDRS_(2,1)*ui+rotDRS_(2,2)*vi+rotDRS_(2,3)*wi + ycm_ - _yshift;
        double zi = rotDRS_(3,1)*ui+rotDRS_(3,2)*vi+rotDRS_(3,3)*wi + zcm_;
        double xo = rotDRS_(1,1)*uo+rotDRS_(1,2)*vo+rotDRS_(1,3)*wo + xcm_ - _xshift;
        double yo = rotDRS_(2,1)*uo+rotDRS_(2,2)*vo+rotDRS_(2,3)*wo + ycm_ - _yshift;
        double zo = rotDRS_(3,1)*uo+rotDRS_(3,2)*vo+rotDRS_(3,3)*wo + zcm_;

        double Ui = iRotM(1,1)*xi+iRotM(1,2)*yi+iRotM(1,3)*zi; // WRS
        double Vi = iRotM(2,1)*xi+iRotM(2,2)*yi+iRotM(2,3)*zi;
        double Wi = iRotM(3,1)*xi+iRotM(3,2)*yi+iRotM(3,3)*zi;
        double Uo = iRotM(1,1)*xo+iRotM(1,2)*yo+iRotM(1,3)*zo; // WRS
        double Vo = iRotM(2,1)*xo+iRotM(2,2)*yo+iRotM(2,3)*zo;
        double Wo = iRotM(3,1)*xo+iRotM(3,2)*yo+iRotM(3,3)*zo;


        if ( eloss <= 0. ) { continue; }

        const double s2pi=2.5066283; // sqrt(2*pi), for whatever this is needed...

        const double nfactor=eGain_; // eGain_ to be given in option file
        const double sigspace=sWidth_; // sWidth_ to be given in option filee
        const double sig2 = sigspace*sigspace*2.;

        double ltrack = fabs(Uo-Ui); // really ok like this... ??
        double aN;
        int dummyCounter = 0;


        //randomize it to get the charge more realistic
        do{
          aN = eloss*(nfactor+eGSig_*CsRandom::gauss()) ; // eGSig_ to be given in option file
        } while( aN<0. && dummyCounter++<10 );

        if(aN<=0.) aN = eloss*nfactor;                  // need to keep it
        double amp      = aN/sigspace/s2pi;             // amplitude
        double norm     = aN+amp*ltrack;                // new integral
        amp = amp*aN/norm;


        double aUo=0.;
        double aUi=0.;
        double aUc=0.;
        double Utmin=0.;
        double Utmax=0.;

        if(  Uo-Ui >0 ) {      // take into account diffusion, sort aUo ">" aUi
          aUo = Uo+sigspace ;
          aUi = Ui-sigspace ;
          Utmin=Ui;
          Utmax=Uo;
          aUc = (Ui + Uo)/2;
        } else {
          aUo = Ui+sigspace ;
          aUi = Uo-sigspace ;
          Utmin=Uo;
          Utmax=Ui;
          aUc = (Ui + Uo)/2;
        }

        int nwires = (isU || isX) ? 1280 : 1024;
        int wireM=-1; // "Main" wire (closest to the hit)
        int wireN=-1; // "Neighbour" wire (closest on the other side of the hit)

        // define "Main" wire number
        if( (aUc-wirD_) <0) {
          wireM = int( (aUc-wirD_)/wirP_ - 0.5);
        } else {
          wireM = int( (aUc-wirD_)/wirP_ + 0.5);
        }

        h(name+"mc_wire", wireM);

        double fdist = ( aUc-wirD_ - wireM*wirP_ +wirP_/2.)*1000; // distance from main wire in um
        double fdist2 =  ( aUc-wirD_ - wireM*wirP_)*1000;

        double eval = cls( &fdist2 , cls_);

        // get the samples of the normalised amplitude, depending on the time (tdc)
        double amps[3];
        getMCAmp( tdc, amps );

        // determine clustersize. for the moment cls 1,2,3 possible
        int clustersize=-1;
        double amp_tmp = amp*amps[2];


        if (CsRandom::flat()*5 < cls3(&amp_tmp, cls3_)) {
          clustersize=3;
          if (fdist2>0.) {
            wireN=wireM+1;
          } else {
            wireN=wireM-1;
          }
        } else if( CsRandom::flat() > eval ) {
          clustersize=1;
          wireN=-1;
        } else {
          clustersize=2;
          if(fdist2>0.) {
            wireN=wireM+1;
          } else {
            wireN=wireM-1;
          }
        }



        h(name+"mc_dist",fdist);
        if(clustersize==2) h(name+"mc_dist1",fdist);
        if(clustersize==1) h(name+"mc_dist2",fdist);
        h(name+"mc_cls", clustersize);

        if (clustersize>1 && wireN==-1 && wireM) { cout << endl << endl << " wireN    is wrong" << endl << endl;}

        const float acut=6;  // to remove meaningless hits

        // create amplitudes
        if(clustersize==1) {
          if( (amps[0]+amps[1]+amps[2])*amp > acut && 0 <= wireM && wireM < nwires ) {
            SIdig sd;
            sd.wire = wireM;
            sd.amp1 = amp*amps[0] + params_noise[wireM]*CsRandom::gauss(); // 3 amplitude samples, scaled with the wire amplitude
            sd.amp2 = amp*amps[1] + params_noise[wireM]*CsRandom::gauss(); // and smeared to simulate noise behaviour
            sd.amp3 = amp*amps[2] + params_noise[wireM]*CsRandom::gauss();

            h(name+"mc_clusteramp",sd.amp3);

            sd.ref  =  (*Ih); // pointer to CsMCHit
            ld.push_back(sd); // save SI digit for main wire (referred to as "amplitude set" in the next lines)
          }
        } else if(clustersize==2) {
          // charge sharing
          const double feta(get_charge_sharing(fdist));
          const double ampM((fdist2>0) ? (amp*(1-feta)) : (amp*feta)); 
          const double ampN(amp - ampM);
          if( (amps[0]+amps[1]+amps[2])*ampM > acut && 0 <= wireM && wireM < nwires ) {
            SIdig sd;
            sd.wire = wireM;
            sd.amp1 = ampM*amps[0] + params_noise[wireM]*CsRandom::gauss(); // 3 amplitude samples, scaled with the wire amplitude
            sd.amp2 = ampM*amps[1] + params_noise[wireM]*CsRandom::gauss(); // and smeared to simulate noise behaviour
            sd.amp3 = ampM*amps[2] + params_noise[wireM]*CsRandom::gauss();

            h(name+"mc_clusteramp",sd.amp3);

            sd.ref  =  (*Ih); // pointer to CsMCHit
            ld.push_back(sd); // save SI digit for main wire (referred to as "amplitude set" in the next lines)
          }

          if( (amps[0]+amps[1]+amps[2])*ampN > acut && 0 <= wireN && wireN < nwires ) {
            SIdig sd;
            sd.wire = wireN;
            sd.amp1 = ampN*amps[0] + params_noise[wireN]*CsRandom::gauss();
            sd.amp2 = ampN*amps[1] + params_noise[wireN]*CsRandom::gauss();
            sd.amp3 = ampN*amps[2] + params_noise[wireN]*CsRandom::gauss();

            h(name+"mc_clusteramp",sd.amp3);

            sd.ref  =  (*Ih); // should already be set correctly, but you never know...
            ld.push_back(sd); // save SI digit for neighbour
          }
        } else if (clustersize==3) {
          // define next-to-next-to-main wire
          const int wireNN((wireN>wireM) ? (wireN+1) : (wireN-1));

          // calculate amplitudes;
          const double ampM (    amp    * (1./3. * (CsRandom::gauss()*0.1+1)));
          const double ampN ((amp-ampM) * (1./2. * (CsRandom::gauss()*0.1+1)));
          const double ampNN(amp-ampM-ampN);

          // first strip (1/3)
          if( (amps[0]+amps[1]+amps[2])*ampM > acut && 0 <= wireM && wireM < nwires ){
            SIdig sd;
            sd.wire = wireM;
            sd.amp1 = ampM*amps[0] + params_noise[wireM]*CsRandom::gauss(); // 3 amplitude samples, scaled with the wire amplitude
            sd.amp2 = ampM*amps[1] + params_noise[wireM]*CsRandom::gauss(); // and smeared to simulate noise behaviour
            sd.amp3 = ampM*amps[2] + params_noise[wireM]*CsRandom::gauss();

            sd.ref  =  (*Ih); // pointer to CsMCHit
            ld.push_back(sd); // save SI digit for main wire (referred to as "amplitude set" in the next lines)
          }
          // second strip (2/3)
          if( (amps[0]+amps[1]+amps[2])*ampN > acut && 0 <= wireN && wireN < nwires ){
            SIdig sd;
            sd.wire = wireN;
            sd.amp1 = ampN*amps[0] + params_noise[wireN]*CsRandom::gauss(); // 3 amplitude samples, scaled with the wire amplitude
            sd.amp2 = ampN*amps[1] + params_noise[wireN]*CsRandom::gauss(); // and smeared to simulate noise behaviour
            sd.amp3 = ampN*amps[2] + params_noise[wireN]*CsRandom::gauss();

            sd.ref  =  (*Ih); // pointer to CsMCHit
            ld.push_back(sd); // save SI digit for main wire (referred to as "amplitude set" in the next lines)
          }
          // third strip (3/3)
          if( (amps[0]+amps[1]+amps[2])*ampNN > acut && 0 <= wireNN && wireNN < nwires ){
            SIdig sd;
            sd.wire = wireNN;
            sd.amp1 = ampNN*amps[0] + params_noise[wireNN]*CsRandom::gauss(); // 3 amplitude samples, scaled with the wire amplitude
            sd.amp2 = ampNN*amps[1] + params_noise[wireNN]*CsRandom::gauss(); // and smeared to simulate noise behaviour
            sd.amp3 = ampNN*amps[2] + params_noise[wireNN]*CsRandom::gauss();

            sd.ref  =  (*Ih); // pointer to CsMCHit
            ld.push_back(sd); // save SI digit for main wire (referred to as "amplitude set" in the next lines)
          }
        }
      }      //end of MCHits loop


      ld.sort();            // sort by wire#

      // Here we have all what's needed for building digits:
      // amplitude sets as they were given from all single hits,
      // now sum up the signals in identical strips from different hits


      list<SIdig>::iterator id,idnext;
      bool newdig=false;
      vector<CsMCHit*> dighits;
      double ampl[12] = { 0.,0.,0.,0., 0.,0.,0.,0., 0.,0.,0.,0.};
      int iwir;
      for(id = ld.begin(); id != ld.end(); id++) { // loop over amplitude sets

        idnext=id; idnext++;          // look on the previous digit
        if(! newdig ) {               // even though it says "!newdig", this is actually when we DO start a new digit
          dighits.clear();            // reset all information the digit will receive
          iwir=-1;
          ampl[0]=0.; ampl[1]=0.; ampl[2]=0.;
          newdig = true;
        }

        iwir=(*id).wire;
        ampl[0]+=(*id).amp1;          // sum up the 3 amplitudes
        ampl[1]+=(*id).amp2;
        ampl[2]+=(*id).amp3;
        dighits.push_back((*id).ref);


        if(idnext == ld.end()) {   // last of all amplitude sets?
          newdig=false;
        } else if ( (*id).wire != (*idnext).wire ) {   //last amplitude set for this wire?
          newdig=false;
        } else {
          continue;
        }

        if( ! newdig ) {
          //if we added up all single signals from one strip, write the digit:
          //this CsDigit goes into clustering method

          CsDigit* digit = new CsMCDigit(*this, iwir, ampl, 12); //!!!!
          vector<CsMCHit*>::iterator ih;

          for(ih = dighits.begin(); ih != dighits.end(); ih++){
            dynamic_cast<CsMCDigit*>(digit)->addHit( *(*ih) );
          }
          myDigits_.push_back( digit );
        }
      }
    }// end of amplitude MC decoding

  } // end of usage of zebra file


  // if you have not even a zebra file...
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
      if( ip1 == 0 ) ip1 == ip2;
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


// --------------------------------------------------
class SiDigit {
// --------------------------------------------------
  // This class is needed to provide a sorting mechanism
  // that ONLY sorts the wire numbers.
  // CsDigit sorts first chip, than channel. For mapped chips
  // this gives a wrong results

  SiDigit(){}

 public:

  SiDigit(CsDigit*dig) {
    wire=dig->getAddress();
    if(dig->getDataSize()>2) {
      A0=dig->getData()[0];
      A1=dig->getData()[1];
      A2=dig->getData()[2];
    } else {
      A0=0;
      A1=20;
      A2=50;
    }
    d=dig;
  }
  int wire;
  double A0,A1,A2;
  CsDigit *d;
  bool operator<(const SiDigit& dig) {return wire<dig.wire;}
};


//output_function which is called in Thiemos clustering algorithm
// --------------------------------------------------

struct cluster_timing {
  double time;
  double el;
  double eh;
  double position;
  double ep; //not yet right implemented in cinderella, changed later
  double a2sum;
  unsigned int cluster_size;
  unsigned int first_strip;
  cluster_timing *next_cl;
};

void output_function(void *reserved, const double time,
		     const double el, const double eh,
		     const double position, const double ep,
		     const unsigned int cluster_size, const unsigned int first_strip,
             const double a2sum) {

  cluster_timing *new_cl = new cluster_timing;

  cluster_timing *r = (cluster_timing*)reserved;
  while (r->next_cl) r=r->next_cl;
  r->next_cl=new_cl;

  new_cl->time  = time;
  new_cl->el    = el;
  new_cl->eh    = eh;
  new_cl->position = position;
  new_cl->ep = ep; //not yet right implemented in cinderella, changed later
  new_cl->a2sum = a2sum;
  new_cl->cluster_size = cluster_size;
  new_cl->first_strip = first_strip;
  new_cl->next_cl=NULL;

  if (getenv("COND_DB_DEBUG")!=0 || getenv("COND_DB_DEBUG_SI")!=0) {
    cout << setprecision(5) << "output_function: "
         << "time " << time
         << "  el " << el
         << "  eh " << eh
         << "  pos " << position
         << "  ep " << ep
         << "  a2sum " << a2sum
         << "  csize " << cluster_size
         << "  first_strip " << first_strip << endl;
  }
}
// --------------------------------------------------

void CsSiTrackerDetector::clusterize() {
    if (cluster_cind_opt) {
        // in case of MC we can only use cinderalla/mate clustering, if we also use
        // amplitude decoding
        if (CsInit::Instance()->IsAMonteCarloJob())
            if (!amplitudeMCDecoding_)
                CsErrLog::msg(elFatal, __FILE__, __LINE__,
                              "CsSiTrackerDetector::clusterize: %s: can only use cinderella/mate clustering with amplitude simulation.", GetTBName().c_str());	

        // this call also takes care of the mate clustering
        clusterizeCind();
    } else {
        clusterizeRatio();
    }
}

// --------------------------------------------------
void CsSiTrackerDetector::clusterizeCind() {
  // Cinderella clustering.
  // Also used in case of Mate Clustering as basic and as default fallback.
  // Historically:
  // Thiemo's clustering algorithm used in CINDERELLA based on a correct
  // error handling of the calculated time
  if (getenv("COND_DB_DEBUG")!=0 || getenv("COND_DB_DEBUG_SI2")!=0)
    cout << "CsSiTrackerDetector[" << GetTBName() << "]::clusterizeCind()" << endl;

  // the planes assignment of "main" planes and "mate" planes is choosen in a way that
  // the orignial sequenze of clustering is not changed
  // skip U/Y planes here, clustering of those planes is done from their master V/X
  // plane
  if (GetTBName().find('Y')!=string::npos || GetTBName().find('U')!=string::npos)
    return;

  // Find pointer to  other plane on same wafer.
  // Each CsSiTrackerDetector has now a member "mateDet" which is the pointer to the other plane
  // on the same waver. The pointer is NULL by default and initialized only once.
  if (mateDet==NULL) {
    string dummy(GetTBName());

    if (dummy.find('V')!=string::npos) {
      dummy[4] = 'U';
      mateDet  = dynamic_cast<CsSiTrackerDetector*>(FindDetector(dummy));
    }

    if (dummy.find('X')!=string::npos) {
      dummy[4] = 'Y';
      mateDet  = dynamic_cast<CsSiTrackerDetector*>(FindDetector(dummy));
    }
  }
  
  if (mateDet!=NULL)
    cluster_mate=true;

  // prepare clusters and digits
  clearClusterList();

  list<SiDigit>     siDigits = sort_digits(getMyDigits());
  list<SiDigit> mateSiDigits;

  if (cluster_mate) {
    mateDet->clearClusterList();
    mateSiDigits = mateDet->sort_digits(mateDet->getMyDigits());
  }

  // flag if this detector has digits for this event
  bool dig_det(true);
  if (siDigits.empty())
    dig_det = false;
  // flag if mate detector has digits for this event
  if (cluster_mate && mateSiDigits.empty())
    cluster_mate = false;

  // Setup Sergei's alignment patch.  Needs to be called inside clusterize()
  // since the run number is only avaible during processing of the events, not
  // at initialisation stage.
  if (dig_det)
    SetupAlignmentPatch();
  if (cluster_mate)
    mateDet->SetupAlignmentPatch();

  // both planes are empty: return
  if ( (!dig_det) && (!cluster_mate))
    return;

  // errors and covariance matrices
  HepMatrix iRotM(3, 3);
  HepMatrix cov  (3, 3, 0.);
  double    wireDCorr(0.);
  comp_cov(iRotM, cov, wireDCorr);

  HepMatrix mateIRotM(3, 3);
  HepMatrix mateCov  (3, 3, 0.);
  double    mateWireDCorr(0.);
  if (cluster_mate) 
    mateDet->comp_cov(mateIRotM, mateCov, mateWireDCorr);

  // correct for time jumps
  // load calibration data
  if (dig_det)
    get_calib(find_jump());
  if (cluster_mate)
    mateDet->get_calib(mateDet->find_jump());

  // Setup parameter arrays for Cinderella clustering.  Needs to be called
  // inside clusterize() because it needs detector mappings which only are
  // available after the first event has passed through the decoding
  // library.
  if (dig_det)
    SetupCinderellaClustering();
  if (cluster_mate)
    mateDet->SetupCinderellaClustering();

  //clustering parameter
  clusterization_option_t params_clust;

  /* if strip hit times differ by more than tear_point_sigma they are not joined to clusters
   * TODO put in option file, depends on run type
   * higher intensity beam needs tighter cut
   * 2009 primakoff setting: 20*/
  params_clust.tear_point_sigmas = 20; //main important cut
  params_clust.min_cluster_a2sum = 0;
  params_clust.max_error = 10000;     // important, in cinderella 9, now disabled
  params_clust.enable_position = 1;   //swich, that cluster position is calculated

  // do the clustering of the current detector
  // in this vector the clusters are stored, and (depending on the chosen options)
  // passed to mate clustering
  std::list<CsCluster*> tentativeClusters;
  if (dig_det) {
    strip_t* paramsStrip = new strip_t[siDigits.size()];
    fill_strips(siDigits, paramsStrip);

    cluster_timing firstCluster;
    firstCluster.next_cl=NULL;

    silicon_clusterize_plane(paramsStrip, params_strip_count,
                             params_noise, &params_time,
                             &params_clust, &output_function,
                             &firstCluster, NULL);

    // apply eta correction, obtain spatial error...
    finalize_clusters(firstCluster, cov, wireDCorr, tentativeClusters, siDigits);
    delete[] paramsStrip;
  }

  // do the same for the mate detector
  std::list<CsCluster*> mateTentativeClusters;
  if (cluster_mate) {
    strip_t* mateParamsStrip = new strip_t[mateSiDigits.size()];
    mateDet->fill_strips(mateSiDigits, mateParamsStrip);

    cluster_timing mateFirstCluster;
    mateFirstCluster.next_cl=NULL;

    silicon_clusterize_plane(mateParamsStrip, mateDet->params_strip_count,
                             mateDet->params_noise, &(mateDet->params_time),
                             &params_clust, &output_function,
                             &mateFirstCluster, NULL);

    // apply eta correction, obtain spatial error...
    mateDet->finalize_clusters(mateFirstCluster, mateCov, mateWireDCorr, mateTentativeClusters, mateSiDigits);
    delete[] mateParamsStrip;
  }

  // if mate option booked in CORAL's options file and all prerequisits fulfilled, start mate clustering
  // else the output of the Cinderella clustering is used
  if (dig_det && cluster_mate && cluster_mate_opt) {
      this->mate_clustering(tentativeClusters, mateTentativeClusters, wireDCorr);
      mateDet->mate_clustering(mateTentativeClusters, tentativeClusters, mateWireDCorr);
  }

  for (std::list<CsCluster*>::iterator icl=tentativeClusters.begin(); icl!=tentativeClusters.end(); ++icl)
    addCluster(**icl);

  for (std::list<CsCluster*>::iterator icl=mateTentativeClusters.begin(); icl!=mateTentativeClusters.end(); icl++)
    mateDet->addCluster(**icl);

  sortClusters();
  setClusteringDone();

  if (mateDet) {
      mateDet->sortClusters();
      mateDet->setClusteringDone();
  }

  cind_clus_setup=false;
  mateDet->cind_clus_setup=false;
}


// --------------------------------------------------
void CsSiTrackerDetector::clusterizeRatio() {
  if (getenv("COND_DB_DEBUG")!=0 || getenv("COND_DB_DEBUG_SI2")!=0)
    cout << "CsSiTrackerDetector[" << GetTBName() << "]::clusterizeRatio()" << endl;

  clearClusterList();

  const list<CsDigit*>& digits = getMyDigits();

  // if there are no digits for this plane leave this method
  if (digits.empty())
    return;

  // sort digits
  std::list<SiDigit> siDigits = sort_digits(digits);

  // prepare position of first wire and covariance matrix
  HepMatrix iRotM(3,3);
  HepMatrix cov(3,3,0);
  double wireDCorr;
  comp_cov(iRotM, cov, wireDCorr);

  // prepare for clustering
  int lastwire;
  double cog(0.);
  double amp0(0.), amp1(0.), amp2(0.);
  double sTime(0.); // sum of digit times, used only in MC case
  unsigned int nClusters(0);
  unsigned int nClustersAccepted(0);

  std::vector<CsDigit*> clDigits;

  for (std::list<SiDigit>::const_iterator siIter=siDigits.begin(); siIter!=siDigits.end(); ) {
    // In case no real amplitudes are existing, some default values are used
    // (0, 20, 50). Though this info is then not helpful, the following is
    // still correct, even the CoG calculation
    amp0 += siIter->A0;
    amp1 += siIter->A1;
    amp2 += siIter->A2;

    lastwire  = siIter->wire;
    cog      += lastwire*siIter->A2;
    sTime    += siIter->d->getDatum();
    clDigits.push_back(siIter->d);

    // look at next digit
    ++siIter;

    // extra conditions to close the current cluster
    bool extraCond(false);
    if (siIter!=siDigits.end()) {
      if (CsInit::Instance()->IsAMonteCarloJob()) {
        // this cluster is also finished if
        // - The time estimate differs from that of the current cluster. A 3 sigma
        //   cut is taken. Which is approximated by "4*tRes", by making a tradeoff
        //   between a single-hit cluster (w/ 3 sigma=3*sqrt(2)*tRes") and an
        //   already multi-hit cluster.

        // mean time of digits seen so far
        const double cTime = sTime/clDigits.size();

        // check on compatibility
        const double tRes=3.; // 3 ns of resolution built-in (for so rudimentary MC digits don't deserve better)
        if (fabs(siIter->d->getDatum() - cTime) > 4.*tRes)
          extraCond = true;
      } else {
        // this cluster is also finished if
        // - the amplitudes are too small to make a time estimate
        // - the time estimate from the amplitude ratios differ

        if ((siIter->A0+siIter->A1+siIter->A2) >= 8 && (amp0+amp1+amp2) >= 8)
          if (fabs(siIter->A0/siIter->A2-amp0/amp2) >= 0.5 || fabs(siIter->A1/siIter->A2-amp1/amp2) >= 0.4)
            extraCond = true;
      }
    }

    // does it belong to the next cluster?
    // current cluster is finished, when
    // - no futher wires anymore
    // - next wire not neighbour (trivial)
    // - extra condition (see above)
    if (siIter==siDigits.end() || siIter->wire!=lastwire+1 || extraCond) {
      // compute center of gravity
      cog/=amp2;

      // calculate final cluster position
      const double u = wireDCorr+cog*wirP_;
      const double v = 0;
      const double w = zcm_;

      // for real data update the first entry in the covariance matrix
      // depending on the cluster size
      if (CsInit::Instance()->IsADataJob()) {
        // TODO: implement resolution in a calibration file
        //spatial resolution depending on cluster size
        if(clDigits.size()==1)
          cov(1,1)= pow(wirP_*0.95, 2)/12 ; //error [mm^2]
        if(clDigits.size()==2)
          cov(1,1) = pow(wirP_*0.70, 2)/12;
        if(clDigits.size()>=3)
          cov(1,1) = pow(wirP_*0.7*(clDigits.size()-1), 2)/12;
      }

      // create cluster
      CsCluster *cluster=new CsCluster(u,v,w,cov);
      nClusters++;

      // add digits to cluster object
      for(std::vector<CsDigit*>::iterator it=clDigits.begin(); it!=clDigits.end(); ++it)
        cluster->addDigit(**it);

      // add backreference
      cluster->addDet(*this);

      // finalize cluster
      bool accept(true);
      if (CsInit::Instance()->IsAMonteCarloJob()) {
        // mean time of digits
        const double cTime = sTime/clDigits.size();
        cluster->setTime(cTime);
      } else {
        cluster->setAnalog(amp2);
        cluster->addAnalogData(amp0);
        cluster->addAnalogData(amp1);
        cluster->addAnalogData(amp2);

        // decode time
        const double cTime = decodeTime(cluster);

        if (hLevel_>=Normal)
          h(GetTBName()+"_cAmp", amp2);

        if (   ( !fRatioCut || (amp0/amp2<1.2 && amp1/amp2>0.05) )     // cut on the ratios
            && ( !fTimeCut  || (fLCut<cTime && cTime<fRCut)        ) ) { // cut on the calculated time
          if (hLevel_>=Normal)
            h(GetTBName()+"_cAmpCut", amp2);
        } else {
          accept = false;
        }
      }

      if (accept) {
        // cluster accepted
        addCluster(*cluster);
        nClustersAccepted++;

        if (hLevel_>=Normal) { // Fill histograms
          h(GetTBName()+"_cPos", cog);
          h(GetTBName()+"_cSize", clDigits.size());
        }
      } else {
        // cluster rejected
        delete cluster;
      }

      // prepare next cluster
      cog            = 0.;
      amp0=amp1=amp2 = 0.;
      sTime          = 0.;
      clDigits.clear();
    } // end of finishing current cluster
  } // end of loop over digit

  // multiplicity histograms
  if (hLevel_>=Normal) {
    h(GetTBName()+"_cMult",nClusters);
    h(GetTBName()+"_cMultCut",nClustersAccepted);
  }

  sortClusters();
  setClusteringDone();
}

// ------------------------------------------------------------
double CsSiTrackerDetector::decodeTime(CsCluster* cl) {
  // ------------------------------------------------------------
  /*
    function uses calibration values to calculate a time estimate of the
    clusters from the 3 sampled amplitudes
  */

  // Monte Carlo
  // ###########
  //Just mean time of MC hits.
  if (CsInit::Instance()->IsAMonteCarloJob()) {
    double cltime=0; int ntimes=0;
    list<CsDigit*> DigList = cl->getDigitsList();
    list<CsDigit*>::iterator DigIter;
    for(DigIter=DigList.begin();DigIter!=DigList.end();++DigIter) { // loop over digits
      CsMCDigit* d = dynamic_cast<CsMCDigit*>((*DigIter));
      list<CsMCHit*> lh = d->getHits();
      list<CsMCHit*>::iterator ih;
      for(ih = lh.begin(); ih != lh.end(); ih++){
	cltime+= (*ih)->getDTime(); ntimes++;
      }
    }
    cltime/=ntimes;
    double cltime_err=3; //[ns]
    cl->setTime(cltime, cltime_err);
    return cltime;
  }

  // Real Data
  // #########
  const string &name=GetTBName();


  // Unpack time calibrations
  // ========================

  // range of validity of linear region of r0=A0/A2
  const double& r0_min = this->timing_calib[0];
  const double& r0_max = this->timing_calib[1];
  // range of validity of linear region of r2=A2/A3
  const double& r1_min = this->timing_calib[2];
  const double& r1_max = this->timing_calib[3];
  // T0s, taken from the fit of "time= f(r)" at r = 0
  // (r is A0/A2 or A1/A2, time is track time from SciFis)
  // Currently is common for all SI planes. To be tuned to compensate some global timing changes.
  const double& t0_0 = this->timing_calib[4];
  const double& t0_1 = this->timing_calib[5];
  // fitted slopes of linear regions of "time = f(r)" dependence
  const double& sl_0 = this->timing_calib[6];
  const double& sl_1 = this->timing_calib[7];
  // fine time corrections, obtained by "Track time - Calculated Time" distribution fit
  // (as fitted T0s are not precise enought).
  const double& tc_0 = this->timing_calib[8];
  const double& tc_1 = this->timing_calib[9];

  const double TCS_T0 = 40.0;     // Calibration constant (Maximum of TCS distr hist.)
  const double time_resol = 3.0; // for weighted mean calculation. (has to be also subject for calibrations)


  // Get TCS Phase
  // =============

  CsEvent* event = CsEvent::Instance();
  tcs_cor = event->getTCSPhaseTime()-TCS_T0;
  if(hLevel_==High)
    h(name+"_tcsPhase",tcs_cor);


  // Calculate amplitude ratios
  // ==========================

  vector<double> analog = cl->getAllAnalogData();
  amp0 = analog[1];
  amp1 = analog[2];
  amp2 = analog[3];

  r0 = amp0/amp2;
  r1 = amp1/amp2;

  hit_time_0 =  1000;  hit_time_err_0 = 1.E+10;
  hit_time_1 =  1000;  hit_time_err_1 = 1.E+10;
  hit_time   =  1000;  hit_time_err   = 1.E+10;

  if(r0_min < r0 && r0 < r0_max)  {
    hit_time_0 = tcs_cor+t0_0 + sl_0*r0 - tc_0;
    hit_time_err_0 = time_resol;
  }
  if(r1_min < r1 && r1 < r1_max)  {
    hit_time_1 = tcs_cor+t0_1 + sl_1*r1 - tc_1;
    hit_time_err_1 = time_resol;
  }

  double w0= 1./(hit_time_err_0* hit_time_err_0);
  double w1= 1./(hit_time_err_1* hit_time_err_1);
  hit_time = (hit_time_0*w0 + hit_time_1*w1)/(w0+w1);
  hit_time_err = sqrt(1./(w0+w1));

  if(hLevel_>=Normal) {
    h(name+"_cTime",hit_time);
    if(hLevel_==High){
      h(name+"_r0",tcs_cor,r0);
      h(name+"_r1",tcs_cor,r1);
      h(name+"_time0",hit_time_0);
      h(name+"_time1",hit_time_1);

    }
  }
  // Set time for Silicon cluster
  cl->setTime(hit_time, hit_time_err);

  return hit_time;

}

// ------------------------------------------------------------
void CsSiTrackerDetector::comp_cov(HepMatrix& iRotM, HepMatrix& cov, double & wireDCorr) {
    // Corrections and Errors
    // ======================

    // compute alignment correction and add 1st wire position
    //!!!OBSOLETE _deltax _deltay == 0 for data!!
    //but neccessary to simulate detectors misalignment with MC data


    if (iRotM.num_col() != 3 || iRotM.num_row() != 3)
        CsErrLog::msg(elFatal, __FILE__, __LINE__,
                      "CsSiTrackerDetector::comp_cov: %s: wrong Matrix Dimension in iRotM.", GetTBName().c_str());	
	if (cov.num_col() != 3 || cov.num_row() != 3)
        CsErrLog::msg(elFatal, __FILE__, __LINE__,
                      "CsSiTrackerDetector::comp_cov: %s: wrong Matrix Dimension in cov.", GetTBName().c_str());	

    int err;
    iRotM = rotWRS_.inverse(err);
    if (err!=0)
        CsErrLog::msg(elError, __FILE__, __LINE__,
                      "CsSiTrackerDetector::comp_cov: %s: error while inverting rotWRS_.", GetTBName().c_str());	

    wireDCorr = wirD_ + iRotM(1,1) * _deltax + iRotM(1,2) *_deltay; 			
    wireDCorr_ = wirD_ + iRotM(1,1) * _deltax + iRotM(1,2) *_deltay;
			
    cov(1,1) = pow(wirP_,2)/12;                      // if cluster error is unknown
    cov(2,2) = pow(max(getXsiz(), getYsiz()),2);     // just very big error
    cov(3,3) = 10.;

    return;
}


// ------------------------------------------------------------
int CsSiTrackerDetector::find_jump() {
// ------------------------------------------------------------
  // Check for timejumps during this spill
  // =====================================
  // Timejumps are expressed in multiples of half clockcycles.
  // This value is saved in tjump, first argument of function

  int tjump(0);

  if (time_12nsjumps){
    unsigned int runnb = CsEvent::Instance()->getRunNumber();
    unsigned int spillnb = CsEvent::Instance()->getBurstNumber();

    for (int i=(timing_spill_jumps.size()/3); i>0; i--){
      if ((int)runnb > timing_spill_jumps[3*i-3] || ((int)runnb == timing_spill_jumps[3*i-3] && (int)spillnb >= timing_spill_jumps[3*i-2])) {
        tjump = timing_spill_jumps[3*i-1];
        break;
      }
    }
    CsErrLog::msg(elInfo,__FILE__,__LINE__,"%s TimeJump value for spill #%u: %i",
                  this->GetTBName().data(),spillnb, tjump);
  }

  return tjump;
}

// ------------------------------------------------------------
void CsSiTrackerDetector::get_calib(int tjump) {
// ------------------------------------------------------------
      // Unpack time calibrations for Jan's approximation function for ratio 0 and ratio 1
      // =================================================================================
      // TODO check if detour with 'const double &' is necessary


      // hard coded values of time calibration and of eta function for the MC case
      if(CsInit::Instance()->IsAMonteCarloJob()){
		
        p_ratio1.t0 = 25.8;
        p_ratio1.r0 = 1.649;
        p_ratio1.a  = 0.05034;
        p_ratio1.c  = 0.005775;
        p_ratio1.b  = 64.68;
        p_ratio1.d  = -0.06434;
        p_ratio1.precal = 0;
		p_ratio1.width  = 0;

        p_ratio2.t0 = 0.1;
        p_ratio2.r0 = 1.227;
        p_ratio2.a  = 0.052;
        p_ratio2.c  = 0.0024;
        p_ratio2.b  = 135.;
        p_ratio2.d  = 0.3728;
        p_ratio2.precal = 0;
		p_ratio2.width  = 0;
	  	
	    if(isU || isX )
	    {
	      p_eta.eta0 = (18.1812);
          p_eta.eta1 = (-342.077);
          p_eta.eta2 = (1762.25);
          p_eta.eta3 = (-3997.68);
          p_eta.eta4 = (4248.65);
          p_eta.eta5 = (-1707.42);
	    }
	    else if (isV || isY )
	    {	 	
	      p_eta.eta0 = 5.52724;
	      p_eta.eta1 = -141.362;
	      p_eta.eta2 = 914.437;
	      p_eta.eta3 = -2376.42;
	      p_eta.eta4 = 2700.69;
	      p_eta.eta5 = -1110.55;
	    }	
		
	  }	  // end if MC job

      else{
		
      shft0= this->timing_calib[8];  // timing shift depending on cluster size
      shft1= this->timing_calib[9];  // timing shift depending on cluster size

      const double& jr0  = this->timing_calib[12];
      const double& jt0  = this->timing_calib[13];
      const double& ja   = this->timing_calib[14];
      const double& jb   = this->timing_calib[16]; //order changed like in CINDERELLA
      const double& jc   = this->timing_calib[15]; //order changed like in CINDERELLA
      const double& jd   = this->timing_calib[17];

      const double& j2r0 = this->timing_calib[20];
      const double& j2t0 = this->timing_calib[21];
      const double& j2a  = this->timing_calib[22];
      const double& j2b  = this->timing_calib[24]; //order changed like in CINDERELLA
      const double& j2c  = this->timing_calib[23]; //order changed like in CINDERELLA
      const double& j2d  = this->timing_calib[25];

      const double& eta0 = this->eta_corr[0];
      const double& eta1 = this->eta_corr[1];
      const double& eta2 = this->eta_corr[2];
      const double& eta3 = this->eta_corr[3];
      const double& eta4 = this->eta_corr[4];
      const double& eta5 = this->eta_corr[5];

      p_ratio1.r0     = jr0;
      p_ratio1.t0     = jt0;
      p_ratio1.a      = ja;
      p_ratio1.b      = jb;
      p_ratio1.c      = jc;
      p_ratio1.d      = jd;
      p_ratio1.precal = 0;
      p_ratio1.width  = 0;

      p_ratio2.r0     = j2r0;
      p_ratio2.t0     = j2t0;
      p_ratio2.a      = j2a;
      p_ratio2.b      = j2b;
      p_ratio2.c      = j2c;
      p_ratio2.d      = j2d;
      p_ratio2.precal = 0;
      p_ratio2.width  = 0;

      p_eta.eta0 =eta0;
      p_eta.eta1 =eta1;
      p_eta.eta2 =eta2;
      p_eta.eta3 =eta3;
      p_eta.eta4 =eta4;
      p_eta.eta5 =eta5;

      p_ratio1.t0  += tjump * 12.86;  //just add jump values to t0's
	  p_ratio2.t0 += tjump * 12.86;
	
	  }

	
	
	
}

// ------------------------------------------------------------
std::list<SiDigit> CsSiTrackerDetector::sort_digits(const std::list<CsDigit*>& digits) {
    // Fetch and sort all digits
    // =========================
    // CsDigits are copied from 'digits' to 'SiDigList', which is a list of SiDigit's
    // These SiDigit's are sorted w. r. t. the wire number and SiDigit is the
    // data type which is used in the following

    std::list<SiDigit> siDigits;

    int counter(0);
    // loop over all digits and fill them into siDigits
    for(list<CsDigit*>::const_iterator Id=digits.begin();Id!=digits.end();++Id) {
        siDigits.push_back(SiDigit(*Id));
        counter++;
    }

    // sort by wire number
    siDigits.sort();

    //multiplicity cut for detector analysis mode
    if (counter>fMultCut){
        CsErrLog::msg(elInfo, __FILE__, __LINE__,
                      "\"%s\" multiplicity cut on event %d",
                      GetTBName().c_str(), CsEvent::Instance()->getEventNumberInRun());
        siDigits.clear();
        return siDigits;
    }

    // record number of strips
    if (hLevel_>=Normal){
        h(GetTBName()+"_strips", siDigits.size());
    }

    params_strip_count = siDigits.size();

    return siDigits;
}

// ------------------------------------------------------------
inline void CsSiTrackerDetector::fill_strips(const std::list<SiDigit>& siDigits, strip_t* params_strip) {
    // SiDigits are converted to strip_t
    // ===================================
    // this step is needed, to make use of the cinderella functions

    size_t index = 0;
    // main important loop in this funtion, rest is for debug purpose...
    for(list<SiDigit>::const_iterator siIter=siDigits.begin(); siIter!=siDigits.end(); ++siIter, ++index) {
        // variables for struct strip_t for Thiemo's timing function
        // division of amplitude for comparison with CINDERELLA, there
        // the digits have lowest-bit truncated ADC units
        params_strip[index].wire =  siIter->wire;
        params_strip[index].a0f  =  (siIter->A0)/2;
        params_strip[index].a1f  =  (siIter->A1)/2;
        params_strip[index].a2f  =  (siIter->A2)/2;

        if (getenv("COND_DB_DEBUG")!=0 || getenv("COND_DB_DEBUG_SI2")!=0) {
            cout << setprecision(7)<<"$$$$$$$$ " << index << "  params_strip.wire        "
                 <<  params_strip[index].wire << endl;
            cout << setprecision(7)<<"$$$$$$$$ " << index << "  params_strip.a0f         "
                 <<  params_strip[index].a0f  << endl;
            cout << setprecision(7)<<"$$$$$$$$ " << index << "  params_strip.a1f         "
                 <<  params_strip[index].a1f  << endl;
            cout << setprecision(7)<<"$$$$$$$$ " << index << "  params_strip.a2f         "
                 <<  params_strip[index].a2f  << endl;
        }
    }

    if (getenv("COND_DB_DEBUG")!=0 || getenv("COND_DB_DEBUG_SI2")!=0) {
        for (int i=0; i < 1300; i++)
            cout << "noise[" << i << "]: " << params_noise[i] << endl;
    }

    //Thiemo's timing algorithm based on Jan's time function. Correct handling and calculation
    //of the errors on the silicon time
    silicon_calculate_timings( params_strip, params_strip_count, params_noise, &params_time, NULL);

    if (getenv("COND_DB_DEBUG")!=0 || getenv("COND_DB_DEBUG_SI2")!=0) {
        for(unsigned int i=0; i < params_strip_count; i++) {
            cout << setprecision(7)<<"CsSiTrackerDetector::clusterize:calc_timing "
                 << i << "  params_strip_count  " <<  params_strip_count    << endl;
            cout << setprecision(7)<<"CsSiTrackerDetector::clusterize:calc_timing "
                 << i << "  params_strip.wire   " <<  params_strip[i].wire  << endl;
            cout << setprecision(7)<<"CsSiTrackerDetector::clusterize:calc_timing "
                 << i << "  params_strip.a0f    " <<  params_strip[i].a0f   << endl;
            cout << setprecision(7)<<"CsSiTrackerDetector::clusterize:calc_timing "
                 << i << "  params_strip.a1f    " <<  params_strip[i].a1f   << endl;
            cout << setprecision(7)<<"CsSiTrackerDetector::clusterize:calc_timing "
                 << i << "  params_strip.a2f    " <<  params_strip[i].a2f   << endl;
            cout << setprecision(7)<<"CsSiTrackerDetector::clusterize:calc_timing "
                 << i << "  params_strip.time   " <<  params_strip[i].time  << endl;
            cout << setprecision(7)<<"CsSiTrackerDetector::clusterize:calc_timing "
                 << i << "  params_strip.el     " <<  params_strip[i].el    << endl;
            cout << setprecision(7)<<"CsSiTrackerDetector::clusterize:calc_timing "
                 << i << "  params_strip.eh     " <<  params_strip[i].eh    << endl;
            cout << setprecision(7)<<"CsSiTrackerDetector::clusterize:calc_timing "
                 << i << "  params_strip.time0  " <<  params_strip[i].time0 << endl;
            cout << setprecision(7)<<"CsSiTrackerDetector::clusterize:calc_timing "
                 << i << "  params_strip.time1  " <<  params_strip[i].time1 << endl;
            cout << setprecision(7)<<"CsSiTrackerDetector::clusterize:calc_timing "
                 << i << "  params_strip.el0    " <<  params_strip[i].el0   << endl;
            cout << setprecision(7)<<"CsSiTrackerDetector::clusterize:calc_timing "
                 << i << "  params_strip.eh0    " <<  params_strip[i].eh0   << endl;
            cout << setprecision(7)<<"CsSiTrackerDetector::clusterize:calc_timing "
                 << i << "  params_strip.el1    " <<  params_strip[i].el1   << endl;
            cout << setprecision(7)<<"CsSiTrackerDetector::clusterize:calc_timing "
                 << i << "  params_strip.eh1    " <<  params_strip[i].eh1   << endl;
            cout << setprecision(7)<<"CsSiTrackerDetector::clusterize:calc_timing "
                 << i << "  params_strip.consistency " <<  params_strip[i].consistency << endl;

            const string &name=GetTBName();

            if(isfinite(params_strip[i].time0)){
                cout << name  << " wire " <<  params_strip[i].wire
                     << " t0 "  << params_strip[i].time0
                     << " el0 " << params_strip[i].el0
                     << " eh0 " << params_strip[i].eh0 << endl;
            }
            if (hLevel_>=High) {
                time_deb_t0 = params_strip[i].time0;
                h(name+"_time_deb_t0", time_deb_t0);
                time_deb_t1 = params_strip[i].time1;
                h(name+"_time_deb_t1", time_deb_t1);
                time_deb = params_strip[i].time;
                h(name+"_time_deb", time_deb);
            }
        }
    }
}


// ------------------------------------------------------------
void CsSiTrackerDetector::finalize_clusters(cluster_timing& first_cluster, HepMatrix& cov, double wireDCorr, std::list<CsCluster*> &vtentative_clusters, const std::list<SiDigit>& siDigits) {
// ------------------------------------------------------------
	  // Clusters are finalized
      // =========================
	  // This includes following main steps
	  // - add time errors to clusters
	  // - add pedestals and noise values to clusters
	  // - add spatial error depending on clustersize
	  // - apply etacorrection
	  // - finally: add cluster to CsDetector


      const string &name=GetTBName();
      i_cl=0; //counter for number of clusters per event
      cluster_timing *temp;
      for (cluster_timing *next_cluster = first_cluster.next_cl;
	  next_cluster;
	  temp = next_cluster, next_cluster = temp->next_cl, delete temp, i_cl++) {

			if (getenv("COND_DB_DEBUG")!=0 || getenv("COND_DB_DEBUG_SI2")!=0) {
			cout << " cluster "  << i_cl
				<< ": time = "  << next_cluster->time
				<< ", pos = "   << next_cluster->position
				<< ", cluster_size = "  << next_cluster->cluster_size
				<< ", first_strip =  " <<  next_cluster->first_strip
				<< endl;
			}

	cluster_size_out = next_cluster->cluster_size;
	first_wire_out = next_cluster->first_strip;
	
			if (hLevel_>=High) {
			time_out = next_cluster->time;
			h(name+"_time_out", time_out);
			pos_out = next_cluster->position;
			h(name+"_pos_out", pos_out);
			h(name+"_cluster_size_out", cluster_size_out);
			}
      	 // h(name+"_cSize",cluster_size_out);

	////calculate position in coordinates [mm]
	double u=wireDCorr_+((double&)next_cluster->position)*wirP_ ;
	
	double v=0;
	
	// error 666 is used as flag, to identify "misidentified 3 strip clusters"
	// to identify this clusters also on phast level, a v coordinate unequal 0 is written
	if(next_cluster->ep == 666) {  v=-1.0;     h(name+"3strip",1);}
	double w=zcm_;


			if (getenv("COND_DB_DEBUG")!=0 || getenv("COND_DB_DEBUG_SI2")!=0) {
			cout << " first_wire_out " << first_wire_out
				<< " cluster_size_out " <<  cluster_size_out
				<< " sum " << cluster_size_out+first_wire_out << endl;
				
			}
				

	// -------------------------------------------------
	// alignment corrections for hadron runs 42320-43350
	// -------------------------------------------------

	// unit of resco in mm!
	//  cout <<  " alignpatch " <<align_patch << "  u before " << u << " ";

	u = u - rescor;

	// cout << GetTBName()<< " " << u << " + " << rescor << " = " << u+rescor<< endl;

	//end of alignment correction
	// ----------------------------------------------------------

	unsigned int cl_count = cluster_size_out+first_wire_out;
	unsigned int w_out    = next_cluster->first_strip;

	double amp2_out_first = 0;
	amp2_out=0;
	amp1_out=0;
	amp0_out=0;
	r0_out=0;
	r1_out=0;

	list<SiDigit>::const_iterator SiIter;
	int maxamp=0;

    std::list<CsDigit*> DigList;

	// DigList fill
	for(SiIter=siDigits.begin(); SiIter!=siDigits.end(); ++SiIter ) {
	  if((unsigned int)SiIter->wire==w_out){
	  // TODO: check if this is the right place???
	  DigList.clear();


	    for(unsigned int iw=w_out; iw< cl_count; iw++){
	
	      if (iw != (unsigned int) SiIter->wire) {
		  CsErrLog::msg(elInfo,__FILE__,__LINE__," Alarm!!! Inconsistency in SI clustering!" );
	      }
	      SiIter->d->getData()[3] = next_cluster->el;
	      SiIter->d->getData()[4] = next_cluster->eh;
	      SiIter->d->getData()[5] = next_cluster->time;
	      SiIter->d->getData()[9] = params_pedestal[iw];
	      SiIter->d->getData()[10] = params_noise[iw];
	      SiIter->d->getData()[11] = pedestal_mean_chip[iw/128];
		
	   	  DigList.push_back(SiIter->d);
	   	
	   	  if(SiIter->A2 > maxamp) maxamp= int(SiIter->A2);
	   	
	      amp2_out+=SiIter->A2;
	      amp1_out+=SiIter->A1;
	      amp0_out+=SiIter->A0;
	
				/* cout << SiIter->d->getData()[0] << " "
				<< SiIter->d->getData()[1] << " "
				<< SiIter->d->getData()[2] << " "
				<< SiIter->d->getData()[6] << " "
				<< SiIter->d->getData()[7] << " "
				<< SiIter->d->getData()[8] << " "
				<< SiIter->d->getData()[9] << " "
				<< endl;*/

	      if(iw == w_out)amp2_out_first=SiIter->A2;

	      if (iw<cl_count-1) SiIter++;
	    }
	  }
	}

	// Spatial erros are written in covariance matrix, depending on clustersize and detector
	get_spatialerror(cov);

	
	// apply eta correction to cluster position
	apply_etacorr(u, amp2_out_first/amp2_out);
		
			if (getenv("COND_DB_DEBUG")!=0 || getenv("COND_DB_DEBUG_SI2")!=0) {
			cout << setprecision(5) << " U " << u
				<< " cov(1,1) " <<  cov(1,1)
				<< " cluster_size_out " <<  cluster_size_out << endl;
			}
    h(name+"_finalizePos",u);

	CsCluster* cluster_thiemo = new CsCluster(u,v,w,cov);
      	
	//cluster ratio calculating
	if(amp2_out>0){
	  r0_out=amp0_out/amp2_out;
	  r1_out=amp1_out/amp2_out;
	}


    list<CsDigit*>::const_iterator DigIter;
	// add digits to cluster object
	for(DigIter=DigList.begin();DigIter!=DigList.end();++DigIter) {
	  cluster_thiemo->addDigit(**DigIter);
	}
	// add backreference
	cluster_thiemo->addDet(*this);

	// set time for Silicon cluster
	double clustertime = next_cluster->time;
	double clustertime_err = (next_cluster->eh+next_cluster->el)/2;

	clustertime_err_out = (next_cluster->eh+next_cluster->el)/2;
	cluster_thiemo->setTime(clustertime,clustertime_err);
	cluster_thiemo->setTimeError(clustertime_err);
	cluster_thiemo->setAnalog(amp2_out);
	cluster_thiemo->addAnalogData(amp2_out);
	cluster_thiemo->addAnalogData(amp1_out);
	cluster_thiemo->addAnalogData(amp0_out);
	cluster_thiemo->addAnalogData(next_cluster->eh);
	cluster_thiemo->addAnalogData(next_cluster->el);

			if (getenv("COND_DB_DEBUG")!=0 || getenv("COND_DB_DEBUG_SI2")!=0) {
			cout << " u " << u << " v " << v << " w " << w << " clustertime_err " << clustertime_err << endl;
				
			CsCluster *c = cluster_thiemo;
			cout << setprecision(5) << "u: " << c->getU() << "   v: " << c->getV()
				<< "   t: " << clustertime << "   dt: " << clustertime_err << endl;
			}

	// store cluster
    vtentative_clusters.push_back(cluster_thiemo);



	if (hLevel_>=High) {
	  h(name+"_cMultCut_thiemo",i_cl);
	  T_SIdeb_Ccl->Fill();
	}
      }


} // end finalize_clusters




// ------------------------------------------------------------
inline void CsSiTrackerDetector::get_spatialerror(HepMatrix& cov) {
// ------------------------------------------------------------
	//!!TODO implement resolution in a calibration file!!
	// spatial resolution depending on cluster size
	// also dependence on amplitude!!
	// also dependence on wafer
	// IMPORTANT: Unit of error is MM !!!

	  const string &name=GetTBName();

	//different resolutions for XU/ and VY planes
	if (isX || isU) {
	  if(cluster_size_out==1) {
			if ( name.find("SI01X")==0 ) cov(1,1)= 1.1e-04; // 10.5 res
			if ( name.find("SI01U")==0 ) cov(1,1)= 8.8e-05; // 9.4 res
		
			if ( name.find("SI02X")==0 ) cov(1,1)= 1.2e-04; // since SI02X was dead since 2009, a high default value is given here
			if ( name.find("SI02U")==0 ) cov(1,1)= 1.0e-04; // since SI02X was dead since 2009, no real resolution for SI02U could determined. the value is extrapolated from other U detectors in front of target
		
			if ( name.find("SI03X")==0 ) cov(1,1)= 1.1e-04; // 10.6 res
			if ( name.find("SI03U")==0 ) cov(1,1)= 9.5e-05; // 9.8 res
		
			if ( name.find("SI04X")==0 ) cov(1,1)= 4.5e-05; // 6.7 res
			if ( name.find("SI04U")==0 ) cov(1,1)= 4.7e-05; // 6.9 res
		
			if ( name.find("SI05X")==0 ) cov(1,1)= 3.9e-05; // 6.2 res
			if ( name.find("SI05U")==0 ) cov(1,1)= 4.9e-05; // 7.0 res
	  }
	
	  if(cluster_size_out==2){
			if ( name.find("SI01X")==0 ) cov(1,1)= 2.7e-05; // res 5.2 um
			if ( name.find("SI01U")==0 ) cov(1,1)= 2.7e-05; // res 5.2 um
		
			if ( name.find("SI02X")==0 ) cov(1,1)= 3.6e-05; // res 6.0 um see comment above
			if ( name.find("SI02U")==0 ) cov(1,1)= 3.6e-05; // res 6.0 um see comment above
		
			if ( name.find("SI03X")==0 ) cov(1,1)= 1.6e-05; // res 4.0 um
			if ( name.find("SI03U")==0 ) cov(1,1)= 1.6e-05; // res 4.0 um
		
			if ( name.find("SI04X")==0 ) cov(1,1)= 9.8e-06; // res 3.1 um
			if ( name.find("SI04U")==0 ) cov(1,1)= 1.0e-05; // res 3.2 um
		
			if ( name.find("SI05X")==0 ) cov(1,1)= 1.0e-05; // res 3.2um
			if ( name.find("SI05U")==0 ) cov(1,1)= 1.6e-05; // res 4.0um  // to be UPDATED
	  }
}


if (isY || isV) {
	  if(cluster_size_out==1) {
			if ( name.find("SI01Y")==0 ) cov(1,1)= 4.8e-05; // res 6.9 um
			if ( name.find("SI01V")==0 ) cov(1,1)= 4.9e-05; // res 7.0 um
		
			if ( name.find("SI02Y")==0 ) cov(1,1)= 2.8e-05; // res 5.3 um
			if ( name.find("SI02V")==0 ) cov(1,1)= 2.9e-05; // res 5.4 um
		
			if ( name.find("SI03Y")==0 ) cov(1,1)= 4.8e-05; //  res 6.9 um
			if ( name.find("SI03V")==0 ) cov(1,1)= 4.8e-05; //  res 6.9 um
		
			if ( name.find("SI04Y")==0 ) cov(1,1)= 5.8e-05; //  res 7.6 um
			if ( name.find("SI04V")==0 ) cov(1,1)= 6.1e-05; //  res 7.8 um
		
			if ( name.find("SI05Y")==0 ) cov(1,1)= 4.8e-05; //  res 6.9 um
			if ( name.find("SI05V")==0 ) cov(1,1)= 4.9e-05; //  res 7.0 um // to be UPDATED
 	}

	  if(cluster_size_out==2){
      //error IN MM !!!!!
			if ( name.find("SI01Y")==0 ) cov(1,1)= 8.4e-06; // res 2.9 um
			if ( name.find("SI01V")==0 ) cov(1,1)= 9.6e-06; // res 3.1 um
		
			if ( name.find("SI02Y")==0 ) cov(1,1)= 9.6e-06; // res 3.1 um
			if ( name.find("SI02V")==0 ) cov(1,1)= 9.6e-06; // res 3.1 um
		
			if ( name.find("SI03Y")==0 ) cov(1,1)= 9.0e-06; //  res 3.0 um
			if ( name.find("SI03V")==0 ) cov(1,1)= 8.4e-06; //  res 2.9 um
		
			if ( name.find("SI04Y")==0 ) cov(1,1)= 1.3e-05; //  res 3.6 um
			if ( name.find("SI04V")==0 ) cov(1,1)= 1.3e-05; //  res 3.6 um
		
			if ( name.find("SI05Y")==0 ) cov(1,1)= 1.3e-05; //  res 3.6 um
			if ( name.find("SI05V")==0 ) cov(1,1)= 1.6e-05; //  res 4.0 um
	  }
	}

    // since clusters with clusersize >2 arise in principal from 2 main cases:
    // - delta electrons
    // - two tracks overlapping in a way, that clustering algorythm cannot distinguish
    // both cases are not handled at the moment, so the error estimate is both large
    // and justified. For residuals of multistrip clusters see phd thesis of A.M. Dinkelbach

	if(cluster_size_out==3) {  cov(1,1) = 6.25e-04; }
	if(cluster_size_out==4) {  cov(1,1) = 2.50e-03; }
	if(cluster_size_out>4)	{  cov(1,1) = 5.63e-03; }

return;
} // end get_spatialerror




//------------------------------------------------------------
inline void CsSiTrackerDetector::apply_etacorr(double& u, double r) {
//------------------------------------------------------------

	// eta correction is applied to cluster position
	// u is the position which is corrected
	// r is the amplitude/charge ratio between 1. (left) strip and total amplitude/charge of cluster

if(cluster_size_out==2){

	if (isX || isU) {
	  //  double r= amp2_out_first/amp2_out;
	    double fr=-(
			((((p_eta.eta5*r + p_eta.eta4) * r + p_eta.eta3) * r + p_eta.eta2) * r + p_eta.eta1) * r + p_eta.eta0
			) * 1e-3;
	    u -= fr;
	}
	
	if (isY || isV) {
    //    double r= amp2_out_first/amp2_out;
	    double fr =-(
			 ((((p_eta.eta5*r + p_eta.eta4) * r + p_eta.eta3) * r + p_eta.eta2) * r + p_eta.eta1) * r + p_eta.eta0
			 ) * 1e-3;
         u -= fr;
	  }
	}

	else {
		// clustersize != 2 , i. e. no eta correction is applied
		return;
		}

}


//---------------------------------------------------------------------------------------
void CsSiTrackerDetector::mate_clustering(std::list<CsCluster*> & vtent_cl, std::list<CsCluster*> & vtent_cl_mate,  double wireDCorr){
//---------------------------------------------------------------------------------------
	// Main function of mate clustering
	// ================================================================================
	// search in clusters of Cinderella output for clusters with size >2
	// For these clusters splitting into two is attempted.
	// If inside find_best two suitable clusters on the same waver
	// are found, this is information is used to split the clusters.
	// Ohterwise 3 strip clusters are split on the middle strip with
	// a factor of 0.5. For clusters of size 4 or lager two new
	// clusters of size 2 are created from the old first two and last two
	// strips.

	const string & name = GetTBName();

	// look for clusters with more than 2 digits
	
	std::list<list<CsCluster*>::iterator> Iwide;		// Indices of clusters with cls > 2
	std::list<list<CsCluster*>::iterator> Iwide_m;
	
	// threshold for a2sum of clusters, for smaller clusters there is no splitting attempt
    const int a2threshold = 20;

    // loop over cluster list, for clusters where splitting should be tried
	int c=0;
	for( list<CsCluster*>::iterator itc = vtent_cl.begin(); itc != vtent_cl.end(); itc++, c++)
		{
		double analog=0;
		const std::list<CsDigit*> & diglist = (*itc)->getDigitsList();
		
		if( (*itc)->getDigitsList().size()>2 && (*itc)->getAnalog(analog) && analog > a2threshold) { Iwide.push_back(itc); }
		}


	c=0;
	for( list<CsCluster*>::iterator itc_m=vtent_cl_mate.begin(); itc_m != vtent_cl_mate.end(); itc_m++, c++)
		{
		double analog=0;
		const std::list<CsDigit*> & diglist = (*itc_m)->getDigitsList();
			
		if( (*itc_m)->getDigitsList().size()>2 && (*itc_m)->getAnalog(analog) && analog > a2threshold) { Iwide_m.push_back(itc_m); }
		}

	// nothing to be done, no cluster fullfilled criteria for being split up
	if(!Iwide.size()) {	return; }
	

	// requirements should be fullfilled, get to start!
	
	for(list<list<CsCluster*>::iterator>::iterator Iiter = Iwide.begin(); Iiter !=Iwide.end(); Iiter++)
	{	int a2sum1, a2sum2;
	
		// fill controll histo
	    h(name+"findbestcalled",1);
	
 		if( find_best_candidates(vtent_cl_mate, Iiter, a2sum1, a2sum2) )
			{
			split_cluster(vtent_cl, *Iiter, a2sum1, a2sum2, wireDCorr);
			// fill histo whenever mate splitting was successfull
			h(name+"findbestsuccess",1);
			}
																		
		else
			{
			// If mate splitting failed and clustersize is smaller than 4,
			// i.e. cls==3 then split cluster into two. Middle strip is
			// halfed	
			if( (*(*Iiter))->getDigitsList().size()<4)
				  {
				  split_cluster(vtent_cl, *Iiter, wireDCorr);
				  }
				
			// Mate splitting failed and clustersize is 4 or lager
			// This clusters are supposed to stem from delta electrons
			// Two new Clusters are generated from the first and the
			// last two strips of the old clusters
			else
			      {
				  split_delta(  vtent_cl, *Iiter, wireDCorr);
				  }
			}
		
	
	}

	   if (getenv("COND_DB_DEBUG")!=0 || getenv("COND_DB_DEBUG_SI2")!=0  || getenv("COND_DB_DEBUG_SI_MATE")!=0   )
			{
			cout << " at end of function mate clustering " << endl;
			}
}

//------------------------------------------------------------
void CsSiTrackerDetector::split_cluster(std::list<CsCluster*> & vtent_cl, std::list<CsCluster*>::iterator Iiter, int a2sum1, int a2sum2, double wireDCorr) {
//------------------------------------------------------------
//	cluster specified by Iiter is split up in relation to a2sum1 and a2sum2
//	2 new clusters are created and stored in the cluster list
//	old cluster is deleted
	
	if (getenv("COND_DB_DEBUG")!=0 || getenv("COND_DB_DEBUG_SI2")!=0 || getenv("COND_DB_DEBUG_SI_MATE")!=0) {
	cout << "CsSiTrackerDetector[" << GetTBName() << "]::split_cluster()" << endl;
    }

	double amp_big, amp_small;
	if(a2sum1>a2sum2) { amp_small=a2sum1; amp_big=a2sum2; }
	else 			  { amp_big=a2sum1; amp_small=a2sum2; }
		
	CsCluster* cluster = *Iiter;	
		
	std::list<CsDigit*> Digits =cluster->getDigitsList();
	
	double a2sum; cluster->getAnalog(a2sum);

	double a2sum_l=0;
	double a2sum_r=0;
	double a2wsum_l=0;
	double a2wsum_r=0;
	
	double amp2out_l=0;
	double amp1out_l=0;
	double amp0out_l=0;
	
	double amp2out_r=0;
	double amp1out_r=0;
	double amp0out_r=0;
	
	bool forward=false, backward =false;
	
	std::list<CsDigit*>::iterator digiter;
	
	if( (*Digits.begin())->getData()[2] < amp_big )     { digiter = Digits.begin(); forward=true; }	
	else if ( (*Digits.back()).getData()[2] < amp_big ) { digiter = (--Digits.end());  backward=true;}
	
	if( forward || backward ) {
	if (getenv("COND_DB_DEBUG")!=0 || getenv("COND_DB_DEBUG_SI2")!=0 || getenv("COND_DB_DEBUG_SI_MATE")!=0) {
	cout << "CsSiTrackerDetector[" << GetTBName() << "]::split_cluster forward || backward" << endl;
    }
		
	int splitpoint =-1;
	
			for(int i=0; i< int(Digits.size()); i++)
				{	
						a2sum_l+=(*digiter)->getData()[2];
											
						if(a2sum_l < amp_big && !a2sum_r)	// here is the strip to split
						{	// left side cluster
							a2wsum_l+=(*digiter)->getAddress()*(*digiter)->getData()[2];
							
							if (getenv("COND_DB_DEBUG")!=0 || getenv("COND_DB_DEBUG_SI2")!=0 || getenv("COND_DB_DEBUG_SI_MATE")!=0) {
							cout << "wire number / digit address   " << (*digiter)->getAddress() << endl;
						    }
							
							amp2out_l+=(*digiter)->getData()[2];
							amp1out_l+=(*digiter)->getData()[1];
							amp0out_l+=(*digiter)->getData()[0];	
						}
						else if( !a2sum_r)
						{   // here the splitting takes place
						splitpoint=i;

							double quotient = double(amp_small)/double(amp_big+amp_small);
							if(  fabs(int( (1.-quotient)*double((*digiter)->getData()[2])) + amp2out_l -amp_small) < fabs(int( (quotient)*double((*digiter)->getData()[2])) + amp2out_l - amp_small) )
							{ quotient=1.-quotient;}
							if( (int (quotient *double((*digiter)->getData()[2])) < 14 ) )
								{
									if( a2sum_l < a2sum -a2sum_l) quotient =1.0;
									else quotient =0.0;
								}
								
							a2wsum_l+= int (quotient * (*digiter)->getAddress() * (*digiter)->getData()[2] );
							amp2out_l+=int (quotient * (*digiter)->getData()[2] );
							amp1out_l+=int (quotient * (*digiter)->getData()[1] );
							amp0out_l+=int (quotient * (*digiter)->getData()[0] );

							
							quotient = 1. -quotient;
							
							a2wsum_r += int (quotient * (*digiter)->getAddress() * (*digiter)->getData()[2]);
							a2sum_r  += int (quotient * (*digiter)->getData()[2]);
							amp2out_r+= int (quotient * (*digiter)->getData()[2]);
							amp1out_r+= int (quotient * (*digiter)->getData()[1]);
							amp0out_r+= int (quotient * (*digiter)->getData()[0]);
						}
						
						else
						{	// right side cluster
	
							a2wsum_r+=(*digiter)->getAddress()*(*digiter)->getData()[2];
							amp2out_r+=(*digiter)->getData()[2];
							amp1out_r+=(*digiter)->getData()[1];
							amp0out_r+=(*digiter)->getData()[0];		
						
							if (getenv("COND_DB_DEBUG")!=0 || getenv("COND_DB_DEBUG_SI2")!=0 || getenv("COND_DB_DEBUG_SI_MATE")!=0){
							cout << "wire number / digit address   " << (*digiter)->getAddress() << endl;
							}
						}
	
						if(forward)  { digiter++; }
						if(backward) { digiter--; }
				}
							
		HepMatrix cov =cluster->getCov();
		cov(1,1)= 2.3e-04;  // res= 15 um

		double u_corr_l=wireDCorr_;
		double u_corr_r=wireDCorr_;
		
		
		double u_l = double(a2wsum_l)/double(amp2out_l)*double(wirP_)+u_corr_l;
		double u_r = double(a2wsum_r)/double(amp2out_r)*double(wirP_)+u_corr_r;
		u_l = u_l -rescor;
		u_r = u_r -rescor;
		
		cout.precision(6);
		
		double w_l=zcm_;
		double w_r=zcm_;
		cout.precision(3);
		
		CsCluster *cluster_l=new CsCluster(u_l,1.0,w_l,cov);
		CsCluster *cluster_r=new CsCluster(u_r,1.0,w_r,cov);
		
		const string & name = GetTBName();
		h(name+"csplitPos",u_l); h(name+"csplitPos",u_r);

		
		  list<CsDigit*>::const_iterator DigIter;
		  int k=0;
		  for(DigIter=Digits.begin(); DigIter!=Digits.end(); ++DigIter,k++) {
		
		  if      ( k< splitpoint-1 ) { cluster_l->addDigit(*(*DigIter)); }
		  else if ( k==splitpoint-1 ) { cluster_l->addDigit(*(*DigIter)); cluster_r->addDigit(*(*DigIter)); }							
		  else if ( k> splitpoint-1 ) { cluster_r->addDigit(*(*DigIter));  }
		
		  }
	
	
		// TODO for the moment old times and errors are usde...
		// this should be chaned, because time and errors have to recalculated
		// idea: use time of clusters used for splitting
		double clustertime;
		double clustertime_err;
		cluster->getTime(clustertime);
		cluster->getTimeError(clustertime_err);
		
	    h(name+"_csplitTime ",clustertime);			

		
		
		
		cluster_l->setTime(clustertime,clustertime_err);
		cluster_l->setTimeError(clustertime_err);
		
		cluster_r->setTime(clustertime,clustertime_err);
		cluster_r->setTimeError(clustertime_err);
				
		cluster_l->setAnalog(amp2out_l);
		cluster_l->addAnalogData(amp2out_l);
		cluster_l->addAnalogData(amp1out_l);
		cluster_l->addAnalogData(amp0out_l);
		
		cluster_r->setAnalog(amp2out_r);
		cluster_r->addAnalogData(amp2out_r);
		cluster_r->addAnalogData(amp1out_r);
		cluster_r->addAnalogData(amp0out_r);
		
	    h(name+"_cAmptsplit",amp2out_r);
        h(name+"_cAmptsplit",amp2out_l);
	
		// add backreference
		cluster_l->addDet(*this);
		cluster_r->addDet(*this); 	

		// make sure cluster digits are not empty and postition is non NaN
		const int mina2sum=1;

		bool add_l=false;
		bool add_r=false;


			

		if( cluster_l->getDigitsList().size() && cluster_l->getU()==cluster_l->getU() && amp2out_l > mina2sum)
			{
			if (getenv("COND_DB_DEBUG")!=0 || getenv("COND_DB_DEBUG_SI2")!=0 || getenv("COND_DB_DEBUG_SI_MATE")!=0)
				{
				cout.precision(6);		
				cout << " cluster left  :  a2sum_l=" << amp2out_l << "   a2wsum_l=" << a2wsum_l << "   qotient="  << fixed << double(a2wsum_l)/double(amp2out_l) << endl;
				cout << "CsSiTrackerDetector[" << GetTBName() << "] old u was "  << cluster->getU()<< endl;
				cout << "CsSiTrackerDetector[" << GetTBName() << "] new u is "  << cluster_l->getU() << endl;
			}
						
			add_l=true;
			h(name+"cSizesplit", cluster_l->getDigitsList().size());			
			vtent_cl.push_back(cluster_l);
			}
			else{
				delete cluster_l;
			}

		if( cluster_r->getDigitsList().size() && cluster_r->getU()==cluster_r->getU() && amp2out_r > mina2sum)
		{
			if (getenv("COND_DB_DEBUG")!=0 || getenv("COND_DB_DEBUG_SI2")!=0 || getenv("COND_DB_DEBUG_SI_MATE")!=0)
			{
			cout.precision(6);		
			cout << " cluster right :  a2sum_r=" << amp2out_r << "   a2wsum_r=" << a2wsum_r << "   qotient="  << fixed << double(a2wsum_r)/double(amp2out_r) << endl;
	        cout << "CsSiTrackerDetector[" << GetTBName() << "] old u was "  << cluster->getU()<< endl;
	        cout << "CsSiTrackerDetector[" << GetTBName() << "] new u is "  << cluster_r->getU() << endl;
			}
			add_r=true;	
			h(name+"cSizesplit", cluster_r->getDigitsList().size());
			vtent_cl.push_back(cluster_r);
		}
		else{
			 delete cluster_r;
		}

		if (getenv("COND_DB_DEBUG")!=0 || getenv("COND_DB_DEBUG_SI2")!=0 || getenv("COND_DB_DEBUG_SI_MATE")!=0)
		{
			double a2sum; cluster->getAnalog(a2sum);
			cout.precision(1);
	        cout << "CsSiTrackerDetector[" << GetTBName() << "] old a2sum was "  << a2sum<< endl;
	        cout << "CsSiTrackerDetector[" << GetTBName() << "] new a2sum is "  << amp2out_l << " + " << amp2out_r << endl;
		}
		if(add_l)    { h(name+"Nsplit",1);	 h(name+"cAmpsplit",amp2out_l);}
		if(add_r) 	 { h(name+"Nsplit",1);	 h(name+"cAmpsplit",amp2out_r);}
		
		
		// if at least one new cluster was created and added, delete the old one
		if( add_l || add_r ) {
							  delete *(Iiter);
							  vtent_cl.erase(Iiter);
							 }
	
	}

}


//------------------------------------------------------------
void CsSiTrackerDetector::split_cluster(std::list<CsCluster*> & vtent_cl, std::list<CsCluster*>::iterator Iiter, double wireDCorr) {
//------------------------------------------------------------
//	cluster specified by Iiter is split up by factor 0.5
//	this version is only called when CsSiTrackerDetector::find_best_candidates() was not sucessful
//	2 new clusters are created and stored in the cluster list
//	old cluster is deleted
	
	if (getenv("COND_DB_DEBUG")!=0 || getenv("COND_DB_DEBUG_SI2")!=0 || getenv("COND_DB_DEBUG_SI_MATE")!=0){
	cout << "CsSiTrackerDetector[" << GetTBName() << "]::split_cluster()" << endl;
	}
	
	CsCluster* cluster = *Iiter;	
		
	std::list<CsDigit*> Digits =cluster->getDigitsList();
	
	double a2sum; cluster->getAnalog(a2sum);

	// a2sums and a2-wire sums.
	// for position of new clusters
	double a2wsum_l=0;
	double a2wsum_r=0;
	
	double amp2out_l=0;
	double amp1out_l=0;
	double amp0out_l=0;
	
	double amp2out_r=0;
	double amp1out_r=0;
	double amp0out_r=0;
	
	std::list<CsDigit*>::iterator digiter;
	digiter = Digits.begin();
	
	if (getenv("COND_DB_DEBUG")!=0 || getenv("COND_DB_DEBUG_SI2")!=0 || getenv("COND_DB_DEBUG_SI_MATE")!=0){
	cout << "CsSiTrackerDetector[" << GetTBName() << "]::split_cluster"<< endl;
	}
		// the point to split is defined.
		// for clusters of odd size, the middle strip is split with factor 0.5
		// for clusters of even side, two clusters with the identical number of strips are created
		
		int dummy_splitpoint;
		if(Digits.size()%2)	{dummy_splitpoint=(Digits.size()-1)/2 +1;}
		else {dummy_splitpoint=(Digits.size())/2;}
		
		const int splitpoint=dummy_splitpoint;
		
			// loop over Digits of old cluster (the one which is to split)
			for(int i=0; i< int(Digits.size()); i++, digiter++)
				{	
						if( i < splitpoint-1)	
							{	
							// in this case the digits are added to the left cluster	
								
							a2wsum_l+=(*digiter)->getAddress()*(*digiter)->getData()[2];
							
								if (getenv("COND_DB_DEBUG")!=0 || getenv("COND_DB_DEBUG_SI2")!=0 || getenv("COND_DB_DEBUG_SI_MATE")!=0)
								{
								cout << "wire number / digit address   " << (*digiter)->getAddress() << endl;
								}

							amp2out_l+=(*digiter)->getData()[2];
							amp1out_l+=(*digiter)->getData()[1];
							amp0out_l+=(*digiter)->getData()[0];	
							}
						
						else if( i == splitpoint-1)
							{
							// this is the strip where splitting takes place	
							
							double quotient;
							if(Digits.size()%2) { quotient = 0.5; } // cluster size is odd, middle strip will be split with 0.5
							else 				{ quotient = 1.0; }	// cluster size is even

							// update wire sums of left cluster
							a2wsum_l+= (*digiter)->getAddress() * int ( quotient*(*digiter)->getData()[2]);
							amp2out_l+= int ( quotient*(*digiter)->getData()[2]);
							amp1out_l+= int ( quotient*(*digiter)->getData()[1]);
							amp0out_l+= int ( quotient*(*digiter)->getData()[0]);
							
							
							// update wire sum of right cluster
							// "rest" of a2 is used
							quotient = 1. -quotient;

							a2wsum_r+= (*digiter)->getAddress() * int ( quotient*(*digiter)->getData()[2]);
							amp2out_r+= int ( quotient*(*digiter)->getData()[2]);
							amp1out_r+= int ( quotient*(*digiter)->getData()[1]);
							amp0out_r+= int ( quotient*(*digiter)->getData()[0]);

							}
						
						else
						{	
							//this digits are added to the right cluster

							a2wsum_r+=(*digiter)->getAddress()*(*digiter)->getData()[2];
							amp2out_r+=(*digiter)->getData()[2];
							amp1out_r+=(*digiter)->getData()[1];
							amp0out_r+=(*digiter)->getData()[0];		
						
							if (getenv("COND_DB_DEBUG")!=0 || getenv("COND_DB_DEBUG_SI2")!=0 || getenv("COND_DB_DEBUG_SI_MATE")!=0) {
							cout << "wire number / digit address   " << (*digiter)->getAddress() << endl;
							}	
						}
				}
							
		HepMatrix cov =cluster->getCov();
		cov(1,1)= 2.3e-04;
	
		double u_corr_l=wireDCorr_;
		double u_corr_r=wireDCorr_;
		
		
		double u_l = double(a2wsum_l)/double(amp2out_l)*double(wirP_)+u_corr_l;
		double u_r = double(a2wsum_r)/double(amp2out_r)*double(wirP_)+u_corr_r;
		u_l = u_l -rescor;
		u_r = u_r -rescor;
		

		double w_l=zcm_;
		double w_r=zcm_;
		cout.precision(3);
		
		CsCluster *cluster_l=new CsCluster(u_l,2.0,w_l,cov);
		CsCluster *cluster_r=new CsCluster(u_r,2.0,w_r,cov);
		
		const string & name = GetTBName();
		h(name+"csplitPos",u_l); h(name+"csplitPos",u_r);

		
	
		list<CsDigit*>::const_iterator DigIter;
		  int k=0;
		  for(DigIter=Digits.begin(); DigIter!=Digits.end(); ++DigIter,k++) {
		
		  if      ( k< splitpoint-1 ) { cluster_l->addDigit(*(*DigIter)); }
		  else if ( k==splitpoint-1 ) { cluster_l->addDigit(*(*DigIter)); cluster_r->addDigit(*(*DigIter)); }							
		  else if ( k> splitpoint-1 ) { cluster_r->addDigit(*(*DigIter));  }
		
		  }
		
		// TODO for the moment old times and errors are usde...
		// this should be chaned, because time and errors have to recalculated
		// idea: use time of clusters used for splitting
		double clustertime;
		double clustertime_err;
		cluster->getTime(clustertime);
		cluster->getTimeError(clustertime_err);
		
	    h(name+"_csplitTime ",clustertime);			

		
		
		
		cluster_l->setTime(clustertime,clustertime_err);
		cluster_l->setTimeError(clustertime_err);
		
		cluster_r->setTime(clustertime,clustertime_err);
		cluster_r->setTimeError(clustertime_err);
				
		cluster_l->setAnalog(amp2out_l);
		cluster_l->addAnalogData(amp2out_l);
		cluster_l->addAnalogData(amp1out_l);
		cluster_l->addAnalogData(amp0out_l);
		
		cluster_r->setAnalog(amp2out_r);
		cluster_r->addAnalogData(amp2out_r);
		cluster_r->addAnalogData(amp1out_r);
		cluster_r->addAnalogData(amp0out_r);
		
	    h(name+"_cAmptsplit",amp2out_r);
        h(name+"_cAmptsplit",amp2out_l);
	
		// add backreference
		cluster_l->addDet(*this);
		cluster_r->addDet(*this); 	

		// make sure cluster digits are not empty and postition is non NaN
		const int mina2sum=1;

		bool add_l=false;
		bool add_r=false;

	
	
		if( cluster_l->getDigitsList().size() && cluster_l->getU()==cluster_l->getU() && amp2out_l > mina2sum)
			{
			if (getenv("COND_DB_DEBUG")!=0 || getenv("COND_DB_DEBUG_SI2")!=0 || getenv("COND_DB_DEBUG_SI_MATE")!=0)
				{
				cout.precision(6);		
				cout << " cluster left  :  a2sum_l=" << amp2out_l << "   a2wsum_l=" << a2wsum_l << "   qotient="  << fixed << double(a2wsum_l)/double(amp2out_l) << endl;
				cout << "CsSiTrackerDetector[" << GetTBName() << "] old u was "  << cluster->getU()<< endl;
				cout << "CsSiTrackerDetector[" << GetTBName() << "] new u is "  << cluster_l->getU() << endl;
				}
						
			add_l=true;
			h(name+"cSizesplit", cluster_l->getDigitsList().size());			
			vtent_cl.push_back(cluster_l);
			}
		else{	
			 delete cluster_l;
		}

		if( cluster_r->getDigitsList().size() && cluster_r->getU()==cluster_r->getU() && amp2out_r > mina2sum)
		{
			if (getenv("COND_DB_DEBUG")!=0 || getenv("COND_DB_DEBUG_SI2")!=0 || getenv("COND_DB_DEBUG_SI_MATE")!=0)
			  {
			  cout.precision(6);		
			  cout << " cluster right :  a2sum_r=" << amp2out_r << "   a2wsum_r=" << a2wsum_r << "   qotient="  << fixed << double(a2wsum_r)/double(amp2out_r) << endl;
	          cout << "CsSiTrackerDetector[" << GetTBName() << "] old u was "  << cluster->getU()<< endl;
	          cout << "CsSiTrackerDetector[" << GetTBName() << "] new u is "  << cluster_r->getU() << endl;
			  }
			
			add_r=true;	
			h(name+"cSizesplit", cluster_r->getDigitsList().size());
			vtent_cl.push_back(cluster_r);
		}
		else{
			 delete cluster_r;
		}	
			
		if(add_l)     h(name+"Nsplit",1);	
		if(add_r) 	  h(name+"Nsplit",1);	
		
		// if at least one new cluster was created and added, delete the old one
		if( add_l || add_r ) {
							 delete *(Iiter);	
							 vtent_cl.erase(Iiter);
						  	 }

		
}

	




//------------------------------------------------------------
bool CsSiTrackerDetector::find_best_candidates(std::list<CsCluster*> & vtent_cl_mate,list<list<CsCluster*>::iterator>::iterator Iiter, int & a2sum1, int& a2sum2) {
//------------------------------------------------------------
// finds best candidates in vtent_cl_mate for splitting cluster Iiter


    if (getenv("COND_DB_DEBUG")!=0 || getenv("COND_DB_DEBUG_SI2")!=0 || getenv("COND_DB_DEBUG_SI_MATE")!=0 )
      {
	  cout << "CsSiTrackerDetector[" << GetTBName() << "]::find_best_candidate()" << endl;
	  }
	
	
	const string & name = GetTBName();
	a2sum1=0;
	a2sum2=0;

	const int min_a2sum=20;

	CsCluster* cluster = *(*Iiter);
	double a2sum;
	cluster->getAnalog(a2sum);
	
	double analog1=-1.;
	double analog2=-1.;
	
	for(std::list<CsCluster*>::iterator Iclus=vtent_cl_mate.begin(); Iclus!= vtent_cl_mate.end(); Iclus++) {
		
			if(!(*Iclus)->getAnalog(analog1))							{ continue; }
			if(analog1< min_a2sum || (*Iclus)->getDigitsList().size()>2) { continue; }
			
			for(std::list<CsCluster*>::iterator Iclus2=vtent_cl_mate.begin(); Iclus2!= vtent_cl_mate.end(); Iclus2++)
			{
				if( (*Iclus) == (*Iclus2) ) 		{ continue; }
				if(!(*Iclus2)->getAnalog(analog2)) 	{ continue; }
				
				if(analog2< min_a2sum || (*Iclus2)->getDigitsList().size()>2) { continue; }
			
					if( abs( (a2sum1+a2sum2) -a2sum) > abs ( (analog1+analog2) -a2sum))
						{
						a2sum1 = (int)analog1;
						a2sum2 = (int)analog2;
						}
				
			}
	}

	
    if (getenv("COND_DB_DEBUG")!=0 || getenv("COND_DB_DEBUG_SI2")!=0 || getenv("COND_DB_DEBUG_SI_MATE")!=0 )
      {
	  cout << " in find best   a2sum1  " << a2sum1 << "    a2sum2  " << a2sum2 << endl;
	  cout << " in find best    sum:  " <<a2sum1+a2sum2 << "   and   " << a2sum  << endl;
	  }
	
	h(name+"cA2_1",analog1);
   	h(name+"cA2_2",analog2);
    h(name+"cA2_wide",a2sum);
    h(name+"cA2_diff",a2sum-(analog1+analog2));
    h(name+"a21:a22", analog1, analog2);
    h(name+"diff_vs_a2sum", a2sum-(analog1+analog2), a2sum);
    h(name+"large_vs_a2sum", a2sum, (analog1+analog2) );

	if(analog1 <30 || analog2<30 ) return false;
	if( fabs( a2sum-(analog1+analog2) ) > 29) return false;



return true;
}





//------------------------------------------------------------
void CsSiTrackerDetector::split_delta(std::list<CsCluster*> & vtent_cl, std::list<CsCluster*>::iterator Iiter, double wireDCorr) {
//------------------------------------------------------------
//	In case mate splitting fails, clusters are nevertheless split up.
//	For clusters of size 4 or lager the cluster is most probably
//	due to a delta electron, meaning that the track passed either
//	on the left or the right side while the strips in the middle
//  do not contain much information. 2 new clusters are created
//  each with size 2, consisting of the first two and the last two
//  strips of the old cluster

	const string & name = GetTBName();

	h(name+"_delta_called",1);
	
	CsCluster* cluster = *Iiter;	
	std::list<CsDigit*> Digits =cluster->getDigitsList();
	
	// size must be larger than 4, otherwise another function is used
	if(Digits.size()<4) { return; }

	const int mina2sum=1;
		
	double a2sum; cluster->getAnalog(a2sum);
	
	h(name+"_delta_pos_in", cluster->getU());

 	// amplitude sums and wire sums
	double amp2out_l = 0;
	double amp1out_l = 0;
	double amp0out_l = 0;	
	double a2wsum_l= 0.;
	
	double amp2out_r = 0;
	double amp1out_r = 0;
	double amp0out_r = 0;	
	double a2wsum_r= 0.;

	std::list<CsDigit*>::iterator digiter;
	
	bool print =true;
		
		
		// -------------------------------------------------------------------------------------------------------------------------
		// left cluster
		//
	   digiter=Digits.begin();
	   a2wsum_l+= (*digiter)->getAddress() * (*digiter)->getData()[2];
	
	   amp2out_l+=(*digiter)->getData()[2];
	   amp1out_l+=(*digiter)->getData()[1];
	   amp0out_l+=(*digiter)->getData()[0];
	   						

	   digiter++;
	   a2wsum_l+= (*digiter)->getAddress() * (*digiter)->getData()[2];
	
	   amp2out_l+=(*digiter)->getData()[2];
	   amp1out_l+=(*digiter)->getData()[1];
	   amp0out_l+=(*digiter)->getData()[0];
	   						


		// -------------------------------------------------------------------------------------------------------------------------
		// right cluster
		//

       // to have the digits still sorted
	   digiter=Digits.end();
	   digiter--;
	   digiter--;

	   a2wsum_r+= (*digiter)->getAddress() * (*digiter)->getData()[2];
	
	   amp2out_r+=(*digiter)->getData()[2];
	   amp1out_r+=(*digiter)->getData()[1];
	   amp0out_r+=(*digiter)->getData()[0];
	   			
	   						
	   digiter++;
	   a2wsum_r+= (*digiter)->getAddress() * (*digiter)->getData()[2];
	
	   amp2out_r+=(*digiter)->getData()[2];
	   amp1out_r+=(*digiter)->getData()[1];
	   amp0out_r+=(*digiter)->getData()[0];
	   						
							
		// define suitable error					
		HepMatrix cov =cluster->getCov();
		cov(1,1)= 2.3e-04;

		double u_corr_l=wireDCorr_;
		double u_corr_r=wireDCorr_;
		
		
		double u_l = double(a2wsum_l)/double(amp2out_l)*double(wirP_)+u_corr_l;
		double u_r = double(a2wsum_r)/double(amp2out_r)*double(wirP_)+u_corr_r;
		u_l = u_l -rescor;
		u_r = u_r -rescor;
		
		// histograms for debug
		h(name+"_delta_pos_1",u_l);
		h(name+"_delta_pos_2",u_r);

				
		double w_l=zcm_;
		double w_r=zcm_;
		cout.precision(3);
		
		CsCluster *cluster_l=new CsCluster(u_l,3.0,w_l,cov);
		CsCluster *cluster_r=new CsCluster(u_r,3.0,w_r,cov);
		
		list<CsDigit*>::const_iterator DigIter;
	
		DigIter=Digits.begin();
		cluster_l->addDigit(*(*DigIter));
		DigIter++;
		cluster_l->addDigit(*(*DigIter));


		DigIter=Digits.end();
		DigIter--;
		DigIter--;
		cluster_r->addDigit(*(*DigIter));
		DigIter++;
		cluster_r->addDigit(*(*DigIter));
					
					
		// TODO for the moment old times and errors are usde...
		// this should be chaned, because time and errors have to recalculated
		// idea: use time of clusters used for splitting
		double clustertime;
		double clustertime_err;
		cluster->getTime(clustertime);
		cluster->getTimeError(clustertime_err);
		
	
	
		// Set Time, and Amplitude Information of new Clusters
		
		cluster_l->setTime(clustertime,clustertime_err);
		cluster_l->setTimeError(clustertime_err);
		
		cluster_r->setTime(clustertime,clustertime_err);
		cluster_r->setTimeError(clustertime_err);
				
		cluster_l->setAnalog(amp2out_l);
		cluster_l->addAnalogData(amp2out_l);
		cluster_l->addAnalogData(amp1out_l);
		cluster_l->addAnalogData(amp0out_l);
		
		cluster_r->setAnalog(amp2out_r);
		cluster_r->addAnalogData(amp2out_r);
		cluster_r->addAnalogData(amp1out_r);
		cluster_r->addAnalogData(amp0out_r);

		h(name+"_delta_A2_1",amp2out_r);
		h(name+"_delta_A2_2",amp2out_l);
		h(name+"_delta_A2_diff", a2sum-(amp2out_r+amp2out_l));

		//h(name+"_delta_cls_1", new_Digits.size());
		//h(name+"_delta_cls_2", new_Digits2.size());
		h(name+"_delta_cls_in",Digits.size());


		
		// Add backreference to detector
		cluster_l->addDet(*this);
		cluster_r->addDet(*this); 	


		bool add_l=false;
		bool add_r=false;
		
		
		// make sure cluster digits are not empty and postition is non NaN
		if( cluster_l->getDigitsList().size() && cluster_l->getU()==cluster_l->getU() && amp2out_l > mina2sum)
			{
			add_l=true;
			vtent_cl.push_back(cluster_l);
			}
		else{	
			delete cluster_l;
		}

		if( cluster_r->getDigitsList().size() && cluster_r->getU()==cluster_r->getU() && amp2out_r > mina2sum)
			{
			add_r=true;	
			vtent_cl.push_back(cluster_r);
			}
		else{		
			delete cluster_r;
		}
		
		
		// if at least one new cluster was created and added, delete the old one
		if( add_l || add_r ) {
						     delete *(Iiter);
						     vtent_cl.erase(Iiter);
						     }
		
} // end of delta splitting


//------------------------------------------------------------
double CsSiTrackerDetector::cls( double *x, double *par){
//------------------------------------------------------------
return par[2]*(TMath::Erf((-x[0]+par[0])/par[1]))  + par[5]*(TMath::Erf((x[0]-par[3])/par[4])) +par[6] ;
}

//------------------------------------------------------------
double CsSiTrackerDetector::cls3( double *x, double *par){
//------------------------------------------------------------

return (par[0]*exp( -1.*exp(- (  (par[1]+par[2])/2.*(x[0]-par[3])+(par[1]-par[2])/2.*sqrt( ((x[0]-par[3])*(x[0]-par[3])+par[4]*par[4])-sqrt(par[4]*par[4]))+par[5]))));
}



//------------------------------------------------------------
double CsSiTrackerDetector::eta( double *x, double *par){
//------------------------------------------------------------

  if( isU || isX )
  {
  return par[0]*(1.-2./(exp(par[1]*(x[0]-par[2]))+1.))+par[3];
  }

  else if( isV || isY )
  {
  return (par[0]*(1.-2./(exp(2.*(par[1]-x[0])/par[2])+1.))+par[3])* (par[4]*(1.-2./(exp(2.*(par[5]-x[0])/par[6])+1.))+par[7]);
  }

  else
  {
  cout << "return 0 " << endl;
  return 0;
  }
}

//------------------------------------------------------------
double CsSiTrackerDetector::get_charge_sharing(double fdist){
//------------------------------------------------------------

	    // charge sharing depends on dectector type, namely on the presence
  	    // of intermediate strips
	    // U, X : NO intermediate strips present
        // V, Y :    intermediate strips present

        double feta=-1;
        const string& name = GetTBName();


	    if( isU || isX )
	      {
	      // randomize distance to main wire
          double shift = 20*CsRandom::gauss();
          double fweight = (-0.25*fabs(fdist-26.05)+26.05)/26.05-0.7;

          shift = fweight*shift;
	
	      // obtain mean value of eta according to fitted function	
          double x_new = fdist+shift;
          feta = eta(&x_new, eta_);

	      feta= feta + 0.015 *CsRandom::gauss();

  	      h(name+"mc_eta_pos", fdist, feta);
  	      h(name+"mc_eta",feta);
 	      }

	    if( isV || isY )
          {
		  double shift= 2.5*CsRandom::gauss();
		
		  // randomize distance to next main wire, but not smaller than 0
		  // or larger than pitch
		  	if( fdist+shift<0.0) { fdist= fdist-shift;}
		  	else if (fdist+shift >wirP_*10000) { fdist=fdist-shift;}
		
		  // get mean of eta distribution according to distance to main wire
		  double eta_1 = eta(&fdist,eta_);
		
		  // randomize result
		  double eta_shift=0.06*CsRandom::gauss() ;
		
		  // shift outliers of distribution, since these cause unphysical structures
		  	if(fdist<20. && eta_shift<0.0) { eta_shift = eta_shift * (fabs(fdist-20.)+9)/9;         }
		  	if(fdist>30. && eta_shift>0.0) { eta_shift = eta_shift * (fabs(fdist-30.)+9)/9;         }
		
		  // increase small shitfs only, otherwise spikes migth occur
		  	if(fabs(eta_shift)<0.01) {eta_shift = eta_shift*1.2;}
		
		  // finally, the value of feta
		  feta=eta_1 + eta_shift;
		
		  // in the data is only between 0.04 and 0.96
		  // this should in principle be displayed in MC, but part of it is due to zero suppresion
		  // where the remaining amplitude is in the order of noise. So a slightly larger range
		  // is simulated to accout for this.
		  if(feta>0.98) { feta = eta_1 -2*eta_shift;}
		  if(feta<0.02) { feta = eta_1 -2*eta_shift;}
		
		  // fill histograms
		  h(name+"mc_eta_pos", fdist, feta);
		  h(name+"mc_eta",feta);
		  }
			
			
			return feta;
	
}

