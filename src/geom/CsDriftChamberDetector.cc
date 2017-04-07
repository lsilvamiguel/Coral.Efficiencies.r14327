// $Id: CsDriftChamberDetector.cc,v 1.117 2010/01/28 12:51:25 tnagel Exp $

/*!
   \file    CsDriftChamberDetector.cc
   \brief   Compass Drift Chamber like detector Class.
   \author  Benigno Gobbo
   \version $Revision: 1.117 $
   \date    $Date: 2010/01/28 12:51:25 $
*/

#include <math.h>
#include <string.h>
#include <strings.h>
#include <functional>
#include "CsDriftChamberDetector.h"
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
//CY
#include "DaqDataDecoding/ChipSinica.h"
//LS Eff
#include <stdexcept>
//#define DEBUG_EFF
#include "CsRandom.h"

using namespace std;
using namespace CLHEP;

extern QhitType Qhit;

//____________________________________________________________________________
CsDriftChamberDetector::CsDriftChamberDetector( const int    row,
      const int    id,    const char* name,   const char *TBname, 
      const int    unit,  const int    type,
      const double rdLen, const double xsiz,  
      const double ysiz,  const double zsiz,
      const double xcm,   const double ycm,   
      const double zcm,   const HepMatrix rotDRS,
      const HepMatrix rotWRS,
      const double wirD,  const double ang,   
      const int    nWir,  const double wirP, 
      const double eff,   const double bkg,
      const double tGate,
      const double vel,   const double t0,
      const double thRes, const double spSli, 
      const double tiSli ) :
  CsDetector( row, id, name, TBname, unit, type, rdLen, xsiz, ysiz, zsiz,
	      xcm, ycm, zcm, rotDRS, rotWRS, wirD, ang,
	      nWir, wirP, eff, bkg, tGate ),
  isMWPC_( false ),
  alwaysTwo_( false ),
  vel_(vel), t0_(t0), thRes_(thRes), spSli_(spSli), tiSli_(tiSli),
  hasRT_(false), rt_(NULL),
  hasAssociateDet_(false), associateDet_(NULL),
  associationCut_(-1.0), LRProbCut_(-1.0), LRMode_(-1),
  propVel_inv_( 0.005 ),   //Unit is ns/mm
  feOrientation_( 0 ),

  //=== histogram initialisation to NULL (0)
  H_t_( 0 ),  H_ch_( 0 ),
  H_t_cor_( 0 ), H_ch_cor_( 0 ),
  H_n_cor_( 0 ),
  H_LR_( 0 ), H_LR_duMin_( 0 ), H_LR_eff_( 0 ),
  H_t_vs_ch_( 0 )
  
{
  CsOpt *opt = CsOpt::Instance( );

  decodeCard_ = false;                            // ***** DECODING METHOD *****
  bool status = opt->getOpt( "" , "make decoding" );
  if (status) {
    list<string> options; list<string>::iterator Is;
    if ((status = opt->getOpt("","make decoding",options))) {
      for (Is = options.begin(); Is!=options.end(); Is++ ) 
	if (*Is=="DC" || *Is=="all") decodeCard_ = true;
    }
    else                             decodeCard_ = true;
  }

  //                                          ***** MULTI-HIT MC DECODING? *****
  multiHitMCDecoding_ = false; 
  if (opt->getOpt(GetTBName(),            "multiHitMCDecoding") ||
      opt->getOpt(GetTBName().substr(0,2),"multiHitMCDecoding")) {
    multiHitMCDecoding_ = true;
  }
  else CsErrLog::msg(elWarning,__FILE__,__LINE__,
		     "%s: multi-hit MC decoding will NOT be used",TBname);
  //LS Eff
  //                                          ***** eff. MAP ENABLED? *****
  mcEffMapsEnabled_ = false;
  if (opt->getOpt("ALL"      ,            "mcEffMapsEnabled") ||
      opt->getOpt(GetTBName(),            "mcEffMapsEnabled") ||
      opt->getOpt(GetTBName().substr(0,2),"mcEffMapsEnabled")) {
    mcEffMapsEnabled_ = true;
#ifdef DEBUG_EFF
    cout<<"**** LS *** : CsDriftChamberDetector:: "<< mcEffMapsEnabled_<<endl;  
#endif
  }
  else CsErrLog::msg(elWarning,__FILE__,__LINE__,
                     "%s: eff. Map will NOT be used",TBname);

  //                                         ***** DC TREATED LIKE a MWPC? *****
  isMWPC_ = opt->getOpt(GetTBName(),            "IS_MWPC") ||
            opt->getOpt(GetTBName().substr(0,2),"IS_MWPC");
  if (isMWPC_)
    CsErrLog::msg(elWarning,__FILE__,__LINE__,"%s: treated as MWPC.",TBname);

  badChannels_.clear();                    // ***** INITIALISE badChannels *****
  
  //                         ********** RT RELATION **********
  if (!isMWPC_ && !_readRTRelation()) {      // ***** RT from OPTIONS FILE *****
    if (CsInit::Instance()->IsADataJob())
      CsErrLog::msg(elInfo,__FILE__,__LINE__,
	"%s: No RT relation in options file. Will search DB for it.",TBname);
    else {                                            // ***** MC CASE ... *****
      rt_ = new CsRTRelation(*this); hasRT_ = true; // ***** ...DEFAULT RT *****
    }
  }
  //=== Check alwaysTwo options ===
  alwaysTwo_ = opt->getOpt("DC","make always two clusters");

  //   ********** INITIALISATION for PROPAGATION TIME CORRECTION **********
  if (CsInit::Instance()->IsADataJob() // ***** MC...
      && !isMWPC_)    // ...or IS MWPC: PROPAGATION TIME = 0...
    feOrientation_ =  // ... ELSE: FE ORIENTATION    *****
      _getFEOrientation();
  //                                            ***** PROPAGATION VELOCITY *****
  if (opt->getOpt(TBname,                 "PropVelInv",propVel_inv_) ||
      opt->getOpt(GetTBName().substr(0,2),"PropVelInv",propVel_inv_)) {
    CsErrLog::msg(elWarning,__FILE__,__LINE__,
		  "%s: inverse Propagation Velocity = %f",TBname,propVel_inv_);
  } 
  //                           ***** DERIVATIVES w.r.t. MRS @ U = constant *****
  double dTdV = feOrientation_*propVel_inv_;  // Derivative w.r.t. coord+pi/2
  if (fabs(ang_-90)<1) { dTdX_ = -dTdV; dTdY_ = 0;  }  // Y coordinate Chambers
  else { dTdX_ = 0    ; dTdY_ = dTdV/cosAng_; }        // Other coordinates

  {
    //            ********** CLUSTERISATION OPTIONS **********

    // ***** Multi-Hit: "DeadTime" key *****
    // This cut can be made dependent upon detector's TB name, e.g:
    // DC       DeadTime 50
    // DC01V1__ DeadTime  0

    deadTime_ = 0;                 // Default is 0: No (software) Dead Time
    if (opt->getOpt(TBname,                 "DeadTime",deadTime_) ||
	opt->getOpt(GetTBName().substr(0,2),"DeadTime",deadTime_)) {
      CsErrLog::msg(elWarning,__FILE__,__LINE__,
		    "%s: DeadTime = %f",TBname,deadTime_);
    }
  }  
}

//____________________________________________________________________________
bool CsDriftChamberDetector::operator==( const CsDriftChamberDetector& det ) const {
  return( CsDetector::operator==(det) );
}

//____________________________________________________________________________
bool CsDriftChamberDetector::operator<( const CsDriftChamberDetector& det ) const {
  return( CsDetector::operator<(det) );
}

//____________________________________________________________________________
void CsDriftChamberDetector::BookHistograms() {

  //=== check if histograms are to be booked ===
  CsDetector::ReadHistLevel();
  if( hLevel_ == None ) return;

  CsHistograms::SetCurrentPath( "/CsDriftChamber" );
  
  if( hLevel_ >= Normal ) {

    char name[100];
    sprintf( name, "%s_t_", GetTBName().c_str() );
    H_t_ = new CsHist1D( name, name, (int) (3*tGate_) , -tGate_, 2*tGate_ );

    sprintf( name, "%s_t_cor_", GetTBName().c_str() );
    H_t_cor_ = new CsHist1D( name, name, (int) (3*tGate_) , -tGate_, 2*tGate_ );

    sprintf( name, "%s_ch_", GetTBName().c_str() );
    H_ch_ = new CsHist1D( name, name, nWir_, 0, nWir_ );
    sprintf( name, "%s_ch_cor_", GetTBName().c_str() );
    H_ch_cor_ = new CsHist1D( name, name, nWir_, 0, nWir_ );
    
    sprintf( name, "%s_n_cor_", GetTBName().c_str() );
    H_n_cor_ = new CsHist1D( name, name, 20, 0, 20 );
  
    sprintf( name, "%s_r_", GetTBName().c_str() );
    H_r_ = new CsHist1D( name, name, (int) (20*wirP_) , -0.05*wirP_, 0.6*wirP_ );
    
    sprintf( name, "%s_LR_", GetTBName().c_str() );
    H_LR_ = new CsHist1D( name, name, (int) (20*wirP_) , 0, wirP_ );
     
    sprintf( name, "%s_LR_duMin_", GetTBName().c_str() );
    H_LR_duMin_ = new CsHist1D( name, name, (int) (20*wirP_) , 0, wirP_ );
    
    sprintf( name, "%s_LR_eff_", GetTBName().c_str() );
    H_LR_eff_ = new CsHist1D( name, name, (int) 100 , 0, 1 );
    
    sprintf( name, "%s_ch_vs_r_", GetTBName().c_str() );
    H_ch_vs_r_ = new CsHist2D( name, name, nWir_, 0, nWir_, (int) 200 , 0.42*wirP_, 0.52*wirP_ );
  
  }

  if( hLevel_ >= High ) {
    
    char name[100];
    sprintf( name, "%s_t_vs_ch_", GetTBName().c_str() );
    H_t_vs_ch_ = new CsHist2D( name, name, nWir_, 0, nWir_,  (int) (3*tGate_) , -tGate_, 2*tGate_ );
  
  }

  CsHistograms::SetCurrentPath("/");
  return;
}

// Local class for tempopary digits in the detector response simulation 
// (Used in CsDriftChamberDetector::makeMCDecoding())
class DCdig
{
public:
  DCdig(int w, double t, CsMCHit *r, int l) {
    wire = w; time = t; ref = r; lr = l;
  }
  int    wire;  // Wire #
  double time;  // Drift time
  CsMCHit *ref; // Reference to MCHit 
  int    lr;    // Left/right
  bool operator < (const DCdig &gd) const 
  { 
    if(wire == gd.wire) return (time < gd.time);
    else                return (wire < gd.wire);
  };
};

//____________________________________________________________________________
void CsDriftChamberDetector::makeMCDecoding() {

  if (!decode_ && !decodeCard_) return;   // Should I proceed?
  if (decodingDone_) return;              // Already done?
  myDigits_.clear();                      // Clear

  double vel = vel_ / tiSli_;           // cast velocity in mm/ns
  double dtres = tiSli_ * thRes_;       // double hit resolution, ns  

  // The following can be done only with zebra (as opposed to NTuple) files... 
  if (!CsGeant3::Instance()->isAnNtFile()) {

    list<DCdig> ld; // ***** LIST OF TEMPORARY DIGITS... *****
    // (...i.e., digits before double hit resolution is applied.) 
    list<CsMCHit*>::iterator Ih;
    double triggerOffset = // Time offset of the trigger w.r.t. the event
      CsEvent::Instance()->getTriggerMCOffset();

    if (!multiHitMCDecoding()) {  // ********** SIMPLISTIC DECODING **********
      // (N.B.: This decoding method is deprecated (having been replaced by
      // "multiHitMCDecoding", cf. infra) and hasn't been tested for some time.)

      for (Ih = myMCHits_.begin(); Ih!=myMCHits_.end(); Ih++) {
	if ((((*Ih)->getMCTrack())->getParticle())->getCharge() ||
	    (*Ih)->getOrigin()  ) {

	  // ***** LOOP ON HITS from charged particles or charged products *****

	  double t  = (*Ih)->getDTime();            // Delay time (ns)
	  double x  = (*Ih)->getX() - _xshift;      // Hit coordinates in (MRS)
	  double y  = (*Ih)->getY() - _yshift;
	  double z  = (*Ih)->getZ();

	  double u, v, w;
	  rotateVectMRS2WRS (x, y, z, u, v, w);

	  int wire;      // Nearest wire
	  if ((u-wirD_)/wirP_<0) wire = int((u-wirD_)/wirP_-0.5);
	  else                   wire = int((u-wirD_)/wirP_+0.5);
	  if (wire<0 || wire>=nWir_) {
	    CsErrLog::msg(elAnomaly,__FILE__,__LINE__,
	   "%s: Wire(=%d) outside range [0,%d[",GetTBName().c_str(),wire,nWir_);
	    continue;
	  }

	  // Drift time in ns: distance/velocity + TOF - trigger offset
	  double driftExact = fabs(u-(wire*wirP_+wirD_))/vel + t -
	    // Trigger time (and hence offset) is subtracted from hit times
	    triggerOffset;

	  // Is the hit inside the detector gate length?
	  if (driftExact>tGate_ || driftExact<0) continue;

	  //                                          ***** SAVE TEMPORARY DIGIT
	  int lr = u>wire*wirP_+wirD_ ? 1 : -1; // Left/right
	  ld.push_back(DCdig(wire,driftExact,*Ih,lr));
	}
      }
    }
    else {                   // ********** MULTI-HIT MC DECODING **********
      // (N.B: This is the recommanded decoding. It allows for a given MC track
      // to create several CsDigits, one for each drift cell traversed by the
      // track. Provided that:
      //  - path length through the cell is long enough: a built-in cut (path >
      //   fraction of thickness) is applied (path being the projection of the
      //   track's path in the (u,z) plane, which is an oversimplification),
      //  - double hit resolution is fulfilled (just as for the simplistic
      //   decoding supra).
      // These so many cells are determined from the info stored in the
      // CsMCTrkHit sub-class. )

      for (Ih = myMCHits_.begin(); Ih!=myMCHits_.end(); Ih++) {

	CsMCTrkHit *thit = dynamic_cast<CsMCTrkHit*>(*Ih);
	if( thit == 0 ) return;
	
	if (((thit->getMCTrack())->getParticle())->getCharge() ||
	    thit->getOrigin()) {
	  // From: Vadim.Alexakhine@cern.ch
	  // Date: 01/30/2003 11:26 PM
	  // I fixed some strange feature in makeMCdecoding for all detectors,
	  // namely hits from secondary charged particles, produced by
	  // originally neuteral particles (gammas, neutrons) via secondary
	  // interactions were skipped on the basis that original particle was
	  // neuteral.
	  // This gives some ~5% more clusters.
	  // Note(Y.B): Doesn't this also include neutral secondaries?
	  // ***** LOOP ON HITS from charged particles or charged products *****

	  double t  = thit->getDTime();            // Delay time (ns)
	  double ui = thit->getUin();    // Hit in point (DRS)
	  double vi = thit->getVin();
	  double wi = thit->getWin();
	  double uo = thit->getUout();   // Hit out point (DRS)
	  double vo = thit->getVout();
	  double wo = thit->getWout();

	  HepMatrix iRotM(3,3); iRotM = getRotWRSInv();
	  double xi = rotDRS_(1,1)*ui+rotDRS_(1,2)*vi+rotDRS_(1,3)*wi+
	    xcm_ - _xshift;
	  double yi = rotDRS_(2,1)*ui+rotDRS_(2,2)*vi+rotDRS_(2,3)*wi+
	    ycm_ - _yshift;
	  double zi = rotDRS_(3,1)*ui+rotDRS_(3,2)*vi+rotDRS_(3,3)*wi+
	    zcm_;
	  double xo = rotDRS_(1,1)*uo+rotDRS_(1,2)*vo+rotDRS_(1,3)*wo+
	    xcm_ - _xshift;
	  double yo = rotDRS_(2,1)*uo+rotDRS_(2,2)*vo+rotDRS_(2,3)*wo+
	    ycm_ - _yshift;
	  double zo = rotDRS_(3,1)*uo+rotDRS_(3,2)*vo+rotDRS_(3,3)*wo+
	    zcm_;
	  double Ui = iRotM(1,1)*xi+iRotM(1,2)*yi+iRotM(1,3)*zi; // WRS
	  //double Vi = iRotM(2,1)*xi+iRotM(2,2)*yi+iRotM(2,3)*zi;
	  double Wi = iRotM(3,1)*xi+iRotM(3,2)*yi+iRotM(3,3)*zi;
	  double Uo = iRotM(1,1)*xo+iRotM(1,2)*yo+iRotM(1,3)*zo; // WRS
	  //double Vo = iRotM(2,1)*xo+iRotM(2,2)*yo+iRotM(2,3)*zo;
	  double Wo = iRotM(3,1)*xo+iRotM(3,2)*yo+iRotM(3,3)*zo;

	  int wireF, wireL;        // First and las wire for this hit
	  if ((Ui-wirD_)/wirP_<0) wireF = int((Ui-wirD_)/wirP_-0.5);
	  else                    wireF = int((Ui-wirD_)/wirP_+0.5);
	  if ((Uo-wirD_)/wirP_<0) wireL = int((Uo-wirD_)/wirP_-0.5);
	  else                    wireL = int((Uo-wirD_)/wirP_+0.5);
	  int incr = wireL<wireF ? -1 : 1;

	  for (int i = wireF; i!=wireL+incr; i += incr) { 
            //LS Eff
            if( mcEffMapsEnabled_ && !MCEff_.empty() ){
              float mceff = MCEff_[i];
	      float random = (float) CsRandom::flat();
#ifdef DEBUG_EFF
	      cout<<"**** LS *** : CsDriftChamberDetector::makeMCDecoding : "<< GetTBName() <<" - "<< "MCEff_[" <<i<<"]= "<< mceff <<", random "<<random<<endl;
#endif
              if ( random > mceff ) break;
            }

	    if (i<0 || i>=nWir_) {
	      CsErrLog::msg(elAnomaly,__FILE__,__LINE__,
	      "%s: Wire(=%d) outside range [0,%d]",GetTBName().c_str(),i,nWir_);
	      continue;
	    }
	    double Uw = i*wirP_+wirD_, dU = Uo-Ui, dW = Wo-Wi;
	    double Ue, Us, We, Ws; // e/s = in/out
	    if (i==wireF && i==wireL) { // E=I, S=O
	      Ue = Ui; We = Wi; Us = Uo; Ws = Wo;
	    }
	    else if (i==wireF) {        // E=I, S=intersection w/ cell boundary
	      Ue = Ui; We = Wi; Us = Uw+incr*wirP_/2; Ws = (Us-Ui)*dW/dU+Wi;
	    }
	    else if (i==wireL) {        // E=intersection w/ cell boundary, S=O
	      Ue = Uw-incr*wirP_/2; We = (Ue-Ui)*dW/dU+Wi; Us = Uo; Ws = Wo;
	    }
	    else {
	      Ue = Uw-incr*wirP_/2; We = (Ue-Ui)*dW/dU+Wi;
	      Us = Uw+incr*wirP_/2; Ws = (Us-Ui)*dW/dU+Wi;
	    }
	    double ltrack = sqrt((Us-Ue)*(Us-Ue)+(Ws-We)*(Ws-We));
	    if (ltrack<zsiz_/20) // Efficiency: Require path length>thickness/20
	      continue;
	    double driftd = ((Uw-Ue)*(Ws-We)-(zcm_-We)*(Us-Ue))/ltrack;

	    // Drift time in ns: distance/velocity + TOF - trigger offset
	    double driftExact = fabs(driftd)/vel + t -
	      // Trigger time (and hence offset) is subtracted from hit times
	      triggerOffset;
	    // Is the hit inside the detector gate length?
	    if (driftExact>tGate_ || driftExact<0) continue;	  

	    //                                        ***** SAVE TEMPORARY DIGIT
	    int lr = // Left/right = sign of ordinate of mid-point of track path
	      (Ue+Us)/2>i*wirP_+wirD_ ? 1 : -1;
	    ld.push_back(DCdig(i,driftExact,*Ih,lr));
	  } // End of loop over wires
	} // End of hit origin check
      } // End of loop over hits
    }

    ld.sort();                  // ***** SORT TEMPORARY DIGITS BY WIRE# AND TIME

    // Here we've all what's needed for build digits...

    list<DCdig>::iterator id,idprev;
    for (id = ld.begin(); id!=ld.end(); id++) {
      //                                          ***** LOOP ON TEMPORARY DIGITS
      int wire = (*id).wire;

      //                  ***** DOUBLE HIT RESOLUTION: COMPARE w/ PREVIOUS DIGIT
      idprev = id; idprev--;
      if (id!=ld.begin() && 
	  wire==(*idprev).wire && (*id).time<(*idprev).time+dtres) continue;

      double drift; int dummyCounter = 0;          // ***** SMEAR THE DRIFT TIME
      do {
	drift = (*id).time + spSli_ * tiSli_ * CsRandom::gauss();
      } while((drift<0 || drift>tGate_) && dummyCounter++<200 );
      if (drift<0 || drift>tGate_) {
	CsErrLog::msg(elAnomaly,__FILE__,__LINE__,
        "%s: Drift time randomisation error: %.2f -> %.2f ns outside gate %.2f",
		      GetTBName().c_str(),(*id).time,drift,tGate_);
	continue;
      }

      //                       ***** MODULO "tiSli_" (here "drift" is always >0)
      drift = int( (drift)/tiSli_ + 0.5 ) * tiSli_;

      CsDigit *digit = new CsMCDigit(*this,wire,&drift);    // ***** NEW CsDigit
      CsMCDigit *mcDigit = dynamic_cast<CsMCDigit*>(digit);
      mcDigit->addHit(*((*id).ref)); mcDigit->setLR((*id).lr);
      myDigits_.push_back( digit );                      // ***** ADD IT TO LIST
    }   // End loop on temporary digits
  } // End of FZ file decoding

  else { // This is (or rather was: it may not work any longer) for old COMGEANT NTUPLES....
    if( Qhit.ndig == 0 && Qhit.nhit != 0 ) {
      // No way to make Digits from ntuple files...
      string str = "Decoding not possible on MC Ntuple and no Digits available on this file";
      CsErrLog::mes( elFatal, str );
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

      double drift = fabs( dig2 - t0_ ); // Is this correct??? Who knows...
      CsDigit* digit = new CsMCDigit( *this, dig1-1, &drift );
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

//____________________________________________________________________________
void CsDriftChamberDetector::DecodeChipDigit(const CS::Chip::Digit &digit)
{
  int ch;        // Channel (i.e. wire) number
  int whichSide; // For "split channels", as in straws: not useful in DCs
  double time;   // Time relative to trigger time

  //  ***** SWITCH BETWEEN SACLAY DCs and DC5 *****
  bool isDC5 = strncmp(GetTBName().c_str(),"DC05",4)==0; if (isDC5) {
    const CS::ChipSinica::Digit *d = dynamic_cast<const CS::ChipSinica::Digit *>(&digit);
    if (d==NULL)
      throw CS::Exception("CsDriftChamberDetector::DecodeRawData(): Wrong digit type");
    ch = d->GetWire();    time = d->GetTimeDecoded(); tiSli_ = d->GetTimeUnit();
    whichSide = d->GetWirePos();
  }
  else {                                                    //  ***** SACALY DCs
    const CS::ChipF1::Digit *d = dynamic_cast<const CS::ChipF1::Digit *>(&digit);
    if (d==NULL)
    throw CS::Exception("CsDriftChamberDetector::DecodeRawData(): Wrong digit type");
    ch = d->GetChannel(); time = d->GetTimeDecoded(); tiSli_ = d->GetTimeUnit();
    whichSide = d->GetChannelPos();
  }

  if( ch < 0 || ch >= int(getNWir()))
  throw CS::Exception("CsDriftChamberDetector::DecodeRawData() ==> Unexpected wire number %d for SiFi %d %s",
    ch,GetID().GetNumber(), GetTBName().c_str());
      
  // Skip if channel belongs to badChannels_ set
  if( badChannels_.find(ch) != badChannels_.end() ) return;

  double driftTime = time-getT0();    // Correct for t0

  vector<double> data; // Let's document CsDigits w/ some (possibly) useful info
  data.push_back(driftTime);
  data.push_back(whichSide);
  data.push_back(time);  // This is the raw, uncorrected, time 
  if (isDC5) { // Switching between Saclay DCs and DC5 again
    CS::ChipSinica::DataID id(digit.GetDataID());
    data.push_back(id.u.s.src_id);
    data.push_back(id.u.s.port);
    data.push_back(id.u.s.channel);
  }
  else {
    CS::ChipF1::DataID id(digit.GetDataID());
    data.push_back(id.u.s.src_id);
    // Geographical address: for debug mode we multiply it by '-1' factor.
    data.push_back(id.u.s.geoID_or_port * (id.u.s.mode?1:-1) );
    data.push_back(id.u.s.chip_chan);
  }  
  myDigits_.push_back(new CsDigit(*this,ch,&data[0],data.size()) );
}

//____________________________________________________________________________
double CsDriftChamberDetector::getDistToWire( const double t, bool &error ) {
  
  double r = 0;        // distance to wire
  error = false;
  	
  // check if RT
  if( hasRT_ ) r = rt_->getRfromT( t, error );
  else if( t < 0 || t > getTGate( ) ) { error = true; r = 0; } 
  else r = t*getVel()/getTiSli();
  
  if (r<0) { r=0; error = true; }
  else if (r>=0.5*wirP_) {
    r = 0.5*wirP_ - 0.001; // Subtract a small quantity so that hit from ith cell be always larger than hit from (i-1)th cell. Otherwise (hit,mirror) pairs may get intertwined.
    error = true;
  }
  
  return r;
}

//____________________________________________________________________________
double CsDriftChamberDetector::getDistToWire( const double t, const int wire, bool &error ) {
  
  double r = 0;        // distance to wire
  error = false;
  	
  // check if RT
  if( hasRT_ ) r = rt_->getRfromT( t, wire, error );
  else if( t < 0 || t > getTGate( ) ) { error = true; r = 0; } 
  else r = t*getVel()/getTiSli();
  
  double p1 = Pitch(double(wire-1));
  double p2 = Pitch(double(wire));
  double pmax = ( p1 >= p2 ? p1 : p2 );

  if (r<0) { r=0; error = true; }
  if (r>=0.5*pmax) {
    r = 0.5*pmax - 0.001; // Subtract a small quantity so that hit from ith cell be always larger than hit from (i-1)th cell. Otherwise (hit,mirror) pairs may get intertwined.
    error = true;
  }
  
  return r;
}

//____________________________________________________________________________
double CsDriftChamberDetector::getDRDT( const double t, bool &error ) {
  
  double v = 0;        
  error = false;
  
  // Check time is in time gate
  if( t < 0 || t > getTGate( ) ) { error = true; v = 0; }  	
  
  // Check if RT
  else if( hasRT_ ) v = rt_->getDRDT( t, error );
  else v = getVel()/getTiSli();
  return v;
}

//____________________________________________________________________________
void CsDriftChamberDetector::clusterize() {

  clearClusterList();

  if( isMWPC_ ) {

    //=== Process DCs as MWPC. ===
    list<CsDigit*>::iterator Id;
    list<CsDigit*> digits = getMyDigits();

    // protection
    if( digits.empty() ) return;

    HepMatrix iRotM(3,3); iRotM = getRotWRSInv();
    double wireDCorr = wirD_ + iRotM(1,1) * _deltax + iRotM(1,2) *_deltay; 

    for( Id=digits.begin(); Id!=digits.end(); Id++ ) {

      // skip if already clusterized
      if( (*Id)->isClusterized() ) continue;
      (*Id)->setClusterized();

      int wire = (*Id)->getAddress();
      double time = (*Id)->getDatum();     // true drift time (ns)
      
      //               ***** EVEN in MWPC mode, TIME HISTOGRAMS ARE FILLED *****
      if (hLevel_>=Normal) {
        H_t_->Fill(time); H_ch_->Fill(wire);
	if (hLevel_>=High) H_t_vs_ch_->Fill(wire,time);
      }

      //                                      ***** SET COORDINATES in WRS *****
      // u: perpendicular to the detector wires direction
      // v: parallel to the detector wires direction
      // w: = Z coordinate in MRS
      // origin: the origin of the MRS
      double u = wireDCorr + wire * wirP_, v = 0, w = zcm_;

      HepMatrix cov(3,3,0);                            // ***** SET ERRORS *****
      cov(1,1) = pow(wirP_/sqrt(12.0),2);  // wire pitch / sqrt(12)
      if (fabs(ang_-90)<1) cov(2,2) = pow(xsiz_/2,2); // wire length/2
      else                 cov(2,2) = pow(ysiz_/2/cos(ang_/180.*M_PI),2);
      cov(3,3) = 1.;  // assume 1 mm resolution in Z
      
      //                                            ***** SAVE THE CLUSTER *****
      CsCluster* cluster = new CsCluster(u,v,w,cov);
      cluster->addDigit(*(*Id)); cluster->addDet(*this); cluster->setTime(time);
      addCluster(*cluster);
    }

  }
  else {                // ********** PROCESS DCs USING DRIFT TIME **********
    int nCorDigits = 0;
    list<CsDigit*>::iterator Id;
    list<CsDigit*> digits = getMyDigits();

    if (digits.empty()) return;  // check number of digits

    HepMatrix iRotM(3,3); iRotM = getRotWRSInv();
    double wireDCorr = wirD_ + iRotM(1,1) * _deltax + iRotM(1,2) *_deltay; 
    double res = vel_ * spSli_;       // Detector resolution (mm)
    //double vel = vel_ / tiSli_;      // Drift velocity (mm/ns)

    bool isMC = CsEvent::Instance()->isAMonteCarloEvent();

    int prv_wire; static double prv_time; // "static" to avoid warning from compiler
    for (Id = digits.begin(), prv_wire = -1; Id!=digits.end(); Id++) {

      int wire = int((*Id)->getAddress()); // Fired wire
      double time = (*Id)->getDatum();     // Drift time (ns)

      if (hLevel_>=Normal) {                               // ***** BASIC HISTOS
        H_t_->Fill(time); H_ch_->Fill(wire);
	if (hLevel_>=High) H_t_vs_ch_->Fill(wire,time);// ***** HIGH-LEVEL HISTO
      }

      if (time<0 || time>tGate_) continue;  // ***** time W/IN TIME WINDOW *****

      // RETAIN ONLY THE FIRST DIGIT (within time gate)
      // NOTA BENE:
      // - Assumption is that digits are chronologically ordered
      //  (True in principle. I check it up in addition (Y.B))
      // - This may not be the ultimate solution. The idea is that later
      //  digits correspond mostly to retriggering of the ASD8 discriminator by
      //  the arrival of electon clusters belonging to the same track as
      //  the first one but originating from primary clusters away from
      //  mid-plane and as such non interesting (and hence to be discarded).
      //  But they may also arise from a second particle. There may be
      //  instances where it is possible to tell one from the other case:
      //  second comes long after first (long enough that all the signal
      //  from the first is exhausted).
      // - May not be valid for discontinued straws. But this should be
      //  marginal. There's no reason for the digit in the second part of the
      //  straw to come right after the digit in the first part (so that
      //  indeed "prv_wire" == "wire").
      if (deadTime_) {
	if (wire==prv_wire && fabs(time-prv_time)<deadTime_) continue;
	else {
	  prv_wire = wire; prv_time = time;
	}
      }

      //    ********** SET THE COORDINATES IN THE WRS **********
      // u: perpendicular to the detector wires direction
      // v: parallel to the detector wires direction
      // w: = Z coordinate in MRS
      // origin: the origin of the MRS

      bool error; double drift, u;  // ***** DRIFT TIME -> DISTANCE TO WIRE -> U
      if (IsVarPitch()) {
	drift = getDistToWire(time,wire,error); u = wireDCorr + Wire2Pos(wire);
      }
      else {
	drift = getDistToWire( time, error );   u = wireDCorr + wire * wirP_;
      }
      double v = 0, w = zcm_;
      //         ********** SET ERRORS ********** 
      HepMatrix cov(3,3,0);  // Zero matrix
      cov(1,1) = pow( res, 2 );    // Drift detector resolution 
      if (fabs(ang_-90)<1) cov(2,2) = pow(xsiz_/2,2);
      else                 cov(2,2) = pow(ysiz_/2/cos(ang_/180.*M_PI),2);
      cov(3,3) = 1.;               // Assume 1 mm resolution in Z

      if (hLevel_>=Normal) {             // ***** TIME-FILTERED HISTOGRAMS *****
        nCorDigits ++;
        H_t_cor_->Fill(time); H_ch_cor_->Fill(wire); H_r_->Fill(drift);
	H_ch_vs_r_->Fill(wire,drift);
      }

      if (drift!=0 || alwaysTwo_) {  // ***** IF drift!=0: SAVE 2 CLUSTERS *****
	if (drift==0) // +1 micron to ensure the ordering of clusters is stable
	  drift += .001;
        CsCluster *cluster1 = new CsCluster(u+drift,v,w,cov);
        cluster1->addDigit(*(*Id)); cluster1->addDet(*this);
        CsCluster *cluster2 = new CsCluster(u-drift,v,w,cov);
        cluster2->addDigit(*(*Id)); cluster2->addDet(*this);

	if (isMC) {
	  int lr = dynamic_cast<CsMCDigit*>(*Id)->getLR();
	  cluster1->isGenuine(lr>0); cluster2->isGenuine(lr<0);
	}

        cluster1->setMirrorCluster(*cluster2); cluster1->setLRProb( 0.5 );
        cluster2->setMirrorCluster(*cluster1); cluster2->setLRProb( 0.5 );
        cluster1->setTime(time);
        cluster2->setTime(time);
        addCluster(*cluster1); addCluster(*cluster2);
      }
      else {                                   // ***** JUST ONE OTHERWISE *****
        CsCluster *cluster = new CsCluster(u,v,w,cov);
        cluster->addDigit(*(*Id)); cluster->addDet(*this);
	if (isMC) cluster->isGenuine(true);
        cluster->setLRProb(1.0);                        // NO MIRROR!
        cluster->setTime(time);
        addCluster(*cluster);
      }
    }   // Loop over digits
    
    if (hLevel_>=Normal) H_n_cor_->Fill(nCorDigits);
  } // is_MWPC ?
  
  sortClusters(); setClusteringDone();
}

//____________________________________________________________________________
double CsDriftChamberDetector::getCorrU(const CsCluster *c,
					const double x, const double y, // MRS, mm
					const double tt,  // Offset w.r.t. *event* time
					bool &error)
{
  // *** This is a fast method to correct for propagation (and track latency):
  //   - All consistency checks are bypassed (on the ground that when
  //    you've reached such an advance step in the reconstruction that you
  //    want to correct for second order effects like propagation you
  //    must have already made sure that clusters are properly defined)
  //   - Resort to mirror to tell which sign to assign to the correction
  // (N.B.: Even faster woud have been to resort to the derivative of
  // the RT relation (instead of calculating R(T) twice, for t and t_cor)...
  // but it's not (yet as of 05/02) available)

  double Ui = c->getU();              // U initial
  double t; c->getTime(t);            // Uncorrected time
  double t_cor = t + dTdX_* (x-xcm_) + dTdY_*(y-ycm_) - tt;// Corrected time
  if (t_cor<0 || t_cor>tGate_) {      // Correction moved cluster out of gate..
    error = true; return Ui;          // ...return Ui
  }
  if (hasRT_) {
    CsCluster *mirror = c->getMirrorCluster();
    if (mirror) {
      double dR = rt_->getRfromT(t_cor,error);   // Corrected R
      dR -= getDistToWire(t,error);              // Delta R
      error = false;  // Disregard error: has to do w/ negative slope...
      if (Ui>mirror->getU()) return Ui+dR;
      else                   return Ui-dR;
    }
    else {      // No mirror => On wire! => Unable to tell which sign to assign
      error = false; return Ui;  // ... to the correction => Return U unchanged
    }
  }
  else {                           // No RT!...
    CsErrLog::msg(elFatal,__FILE__,__LINE__,"%s: No RT!",GetTBName().c_str());
    return 0; // To prevent "warning: control reaches end of non-void function"
  }
}

//____________________________________________________________________________
double CsDriftChamberDetector::getCorrT(const CsCluster *c,
					const double x, const double y, const double tt,
					bool &error)
{

  double t; c->getTime(t);            // Uncorrected time
  double t_cor = t + dTdX_* (x-xcm_) + dTdY_*(y-ycm_) - tt;// Corrected time
  if (t_cor<0 || t_cor>tGate_) {      // Correction moved cluster out of gate..
    error = true; return 1.e6;        // ...return anything
  } else {
    error = false;
    return t_cor;
  }

}

//____________________________________________________________________________
int CsDriftChamberDetector::_getFEOrientation( void )
{
  char tbn[] = "DC01X1__"; strncpy(tbn,GetTBName().c_str(),8);
  bool error = false;
  if      (strncmp(tbn,"DC00",4)==0 || strncmp(tbn,"DC01",4)==0) {
    if      (tbn[5]=='1' &&  tbn[4]=='Y')
      return -1;
    else if (tbn[5]=='1' && (tbn[4]=='X' || tbn[4]=='U' || tbn[4]=='V'))
      return  1;
    else if (tbn[5]=='2' &&  tbn[4]=='Y')
      return 1;
    else if (tbn[5]=='2' && (tbn[4]=='X' || tbn[4]=='U' || tbn[4]=='V'))
      return -1;
    else error = true;
  }
  else if (strncmp(tbn,"DC02",4)==0 || strncmp(tbn,"DC03",4)==0 ||
	   strncmp(tbn,"DC04",4)==0 || strncmp(tbn,"DC05",4)==0) {
    if      (tbn[5]=='1') return  1;  //!< (v>0 => delta t >0 <=> Front-End is at v>0 )
    else if (tbn[5]=='2') return -1;  //!< (v>0 => delta t <0 <=> Front-End is at v<0  )  
    else error = true;
  }
  else {
    CsErrLog::msg(elFatal,__FILE__,__LINE__,
		  "%s: Not catered for => No propagation time correction.",
		  GetTBName().c_str());
    return 0;
  }
  if (error)
    CsErrLog::msg(elFatal,__FILE__,__LINE__,
		  "%s: Unexpected TB name",GetTBName().c_str() );
  return 0; // Control never reaches there: just to fool compiler
}   

//____________________________________________________________________________
void CsDriftChamberDetector::updateClusters(double time)
{
  HepMatrix iRotM(3,3); iRotM = getRotWRSInv();
  double wireDCorr = wirD_ + iRotM(1,1) * _deltax + iRotM(1,2) *_deltay; 

  bool error;

  list<CsCluster*>::iterator ic; CsCluster *prv;
  for (ic = myClusters_.begin(), prv = 0; ic!=myClusters_.end(); ic++) {
    CsCluster *c = *ic, *m = c->hasMirrorCluster() ? c->getMirrorCluster() : 0;
    if (!m)     // No mirror. May happen if hit precisely on wire? not sure...
      continue; // ...Anyway nothing we can do = skip
    if (m==prv) // Mirror of previous. (Note that if there are several pairs
      // (hit,mirror) per cell, hit and mirror may not follow each other. Will
      // then be processed twice. Which doesn't matter, but for CPU time.) ...
      continue; // ...Already processed => skip
    prv = c;
    double t; c->getTime(t); t -= time;
    int wire = c->getDigitsList().front()->getAddress();
    // ***** DISTANCE TO WIRE (Error is ignored just as in "clusterize") *****
    double drift, u; if (IsVarPitch()) {
      drift = getDistToWire(t,wire,error); u = wireDCorr + Wire2Pos(wire);
    }
    else {
      drift = getDistToWire(t,     error); u = wireDCorr + wire * wirP_;
    }
    if (drift==0)// +1 micron to keep the ordering of (cluster,mirror) unchanged
      drift += .001;
    if (m->getU()>=c->getU()) {   // This is what's expected...
      // ...given that the clusters may have been ordered by now.
      // (Note that the two cluster/mirror can be strictly equal. This, despite
      // the fact that in the "clusterize" method, we care to assign them two
      // distinct U coordinates (differing by at least 0.001, cf. the
      // "alwaysTwo_" block). Because the two can become equal again, over the
      // course of the track reconstruction (cf. the "coalesce" block in
      // "TEv::ImportClusters") and then be used again, in a 2nd attempt at
      // event reconstruction (cf. "CsEvent::_reconstructionSchema001", flag
      // "_reTrackingON"; 2nd attempt when this very "updateClusters" is
      // called upon.)
      c->setU(u-drift); m->setU(u+drift);
    }
    else {                   // ...nevertheless provide also for the alternative
      c->setU(u+drift); m->setU(u-drift);
    }
  }
}

//____________________________________________________________________________
bool CsDriftChamberDetector::_readRTRelation( void )
{
  
  CsRTRelation *rt = new CsRTRelation( *this );
  if( rt->readFromOpt() ) {
    rt_ = rt;
    hasRT_ = true;
  } else {
    delete rt;    
    hasRT_ = false;  
  }
  
  return hasRT_;
}   

//=============================================
//=== operator for reading Calibration Data ===
//=============================================
istream& operator>>(istream& in,CsDriftChamberDetector::CalibrationData &c) {
  in >> c.t0;     // read t0
  in >> *(c.rt);  // read RT Relation
  
  // check Key for Dead Wires
  if( in.eof() ) return in;
  string key;
  in >> key;
  if( in.fail() || key != "BADCHANNELS" ){ in.clear(); return in; }

  // Read Dead Wires
  while( in.good() && ! in.eof() ){
    unsigned int bc;
    in >> bc;
    if( in.good() ) c.badChannels.insert( bc );
  }

  if( in.eof() ) in.clear();
  return in;
}

//------------------------------------------------------------------------------

void CsDriftChamberDetector::readCalibration(time_t timePoint){
  tm *t = localtime(&timePoint);
  CalibrationData* calib = 0;

  // read-in corresponding calibration constants
  CDB::Time tp(timePoint,0);
  string data;
  cdb_->read(GetTBName(),data,tp);
  if ((getenv("COND_DB_DEBUG")!=0) || (getenv("COND_DB_DEBUG_DC")!=0)) {
    cout <<GetTBName()<<endl;
    cout << data << endl;
  }

  istringstream is(data);
  calib = new CalibrationData(*this);
  is>>(*calib);

  // update t0:
  t0_ = calib->t0;
  CsErrLog::msg(elInfo,__FILE__,__LINE__,
		"%s: New T0 %f",GetTBName().c_str(),t0_);
  
  // check dead wires
  if( calib->badChannels.size() ){
    badChannels_ = calib->badChannels;
    CsErrLog::msg(elInfo,__FILE__,__LINE__,
		  "%s: %d bad chanel%s.",GetTBName().c_str(),
		  badChannels_.size(),
		  badChannels_.size()>1 ? "s" : "");
  }

  // check rtRelation
  // Option file RTRelations have priority over DB, for debugging
  if( hasRT_ )
    CsErrLog::msg(elWarning,__FILE__,__LINE__,
		  "%s already has RT (prob. from .opt). DB is skipped.",
		  GetTBName().c_str());
  else if (calib->rt->hasGrid() || calib->rt->hasRegularGrid() || calib->rt->hasRTParameters()) {
    rt_ = calib->rt;
    hasRT_ = true;
    CsErrLog::msg(elInfo,__FILE__,__LINE__,
		  "%s: calibration updated.",GetTBName().c_str());
  } else {
    delete calib;
    hasRT_ = false;
    cout << "CsDriftChamberDetector::readCalibration() ==> " << GetTBName() << " RT relation from database not valid\n";
  }

}

//LS Eff
void CsDriftChamberDetector::readMCEffMaps(time_t timePoint)
{
  if( mccdb_ == NULL ) return;
  CDB::Time tp(timePoint,0);
  tm *t = localtime(&tp.first);
  string data("");
  mccdb_->read(GetTBName(),data,tp);
  istringstream is(data);
  string line;
  while( getline(is,line) )
    {
      if( line.length()==0 || line[0]=='#' )
        continue;
      int c;
      float mcEf, mcEf_err;
      if( 3!=sscanf(line.c_str(),"%d %g %g",&c,&mcEf,&mcEf_err) )
	throw std::logic_error("CsDriftChamberDetector_MCEff::read(): bad line structure.");
      MCEff_.push_back(mcEf);
      MCEff_err_.push_back(mcEf_err);
    }
}
