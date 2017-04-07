// $Id: CsRichWallDetector.cc,v 1.19 2010/01/24 16:10:41 suhl Exp $                              

/*!
   \file    CsRichWallDetector.cc
   \brief   Compass RICH wall detector Class.
   \author  Benigno Gobbo
   \version $Revision: 1.19 $
   \date    $Date: 2010/01/24 16:10:41 $
*/


#include "CsRichWallDetector.h"

#include <math.h>
#include <string.h>
#include <strings.h>
#include <functional>
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
#include "CsRTRelation.h"
#include "CDB.h"
//LS Eff
#include <stdexcept>
#define DEBUG_EFF

#include <iostream>
#include <iomanip>
#include <stdio.h>
#include <ctime>

using namespace std;
using namespace CLHEP;

extern QhitType Qhit;

const float CsRichWallDetector::F1bin_ = 0.12892; // ns 

CsRichWallDetector::CsRichWallDetector( const int    row,
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
  propVel_inv_( 0.0033 ),   //Unit is ns/mm
  deadTime_(0.)
{
  CsOpt *opt = CsOpt::Instance( );

  decodeCard_ = false;                            // ***** DECODING METHOD *****
  bool status = opt->getOpt( "" , "make decoding" );
  if (status) {
    list<string> options; list<string>::iterator Is;
    if ((status = opt->getOpt("","make decoding",options))) {
      for (Is = options.begin(); Is!=options.end(); Is++ ) 
	if (*Is=="DR" || *Is=="all") decodeCard_ = true;
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
  //                                          ***** eff. MAP ENABLED ? *****                    
  mcEffMapsEnabled_ = false;
  if (opt->getOpt("ALL"      ,            "mcEffMapsEnabled") ||
      opt->getOpt(GetTBName(),            "mcEffMapsEnabled") ||
      opt->getOpt(GetTBName().substr(0,2),"mcEffMapsEnabled")) {
    mcEffMapsEnabled_ = true;
#ifdef DEBUG_EFF
    cout<<"**** LS *** : CsRichWallDetector:: "<< mcEffMapsEnabled_<<endl;  
#endif
  }
  else CsErrLog::msg(elWarning,__FILE__,__LINE__,
                     "%s: eff. Map will NOT be used",TBname);

  //                                         ***** DR TREATED LIKE a MWPC? *****
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
  alwaysTwo_ = opt->getOpt("DR","make always two clusters");

  //   ********** INITIALISATION for PROPAGATION TIME CORRECTION **********
  //                                            ***** PROPAGATION VELOCITY *****
  if (opt->getOpt(TBname,                 "PropVelInv",propVel_inv_) ||
      opt->getOpt(GetTBName().substr(0,2),"PropVelInv",propVel_inv_)) {
    CsErrLog::msg(elWarning,__FILE__,__LINE__,
		  "%s: inverse Propagation Velocity = %f",TBname,propVel_inv_);
  } 
  else if (CsInit::Instance()->IsAMonteCarloJob())        // ***** MC CASE...
	   propVel_inv_ = 0; // ...So far, as of 12/06, propagation is NOT simulated in coral.
	   

  // ***** CLUSTERISATION OPTION *****

  // put here detector specific clusterization options   

  // Temporary: only needed to speed up the response of "Wire2Pos", which, temporarily, differs in MC and RD. Temporary, but yet still needed as of 2009/08.
  isMC_ = CsInit::Instance()->IsAMonteCarloJob();

  // Coordinate correction reading

  if (opt->getOpt(GetTBName(),            "Shift_corr",shift_corr_path) ||
      opt->getOpt(GetTBName().substr(0,2),"Shift_corr",shift_corr_path)) {

    CsErrLog::msg(elWarning,__FILE__,__LINE__,
		  "%s: inverse Propagation Velocity = %s",TBname,shift_corr_path.c_str());
    
    string file_name = GetTBName().substr(0,6) ;

    cout << " DR FILE NAME - " << file_name << endl << flush;

    shift_corr_path += file_name + ".shift";

    ifstream in_corr_file(shift_corr_path.c_str());

    double curr_coor_val = 0;

    while(in_corr_file >> curr_coor_val){

      //    for(int i =0; i < 74 ; i++){

      //      curr_coor_val = 0;

      //   in_corr_file >> curr_coor_val;

      //      tube_shift[i] = curr_coor_val;

      cout << "curr_coor_val - " << curr_coor_val << endl;

      tube_shift.push_back(curr_coor_val);

    }

    in_corr_file.close();
    
  }

}

//____________________________________________________________________________
bool CsRichWallDetector::operator==( const CsRichWallDetector& det ) const {
  return( CsDetector::operator==(det) );
}

//____________________________________________________________________________
bool CsRichWallDetector::operator<( const CsRichWallDetector& det ) const {
  return( CsDetector::operator<(det) );
}

//____________________________________________________________________________
void CsRichWallDetector::BookHistograms() {

  //=== check if histograms are to be booked ===
  CsDetector::ReadHistLevel();
  if( hLevel_ == None ) return;

  CsHistograms::SetCurrentPath("/CsDriftTube");  

  if( hLevel_ >= Normal ) {
    char name[100];
    sprintf( name, "%s_t_", GetTBName().c_str() );
    H_t_ = new CsHist1D( name, name, (int) (3*tGate_) , -tGate_, 2*tGate_ );

    sprintf( name, "%s_ch_", GetTBName().c_str() );
    H_ch_ = new CsHist1D( name, name, nWir_, 0, nWir_ );
  }

  CsHistograms::SetCurrentPath("/");  
  return; 
}

// Local class for tempopary digits in the detector response simulation 
// (Used in CsDriftChamberDetector::makeMCDecoding())
class RWdig
{
public:
  RWdig(int w, double t, CsMCHit *r, int l) {
    wire = w; time = t; ref = r; lr = l;
  }
  int    wire;  // Wire #
  double time;  // Drift time
  CsMCHit *ref; // Reference to MCHit 
  int    lr;    // Left/right
  bool operator < (const RWdig &gd) const 
  { 
    if(wire == gd.wire) return (time < gd.time);
    else                return (wire < gd.wire);
  };
};

//____________________________________________________________________________
void CsRichWallDetector::makeMCDecoding() {

  if (!decode_ && !decodeCard_) return;   // Should I proceed?
  if (decodingDone_) return;              // Already done?
  myDigits_.clear();                      // Clear

  double vel = vel_ / tiSli_;           // cast velocity in mm/ns
  double dtres = tiSli_ * thRes_;       // double hit resolution, ns  

  // The following can be done only with zebra (as opposed to NTuple) files... 
  if (!CsGeant3::Instance()->isAnNtFile()) {

    list<RWdig> ld; // ***** LIST OF TEMPORARY DIGITS... *****
    // (...i.e. digits before double hit resolution is applied.) 
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

	  int wire = Pos2Wire(u-wirD_,v);      // Nearest wire
	  if (wire<0 || wire>=nWir_) {
	    CsErrLog::msg(elAnomaly,__FILE__,__LINE__,
	   "%s: Wire(=%d) outside range [0,%d[",GetTBName().c_str(),wire,nWir_);
	    continue;
	  }

	  // Drift time in ns: distance/velocity + TOF - trigger offset
	  double driftExact = fabs(u-Wire2Pos(wire)-wirD_)/vel + t -
	    // Trigger time (and hence offset) is subtracted from hit times
	    triggerOffset;

	  // Is the hit inside the detector gate length?
	  if (driftExact>tGate_ || driftExact<0) continue;

	  //                                          ***** SAVE TEMPORARY DIGIT
	  int lr = u>Wire2Pos(wire)+wirD_ ? 1 : -1; // Left/right  //15.03.09
// 	  int lr = u>wire*wirP_+wirD_ ? 1 : -1; // Left/right
	  ld.push_back(RWdig(wire,driftExact,*Ih,lr));
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
      // CsMCTrkHit sub-class.
      // !!! WARNING !!! The present version of the code does not take into
      // account the irreqular spacing of the RICHwall detectors.)
      for (Ih = myMCHits_.begin(); Ih!=myMCHits_.end(); Ih++) {

	CsMCTrkHit *thit = dynamic_cast<CsMCTrkHit*>(*Ih);
	if( thit == 0 ) return;
	
	if (((thit->getMCTrack())->getParticle())->getCharge() ||
	    thit->getOrigin()) {

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
	  double Vi = iRotM(2,1)*xi+iRotM(2,2)*yi+iRotM(2,3)*zi;
	  double Wi = iRotM(3,1)*xi+iRotM(3,2)*yi+iRotM(3,3)*zi;
	  double Uo = iRotM(1,1)*xo+iRotM(1,2)*yo+iRotM(1,3)*zo; // WRS
	  double Vo = iRotM(2,1)*xo+iRotM(2,2)*yo+iRotM(2,3)*zo;
	  double Wo = iRotM(3,1)*xo+iRotM(3,2)*yo+iRotM(3,3)*zo;

	  int wireF = Pos2Wire(Ui-wirD_,Vi); // First and last wire for this hit
	  int wireL = Pos2Wire(Uo-wirD_,Vo);
	  int incr = wireL<wireF ? -1 : 1;

	  for (int i = wireF; i!=wireL+incr; i += incr) { 
	    
            //LS Eff
            if( mcEffMapsEnabled_ && !MCEff_.empty() ){
              float mceff = MCEff_[i];
	      float random = (float) CsRandom::flat();
#ifdef DEBUG_EFF
              cout<<"**** LS *** : CsRichWallDetector::makeMCDecoding : "<< GetTBName() <<" - "<<
		"MCEff_[" <<i<<"]= "<< mceff <<", random "<<random<<endl;
#endif
              if ( random > mceff ) break;
            }

	    if (i<0 || i>=nWir_) {
	      CsErrLog::msg(elAnomaly,__FILE__,__LINE__,
	      "%s: Wire(=%d) outside range [0,%d]",GetTBName().c_str(),i,nWir_);
	      continue;
	    }
	    double Uw = Wire2Pos(i)+wirD_, dU = Uo-Ui, dW = Wo-Wi;
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
	    int lr = (Ue+Us)/2>Wire2Pos(i)+wirD_ ? 1 : -1; //15.03.09
// 	    int lr = // Left/right = sign of ordinate of mid-point of track path
// 	      (Ue+Us)/2>i*wirP_+wirD_ ? 1 : -1;
	    ld.push_back(RWdig(i,driftExact,*Ih,lr));
	  } // End of loop over wires
	} // End of hit origin check
      } // End of loop over hits
    }

    ld.sort();                  // ***** SORT TEMPORARY DIGITS BY WIRE# AND TIME

    // Here we've all what's needed for build digits...

    list<RWdig>::iterator id,idprev;
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
void CsRichWallDetector::DecodeChipDigit(const CS::Chip::Digit &digit)
{

  const CS::ChipF1::Digit *d = dynamic_cast<const CS::ChipF1::Digit *>(&digit);
  if( d==NULL )
    throw CS::Exception("CsRichWallDetector::DecodeRawData(): Wrong digit type");

  int ch = d->GetChannel();
  if( ch < 0 || ch >= int(getNWir()))
  throw CS::Exception("CsRichWallDetector::DecodeRawData() ==> Unexpected wire number %d for RW %d %s",
    d->GetChannel(),GetID().GetNumber(), GetTBName().c_str());
     

  // ***** TIME RELATIVE to TRIGGER *****
  double time = d->GetTimeDecoded();

  //  if (time<hittMin_ || time>hittMax_) return;  // ***** CUT on Hit Time *****

  time -= getT0();            // correct from time offset

  //  cout << setprecision(10) << "Time = " << d->GetTimeDecoded() << " ,T0 = " << getT0() << endl;

  tiSli_ = d->GetTimeUnit();

  if ( time < 0 || time > getTGate() ) return;

  vector<double> data;
  data.push_back(time);
  data.push_back(d->GetChannelPos());
  data.push_back(d->GetTimeDecoded());

  myDigits_.push_back( new CsDigit(*this, ch, &data[0], 3) );
}

//____________________________________________________________________________
double CsRichWallDetector::getDistToWire( const double t, bool &error ) {
  
  double r = 0;        // distance to wire
  error = false;

  // check if RT
  if( hasRT_ ) r = rt_->getRfromT( t, error );
  else if( t < 0 || t > getTGate() ) { error = true; r = 0; } 
  else r = t*getVel()/getTiSli();
  
  if (r<0) { r=0; error = true; }
  else if (r>=0.5*wirP_) {
    r = 0.5*wirP_ - 0.001; // Subtract a small quantity so that hit from ith cell be always larger than hit from (i-1)th cell. Otherwise (hit,mirror) pairs may get intertwined.
    error = true;
  }

  return r;
}

double CsRichWallDetector::getDRDT( const double t, bool &error ) {
  
  double v = 0;        
  error = false;
  
  // Check time is in time gate
  if( t < 0 || t > getTGate( ) ) { error = true; v = 0; }  	
  
  // Check if RT
  else if( hasRT_ ) v = rt_->getDRDT( t, error );
  else v = getVel()/getTiSli();
  return v;
}

double CsRichWallDetector::Wire2Pos( double wire ) const {
  // Converts distance to first wire from wire units to mm.
  // (Would there be a uniform pitch <P>, it would return: wire*P.)

  if (isMC_) return wire*wirP_;

  double wire_tmp = 0;
  double pos = 0;

  if ( nWir_ == 592 ){ //x coord

    if(wire < 496){
      pos = wire * wirP_ + (wirP_/2)*(int)(wire/8); // *ModifComm*: CsDetector's "Wire2Pos" does not return the position corrected for "wirD_", but the distance to 1st wire. Let's do the same here, in order to avoid confusion.
    }else{
      wire_tmp = wire - 296;
      pos = wire_tmp * wirP_ + (wirP_/2)*(int)(wire_tmp/8); 
    }

  }else{ //y coord

    if(wire < 368){
      pos = wire * wirP_ + (wirP_/2)*(int)(wire/8);
    }else{
      wire_tmp = wire - 208;
      pos = wire_tmp * wirP_ + (wirP_/2)*(int)(wire_tmp/8);
    }
    
  }

  return pos;

}

double CsRichWallDetector::Wire2Pos_cor( double wire ) const {
  
  double wire_tmp = 0;
  double pos = 0;
  
  if ( nWir_ == 592 ){ //x coord
    
    if(wire < 496){
      pos = wire * wirP_ + (wirP_/2)*(int)(wire/8) + tube_shift[(int)(wire/8)];
    }else{
      wire_tmp = wire - 296;
      pos = wire_tmp * wirP_ + (wirP_/2)*(int)(wire_tmp/8) + tube_shift[(int)(wire/8)]; 
    }
    
  }else{ //y coord
    
    if(wire < 368){
      pos = wire * wirP_ + (wirP_/2)*(int)(wire/8) + tube_shift[(int)(wire/8)];
    }else{
      wire_tmp = wire - 208;
      pos = wire_tmp * wirP_ + (wirP_/2)*(int)(wire_tmp/8) + tube_shift[(int)(wire/8)];
    }
    
  }
  
  return pos;
  
}

int CsRichWallDetector::Pos2Wire( double pos , double perp_koord) {
  // Converts distance to first wire to wire number
  // (Would there be a uniform pitch <P>, it would return nearest integer to
  // pos/P.)

  // This function is only useful in the MC case.

  // *ModifComm*: The v1.5 version does not seem to be OK and is in any case awkward => I decide to revert, for the time being, to the simplistic scheme that ignores the superpitch intricacy. Of course, the "Wire2Pos" function has to have a special processing, corresponding to the simplistic scheme, for MC data.
   if (isMC_) {
     if (pos/wirP_<0) return int(pos/wirP_-0.5);
     else             return int(pos/wirP_+0.5);
   }

//   double first_dist = 10 , second_dist = 5;
//   double wire_tmp = 0;
//   double coordinate = 0;

//   int found_wire = -1;

  int tube = 0, rec_wire = 0;

  if ( nWir_ == 592 ){

    if ( pos < -5. )

      return -1;

    else if (pos < 0.){

      return 0;
      
    }else if (pos > (getXsiz()- 5.) ){

      return 600;

    }else{

      tube = (int)((pos + 7.5)/85) ;
      rec_wire = (int)(pos +5  - tube*5 )/10 ;
      
    }

    if (DetectorSection(rec_wire) == 1 and perp_koord < 0)
      return rec_wire;
    else if (DetectorSection(rec_wire) == 1 and perp_koord > 0)
      return rec_wire + 296;
    else
      return rec_wire;

  }else{

    perp_koord = - perp_koord;

    if ( pos < -5. )

      return -1;

    else if (pos < 0.){

      return 0;
      
    }else if (pos > (getYsiz()- 5.) ){

      return 600;

    }else{

      tube = (int)((pos + 7.5)/85) ;
      rec_wire = (int)(pos +5  - tube*5 )/10 ;
      
    }

    if (DetectorSection(rec_wire) == 1 and perp_koord < 0)
      return rec_wire;
    else if (DetectorSection(rec_wire) == 1 and perp_koord > 0)
      return rec_wire + 208;
    else
      return rec_wire;

  }

}

int CsRichWallDetector::Pos2Wire_cor( double pos , double perp_koord) {
  
  double first_dist = 0 , second_dist = 0;
  double coordinate = 0;
  double coordinate_p = 0;
  
  int found_wire = -1, wire = -1, wire_p = -1;
  
  if ( nWir_ == 592 ){
    
    for(int i = 0; i < nWir_ ; i++){
      
      coordinate = 0;
      coordinate_p = 0;
      
      wire = 0;
      wire_p = 0;
      
      first_dist = 0;
      second_dist = 0;

      if ( DetectorSection(i) == 1 && perp_koord > 0 ){
	
	wire = i + 296;
	
      }else{
	
	wire = i;
	
      }
      
      if(i == 0){

	wire_p = i;
	
      }else if( DetectorSection(i-1) == 1 && perp_koord > 0 ){
	
	wire_p = i -1 + 296;
	
      }else{
	
	wire_p = i -1;
	
      }
      
      coordinate = Wire2Pos_cor(wire);
      coordinate_p = Wire2Pos_cor(wire_p);
      
      if(wire == wire_p){
	
	if(coordinate >= pos){
	  
	  if(fabs(coordinate - pos) > 5.0){
	    found_wire = -1;
	    break;
	  }else{
	    found_wire = wire;
	    break;
	  }
	  
	}
	
      }else if(i == 495 && coordinate <= pos){
	
	if( fabs(coordinate - pos) > 5.0){
	  found_wire = -1;
	  break;
	}else{
	  found_wire = wire;
	  break;
	}
	
      }else{
	
	if(coordinate >= pos){
	  
	  first_dist = fabs(coordinate - pos);
	  second_dist = fabs(coordinate_p - pos);
	  
	  if(first_dist < second_dist){
	    found_wire = wire;
	    break;
	  }else if(second_dist < first_dist){
	    found_wire = wire_p;
	    break;
	  }else{
	    
	    if(DetectorSection(i) != DetectorSection(i-1)){
	      
	      if(DetectorSection(i) == 1){
		found_wire = wire_p;
		break;
	      }else{
		found_wire = wire;
		break;
	      }
	      
	    } else {
	      
	      srand(time(0));
	      
	      if( ((double)rand()/(double)RAND_MAX) > 0.5){
	        found_wire = wire;
	        break;
	      }else{
	        found_wire = wire_p;
	        break;
	      }
	      
	    }
	    
	  }
	  
	}
	
      }
      
    }
    
  }else{
    
    perp_koord = -perp_koord;

    for(int i = 0; i < nWir_; i++){
      
      coordinate = 0;
      coordinate_p = 0;
      
      wire = 0;
      wire_p = 0;
      
      first_dist = 0;
      second_dist = 0;
      
      if ( DetectorSection(i) == 1 && perp_koord > 0 ){
	
	wire = i + 208;
	
      }else{
	
	wire = i;
	
      }
      
      if(i == 0){
	
	wire_p = i;
	
      }else if( DetectorSection(i-1) == 1 && perp_koord > 0 ){
	
	wire_p = i -1 + 208;
	
      }else{
	
	wire_p = i -1;
	
      }
      
      coordinate = Wire2Pos_cor(wire);
      coordinate_p = Wire2Pos_cor(wire_p);
      
      if(wire == wire_p){
	
	if(coordinate >= pos){
	  
	  if(fabs(coordinate - pos) > 5.0){
	    found_wire = -1;
	    break;
	  }else{
	    found_wire = wire;
	    break;
	  }
	  
	}
	
      }else if(i == 367 && coordinate <= pos){
	
	if( fabs(coordinate - pos) > 5.0){
	  found_wire = -1;
	  break;
	}else{
	  found_wire = wire;
	  break;
	}
	
      }else{
	
	if(coordinate >= pos){
	  
	  first_dist = fabs(coordinate - pos);
	  second_dist = fabs(coordinate_p - pos);
	  
	  if(first_dist < second_dist){
	    found_wire = wire;
	    break;
	  }else if(second_dist < first_dist){
	    found_wire = wire_p;
	    break;
	  }else{
	    
	    if(DetectorSection(i) != DetectorSection(i-1)){
	      
	      if(DetectorSection(i) == 1){
		found_wire = wire_p;
		break;
	      }else{
		found_wire = wire;
		break;
	      }
	      
	    } else {
	      
	      srand(time(0));
	      
	      if( ((double)rand()/(double)RAND_MAX) > 0.5){
	        found_wire = wire;
	        break;
	      }else{
	        found_wire = wire_p;
	        break;
	      }
	      
	    }
	    
	  }
	  
	}
	
      }
      
    }
    
  }
  
  return found_wire;

}

vector<int> CsRichWallDetector::WiresCloseTo( int range, int current_wire, double h_V ) {
  
  vector<int> out_wires;
  vector<int>::iterator out_w_it;
  
  int chann_range = (2*range + 1);
  int current_det_section = DetectorSection(current_wire);

  if(getNWir() != 592){
    h_V = -h_V; 
  }

  // BEGIN channels check
  
  for(int i=0 ; i < chann_range ; i++){
    out_wires.push_back(-1);
  }
  
  if (current_wire > -1 && current_wire < getNWir() ){
    
    out_wires[range] = current_wire;
    
    for(int j = 1 ; j < range+1 ; j++){  // BEGIN , SMALLER CHANNELS
      
      int current_wire_tmp = current_wire - j;
      
      int current_det_section_tmp = DetectorSection(current_wire_tmp);
      
      if(current_det_section_tmp != current_det_section){
	
	if(current_det_section == 3){
	  
	  if(getNWir() == 592){
	    current_wire_tmp = current_wire_tmp - 296;
	  }else{
	    current_wire_tmp = current_wire_tmp - 208; 
	  }
	  
	  out_wires[range-j] = current_wire_tmp;
	  
	}else if(current_det_section == 2){
	  
	  if(getNWir() == 592){
	    if(h_V > 0){
		current_wire_tmp = current_wire_tmp + 201;
	    }
	  }else{
	    if(h_V > 0){
	      current_wire_tmp = current_wire_tmp + 161;
	    }
	  }
	  
	  out_wires[range-j] = current_wire_tmp;
	  
	}else{
	  
	  out_wires[range-j] = current_wire_tmp;
	  
	}
	
	
      } else {
	
	if (current_wire_tmp > -1){
	  
	  out_wires[range-j] = current_wire_tmp;
	  
	}else{
	  
	  out_wires[range-j] = -2;
	  
	}
	
      }
      
      
    } // END , SMALLER CHANNELS
    
    
    for(int j = 1 ; j < range+1 ; j++){  // BEGIN , BIGGER CHANNELS
      
      int current_wire_tmp = current_wire + j;
      
      int current_det_section_tmp = DetectorSection(current_wire_tmp);
      
      if(current_det_section_tmp != current_det_section){
	
	if(current_det_section == 2){
	  
	  out_wires[range+j] = -2;
	  
	}else if(current_det_section == 0){
	  
	  if(getNWir() == 592){
	    if(h_V > 0){
	      current_wire_tmp = current_wire_tmp + 296;
	    }
	  }else{
	    if(h_V > 0){
	      current_wire_tmp = current_wire_tmp + 208;
	    }
	  }
	  
	  out_wires[range+j] = current_wire_tmp;
	  
	}else{
	  
	  out_wires[range+j] = current_wire_tmp;
	  
	}
	
	} else {
	
	if (current_wire_tmp < getNWir() ){
	  
	  out_wires[range+j] = current_wire_tmp;
	  
	}else{
	  
	  out_wires[range+j] = -2;
	  
	}
	
      }
      
      
    } // END , BIGGER CHANNELS
    
    
  }
    
  // END channels check
  
  return out_wires;
  
}


int CsRichWallDetector::DetectorSection( int wire ) {

  if ( nWir_ == 592 ){ //x coord
   
    if(wire < 200)
      return 0;
    else if(wire >= 200 && wire <= 295)
      return 1;
    else if(wire > 295 && wire < 496)
      return 2;
    else if(wire >= 496)
      return 3;

  }else{ //y coord

    if(wire < 160)
      return 0; 
    else if(wire >= 160 && wire <= 207)
      return 1;
    else if(wire > 207 && wire < 368)
      return 2;
    else if(wire >= 368)
      return 3;
  }

  return -1;
  
}

void CsRichWallDetector::clusterize() {

  clearClusterList();

  if( isMWPC_ ) {

    //=== Process Drift Tubes as MWPC (no drift time). === 
    list<CsDigit*>::iterator Id;
    list<CsDigit*> digits = getMyDigits();

    // protection
    if( digits.empty() ) return;

    HepMatrix iRotM(3,3); iRotM = getRotWRSInv();
    double wireDCorr = wirD_ + iRotM(1,1) * _deltax + iRotM(1,2) *_deltay; 

    int prv_wire; static double prv_time; // "static" to avoid warning from compiler
    for (Id = digits.begin(), prv_wire = -1; Id!=digits.end(); Id++) {

      // skip if already clusterized
      if( (*Id)->isClusterized() ) continue;
      (*Id)->setClusterized();

      int wire = (*Id)->getAddress();
      double time = (*Id)->getDatum();     // true drift time (ns)

      if (deadTime_) {
        if (wire==prv_wire && fabs(time-prv_time)<deadTime_) continue;
	else {
 	  prv_wire = wire; prv_time = time;
 	}
      }

      //               ***** EVEN in MWPC mode, TIME HISTOGRAMS ARE FILLED *****
      if (hLevel_>=Normal) {
        H_t_->Fill(time); H_ch_->Fill(wire);
      }

      //                                      ***** SET COORDINATES in WRS *****
      // u: perpendicular to the detector wires direction
      // v: parallel to the detector wires direction
      // w: = Z coordinate in MRS
      // origin: the origin of the MRS
      double u = wireDCorr + Wire2Pos(wire); // *ModifComm*: CsDetector's "Wire2Pos" does not return the position corrected for "wirD_", but the distance to 1st wire. Let's do the same here, in order to avoid confusion.
      double v = 0, w = zcm_;

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
  else {                // ********** PROCESS DR's USING DRIFT TIME **********
    list<CsDigit*>::iterator Id;
    list<CsDigit*> digits = getMyDigits();

    if (digits.empty()) return;  // check number of digits

    HepMatrix iRotM(3,3); iRotM = getRotWRSInv();
    double wireDCorr = wirD_ + iRotM(1,1) * _deltax + iRotM(1,2) *_deltay; 
    double res  = vel_ * spSli_;       // Detector resolution (mm)
    //double vel = vel_ / tiSli_;      // Drift velocity (mm/ns)

    bool isMC = CsEvent::Instance()->isAMonteCarloEvent();

//     cout << setprecision(10) << "wirD_ = " <<  wirD_ << " wireDCorr = " <<  wireDCorr << " iRotM(1,1) = " << iRotM(1,1) 
// 	 << " _deltax = " << _deltax << " iRotM(1,2) = " << iRotM(1,2)
// 	 << " _deltay = " << _deltay << endl << flush;

    int prv_wire; static double prv_time; // "static" to avoid warning from compiler
    for (Id = digits.begin(), prv_wire = -1, prv_time = -10000000.; Id!=digits.end(); Id++) {

      int wire = int((*Id)->getAddress()); // Fired wire
      double time = (*Id)->getDatum();     // Drift time (ns)

      if (hLevel_>=Normal) {                               // ***** BASIC HISTOS
        H_t_->Fill(time); H_ch_->Fill(wire);
      }

//      cout << "Time Gate = " << tGate_ << ", time = " << time << endl << flush ;

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

      //      if (deadTime_) {
      //        if (wire==prv_wire && fabs(time-prv_time)<deadTime_) continue;
      //        else {
      //          prv_wire = wire; prv_time = time;
      //        }
      //      }   // Yu.Kh. Will apply that  after discarding negative and extra-large time 

      //    ********** SET THE COORDINATES IN THE WRS **********
      // u: perpendicular to the detector wires direction
      // v: parallel to the detector wires direction
      // w: = Z coordinate in MRS
      // origin: the origin of the MRS

      //                                    ***** DRIFT TIME -> DISTANCE TO WIRE
      bool error; double drift = getDistToWire(time,error);

//       cout << setprecision(10) << "drift = " << drift << endl << flush;

      // Yu. Kh. Uncommented this :
      if( error ) continue;
      // End. Yu.Kh.

      if (deadTime_) {
	//	cout << "Dead Time = " << deadTime_ << endl << flush;
	if (wire==prv_wire && fabs(time-prv_time)<deadTime_){
	//	  cout << "Time difference of double hit(" << prv_wire <<"-" << wire << ") = " <<  fabs(time-prv_time) << endl << flush;
	  continue;
        }else {
	  prv_wire = wire; prv_time = time;
	}
      }

      double u = wireDCorr + Wire2Pos(wire); // *ModifComm*: CsDetector's "Wire2Pos" does not return the position corrected for "wirD_", but the distance to 1st wire. Let's do the same here, in order to avoid confusion.
      double v = 0, w = zcm_;
      //         ********** SET ERRORS ********** 
      HepMatrix cov(3,3,0);  // Zero matrix
      cov(1,1) = pow( res, 2 );    // Drift detector resolution 
      if (fabs(ang_-90)<1) cov(2,2) = pow(xsiz_/2,2);
      else                 cov(2,2) = pow(ysiz_/2/cos(ang_/180.*M_PI),2);
      cov(3,3) = 0.4;              // Assume 0.4 mm resolution in Z
    
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
  } // is_MWPC ?
  
  sortClusters(); setClusteringDone();
}

double CsRichWallDetector::getCorrU(const CsCluster *c,
				    const double x, const double y, // MRS, mm
				    const double tt,  // Offset w.r.t. *event* time
				    bool &error)
{

  //  cout << "nWir_ = " << nWir_ << " ,wirP_ = " << wirP_ << endl << flush;

  CsDigit *dg = c->getDigitsList().front();
  int nch=dg->getAddress();

  // Time propagation is a 2nd order effect: no need to go through a precise and time consuming transform
  //double xx,yy,zz; rotatePointMRS2WRS(x,y,0,xx,yy,zz);

  //  cout << "x= " << x << " ,y= " << y << " ,tt= " << tt << endl << flush; 

  double xx,yy,zz;
  rotatePointMRS2WRS(x,y,0,xx,yy,zz);

  if(getNWir() != 592){
    yy = -yy; 
  }

  //  cout << "xx= " << xx << " ,yy= " << yy << " ,zz= " << zz << endl << flush; 

  double Ui = c->getU();

  //  cout << setprecision(10) << "Ui_def = " <<  Ui << ", " ;

  double t = 0;
  c->getTime(t);

  //  cout << setprecision(10) << "Time_def = " <<  t << ", " ;

  // Signal propagation time velocity.
  //static const double v = 333; // *ModifCom*: I made a member, taht can be set by option. Default inverse velocity is 5 ns/m => significantly slower than light

  double propagation_distance = 0;
  // "propagation_distance"
  // - ==0 corresponds to mid-plane. Because fixed by calibration, which in turn
  //  is determined on a tracks distribution which, a priori, is peaked in
  //  mid-plane
  // - It's measured along perp coordinate axis, which equates one of COMPASS
  //  axes in the DR case, DR being either DRX or DRY.
  // - Therefore it's = perp coordinate...
  // - ...Times +/-, depending upon whether signal propagates according or
  //  contrary to COMPASS coordinate system: 2 cases:
  //    - DetectorSection == 3: signal goes to ??? // *ModifCom: Top/Jura or bottom/Saleve
  //    - DetectorSection != 3: signal goes to ???
  // - ...Corrected by offset: mainly useful for the Y planes, where beam (due
  //  to its bending by SM1), and hence detector's center, departs from COMPASS
  //  axis.

  //  cout << setprecision(10) << "getXsiz()/2 = " <<  getXsiz()/2 << ", " << "getYsiz()/2 = " <<  getYsiz()/2 << endl << flush;

//   if (nWir_==592) {  // x coord
//     if (DetectorSection( nch ) == 3)
//       //propagation_distance = ysiz_/2 - yy - 65;
//       propagation_distance = -y+ycm_;
//     else
//       //propagation_distance = ysiz_/2 + yy - 65;
//       propagation_distance = y-ycm_;
//   }
//   else {             // y coord
//     if (DetectorSection( nch ) == 3)
//       //propagation_distance = xsiz_/2 - yy - 65; // *ModifCom*: I guess you meant "xx" instead of "yy" => I replaced "y" b "x", but am not sure about the sign.
//       propagation_distance = -x+xcm_;
//     else
//       //propagation_distance = xsiz_/2 + yy - 65 ;
//       propagation_distance = x-xcm_;
//   }

  if ( nWir_ == 592 ){ //x coord
   
    if (DetectorSection( nch ) == 3) {

      propagation_distance = getYsiz()/2 - yy;

    }else{

      propagation_distance = getYsiz()/2 + yy;
     
    }

  }else{ //y coord

    if (DetectorSection( nch ) == 3) {

      propagation_distance = getXsiz()/2 - yy;

    }else{
      
      propagation_distance = getXsiz()/2 + yy;

    }

  }


  if (propagation_distance < 0) propagation_distance = 0; // *ModifCom*: propagation distance is <0 half of the time

  //                 ***** SPLIT TUBE? *****
  // *ModifCom*: This what the code could like, had the "wireP" feature of the decoding library been used:
  // Check data[1] => If indeed a split tube, correction is always >0.
  //if (dg->getData()[1]!=0) propagation_distance = fabs(propagation_distance);

  double t_cor = t - propagation_distance*propVel_inv_ - tt;

  //  cout << " t_cor = " << t_cor << endl << flush;

  //  if (propagation_distance < 0) propagation_distance = 0;

  if( t_cor<0 || t_cor>tGate_ )
    {
      // Correction has moved cluster out of gate..
      error = true;
    }
  else
    {
      CsCluster *mirror = c->getMirrorCluster();
      if( mirror!=NULL )
	{
	  double dR = rt_->getRfromT(t_cor,error);   // Corrected R
	  dR -= getDistToWire(t,error);              // Delta R

	  error = false;
	  if (Ui>mirror->getU())
	    Ui += dR;
	  else
	    Ui -= dR;
	}
      else
	{
	  // No mirror => On wire! => Unable to tell which sign to assign
	  error = false;
	  // ... to the correction => Return U unchanged
	  //return Ui;
	}
    }
    
  //  cout << " Ui_cor = " << Ui << endl << flush;
  
  return Ui;

}

//____________________________________________________________________________
void CsRichWallDetector::updateClusters(double time)
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
    double drift =       // ***** DISTANCE TO WIRE *****
      // Error is ignored here, whereas it's not in "clusterize" at variance w/
      // what is done in other drift detector classes. Anyway we assume that all
      // very bad cases have already been filtered out. Marginally bad clusters
      // can show up between of the time offset introduced here: never mind...
      getDistToWire(t,error);
    if (drift==0)// +1 micron to keep the ordering of (cluster,mirror) unchanged
      drift += .001;
    double u = wireDCorr + Wire2Pos(wire);
    if (m->getU()>=c->getU()) {    // This is what's expected...
      // ...given that the clusters may have been ordered by now.
      // (Note that the two cluster/mirror can be strictly equal: cf. the
      // "CsDriftChamberDetector" version of this "updateClusters" method for
      // explanations.)
      c->setU(u-drift); m->setU(u+drift);
    }
    else {                   // ...nevertheless provide also for the alternative
      c->setU(u+drift); m->setU(u-drift);
    }
  }
}

//------------------------------------------------------------------------------
void CsRichWallDetector::readCalibration(time_t timePoint){
  tm *t = localtime(&timePoint);
  CalibrationData* calib = 0;

  // read-in corresponding calibration constants
  CDB::Time tp(timePoint,0);
  string data;
  cdb_->read(GetTBName(),data,tp);
  if ((getenv("COND_DB_DEBUG")!=0) || (getenv("COND_DB_DEBUG_DR")!=0)) {
    cout <<GetTBName()<<endl;
    cout << data << endl;
  }

  if(data.size()!=0){
    istringstream is(data);
    calib = new CalibrationData(*this);
    is>>(*calib);

    // update t0:
    t0_ = calib->t0;
#ifdef CsRW_DEBUG
    std::cout << "CRichWallDetector::readCalibration() ==> "
	      << GetTBName().c_str()
	      <<" t0="
	      << std::setw(10) << std::setprecision(2) << t0_
	      <<" ns ";
    cout << "(" 
      << t <<").\n";
#endif

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
      cout << "CsRichWallDetector::readCalibration() ==> " << GetTBName() << " RT relation from database not valid\n";
    }

  }else{
    hasRT_ = false;
    if( calib ) delete calib;
    cout<<"CsRichWallDetector::readCalibration() ==> "<<GetTBName()<<" calibrations, valid for "
	<< t <<", not found in DB"<<endl;
  }
}

//=======================
//=== PRIVATE METHODS ===
//=======================

bool CsRichWallDetector::_readRTRelation( void )
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
istream& operator>>(istream& in,CsRichWallDetector::CalibrationData &c) {
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

//LS Eff                                                                                            
void CsRichWallDetector::readMCEffMaps(time_t timePoint)                                          
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
        throw std::logic_error("CsRichWallDetector_MCEff::read(): bad line structure.");          
      MCEff_.push_back(mcEf);                                                                       
      MCEff_err_.push_back(mcEf_err);                                                               
                                                                                                    
    }                                                                                               
}       
