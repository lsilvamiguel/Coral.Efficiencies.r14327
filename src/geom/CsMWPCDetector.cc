// $Id: CsMWPCDetector.cc 14094 2015-11-06 15:28:48Z lsilva $

/*!
   \file    CsMWPCDetector.cc
   \brief   Compass MWPC like detector Class.
   \author  Benigno Gobbo
   \version $Revision: 14094 $
   \date    $Date: 2015-11-06 16:28:48 +0100 (Fri, 06 Nov 2015) $
*/

#include "CsMWPCDetector.h"

#include <math.h>
#include <string.h>
#include <strings.h>
#include <functional>
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
#include "CsMCTrkHit.h"
#include "CsRandom.h"
#include "CsGeant3.h"
#include "CsMCTrack.h"
#include "CDB.h"
//LS Eff
#include <stdexcept>
//#define DEBUG_EFF
#include "CsRandom.h"
  
#include "DaqDataDecoding/ChipF1.h"

using namespace std;
using namespace CLHEP;

extern QhitType Qhit;

const float CsMWPCDetector::F1bin_ = 0.12892; // ns

CsMWPCDetector::CsMWPCDetector( const int    row,
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
			const double tGate, const double vel,
                        const double t0,    const double thRes,
                        const double spSli, const double tiSli) :
  CsDetector( row, id, name, TBname, unit, type, rdLen, xsiz, ysiz, zsiz,
			xcm, ycm, zcm, rotDRS, rotWRS, wirD, ang,   
		        nWir, wirP, eff, bkg, tGate ),
  vel_(vel),    
  t0_(t0), 
  thRes_(thRes),
  spSli_(spSli),
  tiSli_(tiSli),
  hist_(NULL),
  histCalib_(NULL),
  histCluster_(NULL)
{
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
	if( *Is == "MW" || *Is == "MWPC" || *Is == "all" ) {
	  decodeCard_ = true;
	}
      }
    }      
    else {
      decodeCard_ = true;
    }
  }

  // ***** CLUSTERISATION OPTION *****

  key = "no associated MWPC digits";
  associate_ = true;
  if( CsOpt::Instance()->getOpt( tag, key ) ) {
    associate_ = false;
  }

  // ***** Drift time simulation in MC decoding OPTION *****

  driftTimeMC_ = false; 
  if (CsOpt::Instance()->getOpt(GetTBName(),"driftTimeMCDecoding") ||
      CsOpt::Instance()->getOpt(GetTBName().substr(0,2),"driftTimeMCDecoding"))  {
    driftTimeMC_ = true;
    CsErrLog::msg(elInfo,__FILE__,__LINE__,
		  "%s: drift time simulation in MC decoding will be used",TBname);
  }
  
  //LS Eff
  //                                          ***** eff. MAP ENABLED? *****
  mcEffMapsEnabled_ = false;
  if (CsOpt::Instance()->getOpt("ALL"      ,            "mcEffMapsEnabled") ||
      CsOpt::Instance()->getOpt(GetTBName(),            "mcEffMapsEnabled") ||
      CsOpt::Instance()->getOpt(GetTBName().substr(0,2),"mcEffMapsEnabled")) {
    mcEffMapsEnabled_ = true;
#ifdef DEBUG_EFF
    cout<<"**** LS *** : CsMWPCDetector:: "<< mcEffMapsEnabled_<<endl;  
#endif
  }
  else CsErrLog::msg(elWarning,__FILE__,__LINE__,
                     "%s: eff. Map will NOT be used",TBname);

  {
    // ***** CUT ON Cluster Time *****

    // This cut can be made dependent upon detector's TB name, e.g:
    // PC       HitTime [0-1] -25 25
    // PA01X1__ HitTime [0-1] -50 50
    string PC = "PC";  // tag is "PC" or TB name.
    cltMin_ = -1.e+10; cltMax_ = 1.e+10;      // Default: no time cut
    key = "ClusterTime";
    vector<double> v;
    if (CsOpt::Instance()->getOpt(TBname,key,v) ||
	CsOpt::Instance()->getOpt(PC,key,v)) {
      if (v.size()!=2) {
	CsErrLog::mes(elFatal,
		      "Error syntax in specifying Hit Time cut");
      }
      else {
	cltMin_ = v[0]; cltMax_ = v[1];
	CsErrLog::msg(elInfo,__FILE__,__LINE__,"\"%s\" %f < ClusterTime < %f",
		      TBname,cltMin_,cltMax_);
      }
    }

    // ***** CUT ON Hit Time *****

    if (CsInit::Instance()->IsAMonteCarloJob()) {
      hittMin_ = -tGate_/2/F1bin_; hittMax_ = tGate_/2/F1bin_;
    }
    else {
      hittMin_ = -75/F1bin_; hittMax_ = 175/F1bin_;
    }

    // This cut can be made dependent upon detector's TB name, cf. supra:
    key = "HitTime";
    if (CsOpt::Instance()->getOpt(TBname,key,v) ||
	CsOpt::Instance()->getOpt(PC,key,v)) {
      if (v.size()!=2)
	CsErrLog::mes(elFatal,
		      "Error syntax in specifying Hit Time cut");
      hittMin_ = v[0]; hittMax_ = v[1];
      CsErrLog::msg(elInfo,__FILE__,__LINE__,"\"%s\" %f < HitTime < %f",
		    TBname,hittMin_,hittMax_);
      if (hittMin_*F1bin_<-tGate_/2 || hittMax_*F1bin_>tGate_/2)
	CsErrLog::msg(elWarning,__FILE__,__LINE__,
		      "\"%s\" HitTime window [%f,%f] > tGate %f",
		      TBname,hittMin_*F1bin_,hittMax_*F1bin_,tGate_);
    }
  }


  // ***** CALIBRATION OPTIONS *****

  if (CsInit::Instance()->IsAMonteCarloJob()) useCalib_ = false;
  else useCalib_ = CsInit::Instance()->useCalibration();

  // ***** MONTE CARLO: Drift time simulation parameters *****
  // These pars can be made dependent upon detector's TB name, cf. supra:
  key = "driftParsMC";
  vector<double> v;
  if (CsOpt::Instance()->getOpt(GetTBName(),key,v) ||
      CsOpt::Instance()->getOpt(GetTBName().substr(0,2),key,v)) {

    if (v.size()!=5)
      CsErrLog::msg(elFatal,__FILE__,__LINE__,
		    "%s: Syntax error in specifying MC drift parameters!",TBname);
    vel_ = v[0]; t0_ = v[1]; thRes_ = v[2]; spSli_ = v[3]; tiSli_ = v[4];
    CsErrLog::msg(elInfo,__FILE__,__LINE__,"%s: drift vel. %f mm/tsl, T0 %f, "
                  "tsl, 2hitres %f tsl, spaceres %f tsl, time slice %f ns", 
                  TBname,vel_,t0_,thRes_,spSli_,tiSli_);
  }

  if(driftTimeMC_) {
    if(tiSli_==0 || vel_==0 ) {
      CsErrLog::msg(elFatal,__FILE__,__LINE__,
		    "%s: Drift time simulation impossible, wrong time slice or drift velocity %f %f",TBname,tiSli_,vel_);
    }
  }

}

//_________________________________________________
void CsMWPCDetector::BookHistograms( void )
{
  CsDetector::ReadHistLevel();
  if (hLevel_==None) return;  //=== check if histograms are to be booked ===

  CsHistograms::SetCurrentPath("/MWPC");
  
  if (hLevel_>=Normal) {   
    // normal level histograms
    histCluster_ = new CsHist1D(GetTBName()+"_ct",GetTBName()+" Cluster Time",
				200,-65,200);
    // high level histograms
    if (hLevel_>=High) {
      hist_      = new CsHist2D(GetTBName()+"_ht",GetTBName()+" Hit Time",
				nWir_,0,nWir_,200,-15000,500);
      histCalib_ = new CsHist2D(GetTBName()+"_htc",
				GetTBName()+" Hit Time (calib.)",
				nWir_,0,nWir_,200,-500,1500);

      if (CsInit::Instance()->IsAMonteCarloJob()) {
        mH1[GetTBName()+"_htMC"] = new CsHist1D(GetTBName()+"_htMC",
			           GetTBName()+" Digit Time",60,-100,200);
      }

    }
  }

  CsHistograms::SetCurrentPath("/");
  
}  

//_________________________________________________
bool CsMWPCDetector::operator==( const CsMWPCDetector& det ) const {
  return( CsDetector::operator==( det ) );
}

bool CsMWPCDetector::operator<( const CsMWPCDetector& det ) const {
  return( CsDetector::operator<( det ) );
}

void CsMWPCDetector::makeMCDecoding() {


  // should I proceed?
  if( !decode_ && !decodeCard_ ) return;

  // Already done?
  if( decodingDone_ ) return;

  // clear
  myDigits_.clear();

  // get a link to relevant pieces...
  CsEvent* event = CsEvent::Instance();

  double vel = vel_ / tiSli_;           // cast velocity in mm/ns

  // this can be done only with zebra binary files... 

  double extra = // Extra time width added to compensate for the trigger jitter
    // of the current event being larger than reference trigger jitter.
    CsEvent::Instance()->getExtraTimeWidth();
  double nsmin =  hittMin_*F1bin_-extra, // Time window (meant to be strict)
    nsmax = hittMax_*F1bin_+extra;
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

	t -= // Trigger time (and hence offset) is subtracted from hit times
	  triggerOffset;
	//         ********** TIME CUT **********
	// - Contrary to what's done for scifis or microMegas, this cut is not
	//  to be supplemented by a cut on the actual hit time (that would be
	//  performed in the "clusterize" method later on). Therefore need to be
	//  more rigorous here: cf. infra the computation of the drift time and
	//  the taking into account of the time resolution.
	if (!hasDriftTimeMC() && (t<nsmin || nsmax<t))
	  // Here we cut directly on the MC true time. Which should not be
	  // allowed. Let it be so for the time being (givent that I haven't
	  // understood yet the meaning of "hasDriftTimeMC")...
	  continue;

	int err;
	HepMatrix iRotM(3,3); iRotM = rotWRS_.inverse( err );
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
	
	int wireF, wireL;  // First and last wire for this hit
	if ((Ui-wirD_)/wirP_<0) wireF = int( (Ui-wirD_)/wirP_-0.5 );
	else                    wireF = int( (Ui-wirD_)/wirP_+0.5 );
	if ((Uo-wirD_)/wirP_<0) wireL = int( (Uo-wirD_)/wirP_-0.5 );
	else                    wireL = int( (Uo-wirD_)/wirP_+0.5 );
	if (wireL<wireF) {
	  int tmp = wireL; wireL = wireF; wireF = tmp;
	}

	for (int i = wireF; i<=wireL; i++) {
	  //LS Eff
	  if( mcEffMapsEnabled_ && !MCEff_.empty() ){
	    float mceff = MCEff_[i];
	    float random = (float) CsRandom::flat();
#ifdef DEBUG_EFF
	    cout<<"**** LS *** : CsMWPCDetector::makeMCDecoding : "<< GetTBName() <<" - "<<  "MCEff_[" <<i<<"]= "<< mceff <<", random "<<random<<endl;
#endif
	    if ( random > mceff ) break;
	  }
	  
          double tdc = t; 	  
          if (hasDriftTimeMC()) {
	    double driftd;          
	    double ltrack = sqrt((Uo-Ui)*(Uo-Ui)+(Wo-Wi)*(Wo-Wi));
	    if (ltrack>0) {
	      driftd = ((i*wirP_+wirD_-Ui)*(Wo-Wi)-(zcm_-Wi)*(Uo-Ui))/ltrack;
	    }
	    else {
	      driftd = sqrt((Ui-(i*wirP_+wirD_))*(Ui-(i*wirP_+wirD_)) + 
			    (zcm_-Wi)*(zcm_ - Wi));
	    }            	    	    
	    tdc = fabs(driftd)/vel + t + spSli_ * tiSli_ * CsRandom::gauss();
            if (tdc>0) tdc = int( (tdc)/tiSli_ + 0.5 ) * tiSli_ - t0_*tiSli_;
            else       tdc = int( (tdc)/tiSli_ - 0.5 ) * tiSli_ - t0_*tiSli_;

            if (mH1[GetTBName()+"_htMC"]!=NULL)
	      mH1[GetTBName()+"_htMC"]->Fill(tdc);

	    if (tdc<nsmin || tdc>nsmax) continue;
	  }

	  // Look if a digit on this detector with this wire already exists. 
	  // Boring, but no idea how to do in other vay... :(
	  list<CsDigit*>::iterator Id; bool found;
	  for (Id=myDigits_.begin(), found = false;
	       Id!=myDigits_.end() && !found; Id++) {
	    if (i==(*Id)->getAddress()) {
	      found = true; 	      // Here it is, add this hit to it...
	      dynamic_cast<CsMCDigit*>(*Id)->addHit(*(*Ih));
	      if (tdc<(*Id)->getDatum()) (*Id)->replaceDatum( tdc );
	    }
	  }
	  if (!found) { // No digits found with these wire and detector
	    if (i<0 || nWir_<=i) {
	      CsErrLog::msg(elAnomaly,__FILE__,__LINE__,
		"%s: Wire# =%d outside range [0,%d]",GetTBName().c_str(),i,nWir_);
	    }
	    else {
	      CsDigit* digit = new CsMCDigit( *this, i, &tdc );
	      dynamic_cast<CsMCDigit*>(digit)->addHit( *(*Ih) );
	      // add this digit to my list
	      myDigits_.push_back( digit );
	    }
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


////////////////////////////////////////////////////////////////////////////////

void CsMWPCDetector::DecodeChipDigit(const CS::Chip::Digit &digit)
{
  const CS::ChipF1::Digit *d = dynamic_cast<const CS::ChipF1::Digit *>(&digit);
  if( d==NULL )
    throw CS::Exception("CsMWPCDetector::DecodeRawData(): Wrong digit type");
  
  int    wire = d->GetChannel();

  // ***** TIME RELATIVE to TRIGGER *****
  // MWPC F1 chip assumed to have normal precision
  double time = d->GetTimeDecoded()/d->GetTimeUnit();
  if(hist_!=NULL) hist_->Fill(double(wire),time);

  // ***** APPLY CALIBRATIONS *****

  if (useCalib_) {    
    if (badChannels_.find(wire)!=badChannels_.end()) // Discard bad channel
      return;
    time -= calib_data[0];                           // Apply calibration
  }

  if(histCalib_ != NULL) histCalib_->Fill(double(wire),double(time));

  if (time<hittMin_ || time>hittMax_) return;  // ***** CUT on Hit Time *****

  // ***** CONVERSION F1 UNIT-> ns *****
  // proper conversion: PS planes have different F1bin_, so always ask
  // for correct time unit
  time *= d->GetTimeUnit();

 // Create CORAL digit
  myDigits_.push_back( new CsDigit(*this, wire, &time, 1) );
}

////////////////////////////////////////////////////////////////////////////////

void CsMWPCDetector::clusterize() {

  clearClusterList();

  list<CsDigit*>::iterator Id, Id2;
  list<CsDigit*> digits = getMyDigits();
  
  vector<list<CsDigit*>::iterator> iterators;
  vector<list<CsDigit*>::iterator>::iterator Idi;
  iterators.clear();

  // protection
  if( digits.empty() ) return;

  int err;
  HepMatrix iRotM(3,3); iRotM = rotWRS_.inverse( err );
  double wireDCorr = wirD_ + iRotM(1,1) * _deltax + iRotM(1,2) *_deltay; 

  if (associate_) {
    int firstwire, lastwire, previouswire;
    do {
    // Outer loop on digits.
      for( Id=digits.begin(); Id!=digits.end(); ++Id) {
      // skip if already clusterized
        if( (*Id)->isClusterized() ) continue;
        (*Id)->setClusterized();
      
        iterators.clear();
        firstwire    = int((*Id)->getAddress());
        lastwire     = firstwire;
        previouswire = firstwire;
        int delta = 0;
        iterators.push_back( Id );
        CsDigit* d = (*Id);
      
#ifdef MWPC_CLUSTER_DEBUG
        cout << "CsMWPCDetector::clusterize(): cluster seed is " << firstwire
             << endl;
#endif      
	// Then we loop until no hits can be assigned to that cluster anymore.
        bool added;
	do {
	  added = false;
	  // Inner loop on digits.
          for( Id2=digits.begin(); Id2!=digits.end(); ++Id2) {
	    // Skip if already clusterized
	    if((*Id2)->isClusterized() ) continue;
	    // Loop on cluster digits.
	    for(Idi = iterators.begin(); Idi != iterators.end(); ++Idi) {
	      // Test if digit belongs to cluster
	      CsDigit* dd = *(*Idi);
	      int w1 = (int)(*Id2)->getAddress();
	      int w2 = (int)(*(*Idi))->getAddress();
	      delta = abs(int((*Id2)->getAddress()) - (*(*Idi))->getAddress());
	      if( delta <= 1 ) {
		(*Id2)->setClusterized();
		iterators.push_back( Id2 );
		added = true;
		if((*Id2)->getAddress() < firstwire) 
		  firstwire = (*Id2)->getAddress();
		if((*Id2)->getAddress() > lastwire) 
		  lastwire = (*Id2)->getAddress();
#ifdef MWPC_CLUSTER_DEBUG
                cout << "CsMWPCDetector::clusterize(): next wire is " 
		     << (*Id2)->getAddress() << endl;
#endif      
                break;
              }// if
	    }// for(Idi = iterators.begin();
	    // If a new hit has been added to the cluster, the inner loop on hits
	    // is restarted.
	    if(added) break;
	  }// for( Id2=digits.begin();
	} while(added);

        double wire = double( lastwire + firstwire ) / 2.;
        int nwires = abs( lastwire - firstwire ) + 1;

#ifdef MWPC_CLUSTER_DEBUG
        cout << "CsMWPCDetector::clusterize(): cluster = (" << wire
             << "," << nwires << ")" << endl;
#endif      
	// Set the values in the WRS 
	// u: perpendicular to the detector wires direction
	// v: parallel to the detector wires direction
	// w: = Z coordinate in MRS
	// origin: the origin of the MRS
    
        double u = wireDCorr + wire * wirP_;
        double v = 0;
        double w = zcm_;
      
	// Set errors:
        HepMatrix cov(3,3,0);  // Zero matrix
        cov(1,1) = pow( nwires * wirP_ / sqrt(12.0), 2 ); // wire pitch / sqrt(12)
        if( fabs(ang_ - 90 )<1) {           // wire length/2
	  cov(2,2) = pow( (getXsiz())/2., 2 );
        }
        else {
	  cov(2,2) = pow( (getYsiz())/2./cos(ang_/180.*(M_PI)), 2 );
        }
        cov(3,3) = 1.;               // assume 1 mm resolution in Z
      
	// Save the cluster:
        CsCluster* cluster = new CsCluster( u, v, w, cov );
	static double ClusterTime; // "static" to avoid warning "uninitialized" 
        for( unsigned int i=0; i<iterators.size(); i++ ) {
	  CsDigit* Hit = *(iterators[i]);
	  cluster->addDigit(*Hit);
	  double HitTime = Hit->getDatum();
	  if (i==0) ClusterTime = HitTime;
	  else if (HitTime<ClusterTime) ClusterTime = HitTime;
        }
	if (histCluster_) histCluster_->Fill(ClusterTime);
	if (ClusterTime<cltMin_ || ClusterTime>cltMax_) {
	  delete cluster; continue;
	}
	cluster->setTime(ClusterTime);
        cluster->addDet( *this );
        addCluster( *cluster );
      }// for( Id=digits.begin();
    } while( 0 );
  }
  else { // !!!!!!!!!!! Temporary

    for( Id=digits.begin(); Id!=digits.end(); Id++ ) {
      
      // skip if already clusterized
      if( (*Id)->isClusterized() ) continue;
      (*Id)->setClusterized();
      
      int wire = (*Id)->getAddress();

      // Set the values in the WRS 
      // u: perpendicular to the detector wires direction
      // v: parallel to the detector wires direction
      // w: = Z coordinate in MRS
      // origin: the origin of the MRS

      double u = wireDCorr + wire * wirP_;
      double v = 0;
      double w = zcm_;

      // Set errors:
      HepMatrix cov(3,3,0);  // Zero matrix
      cov(1,1) = pow( wirP_ / sqrt(12.0), 2 ); // wire pitch / sqrt(12)
      if( ang_ == 90 ) {           // wire length/2
	cov(2,2) = pow( (getXsiz())/2., 2 );
      }
      else {
	cov(2,2) = pow( (getYsiz())/2./cos(ang_/180.*(M_PI)), 2 );
      }
      cov(3,3) = 1.;               // assume 1 mm resolution in Z
      
      // Save the cluster:
      CsCluster* cluster = new CsCluster( u, v, w, cov );
      cluster->addDigit( *(*Id) );
      cluster->addDet( *this );
      cluster->setTime( (*Id)->getDatum() );
      addCluster( *cluster );
    }
  } // !!!!!!!!!!!!!
  sortClusters();
  setClusteringDone();
}


//----------------------------------------------------------------------------

void CsMWPCDetector::readCalibration(time_t timePoint){

  CDB::Time tp(timePoint,0);

  string strdata("");
  cdb_->read(GetTBName(),strdata,tp);
  istringstream istrdata(strdata);
  istrdata >> calib_data;

  if(calib_data.size()!=0){
    CsErrLog::msg(elInfo,__FILE__,__LINE__,"%s: Calibrations are found!",
		  GetTBName().c_str());  
    if (0<calib_data.size() && calib_data.size()<3) {
      CsErrLog::mes(elError,"Calibrations data size < 3 !");
      useCalib_=false;   
    } else if (calib_data.size()>3) {
      for (unsigned int idata = 3; idata<calib_data.size(); idata++) 
	badChannels_.insert((int)calib_data[idata]);
    }
  }else{
    tm *t = localtime(&tp.first);
    cout << GetTBName() << ", no calibration for local time "
	 << t <<" in CDB"<< endl;
    useCalib_ = false;
  }

  if ((getenv("COND_DB_DEBUG")!=0) || (getenv("COND_DB_DEBUG_MWPC")!=0)) {
    cout << GetTBName() << endl;
    for(unsigned i=0;i<calib_data.size();i++){
      cout << i << " " << calib_data[i] << endl;
    }
  }
}

//LS Eff
 void CsMWPCDetector::readMCEffMaps(time_t timePoint) {
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
	 throw std::logic_error("CsMWPCDetector_MCEff::read(): bad line structure.");
       MCEff_.push_back(mcEf);
       MCEff_err_.push_back(mcEf_err);
     }
 }
