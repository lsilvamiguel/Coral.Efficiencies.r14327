// $Id: CsJFiberHodoDetector.cc,v 1.32 2010/10/17 20:08:13 ybedfer Exp $

/*!
   \file    CsJFiberHodoDetector.cc
   \brief   Compass Fiber Hodoscope like detector Class.
   \author  Benigno Gobbo
   \version $Revision: 1.32 $
   \date    $Date: 2010/10/17 20:08:13 $
*/

#include "CsJFiberHodoDetector.h"

#include <math.h>
#include <string.h>
#include <strings.h>
#include <set>
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

#include "DaqDataDecoding/ChipF1.h"

using namespace std;
using namespace CLHEP;

extern QhitType Qhit;

CsJFiberHodoDetector::CsJFiberHodoDetector( const int    row,
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
			const double tGate, const double tRes ) :
  CsDetector( row, id, name, TBname, unit, type, rdLen, xsiz, ysiz, zsiz,
			xcm, ycm, zcm, rotDRS, rotWRS, wirD, ang,   
	      nWir, wirP, eff, bkg, tGate ),
  tRes_(tRes)
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
	if( *Is == "FI" || *Is == "JFiberHodo" || *Is == "all" ) {
	  decodeCard_ = true;
	}
      }
    }      
    else {
      decodeCard_ = true;
    }
  }
  
  { //              ***** CUT OPTIONS *****

    // These cuts can be made dependent upon detector's TB name, e.g:
    // FI       HitTime [0-1] -10 10
    // FI01X1__ HitTime [0-1] -25 25

    // ***** CUT ON HIT TIME
    string FI = "FI";  // tag is "FI" or TB name.
    hitTMin_ = -1.e+10; hitTMax_ = +1e+10;         // Default: no timing cut
    if (CsInit::Instance()->IsAMonteCarloJob()) {
      hitTMin_ = -tGate_/2; hitTMax_ = tGate_/2;
    }
    key = "HitTime";
    vector<double> v;
    if (CsOpt::Instance()->getOpt(TBname,key,v) ||
	CsOpt::Instance()->getOpt(FI,key,v)) {
      if (v.size()!=2)
	CsErrLog::mes(elFatal,"Error syntax in specifying Hit Time cut");
      else {
	hitTMin_ = v[0]; hitTMax_ = v[1];
	if (hitTMin_<-tGate_/2 || hitTMax_>tGate_/2)
	  CsErrLog::msg(elWarning,__FILE__,__LINE__,
  "\"%s\" HitTime window [%f,%f] > tGate %f",TBname,hitTMin_,hitTMax_,tGate_);
	else
	  CsErrLog::msg(elInfo,__FILE__,__LINE__,"\"%s\" %f < HitTime < %f ns",
			TBname,hitTMin_,hitTMax_);
      }
    }

    // ***** CUT ON HIT TIME DIFF. FOR CLUSTER
    cluDtMax_ = 5;   // Default = 5ns (as had long been the built-in value)
    // But if the local data member ("cluDtmax_") is in unit of ns, let's adopt
    // for the entries in options file a unit of resolution (sigma). Resolution
    // being itself, as of 2015/09, built-in = 0.6ns)
    float timeResolution = 0.6; 
    key = "ClusterDTime";
    if (CsOpt::Instance()->getOpt(TBname,key,cluDtMax_) ||
	CsOpt::Instance()->getOpt(FI,key,cluDtMax_)) {
      cluDtMax_ *= timeResolution;
      CsErrLog::msg(elInfo,__FILE__,__LINE__,"\"%s\" ClusterDTime < %f ns",
			TBname,cluDtMax_);
    }
  }

  // calibrations
  use_calib = CsInit::Instance()->useCalibration();
  if (CsInit::Instance()->IsAMonteCarloJob()) use_calib = false;

  this->do_clustering=true;
}

//________________________________________________
void CsJFiberHodoDetector::BookHistograms( void )
{
  for (unsigned int i = 0; i<sizeof(mH1)/sizeof(CsHist1D*); i++) mH1[i] = 0;
  for (unsigned int i = 0; i<sizeof(mH2)/sizeof(CsHist2D*); i++) mH2[i] = 0;
  //=== check if histograms are to be booked ===
  CsDetector::ReadHistLevel();
  if( hLevel_ == None ) return;

  CsHistograms::SetCurrentPath("/SciFi");

  string tbn  = GetTBName(); 

  int range_low = -13000;
  int range_high = -9500; // zoom & calibration histogram range

  if (hLevel_>=Normal) {   // ********** NORMAL LEVEL HISTOGRAMS **********

    mH1[1] = new CsHist1D(tbn+"1",tbn+"Time - trig. time (zoom)",(range_high-range_low), range_low, range_high);
    mH1[2] = new CsHist1D(tbn+"2",tbn+"Hit pattern", nWir_, 0, nWir_);

    //                                           ***** CLUSTERIZATION HISTOGRAMS
    mH1[5] = new CsHist1D(tbn+"15",tbn+"Cluster size", 50, 0, 50);

    if (this->use_calib) {     //                   ***** CALIBRATION HISTOGRAMS
      mH1[4] = new CsHist1D(tbn+"9",tbn+"Time - trig. time (calibr,closer zoom),ns",800, -80, 80);
    }
  }

  if (hLevel_>=High) {   //    ********** HIGH LEVEL HISTOGRAMS **********
    mH1[0] = new CsHist1D(tbn+"0",tbn+"Time - trig. time",       128,  -65536, 65536);
    mH1[6] = new CsHist1D(tbn+"16",tbn+"Time diff. of neighbour fibres  [ns]", 120, -30, 30);
    mH1[7] = new CsHist1D(tbn+"17",tbn+"Time diff. of same-fibre digits [ns]", 120,   0, 60);
    mH2[0] = new CsHist2D(tbn+"3",tbn+"Fibre VS time", (range_high-range_low), range_low, range_high, nWir_, 0, nWir_);
    if(this->use_calib){
      mH1[3] = new CsHist1D(tbn+"7",tbn+"Time - trig. time (calibr,closer zoom)",1200, -600, +600);
      mH2[1] = new CsHist2D(tbn+"4",tbn+"Fibre VS time (calibrated)",200, -200, 200, nWir_, 0, nWir_);
    }
  }

  CsHistograms::SetCurrentPath("/");

}


//________________________________________________

void CsJFiberHodoDetector::DecodeChipDigit(const CS::Chip::Digit &digit)
{
  const CS::ChipF1::Digit *d = dynamic_cast<const CS::ChipF1::Digit *>(&digit);
  if( d==NULL )
    throw CS::Exception("CsJFiberHodoDetector::DecodeRawData(): Wrong digit type");

  if(d->GetChannelPos() != 0) {
    // tmp
    //cout<<this->GetTBName()<<"  ChanPos = "<<d->GetChannelPos()<<endl;
    if(d->GetChannelPos() > 0) return; // skip high threshold 
  }
  
  if( d->GetChannel() < 0 || d->GetChannel() >= int(getNWir()))
    throw CS::Exception("CsJFiberHodoDetector::DecodeRawData() ==> Unexpected wire number %d for SiFi %d %s",
                         d->GetChannel(),GetID().GetNumber(), GetTBName().c_str());

  // ***** TIME RELATIVE TO TRIGGER TIME *****
  double time = d->GetTimeDecoded()/d->GetTimeUnit(); 

  if (mH1[0]) mH1[0]->Fill(time);
  if (mH1[1]) mH1[1]->Fill(time);
  if (mH1[2]) mH1[2]->Fill(d->GetChannel());
  if (mH2[0]) mH2[0]->Fill(time,d->GetChannel());

  // ***** APPLY CALIBRATIONS *****
  bool useIt=true;
  if(d->GetChannel()<0 || d->GetChannel() >= int(calib.size()) ) {
    if(calib.size() != 0) {
      CsErrLog::msg(elError,__FILE__,__LINE__,
		    "SciFi %s wire #%d but calib. container size = %d",
		    GetTBName().c_str(),d->GetChannel(),calib.size());
      useIt=false;
    }
  }
  else {
    if(calib[d->GetChannel()].flag!=0 && calib[d->GetChannel()].flag!=2)
      useIt=false;
    time -= calib[d->GetChannel()].data; // apply calibration

    if (hLevel_>=Normal) {
      if (mH2[1]) mH2[1]->Fill(time,d->GetChannel());
      if (mH1[3]) mH1[3]->Fill(time);
      if (mH1[4]) mH1[4]->Fill(time*d->GetTimeUnit());
    }
  }
  
  time *= d->GetTimeUnit(); // convert "bins" to "ns"

  double extra =  // ***** EVENT MAY HAVE LARGE TRIGGER JITTER? ENLARGE t GATE
    // (The expectation for the trigger jitter may be refined later on, after
    // all data including the trigger pattern TDC have been decoded: the present
    // "getExtraTimeWidth" is an upper bound and is used to cut away hits that
    // can't in any case be retained, in order to speed up the processing.)
    CsEvent::Instance()->getExtraTimeWidth();
  if (time<hitTMin_-extra || hitTMax_+extra<time)      // CUT on Hit Time
    return;

   double data[2] = {time,(double)d->GetChannelPos()};

   // Create CORAL digit
   if(useIt) myDigits_.push_back( new CsDigit(*this, d->GetChannel(), data, 2) );
}

//------------------------------------------------------------------------------
//------------------------------------------------------------------------------

//
// Local classes just for interfacing
// with FI clusterization procedure
// (Used in CsJFiberHodoDetector::clusterize())
//

// input digits
class FIdig 
{
public:
  int    wire;  // fibre #
  double time;  // time
  CsDigit* ref; // reference to original CsDigit 
  bool operator < (const FIdig& gd) const 
  { 
    if(wire == gd.wire) return (time < gd.time);
    else                return (wire < gd.wire);
  };
};

// output clusters
class FIcl
{
public:
  double dwire; // cluster center position in the units of "wires"
  double sigma; // error of cluster center position (also in "wires")
  double time;  // cluster time 
  double sigt;  // error of cluster time
  list<CsDigit*> lref; // list of referncies to original CsDigits
  bool operator < (const FIcl& gc) const 
  { 
    if(dwire == gc.dwire) return (time  < gc.time);
    else                  return (dwire < gc.dwire);
  };

};

//
// FI clusterization. 
// (Used in CsJFiberHodoDetector::clusterize())
// 
// input :  list of digits    (FIdig)
// output:  list of clusters  (FIcl)
//
// do_clustering - flag to switch ON/OF clusterization
// mh - array of histograms

void JFIclusterize(bool do_clustering, list<FIdig>& ld, list<FIcl>& lc, CsHist1D **mH, double cluDtMax_)
{
  assert(lc.empty()); // ouput list must be empty

  list<FIdig>::iterator id, idnext;

  if(do_clustering) { // clustering is ON
    // Variables to be reset when opening a new cluster
    // (static to prevent warning: might be used uninitialized)
    static double dwire, time, time2; static int ndig;
    FIcl gc;          // create "working" cluster object
    bool opened=false;
    for(id = ld.begin(); id != ld.end(); id++) { // loop over digits
      if( !opened) { // open new cluster (clear all, reset FIcl obj.)
	dwire = 0; time = 0; time2=0;
	ndig = 0;
	gc.lref.clear();
	opened=true;
      }
      dwire += (*id).wire;
      time  += (*id).time;
      time2 += (*id).time * (*id).time;
      gc.lref.push_back((*id).ref); // add digit ref. to cluster
      ndig++;

      idnext=id; idnext++;          // look on the next digit and close cluster if ...
      if(idnext == ld.end())                     opened=false; // just an end
      if(    ((*idnext).wire - (*id).wire) > 1)  opened=false; // cut on distance (in "wires")
      if(fabs((*idnext).time - (*id).time) > cluDtMax_)  opened=false; // cut on time difference (ns).

      if(!opened) { // cluster had been closed. Save it.
	gc.dwire = dwire/ndig; // simple mean
	//	gc.sigma = -1;         // means "use default"
	gc.sigma = ndig;         // means "use cluster_size"
	gc.time  =  time/ndig;
	/*
	if(ndig > 1) { // time dispersion estimation
	  gc.sigt  = sqrt((time2 - time*time/ndig)/(ndig-1));
	} else {
	  gc.sigt  = -1; // means "use default"
	} 
	*/
	gc.sigt  = -1; // "use default" i.e. detectpr's time resolution
	lc.push_back(gc); // save (copy) found cluster

	if (mH[5]) mH[5]->Fill(gc.lref.size()+0.5); // cluster size histogram
      }
      // just histogramming
      if(idnext != ld.end()) {
	float dt = (*idnext).time -  (*id).time; // neighbor digits' time difference
	if (mH[6]) mH[6]->Fill(dt);
	if((*idnext).wire == (*id).wire) { if (mH[7]) mH[7]->Fill(dt); }
      }
    }// end of loop over digits

  } else { // clustering is OFF. Primitive (1 digit --> 1 cluster) procedure instead.

    for(id = ld.begin(); id != ld.end(); id++) { // loop over digits
      FIcl gc;                       // create new cluster candidate
      gc.dwire = double((*id).wire); // store cluster position in terms of "wires"
      gc.sigma = -1;                 // means: "cluster position error is unknown"
      gc.time  = (*id).time;
      gc.sigt  = -1;
      gc.lref.push_back((*id).ref);  // store into reference to original CsDigit
      
      lc.push_back(gc); // save found cluster
    }// end of loop over digits

  }
}

//------------------------------------------------------------------------------
//------------------------------------------------------------------------------


void CsJFiberHodoDetector::clusterize() {

  clearClusterList();

  list<CsDigit*>::iterator Id;
  list<CsDigit*> digits = getMyDigits();
  
  vector<list<CsDigit*>::iterator> iterators;
  iterators.clear();

  // protection
  if( digits.empty() ) return;
  int err;
  HepMatrix iRotM(3,3); iRotM = rotWRS_.inverse( err );
  double wireDCorr = wirD_ + iRotM(1,1) * _deltax + iRotM(1,2) *_deltay; 

  FIdig gd; 
  list<FIdig> ld; // list of digits (input for clusterization)
  FIcl gc; 
  list<FIcl>  lc; // list of clusters (output of clusterization) 

  //    ********** PREPARE INPUT FOR "FIclusterize" **********
  double extra = // Extra time width added to compensate for the trigger jitter
    // of the current event being larger than reference trigger jitter.
    CsEvent::Instance()->getExtraTimeWidth();
  double t0    =   // Wrong master trigger was selected by Decoding library
    CsEvent::Instance()->getTriggerTimeCorr();
  double tmin = t0+hitTMin_-extra, tmax = t0+hitTMax_+extra;

  for (Id = digits.begin(); Id!=digits.end(); Id++) {  // ***** LOOP OVER DIGITS
    assert( (*Id)->getDet() == this );
    if ((*Id)->isClusterized()) continue;   // ***** SKIP IF ALREADY CLUSTERIZED
    double time = (*Id)->getDatum();
    if (time<tmin || tmax<time) continue;               // ***** CUT ON HIT TIME
    (*Id)->setClusterized();
    gd.wire =  (*Id)->getAddress();
    gd.time =  time;
    gd.ref  =  *Id;   // Pointer to CsDigit
    ld.push_back(gd); // save FI digit
  } // end of loop over digits


  ld.sort();            // sort by wire# and time
  JFIclusterize(this->do_clustering, ld, lc, mH1, cluDtMax_); // do FI clusterizetion <==========

  // store found clusters to CsCluster objects
  list<FIcl>::iterator ic;
  for(ic = lc.begin(); ic != lc.end(); ic++){ // loop over found clusters

    double u;
    if ( IsVarPitch() )
      u = wireDCorr + Wire2Pos( (*ic).dwire );
    else
      u = wireDCorr + (*ic).dwire * wirP_;
    double v = 0;
    double w = zcm_;

    // Set errors
    HepMatrix cov(3,3,0);  // Zero matrix
    if((*ic).sigma < 0 ) {
      if ( IsVarPitch() )
        cov(1,1) = pow( Pitch((*ic).dwire) / sqrt(12.0), 2 ); // if cluster error is unknown
      else
        cov(1,1) = pow( wirP_ / sqrt(12.0), 2 ); // if cluster error is unknown
    } else {
      if ( IsVarPitch() )
        cov(1,1) = (*ic).sigma/sqrt(12.0)*Pitch((*ic).dwire) * (*ic).sigma/sqrt(12.0)*Pitch((*ic).dwire);
      else
        cov(1,1) = (*ic).sigma/sqrt(12.0)*wirP_ * (*ic).sigma/sqrt(12.0)*wirP_;
    }
    cov(2,2) = getXsiz()*getXsiz() + getYsiz()*getYsiz(); // just very big error (diagonale)
    cov(3,3) = 10.;

    // create CsCluster
    CsCluster* cluster = new CsCluster( u, v, w, cov );
    // save CsDigits the CsCluster is made from
    list<CsDigit*>::iterator idig, idig_next;
    for(idig = (*ic).lref.begin(); idig != (*ic).lref.end(); idig++){ 
      cluster->addDigit(*(*idig) ); // add digits to the cluster
    }

    cluster->setTime((*ic).time, (*ic).sigt);
    cluster->addDet( *this );
    addCluster( *cluster );
  } // end of loop over found clusters

  sortClusters();
  setClusteringDone();
}

//------------------------------------------------------------------------------
//------------------------------------------------------------------------------




bool CsJFiberHodoDetector::operator==( const CsJFiberHodoDetector& det ) const {
  return( CsDetector::operator==(det) );
}

bool CsJFiberHodoDetector::operator<( const CsJFiberHodoDetector& det ) const {
  return( CsDetector::operator<(det) );
}

void CsJFiberHodoDetector::makeMCDecoding() {

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

	t -= // Trigger time (and hence offset) is subtracted from hit times
	  triggerOffset;
	//         ********** TIME -> TDC and TIME CUT **********
	// - This time cut is to be supplemented by a cut performed in the
	//  "clusterize" method. (Which, in the MC case, is done against the
	//  very same time gate as the present one: no equivalent here of the
	//  ambiguities affecting the definition of the master time in RD).
	// - Therefore a loose cut, performed on "t" (as opposed to "tdc"),
	//  could have been done here, saving the CPU it takes to randomize all
	//  hit times, including those far outside time gate. We prefer to cut
	//  on tdc already at the present MC decoding stage, so that the random
	//  generation be not affected by the actual setting of time gates: 
	//  "CsRandom" is called in any case, and therefore the sequence of 
	//  random numbers is synchronized on sequence of hits in a unique way.
	// - This implicitly sets = 0 the relative jitter betweeen CsHit's
	//  originating from a same CsMCDigit. Would the randomisation have been
	//  applied at CsHit instantiation time, the other extreme option, i.e.
	//  relative jitter = (in average) absolute jitter, would have been
	//  taken. Truth lies in between...
	double tdc = t + tRes_ * CsRandom::gauss();
	if (tdc<nsmin || nsmax<tdc) continue;

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
	  // Look if a digit on this detector with this wire already exists. 
	  // Boring, but no idea how to do in other vay... :(
	  list<CsDigit*>::iterator Id; bool found;
	  for (Id=myDigits_.begin(), found = false;
	       Id!=myDigits_.end() && !found; Id++) {
	    if (i==(*Id)->getAddress()) {
	      found = true; 	      // Here it is, add this hit to it...
	      dynamic_cast<CsMCDigit*>(*Id)->addHit(*(*Ih));
	      if (tdc<(*Id)->getDatum()) (*Id)->replaceDatum(tdc);
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
	      myDigits_.push_back( digit );   // Add this digit to my list
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
			   
      double bignum = 9999999.;
      CsDigit* digit = new CsMCDigit( *this, dig1-1, &bignum );
      list<CsMCHit*>::iterator Ih;
      for( int k=ip1-1; k<ip2-1; k ++ ) {
	int hitn =  Qhit.jpdig[k];
	if( hitn != 0 ) {
	  for( Ih=myMCHits_.begin(); Ih!=myMCHits_.end(); Ih++ ) {
	    if( ( (*Ih)->getX() == (Qhit.hit[hitn-1][0]*10.) ) &&
		( (*Ih)->getY() == (Qhit.hit[hitn-1][1]*10.) ) ) {
	      dynamic_cast<CsMCDigit*>(digit)->addHit( *(*Ih) );
	      double t = (*Ih)->getDTime();
	      double tdc = t + tRes_ * CsRandom::gauss();
	      if( tdc < digit->getDatum() ) {
		digit->replaceDatum( tdc );
	      }
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


//-------------------------------------------------------------------------

void CsJFiberHodoDetector::readCalibration(time_t timePoint) {

  CDB::Time tp(timePoint,0);

  calib.clear();
  string datastr;
  cdb_->read(GetTBName(),datastr,tp);
  istringstream datais(datastr);
  datais >> calib;
  
  vector<Calib>::iterator it;
  if(calib.size()!=0){
    CsErrLog::msg(elInfo,__FILE__,__LINE__,
		  "%s calibration OK.",GetTBName().c_str());
    for (it=calib.begin(); it!=calib.end(); it++) {
      int chn=(*it).chn;
      if (chn<0 || chn>=getNWir()) {
	cout << "CsJFiberDetector::readCalibration() invalid channel number " << chn << "in " << GetTBName() << endl;
      }
    }
    string debug_level;
    if ((getenv("COND_DB_DEBUG")!=0) || (getenv("COND_DB_DEBUG_JFIBER")!=0)) {
      cout << GetTBName() << endl;
      for (it=calib.begin(); it!=calib.end(); it++) {
        cout << (*it).chn << " ";
        cout << (*it).data << " ";
        cout << (*it).flag <<endl;
      }
    }
  } else {
    tm *t = localtime(&tp.first);
    cout << GetTBName() << ", no calibration for local time "
         << t <<" in CDB"<< endl;
  }
}


//------------------------------------------------------------------------------

istream& operator>>(istream& in, CsJFiberHodoDetector::Calib &c) {
  string line;
  std::getline(in, line);
  istringstream is(line);
  is >> c.chn;
  is >> c.data;
  c.flag = 0;
  is >> c.flag;
  return in;
}


//----------------------------------------------------------------------------

istream& operator>>(istream& in, vector<CsJFiberHodoDetector::Calib> &vc) {
  CsJFiberHodoDetector::Calib c;

//  vc.clear();
  while (in >> c) {
    if (((int)vc.size()) < (c.chn+1)) vc.resize(c.chn+1, CsJFiberHodoDetector::Calib(0, 0.0, 17));
    vc[c.chn] = c;
  }

  if (in.eof()) in.clear();
  return in;
}

