// $Id: CsTriggerHodoDetector.cc,v 1.60 2011/03/01 03:31:27 ybedfer Exp $
/*!
   \file    CsTriggerHodoDetector.cc
   \brief   Compass Trigger Hodoscope like detector Class.
   \author  Benigno Gobbo
   \version $Revision: 1.60 $
   \date    $Date: 2011/03/01 03:31:27 $
*/

#include "CsTriggerHodoDetector.h"

#include <math.h>
#include <string>
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

#include "DaqDataDecoding/ChipF1.h"

using namespace std;
using namespace CLHEP;

extern QhitType Qhit;

const float CsTriggerHodoDetector::F1bin_ = 0.129695; // ns

CsTriggerHodoDetector::CsTriggerHodoDetector( const int    row,
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
  tRes_(tRes), initDecodingDone(false), calcMeanTime(false)
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
	if( *Is == "HO" || *Is == "TriggerHodo" || *Is == "all" ) {
	  decodeCard_ = true;
	}
      }
    }
    else {
      decodeCard_ = true;
    }
  }

  string name1(getName());
  if(name1.compare(3,1,"U")==0) isup_=true;
  else isup_=false;

  
  { //              ***** CUT ON HIT TIME *****

    // This cut can be made dependent upon detector's TB name, e.g:
    // HI       HitTime [0-1] -10 10
    // HI05X1__ HitTime [0-1] -25 25
    string H = "H";  // tag is "H" or TB name.
    hitTMin_ = -1.e+10; hitTMax_ = +1e+10;         // Default: no timing cut
    if (CsInit::Instance()->IsAMonteCarloJob()) {
      hitTMin_ = -tGate_/2; hitTMax_ = tGate_/2;
    }
    key = "HitTime";
    vector<double> v;
    if (CsOpt::Instance()->getOpt(TBname,key,v) ||
	CsOpt::Instance()->getOpt(H,key,v)) {
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
  }

  // add detector IDs to set of ids belonging to this detector
  ids_.clear();
  for (int i=0; i<nWir; i++) {
    ids_.insert(id+i+1);       // slabs detector ID is base detector id + copy number (this starts from 1)
  }

  // check that we did not mess up the set
  if((unsigned) nWir_ != ids_.size()) {
    cout<<nWir_<<" "<<ids_.size()<<endl;
    throw CS::Exception("CsTriggerHodoDetector::CsTriggerHodoDetector() : number of ids and number of channels don't match !");
  }
}

void CsTriggerHodoDetector::AddSubDetector( const int    row,
                                            const int    id,    const char* name,    
                                            const char *TBname,
                                            const int    unit,  const int    type,
                                            const double rdLen, const double xsiz,  
                                            const double ysiz,  const double zsiz,
                                            const double xcm,   const double ycm,   
                                            const double zcm,   const HepMatrix rotDRS,
                                            const HepMatrix rotWRS,
                                            const double wirD,  const double ang,   
                                            const int    nWir,  const double wirP, 
                                            const double eff,   const double bkg,
                                            const double tGate ) {
  // save ids_ set before call to CsDetector::AddSubDetector
  // this method modifies the set, in a way which we do not want to have
  std::set<int> backupIds = ids_;

  // call AddSubDetector of CsDetector parent class
  CsDetector::AddSubDetector( row, id, name, TBname, 
                              unit, type, rdLen,
                              xsiz, ysiz, zsiz, xcm, ycm, zcm,
                              rotDRS, rotWRS, wirD, ang, nWir, wirP,
                              eff, bkg, tGate );
  
  // now put back the ids_ set to what it was before
  ids_ = backupIds;

  // and add the ids of new slabs
  for (int i=0; i<nWir; i++) {
    ids_.insert(id+i+1);       // slabs detector ID is base detector id + copy number (this starts from 1)
  }

  // check that we did not mess up the set
  if((unsigned) nWir_ != ids_.size()) {
    cout<<nWir_<<" "<<ids_.size()<<endl;
    throw CS::Exception("CsTriggerHodoDetector::AddSubDetector() : number of ids and number of channels don't match !");
  }
}

//_____________________________________________________________________________
void CsTriggerHodoDetector::BookHistograms() {

  //=== check if histograms are to be booked ===
  CsDetector::ReadHistLevel();

  // print ids which are associated to this detector in MC
//  cout << GetTBName() << " ";
//  for(IS is=ids_.begin(); is!=ids_.end(); is++) {
//    cout<<*is<<" ";
//  }
//  cout<<endl;

  if( hLevel_ >= Normal ) {
    CsHistograms::SetCurrentPath("/HS");
    string tbn = GetTBName();
    cabs_ = new CsHist1D(tbn+"_cabs",tbn+", Cluster position (LWRS)",
			 100,0,Wire2Pos(nWir_));
    hists1D_.push_back(cabs_);
    hittimes_ = new CsHist1D(tbn+"hittimes",tbn+"Hit Times (ns, calibrated)",
                              400,-0.129695*2000,0.129695*2000);
    hists1D_.push_back(hittimes_);
    hittimesz_ = new CsHist1D(tbn+"hittimesz",tbn+"Hit Times Zoom (ns, calibrated)",
                              400,-0.129695*200,0.129695*200);
    hists1D_.push_back(hittimesz_);

    CsHistograms::SetCurrentPath("/");
  }
}

CsTriggerHodoDetector::~CsTriggerHodoDetector() {

  for(unsigned i=0; i<hists1D_.size(); i++)
    delete hists1D_[i];
  hists1D_.clear();
}

bool CsTriggerHodoDetector::operator==( const CsTriggerHodoDetector& det ) const {
  return( CsDetector::operator==(det) );
}

bool CsTriggerHodoDetector::operator<( const CsTriggerHodoDetector& det ) const {
  return( CsDetector::operator<(det) );
}

void CsTriggerHodoDetector::makeMCDecoding() {

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
//  	double ui = thit->getUin();    // Hit in point (DRS)
//  	double vi = thit->getVin();
//  	double wi = thit->getWin();
//  	double uo = thit->getUout();   // Hit out point (DRS)
//  	double vo = thit->getVout();
//  	double wo = thit->getWout();
	int id = thit->GetID();

	t -= // Trigger time (and hence offset) is subtracted from hit times
	  triggerOffset;
	//         ********** TIME -> TDC and TIME CUT **********
	// - This time cut is to be supplemented by a cut performed in the
	//  "clusterize" method. (Which, in the MC case, is done against the
	//  very same time gate as the present one: no equivalent here of the
	//  ambiguities affecting the definition of the master time in RD).
	// - Therefore a loose cut, performed on "t" (as opposed to "tdc"),
	//  could have been done here, saving th CPU it takes to randomize all
	//  hit times. We prefer to do it already at the MC decoding stage, so
	//  that the random generation be not affected by the actual setting of
	//  time gates: the call to "CsRandom" is done in any case.
	// - This implicitly sets = 0 the relative jitter betweeen CsHit's
	//  originating from a same CsMCDigit. Would the randomisation have been
	//  applied at CsHit instantiation time, the other extreme option, i.e.
	//  relative jitter = (in average) absolute jitter, would have been
	//  taken. Truth lies in between...
	double tdc = t + tRes_ * CsRandom::gauss();
	if (tdc<nsmin || nsmax<tdc) continue;

//  	int err;
//  	HepMatrix iRotM(3,3); iRotM = rotWRS_.inverse( err );
//  	double xi = rotDRS_(1,1)*ui+rotDRS_(1,2)*vi+rotDRS_(1,3)*wi+
//  	  xcm_ - _xshift;
//  	double yi = rotDRS_(2,1)*ui+rotDRS_(2,2)*vi+rotDRS_(2,3)*wi+
//  	  ycm_ - _yshift;
//  	double zi = rotDRS_(3,1)*ui+rotDRS_(3,2)*vi+rotDRS_(3,3)*wi+
//  	  zcm_;
//  	double xo = rotDRS_(1,1)*uo+rotDRS_(1,2)*vo+rotDRS_(1,3)*wo+
//  	  xcm_ - _xshift;
//  	double yo = rotDRS_(2,1)*uo+rotDRS_(2,2)*vo+rotDRS_(2,3)*wo+
//  	  ycm_ - _yshift;
//  	double zo = rotDRS_(3,1)*uo+rotDRS_(3,2)*vo+rotDRS_(3,3)*wo+
//  	  zcm_;
//  	double Ui = iRotM(1,1)*xi+iRotM(1,2)*yi+iRotM(1,3)*zi; // WRS
//  	double Vi = iRotM(2,1)*xi+iRotM(2,2)*yi+iRotM(2,3)*zi;
//  	double Wi = iRotM(3,1)*xi+iRotM(3,2)*yi+iRotM(3,3)*zi;
//  	double Uo = iRotM(1,1)*xo+iRotM(1,2)*yo+iRotM(1,3)*zo; // WRS
//  	double Vo = iRotM(2,1)*xo+iRotM(2,2)*yo+iRotM(2,3)*zo;
//  	double Wo = iRotM(3,1)*xo+iRotM(3,2)*yo+iRotM(3,3)*zo;
	
//  	int wireF;  // first wire for this hit
//  	wireF = Pos2Wire(Ui);
//  	if( (Ui-wirD_)/wirP_ < 0 ) {
//  	  wireF = int( (Ui-wirD_)/wirP_-0.5 );
//  	}
//  	else {
//  	  wireF = int( (Ui-wirD_)/wirP_+0.5 );
//  	}
	//	int wireL; // last  wire for this hit
	//wireL = Pos2Wire(Uo);
//  	if( (Uo-wirD_)/wirP_ < 0 ) {
//  	  wireL = int( (Uo-wirD_)/wirP_-0.5 );
//  	}
//  	else {
//  	  wireL = int( (Uo-wirD_)/wirP_+0.5 );
//  	}
//  	if( wireL < wireF ) {
//  	  int tmp = wireL;
//  	  wireL = wireF;
//  	  wireF = tmp;
//  	}

	//       	for( int i=wireF; i<=wireL; i++ ) {
	
	  // look if a digit on this detector with this wire already exists.
	  // Boring, but no idea how to do in other vay... :(
	int i=id - *(ids_.begin());
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
			
      double bignum = 999999.;
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

////////////////////////////////////////////////////////////////////////////////

void CsTriggerHodoDetector::DecodeChipDigit(const CS::Chip::Digit &digit)
{
  const CS::ChipF1::Digit *d = dynamic_cast<const CS::ChipF1::Digit *>(&digit);
  if( d==NULL )
    throw CS::Exception("CsTriggerHodoDetector::DecodeRawData(): Wrong digit type");

  //        ********** TIME RELATIVE to TRIGGER **********
  double time = d->GetTimeDecoded()/d->GetTimeUnit(); // Time in F1 units
  //                                                     ***** APPLY CALIBRATION
  bool use = true; if(calib.size() != 0) {
    if ((int)calib.size()>d->GetChannel()) {
      use = (calib[d->GetChannel()].flag&1)==0;
      time -= calib[d->GetChannel()].data;
    } else use = false;
  }
  time *= d->GetTimeUnit();                           // Convert to "ns"

  double extra =  // ***** EVENT MAY HAVE LARGE TRIGGER JITTER? ENLARGE t GATE
    // (The expectation for the trigger jitter may be refined later on, after
    // all data including the trigger pattern TDC have been decoded: the present
    // "getExtraTimeWidth" is an upper bound and is used to cut away hits that
    // can't in any case be retained, in order to speed up the processing.)
    CsEvent::Instance()->getExtraTimeWidth();
  if (time<hitTMin_-extra || hitTMax_+extra<time)      // CUT on Hit Time
    return;

  if (use && hLevel_>=Normal) {
    hittimes_->Fill(time); hittimesz_->Fill(time);
  }

  if (use) myDigits_.push_back( new CsDigit(*this, d->GetChannel(), &time, 1) );
}

////////////////////////////////////////////////////////////////////////////////

void CsTriggerHodoDetector::DecodeChipDigits( const CS::Chip::Digits &digits )
{
  if (!initDecodingDone) {
    // All decoding maps
    const CS::Chip::Maps &daq_maps = CsInit::Instance()->getDaqMaps(); 

    // check if this detector is known to DDD
    set<uint16> srcIds;
    daq_maps.GetSrcIDs(GetTBName(), srcIds);
    if (srcIds.empty()) {
      // this detector is not known to DDD directly
      // we will now test for the existance of TDC information on Jura and
      // Saleve side (first seven characters of TBName + "j" or "s". If both
      // are found, the mean time will be calculated here
      const char app[] = { 's', 'j' };
      bool found(true);
      for (unsigned int i=0; i<2; i++) {
        std::string name = GetTBName().substr(0, 7) + app[i];
        daq_maps.GetSrcIDs(name, srcIds);
        found &= !srcIds.empty();
      }

      if (found) {
	    CsErrLog::msg(elWarning, __FILE__, __LINE__,
                      "\"%s\" Using software mean timer",GetTBName().c_str());
        calcMeanTime = true;
      }
    }

    initDecodingDone = true;
  }

  // calculate the mean time here in CORAL, or can we use information from
  // raw data
  if (calcMeanTime) {
    // map from channel to sum of all measured times and number of measurements
    std::map<int32, std::pair<double, unsigned int> > meanTime;

    const char app[] = { 's', 'j' };
    for (unsigned int i=0; i<2; i++) {
      // build the special name
      std::string name = GetTBName().substr(0, 7) + app[i];

      typedef std::multimap<CS::DetID,CS::Chip::Digit*>::const_iterator m_it; // iterator type
      // get all digits for the detector
      std::pair<m_it,m_it> m_range = digits.equal_range(name);

      // loop on all found digits
      for( m_it d_it=m_range.first; d_it!=m_range.second; d_it++ ) {
        const CS::ChipF1::Digit *d = dynamic_cast<const CS::ChipF1::Digit *>(d_it->second);
        if( d==NULL )
          throw CS::Exception("CsTriggerHodoDetector::DecodeChipDigits(): Wrong digit type");

        //        ********** TIME RELATIVE to TRIGGER **********
        double time = d->GetTimeDecoded()/d->GetTimeUnit(); // Time in F1 units
        //                                                     ***** APPLY CALIBRATION
        bool use = true;
        if(calib.size() != 0) {
          if (calib.size()>(d->GetChannel()+i*getNWir())) {
            use = (calib[d->GetChannel()+i*getNWir()].flag&1)==0;
            time -= calib[d->GetChannel()+i*getNWir()].data;
          } else {
            use = false;
            CsErrLog::msg(elError, __FILE__, __LINE__,
                          "\"%s\" No calibration available for channel %d, channel will not be used", GetTBName().c_str(), d->GetChannel()+i*getNWir());
          }
        }

        // channel has been flagged or no calibration available for this
        // specific channel
        if (!use)
          continue;

        time *= d->GetTimeUnit();                           // Convert to "ns"

        // according to J.Barth there is at maximum one measurement per counter
        // (due to the deadtime of the descriminator of 160ns)
        // so we just sum up and calculate the mean from everything
        if (meanTime.count(d->GetChannel())==0) {
          meanTime.insert(std::pair<int32, std::pair<double, unsigned int> >
                                   (d->GetChannel(), std::pair<double, unsigned int>(time, 1)));
        } else {
          meanTime[d->GetChannel()].first  += time;
          meanTime[d->GetChannel()].second++;
        }
      }
    }

    typedef std::map<int32, std::pair<double, unsigned int> >::const_iterator mT_it_t;
    for (mT_it_t it=meanTime.begin(); it!=meanTime.end(); it++) {
      double time = it->second.first / it->second.second;

      double extra =  // ***** EVENT MAY HAVE LARGE TRIGGER JITTER? ENLARGE t GATE
        // (The expectation for the trigger jitter may be refined later on, after
        // all data including the trigger pattern TDC have been decoded: the present
        // "getExtraTimeWidth" is an upper bound and is used to cut away hits that
        // can't in any case be retained, in order to speed up the processing.)
        CsEvent::Instance()->getExtraTimeWidth();

      bool use(true);
      if (time<hitTMin_-extra || hitTMax_+extra<time)      // CUT on Hit Time
        use = false;

      if (use) {
        if (hLevel_>=Normal) {
          hittimes_->Fill(time); hittimesz_->Fill(time);
        }
        myDigits_.push_back( new CsDigit(*this, it->first, &time, 1) );
      }
    }
  } else
    // use the standard routine
    CsDet::DecodeChipDigits(digits);
}

////////////////////////////////////////////////////////////////////////////////

void CsTriggerHodoDetector::clusterize() {

  clearClusterList();

  list<CsDigit*>::iterator Id;
  list<CsDigit*> digits = getMyDigits();

  // protection
  if( digits.empty() ) return;

  //    ***** INITIALISE CLUSTERISATION *****
  int err; HepMatrix iRotM(3,3); iRotM = rotWRS_.inverse( err );
  double wireDCorr = wirD_ + iRotM(1,1) * _deltax + iRotM(1,2) *_deltay;

  double extra = // Extra time width added to compensate for the trigger jitter
    // or the current event being larger than reference trigger jitter.
    CsEvent::Instance()->getExtraTimeWidth();
  double t0    =   // Wrong master trigger was selected by Decoding library
    CsEvent::Instance()->getTriggerTimeCorr();
  double tmin = t0+hitTMin_-extra, tmax = t0+hitTMax_+extra;

  for( Id=digits.begin(); Id!=digits.end(); Id++ ) {

    double time = (*Id)->getDatum();
    if (time<tmin || tmax<time) continue;               // ***** CUT ON HIT TIME

    int wire = int((*Id)->getAddress());

    // Set the values in the WRS
    // u: perpendicular to the detector wires direction
    // v: parallel to the detector wires direction
    // w: = Z coordinate in MRS
    // origin: the origin of the MRS

    double d = Wire2Pos(wire);
    double u = wireDCorr + d;
    //double u = wireDCorr + wire * wirP_;
    double v = 0;
    double w = zcm_;
    if( hLevel_ >= Normal ) cabs_ -> Fill(d);                // cluster pos in local WRS

    // Set errors:
    HepMatrix cov(3,3,0);  // Zero matrix
    cov(1,1) = pow(Pitch(wire) / sqrt(12.0), 2 );
    if( fabs(ang_) == 90 ) {           // wire length/2
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
  sortClusters();
  setClusteringDone();
}


//-------------------------------------------------------------------------

void CsTriggerHodoDetector::readCalibration(time_t timePoint){

  CDB::Time tp(timePoint,0);

  calib.clear();
//   for (int a=0; a<getNWir(); a++) {
//     calib.push_back(Calib(a, 0.0, 17));
//   }

  string datastr;
  cdb_->read(GetTBName(),datastr,tp);
  istringstream datais(datastr);
  datais >> calib;

  vector<Calib>::iterator it;
  if(calib.size()!=0){
    CsErrLog::msg(elInfo,__FILE__,__LINE__,
		  "%s calibration OK.",GetTBName().c_str());
    for (it=calib.begin();it!=calib.end();it++) {
      int chn=(*it).chn;
      if ((int)calib.size()==2*getNWir() && chn>=getNWir()) chn-=getNWir(); // might be the special case of meantimer in software
      if (chn<0 || chn>=getNWir()) {
	cout << "CsTriggerHodoDetector::readCalibration() invalid channel number " << chn << " in " << GetTBName() << endl;
      }
    }
    if (getenv("COND_DB_DEBUG")!=0 || (getenv("COND_DB_DEBUG_HODO")!=0)) {
      cout << GetTBName() << endl;
      for (it=calib.begin(); it!=calib.end(); it++) {
        cout << (*it).chn << " ";
        cout << (*it).data << " ";
        cout << (*it).flag <<endl;
      }
    }
  } else {
    calib.clear();
    tm *t = localtime(&tp.first);
    cout << GetTBName() << ", no calibration for local time "
         << t <<" in CDB"<< endl;
  }
  //cout<<"T0-test: "<<GetTBName().c_str()<<endl;
  //for(int m=0;m<getNWir();m++){
  //cout<<"BMS-ch "<<m<<"t0: "<<calib_data[m]<<" flag: "<<calib_flag[m]<<endl;
  //}
}


//----------------------------------------------------------------------------

istream& operator>>(istream& in, CsTriggerHodoDetector::Calib &c) {
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

istream& operator>>(istream& in, vector<CsTriggerHodoDetector::Calib> &vc) {
  CsTriggerHodoDetector::Calib c;

//  vc.clear();
  while (in >> c) {
    if (((int)vc.size()) < (c.chn+1)) vc.resize(c.chn+1, CsTriggerHodoDetector::Calib(0, 0.0, 17));
    vc[c.chn] = c;
  }

  if (in.eof()) in.clear();
  return in;
}


