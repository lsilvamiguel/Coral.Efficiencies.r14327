// $Id: CsMW1Detector.cc,v 1.20 2010/02/11 15:31:19 suhl Exp $

/*!
   \file    CsMW1Detector.cc
   \brief   Compass MW1 detector Class.
   \author  Benigno Gobbo
   \version $Revision: 1.20 $
   \date    $Date: 2010/02/11 15:31:19 $
*/

#include "CsMW1Detector.h"

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

const float CsMW1Detector::F1bin_ = 0.12892; // ns

CsMW1Detector::CsMW1Detector( const int    row,
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
  CsDetector( row, id, name, TBname, unit, type, rdLen, xsiz, ysiz, zsiz,
			xcm, ycm, zcm, rotDRS, rotWRS, wirD, ang,   
		        nWir, wirP, eff, bkg, tGate ),
  hist_(NULL),
  histCalib_(NULL),
  histCluster_(NULL)
{

  hasVarP_ = true; // MW1 is variable pitch det. (see virtual Wire2Pos() below)
  
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

  key = "yes associated MW1 digits";
  associate_ = false;   // No reason for ``association'' on these MW1s
  if( CsOpt::Instance()->getOpt( tag, key ) ) {
    associate_ = false;
  }
  
  //LS Eff                                                                                          
  //                                          ***** eff. MAP ENABLED ? *****
  mcEffMapsEnabled_ = false;
  if (CsOpt::Instance()->getOpt("ALL"      ,            "mcEffMapsEnabled") ||
      CsOpt::Instance()->getOpt(GetTBName(),            "mcEffMapsEnabled") ||
      CsOpt::Instance()->getOpt(GetTBName().substr(0,2),"mcEffMapsEnabled")) {
    mcEffMapsEnabled_ = true;
#ifdef DEBUG_EFF
    cout<<"**** LS *** : CsMW1Detector:: "<< mcEffMapsEnabled_<<endl;  
#endif
  }
  else CsErrLog::msg(elWarning,__FILE__,__LINE__,
                     "%s: eff. Map will NOT be used",TBname);
  {
    // ***** CUT ON Cluster Time *****

    // This cut can be made dependent upon detector's TB name, e.g:
    // MA       HitTime [0-1] -25 25
    // MA01X1__ HitTime [0-1] -50 50
    string MA = "MA";  // tag is "MA" or TB name.
    cltMin_ = -1.e+10; cltMax_ = +1.e+10;      // Default: no cluster time cut
    key = "ClusterTime";
    vector<double> v;
    if (CsOpt::Instance()->getOpt(TBname,key,v) ||
	CsOpt::Instance()->getOpt(MA,key,v)) {
      if (v.size()!=2) {
	CsErrLog::msg(elFatal,__FILE__,__LINE__,
		      "%s: Error syntax in specifying Hit Time cut!",TBname);
      }
      else {
	cltMin_ = v[0]; cltMax_ = v[1];
	CsErrLog::msg(elInfo,__FILE__,__LINE__,"%s: %f < ClusterTime < %f",
		      TBname,cltMin_,cltMax_);
      }
    }

    // ***** CUT ON Hit Time *****

    if (CsInit::Instance()->IsAMonteCarloJob()) {
      hittMin_ = -tGate_/2/F1bin_; hittMax_ = tGate_/2/F1bin_;
    }
    else {
      hittMin_ = -1.e+10; hittMax_ = +1.e+10;      // Default: no hit time cut
    }

    // This cut can be made dependent upon detector's TB name, cf. supra:
    key = "HitTime";
    if (CsOpt::Instance()->getOpt(TBname,key,v) ||
	CsOpt::Instance()->getOpt(MA,key,v)) {
      if (v.size()!=2)
	CsErrLog::msg(elFatal,__FILE__,__LINE__,
		      "%s: Error syntax in specifying Hit Time cut!",TBname);
      hittMin_ = v[0]; hittMax_ = v[1];
      CsErrLog::msg(elInfo,__FILE__,__LINE__,"\"%s\" %f < HitTime < %f",
		    TBname,hittMin_,hittMax_);
      if (hittMin_*F1bin_<-tGate_/2 || hittMax_*F1bin_>tGate_/2)
	CsErrLog::msg(elWarning,__FILE__,__LINE__,
		      "%s: HitTime window [%f,%f] > tGate %f",
		      TBname,hittMin_*F1bin_,hittMax_*F1bin_,tGate_);
    }
  }


  // ***** CALIBRATION OPTIONS *****

  if (CsInit::Instance()->IsAMonteCarloJob()) useCalib_ = false;
  else useCalib_ = CsInit::Instance()->useCalibration();
}

//_________________________________________________
void CsMW1Detector::BookHistograms( void )
{

  //=== check if histograms are to be booked ===
  CsDetector::ReadHistLevel();
  if( hLevel_ == None ) return;

  CsHistograms::SetCurrentPath("/MWPC");
  
  //book histograms
  unsigned int id = (unsigned int) GetID();
  
  // high level histograms
  if( hLevel_ >= High ) {
    char hname[20];
    char hcname[20];
    sprintf(hname,"%s%d","MWPC_",id);
    sprintf(hcname,"%s%s",hname,"_calib");
    hist_        = new CsHist2D(hname, this->GetTBName(), getNWir(), 0, getNWir(), 200, -15000, 500);
    histCalib_   = new CsHist2D(hcname, this->GetTBName(), getNWir(),0,getNWir(),200,-500,1500);
  }  
  
  // normal level histograms
  if( hLevel_ >= Normal ) {   
    histCluster_ = new CsHist1D(GetTBName()+"_ct",GetTBName()+" Cluster Time",
      200,-65,200);
  }
  
  CsHistograms::SetCurrentPath("/");
  
}  

//_________________________________________________
bool CsMW1Detector::operator==( const CsMW1Detector& det ) const {
  return( CsDetector::operator==( det ) );
}

bool CsMW1Detector::operator<( const CsMW1Detector& det ) const {
  return( CsDetector::operator<( det ) );
}

void CsMW1Detector::makeMCDecoding() {


  // should I proceed?
  if( !decode_ && !decodeCard_ ) return;

  // Already done?
  if( decodingDone_ ) return;

  // clear
  myDigits_.clear();

  // get a link to relevant pieces...
  CsEvent* event = CsEvent::Instance();

  // this can be done only with zebra binary files... 

  double nsmin = hittMin_*F1bin_, nsmax = hittMax_*F1bin_;
  if( !CsGeant3::Instance()->isAnNtFile() ) {
    double triggerOffset = // Time offset of the trigger w.r.t. the event
      CsEvent::Instance()->getTriggerMCOffset();
    list<CsMCHit*>::iterator Ih;
    for( Ih=myMCHits_.begin(); Ih!=myMCHits_.end(); Ih++ ) { // loop on hits

      CsMCTrkHit* thit = dynamic_cast<CsMCTrkHit*>(*Ih);
      if( thit == 0 ) return;

      //only charged particles or charged products
      if(  ((thit->getMCTrack())->getParticle())->getCharge() ||
           thit->getOrigin()  ) {


	double t  = thit->getDTime();  // Delay time (ns)
	double ui = thit->getUin();    // Hit in point (DRS)
	double vi = thit->getVin();
	double wi = thit->getWin();
	double uo = thit->getUout();   // Hit out point (DRS)
	double vo = thit->getVout();
	double wo = thit->getWout();

	// Check if hit inside time window (==time gate by default)
	if (t<nsmin || t>nsmax) continue;

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
	
	int wireF;  // first wire for this hit

	if( (Ui-wirD_)/wirP_ < 0 ) {
	  wireF = int( (Ui-wirD_)/wirP_-0.5 );
	}
	else {
	  wireF = int( (Ui-wirD_)/wirP_+0.5 );
	}
	int wireL; // last  wire for this hit
	if( (Uo-wirD_)/wirP_ < 0 ) {
	  wireL = int( (Uo-wirD_)/wirP_-0.5 );
	}
	else {
	  wireL = int( (Uo-wirD_)/wirP_+0.5 );
	} 
	if( wireL < wireF ) {
	  int tmp = wireL;
	  wireL = wireF;
	  wireF = tmp;
	}

	for( int i=wireF; i<=wireL; i++ ) { 	  
	  //LS Eff
          if( mcEffMapsEnabled_ && !MCEff_.empty() ){
            float mceff = MCEff_[i];
	    float random = (float) CsRandom::flat();
#ifdef DEBUG_EFF
            cout<<"**** LS *** : CsMW1Detector::makeMCDecoding : "<< GetTBName() <<" - "<<      
              "MCEff_[" <<i<<"]= "<< mceff <<", random "<<random<<endl;
#endif
            if ( random > mceff ) break;
          }

	  // look if a digit on this detector with this wire already exists. 
	  // Boring, but no idea how to do it in aother way... :(
	  list<CsDigit*>::iterator Id;
	  bool found = false;
	  for( Id=myDigits_.begin(); (Id!=myDigits_.end())&&(!found); Id++ ) {
	    if( i == (*Id)->getAddress() ) {
	      // Here it is, add this hit to it...
	      found = true;
	      dynamic_cast<CsMCDigit*>(*Id)->addHit( *(*Ih) );
	      //double tdc = t + tRes_ * CsRandom::gauss();
	      double tdc = t -
		// Trigger time (and hence offset) is subtracted from hit times
		triggerOffset;
	      if( tdc < (*Id)->getDatum() ) {
		(*Id)->replaceDatum( tdc );
	      }
	    }
	  }
	  if( !found ) { // No digits found with these wire and detector
	    if( i<0 || i>=nWir_ ) {
	      ostringstream ost;
	      ost << "Unreliable wire number: " << i
		  << " (0," << nWir_ <<"), "
		  << " detector : " << GetID() 
		  << " " << unit_
		  << " " << type_ << ".";
	      CsErrLog::Instance()->mes( elAnomaly, ost.str() );
	    }
	    else {
	      //double tdc = t + tRes_ * CsRandom::gauss();
	      double data[2] = {t,0.};
	      CsDigit* digit = new CsMCDigit( *this, i, data, 2 );
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

void CsMW1Detector::DecodeChipDigit(const CS::Chip::Digit &digit)
{
  const CS::ChipF1::Digit *d = dynamic_cast<const CS::ChipF1::Digit *>(&digit);
  if( d==NULL )
    throw CS::Exception("CsMW1Detector::DecodeRawData(): Wrong digit type");
  
  int    wire = d->GetChannel();
  double time = d->GetTimeDecoded()/d->GetTimeUnit(); // just to have histograms in timeslices

  if(hist_!=NULL) hist_->Fill(double(wire),time);

  if (useCalib_) {    
    if (badChannels_.find(wire)!=badChannels_.end()) return; // Discard bad channel
    time -= calib_data[0];                                   // Apply calibration
  }
  if(histCalib_ != NULL) histCalib_->Fill(double(wire),double(time));

//  if (time<hittMin_ || time>hittMax_) return;  // ***** CUT on Hit Time *****

  time *= d->GetTimeUnit(); // convert "bins" back to "ns"
  double data[2] = {time,(double)d->GetChannelPos()};

 // Create CORAL digit
  myDigits_.push_back( new CsDigit(*this, wire, data, 2) );
}

////////////////////////////////////////////////////////////////////////////////

void CsMW1Detector::clusterize() {

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
        cout << "CsMW1Detector::clusterize(): cluster seed is " << firstwire
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
                cout << "CsMW1Detector::clusterize(): next wire is " 
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
        cout << "CsMW1Detector::clusterize(): cluster = (" << wire
             << "," << nwires << ")" << endl;
#endif      

	// ***** SET THE VALUES IN the WRS *****
	double u;	   // u: perpendicular to the detector wires direction
	if (CsInit::Instance()->IsAMonteCarloJob())
	  u = wireDCorr + wire * wirP_;
        else
	  u = wireDCorr + Wire2Pos(wire);
        double v = 0;	   // v: parallel to the detector wires direction
        double w = zcm_;   // w: = Z coordinate in MRS      

	// Set errors:
        HepMatrix cov(3,3,0);  // Zero matrix
        cov(1,1) = pow( nwires * wirP_ / sqrt(12.0), 2 ); // wire pitch / sqrt(12)
        if (fabs(ang_-90)<1) {           // wire length/2
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
        //MAXIM test
	//if (ClusterTime<cltMin_ || ClusterTime>cltMax_) {
	//  delete cluster; continue;
	//}
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

      // ***** SET THE VALUES IN the WRS *****
      double u;	   // u: perpendicular to the detector wires direction
      if (CsInit::Instance()->IsAMonteCarloJob())
	u = wireDCorr + wire * wirP_;
      else
	u = wireDCorr + Wire2Pos(wire);
      double v = 0;	 // v: parallel to the detector wires direction
      double w = zcm_;   // w: = Z coordinate in MRS      

      // Set errors:
      HepMatrix cov(3,3,0);  // Zero matrix
      cov(1,1) = pow( wirP_ / sqrt(12.0), 2 ); // wire pitch / sqrt(12)
      if (fabs(ang_-90)<1) {           // wire length/2
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

void CsMW1Detector::readCalibration(time_t timePoint) {

 // return; // not constants are needed for this detectors (yet)

  CDB::Time tp(timePoint,0);

  string strdata("");
  cdb_->read(GetTBName(),strdata,tp);
  istringstream istrdata( strdata );
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

  if ((getenv("COND_DB_DEBUG")!=0) || (getenv("COND_DB_DEBUG_MW1")!=0)) {
    cout << GetTBName() << endl;
    for(unsigned i=0;i<calib_data.size();i++){
      cout << i << " " << calib_data[i] << endl;
    }
  }
}



//_________________________________________
// redefine virtual Wire2Pos of CsDetector
double CsMW1Detector::Wire2Pos(double wire) const
{
  int w = int(wire);
  const int ss=8;
  int nmod = w / ss;

//=== Commented by Hugo PEREIRA on Wed Oct 16 2002 to remove hard coded pitch and use det.dat pitch instead
//=== MA are done so that every 8 wires are spaced by the pitch, then there is a 0.5*pitch gap between 8th and 9th wire
//   double x = (wirD_)/10. + double(nmod)*(double(ss)+0.5) + double(w % ss);
//   return 10*x;

  double x = double(nmod)*(double(ss)+0.5) + double(w % ss);
  return wirP_*x;
}


//LS Eff
void CsMW1Detector::readMCEffMaps(time_t timePoint)
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
        throw std::logic_error("CsMW1Detector_MCEff::read(): bad line structure.");          
      MCEff_.push_back(mcEf);                                                                       
      MCEff_err_.push_back(mcEf_err);                                                               
                                                                                                    
    }                                                                                               
}             
