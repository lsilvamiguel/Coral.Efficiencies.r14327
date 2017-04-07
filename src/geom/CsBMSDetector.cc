// $Id: CsBMSDetector.cc,v 1.50 2010/01/28 12:51:25 tnagel Exp $

/*!
   \file    CsBMSDetector.cc
   \brief   Compass Fiber Hodoscope like detector Class.
   \author  Benigno Gobbo
   \version $Revision: 1.50 $
   \date    $Date: 2010/01/28 12:51:25 $
*/

#include "CsBMSDetector.h"

#include <math.h>
#include <string.h>
#include <strings.h>
#include <functional>
#include "CsZone.h"
#include "CsDigit.h"
#include "CsInit.h"
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

using namespace std;
using namespace CLHEP;

extern QhitType Qhit;

CsBMSDetector::CsBMSDetector( const int    row,
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
	      nWir, wirP, eff, bkg, tGate ), tRes_(tRes)
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
	if( *Is == "BM" || *Is == "BMS" || *Is == "all" ) {
	  decodeCard_ = true;
	}
      }
    }
    else {
      decodeCard_ = true;
    }
  }

  this->use_calib = CsInit::Instance()->useCalibration();
  if(CsOpt::Instance()->getOpt( "", "Monte Carlo job")) this->use_calib = false;

}

//______________________________________________________________________
void CsBMSDetector::BookHistograms( void )
{


  //=== check if histograms are to be booked ===
  CsDetector::ReadHistLevel();
  if( hLevel_ == None ) return;

  CsHistograms::SetCurrentPath("/BMS");

  string tbn  = GetTBName();
  unsigned int id = (unsigned int) GetID();

  char cid[5]; sprintf(cid,"%04u",id); // char[] with this det. id
  string hnam = "BMS_" + string(cid) + "_"; // histograms family name

  //int range_low = -3000;
  int range_low = -7500;
  //int range_high = 1000; // zoom & calibration histogram range
  int range_high = 5000;
  int range_low2 = -500;
  int range_high2 = 500;

  // normal level histograms
  if( hLevel_ >= Normal ) {
    //mH1[0]=new CsHist1D(tbn+"0",tbn+"Time - trig. time",       128,  -65536, 65536);
    //mH1[1]=new CsHist1D(tbn+"1",tbn+"Time - trig. time (zoom)",(range_high-range_low), range_low, range_high);
    mH1[1]=new CsHist1D(tbn+"1",tbn+"Time - trig. time (zoom)",(range_high-range_low)/2, range_low, range_high);
    mH1[2]=new CsHist1D(tbn+"2",tbn+"Hit pattern", nWir_, 0, nWir_);

    if(this->use_calib)
    mH1[6]=new CsHist1D(tbn+"6",tbn+"Time - trig. time (calibrated)",500, -500, 500);

  }

  // high level histograms
  if( hLevel_ >= High ) {
    mH1[0]=new CsHist1D(tbn+"0",tbn+"Time - trig. time",       128,  -65536, 65536);
    mH2[3]=new CsHist2D(tbn+"3",tbn+"Fiber VS time", (range_high-range_low)/2, range_low, range_high, nWir_, 0, nWir_);
    mH2[9]=new CsHist2D(tbn+"9",tbn+"Fiber VS time(zoom)", (range_high2-range_low2)/2, range_low2, range_high2, nWir_, 0, nWir_);
    if(this->use_calib)
    mH2[4]=new CsHist2D(tbn+"4",tbn+"Fiber VS time (calibrated)",200, -200, 200, nWir_, 0, nWir_);
  }

  CsHistograms::SetCurrentPath("/");

  return;
}

//______________________________________________________________________
bool CsBMSDetector::operator==( const CsBMSDetector& det ) const {
  return( CsDetector::operator==(det) );
}

bool CsBMSDetector::operator<( const CsBMSDetector& det ) const {
  return( CsDetector::operator<(det) );
}

void CsBMSDetector::makeMCDecoding() {

  // should I proceed?
  if( !decode_ && !decodeCard_ ) return;

  // Already done?
  if( decodingDone_ ) return;

  // clear
  myDigits_.clear();

  // get a link to relevant pieces...
  CsEvent* event = CsEvent::Instance();

  // this can be done only with zebra binary files...

  if( !CsGeant3::Instance()->isAnNtFile() ) {
    list<CsMCHit*>::iterator Ih;
    for( Ih=myMCHits_.begin(); Ih!=myMCHits_.end(); Ih++ ) { // loop on hits

      CsMCTrkHit* thit = dynamic_cast<CsMCTrkHit*>(*Ih);
      if( thit == 0 ) return;

      //only charged particles
      if(  ((thit->getMCTrack())->getParticle())->getCharge() ) {

	double t  = thit->getDTime();  // Delay time (ns)
	double ui = thit->getUin();    // Hit in point (DRS)
	double vi = thit->getVin();
	double wi = thit->getWin();
	double uo = thit->getUout();   // Hit out point (DRS)
	double vo = thit->getVout();
	double wo = thit->getWout();

	// Check if the hit is inside the time gate
	if ( (-tGate_/2) > t ||  t > (tGate_/2) ) continue;

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
	
	  // look if a digit on this detector with this wire already exists.
	  // Boring, but no idea how to do in other vay... :(
	  list<CsDigit*>::iterator Id;
	  bool found = false;
	  for( Id=myDigits_.begin(); (Id!=myDigits_.end())&&(!found); Id++ ) {
	    if( i == (*Id)->getAddress() ) {
	      // Here it is, add this hit to it...
	      found = true;
	      dynamic_cast<CsMCDigit*>(*Id)->addHit( *(*Ih) );
	      double tdc = t + tRes_ * CsRandom::gauss();
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
	      double tdc = t + tRes_ * CsRandom::gauss();
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

void CsBMSDetector::clusterize() {

  return;   //clusterization of BMS planes exclusively in the BMS reconstruction package

  clearClusterList();

  list<CsDigit*>::iterator Id;
  list<CsDigit*> digits = getMyDigits();

  // protection
  if( digits.empty() ) return;

  int err;
  HepMatrix iRotM(3,3); iRotM = rotWRS_.inverse( err );
  double wireDCorr = wirD_ + iRotM(1,1) * _deltax + iRotM(1,2) *_deltay;

  for( Id=digits.begin(); Id!=digits.end(); Id++ ) {

    int wire = int((*Id)->getAddress());

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
  sortClusters();
  setClusteringDone();
}

////////////////////////////////////////////////////////////////////////////////

void CsBMSDetector::DecodeChipDigit(const CS::Chip::Digit &digit)
{
  const CS::ChipF1::Digit *d = dynamic_cast<const CS::ChipF1::Digit *>(&digit);
  if( d==NULL )
    throw CS::Exception("CsBMSDetector::DecodeRawData(): Wrong digit type");

  if( d->GetChannel() < 0 || d->GetChannel() >= int(getNWir()))
    throw CS::Exception("CsBMSDetector::DecodeRawData() ==> Unexpected wire number %d for SiFi %d %s",
                         d->GetChannel(),GetID().GetNumber(),GetID().GetName().c_str());

  // ***** TIME RELATIVE TO TRIGGER TIME *****
  double time = d->GetTimeDecoded()/d->GetTimeUnit();

  // 'time' value now is in "bins", but it can be real number, because
  // trigger time decoding can give non-integer trigger time

  if( mH1[0]!=NULL )  mH1[0]->Fill(time);
  if( mH1[1]!=NULL )  mH1[1]->Fill(time);
  if( mH1[2]!=NULL )  mH1[2]->Fill(d->GetChannel());
  if( mH2[3]!=NULL )  mH2[3]->Fill(time,d->GetChannel());
  if( mH2[9]!=NULL )  mH2[9]->Fill(time,d->GetChannel());

  // ***** APPLY CALIBRATIONS *****
  bool useIt=false;
  if(calib.size() != 0) {
    if(d->GetChannel()<0 || d->GetChannel() >= int(calib.size()) ) {
      CsErrLog::msg(elError,__FILE__,__LINE__,
		    "BMS %s wire #%d but calib. container size = %d",
		    GetTBName().c_str(),d->GetChannel(),calib.size());
    } else {
      if(calib[d->GetChannel()].flag==0 || calib[d->GetChannel()].flag==2) useIt=true;
      time -= calib[d->GetChannel()].data; // apply calibration
    }
    if( mH2[4]!=NULL )  mH2[4]->Fill(time,d->GetChannel());
    if( mH1[6]!=NULL )  mH1[6]->Fill(time);
  }

  time *= d->GetTimeUnit(); // convert "bins" to "ns"
  //
  // Create CORAL digit
  //
  if(useIt) myDigits_.push_back( new CsDigit(*this, d->GetChannel(), &time, 1) );
  //else cout<<"Rejecting digit!"<<endl;
}


//-------------------------------------------------------------------------

void CsBMSDetector::readCalibration(time_t timePoint){

  CDB::Time tp(timePoint,0);

  calib.clear();
  for (int a=0; a<getNWir(); a++) {
    calib.push_back(Calib(a, 0.0, 17));
  }

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
      if (chn<0 || chn>=getNWir()) {
	cout << "CsBMSDetector::readCalibration() invalid channel number " << chn << "in " << GetTBName() << endl;
      }
    }
  } else {
    tm *t = localtime(&tp.first);
    cout << GetTBName() << ", no calibration for local time "
         << t <<" in CDB"<< endl;
  }
  //cout<<"T0-test: "<<GetTBName().c_str()<<endl;
  //for(int m=0;m<getNWir();m++){
  //cout<<"BMS-ch "<<m<<"t0: "<<calib_data[m]<<" flag: "<<calib_flag[m]<<endl;
  //}
  string debug_level;
  if (getenv("COND_DB_DEBUG")!=0) debug_level = getenv("COND_DB_DEBUG");
  if ((debug_level == "devel") || (getenv("COND_DB_DEBUG_BM")!=0)) {
    cout << GetTBName() << endl;
    for (it=calib.begin(); it!=calib.end(); it++) {
      cout << (*it).chn << " ";
      cout << (*it).data << " ";
      cout << (*it).flag <<endl;
    }
  }
}


//----------------------------------------------------------------------------

istream& operator>>(istream& in, CsBMSDetector::Calib &c) {
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

istream& operator>>(istream& in, vector<CsBMSDetector::Calib> &vc) {
  CsBMSDetector::Calib c;

//  vc.clear();
  while (in >> c) {
    if (((int)vc.size()) < (c.chn+1)) vc.resize(c.chn+1, CsBMSDetector::Calib(0, 0.0, 17));
    vc[c.chn] = c;
  }

  if (in.eof()) in.clear();
  return in;
}


