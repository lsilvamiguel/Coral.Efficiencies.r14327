// $Id: CsRICH1Detector.cc,v 1.45 2010/05/04 12:45:59 tnagel Exp $

/*!
   \file    CsRICH1Detector.cc
   \brief   Compass Fiber Hodoscope like detector Class.
   \author  Benigno Gobbo
   \version $Revision: 1.45 $
   \date    $Date: 2010/05/04 12:45:59 $
*/

#include "CsStopwatch.h"
#include "CsRICH1Detector.h"
#include "CsMCRICH1Hit.h"
#include "CsMCTrack.h"

#include "CsRCDetectors.h"
#include "CsRCPhotonDet.h"
#include "CsRCCathode.h"
#include "CsRCMirrors.h"
#include "../rich/CsRCRecConst.h"

//#include "CsRCEventAnalysis.h"

#include <cmath>
#include <cstring>
#include <strings.h>
#include <functional>
#include "CsZone.h"
#include "CsDigit.h"
#include "CsMCDigit.h"
#include "CsOpt.h"
#include "CsEvent.h"
#include "CsDigit.h"
#include "CsComgNtCommons.h"
#include "CsGeom.h"
#include "CsGeant3.h"
#include <cstdlib>
#include "DaqDataDecoding/ChipGassiplex.h"
#include "CDB.h"
#include <sstream>
#include <fstream>
#include <cstdlib>

using namespace std;
using namespace CLHEP;
using CS::DetID;

extern QhitType Qhit;


CsRICH1Detector::CsRICH1Detector( const int id, const string &TBname ) : CsDet(DetID("RIC1",id), TBname) {

  _borasOK = false;

  nPhotChmb_ = 0;

  decodingDone_ = false;
  decode_ = false;

  calib_flag = false;  

/* Andrea ... setting flag for database reading of PMT T0 Calibrations */

  _calibT0_flag = false;

/* Andrea ... end setting flag for database reading of PMT T0 Calibrations */

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
	if( *Is == "BM" || *Is == "RICH1" || *Is == "all" ) {
	  decodeCard_ = true;
	}
      }
    }      
    else {
      decodeCard_ = true;
    }
  }
}

bool CsRICH1Detector::operator==( const CsRICH1Detector& det ) const
{
  return GetID()==det.GetID();
}

double CsRICH1Detector::getZcm() const {
  list<CsRCPhotonDet*> lPhoDet = CsRCDetectors::Instance()->lPhoDet();
  list<CsRCPhotonDet*>::iterator id;
  id = lPhoDet.begin();
  Hep3Vector vDet = (*id)->vDet0();
  return( vDet.z()+1500. );            //   provisional !!!
}

void CsRICH1Detector::print() {};

void CsRICH1Detector::addMCHit( CsMCHit& hit ) {
  myMCHits_.push_back( &hit );
}

void CsRICH1Detector::clearMCHitList() {
  myMCHits_.clear();
  if( !myDigits_.empty() ) {
    list<CsDigit*>::iterator id;
    //    for( id=myDigits_.begin(); id!=myDigits_.end(); id++ ) delete *id;
    myDigits_.clear();
  }
  decodingDone_ = false;
  decode_ = false;
}

//===========================================================================
void CsRICH1Detector::makeMCDecoding() {
//--------------------------------------


  // should I proceed?
  if( !decode_ && !decodeCard_ ) return;

  // Already done?
  if( decodingDone_ ) return;

  // clear
  myDigits_.clear();

  // get a link to relevant pieces...
  CsEvent* event = CsEvent::Instance();

  // this can be done only with zebra binary files... 

  // loop on hits
  if( !CsGeant3::Instance()->isAnNtFile() ) {

    getMCDigits( myMCHits_ );
//  ------------------------

  }
  decodingDone_ = true;
}

//===========================================================================
  void  CsRICH1Detector::getPadADR( CsDigit& digit, int& Cathode,
                                    int& ix, int& iy ) const             {
//------------------------------------------------------------------------

    int address = digit.getAddress();
    ix = address & 0x3FF;
    iy = (address >>10)& 0x3FF;
    Cathode =  address >>20;
    return;
  }


  extern "C"  {
    void digits_( const int&, const int*,
                  const float*, const float*, const float*,
// 1900                  int&, int*, int*, int*, float* );
// 0104                  int&, int*, int*, int*, float*, int* );
                  int&, int*, int*, int*, float*, long* );
    void initall_( const int&, const float&,
                   const int&, const float&, const int&,
                   const int&, const int&,
                   const float*, const float*, const float*, const float*,
                   const int&, const int&, const int&,
                   const float&, const float&, const float&, const float&,
                   const float&, const int&);
    void geomet_();
    void pedsigma_();
  }

//===========================================================================
  void CsRICH1Detector::getMCDigits( list<CsMCHit*>& richHits )   {
//-----------------------------------------------------------------

    //cout << " Enter CsRICH1Detector::getMCDigits " << endl;


//--- Digitization of RICH MC hits.
//    -----------------------------
//--- Paolo  -  16 August 2000

    static double richTime;
    int chr1;

    static int nPadx;
    static float padx;
    static int nPady;
    static float pady;
    static int nCathode;
    static int kDoPhoGen = 0;
    static int kNoise = 0;
    static bool firstCall = true;
    static float timegate[16];
    if( firstCall ) {
      firstCall = false;

      richTime = 0.;

      list<CsRCCathode*> lCathodes = CsRCDetectors::Instance()->lCathodes();
      if( lCathodes.size() == 0 ) {
        string str = "CsRICH1Detector::getMCDigits() : wrong detectors.dat for RICH1";
        CsErrLog::Instance()->mes( elFatal, str );
      }
      nPadx = lCathodes.front()->nPadx();
      padx = lCathodes.front()->padx();
      nPady = lCathodes.front()->nPady();
      pady = lCathodes.front()->pady();
      nCathode = CsRCDetectors::Instance()->nCathode();

      CsOpt* opt = CsOpt::Instance();
      bool boo;
      string str;
      vector<float> vec;
      int kin = 0;
      float flo;
      int k = 0;

      bool lprint = false;
      boo = opt->CsOpt::getOpt( "RICHONE", "PrintMCDig", kin );
      if( boo ) { if( kin > 0 ) lprint = true; }

      int kboo = 0;
      int kDoPhoGen = 1;
      boo = opt->CsOpt::getOpt( "RICHONE", "DoPhoGen", str );
      if( boo ) { if( str[0] == 'Y' ) kDoPhoGen = 0; kboo++; }
      int kNoise = 1;
      boo = opt->CsOpt::getOpt( "RICHONE", "ElecNoise", str );
      if( boo ) { 
        if( str[0] == 'Y' ) kNoise =  0;
        if( str[0] == 'O' ) kNoise = -1;
	kboo++;
      }
      float threshin = 0.;
      boo = opt->CsOpt::getOpt( "RICHONE", "Threshold", flo );
      if( boo ) { threshin = flo; kboo++; }

      int dofeedbackin = 0;
      boo = opt->CsOpt::getOpt( "RICHONE", "DoFeedback", str );
      if( boo ) { if( str[0] == 'Y' ) dofeedbackin = 1; kboo++; }
      int doreflin = 0;
      boo = opt->CsOpt::getOpt( "RICHONE", "DoReflection", str );
      if( boo ) { if( str[0] == 'Y' ) doreflin = 1; kboo++; }
      int doQEtestin = 0;
      boo = opt->CsOpt::getOpt( "RICHONE", "DoQEtest", str );
      if( boo ) { if( str[0] == 'Y' ) doQEtestin= 1; kboo++; }
      ;
      float pedmin = 0.;
      float sigpin = 0.;
      float sigmin = 0.;
      float sigsin = 0.;
      boo = opt->CsOpt::getOpt( "RICHONE", "PedSigma", vec );
      if( boo ) {
        pedmin = vec[0];
        sigpin = vec[1];
        sigmin = vec[2];
        sigsin = vec[3];
	kboo++;
      }
      int dosigcompin = 1;
      boo = opt->CsOpt::getOpt( "RICHONE", "DoNSigComp", str ); //031005
      if( boo ) { if( str[0] == 'N' ) dosigcompin = 0; kboo++; }

      float csiqein[25];
      for( int k=0; k<25; k++ ) csiqein[k] = 0.;           //   020619
      boo = opt->CsOpt::getOpt( "RICHONE", "CsIQuantumEff", vec );
      if( boo ) {
        for( k=0; k<25; k++ ) {
          csiqein[k] = vec[k];
        }
	kboo++;
      }
      float alfain[16];
      for( int k=0; k<16; k++ ) alfain[k] = 0.03;          //   020619
      boo = opt->CsOpt::getOpt( "RICHONE", "Alfa", vec );
      if( boo ) {
        for( k=0; k<16; k++ ) {
          alfain[k] = vec[k];
        }
	kboo++;
      }
      float phsloin[16];
      for( int k=0; k<16; k++ ) phsloin[k] = 30.;          //   020619
      boo = opt->CsOpt::getOpt( "RICHONE", "PHslo", vec );
      if( boo ) {
        for( k=0; k<16; k++ ) {
          phsloin[k] = vec[k];
        }
	kboo++;
      }
      float ridfactin[16];
      for( int k=0; k<16; k++ ) ridfactin[k] = 1.;          //   020619
      boo = opt->CsOpt::getOpt( "RICHONE", "RedFactor", vec );
      if( boo ) {
        for( k=0; k<16; k++ ) {
          ridfactin[k] = vec[k];
        }
	kboo++;
      }
      for( int k=0; k<16; k++ ) timegate[k] = 2000.;     //031005
      boo = opt->CsOpt::getOpt( "RICHONE", "TimeGate", vec );
      if( boo ) {
	for( k=0; k<16; k++ ) {
	  timegate[k] = vec[k];
	}
      	kboo++;
      }     
      //      cout << "kboo = " << kboo << endl;

      if( kboo < 13 ) {
	//CsErrLog::mes( elFatal,
	//  "CsRICH1Detector::getMCDigits() : NOT initialized !" );
	CsErrLog::mes( elError,
	  "CsRICH1Detector::getMCDigits() : default initialization !" );
      }

      if( lprint ) {
        cout << endl;
        cout << "CsRICH1Detector::getMCDigits() : MC Digitization Constants"
             << endl
             << "----------------------------------------------------------"
             << endl;
        cout << "kDoPhoGen = " << kDoPhoGen 
             << ", kNoise = " << kNoise
             << ",  theshold = " << threshin << endl;
        cout << "dofeedback = " << dofeedbackin
             << ",  dorefl = " << doreflin 
             << ",  doQEtest = " << doQEtestin << endl;
        cout << "Peds & Sigmas = " << pedmin << ",  " << sigpin
             << ",  " << sigmin << ",  " << sigsin << endl;
        cout << "dosigcomp = " << dosigcompin<<endl;
        cout << "CsIQEff = ";
        for( k=0; k<25; k++ ) { cout << csiqein[k] << "  "; }
        cout << endl;
        cout << "Alfa = ";
        for( k=0; k<16; k++ ) { cout << alfain[k] << "  "; }
        cout << endl;
        cout << "PHslo = ";
        for( k=0; k<16; k++ ) { cout << phsloin[k] << "  "; }
        cout << endl;
        cout << "RedFactor = ";
        for( k=0; k<16; k++ ) { cout << ridfactin[k] << "  "; }
        cout << endl;
        cout << "TimeGate = ";
        for( k=0; k<16; k++ ) { cout << timegate[k] << "  "; }
        cout << endl;
      }

      initall_( nPadx, padx, nPady, pady, nCathode,
//-------------------------------------------------
		kDoPhoGen, kNoise,
		csiqein, alfain, phsloin, ridfactin,
		dofeedbackin, doreflin, doQEtestin,
      //                nwiresin, nwirepadin, pitchin,
		pedmin, sigpin, sigmin, sigsin, threshin,
                dosigcompin );

      geomet_();
//-------------
      pedsigma_();
//---------------
    }
    CsStopwatch chronos( 1 );
    chr1 = chronos.start();

    int nHit = richHits.size();
    int  icHit[nHit];
    float xHit[nHit];
    float yHit[nHit];
    float eHit[nHit];
    CsMCHit* pHit[nHit];
//- coordinates in DRS
//  ------------------
//- WARNING : hit cathodes from 1 to 16, as from COMGeant
//  -----------------------------------------------------
    //    list<CsMCHit*> lHits;

    CsMCTrack* mctrackSel = NULL;
    int kHit = 0;
    list<CsMCHit*>::iterator ih;
    for( ih=richHits.begin(); ih!=richHits.end(); ih++ ) {

      bool bOnePart = false;
//-------------------------------------------------------------------
//--- select the photons (=digits) of only one particle (pion)
//    for background studies (9/01) :
//    use together with the corresponding CsRCEventPads.cc code
//    ---------------------------------------------------------
      if( bOnePart ) {
	CsMCTrack* mctrack = (*ih)->getMCTrack();
	int kType = mctrack->getParticle()->getGeantNumber();
	if( kType >= 200 ) continue;
	if( (*ih)->getOrigin() != 0 ) continue;
	if( kType != 8  &&  kType != 9 ) continue;
	if( mctrackSel == NULL ) mctrackSel = mctrack;
	if( mctrack != mctrackSel ) continue;
      }
//-------------------------------------------------------------------

      CsMCRICH1Hit* hit = (CsMCRICH1Hit*)(*ih);
      int icth=hit->getCathode();
      float htim = hit->getDTime();      
      if( fabs( htim ) > timegate[icth] ) continue;

      icHit[kHit] = icth;
      xHit[kHit] = hit->getYdetDRS();
      yHit[kHit] = hit->getZdetDRS();
      eHit[kHit] = hit->getPhotEnergy();
      pHit[kHit] = (*ih);
      kHit++;
      //      lHits.push_back( hit );

//-------------------------------------------------------------------
      //if( bOnePart ) break;         // NOISE ONLY!; nHit must be > 0!
//-------------------------------------------------------------------
    }
    //    CsRCEventAnalysis::Instance()->hitDisplay( 5000, lHits );
    nHit = kHit;

    static int mDigit = 5000;
    int icto[mDigit];
    int ixto[mDigit];
    int iyto[mDigit];
    float phto[mDigit];
    // 0104    int ihcco[mDigit];         //   1900
    long ihcco[mDigit];         //   0104

    int nDigit = 0;
    //    nDigit=1; icto[0]=1; ixto[0]=1; iyto[0]=1; phto[0]=1.; ihcco[0]=1;

    if( nHit > 0 ) {
      digits_( nHit, icHit, xHit, yHit, eHit,
//-------------------------------------------
	       nDigit, icto, ixto, iyto, phto, ihcco );
    //cout << "getMCDigits, Hits : " << nHit << "  Digits " << nDigit << endl;
    }


    int k;
    CsEvent* event = CsEvent::Instance();
    int kDig = 0;
    for ( kDig=0; kDig < nDigit; kDig++ ) {

      int ic = icto[kDig];
      int ix = ixto[kDig];
      int iy = iyto[kDig];
      double ph = phto[kDig];

      if( ix < 1 || ix > nPadx || iy < 1 || iy > nPady || 
          ic < 1 || ic > nCathode ) {
	ostringstream ost;
        ost << "Unreliable pad number (cathode,ix,iy): " 
	    << ic <<"  "<< ix <<"  "<< iy 
	    << " detector : " << GetID();
	CsErrLog::Instance()->mes( elWarning, ost.str() );
      }
      else {
        CsDigit* digit = new CsMCDigit( *this, setPadADR( ic, ix, iy ), &ph );
        myDigits_.push_back( digit );

        // added 1900
        int MAXPAD = 5000;                 //   provisional !!!
        int maxw = MAXPAD*MAXPAD;
        int iHitco = ihcco[kDig];
	//cout << iHitco << " --- " << nDigit << "  " << nHit << endl;
	//if( iHitco < 0. ) cout << iHitco << " --------------------------- " 
	//	       << nDigit << "  " << kDig << endl;
        if( iHitco < 0. ) continue;        //   provisional !!! 0104
        do {
	  CsMCHit* pHitsv = NULL;
	  //cout << iHitco << "   " << maxw << endl; 
          int iHit = iHitco / maxw;
          iHitco = iHitco - iHit*maxw;
          maxw /= MAXPAD; 
          if( iHit == 0 ) continue;
          if( iHit > nHit ) {
            ostringstream ost;
            ost << "Wrong hit reference : Digit " << kDig
	        << "  nHit  " << nHit << "  Hit  " << iHit;
	    CsErrLog::Instance()->mes( elWarning, ost.str() );
	    //cout << "@@@@@@  nHit  " << nHit << "  Hit  " << iHit
	    //     << "  Digit  " << kDig << endl;
	    continue;
	  }
	  iHit--;
	  //cout << iHitco << "  Digit  " << kDig 
	  //     << "  iHit  " << iHit << endl;

          //int kHit = 1;                      //   kHit from 1 to n !
          //list<CsMCHit*>::iterator ihsv = 0;
          //list<CsMCHit*>::iterator ih;
          //for( ih=richHits.begin(); ih!=richHits.end(); ih++ ) {
          //  if( kHit == iHit ) { ihsv = ih;  break; }
          //  kHit++;
          //}
          //if( ihsv != 0 ) {
          //  dynamic_cast<CsMCDigit*>(digit)->addHit( *(*ihsv) );
          //}

	  pHitsv = pHit[iHit];
	  //cout << iHit << "  ";
          if( pHitsv != NULL ) {
            dynamic_cast<CsMCDigit*>(digit)->addHit( *(pHitsv) );
	    //}
	    //cout << pHitsv->getMCTrack()->getParticle()->getGeantNumber()
	    //     << endl;
	    //cout << pHitsv->getDTime() << endl;
	  }
        } while( maxw >= 1. );
	//cout << dynamic_cast<CsMCDigit*>(digit)->getHits().size() << endl;
	//list<CsMCHit*> lDigHits= dynamic_cast<CsMCDigit*>(digit)->getHits();
	//list<CsMCHit*>::iterator ih;
	//for( ih=lDigHits.begin(); ih!=lDigHits.end(); ih++ ) {
	//std::cout << (*ih)->getDTime() << "  ";
	//}
	//std::cout << std::endl;
      }
    }


    double time2 = chronos.stop( chr1 ) * 1000.;
    richTime += time2;
    //cout << "rich digit time = " << richTime << endl;

    //cout << " Exit CsRICH1Detector::getMCDigits " << endl;

  }

//===========================================================================
  void CsRICH1Detector::setPhotonDet( string name, 
              double xDet0, double yDet0, double zDet0, 
	      HepMatrix rotMatrix ) {
//-----------------------------------

    CsRCDetectors* det = CsRCDetectors::Instance();
    det->setPhotonDet( name, xDet0, yDet0, zDet0, rotMatrix );

  }

//===========================================================================
  void CsRICH1Detector::setCathode( int id, string TBname, string name,
              double xOffCat0, double yOffCat0, double zOffCat0,
              int nPadx, int nPady, double padx, double pady,
	      double ddQzW, double ddGap, HepMatrix corrRotMx ) {
//---------------------------------------------------------------

    CsRCDetectors* det = CsRCDetectors::Instance();
    det->setCathode( id, TBname, name, xOffCat0, yOffCat0, zOffCat0,
                     nPadx, nPady, padx, pady, ddQzW, ddGap,
		     corrRotMx );
  }

//===========================================================================
  void CsRICH1Detector::setMirrNom( string name,
                                    double xC0, double yC0, double zC0,
                                    double RR ) {
//-----------------------------------------------

    CsRCMirrors* mirr = CsRCMirrors::Instance();
    mirr->setMirrNom( name, xC0, yC0, zC0, RR );
  }

//===========================================================================
  void CsRICH1Detector::setMirrEle( string name, 
                                    double theta, double phi, double RR,
                                    double deTheta, double dePhi,
                                    double delta, double qfact, int align ) {
//---------------------------------------------------------------------------

    CsRCMirrors* mirr = CsRCMirrors::Instance();
    mirr->setMirrEle( name, theta, phi, RR, deTheta, dePhi,
                      delta, qfact, align );
  }

//===========================================================================
  void CsRICH1Detector::setMCCFRefInd( double index ) {
//-----------------------------------------------------

      CsRCRecConst::Instance()->setMCCFRefInd( index );
  }

//===========================================================================
  void CsRICH1Detector::setMCCFRefIndVS( double index ) {
//-----------------------------------------------------

      CsRCRecConst::Instance()->setMCCFRefIndVS( index );
  }


//===========================================================================
  void CsRICH1Detector::DecodeChipDigits(const CS::Chip::Digits &digits)
{
  //^cout << "CsRICH1Detector::DecodeChipDigits " << digits.size() << endl;
  CsDet::DecodeChipDigits( digits );
  setDecodingDone();
}
//===========================================================================
  void CsRICH1Detector::DecodeChipDigit(const CS::Chip::Digit &digit)
{
  //^cout << "CsRICH1Detector::DecodeChipDigit" << endl;
  const CS::ChipGassiplex::Digit *d = dynamic_cast<const CS::ChipGassiplex::Digit *>(&digit);
  if( d==NULL )
    throw CS::Exception("CsRICH1Detector::DecodeRawData(): Wrong digit type");

//   //
//   // Apply calibrations
//   //
//   if(calib_data.size() != 0)
//   {
//     if(d->GetChannel()<0 || d->GetChannel() >= int(calib_data.size()) ) {
// 	cout<<"CsFiberHodoDetector::DecodeRawData() ==> SiFi Id = "
// 	    <<GetID()<<" 'wire #' "<<d->GetChannel()<<"  \t but calib. container size = "
// 	    <<calib_data.size()<< endl;
//       } else {
// 	time -= calib_data[d->GetChannel()]; // apply calibration
//       }
// 
//       if(mH2[4]!=NULL)
//         mH2[4]->Fill(time,d->GetChannel());
//   }

    //
    // Create CORAL digit
    //

    //------- CHECK THIS ------
    //    d->Print(cout,"RICH1 digit: ");

  //  if( d->GetCathode() == 13  &&
  //      d->GetX() < 20  &&  d->GetY() < 20 ) {

  //  OLd code
//     double data[4]={d->GetCathode(),d->GetX(),d->GetY(),d->GetAmplitude()};
//     myDigits_.push_back( new CsDigit(*this, 0, data, 4) );
  //  New code

    double amp = (double)d->GetAmplitude();
    int adr = setPadADR( d->GetCathode()+1, d->GetX()+1, d->GetY()+1 );
    CsDigit* rich_digit = new CsDigit( *this, adr, &amp );
    myDigits_.push_back( rich_digit );
    //^cout << "CsRICH1Detector::DecodeChipDigit" << myDigits_.size() << endl;
    //------- CHECK THIS ------

    //  }
    //cout<< this->getDTime() << endl;
    //cout<<"CsRICH1Detector::ReadCalib() ==> "<<GetTBName()<<" calibrations, valid for ";
    //cout<<t.tm_mday<<"."<<t.tm_mon+1<<"."<<t.tm_year+1900<<" "
    //  <<t.tm_hour<<":"<<t.tm_min<<":"<<t.tm_sec
    //  <<", not found in DB"<<endl;

}


//===========================================================================
void CsRICH1Detector::readCalibration(time_t timePoint) {

//append three underscore to rich1 TBName
//   size_t npo = folderSet_.find( GetTBName(), 0 );
//   folderSet_ = folderSet_.erase( npo+5, 3 );
//   folderSet_ = folderSet_.append( "___", 0, 3 );
//   string folder = folderSet_+"/calibration";
  //cout << "======= " << folderSet_ << endl;

  string tbname_(GetTBName());

  tbname_.erase( 5, 3 );
  tbname_.append( "___", 0, 3 );

  CDB::Time tp(timePoint,0);

  string strdata("");
  cdb_->read(tbname_,strdata,tp);
  istringstream istrdata(strdata);
  istrdata >> calib_data;

  if (calib_data.size()!=0) {
    CsErrLog::msg(elInfo,__FILE__,__LINE__,
		  "%s calibration OK.",GetTBName().c_str());

    calib_flag = true;
    index_ = calib_data[0];
    indexUV_ = index_;
    thresh_ = calib_data[1];
    indexVS_ = 0.;
    if (calib_data.size()>2) indexVS_ = calib_data[2];
    std::cout << setprecision( 6 );
    cout << "RICH1 refractive index from DB " << index_ << " ,  ";
    if (calib_data.size()>2) cout << indexVS_ << " ,  ";
    std::cout << setprecision( 2 );
    cout << "pad global software threshold " << thresh_ << endl;

  }else{
    tm *t = localtime(&tp.first);
    cout << GetTBName() << ", no calibration for local time "
	 << t <<" in CDB"<< endl;
  }

  std::cout << "Read RICH1 event1 data. Temporarly from a file." << std::endl;
  int run_nb = CsEvent::Instance()->getRunNumber();
  bool status = readBorasTables( run_nb );

  /* Andrea ... reading of PMT T0 Calibrations */
  // Initialization
  for ( int ge=0;ge<160;ge++ ) {
    for ( int ch=0;ch<  8;ch++ ) {
      for ( int an=0;an<  8;an++ ) {
	_pmt_t0calarray[ge][ch][an] = 0;
      }
    }  
  }

  /*
  CsRCDetectors   * dets = CsRCDetectors::Instance();
  list<CsRCCathode*> lCathodes = dets->lCathodes();
  int nCathode = lCathodes.size();
   for( int ic=0; ic< lCathodes.size(); ic++ ) {
     CsRCCathode* catP = dets->ptrToCat( ic );
     cout << " " << catP->TBname() << endl;
     string pmtbname_(GetTBName());
     pmtbname_.replace(1,1,"M");
     pmtbname_.erase( 5, 3 );
     pmtbname_.append( "___", 0, 3 );
   }
  */

  string pmtbname_(GetTBName());
  pmtbname_.replace(1,1,"M");
  pmtbname_.erase( 5, 3 );
  pmtbname_.append( "___", 0, 3 );

  CDB::Time pmtp(timePoint,0);

  string pmstrdata("");
  cdb_->read(pmtbname_,pmstrdata,pmtp,"T0");
  istringstream pmistrdata(pmstrdata);

  _calibT0_flag = true;

  std::string aa,bb,cc,dd;
  pmistrdata >> aa >> bb >> cc >> dd >> dd >> dd >> dd >> dd >> calibT0_data;
  /* SKIPPING FIRST VALUES (COMMENTS FROM INPUT FILE */
  //  cout << " AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA " <<  dd << " AAAA " << calibT0_data[0] << endl;

  if (calibT0_data.size()!=0) {
    CsErrLog::msg(elInfo,__FILE__,__LINE__,
		  "%s calibration OK.",GetTBName().c_str());

    _calibT0_flag = true;
    for (unsigned int kk=0; kk < calibT0_data.size(); kk=kk+5 ) {
      int srcID = calibT0_data[0+kk];
      int port  = calibT0_data[1+kk];
      int GeoID = (srcID-500)*16+port;
      int ChiID = calibT0_data[2+kk];
      int ChaID = calibT0_data[3+kk];
      if ( GeoID >= 16 && GeoID < 260 && ChiID >= 0 && ChiID <8 && ChaID >=0 && ChaID < 8 ) { 
	_pmt_t0calarray[GeoID][ChiID][ChaID] = calibT0_data[4+kk];
      } else {
	std::cout << " CsRICH1Detector: WARNING xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx " << endl;
	std::cout << " CsRICH1Detector: bad reading of PMT T0 calibrations - GeoID "  << GeoID << " ChipID " << ChiID 
		  << " ChannelID " << ChaID << endl ; 
      }

      CsOpt* opt = CsOpt::Instance();
      int kin = 0;
      bool lprint = false;
      bool boo = opt->CsOpt::getOpt( "RICHONE", "kPrintRichRec", kin );
      if( boo ) { if( kin > 0 ) lprint = true; }

      if(lprint){
	std::cout << setprecision( 6 );
	std::cout << "MaPMT T0 calib from DB: srcID " << srcID << " Port " << port << " GeoID " << GeoID << " ChipID " << ChiID 
		  << " ChannelID " << ChaID << " T0 " << _pmt_t0calarray[GeoID][ChiID][ChaID] << endl;
      }
      
    }

  }else{
    tm *t = localtime(&pmtp.first);
    cout << GetTBName() << ", ReadT0calib - no calibration for local time "
	 << t <<" in CDB"<< endl;
  }

  /* Andrea ... end of reading of PMT T0 Calibrations */

}

 /* Andrea ... way to access private array for calibrations */

int CsRICH1Detector::getPMT_T0( int geoid, int chipid, int chanid ) {
  if( geoid >=0 && geoid <160 && chipid>=0 && chipid<8 && chanid>=0 && chanid <8 ) {
    return( _pmt_t0calarray[geoid][chipid][chanid] );
  }
  else {
    return( -1 );
  }
}

 /* Andrea ... end of way to access private array for calibrations */

int CsRICH1Detector::getBoraVoltage( int geoid, int i ) {
  if( geoid >=0 && geoid <248 && i>=0 && i<5 ) {
    return( _voltage[geoid][i] );
  }
  else {
    return( -1 );
  }
}


int CsRICH1Detector::getBoraTemperature( int geoid, int i ) {
  if( geoid >=0 && geoid <248 && i>=0 && i<4 ) {
    return( _temperature[geoid][i] );
  }
  else {
    return( -1 );
  }
}

bool CsRICH1Detector::readBorasTables( int run ) {

  _borasOK = false;

  // Preset to -1
  for( int i=0; i<16; i++ ) {
    for( int j=0; j<72; j++ ) {
      for( int k=0; k<72; k++ ) {
	_threshold[i][j][k] = -1;
	_pedestal[i][j][k]  = -1.;
	_sigma[i][j][k]     = -1.;
      }
    }
  }
  for( int i=0; i<248; i++ ) {
    for( int j=0; j<5; j++ ) {
      _voltage[i][j] = -1;
    }
    for( int j=0; j<4; j++ ) {
      _temperature[i][j] = -1;
    }
  }

  // At the moment data is read from a file...

  std::ostringstream filename;
  char* filepath = getenv( "RICH1_EVENT1_DATA_PATH" );
  if ( filepath == NULL ) {
    filename << "run" << run << ".rich1_event1_data" << std::ends;
  }
  else {
    filename << filepath << "/run" << run << ".rich1_event1_data" << std::ends;
  }
  std::cout << "Trying to open " << filename.str().c_str() << "..." << std::endl;
  std::fstream f( filename.str().c_str(), std::ios::in );

  if( !f ) {
    std::cerr << "+-------------------------------------------------+" << std::endl
	      << " RICH1 EVENT1 Data not found. RICH1 multi-array" << std::endl
	      << " set to (-1). Ignore this message if you do not" << std::endl
	      << " understand what it means..." << std::endl
	      << "+-------------------------------------------------+" << std::endl;
    return false;
  }

  const int linesize = 256;
  char line[linesize];
  // check:
  f.getline( line, linesize, '\n' );
  if( strncmp( line, "#GId ", 5 ) != 0 ) {
    std::cerr << "+-------------------------------------------------+" << std::endl
	      << " RICH1 EVENT1 Data file wrong. RICH1 multi-array" << std::endl
	      << " set to (-1). Ignore this message if you do not" << std::endl
	      << " understand what it means..." << std::endl
	      << "+-------------------------------------------------+" << std::endl;
    return false;
  }

  _borasOK = true;

  do {
    f.getline( line, linesize, '\n' );
    if( f.eof() || strncmp( line, "# c ", 4 ) == 0 ) continue;
    std::istringstream s(line);
    int geoid, i;
    s >> geoid;
    if( geoid < 0 || geoid >= 248 ) { 
      std::cerr << "+----------------------------------------------------+" << std::endl
		<< " RICH1 EVENT1 Data file corrupted. RICH1 multi-array" << std::endl
		<< " partially set to (-1). Ignore this message if you" << std::endl
		<< " do not understand what it means..." << std::endl
		<< "+----------------------------------------------------+" << std::endl;
      return false;
    }
    i = 0;
    s >> _voltage[geoid][i++];
    s >> _voltage[geoid][i++];
    s >> _voltage[geoid][i++];
    s >> _voltage[geoid][i++];
    s >> _voltage[geoid][i++];
    i = 0;
    s >> _temperature[geoid][i++];
    s >> _temperature[geoid][i++];
    s >> _temperature[geoid][i++];
    s >> _temperature[geoid][i++];
  } while( !f.eof() && strncmp( line, "# c ", 4 ) != 0 );

  if( strncmp( line, "# c ", 4 ) != 0 ) {
    std::cerr << "+-------------------------------------------------+" << std::endl
	      << " RICH1 EVENT1 Data file wrong. RICH1 multi-array" << std::endl
	      << " set to (-1). Ignore this message if you do not" << std::endl
	      << " understand what it means..." << std::endl
	      << "+-------------------------------------------------+" << std::endl;
    return false;
  }

  do {
    f.getline( line, linesize, '\n' );
    if( f.eof() ) continue;
    std::istringstream s(line);
    int i, j, k;
    s >> i >> j >> k;
    if( i<-1 || i>=16 || j<-1 || j>=72 || k<-1  || k>= 72 ) {
      std::cerr << "+----------------------------------------------------+" << std::endl
		<< " RICH1 EVENT1 Data file corrupted. RICH1 multi-array" << std::endl
		<< " partially set to (-1). Ignore this message if you" << std::endl
		<< " do not understand what it means..." << std::endl
		<< "+----------------------------------------------------+" << std::endl;
      return false;
    }
    if( i>=0 && j>=0 && k>=0 ) s >> _threshold[i][j][k] >> _pedestal[i][j][k] >> _sigma[i][j][k];
  } while( !f.eof() );

  std::cout << "+----------------------------------------------------+" << std::endl
	    << " RICH1 EVENT1 Data file for run " << run << std::endl
	    << " successfully read and stored to multiarrays." << std::endl
	    << "+----------------------------------------------------+" << std::endl;
  return true;
}

//===========================================================================

