// $Id: CsRPDetector.cc,v 1.8 2010/01/24 16:10:41 suhl Exp $
/*!
   \file   CsRPDetector.cc
   \brief   Compass Recoil Proton detector Class.
   \author  Benigno Gobbo
   \version $Revision: 1.8 $
   \date    $Date: 2010/01/24 16:10:41 $
*/

#include "CsRPDetector.h"

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

#include "TString.h"
#include "TLorentzVector.h"

using namespace std;
using namespace CLHEP;

// Dymmy func6tions for time being.
// MC digitization and clusterization to be added.

////////////////////////////////////////////////////////////////////////////////

CsRPDetector::CsRPDetector( const int    row,
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
              nWir, wirP, eff, bkg, tGate ) {
  InitRPD();

  string TB = GetTBName().substr(0,2);  // First 2 letter of TB name

  // should this detector be decoded
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
        if( *Is==TB || *Is == "RPD" || *Is == "RP" || *Is == "all" ) {
          decodeCard_ = true;
        }
      }
    } else {
      decodeCard_ = true;
    }
  }

  // add detector IDs to set of ids belonging to this detector
  ids_.clear();
  for (int i=0; i<nWir; i++) {
    ids_.insert(id+i+1);       // slabs detector ID is base detector id + copy number (this starts from 1)
  }

  // check that we did not mess up the set
  if((unsigned) nWir_ != ids_.size())
    CsErrLog::msg(elFatal, __FILE__, __LINE__,
                  "%s: Set of IDs contains more elements than there are slabs on the detector.", GetTBName().c_str());
}

////////////////////////////////////////////////////////////////////////////////

// we need a correct ids_ set for identification of which detector IDs belong
// to this detector in Monte Carlo, therefore we have to make sure the
// AddSubDetector method from CsDetector is not called. Have a look at
// CsTriggerHodoDetector when you want to implement this method for RPD detectors
void CsRPDetector::AddSubDetector( const int    row,
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
  CsErrLog::msg( elFatal, __FILE__, __LINE__,
                 "CsRPDetector::AddSubDetector not implemented!" );
}

////////////////////////////////////////////////////////////////////////////////

void CsRPDetector::BookHistograms() {
  // get level of histograms to be filled from CORAL options file
  ReadHistLevel();
  if( hLevel_ == None ) return;

  string tbn  = GetTBName();
  CsHistograms::SetCurrentPath("/RPDs");

  if( hLevel_ >= High ) {
    // only create these histograms for Monte Carlo data
    if( CsInit::Instance()->IsAMonteCarloJob() )  {
      mH1[tbn+"_E"]=new CsHist1D(tbn+"_E", tbn+" energy deposition (GeV)", 100, 0., 0.1);
      mH1[tbn+"_phi"]=new CsHist1D(tbn+"_phi", tbn+" phi angle", 100, -M_PI, M_PI);

      // there could also be 2d histograms
      // mH2[tbn+"_bla"]=new CsHist2D(tbn+"_bla", tbn+" bla", 100, -0.1, 2.4, 100, -0.1, 2.4);
    }
  }
}

////////////////////////////////////////////////////////////////////////////////

CsRPDetector::~CsRPDetector() {
}

////////////////////////////////////////////////////////////////////////////////

// create CsDigits from Monte Carlo hits
void CsRPDetector::makeMCDecoding() {
  if (!decode_ && !decodeCard_) return;   // Should I proceed?
  if (decodingDone_) return;              // Already done?
  myDigits_.clear();                      // Clear

  // check if we are working on a zebra file
  if( !CsGeant3::Instance()->isAnNtFile() ) {
    list<CsMCHit*>::iterator Ih;
    for( Ih=myMCHits_.begin(); Ih!=myMCHits_.end(); Ih++ ) { // loop on hits
      // get the track hit of current hit
      CsMCTrkHit* hit = dynamic_cast<CsMCTrkHit*>(*Ih);

      if (hit==0) {
        CsErrLog::msg(elError, __FILE__, __LINE__,
                      "CsRPDetector::makeMCDecoding: Monte Carlo hit in RPD is no track hit, ignoring this hit.");
      } else {
        // here the creation of CsMCDigits should happen ...

        // fill the histograms
        string tbn  = GetTBName();
        if (hLevel_ >= High) {
          if (mH1[tbn+"_E"])   mH1[tbn+"_E"]  ->Fill(hit->getELos());
          if (mH1[tbn+"_phi"]) mH1[tbn+"_phi"]->Fill(atan2(hit->getX(), hit->getY()));
        }
      }
    }
  } else {
    CsErrLog::msg(elFatal, __FILE__, __LINE__,
                  "CsRPDetector::makeMCDecoding cannot work on n-tuple input file.");
  }

  decodingDone_ = true;
}

////////////////////////////////////////////////////////////////////////////////

void CsRPDetector::DecodeChipDigits(const CS::Chip::Digits &digits) {
//   bool debug = true;
  bool debug = false;
  if( debug ) cout <<" This is CsRPDetector::DecodeChipDigits function and it is really called in Coral!! " << endl;
  typedef std::multimap<CS::DetID,CS::Chip::Digit*>::const_iterator m_it; // iterator type

  {
    std::string dname("RP01TA__");
    std::pair<m_it,m_it> m_range = digits.equal_range( dname );
    for( m_it d_it=m_range.first; d_it!=m_range.second; d_it++ )
    {
      if( debug ) cout << " decode RP01TA__ ChipF1 digit " << endl;
      const CS::ChipF1::Digit *dt = reinterpret_cast<const CS::ChipF1::Digit *>(d_it->second);
      int32 ch_digit = dt->GetChannel();
      int32 ch_pos_digit = dt->GetChannelPos();
      int32 t_digit = dt->GetTime();
      int32 x_digit = dt->GetX();
      int32 y_digit = dt->GetY();
      int ctr=0;
      TString planeName = ChannelToPMT(dname.data(),int(x_digit),ctr);

      if( debug ) cout << " Store RP01TA__ digit x=" << x_digit << " y=" << y_digit <<
                       " ch_pos_digit =" << ch_pos_digit << " ch_digit=" << ch_digit <<
                       " planeName " << planeName << " ctr " << ctr << endl;
      if (planeName.Contains("Au")&&fDataSentAupt[ctr]<4 && TMath::Abs(t_digit+970.)<40){ // 
      	fTdcAup[ctr][fDataSentAupt[ctr]]=t_digit+970.;
      	fDataSentAupt[ctr]++;
      } 
      if (planeName.Contains("Ad")&&fDataSentAdot[ctr]<4 && TMath::Abs(t_digit+970.)<40){ // 
      	fTdcAdo[ctr][fDataSentAdot[ctr]]=t_digit+970.;
      	fDataSentAdot[ctr]++;
      } 
    }
  }

  {
    std::string dname("RP01TBl__");
    std::pair<m_it,m_it> m_range = digits.equal_range( dname );
    for( m_it d_it=m_range.first; d_it!=m_range.second; d_it++ )
    {
      if( debug ) cout << " decode RP01TBl__ ChipF1 digit " << endl;
      const CS::ChipF1::Digit *dt = reinterpret_cast<const CS::ChipF1::Digit *>(d_it->second);
//      const float F1bin = 0.1289231; // ns
      int32 ch_digit = dt->GetChannel();
      int32 ch_pos_digit = dt->GetChannelPos();
      int32 t_digit = dt->GetTime();
      int32 x_digit = dt->GetX();
      int32 y_digit = dt->GetY();

      int ctr=0;
      TString planeName = ChannelToPMT(dname.data(),int(x_digit),ctr);
      if( debug ) cout << " Store RP01TBl__ digit x=" << x_digit << " y=" << y_digit <<
                       " ch_pos_digit =" << ch_pos_digit << " ch_digit=" << ch_digit <<
                       " planeName " << planeName << " ctr " << ctr << endl;
      if (planeName.Contains("Bu")&&fDataSentBupt[ctr]<4 && TMath::Abs(t_digit+970.)<40){ // 
      	fTdcBup[ctr][fDataSentBupt[ctr]]=t_digit+970.;
      	fDataSentBupt[ctr]++;
      } 
      if (planeName.Contains("Bd")&&fDataSentBdot[ctr]<4 && TMath::Abs(t_digit+970.)<40){ // 
      	fTdcBdo[ctr][fDataSentBdot[ctr]]=t_digit+970.;
      	fDataSentBdot[ctr]++;
      } 
    }
  }

  {
    std::string dname("RP01TBh__");
    std::pair<m_it,m_it> m_range = digits.equal_range( dname );
    for( m_it d_it=m_range.first; d_it!=m_range.second; d_it++ )
    {
      if( debug ) cout << " decode RP01TBh__ ChipF1 digit " << endl;
      const CS::ChipF1::Digit *dt = reinterpret_cast<const CS::ChipF1::Digit *>(d_it->second);
      int32 ch_digit = dt->GetChannel();
      int32 ch_pos_digit = dt->GetChannelPos();
      int32 t_digit = dt->GetTime();
      int32 x_digit = dt->GetX();
      int32 y_digit = dt->GetY();

      int ctr=0;
      TString planeName = ChannelToPMT(dname.data(),int(x_digit),ctr);
      if( debug ) cout << " Store RP01TBh__ digit x=" << x_digit << " y=" << y_digit <<
                       " ch_pos_digit =" << ch_pos_digit << " ch_digit=" << ch_digit <<
                       " planeName " << planeName << " ctr " << ctr << endl;
    }
  }

  {
    std::string dname("RP01Ql__");
    std::pair<m_it,m_it> m_range = digits.equal_range( dname );
    for( m_it d_it=m_range.first; d_it!=m_range.second; d_it++ )
    {
      if( debug ) cout << " decode RP01Ql__ SADC digit " << endl;
      const CS::ChipSADC::Digit *ds = reinterpret_cast<const CS::ChipSADC::Digit *>(d_it->second);
      int32 x_digit = ds->GetX();
      int32 y_digit = ds->GetY();

      int ctr=0;
      TString planeName = ChannelToPMT(dname.data(),int(x_digit),ctr);
      if( debug ) cout << " Store RP01Ql__ SADC digit x=" << x_digit << " y=" << y_digit <<
                       " planeName " << planeName << " ctr " << ctr << endl;
      if (planeName.Contains("Au")){
      	fDataSentAupq[ctr]=true;
      	fAdcAup[ctr]=treatADC(ds);
      }
      if (planeName.Contains("Ad")){
      	fDataSentAdoq[ctr]=true;
      	fAdcAdo[ctr]=treatADC(ds);
      }
      if (planeName.Contains("Bu")){
      	fDataSentBupq[ctr]=true;
      	fAdcBup[ctr]=treatADC(ds);
      }
      if (planeName.Contains("Bd")){
      	fDataSentBdoq[ctr]=true;
      	fAdcBdo[ctr]=treatADC(ds);
      }
    }
  }

}

/////////////////////////////////////////////////////////////////////////////////

TString CsRPDetector::ChannelToPMT(const char *RpdPlanes, int chId, int &thisPMT)
{

    TString nameIs;
    thisPMT=-1;

    if ( strcmp(RpdPlanes,"RP01TA__") ==0){
       if ( chId >=0 && chId< 12) { // Au0 -Au11
	nameIs="TAu";   
	nameIs+=chId;
	thisPMT=chId;
      }
      else if(chId >= 32 && chId < 44)  { // Ad0 -Ad11
	nameIs="TAd";  
	nameIs+=chId-32;
	thisPMT=chId-32;
      }
      else {  
	nameIs="Unkn"; // unknown A channel
	thisPMT=-1;
      }
      
    }
    else if ( strcmp(RpdPlanes,"RP01TBl_")==0){
     
      if ( chId >=0 && chId< 12) { // Bu0 -Bu11
	nameIs="TBu";
	nameIs+=chId;
	thisPMT=chId;
      }
      else if(chId >= 16 && chId < 28)  { // Bu12 -Bu23
	nameIs="TBu";  
	nameIs+=12+chId-16;
	thisPMT=12+chId-16;
      }
      else if(chId >= 32 && chId < 44)  { // Bd0 -Bd11
	nameIs="TBd";  
	nameIs+=chId-32;
	thisPMT=chId-32;
      }
       else if(chId >= 48 && chId < 60)  { // Bd12 -Bd23
	 nameIs="TBd";  
	 nameIs+=12+chId-48;
	 thisPMT=12+chId-48;
       }
      else {  
	nameIs="Unkn"; // unknown B channel
	thisPMT=-1;
      }
      
    }
    else if ( strcmp(RpdPlanes,"RP01Ql__") ==0){
      

       if ( chId >=0 && chId< 12) { // Au0- Au11
	nameIs="QAu";
	nameIs+=chId;
	thisPMT=chId;
      }
       else if ( chId >=12 && chId< 24) { // Ad0- Ad11
	nameIs="QAd";
	nameIs+=chId-12;
	thisPMT=chId-12;
      }
        else if ( chId >=24 && chId< 48) { // Bu0- Bu23
	nameIs="QBu";
	nameIs+=chId-24;
	thisPMT=chId-24;
      }
      else if ( chId >=48 && chId< 72) { // Bd0- Bd23
	nameIs="QBd";
	nameIs+=chId-24-24;
	thisPMT=chId-24-24;
      } 
    }
   
    else {   
      nameIs="Unkn";  // function works only for RP01TA and RP01TBl and RP01Ql 
      thisPMT=-1;
    }
   
    return nameIs;
}

////////////////////////////////////////////////////////////////////////////////

void CsRPDetector::DecodeChipDigit(const CS::Chip::Digit &digit) {
}

////////////////////////////////////////////////////////////////////////////////

void CsRPDetector::clusterize() {
}

////////////////////////////////////////////////////////////////////////////////

void  CsRPDetector::readCalibration(time_t timePoint) {
}

/////////////////////////////////////////////////////////////////////////////////

void CsRPDetector::Clear(void) {
  CsDetector::Clear();
  ClearRPD();
}

/////////////////////////////////////////////////////////////////////////////////
// /*
// RPD Helper for Phast
// 
// comments & questions to 
// Etienne Burtin <Etienne.Burtin@cea.fr>
// Johannes Bernhard <Johannes.Bernhard@cern.ch>
// 
// Enjoy!
// 
// 09-02-12 18:50 jbernhar : Changed defintion of phi according to root definition (i.e. -pi < phi <= pi)
// 09-02-18 00:18 jbernhar : Implemented Etienne's corrections (energy loss in the target, calibration constants, vertex correction) 
// 09-02-21 18:30 jbernhar : Introduced multiple RPD Tracks
// 09-02-23 18:50 jbernhar : guarded several variables against division by 0
// 09-02-26 23:53 jbernhar : changed class structure to be according to Phast
// 
// */

void CsRPDetector::ClearRPD(void) {

fRPDTrack.clear(); // RPD Tracks are stored here
fzT.clear(); // primary vertex seen by the RPD for every track
dEA.clear(); // energy loss in A for every track
dEB.clear();// energy loss in B for every track
zA_vec.clear(); //z Position of Ring A hit
zB_vec.clear(); //z Position of Ring B hit
fCalibratedEvent.clear(); // Helper to see whether beta is between 0 and 1 (for every track)
hits.clear(); // pair of A and B hits for every track

nTracks = 0;
fMomentum = 0.;
fMom_x = 0.;
fMom_y = 0.;
fMom_z = 0.;
fTheta = 0.;
fPhi   = 0.;
fCosTheta = 0.;

memset(&fDataSentAupt, 0, sizeof(fDataSentAupt));
memset(&fDataSentAdot, 0, sizeof(fDataSentAdot));
memset(&fDataSentBupt, 0, sizeof(fDataSentBupt));
memset(&fDataSentBdot, 0, sizeof(fDataSentBdot));

memset(&fDataSentAupq, 0, sizeof(fDataSentAupq));
memset(&fDataSentAdoq, 0, sizeof(fDataSentAdoq));
memset(&fDataSentBupq, 0, sizeof(fDataSentBupq));
memset(&fDataSentBdoq, 0, sizeof(fDataSentBdoq));

memset(&fAdcAup, 0, sizeof(fAdcAup));
memset(&fAdcAdo, 0, sizeof(fAdcAdo));
memset(&fAdcBup, 0, sizeof(fAdcBup));
memset(&fAdcBdo, 0, sizeof(fAdcBdo));

memset(&fTdcAup, 0, sizeof(fTdcAup));
memset(&fTdcAdo, 0, sizeof(fTdcAdo));
memset(&fTdcBup, 0, sizeof(fTdcBup));
memset(&fTdcBdo, 0, sizeof(fTdcBdo));

zA    = -1000.;
tAoff = -1000.;
dtA   = -1000.;
tA    = -1000.;

zB    = -1000.;
tBoff = -1000.;
dtB   = -1000.;
tB    = -1000.;

beta = -1000.;
fdEA  = -1000.; 
fdEB  = -1000.;
dof  = -1000.;

something=false;

}

/////////////////////////////////////////////////////////////////////////////////

float CsRPDetector::treatADC(const CS::ChipSADC::Digit *ds)
{
  const vector<CS::uint16>& fsample = ds->GetSamples();
//calculate pedestal
float fPedestal=0.;
float maxAmplitude=0.;
float sumAmplitude=0.;
//int nSamples=32;
int nSamples=fsample.size();
int begSignal=0;
int samples=0;
for (int i=0;i<4;i++){
   fPedestal+=fsample[i];
}
fPedestal=fPedestal/4.;
//analyze waveform
for (int i=4;i<nSamples;i++){
   if (fsample[i]>maxAmplitude) maxAmplitude=fsample[i]; 
   if (fsample[i]>fsample[i-1]+3) begSignal=i;
   sumAmplitude+=fsample[i];
   samples++;
}
return sumAmplitude-samples*fPedestal;
  
}

/////////////////////////////////////////////////////////////////////////////////

void CsRPDetector::InitRPD( void )
{
initialized = true;

mP=0.938272;

tofOffset[0][23] = 4.44;     positionOffsetForA[0]=5.59;
tofOffset[0][0] = -0.64;     positionOffsetForB[0]=26.52;
tofOffset[0][1] = -0.78;     positionOffsetForB[1]=13.06;

tofOffset[1][1] = -1.97;     positionOffsetForA[1]=17.62;
tofOffset[1][2] = -2.55;     positionOffsetForB[2]=24.81;
tofOffset[1][3] = -2.37;     positionOffsetForB[3]=28.24;

tofOffset[2][3] = -0.10;     positionOffsetForA[2]=48.34;
tofOffset[2][4] = -0.86;     positionOffsetForB[4]=39.47;
tofOffset[2][5] = 1.35;      positionOffsetForB[5]=35.72;

tofOffset[3][5] = -3.74;     positionOffsetForA[3]=30.46;
tofOffset[3][6] = -4.97;     positionOffsetForB[6]=11.37;
tofOffset[3][7] = -5.70;     positionOffsetForB[7]=15.20;

tofOffset[4][7] = -4.92;     positionOffsetForA[4]=41.52;
tofOffset[4][8] = -4.09;     positionOffsetForB[8]=49.74;
tofOffset[4][9] = -0.68;     positionOffsetForB[9]=14.60;

tofOffset[5][9] = -0.44;     positionOffsetForA[5]=52.19;
tofOffset[5][10] = -2.27;    positionOffsetForB[10]=20.41;
tofOffset[5][11] = 0.40;     positionOffsetForB[11]=3.78;

tofOffset[6][11] = 0.72;     positionOffsetForA[6]=26.58;
tofOffset[6][12] = -4.99;    positionOffsetForB[12]=25.08;
tofOffset[6][13] = -0.82;    positionOffsetForB[13]=42.70;

tofOffset[7][13] = -2.75;    positionOffsetForA[7]=24.33;
tofOffset[7][14] = -4.92;    positionOffsetForB[14]=-0.18;
tofOffset[7][15] = -3.33;    positionOffsetForB[15]=17.47;

tofOffset[8][15] = 2.39;     positionOffsetForA[8]=21.59;
tofOffset[8][16] = -0.23;    positionOffsetForB[16]=-10.47;
tofOffset[8][17] = 5.51;     positionOffsetForB[17]=27.07;

tofOffset[9][17] = 0.92;     positionOffsetForA[9]=4.85;
tofOffset[9][18] = -3.50;    positionOffsetForB[18]=17.68;
tofOffset[9][19] = -2.11;    positionOffsetForB[19]=18.48;

tofOffset[10][19] = -2.16;   positionOffsetForA[10]=13.08;
tofOffset[10][20] = -2.62;   positionOffsetForB[20]=1.65;
tofOffset[10][21] = 0.32;    positionOffsetForB[21]=9.61;

//tofOffset[11][21] = -1.3;  positionOffsetForA[11]=7.31;
tofOffset[11][21] = -1.01;   positionOffsetForA[11]=7.31;
tofOffset[11][22] = -5.85; positionOffsetForB[22]=-0.17;
tofOffset[11][23] = -0.78; positionOffsetForB[23]=-17.07;

lightSpeedInA[0]=12.;
lightSpeedInA[1]=12.;
lightSpeedInA[2]=11.5;
lightSpeedInA[3]=11.5;
lightSpeedInA[4]=12.;
lightSpeedInA[5]=12.;
lightSpeedInA[6]=12.;
lightSpeedInA[7]=11.;
lightSpeedInA[8]=11.;
lightSpeedInA[9]=12.;
lightSpeedInA[10]=12.;
lightSpeedInA[11]=12.;

memset(&lightSpeedInB, 13, sizeof(lightSpeedInB));

// dEBcusp[ 0]= 2280.; //old calibration
// dEBcusp[ 1]=  2150.;
// dEBcusp[ 2]= 1830.;
// dEBcusp[ 3]=  2700.;
// dEBcusp[ 4]= 2180.;
// dEBcusp[ 5]=  2100.;
// dEBcusp[ 6]= 2340.;
// dEBcusp[ 7]=  2000.;
// dEBcusp[ 8]= 2270.;
// dEBcusp[ 9]=  2300.;
// dEBcusp[10]= 2280.;
// dEBcusp[11]=  2050.;
// dEBcusp[12]= 2165.;
// dEBcusp[13]=  2200.;
// dEBcusp[14]= 2250.;
// dEBcusp[15]=  2100.;
// dEBcusp[16]= 2430.;
// dEBcusp[17]=  2000.;
// dEBcusp[18]= 2440.;
// dEBcusp[19]=  2000.;
// dEBcusp[20]= 2360.;
// dEBcusp[21]=  2400.;
// dEBcusp[22]= 2340.;
// dEBcusp[23]=  2000.;

dEBcusp[ 0]= 2403.;
dEBcusp[ 1]=2117.;
dEBcusp[ 2]= 1886.;
dEBcusp[ 3]=2690.2;
dEBcusp[ 4]= 2244.;
dEBcusp[ 5]=2072.;
dEBcusp[ 6]= 2409.;
dEBcusp[ 7]=2022.;
dEBcusp[ 8]= 2318.;
dEBcusp[ 9]=2330.;
dEBcusp[10]= 2317.;
dEBcusp[11]=2047.;
dEBcusp[12]= 2238.;
dEBcusp[13]=2247.;
dEBcusp[14]= 2332.;
dEBcusp[15]=2242.;
dEBcusp[16]= 2490.;
dEBcusp[17]=2052.;
dEBcusp[18]= 2549.;
dEBcusp[19]=2050.;
dEBcusp[20]= 2495.;
dEBcusp[21]=2368.;
dEBcusp[22]= 2371.;
dEBcusp[23]=2048.;

dEAcusp[ 0]= 1470.;
dEAcusp[ 1]= 1280.;
dEAcusp[ 2]= 1100.;
dEAcusp[ 3]= 1280.;
dEAcusp[ 4]= 1190.;
dEAcusp[ 5]= 1260.;
dEAcusp[ 6]= 1260.;
dEAcusp[ 7]= 1480.;
dEAcusp[ 8]=  800.;
dEAcusp[ 9]= 1370.;
dEAcusp[10]= 1230.;
dEAcusp[11]= 1240.;

zVertexOffset[0][23] = -1.80; 
zVertexOffset[0][0] = -2.40; 
zVertexOffset[0][1] = -3.12; 
zVertexOffset[1][1] = -3.31; 
zVertexOffset[1][2] = -3.59; 
zVertexOffset[1][3] = -3.84; 
zVertexOffset[2][3] = -4.26; 
zVertexOffset[2][4] = -4.43; 
zVertexOffset[2][5] = -4.71; 
zVertexOffset[3][5] = -3.83; 
zVertexOffset[3][6] = -3.75; 
zVertexOffset[3][7] = -3.36; 
zVertexOffset[4][7] = -1.88; 
zVertexOffset[4][8] = -1.56; 
zVertexOffset[4][9] = -1.04; 
zVertexOffset[5][9] = 0.53; 
zVertexOffset[5][10] = 1.07; 
zVertexOffset[5][11] = 2.13; 
zVertexOffset[6][11] = 2.49; 
zVertexOffset[6][12] = 2.94; 
zVertexOffset[6][13] = 3.44; 
zVertexOffset[7][13] = 4.78; 
zVertexOffset[7][14] = 4.69; 
zVertexOffset[7][15] = 4.67; 
zVertexOffset[8][15] = 4.62; 
zVertexOffset[8][16] = 4.73; 
zVertexOffset[8][17] = 4.47; 
zVertexOffset[9][17] = 5.00; 
zVertexOffset[9][18] = 4.58; 
zVertexOffset[9][19] = 4.35; 
zVertexOffset[10][19] = 3.63; 
zVertexOffset[10][20] = 3.07; 
zVertexOffset[10][21] = 2.56; 
zVertexOffset[11][21] = 0.77; 
zVertexOffset[11][22] = 0.10; 
zVertexOffset[11][23] = -0.62;

rA=12.5;
rB=75.;
deltaR=rB-rA;
lightSpeed=30.; // cm/ns
dEAcuspMeV=8.5;
dEBcuspMeV=24.;
tofGlobalOffset=2.2;
zTargetCenter=-48.;
}

/////////////////////////////////////////////////////////////////////////////////

float CsRPDetector::correctEnergyLoss(float pRpd, float sinTheta, float phi, float xVertex, float yVertex, int debug){
  const float rTarget=17.5;
  const float chaussetteThickness=0.2;
  const float aluminiumThickness=1.8;
  const float ringAThickness=5.;
  const float cellThickness=0.125;
  float lengthInH2;
  pRpd=correction_energy(pRpd,4,ringAThickness/sinTheta,debug); // scint A
  pRpd=correction_energy(pRpd,3,chaussetteThickness/sinTheta,debug); // mylar
  pRpd=correction_energy(pRpd,2,aluminiumThickness/sinTheta,debug); // aluminium
  pRpd=correction_energy(pRpd,3,cellThickness/sinTheta,debug); // mylar
  lengthInH2=TMath::Sqrt(pow(rTarget*TMath::Cos(phi)-xVertex,2)+pow(rTarget*TMath::Sin(phi)-yVertex,2));
  lengthInH2=lengthInH2/sinTheta;
  pRpd=correction_energy(pRpd,1,lengthInH2,debug);// hydrogen
  return pRpd;
}

/////////////////////////////////////////////////////////////////////////////////

float CsRPDetector::correction_energy(float p, int mat, float thickness, int debug){
        const int nData=151;                                                          

// --- kinetic energy in MeV  						 

        float Ek[nData]={                                                        
     	  4.0,  8.0, 12.0, 16.0, 20.0, 24.0, 28.0, 32.0, 36.0, 40.0,	       
     	 44.0, 48.0, 52.0, 56.0, 60.0, 64.0, 68.0, 72.0, 76.0, 80.0,	       
     	 84.0, 88.0, 92.0, 96.0,100.0,104.0,108.0,112.0,116.0,120.0,	       
     	124.0,128.0,132.0,136.0,140.0,144.0,148.0,152.0,156.0,160.0,	       
     	164.0,168.0,172.0,176.0,180.0,184.0,188.0,192.0,196.0,200.0,	       
     	204.0,208.0,212.0,216.0,220.0,224.0,228.0,232.0,236.0,240.0,	       
     	244.0,248.0,252.0,256.0,260.0,264.0,268.0,272.0,276.0,280.0,	       
     	284.0,288.0,292.0,296.0,300.0,304.0,308.0,312.0,316.0,320.0,	       
     	324.0,328.0,332.0,336.0,340.0,344.0,348.0,352.0,356.0,360.0,	       
     	364.0,368.0,372.0,376.0,380.0,384.0,388.0,392.0,396.0,400.0,	       
     	404.0,408.0,412.0,416.0,420.0,424.0,428.0,432.0,436.0,440.0,	       
     	444.0,448.0,452.0,456.0,460.0,464.0,468.0,472.0,476.0,480.0,	       
     	484.0,488.0,492.0,496.0,500.0,504.0,508.0,512.0,516.0,520.0,	       
     	524.0,528.0,532.0,536.0,540.0,544.0,548.0,552.0,556.0,560.0,	       
     	564.0,568.0,572.0,576.0,580.0,584.0,588.0,592.0,596.0,600.0};	       

// --- protons range in H2 mm  						 										 

      float range_H2[nData]={                                                       
     	0.1373E+01, 0.4879E+01, 0.1030E+02, 0.1751E+02, 0.2638E+02,
     	0.3700E+02, 0.4917E+02, 0.6288E+02, 0.7799E+02, 0.9461E+02,
     	0.1130E+03, 0.1323E+03, 0.1534E+03, 0.1756E+03, 0.1994E+03,
     	0.2239E+03, 0.2507E+03, 0.2774E+03, 0.3067E+03, 0.3368E+03,
     	0.3668E+03, 0.4001E+03, 0.4338E+03, 0.4674E+03, 0.5035E+03,
     	0.5411E+03, 0.5787E+03, 0.6164E+03, 0.6574E+03, 0.6993E+03,
     	0.7413E+03, 0.7832E+03, 0.8269E+03, 0.8735E+03, 0.9202E+03,
     	0.9669E+03, 0.1014E+04, 0.1061E+04, 0.1113E+04, 0.1165E+04,
     	0.1217E+04, 0.1268E+04, 0.1320E+04, 0.1373E+04, 0.1431E+04,
     	0.1488E+04, 0.1545E+04, 0.1602E+04, 0.1659E+04, 0.1716E+04,
     	0.1776E+04, 0.1839E+04, 0.1902E+04, 0.1965E+04, 0.2028E+04,
     	0.2091E+04, 0.2154E+04, 0.2216E+04, 0.2282E+04, 0.2351E+04,
     	0.2420E+04, 0.2489E+04, 0.2558E+04, 0.2627E+04, 0.2696E+04,
     	0.2765E+04, 0.2834E+04, 0.2905E+04, 0.2980E+04, 0.3055E+04,
     	0.3130E+04, 0.3205E+04, 0.3281E+04, 0.3356E+04, 0.3431E+04,
     	0.3506E+04, 0.3581E+04, 0.3657E+04, 0.3736E+04, 0.3818E+04,
     	0.3899E+04, 0.3981E+04, 0.4062E+04, 0.4144E+04, 0.4226E+04,
     	0.4307E+04, 0.4389E+04, 0.4470E+04, 0.4552E+04, 0.4633E+04,
     	0.4718E+04, 0.4806E+04, 0.4894E+04, 0.4982E+04, 0.5070E+04,
     	0.5158E+04, 0.5246E+04, 0.5334E+04, 0.5422E+04, 0.5510E+04,
     	0.5598E+04, 0.5686E+04, 0.5774E+04, 0.5862E+04, 0.5951E+04,
     	0.6045E+04, 0.6140E+04, 0.6234E+04, 0.6329E+04, 0.6423E+04,
     	0.6517E+04, 0.6612E+04, 0.6706E+04, 0.6801E+04, 0.6895E+04,
     	0.6989E+04, 0.7084E+04, 0.7178E+04, 0.7273E+04, 0.7367E+04,
     	0.7462E+04, 0.7560E+04, 0.7661E+04, 0.7762E+04, 0.7862E+04,
     	0.7963E+04, 0.8064E+04, 0.8164E+04, 0.8265E+04, 0.8366E+04,
     	0.8466E+04, 0.8567E+04, 0.8668E+04, 0.8768E+04, 0.8869E+04,
     	0.8970E+04, 0.9070E+04, 0.9171E+04, 0.9272E+04, 0.9372E+04,
     	0.9477E+04, 0.9584E+04, 0.9690E+04, 0.9797E+04, 0.9904E+04,
     	0.1001E+05, 0.1012E+05, 0.1022E+05, 0.1033E+05, 0.1044E+05};
                                                                               

      float range_Aluminium[nData]={                                                      
     	 0.1290E+00, 0.4173E+00, 0.8474E+00, 0.1408E+01, 0.2087E+01,
     	 0.2893E+01, 0.3809E+01, 0.4833E+01, 0.5958E+01, 0.7188E+01,
     	 0.8542E+01, 0.9960E+01, 0.1151E+02, 0.1312E+02, 0.1486E+02,
     	 0.1664E+02, 0.1857E+02, 0.2051E+02, 0.2262E+02, 0.2479E+02,
     	 0.2695E+02, 0.2934E+02, 0.3175E+02, 0.3416E+02, 0.3674E+02,
     	 0.3943E+02, 0.4211E+02, 0.4480E+02, 0.4771E+02, 0.5069E+02,
     	 0.5367E+02, 0.5665E+02, 0.5975E+02, 0.6305E+02, 0.6635E+02,
     	 0.6965E+02, 0.7295E+02, 0.7633E+02, 0.7998E+02, 0.8362E+02,
     	 0.8727E+02, 0.9092E+02, 0.9456E+02, 0.9832E+02, 0.1023E+03,
     	 0.1063E+03, 0.1104E+03, 0.1144E+03, 0.1184E+03, 0.1224E+03,
     	 0.1266E+03, 0.1310E+03, 0.1354E+03, 0.1398E+03, 0.1442E+03,
     	 0.1486E+03, 0.1530E+03, 0.1574E+03, 0.1620E+03, 0.1668E+03,
     	 0.1716E+03, 0.1764E+03, 0.1812E+03, 0.1860E+03, 0.1908E+03,
     	 0.1956E+03, 0.2004E+03, 0.2053E+03, 0.2105E+03, 0.2158E+03,
     	 0.2210E+03, 0.2262E+03, 0.2314E+03, 0.2366E+03, 0.2419E+03,
     	 0.2471E+03, 0.2523E+03, 0.2575E+03, 0.2630E+03, 0.2687E+03,
     	 0.2743E+03, 0.2800E+03, 0.2856E+03, 0.2912E+03, 0.2969E+03,
     	 0.3025E+03, 0.3082E+03, 0.3138E+03, 0.3195E+03, 0.3251E+03,
     	 0.3309E+03, 0.3370E+03, 0.3431E+03, 0.3491E+03, 0.3552E+03,
     	 0.3613E+03, 0.3673E+03, 0.3734E+03, 0.3795E+03, 0.3855E+03,
     	 0.3916E+03, 0.3977E+03, 0.4037E+03, 0.4098E+03, 0.4159E+03,
     	 0.4224E+03, 0.4289E+03, 0.4354E+03, 0.4419E+03, 0.4484E+03,
     	 0.4549E+03, 0.4614E+03, 0.4678E+03, 0.4743E+03, 0.4808E+03,
     	 0.4873E+03, 0.4938E+03, 0.5003E+03, 0.5068E+03, 0.5133E+03,
     	 0.5197E+03, 0.5265E+03, 0.5334E+03, 0.5403E+03, 0.5472E+03,
     	 0.5541E+03, 0.5610E+03, 0.5679E+03, 0.5748E+03, 0.5817E+03,
     	 0.5886E+03, 0.5955E+03, 0.6024E+03, 0.6093E+03, 0.6161E+03,
     	 0.6230E+03, 0.6299E+03, 0.6368E+03, 0.6437E+03, 0.6506E+03,
     	 0.6578E+03, 0.6651E+03, 0.6723E+03, 0.6796E+03, 0.6869E+03,
     	 0.6942E+03, 0.7014E+03, 0.7087E+03, 0.7160E+03, 0.7233E+03};
                                                                               

// --- protons range in Mylar mm  						 										 

      float range_Mylar[nData]={                                                        
     	 0.1851E+00, 0.6286E+00, 0.1300E+01, 0.2182E+01, 0.3259E+01,
     	 0.4540E+01, 0.6001E+01, 0.7640E+01, 0.9442E+01, 0.1142E+02,
     	 0.1360E+02, 0.1588E+02, 0.1837E+02, 0.2099E+02, 0.2379E+02,
     	 0.2667E+02, 0.2981E+02, 0.3294E+02, 0.3637E+02, 0.3988E+02,
     	 0.4339E+02, 0.4728E+02, 0.5120E+02, 0.5512E+02, 0.5932E+02,
     	 0.6369E+02, 0.6806E+02, 0.7244E+02, 0.7719E+02, 0.8206E+02,
     	 0.8692E+02, 0.9178E+02, 0.9683E+02, 0.1022E+03, 0.1076E+03,
     	 0.1130E+03, 0.1184E+03, 0.1239E+03, 0.1299E+03, 0.1359E+03,
     	 0.1418E+03, 0.1478E+03, 0.1538E+03, 0.1599E+03, 0.1665E+03,
     	 0.1731E+03, 0.1797E+03, 0.1862E+03, 0.1928E+03, 0.1994E+03,
     	 0.2062E+03, 0.2135E+03, 0.2207E+03, 0.2279E+03, 0.2351E+03,
     	 0.2424E+03, 0.2496E+03, 0.2568E+03, 0.2644E+03, 0.2723E+03,
     	 0.2801E+03, 0.2880E+03, 0.2959E+03, 0.3038E+03, 0.3117E+03,
     	 0.3196E+03, 0.3275E+03, 0.3356E+03, 0.3442E+03, 0.3528E+03,
     	 0.3614E+03, 0.3700E+03, 0.3786E+03, 0.3872E+03, 0.3958E+03,
     	 0.4044E+03, 0.4130E+03, 0.4215E+03, 0.4306E+03, 0.4399E+03,
     	 0.4492E+03, 0.4585E+03, 0.4678E+03, 0.4771E+03, 0.4864E+03,
     	 0.4957E+03, 0.5050E+03, 0.5143E+03, 0.5236E+03, 0.5329E+03,
     	 0.5425E+03, 0.5525E+03, 0.5626E+03, 0.5726E+03, 0.5826E+03,
     	 0.5926E+03, 0.6026E+03, 0.6126E+03, 0.6226E+03, 0.6327E+03,
     	 0.6427E+03, 0.6527E+03, 0.6627E+03, 0.6727E+03, 0.6829E+03,
     	 0.6936E+03, 0.7043E+03, 0.7150E+03, 0.7257E+03, 0.7365E+03,
     	 0.7472E+03, 0.7579E+03, 0.7686E+03, 0.7794E+03, 0.7901E+03,
     	 0.8008E+03, 0.8115E+03, 0.8222E+03, 0.8330E+03, 0.8437E+03,
     	 0.8544E+03, 0.8656E+03, 0.8770E+03, 0.8884E+03, 0.8998E+03,
     	 0.9112E+03, 0.9226E+03, 0.9341E+03, 0.9455E+03, 0.9569E+03,
     	 0.9683E+03, 0.9797E+03, 0.9911E+03, 0.1002E+04, 0.1014E+04,
     	 0.1025E+04, 0.1037E+04, 0.1048E+04, 0.1060E+04, 0.1071E+04,
     	 0.1083E+04, 0.1095E+04, 0.1107E+04, 0.1119E+04, 0.1131E+04,
     	 0.1143E+04, 0.1155E+04, 0.1167E+04, 0.1179E+04, 0.1191E+04};
                                                                               

// --- protons range in Scintillator mm  						 										 

      float range_Scintillator[nData]={                                                       
     	 0.2301E+00, 0.7883E+00, 0.1636E+01, 0.2752E+01, 0.4114E+01,
     	 0.5737E+01, 0.7590E+01, 0.9670E+01, 0.1196E+02, 0.1446E+02,
     	 0.1723E+02, 0.2014E+02, 0.2331E+02, 0.2663E+02, 0.3019E+02,
     	 0.3385E+02, 0.3784E+02, 0.4183E+02, 0.4620E+02, 0.5066E+02,
     	 0.5513E+02, 0.6008E+02, 0.6508E+02, 0.7007E+02, 0.7541E+02,
     	 0.8098E+02, 0.8655E+02, 0.9212E+02, 0.9818E+02, 0.1044E+03,
     	 0.1106E+03, 0.1168E+03, 0.1232E+03, 0.1301E+03, 0.1370E+03,
     	 0.1438E+03, 0.1507E+03, 0.1578E+03, 0.1654E+03, 0.1730E+03,
     	 0.1806E+03, 0.1882E+03, 0.1958E+03, 0.2036E+03, 0.2120E+03,
     	 0.2204E+03, 0.2288E+03, 0.2372E+03, 0.2456E+03, 0.2540E+03,
     	 0.2627E+03, 0.2719E+03, 0.2812E+03, 0.2904E+03, 0.2996E+03,
     	 0.3088E+03, 0.3180E+03, 0.3272E+03, 0.3369E+03, 0.3470E+03,
     	 0.3570E+03, 0.3671E+03, 0.3772E+03, 0.3873E+03, 0.3973E+03,
     	 0.4074E+03, 0.4175E+03, 0.4278E+03, 0.4388E+03, 0.4498E+03,
     	 0.4608E+03, 0.4717E+03, 0.4827E+03, 0.4937E+03, 0.5046E+03,
     	 0.5156E+03, 0.5266E+03, 0.5375E+03, 0.5491E+03, 0.5610E+03,
     	 0.5729E+03, 0.5848E+03, 0.5966E+03, 0.6085E+03, 0.6204E+03,
     	 0.6323E+03, 0.6442E+03, 0.6560E+03, 0.6679E+03, 0.6798E+03,
     	 0.6920E+03, 0.7048E+03, 0.7176E+03, 0.7304E+03, 0.7432E+03,
     	 0.7560E+03, 0.7688E+03, 0.7816E+03, 0.7944E+03, 0.8072E+03,
     	 0.8200E+03, 0.8328E+03, 0.8456E+03, 0.8584E+03, 0.8713E+03,
     	 0.8850E+03, 0.8987E+03, 0.9124E+03, 0.9261E+03, 0.9398E+03,
     	 0.9535E+03, 0.9672E+03, 0.9809E+03, 0.9946E+03, 0.1008E+04,
     	 0.1022E+04, 0.1036E+04, 0.1049E+04, 0.1063E+04, 0.1077E+04,
     	 0.1090E+04, 0.1105E+04, 0.1119E+04, 0.1134E+04, 0.1149E+04,
     	 0.1163E+04, 0.1178E+04, 0.1192E+04, 0.1207E+04, 0.1221E+04,
     	 0.1236E+04, 0.1251E+04, 0.1265E+04, 0.1280E+04, 0.1294E+04,
     	 0.1309E+04, 0.1324E+04, 0.1338E+04, 0.1353E+04, 0.1367E+04,
     	 0.1382E+04, 0.1398E+04, 0.1413E+04, 0.1429E+04, 0.1444E+04,
     	 0.1459E+04, 0.1475E+04, 0.1490E+04, 0.1506E+04, 0.1521E+04};

p=p*1000.;
float E=TMath::Sqrt(p*p+mP*mP)-mP;
float beta=p/TMath::Sqrt(p*p+mP*mP);

bool found=false;
float range=0.;
float energy=0.;
float a=0.,b=0.,c=0.;
//int mat=3;

for (int i=0;i<nData;i++){
   if(E<Ek[i] && !found) {
      found=true;
      if (mat==1){
        if (debug) printf("***********************\n in hydrogen\n");
        a=range_H2[i-1];
        b=range_H2[i];
        c=range_H2[i+1];
      }else if (mat==2){
        if (debug) printf("***********************\n in Aluminium\n");
        a=range_Aluminium[i-1];
        b=range_Aluminium[i];
        c=range_Aluminium[i+1];
      }else if (mat==3){
        if (debug) printf("***********************\n in Mylar\n");
        a=range_Mylar[i-1];
        b=range_Mylar[i];
        c=range_Mylar[i+1];
      }else if (mat==4){
        if (debug) printf("***********************\n in Scintillator\n");
        a=range_Scintillator[i-1];
        b=range_Scintillator[i];
        c=range_Scintillator[i+1];
      }
      if (debug) printf("p=%f beta=%f E=%3.2f Ekin=%3.2f,%3.2f,%3.2f range=%4.2f,%4.2f,%4.2f\n",p,beta,E,Ek[i-1],Ek[i],Ek[i+1],a,b,c);
      range=terpol3(Ek[i-1],Ek[i],Ek[i+1],a,b,c,E);
      if (debug) printf("interpolated range = %4.2f \n",range);
      range+=thickness;
      energy=terpol3(a,b,c,Ek[i-1],Ek[i],Ek[i+1],range);
      if (debug) printf("energy:%f p=%f\n",energy, TMath::Sqrt(energy*energy+2*energy*mP));
   }
}
return (isnan(energy) || isinf(energy)) ? 0 : TMath::Sqrt(energy*energy+2*energy*mP)/1000.;
}

/////////////////////////////////////////////////////////////////////////////////

float CsRPDetector::terpol3(float xx1,float xx2,float xx3,float yy1,float yy2,float yy3,float x){
  float val  = yy1*(x-xx2)*(x-xx3)/((xx1-xx2)*(xx1-xx3))
             + yy2*(x-xx1)*(x-xx3)/((xx2-xx1)*(xx2-xx3))
             + yy3*(x-xx1)*(x-xx2)/((xx3-xx1)*(xx3-xx2));
  return val;
}

/////////////////////////////////////////////////////////////////////////////////

//void RPD::Search(const PaEvent& e, bool correct_momentum = false, double PV_x = 0., double PV_y = 0., bool correct_overlaps = false) {
//void RPD::Search(const PaEvent& e) {
void CsRPDetector::Search( void) {

//  bool debug = true;
//  bool correct_momentum = false;
//  double PV_x = 0.;
//  double PV_y = 0.;
//  bool correct_overlaps = false;
// // Searches the data for RPD infos and calculates the tracks
// 
// Clear();
//  if( debug )
//  {
//    cout <<" Sorry RPD::Search is not debuged in coral exit " << endl;
//    exit(0);
//  }
// //DecodePhast(e);
// DecodePhast();
// 
// int iB=0;
// for (int iA=0;iA<12;iA++){
// for (int ilB=-1;ilB<2;ilB++){
// 
// if (iA==0 && ilB==-1) {
// iB=23;
// }else{ 
// iB=2*iA+ilB;
// }
// 
// for (int iHitAup=0;iHitAup<fDataSentAupt[iA];iHitAup++){
// for (int iHitAdo=0;iHitAdo<fDataSentAdot[iA];iHitAdo++){
// for (int iHitBup=0;iHitBup<fDataSentBupt[iB];iHitBup++){
// for (int iHitBdo=0;iHitBdo<fDataSentBdot[iB];iHitBdo++){
// 
// if (fAdcAdo[iA]<10 || fAdcAup[iA]<10) continue;
// if (fAdcBdo[iB]<10 || fAdcBup[iB]<10) continue;
// 
// dtA = (float)(fTdcAup[iA][iHitAup]-fTdcAdo[iA][iHitAdo]);
// zA = dtA*lightSpeedInA[iA]/2.+positionOffsetForA[iA];
// tA = (fTdcAup[iA][iHitAup]+fTdcAdo[iA][iHitAdo])/2.;
// 
// dtB = (float)(fTdcBup[iB][iHitBup]-fTdcBdo[iB][iHitBdo]);
// zB = dtB*lightSpeedInB[iB]/2. + positionOffsetForB[iB];
// tB = (fTdcBup[iB][iHitBup]+fTdcBdo[iB][iHitBdo])/2.;
// 
// zA = zTargetCenter+zA-13.; 
// zB = zTargetCenter+zB-8.;
// 
// // if (zB - zA < -50. || zA > 75. || zB > 200.) continue;
// 
// dof = TMath::Sqrt((zA-zB)*(zA-zB)+deltaR*deltaR);
// beta = (tB-tA+tofOffset[iA][iB]+tofGlobalOffset > 1e-10) ? dof/(tB-tA+tofOffset[iA][iB]+tofGlobalOffset)/lightSpeed : 0.995; // was 0
// beta = (beta > 0.9995) ? 0 : beta; // introduced cutoff for not physical momenta 
// beta = (beta < 0.) ? 0. : beta;
// fdEA = TMath::Sqrt(fAdcAdo[iA]*fAdcAup[iA]); 
// fdEB = TMath::Sqrt(fAdcBdo[iB]*fAdcBup[iB]); 
// fdEA = fdEA*dEAcuspMeV/dEAcusp[iA];
// fdEB = fdEB*dEBcuspMeV/dEBcusp[iB];
// 
// something=true;
// 
// fzT.push_back(zA-rA*(zB-zA)/deltaR+zVertexOffset[iA][iB]);
// fCosTheta = (zB-zA)/dof;
// fSinTheta = deltaR/dof;
// fTheta = TMath::ACos((zB-zA)/dof);
// fPhi   = 2. * TMath::Pi() * 15. / 360. * (iB + 1e-10);
// if (fPhi > TMath::Pi()) fPhi -= 2. * TMath::Pi();
// if (ilB) fPhi -= 2. * TMath::Pi() * 3.75 / 360. * ilB;
// fMomentum = beta*mP/TMath::Sqrt(1.-beta*beta);
// if (correct_momentum && fMomentum > 1e-6 && fMomentum < 0.596) fMomentum = correctEnergyLoss(fMomentum,fSinTheta,fPhi,PV_x,PV_y,0);// correct for energy loss in target naterial
// fMom_x = fMomentum * fSinTheta * TMath::Cos(fPhi);
// fMom_y = fMomentum * fSinTheta * TMath::Sin(fPhi);
// fMom_z = fMomentum * fCosTheta;
// hits.push_back(make_pair(iA,iB));
// fRPD4Vector.SetXYZM(fMom_x,fMom_y,fMom_z,mP); // assigning proton mass to the 4-vector
// nTracks++;
// fRPDTrack.push_back(fRPD4Vector);
// dEA.push_back(fdEA);
// dEB.push_back(fdEB);
// zA_vec.push_back(zA);
// zB_vec.push_back(zB);
// fCalibratedEvent.push_back(true);
// }
// }
// }
// }
// }
// }
// 
// if (!something) {
// fRPD4Vector.SetXYZT(0.,0.,0.,0.);
// fzT.push_back(-777);
// hits.push_back(make_pair(-1,-1));
// fCalibratedEvent.push_back(false);
// fRPDTrack.push_back(fRPD4Vector);
// dEA.push_back(0);
// dEB.push_back(0);
// zA_vec.push_back(-777);
// zB_vec.push_back(-777);
// }
// 
// // start correction for possible Ring B overlaps (also possible target excitation killer, delta would be seen as proton, careful!)
// // to be added
// 
// // guess best proton track
// fBestProtonTrack = 0;
// for (unsigned int i=0; i<fRPDTrack.size(); i++) if ( (fRPDTrack[i].Vect().Mag() > 0.) && (fabs(fzT[i] + 48.) < 20.) && fCalibratedEvent[i] ) fBestProtonTrack = i;
// 
}

////////////////////////////////////////////////////////////////////////////////


