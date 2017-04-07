/*!
   \file    
   \brief   
   \author  Kolosov Vladimir            Vladimir.Kolosov@cern.ch   
   \version 
   \date    
*/

//////////////////////////////// CsRICH1UpGrade FUNCTIONS //////////////////////////////////////////////////

#include "DaqDataDecoding/ChipGassiplex.h"
#include "DaqDataDecoding/ChipAPVRICH.h"
#include "DaqDataDecoding/ChipF1.h"       //^

#include "CsRICH1UpGrade.h"
#include "CsRICH1Detector.h"
#include "CsMCRICH1Hit.h"
#include "CsOpt.h"
#include "CsDigit.h"
#include "CsMCDigit.h"
#include "CsRCDetectors.h" // Problematic header !!!
#include "CathodeMAPMT.h"
#include "CathodeAPV.h"
#include <cassert>

// #include "RICH1Analyser.h"

/* Andrea ... include of way CsRICH1Detector.h to have access to getPMT_T0 */

#include "CsRICH1Detector.h"

/* Andrea ... end of include of way CsRICH1Detector.h to have access to getPMT_T0 */

using namespace std;
using namespace CLHEP;


//////////////////////////////////////////////////////////////////////////////////

CsRICH1UpGrade::CsRICH1UpGrade ( const int id, const std::string &TBname ) : CsRICH1Detector(id,TBname) 
{
  bool debug = false;
  if( debug ) cout << " Construct CsRICH1UpGrade " << endl;
  CsOpt* opt = CsOpt::Instance();
  bool boo;
  std::string tag = "CsRICH1UpGrade";
  std::string key = "CONFIG";
  std::string value = "";
  if( opt->getOpt( tag, key, value ) )
// Warning : if 'value' is left blank,  opt->getOpt(... fails
  {
    cout << " CsRICH1UpGrade CONFIG " << value << endl;
    if( value == "Gassiplex" ) 
      cout << " Configure RICH1 with Gassiplex read-out " << endl;
    else if( value == "APV" ) 
      cout << " Configure RICH1 with APV read-out " << endl;
    else if( value == "TEST" ) 
      cout << " Configure RICH1 with Gassiplex and one cathode with APV read-out " << endl;
    else  
      cout << " Configure RICH1 with MAPMT and APV read-out " << endl;
  }

  double xsize = 600.;
  double ysize = 600.;
  char name[132];
  int ic = 0;
  for( unsigned updown=0; updown < 2; updown++ )
  {
    for( unsigned col=0; col < 4; col++ )
    {
      for( unsigned row=0; row < 2; row++ )
      {
        if( debug ) cout << " Construct CsRICH1UpGrade make new cathode ic " << ic << endl;
        CathodePlane *ch = NULL;
        if( value == "IDEAL" )
        {
          ch = new CathodePlane( ic);
          if( ic < 10 )  
            sprintf(name,"RI01P0%d",ic);
          else
            sprintf(name,"RI01P%d",ic);
          ch->SetTBname( name );
        }
        if( value == "GASSIPLEX" )
        {
          ch = new CathodeAPV( ic );
          ch->SetReadOutType(CsRICH1UpGrade::CathodePlane::GASSIPLEX);
          if( ic < 10 )  
            sprintf(name,"RI01P0%d",ic);
          else
            sprintf(name,"RI01P%d",ic);
          ch->SetTBname( name );
        }
        else if( value == "APV" )
        {
          ch = new CathodeAPV( ic );
          if( ic < 10 )  
            sprintf(name,"RA01P0%d",ic);
          else
            sprintf(name,"RA01P%d",ic);
          ch->SetTBname( name );
        } 
        else if( value == "TEST" )
        {
          ch = new CathodeAPV( ic );
          if( ic != 12)
          { 
            ch->SetReadOutType(CsRICH1UpGrade::CathodePlane::GASSIPLEX);
            if( ic < 10 )  
              sprintf(name,"RI01P0%d",ic);
            else
              sprintf(name,"RI01P%d",ic);
            ch->SetTBname( name );
          }
          else
          {
            if( ic < 10 )  
              sprintf(name,"RA01P0%d",ic);
            else
              sprintf(name,"RA01P%d",ic);
            ch->SetTBname( name );
          }
        }
        else
        {
          if( ic == 3 || ic == 5 || ic == 10 || ic == 12 )
          {
            ch = new CathodeMAPMT( ic );
            if( ic < 10 )  
              sprintf(name,"RM01P0%d",ic);
	    //^sprintf(name,"RP01P0%d",ic);
            else
              sprintf(name,"RM01P%d",ic);
	    //^sprintf(name,"RP01P%d",ic);
            ch->SetTBname( name );
          }
          else
          {  
            ch = new CathodeAPV( ic );
            if( ic < 10 )  
              sprintf(name,"RA01P0%d",ic);
            else
              sprintf(name,"RA01P%d",ic);
            ch->SetTBname( name );
          }
        }
        if( ch == NULL ) assert( false );

        ch->x0_image_ = xsize + xsize/2. - xsize*( (double)col );
        ch->y0_image_ = ysize + ysize/2. - ysize*2*( (double)updown ) - ysize*( (double)row );
        cathodes_.push_back( ch );
        ic++;
      }
    }
  }
  
  Clear();
  InitCsOpt();  
  
  if( debug ) cout << " Construct CsRICH1UpGrade OK " << endl;
  if( options_.printinfo_ )
    std::cout << " Construct CsRICH1UpGrade detector " << std::endl;
}

//////////////////////////////////////////////////////////////////////////////////

void CsRICH1UpGrade::DecodeChipDigits(const CS::Chip::Digits &digits)
{
  bool debug = false;
//   bool debug = true;
  typedef std::multimap<CS::DetID,CS::Chip::Digit*>::const_iterator m_it; // iterator type
// Decode Gassiplex digits
  const static std::string richname("RI01P");
  std::pair<m_it,m_it> m_range = digits.equal_range( richname );
  int ndig=0;
  for( m_it d_it=m_range.first; d_it!=m_range.second; d_it++ )
  {
    CsRICH1Detector::DecodeChipDigit(*d_it->second);
    ndig++;
    //^std::cout << "DecodeChipDigits " << ndig << endl;
  }  
  if( debug ) cout << " RI01P digits " << ndig << endl;

  //^std::cout << "CsRICH1UpGrade::DecodeChipDigits " << digits.size() << endl;
  for( unsigned i=0; i < NCathodes(); i++ )
  {
    //^std::cout << "DecodeChipDigits " << cathodes_[i]->GetTBname() << endl;
    m_range = digits.equal_range( cathodes_[i]->GetTBname() );
    ndig = 0;
    for( m_it d_it=m_range.first; d_it!=m_range.second; d_it++ )
    {
      DecodeChipDigitPD( i, *d_it->second);
      ndig++;
      //^std::cout << "DecodeChipDigits " << cathodes_[i]->GetTBname() << " " << ndig << endl;
    }  
//     if( i == 12 ) cout << cathodes_[i]->GetTBname() << "  digits " << ndig << endl;
  }

  setDecodingDone();
  Reconstruction();
}

//////////////////////////////////////////////////////////////////////////////////

void CsRICH1UpGrade::DecodeChipDigitPD(int ic, const CS::Chip::Digit &digit)
{
  //^std::cout << " DecodeChipDigitPD" << std::endl;
  bool debug = false;
//   bool debug = true;
//   const CS::ChipGassiplex::Digit *d = dynamic_cast<const CS::ChipGassiplex::Digit *>(&digit);
  const CS::ChipAPVRICH::Digit *d = dynamic_cast<const CS::ChipAPVRICH::Digit *>(&digit);
//   if( d==NULL )
//     throw CS::Exception("CsRICH1Detector::DecodeRawData(): Wrong digit type");
  if( d!=NULL )
  {
    const uint32*  amplitude = d->GetAmplitude();
    CathodeAPV * c = reinterpret_cast<CathodeAPV *>(cathodes_[ic]);
    if( debug )
    { 
      cout << " Is sparsifed " << d->IsSparsifed() << endl;
      cout << " a0 =" << amplitude[0] <<" a1 =" << amplitude[1] << " a2 =" << amplitude[2] << endl;
    }


    if( debug )
    { 
      cout << " PDx " << d->GetPDx()  << " PDy " << d->GetPDy()  <<
              " PixelX " << d->GetPixelX()  << " PixelY " << d->GetPixelY()  <<endl;
    }
    double amp[10];
    amp[0] =  amplitude[2];
    amp[1] =  1.;
    for( unsigned itdat=0; itdat < 3; itdat++)
    {
      amp[2+itdat]=(double)( amplitude[itdat] );
    }
    int ix = d->GetPixelX();
    int iy = d->GetPixelY();
    
//     if( cathodes_[i]->GetReadOutType() == CsRICH1UpGrade::CathodePlane::GASSIPLEX ) data_size = 1;
//     if( cathodes_[i]->GetReadOutType() == CsRICH1UpGrade::CathodePlane::MAPMT ) data_size = 2;
//     if( cathodes_[ic]->GetReadOutType() == CsRICH1UpGrade::CathodePlane::APV ) data_size = 5;
    int data_size = 5;
    
//     cout << " Digit " << myDigits_.size() << " Amp " << amp[0] << " ic " << ic << " ix " << ix << " iy " << iy << endl;
    if( ix< 0 || iy < 0 )
    {
//       cerr << " BAD Digit Amp= " << amp[0] << " ic " << ic << " ix " << ix << " iy " << iy << endl;
    }
    else
    {
      if( c->GetRDreadoutDirX() < 0 ) ix = 71-ix;
      if( c->GetRDreadoutDirY() < 0 ) iy = 71-iy;
    
      CsDigit* digit = new CsDigit( *this, setPadADR( ic+1, ix+1, iy+1 ), amp, data_size );
      myDigits_.push_back( digit );
      //^std::cout << " DecodeChipDigitPD APV" << myDigits_.size() << std::endl;
    }
  }
  else
  {
// other RICH digits...
  }
  /* Andrea ... Modifications to include T0 */
  const CS::ChipF1::Digit *p = dynamic_cast<const CS::ChipF1::Digit*>(&digit);
  if( p!=NULL ) {
    //^int32 time = p->GetTime();
    string name = p->GetDetID().GetName();
    double gtu  = p->GetTimeUnit();
    double time = p->GetTimeDecoded();
    int32 ampli = p->GetAmplitude();
    int32 ix = p->GetX();
    int32 iy = p->GetY();
    const CS::ChipF1::DataID  &id = p->GetDataID();


    bool t0Cal= CsRICH1Detector::T0_flag();

    if( t0Cal ) {
      unsigned int gid = id.u.s.geoID_or_port;
      unsigned int gpo = id.u.s.src_id      ;
      unsigned int gcp = id.u.s.chip_chan>>3;
      unsigned int gch = id.u.s.chip_chan&7;
      bool         gmm = id.u.s.mode;

      unsigned int geo = gid;
      if(!gmm) geo = (gpo-500)*16+gid;

      double T0_cal = CsRICH1Detector::getPMT_T0( geo, gcp, gch )*gtu;

      if( debug ) {
	std::cout.setf(ios::fixed);
	std::cout << name 
		  << "  mode " << setprecision (0) << setw( 3 ) << gmm
		  << " gid-p " << setprecision (0) << setw( 4 ) << gid
		  << "  src  " << setprecision (0) << setw( 4 ) << gpo
		  << "  geo  " << setprecision (0) << setw( 4 ) << geo
		  << "  CHIP " << setprecision (0) << setw( 4 ) << gcp
		  << "  Chan." << setprecision (0) << setw( 4 ) << gch
		  << "  ix   " << setprecision (0) << setw( 4 ) << ix 
		  << "  iy   " << setprecision (0) << setw( 4 ) << iy  
		  << "  time " << setprecision (3) << setw( 8 ) << time 
		  << "  T0   " << setprecision (0) << setw( 6 ) << T0_cal
		  << "  amp  " << setprecision (0) << setw( 6 ) << ampli 
		  << "  t-T0 " << setprecision (0) << setw( 6 ) << time-T0_cal 
		  <<  std::endl;
      }
      
      double amp[3];
      int data_size = 3;
      amp[0] = ampli;
      amp[1] = time - T0_cal;
      amp[2] = T0_cal;
      CsDigit* digit = 
	new CsDigit( *this, setPadADR( ic+1, ix+1, iy+1 ), amp, data_size );
      myDigits_.push_back( digit );

    } else {

      double amp[2];
      int data_size = 2;
      amp[0] = ampli;
      amp[1] = time;
      CsDigit* digit = 
	new CsDigit( *this, setPadADR( ic+1, ix+1, iy+1 ), amp, data_size );
      myDigits_.push_back( digit );

    }
    //^std::cout << " DecodeChipDigitPD other" << myDigits_.size() << std::endl;
  /* Andrea ... end of modifications to include T0 */
  }

}

//////////////////////////////////////////////////////////////////////////////////

void CsRICH1UpGrade::getMCDigits( list<CsMCHit*>& richHits )   
{
  if( !flags_.initialization_done_ )
  {
    if( CoralIsReady() )
    { 
      InitFromCoral();
    }
    else
    {
      CsErrLog::Instance()->mes(elFatal, "CsRICH1UpGrade Initialization failed." );
    }   
  }
  MakeMCResponse( richHits );
}

//////////////////////////////////////////////////////////////////////////////////

bool    CsRICH1UpGrade::CoralIsReady  ( void )
{
  CsRCDetectors* det = CsRCDetectors::Instance();
  if( det == NULL) return false; 
  return ( det->lCathodes().size() != 0);
}

//////////////////////////////////////////////////////////////////////////////////

void CsRICH1UpGrade::InitCsOpt( void )   
{
   CsOpt* opt = CsOpt::Instance();
   string str;
   if ( opt->CsOpt::getOpt( "CsRICH1UpGrade", "SelectInTimeMCHits", str ) )
   {
     if( str[0] == 'Y' )
     {
       options_.mchits_time_min_=-10.;
       options_.mchits_time_max_= 10.;
     }
   }

   if ( opt->CsOpt::getOpt( "CsRICH1UpGrade", "IgnoreMCHits", str ) )
   {
     options_.ignore_mc_hits_=true;
   }

   if ( opt->CsOpt::getOpt( "CsRICH1UpGrade", "InvertMCTimeGate", str ) )
   {
     options_.invert_mc_time_gate_=true;
   }

   if ( opt->CsOpt::getOpt( "CsRICH1UpGrade", "FillMCHito", str ) )
   {
     if( str[0] == 'Y' )
     {
       options_.fill_histo_=true;
     }
   }
   
   if ( opt->CsOpt::getOpt( "CsRICH1UpGrade", "FillHito", str ) )
   {
     if( str[0] == 'Y' )
     {
       options_.fill_histo_=true;
     }
   }

   vector<float> vec;
   if( opt->CsOpt::getOpt( "CsRICH1UpGrade", "RedFactor", vec ) ) 
   {
     for( unsigned i=0; i< NCathodes(); i++ )
     {
       if( vec[i] >= 0 ) cathodes_[i]->SetReductionFactor( double(vec[i]) );
     }
   }

//     CsOpt* opt = CsOpt::Instance();
//     bool boo;
//     string str;
//     vector<float> vec;
//     int kin = 0;
//     float flo;
//     int k = 0;
// 
//     bool lprint = false;
//     boo = opt->CsOpt::getOpt( "RICHONE", "PrintMCDig", kin );
//     if( boo ) { if( kin > 0 ) lprint = true; }
// 
//     int kboo = 0;
//     int kDoPhoGen = 1;
//     boo = opt->CsOpt::getOpt( "RICHONE", "DoPhoGen", str );
//     if( boo ) { if( str[0] == 'Y' ) kDoPhoGen = 0; kboo++; }
//     int kNoise = 1;
//     boo = opt->CsOpt::getOpt( "RICHONE", "ElecNoise", str );
//     if( boo ) { 
//       if( str[0] == 'Y' ) kNoise =  0;
//       if( str[0] == 'O' ) kNoise = -1;
//       kboo++;
//     }
//     float threshin = 0.;
//     boo = opt->CsOpt::getOpt( "RICHONE", "Threshold", flo );
//     if( boo ) { threshin = flo; kboo++; }
// 
//     int dofeedbackin = 0;
//     boo = opt->CsOpt::getOpt( "RICHONE", "DoFeedback", str );
//     if( boo ) { if( str[0] == 'Y' ) dofeedbackin = 1; kboo++; }
//     int doreflin = 0;
//     boo = opt->CsOpt::getOpt( "RICHONE", "DoReflection", str );
//     if( boo ) { if( str[0] == 'Y' ) doreflin = 1; kboo++; }
//     int doQEtestin = 0;
//     boo = opt->CsOpt::getOpt( "RICHONE", "DoQEtest", str );
//     if( boo ) { if( str[0] == 'Y' ) doQEtestin= 1; kboo++; }
//     ;
//     float pedmin = 0.;
//     float sigpin = 0.;
//     float sigmin = 0.;
//     float sigsin = 0.;
//     boo = opt->CsOpt::getOpt( "RICHONE", "PedSigma", vec );
//     if( boo ) {
//       pedmin = vec[0];
//       sigpin = vec[1];
//       sigmin = vec[2];
//       sigsin = vec[3];
//       kboo++;
//     }
//     int dosigcompin = 1;
//     boo = opt->CsOpt::getOpt( "RICHONE", "DoNSigComp", str ); //031005
//     if( boo ) { if( str[0] == 'N' ) dosigcompin = 0; kboo++; }
// 
//     float csiqein[25];
//     for( int k=0; k<25; k++ ) csiqein[k] = 0.;           //   020619
//     boo = opt->CsOpt::getOpt( "RICHONE", "CsIQuantumEff", vec );
//     if( boo ) {
//       for( k=0; k<25; k++ ) {
//         csiqein[k] = vec[k];
//       }
//       kboo++;
//     }
//     float alfain[16];
//     for( int k=0; k<16; k++ ) alfain[k] = 0.03;          //   020619
//     boo = opt->CsOpt::getOpt( "RICHONE", "Alfa", vec );
//     if( boo ) {
//       for( k=0; k<16; k++ ) {
//         alfain[k] = vec[k];
//       }
//       kboo++;
//     }
//     float phsloin[16];
//     for( int k=0; k<16; k++ ) phsloin[k] = 30.;          //   020619
//     boo = opt->CsOpt::getOpt( "RICHONE", "PHslo", vec );
//     if( boo ) {
//       for( k=0; k<16; k++ ) {
//         phsloin[k] = vec[k];
//       }
//       kboo++;
//     }
//     float ridfactin[16];
//     for( int k=0; k<16; k++ ) ridfactin[k] = 1.;          //   020619
//     boo = opt->CsOpt::getOpt( "RICHONE", "RedFactor", vec );
//     if( boo ) {
//       for( k=0; k<16; k++ ) {
//         ridfactin[k] = vec[k];
//       }
//       kboo++;
//     }
//     float timegate[16];
//     for( int k=0; k<16; k++ ) timegate[k] = 2000.;     //031005
//     boo = opt->CsOpt::getOpt( "RICHONE", "TimeGate", vec );
//     if( boo ) {
//       for( k=0; k<16; k++ ) {
//         timegate[k] = vec[k];
//       }
//       kboo++;
//     }     
//     //      cout << "kboo = " << kboo << endl;
// 
//     if( kboo < 13 ) {
//       //CsErrLog::mes( elFatal,
//       //  "CsRICH1Detector::getMCDigits() : NOT initialized !" );
//       CsErrLog::mes( elError,
//         "CsRICH1Detector::getMCDigits() : default initialization !" );
//     }
// 
//     if( lprint ) {
//       cout << endl;
//       cout << "CsRICH1Detector::getMCDigits() : MC Digitization Constants"
//            << endl
//            << "----------------------------------------------------------"
//            << endl;
//       cout << "kDoPhoGen = " << kDoPhoGen 
//            << ", kNoise = " << kNoise
//            << ",  theshold = " << threshin << endl;
//       cout << "dofeedback = " << dofeedbackin
//            << ",  dorefl = " << doreflin 
//            << ",  doQEtest = " << doQEtestin << endl;
//       cout << "Peds & Sigmas = " << pedmin << ",  " << sigpin
//            << ",  " << sigmin << ",  " << sigsin << endl;
//       cout << "dosigcomp = " << dosigcompin<<endl;
//       cout << "CsIQEff = ";
//       for( k=0; k<25; k++ ) { cout << csiqein[k] << "  "; }
//       cout << endl;
//       cout << "Alfa = ";
//       for( k=0; k<16; k++ ) { cout << alfain[k] << "  "; }
//       cout << endl;
//       cout << "PHslo = ";
//       for( k=0; k<16; k++ ) { cout << phsloin[k] << "  "; }
//       cout << endl;
//       cout << "RedFactor = ";
//       for( k=0; k<16; k++ ) { cout << ridfactin[k] << "  "; }
//       cout << endl;
//       cout << "TimeGate = ";
//       for( k=0; k<16; k++ ) { cout << timegate[k] << "  "; }
//       cout << endl;
//     }

}

//////////////////////////////////////////////////////////////////////////////////

void    CsRICH1UpGrade::InitFromCoral ( void )
{
//   bool debug = true;
  bool debug = false;
  if( !CoralIsReady() )
  {
    string str = "CsRICH1UpGrade::getMCDigits() : wrong detectors.dat for RICH1";
    CsErrLog::Instance()->mes( elFatal, str );
  }
  if( debug ) cout << " CsRICH1UpGrade::InitFromCoral  NCathodes()= " << NCathodes() << endl;
  for( unsigned i=0; i < NCathodes(); i++ )
  {
    if( debug ) cout << " InitFromCoral cathode plane " << i << endl; 
    cathodes_[i]->InitFromCoral();
  }
  flags_.initialization_done_=true;
}

//////////////////////////////////////////////////////////////////////////////////

void    CsRICH1UpGrade::StoreMCHits ( list<CsMCHit*>& richHits, double tmin,  double tmax )
{
//   bool debug = true;
  bool debug = false;
  if( debug ) cout << " CsRICH1UpGrade::StoreMCHits tmin = " << tmin << " tmax = " << tmax << endl;
  int kHit = 0;
  list<CsMCHit*>::const_iterator ih;
  for( ih=richHits.begin(); ih!=richHits.end(); ih++ ) 
  {
    CsMCRICH1Hit* hit = (CsMCRICH1Hit*)(*ih);
    int icth=hit->getCathode();
    double htim = hit->getDTime();
    if( options_.printinfo_ )
    {
      cout << " kHit " << kHit << " htim " << htim << endl;
      cout << " Cathod " << icth  <<
       " xHit " << hit->getYdetDRS() << " yHit " << hit->getZdetDRS() <<
       " eHit " << hit->getPhotEnergy() << endl;
    }      
    kHit++;
//     if( htim  < tmax &&  htim  > tmin && fabs( htim ) > 1.e-05 )
    if( options_.invert_mc_time_gate_ ) //  to see effect from pile-up hits
    {
      if( htim  <= tmin ||  htim  >= tmax  )
      {
        mchits_.push_back( hit );
      }
    }
    else               //   select "in-time" hits
    {
      if( htim  < tmax &&  htim  > tmin  )
      {
        mchits_.push_back( hit );
      }
    }
  }

  for( unsigned ih=0; ih < NMCHits(); ih++ )
  {
    int ic = mchits_[ih]->getCathode()-1;
    if( ic >= 0 && ic < (int)NCathodes() )
    {
      cathodes_[ic]->AddMCHit(mchits_[ih]);
    }
    else
    {
      cerr << " Bad cathode index " << ic << endl;
      exit(1);
    }
  }

  if( options_.printinfo_ )
  {
    std::cout << " CsRICH1UpGrade MakeMCDecoding Selected N MC hits = " << NMCHits() << std::endl;
  }
  if( options_.printinfo_ )
    std::cout << " CsRICH1UpGrade MakeMCDecoding done " << std::endl;
}

//////////////////////////////////////////////////////////////////////////////////
  
void  CsRICH1UpGrade::Clear ( void )
{

  mchits_.clear();
  for( unsigned i=0; i < NCathodes(); i++ )
  {
    cathodes_[i]->Clear();
  }

  if( options_.printinfo_ )
    std::cout << " CsRICH1UpGrade cleaning done " << std::endl;
}

//////////////////////////////////////////////////////////////////////////////////

void    CsRICH1UpGrade::MakeMCResponse    ( list<CsMCHit*>& richHits )
{
  Clear(); 
  if( !options_.ignore_mc_hits_ ) 
    StoreMCHits( richHits, options_.mchits_time_min_, options_.mchits_time_max_ );
  MakeMCResponse();
  FillHisto();
}

//////////////////////////////////////////////////////////////////////////////////

void    CsRICH1UpGrade::MakeMCResponse    ( void )
{
//   bool debug = true;
  bool debug = false;
  bool fasttrack = false;
  if( fasttrack )
  {
    for( unsigned ih=0; ih < NMCHits(); ih++ )
    {
      int ic = mchits_[ih]->getCathode()-1;
      double x = mchits_[ih]->getYdetDRS();
      double y = mchits_[ih]->getZdetDRS();
      double time = mchits_[ih]->getDTime();
      pair < int, int> pad_indexes = cathodes_[ic]->GetPadIndexes( x, y );
      if( pad_indexes.first < 0 ) continue;
      int ix = pad_indexes.first;
      int iy = pad_indexes.second;
      double amp[10];
      int data_size = 1;
      amp[0] = 10.;
      CsMCDigit* digit = new CsMCDigit( *this, setPadADR( ic+1, ix+1, iy+1 ), amp, data_size );
      myDigits_.push_back( reinterpret_cast<CsDigit*>(digit) );
      digit->addHit(*mchits_[ih]);
    }
  }
  else
  {
    for( unsigned i=0; i < NCathodes(); i++ )
    {
      if( debug ) cout << " CsRICH1UpGrade::MakeMCResponse cathode_[" << i << "]->MakeMCResponse()" << endl;
      cathodes_[i]->MakeMCResponse();
      if( debug ) cout << " Got NPadsFired() in " << i << " = " << cathodes_[i]->NPadsFired() << endl;
      for( unsigned ip=0; ip < cathodes_[i]->NPadsFired(); ip++ )
      {
        CathodePAD &pad = cathodes_[i]->GetPads()[ip];
        int ic =  pad.cathode_id_;
        int ix =  pad.ix_;
        int iy =  pad.iy_;
        double amp[10];
        int data_size = 1;
        amp[0] =  pad.amp_;
        amp[1] =  pad.time_;
        vector<int> &data = pad.el_digits_;
        for( unsigned itdat=0; itdat < pad.el_digits_.size(); itdat++)
        {
          amp[2+itdat]=(double)( pad.el_digits_[itdat] );
        }
        if( cathodes_[i]->GetReadOutType() == CsRICH1UpGrade::CathodePlane::GASSIPLEX ) data_size = 1;
        if( cathodes_[i]->GetReadOutType() == CsRICH1UpGrade::CathodePlane::MAPMT ) data_size = 2;
        if( cathodes_[i]->GetReadOutType() == CsRICH1UpGrade::CathodePlane::APV ) data_size = 5;

        CsMCDigit* digit = new CsMCDigit( *this, setPadADR( ic+1, ix+1, iy+1 ), amp, data_size );
        myDigits_.push_back( reinterpret_cast<CsDigit*>(digit) );
        list < CsMCRICH1Hit* > &hitsass = cathodes_[i]->GetPadMCHits( ix, iy);
        for( list< CsMCRICH1Hit*>::iterator it=hitsass.begin(); it!=hitsass.end(); it++ )
        {
          (digit)->addHit(**it);
        }
      }      
    }
  }

}

//////////////////////////////////////////////////////////////////////////////////

void    CsRICH1UpGrade::FillHisto    ( void )
{
  FillMCDecodingHisto();
}

//////////////////////////////////////////////////////////////////////////////////

void    CsRICH1UpGrade::Analysis    ( void )
{
  FillMCAnalysisHisto();
}

// //////////////////////////////////////////////////////////////////////////////////
// 
// void    CsRICH1UpGrade::BookHisto    ( void )
// {
//   bool debug = false;
//   if( flags_.histo_booked_) return;
//   if( debug ) cout << " CsRICH1UpGrade::BookHisto " << endl;
//   CsHistograms::SetCurrentPath("/CsRICH1UpGrade/MCHits");
// 
//   histo_.h1_MC_Chambers = new  CsHist1D("MC_Chambers"," MC Chambers ",32,0.,32.);
//   histo_.h1_MC_Time = new  CsHist1D("MC_Time"," MC Time ",10000,-2500.,2500.);
//   histo_.h1_MC_PreciseTime = new  CsHist1D("MC_PreciseTime"," MC Precise Time ",20000,-10.,10.);
//   histo_.h1_MC_PhotonE = new  CsHist1D("MC_PhotonE"," MC Photon E ",1000,0.,10.);
//   histo_.h2_MC_ChambersXY = new  CsHist2D("MC_ChambersXY"," MC Chambers XY ",800,-1600.,1600.,800,-1600.,1600.);
//   histo_.h2_MC_ChambersXYProblems = new  CsHist2D("MC_ChambersXYProblems"," MC Chambers XY with pad problems ",800,-1600.,1600.,800,-1600.,1600.);
//   char hist_name[132];
//   for( unsigned ich=0; ich< NCathodes(); ich++ )
//   {
//     sprintf(hist_name,"Chamber_%d",ich);
//     histo_.h2_MC_ChamberXY[ich] = new  CsHist2D(hist_name,hist_name,560,-280.,280.,560,-280.,280.);
//     sprintf(hist_name,"ChamberTime_%d",ich);
//     histo_.h1_MC_ChamberTime[ich] = new  CsHist1D(hist_name,hist_name,10000,-2500.,2500.);
//     sprintf(hist_name,"ChamberPreciseTime_%d",ich);
//     histo_.h1_MC_ChamberPreciseTime[ich] = new  CsHist1D(hist_name,hist_name,20000,-10.,10.);
//   } 
// 
//   histo_.h1_RD_Amp = new  CsHist1D("RD_Amp"," MC Photon E ",1000,0.,1000.);;
//   histo_.h2_RD_ChambersXY = new  CsHist2D("RD_ChambersXY"," RD Chambers XY ",800,-1600.,1600.,800,-1600.,1600.);;
// 
//   if( debug ) cout << " CsRICH1UpGrade::BookHisto OK " << endl;
//   flags_.histo_booked_= true;
// }
// 
//////////////////////////////////////////////////////////////////////////////////

void    CsRICH1UpGrade::BookMCHisto    ( void )
{
  bool debug = false;
  if( debug ) cout << " CsRICH1UpGrade::BookHisto " << endl;
  CsHistograms::SetCurrentPath("/CsRICH1UpGrade/MCHits");

  histo_.h1_MC_Chambers = new  CsHist1D("MC_Chambers"," MC Chambers ",32,0.,32.);
  histo_.h1_MC_Time = new  CsHist1D("MC_Time"," MC Time ",10000,-2500.,2500.);
  histo_.h1_MC_PreciseTime = new  CsHist1D("MC_PreciseTime"," MC Precise Time ",20000,-10.,10.);
  histo_.h1_MC_PhotonE = new  CsHist1D("MC_PhotonE"," MC Photon E ",1000,0.,10.);
  histo_.h2_MC_ChambersXY = new  CsHist2D("MC_ChambersXY"," MC Chambers XY ",800,-1600.,1600.,800,-1600.,1600.);
  histo_.h2_MC_ChambersXYProblems = new  CsHist2D("MC_ChambersXYProblems"," MC Chambers XY with pad problems ",800,-1600.,1600.,800,-1600.,1600.);
  char hist_name[132];
  for( unsigned ich=0; ich< NCathodes(); ich++ )
  {
    sprintf(hist_name,"Chamber_%d",ich);
    histo_.h2_MC_ChamberXY[ich] = new  CsHist2D(hist_name,hist_name,560,-280.,280.,560,-280.,280.);
    sprintf(hist_name,"ChamberTime_%d",ich);
    histo_.h1_MC_ChamberTime[ich] = new  CsHist1D(hist_name,hist_name,10000,-2500.,2500.);
    sprintf(hist_name,"ChamberPreciseTime_%d",ich);
    histo_.h1_MC_ChamberPreciseTime[ich] = new  CsHist1D(hist_name,hist_name,20000,-10.,10.);
  } 

  histo_.h1_RD_Amp = new  CsHist1D("RD_Amp"," MC Photon E ",1000,0.,1000.);;
  histo_.h2_RD_ChambersXY = new  CsHist2D("RD_ChambersXY"," RD Chambers XY ",800,-1600.,1600.,800,-1600.,1600.);;

  if( debug ) cout << " CsRICH1UpGrade::BookHisto OK " << endl;
}

//////////////////////////////////////////////////////////////////////////////////

void    CsRICH1UpGrade::BookHisto    ( void )
{
  if( flags_.histo_booked_) return;
  BookMCHisto();
  BookRDHisto();
  flags_.histo_booked_= true;
}

//////////////////////////////////////////////////////////////////////////////////

void    CsRICH1UpGrade::BookRDHisto    ( void )
{
  bool debug = false;
  if( debug ) cout << " CsRICH1UpGrade::BookHisto " << endl;
  CsHistograms::SetCurrentPath("/CsRICH1UpGrade/Digits");

  histo_.h1_Digit_Amp = new  CsHist1D("Digit_Amp"," MC Photon E ",1000,0.,1000.);;
  histo_.h2_Digit_ChambersXY = new  CsHist2D("Digit_ChambersXY"," RD Chambers XY ",800,-1600.,1600.,800,-1600.,1600.);;

  if( debug ) cout << " CsRICH1UpGrade::BookHisto OK " << endl;
}

//////////////////////////////////////////////////////////////////////////////////

void    CsRICH1UpGrade::FillMCAnalysisHisto    ( void )
{
//   if( analyser_ != NULL )
//   { 
//     analyser_->FillMCHisto();
//     analyser_->FillMCHistoSpecial();
//   }
}

// //////////////////////////////////////////////////////////////////////////////////
// 
// void    CsRICH1UpGrade::FillMCDecodingHisto    ( void )
// {
//   bool debug = false;
// //   bool debug = true;
//   if( debug ) cout << " CsRICH1UpGrade::FillMCDecodingHisto " << endl;
//   if( ! options_.fill_histo_ ) return;
//   static double cathodes_x0[16];
//   static double cathodes_y0[16];
//   if( !flags_.histo_booked_ )
//   {
//     BookHisto();
//     double xsize = 600.;
//     double ysize = 600.;
// 
//     int ic = 0;
//     for( unsigned updown=0; updown < 2; updown++ )
//     {
//       for( unsigned col=0; col < 4; col++ )
//       {
//         double x = xsize + xsize/2. - xsize*( (double)col );
//         for( unsigned row=0; row < 2; row++ )
//         {        
//           double y = ysize + ysize/2. - ysize*2*( (double)updown ) - ysize*( (double)row );
//           cathodes_x0[ic]=x;
//           cathodes_y0[ic]=y;
//           ic++;
//         }
//       }
//     }
//   }
// 
//   for( unsigned ich=0; ich < NCathodes(); ich++ )
//   {
//     cathodes_[ich]->FillMCDecodingHisto();
//   }
// 
//   if( debug ) cout << " Start filling  " << NMCHits() << endl;
//   for( unsigned ih=0; ih < NMCHits(); ih++ )
//   {
//     int icatod = mchits_[ih]->getCathode()-1;
//     histo_.h1_MC_Chambers->Fill( double(icatod) );
//     if( icatod < 0 || icatod > 15 )
//     {
//       CsErrLog::Instance()->mes(elFatal, "CsRICH1UpGrade::FillMCDecodingHisto Bad cathod index" );
//     }  
//     
//     double xcathod = cathodes_x0[icatod];
//     double ycathod = cathodes_y0[icatod];
//     
//     double xhit = mchits_[ih]->getYdetDRS();
//     double yhit = mchits_[ih]->getZdetDRS();
//     double time = mchits_[ih]->getDTime();
//     double phe = mchits_[ih]->getPhotEnergy();
// 
// // the following is a quick and dirty fix to fill the top level histogram for MC...
// //    double xhit_global = xhit + xcathod;
// //    double yhit_global = yhit + ycathod;
//     double xhit_global = mchits_[ih]->getX() ;
//     double yhit_global = mchits_[ih]->getY();
//     double chadist = 1300;  // just a guess...
//     if(yhit_global < 0) yhit_global = yhit_global + chadist;
//     if(yhit_global > 0) yhit_global = yhit_global - chadist;
// // end of my dirty changes, Andreas Mutter, 3.3.06
//     
//     pair < int, int> padindex =  cathodes_[icatod]->GetPadIndexes( xhit, yhit );
//     if( padindex.first < 0 )
//     {
//       histo_.h2_MC_ChambersXYProblems -> Fill( xhit_global, yhit_global);
//     }
// 
//     histo_.h1_MC_Time->Fill( time );
//     histo_.h1_MC_PreciseTime->Fill( time );
//     histo_.h1_MC_PhotonE->Fill( phe );
// 
//     histo_.h2_MC_ChambersXY -> Fill( xhit_global, yhit_global);
// 
//     if( icatod >=0 && (unsigned)icatod < NCathodes() )
//     {
//       histo_.h1_MC_ChamberTime[icatod]->Fill( time );
//       histo_.h1_MC_ChamberPreciseTime[icatod]->Fill( time );
//       histo_.h2_MC_ChamberXY[icatod]->Fill( xhit, yhit ); 
//     }
//   }
// 
//   for( unsigned ich=0; ich< NCathodes(); ich++ )
//   {
//     for( unsigned it=0; it< cathodes_[ich]->NPadsFired(); it++ )
//     {
//       CathodePAD &pad = cathodes_[ich]->GetPads()[it];
//       int icatod =  pad.cathode_id_; // Must be the same as ich
//       int ixpad =  pad.ix_;
//       int iypad =  pad.iy_;
//       double amp =  pad.amp_;
//       double time =  pad.time_;
//       pair < double, double> pad_position = cathodes_[icatod]->GetPadPosition( ixpad, iypad );
//      
//       double xpad = pad_position.first;  
//       double ypad = pad_position.second;  
//      
//       double xcathod = cathodes_x0[icatod];
//       double ycathod = cathodes_y0[icatod];
//      
//       double xpad_global = xpad + xcathod; 
//       double ypad_global = ypad + ycathod; 
//  
//       histo_.h1_RD_Amp->Fill( amp );
//       histo_.h2_RD_ChambersXY ->Fill( xpad_global, ypad_global);
//     }
//   }
// 
//   if( debug ) cout << " CsRICH1UpGrade::FillMCDecodingHisto OK " << NMCHits() << endl;
// 
// }
// 
//////////////////////////////////////////////////////////////////////////////////

void    CsRICH1UpGrade::FillMCDecodingHisto    ( void )
{
  bool debug = false;
//   bool debug = true;
  if( debug ) cout << " CsRICH1UpGrade::FillMCDecodingHisto " << endl;
  if( ! options_.fill_histo_ ) return;
  BookHisto();

  for( unsigned ich=0; ich < NCathodes(); ich++ )
  {
    cathodes_[ich]->FillMCDecodingHisto();
  }

  if( debug ) cout << " Start filling  " << NMCHits() << endl;
  for( unsigned ih=0; ih < NMCHits(); ih++ )
  {
    int icatod = mchits_[ih]->getCathode()-1;
    histo_.h1_MC_Chambers->Fill( double(icatod) );
    if( icatod < 0 || icatod > 15 )
    {
      CsErrLog::Instance()->mes(elFatal, "CsRICH1UpGrade::FillMCDecodingHisto Bad cathod index" );
    }  
    
    double xcathod = cathodes_[icatod]->x0_image_;
    double ycathod = cathodes_[icatod]->y0_image_;
    
    double xhit = mchits_[ih]->getYdetDRS();
    double yhit = mchits_[ih]->getZdetDRS();
    double time = mchits_[ih]->getDTime();
    double phe = mchits_[ih]->getPhotEnergy();
    double xhit_global = xhit + xcathod;
    double yhit_global = yhit + ycathod;
    
    
    pair < int, int> padindex =  cathodes_[icatod]->GetPadIndexes( xhit, yhit );
    if( padindex.first < 0 )
    {
      histo_.h2_MC_ChambersXYProblems -> Fill( xhit_global, yhit_global);
    }

    histo_.h1_MC_Time->Fill( time );
    histo_.h1_MC_PreciseTime->Fill( time );
    histo_.h1_MC_PhotonE->Fill( phe );

    histo_.h2_MC_ChambersXY -> Fill( xhit_global, yhit_global);

    if( icatod >=0 && (unsigned)icatod < NCathodes() )
    {
      histo_.h1_MC_ChamberTime[icatod]->Fill( time );
      histo_.h1_MC_ChamberPreciseTime[icatod]->Fill( time );
      histo_.h2_MC_ChamberXY[icatod]->Fill( xhit, yhit ); 
    }
  }

  for( unsigned ich=0; ich< NCathodes(); ich++ )
  {
    for( unsigned it=0; it< cathodes_[ich]->NPadsFired(); it++ )
    {
      CathodePAD &pad = cathodes_[ich]->GetPads()[it];
      int icatod =  pad.cathode_id_; // Must be the same as ich
      int ixpad =  pad.ix_;
      int iypad =  pad.iy_;
      double amp =  pad.amp_;
      double time =  pad.time_;
      pair < double, double> pad_position = cathodes_[icatod]->GetPadPosition( ixpad, iypad );
     
      double xpad = pad_position.first;  
      double ypad = pad_position.second;  
     
      double xpad_global = xpad + cathodes_[ich]->x0_image_; 
      double ypad_global = ypad + cathodes_[ich]->y0_image_; 
 
      histo_.h1_RD_Amp->Fill( amp );
      histo_.h2_RD_ChambersXY ->Fill( xpad_global, ypad_global);
    }
  }

  if( debug ) cout << " CsRICH1UpGrade::FillMCDecodingHisto OK " << NMCHits() << endl;

}

//////////////////////////////////////////////////////////////////////////////////

void    CsRICH1UpGrade::Reconstruction    ( void )
{
  FillDecodingHisto();
}

//////////////////////////////////////////////////////////////////////////////////

void    CsRICH1UpGrade::FillDecodingHisto    ( void )
{
  bool debug = false;
//   bool debug = true;
  if( debug ) cout << " CsRICH1UpGrade::FillDecodingHisto options_.fill_histo_ = " << options_.fill_histo_ << endl;
  if( ! options_.fill_histo_ ) return;
  BookHisto();

// Try with external access
  const std::list<CsDigit*> &mydigs = (this)->getMyDigits();
  for( list<CsDigit*>::const_iterator id=mydigs.begin(); id!=mydigs.end(); id++ )
  {
    int adr = (*id)->getAddress();
    int icatod = CathodePlane::DecodeIndxCathode( adr );
    int ixpad = CathodePlane::DecodeIndxPadX( adr );
    int iypad = CathodePlane::DecodeIndxPadY( adr );
    if( icatod <0 || icatod > 15 )
    {
       cerr << " Bad cathode number " << icatod << endl;
       exit(1);
    }
    if( ixpad <0 || ixpad > 71 )
    {
       cerr << " Bad ixpad number " << ixpad << endl;
       exit(1);
    }
    if( iypad <0 || iypad > 71 )
    {
       cerr << " Bad iypad number " << iypad << endl;
       exit(1);
    }
    
    double* data = (*id)->getData();
    int data_size = (*id)->getDataSize();
    double amp = data[0];

    pair < double, double> pad_position = cathodes_[icatod]->GetPadPosition( ixpad, iypad );     
    double xpad = pad_position.first;  
    double ypad = pad_position.second;  
    double xpad_global = xpad + cathodes_[icatod]->x0_image_; 
    double ypad_global = ypad + cathodes_[icatod]->y0_image_; 

    histo_.h1_Digit_Amp->Fill( amp );
    histo_.h2_Digit_ChambersXY ->Fill( xpad_global, ypad_global);
  }

  if( debug ) cout << " CsRICH1UpGrade::FillDecodingHisto OK " << NMCHits() << endl;

}

//////////////////////////////////////////////////////////////////////////////////

CsRICH1UpGrade::CathodePlane::CathodePlane ( int id ) : 
  id_(id),
  read_out_type_(IDEAL)
{
  bool debug = false;
//   bool debug = true;
  if( debug ) cout << " Construct  CathodePlane id " << id << endl;
  
  nxpads_ = 72;
  nypads_ = 72;
  npads_= nxpads_*nypads_;

  pad_size_=8.;

  size_x_= nxpads_*pad_size_;
  size_y_= nypads_*pad_size_;

  pad0_x_position_=-size_x_/2.+pad_size_/2.;
  pad0_y_position_=-size_y_/2.+pad_size_/2.;
  
  x0_image_ = 0.;
  y0_image_ = 0.;  
  
  if( CoralIsReady() ) InitFromCoral(); 
}

//////////////////////////////////////////////////////////////////////////////////

CsRICH1UpGrade::CathodePlane::CathodePlane ( int id, double center_position[], double dir_normal[]) : 
  id_(id),
  efficiency_reduction_factor_(1.)
{
  bool debug = false;
//   bool debug = true;
  if( debug ) cout << " Construct  CathodePlane id " << id << endl;
  
  nxpads_ = 72;
  nypads_ = 72;
  npads_= nxpads_*nypads_;

  pad_size_=8.;

  center_position_[0]=center_position[0];
  center_position_[1]=center_position[1];
  center_position_[2]=center_position[2];
  
  size_x_= nxpads_*pad_size_;
  size_y_= nypads_*pad_size_;

  pad0_x_position_=-size_x_/2.+pad_size_/2.;
  pad0_y_position_=-size_y_/2.+pad_size_/2.;
  
  dir_normal_[0]=dir_normal[0];
  dir_normal_[1]=dir_normal[1];
  dir_normal_[2]=dir_normal[2];

  dir_x_[0]=1.;
  dir_x_[1]=0.;
  dir_x_[2]=0.;

  dir_y_[0]=0.;
  dir_y_[1]=-dir_normal[2];
  dir_y_[2]= dir_normal[1];

}

//////////////////////////////////////////////////////////////////////////////////

void CsRICH1UpGrade::CathodePlane::InitCsOpt( void )
{
   CsOpt* opt = CsOpt::Instance();
}

//////////////////////////////////////////////////////////////////////////////////

void    CsRICH1UpGrade::CathodePlane::InitFromCoral( void )
{
//   bool debug = true;
  bool debug = false;
  if( !CsRICH1UpGrade::CoralIsReady() )
  {
    CsErrLog::Instance()->mes(elFatal, "No cathodes in CsRCDetectors::Instance() CsRICH1UpGrade::CathodePlane::InitFromCoral FAILED" );
  }
  CsRCDetectors* det = CsRCDetectors::Instance();
  CsRCCathode* cc = det->ptrToCat( id_ );

  if( debug ) cout << " xcMRS " << cc->vCat0()[0] << " ycMRS " << cc->vCat0()[1] << " zcMRS " << cc->vCat0()[2] <<endl;

  center_position_[0] = cc->vCat0()[0];
  center_position_[1] = cc->vCat0()[1];
  center_position_[2] = cc->vCat0()[2];

  HepMatrix rotMRC = cc->rotMatrix();
  if( debug)
  { 
    cout << " Rotation matrix for cathod id=" << id_ << endl;
    for( unsigned i=0; i<3; i++ ) { 
      for( unsigned j=0; j<3; j++ ) { 
        cout << " ("<<i<<","<<j<<")="<< rotMRC(i+1,j+1);
      }
      cout << endl;
    }
  }
  dir_x_[0]=rotMRC(1,1);
  dir_x_[1]=rotMRC(1,2);
  dir_x_[2]=rotMRC(1,3);
  
  dir_y_[0]=rotMRC(2,1);
  dir_y_[1]=rotMRC(2,2);
  dir_y_[2]=rotMRC(2,3);
  
  dir_normal_[0]=rotMRC(3,1);
  dir_normal_[1]=rotMRC(3,2);
  dir_normal_[2]=rotMRC(3,3);

// ????????????????????????????????
//   double vcorr = 11.902;
//   
//   center_position_[0] += vcorr*dir_normal_[0];
//   center_position_[1] += vcorr*dir_normal_[1];
//   center_position_[2] += vcorr*dir_normal_[2];
  
}

//////////////////////////////////////////////////////////////////////////////////

list < CsMCRICH1Hit* > &CsRICH1UpGrade::CathodePlane::GetPadMCHits( int ixpad, int iypad )
{
  int pad_index = iypad*72 + ixpad;
  assert( (pad_index >= 0 && pad_index < (int)npads_ ));
  return padmc_hits_[pad_index]; 
}

//////////////////////////////////////////////////////////////////////////////////

void    CsRICH1UpGrade::CathodePlane::AddMCHit( CsMCRICH1Hit* hit )
{
  if( hit->getCathode()-1 == id_ )
  { 
    mchits_.push_back( hit );
  }
  else
  {
    cerr << " CsRICH1UpGrade::CathodePlane::AddMCHit wrong cathode " <<
                                   hit->getCathode()-1 << " but id " << id_ << endl;
  }
}

//////////////////////////////////////////////////////////////////////////////////

void    CsRICH1UpGrade::CathodePlane::SortMCHits( void )
{
  for ( unsigned ih=0; ih< NMCHits(); ih++ )
  {
    pair < int, int> pad_indexes = GetPadIndexes( mchits_[ih]->getYdetDRS(), mchits_[ih]->getZdetDRS() );
    if( pad_indexes.first < 0 ) continue;
    int ixpad = pad_indexes.first;
    int iypad = pad_indexes.second;
    list < CsMCRICH1Hit* > &lmcp = GetPadMCHits( ixpad, iypad); 
    lmcp.push_back(mchits_[ih]);
  }
}

//////////////////////////////////////////////////////////////////////////////////

void    CsRICH1UpGrade::CathodePlane::MakeMCResponse( void )
{
  SortMCHits();
  for( unsigned ih=0; ih < NMCHits(); ih++ )
  {
    int ic = mchits_[ih]->getCathode()-1;
    double x = mchits_[ih]->getYdetDRS();
    double y = mchits_[ih]->getZdetDRS();
    double time = mchits_[ih]->getDTime();
    pair < int, int> pad_indexes = GetPadIndexes( x, y );
    if( pad_indexes.first < 0 ) continue;
    pads_.push_back( CathodePAD( ic, pad_indexes.first, pad_indexes.second, 10., time) );
  }
}

//////////////////////////////////////////////////////////////////////////////////

void    CsRICH1UpGrade::CathodePlane::Clear( void )
{
  mchits_.clear(); 
  pads_.clear();
  padclusters_.clear();
  for( unsigned i=0; i < npads_; i++ )
  {
    padmc_hits_[i].clear();
  }
}

//////////////////////////////////////////////////////////////////////////////////

void    CsRICH1UpGrade::CathodePlane::MakeClusters( void )
{
  for( unsigned i=0; i < NPadsFired(); i++ )
  {
    CathodePAD &pad = pads_[i];
//    int icatod =  pad.cathode_id_;
    int ixpad =  pad.ix_;
    int iypad =  pad.iy_;
    double amp =  pad.amp_;
    double time =  pad.time_;
    pair < double, double> pad_position = GetPadPosition( ixpad, iypad );
     
    double xpad = pad_position.first;  
    double ypad = pad_position.second;  
    padclusters_.push_back( CsRICH1UpGrade::Cluster( xpad, ypad, amp, time ) );
  }
}

//////////////////////////////////////////////////////////////////////////////////

int  CsRICH1UpGrade::CathodePlane::SetPadAdr( int ix, int iy ) const
{
  return  ((ix+1) + ((iy+1)<<10)+ ((id_+1)<<20));
}

//////////////////////////////////////////////////////////////////////////////////

int  CsRICH1UpGrade::CathodePlane::DecodeIndxCathode( int iadr )
{
  return  ( ((iadr)>>20)-1 );
}

//////////////////////////////////////////////////////////////////////////////////

int  CsRICH1UpGrade::CathodePlane::DecodeIndxPadX( int iadr )
{
  return  ((iadr&1023)-1);
}

//////////////////////////////////////////////////////////////////////////////////

int  CsRICH1UpGrade::CathodePlane::DecodeIndxPadY( int iadr )
{
  return  (( (iadr>>10)&1023 )-1);
}

//////////////////////////////////////////////////////////////////////////////////

std::pair < int, int>  CsRICH1UpGrade::CathodePlane::GetPadIndexes( double x, double y )
{
  if( x < -size_x_/2. || x > size_x_/2. || y < -size_y_/2. || y > size_y_/2. )
  {
//     cerr << " WARNING! Bad coordinates in CsRICH1UpGrade::CathodePlane::GetPadIndexes x =" <<
//                                                                        x << " y =" << y << endl;
    return   std::pair < int, int> ( -1, -1);                                                                 
  }                                                                     
  int ixpad = (int)(( x - pad0_x_position_ + pad_size_/2.)/pad_size_); 
  int iypad = (int)(( y - pad0_y_position_+ pad_size_/2. )/pad_size_);
  if( ixpad < 0 || ixpad >= (int)nxpads_ ||iypad < 0 || iypad >= (int)nypads_ )
  {
    if( ixpad == -1 ) ixpad=0; 
    if( ixpad == 72 ) ixpad=71; 
    if( iypad == -1 ) iypad=0; 
    if( iypad == 72 ) iypad=71; 

    if( ixpad < 0 || ixpad >= (int)nxpads_ ||iypad < 0 || iypad >= (int)nypads_ )
    {
      cerr << " !!!!!!!!!!!!!!!!!!!!!!!!! ERROR Bad indexes in CsRICH1UpGrade::CathodePlane::GetPadIndexes x =" <<
                       x << " y =" << y <<" ixpad =" << ixpad << " iypad =" << iypad << endl;
      string str = " CsRICH1UpGrade::CathodePlane::GetPadIndexes : wrong x, y ";
      return std::pair < int, int> (71,71);
//      CsErrLog::Instance()->mes( elFatal, str );
    }
  }
  return std::pair < int, int> ( ixpad, iypad);
}

//////////////////////////////////////////////////////////////////////////////////

std::pair < double, double>  CsRICH1UpGrade::CathodePlane::GetPadPosition( int ixpad, int iypad )
{
  double xpad = pad0_x_position_ + ixpad*pad_size_;
  double ypad = pad0_y_position_ + iypad*pad_size_;
  return std::pair < double, double> ( xpad, ypad );
}

//////////////////////////////////////////////////////////////////////////////////

pair<double,double> CsRICH1UpGrade::RanGau( void )
{
  double u1,u2,r1,r2;
  u1=2*M_PI*drand48(), u2=sqrt(-2*log(drand48())), r1=sin(u1)*u2, r2=cos(u1)*u2;
  return pair<double,double> (r1,r2);
}

//////////////////////////////////////////////////////////////////////////////////

double CsRICH1UpGrade::Randm( void )
{
  return drand48();
}

//////////////////////////////////////////////////////////////////////////////////

double CsRICH1UpGrade::RanExp( double b )
{
  return log(1./drand48())*b;
}

//////////////////////////////////////////////////////////////////////////////////
