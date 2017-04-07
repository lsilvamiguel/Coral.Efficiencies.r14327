// $Id: DetFileInfo.cc,v 1.12 2002/09/16 17:23:30 hpereira Exp $

#include "DetFileInfo.h"
#include "Utils.h"
#include <strstream.h>
#include <fstream.h>
#include <TMath.h>

//________________________________________________________________
ClassImp(DetFileInfo)

//________________________________________________________________
DetFileInfo::DetFileInfo( string detFileName, unsigned int lNumber, string detInfoLine  ):
  TObject(),
  lineNumber_( lNumber ),
  Comment_(""), id_( 0 ), 
  TBName_(""), name_(""), unit_( 0 ),
  type_(-1),
  rdlen_( 0 ),
  xsiz_(0), ysiz_(0), zsiz_(0),
  xcm_(0), ycm_(0), zcm_(0),
  rotMNum_(-1),
  wirP_(0), ang_(0), nWir_(0), wirD_(0),
  eff_(1), bkg_(0),
  tGate_(0), vel_(0),  
  t0_(0), thRes_(0), spSli_(0), tiSli_(0), 
  rotM_(TMatrix(3,3)),
  irotM_(TMatrix(3,3)),
  res_(0),
  isMWPC_( false )

{
  //=== check detFileName
  if( !detFileName.size() ) return;

  //=== first get DRSRotM_ from DetectorFile
  ifstream in( detFileName.c_str(), ios::in );
  if( !in ) {
    cout << "DetFileInfo::DetFileInfo - FATAL: cannot read file \"" 
         << detFileName  << "\"." 
         << endl;
    return;
  }
  
  char line[512];
  TMatrix geantToCoral(3,3);
  double set[] = { 0, 1, 0, 0, 0, 1, 1, 0, 0 };
  for( int i=0; i<9; i++ ) geantToCoral( i/3, i%3 ) = set[i];
  TMatrix geantToCoralT(TMatrix::kTransposed, geantToCoral);
  
  while( !in.eof() ) { 
    in.getline( line, sizeof( line ), '\n');
    if( in.eof() || !in.good() ) continue;
    
    if( line[0] == '\0' || line[0] != ' ' ) continue;    
    if( line[1] == ' ' || line[1] == '-' || line[1] == '=' ) continue;
    
    istrstream s(line);
    string opt;
    s >> opt;

    if( opt != "rot" ) continue;
    
    TMatrix m(3,3);
    int dummy; s >> dummy;
		for( unsigned int i=0; i<3; i++ )
		for( unsigned int j=0; j<3; j++ )
		s >> m(i,j);
    TMatrix tmp = TMatrix(geantToCoral);
    tmp*=m;
    tmp*=geantToCoralT;
    m = tmp;
		DRSrotM_.push_back( tmp );
  }
  
  in.close();
  
  //=== second, parse detInfoLine to set parameters
  istrstream s( detInfoLine.c_str() );
  string opt;
  s >> opt;

  if ( opt != "det" ) {
    Comment_ = opt;
    s >> opt;
  }  

  if ( opt != "det" ) {
    cout << "DetFileInfo::DetFileInfo - FATAL: wrong line type\"" 
         << opt  << "\"." 
         << endl;
    Add( this );
    return;
  }
  
  s >> id_;
  s >> TBName_;
              
  s >> name_;     // detector name
  s >> unit_;     // detector number in station 
  s >> type_;     // detector type
  s >> rdlen_;    // radiation length
  
  //=== detector size (mm)
  s >> zsiz_;  zsiz_ *= 10;    // convert [cm] to [mm]       
  s >> xsiz_;  xsiz_ *= 10;
  s >> ysiz_;  ysiz_ *= 10;  
  
  //=== detector centre (MRS) (cm)   
  s >> zcm_; zcm_*=10;  // convert [cm] to [mm]       
  s >> xcm_; xcm_*=10; 
  s >> ycm_; ycm_*=10;    
  
  s >> rotMNum_;                          // rotation matrix number
  s >> wirD_; wirD_*=10;                  // 1st wire offset [cm]->[mm]
  s >> ang_;                              // angle of wires in DRS [deg]
  
  //=== Set rotM_ and irotM_       TMatrix rotWire(3,3);
  if( rotMNum_ < 0 || rotMNum_ > (int) DRSrotM_.size() ) {
    cout << "DetFileInfo::DetFileInfo -FATAL: unexpected matrix number " << rotMNum_ << ".\n";
    exit(0);
  }
  
  TMatrix rotWire(3,3);
  double ang_rad = ang_*TMath::Pi()/180.;
  rotWire(0,0) = cos(ang_rad); rotWire(0,1) = -sin(ang_rad); rotWire(0,2) = 0;
  rotWire(1,0) = sin(ang_rad); rotWire(1,1) =  cos(ang_rad); rotWire(1,2) = 0;
  rotWire(2,0) = 0;            rotWire(2,1) = 0;             rotWire(2,2) = 1;
            
  rotM_  = TMatrix( DRSrotM_[rotMNum_-1], TMatrix::kMult, rotWire );
  irotM_ = TMatrix( TMatrix::kInverted, rotM_ );   
  
  s >> nWir_;               // number of wires
  s >> wirP_; wirP_*=10;     // wires pitch [cm]->[mm]
  
  //=== Extra keys for drift-like detectors
  s >> eff_;      // detector efficiency
  s >> bkg_;      // detector background
  
  s >> tGate_;    // detector time gate
  
  //=== parsing ends here for non drift-like detectors
  if( type_ != 11 ) {
    res_ = wirP_/sqrt(12.0);
    Add( this );
    return;
  }

  double vel;          // velocity in mm/time_slices
  if( !( s >> vel ) ) {
    cout << "DetFileInfo::DetFileInfo - FATAL: bad line for drift like detector\"" 
         << TBName_  << "\"." 
         << endl;
    Add( this );
    return;
  }
  
  vel_ = vel*10;  // velocity mm/F1unit  
  s >> t0_;       // t0 [ns] 
  s >> thRes_;     
  s >> spSli_;     
  s >> tiSli_;     
    
  res_ = vel_*spSli_;
  Add( this );
  return;
}  

//________________________________________________________
void DetFileInfo::UpdateRotMatrices( void ) 
{
  TMatrix rotWire(3,3);
  double ang_rad = ang_*TMath::Pi()/180.;
  rotWire(0,0) = cos(ang_rad); rotWire(0,1) = -sin(ang_rad); rotWire(0,2) = 0;
  rotWire(1,0) = sin(ang_rad); rotWire(1,1) =  cos(ang_rad); rotWire(1,2) = 0;
  rotWire(2,0) = 0;            rotWire(2,1) = 0;             rotWire(2,2) = 1;
            
  rotM_  = TMatrix( DRSrotM_[rotMNum_-1], TMatrix::kMult, rotWire );
  irotM_ = TMatrix( TMatrix::kInverted, rotM_ );   

  return;
}

//________________________________________________________
unsigned int DetFileInfo::Add( DetFileInfo* dfi )
{
  sub_.push_back( dfi );
  return sub_.size();
}
  
//________________________________________________________
DetFileInfo* DetFileInfo::GetMain( void )
{
  //=== no subdetectors
  if( sub_.size() <= 1 ) return this;

  //=== Special case for straws
  DetFileInfo * out = 0;
  if( sub_[0]->TBName_.substr(0,2) == "ST" ) {
    for( unsigned int i=0; i< sub_.size(); i++ ) 
    if( sub_[i]->TBName_[7] == 'b' ) {
      out = sub_[i];
      break;
    }
  }  
  
  //=== Standard case
  else {
    char digit;
   for( unsigned int i=0; i< sub_.size(); i++ ) 
    if( i == 0 ||  sub_[i]->name_[sub_[i]->name_.size()-1] < digit )  { 
      digit = sub_[i]->name_[sub_[i]->name_.size()-1];
      out = sub_[i];
    } 
  }
  
  return out;
}

//____________________________________________________________________
// Warning: this is just a printout method. Does not modify any member
#define resout "res.out"
void DetFileInfo::DUtoCenter_CM( double du )
{
  ofstream out( resout, ios::app );
  out.form("DetFileInfo::DUtoCenter - \"%s\" [%10.3f] [ %10.4f %10.4f ] -> [ %10.4f %10.4f ].\n",
    TBName_.c_str(),
    du,
    xcm_/10,
    ycm_/10,
    (xcm_+rotM_(0,0)*du)/10,
    (ycm_+rotM_(1,0)*du)/10);
  out.close();
}
  
//________________________________________________________
string DetFileInfo::Dump( void )
{
  char* line = new char[512];
  string out("");
  sprintf( line,
    "%s det %4i  %s  %s %2i %4i %8.2f %7.3f %7.3f %7.3f %11.4f %10.4f %10.4f ",
    Comment_.c_str(),
    id_, TBName_.c_str(), name_.c_str(), unit_, type_,
    rdlen_, 
    zsiz_/10, xsiz_/10, ysiz_/10,
    zcm_/10,  xcm_/10,  ycm_/10);
  out += string(line);
  
  sprintf( line,
    "%6i %13.4f %8.3f %5i %10.6f %5.3f %5.1f %6.1f ",
    rotMNum_,
    wirD_/10, ang_, nWir_, wirP_/10,
    eff_, bkg_, tGate_ );
  out += string(line);

  if( type_ != 11 ) { delete line; return out; }
  
  sprintf( line, 
    "%7.5f %7.1f %5.1f %5.2f %8.5f",
    vel_/10, t0_, thRes_, spSli_, tiSli_ );
  out += string(line);
  delete line;
  return out;

}
