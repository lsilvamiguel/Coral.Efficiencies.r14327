// $Id: DetectorInfo.cc,v 1.30 2009/08/31 01:38:57 ybedfer Exp $

/*!
   \file    DetectorInfo.cc
   \brief   store detector.dat detector informations
   \author  Hugo Pereira
   \version $Revision: 1.30 $
   \date    $Date: 2009/08/31 01:38:57 $
*/

#include "DetectorInfo.h"
#include "Utils.h"
#include "Opt.h"
#include "Defs.h"
#include "Point.h"
#include "Obj3D.h"
#include <sstream>
#include <fstream>
#include <TMath.h>
char* operator+( std::streampos&, char* );
using namespace std;

//________________________________________________________________
ClassImp(DetectorInfo)
DetectorInfo::DetectorInfo( string detFileName, unsigned int lNumber, string detInfoLine  ):
  TObject(),
  lineNumber_( lNumber ),
  id_( 0 ), 
  TBName_(""), name_(""), unit_( 0 ),
  type_(-1),
  rdlen_( 0 ),
  zsiz_(0), xsiz_(0), ysiz_(0),
  xcm_(0), ycm_(0), zcm_(0),
  rotMNum_(-1),
  wirD_(0), ang_(0), nWir_(0), wirP_(0),
  rotM_(TMatrix(3,3)),
  irotM_(TMatrix(3,3)),
  cosTheta_(0), sinTheta_(0),
  eff_(1), bkg_(0),
  tGate_(0), vel_(0),  
  t0_(0), thRes_(0), spSli_(0), tiSli_(0), 
  res_(0),
  isMWPC_( false ),
  projection_( -1 ),
  
  dead_(0),

  //! alignment bias
  biasU_(0), 
  biasV_(0), 
  biasZ_(0),
  biasT_(0),
  biasP_(0),
  biasR_(0),
  biasL_(0),
    
  area_( 0 )
{
  // check detFileName
  if( !detFileName.size() ) return;

  // first get DRSRotM_ from DetectorFile
  ifstream in( detFileName.c_str(), ios::in );
  if( !in ) {
    cout << "DetectorInfo::DetectorInfo - FATAL: cannot read file \"" 
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
    
    istringstream s(line);
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
  
  // second, parse detInfoLine to set parameters
  istringstream s( detInfoLine.c_str() );
  string opt;
  s >> opt;

  if ( opt != "det" ) {
    cout << "DetectorInfo::DetectorInfo - FATAL: wrong line type\"" 
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
  
  // detector size (mm)
  s >> zsiz_;  zsiz_ *= 10;    // convert [cm] to [mm]       
  s >> xsiz_;  xsiz_ *= 10;
  s >> ysiz_;  ysiz_ *= 10;  
  
  // detector centre (MRS) (cm)   
  s >> zcm_; zcm_*=10;  // convert [cm] to [mm]       
  s >> xcm_; xcm_*=10; 
  s >> ycm_; ycm_*=10;    
  
  s >> rotMNum_;                          // rotation matrix number
  s >> wirD_; wirD_*=10;                  // 1st wire offset [cm]->[mm]
  s >> ang_;                              // angle of wires in DRS [deg]
  
  // Set rotM_ and irotM_       TMatrix rotWire(3,3);
  if( rotMNum_ < 0 || rotMNum_ > (int) DRSrotM_.size() ) {
    cout << "DetectorInfo::DetectorInfo -FATAL: unexpected matrix number " << rotMNum_ << ".\n";
    exit(0);
  }
  
  TMatrix rotWire(3,3);
  double ang_rad = ang_*TMath::Pi()/180.;
  rotWire(0,0) = cos(ang_rad); rotWire(0,1) = -sin(ang_rad); rotWire(0,2) = 0;
  rotWire(1,0) = sin(ang_rad); rotWire(1,1) =  cos(ang_rad); rotWire(1,2) = 0;
  rotWire(2,0) = 0;            rotWire(2,1) = 0;             rotWire(2,2) = 1;
            
  rotM_  = TMatrix( DRSrotM_[rotMNum_-1], TMatrix::kMult, rotWire );

  // irotM_ = TMatrix( TMatrix::kInverted, rotM_ );   
  // do inversion by hand (valid only for dets perp to the beam)
  irotM_ = TMatrix( rotM_ );
  irotM_(0,1) = rotM_(1,0);
  irotM_(1,0) = rotM_(0,1);
  
  
  cosTheta_ = rotM_(0,0);
  sinTheta_ = rotM_(1,0);    
  
  s >> nWir_;               // number of wires
  s >> wirP_; wirP_*=10;     // wires pitch [cm]->[mm]
  
  // Extra keys for drift-like detectors
  s >> eff_;      // detector efficiency
  s >> bkg_;      // detector background
  
  s >> tGate_;    // detector time gate
  
  // parsing ends here for non drift-like detectors
  if( type_ != 11 ) {
    res_ = wirP_/sqrt(12.0);
    Add( this );
    return;
  }

  double vel;          // velocity in mm/time_slices
  if( !( s >> vel ) ) {
    cout << "DetectorInfo::DetectorInfo - FATAL: bad line for drift like detector\"" 
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
unsigned int DetectorInfo::Add( DetectorInfo* dfi )
{
  sub_.push_back( dfi );
  return sub_.size();
}

//________________________________________________________
void DetectorInfo::SetBias( const char* biasL ) 
{
  bool isPixel=false;
  // special treatment of pixel detectors
  if ( (TBName_.substr(0,2)=="GP" ||
	TBName_.substr(0,2)=="MP") &&
       (TBName_.substr(4,1)=="P"  ||
	TBName_.substr(4,1)=="M" ) )
    isPixel=true;

  istringstream in( biasL );
  unsigned int ib=0;
  while( in ) {
    if( in.eof() ) break;
    double value; in >> value;
      
    switch ( ib ) {
      case 0: biasU_ = value;
      case 1: biasZ_ = value; break;
      case 2: biasT_ = value; break;
      case 3: biasP_ = value; break;
      case 4: 
        //! check if det ist drift-like
        if( type_ == 11 ) biasR_ = value; 
        else if ( isPixel ) biasV_ = value;
        else if( value != 0 ) cout << "DetectorInfo::SetBias \"" << TBName_ << "\". Cannot add R bias to non drift-like detectors.\n";
        break;
        
      case 5: 
        //! check if det ist drift-like
        if( type_ == 11 ) biasL_ = value; 
        else if( value != 0 ) cout << "DetectorInfo::SetBias \"" << TBName_ << "\". Cannot add L bias to non drift-like detectors.\n";
        break;
        
      default: break;
    }
    ib++;
  
  }
  
  // dump
  cout << "DetectorInfo - adding bias for \"" 
       << TBName_ 
       << "\": U+= " <<  biasU_ << "mm";
  if ( isPixel )
    cout << " V+= " << biasV_ << "mm";
  cout << " Z += " <<    biasZ_ << "mm";
  cout << " T += " <<    biasT_ << "rad";
  cout << " P *= (1+" << biasP_ << ")";
  if( type_ == 11 ) cout << " R += " <<    biasR_ << "mm";
  else cout << " R (ignored) ";
  if( type_ == 11 ) cout << " R *= (1+" <<    biasL_ << ")";
  else cout << " L (ignored) ";
  cout << ".\n";
}  

//________________________________________________________
void DetectorInfo::UpdateRotMatrices( void ) 
{
  TMatrix rotWire(3,3);
  double ang_rad = ang_*TMath::Pi()/180.;
  rotWire(0,0) = cos(ang_rad); rotWire(0,1) = -sin(ang_rad); rotWire(0,2) = 0;
  rotWire(1,0) = sin(ang_rad); rotWire(1,1) =  cos(ang_rad); rotWire(1,2) = 0;
  rotWire(2,0) = 0;            rotWire(2,1) = 0;             rotWire(2,2) = 1;
            
  rotM_  = TMatrix( DRSrotM_[rotMNum_-1], TMatrix::kMult, rotWire );

  // irotM_ = TMatrix( TMatrix::kInverted, rotM_ );   
  // do inversion by hand (valid only for dets perp to the beam)
  irotM_ = TMatrix( rotM_ );
  irotM_(0,1) = rotM_(1,0);
  irotM_(1,0) = rotM_(0,1);

  cosTheta_ = rotM_(0,0);
  sinTheta_ = rotM_(1,0);

  return;
}
  
//________________________________________________________
DetectorInfo* DetectorInfo::GetMain( void )
{
  // no subdetectors
  if( sub_.size() <= 1 ) return this;

  DetectorInfo * out = 0;

  // Special case for straws 
  // (jj)
  //  if( sub_[0]->TBName_.substr(0,2) == "ST" ) {
  //   for( unsigned int i=0; i< sub_.size(); i++ ) 
  //   if( sub_[i]->TBName_[7] == 'b' ) {
  //     out = sub_[i];
  //     break;
  //   }
  // }  
  
  // Special case for Muon wall 2
  // jj
  if( sub_[0]->TBName_.substr(0,2) == "MB" ) {
    for( unsigned int i=0; i< sub_.size(); i++ ) 
    if( sub_[i]->TBName_[7] == 'b' ||
        sub_[i]->TBName_[7] == 'r' ) {
      out = sub_[i];
      break;
    }
  }  
  
  // Standard case (use name last digit - and not TBName)
  else {
    char digit ='\0';
    for( unsigned int i=0; i< sub_.size(); i++ ) 
    if( i == 0 ||  sub_[i]->name_[sub_[i]->name_.size()-1] < digit )  { 
      digit = sub_[i]->name_[sub_[i]->name_.size()-1];
      out = sub_[i];
    } 
  }

  return out;
}

//____________________________________________________________________
void DetectorInfo::DUtoCenter( double du )
{
  xcm_ += cosTheta_*du;
  ycm_ += sinTheta_*du;
}

//_______for_pixel_detector___________________________________________
void DetectorInfo::DVtoCenter( double du )
{
  xcm_ += sinTheta_*du;
  ycm_ -= cosTheta_*du;
}
 
//________________________________________________________
int DetectorInfo::UtoWire( double u ) 
{
  // position of the first wire
  double uFirst = wirD_ + cosTheta_*xcm_ + sinTheta_*ycm_;
  if( u < uFirst-0.5*wirP_ || u > uFirst + wirP_*(double(nWir_)-0.5) ) return -1;
  return ( ( u-uFirst )/wirP_ < 0 ) ?
    int( (u - uFirst)/wirP_ - 0.5 ):
    int( (u - uFirst)/wirP_ + 0.5 );
}     
  
//________________________________________________________
string DetectorInfo::Dump( void )
{
  char* line = new char[512];
  string out("");
  sprintf( line,
    " det %4i  %s  %s %2i %4i %8.2f %7.3f %7.3f %8.3f %11.4f %10.5f %10.5f ",
    id_, TBName_.c_str(), name_.c_str(), unit_, type_,
    rdlen_, 
    zsiz_/10, xsiz_/10, ysiz_/10,
    zcm_/10,  xcm_/10,  ycm_/10);
  out += string(line);
  
  sprintf( line,
    "%4i %11.4f %8.3f %5i %10.7f %5.3f %5.1f %6.1f ",
    rotMNum_,
    wirD_/10, ang_, nWir_, wirP_/10,
    eff_, bkg_, tGate_ );
  out += string(line);

  //  if( type_ != 11 ) { delete line; return out; } // jj - obvious bug!
  if( type_ != 11 ) { delete[] line; return out; }
  
  sprintf( line, 
    "%7.5f %7.1f %5.1f %6.2f %8.5f",
    vel_/10, t0_, thRes_, spSli_, tiSli_ );
  out += string(line);
  // delete line;
  delete[] line; //jj - obvious bug!
  return out;

}

#ifndef __CINT__
//_____________________________________________________
Obj3D* DetectorInfo::GetArea3D( void )
{
  if( !area_ ) 
  area_ = (Obj3D*) new Pave3D( TBName_, Point(xcm_, ycm_, zcm_), DRSrotM_[rotMNum_-1], 0.5*xsiz_, 0.5*ysiz_, 0.5*zsiz_);
  return (Obj3D*) area_;
}

//____________________________________________________
Obj3D* DetectorInfo::GetWire3D( int iw )
{
  
  if( !wires_.size() ) 
  for( int i=0; i< nWir_; i++ ) wires_.push_back( 0 );
  if( iw < 0 || iw > nWir_ ) return 0;
  if( wires_[iw] ) return wires_[iw];
      
  // loop over wires
  double wang = ang_+90;
  double sinw = sin(wang*PI/180);
  double cosw = cos(wang*PI/180);
  Point Mid_WRS( wirD_+iw*wirP_, 0 );
  Point Mid_DRS = Mid_WRS.Rotate( ang_ );
  
  // Extrapolate to active area
  Point W1_DRS, W2_DRS;
  
  if( fabs(wang)<45 || fabs(wang)>135 ) {
    W1_DRS = Point( -xsiz_/2, Mid_DRS.y + sinw/cosw*( -xsiz_/2 - Mid_DRS.x ) );
    W2_DRS = Point(  xsiz_/2, Mid_DRS.y + sinw/cosw*(  xsiz_/2 - Mid_DRS.x ) );
  } else {
    W1_DRS = Point( Mid_DRS.x + cosw/sinw*( -ysiz_/2 - Mid_DRS.y ), -ysiz_/2 ); 
    W2_DRS = Point( Mid_DRS.x + cosw/sinw*(  ysiz_/2 - Mid_DRS.y ),  ysiz_/2 ); 
  }
  
  Point W1_MRS = W1_DRS.Rotate( DRSrotM_[rotMNum_-1] )+Point( xcm_, ycm_ );
  Point W2_MRS = W2_DRS.Rotate( DRSrotM_[rotMNum_-1] )+Point( xcm_, ycm_ ); 
  
  // Recalculate wire center position according to limits
  double xcm = 0.5*(W1_MRS.x+W2_MRS.x);
  double ycm = 0.5*(W1_MRS.y+W2_MRS.y);
  double dv  = 0.5*sqrt(pow(W1_MRS.x-W2_MRS.x,2)+pow(W1_MRS.y-W2_MRS.y,2) );
  
  char* buf = new char[64];
  sprintf( buf, "%s_w_%i", TBName_.c_str(), iw );
  wires_[iw] = new Pave3D( string( buf ), Point( xcm, ycm, zcm_), rotM_, 0.5*wirP_, dv, 0.5*zsiz_ ); 
  SafeDelete( buf );
   
  return wires_[iw];
}
#endif 
      
//______________________________________________________
//! Used to dump informations needed for the alignment
ostream &operator << (ostream &out,const DetectorInfo &di) 
{
  printf("%13.4f %10.4f %10.4f %10.4f %10.4f %10.4f %10.4f %10.4f %10.4f %4d %s\n",
    di.zcm_,
    di.res_,
    di.cosTheta_,
    di.sinTheta_,
    di.biasU_,
    di.biasZ_,
    di.biasT_,
    di.biasP_,
    di.biasR_,
    di.id_,
    di.TBName_.c_str()); 
  return out;
}
