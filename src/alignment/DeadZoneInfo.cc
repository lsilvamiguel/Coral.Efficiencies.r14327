// $Id: DeadZoneInfo.cc,v 1.15 2008/06/05 12:50:17 rgazda Exp $

/*!
   \file    DeadZoneInfo.cc
   \brief   store detector.dat dead zone informations
   \author  Hugo Pereira
   \version $Revision: 1.15 $
   \date    $Date: 2008/06/05 12:50:17 $
*/

#include "DeadZoneInfo.h"
#include "DetectorInfo.h"
#include "Utils.h"
#include "Obj3D.h"
#include "Point.h"

#include <sstream>
#include <fstream>
#include <TMath.h>
char* operator+( std::streampos&, char* );
using namespace std;

//________________________________________________________________
ClassImp(DeadZoneInfo)
DeadZoneInfo::DeadZoneInfo( string detFileName, unsigned int lNumber, string deadInfoLine  ):
  TObject(),
  lineNumber_( lNumber ),
  id_( 0 ),
  TBName_(""), name_(""), unit_( 0 ),
  shape_(0),
  zsiz_(0), xsiz_(0), ysiz_(0),
  xcm_(0),  ycm_(0),  zcm_(0),
  rotMNum_(-1),
  
  rotM_(TMatrix(3,3)),
  irotM_(TMatrix(3,3)),

  det_(0),   
  area_(0)
  
{
  //=== check detFileName
  if( !detFileName.size() ) return;

  //=== first get DRSRotM_ from DetectorFile
  ifstream in( detFileName.c_str(), ios::in );
  if( !in ) {
    cout << "DeadZoneInfo::DeadZoneInfo - FATAL: cannot read file \"" 
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
  
  //=== second, parse detInfoLine to set parameters
  istringstream s( deadInfoLine.c_str() );
  string opt;
  s >> opt;

  if ( opt != "dead" ) {
    cout << "DeadZoneInfo::DeadZoneInfo - FATAL: wrong line type\"" 
         << opt  << "\"." 
         << endl;
    return;
  }
  
  
  s >> id_;
  s >> TBName_;
              
  s >> name_;     // detector name
  s >> unit_;     // detector number in station 
  s >> shape_;
  
  //=== detector size (mm)
  s >> zsiz_;  zsiz_ *= 10;     // convert [cm] to [mm]      
  s >> xsiz_;  xsiz_ *= 10;
  s >> ysiz_;  ysiz_ *= 10;  
  
  //=== detector centre (MRS) (cm)   
  s >> zcm_; zcm_*=10;  // convert [cm] to [mm]       
  s >> xcm_; xcm_*=10; 
  s >> ycm_; ycm_*=10;    
  
  s >> rotMNum_;  // rotation matrix number
  if( rotMNum_ < 0 || rotMNum_ > (int) DRSrotM_.size() ) {
    cout << "DeadFileInfo::DeadFileInfo -FATAL: unexpected matrix number " << rotMNum_ << " for det: " << TBName_ << ".\n";
    exit(0);
  }
  
  rotM_  = TMatrix( DRSrotM_[rotMNum_-1] );
  
  // irotM_ = TMatrix( TMatrix::kInverted, rotM_ );   
  //=== do inversion by hand (valid only for dets perp to the beam)
  irotM_ = TMatrix( rotM_ );
  irotM_(0,1) = rotM_(1,0);
  irotM_(1,0) = rotM_(0,1);

  return;
}  
  
      
//________________________________________________________
string DeadZoneInfo::Dump( void )
{
  char* line = new char[512];
  sprintf( line,
    " dead %4i  %s  %s %2i %2i %7.3f %7.3f %7.3f %10.4f %10.4f %10.4f %2i",
    id_, TBName_.c_str(), name_.c_str(), unit_, shape_,
    zsiz_/10, xsiz_/10, ysiz_/10,
    zcm_/10,  xcm_/10,  ycm_/10,
    rotMNum_);
  string out( line );
  //delete line; // jj - obvious bug!
  delete[] line;
  return out;
}

#ifndef __CINT__
//_____________________________________________________
Obj3D* DeadZoneInfo::GetArea3D( void )
{
  if( area_ ) return area_;
  switch ( shape_ ) {
  
    case 1:  
    area_ = (Obj3D*) new Pave3D( TBName_, Point(xcm_, ycm_, zcm_), DRSrotM_[rotMNum_-1], 0.5*xsiz_, 0.5*ysiz_, 0.5*zsiz_);
    break;
        
    case 5: {
      TMatrix m(3,3); m = TMatrix( TMatrix::kUnit, m );
      area_ = (Obj3D*) new Circle3D( TBName_, Point(xcm_, ycm_, zcm_), m, 0.5*xsiz_, 0.5*zsiz_);
    }
    break;
    
    default:
    cout << "DeadZoneInfo::GetArea - \"" << TBName_ << "\" unrecognized shape " << shape_ << ".\n";
    area_ = 0;
    break;
  }
  
  return area_;
  
}
#endif 
