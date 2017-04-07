// $Id: DeadZoneInfo.cc,v 1.1 2002/07/08 20:34:22 hpereira Exp $

#include "DeadZoneInfo.h"
#include "Utils.h"

#include "Utils.h"
#include <strstream.h>
#include <fstream.h>
#include <TMath.h>

//________________________________________________________________
ClassImp(DeadZoneInfo)
DeadZoneInfo::DeadZoneInfo( string detFileName, unsigned int lNumber, string deadInfoLine  ):
  TObject(),
  Comment_( "" ),
  lineNumber_( lNumber ),
  id_( 0 ),
  TBName_(""), name_(""), unit_( 0 ),
  shape_(0),
  xsiz_(0), ysiz_(0), zsiz_(0),
  xcm_(0), ycm_(0), zcm_(0),
  rotMNum_(-1),
  rotM_(TMatrix(3,3)),
  irotM_(TMatrix(3,3))
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
  istrstream s( deadInfoLine.c_str() );
  string opt;
  s >> opt;
  if ( opt != "dead" ) {
    Comment_ = opt;
    s >> opt;
  }
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
    cout << "DeadFileInfo::DeadFileInfo -FATAL: unexpected matrix number " << rotMNum_ << ".\n";
    exit(0);
  }
  
  rotM_  = TMatrix( DRSrotM_[rotMNum_-1] );
  irotM_ = TMatrix( TMatrix::kInverted, rotM_ );   

  return;
}  
  
      
//________________________________________________________
string DeadZoneInfo::Dump( void )
{
  char* line = new char[512];
  sprintf( line,
    "%s dead %4i  %s  %s %2i %2i %7.3f %7.3f %7.3f %10.4f %10.4f %10.4f %2i",
    Comment_.c_str(), id_, TBName_.c_str(), name_.c_str(), unit_, shape_,
    zsiz_/10, xsiz_/10, ysiz_/10,
    zcm_/10,  xcm_/10,  ycm_/10,
    rotMNum_);
  string out( line );
  delete line;
  return out;
}

