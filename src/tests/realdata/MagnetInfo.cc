// $Id: MagnetInfo.cc,v 1.2 2007/06/01 09:13:52 conrad Exp $


/*!
   \file    MagnetInfo.cc
   \brief   store detector.dat magnet informations
   \author  Hugo Pereira
   \version $Revision: 1.2 $
   \date    $Date: 2007/06/01 09:13:52 $
*/

#include "MagnetInfo.h"
#include "Obj3D.h"
#include "Point.h"

#include <sstream>
#include <fstream>

#include <TMatrix.h>

using namespace std;

//______________________________________________________________________________
const MagnetSize MagnetInfo::SM1Size_( 1150, 1100, 2100, 2000, 850 );
const MagnetSize MagnetInfo::SM2Size_( 1005, 495,  3025, 2000, 1500 );

//________________________________________________________________
ClassImp(MagnetInfo)
MagnetInfo::MagnetInfo(
  string detFileName, 
  unsigned int lNumber, 
  string detInfoLine  ):
  TObject(),
  lineNumber_( lNumber ),
  id_( 0 ), 
  xcm_( 0 ),
  ycm_( 0 ),
  zcm_( 0 ),
  rot_( 0 ),
  scale_( 0 ),
  flag1_( 0 ),
  flag2_( 0 ),
  current_( 0 ),
  gotCurrent_( false ),
  area_( 0 )
{
  //=== second, parse detInfoLine to set parameters
  istringstream s( detInfoLine.c_str() );
  string opt;
  s >> opt;

  if ( opt != "mag" ) {
    cout << "MagnetInfo::MagnetInfo - FATAL: wrong line type\"" 
         << opt  << "\"." 
         << endl;
    return;
  }
  
  s >> id_;
  s >> zcm_; zcm_ *= 10;  // convert [cm] to [mm]
  s >> xcm_; xcm_ *= 10;
  s >> ycm_; ycm_ *= 10;
  
  s >> rot_;
  s >> scale_;
  s >> flag1_;
  s >> flag2_;
  
  s >> current_ ;
  if( s.good() || s.eof() ) gotCurrent_ = true;
  return;
}

#ifndef __CINT__
//_____________________________________________________
Obj3D* MagnetInfo::GetArea3D( void )
{
  if( !area_ ) {
    TMatrix I(3,3); I = TMatrix( TMatrix::kUnit, I );
    switch( id_ ) {
     
      case 4: // SM1
      area_ = (Obj3D*) new Tunnel3D( "SM1", Point(xcm_, ycm_, zcm_), I, 
        SM1Size_.xsizIn_, SM1Size_.ysizIn_, 
        SM1Size_.xsizOut_, SM1Size_.ysizOut_, 
        SM1Size_.zsiz_);
      break;
      
      case 5: // SM2
      area_ = (Obj3D*) new Tunnel3D( "SM2", Point(xcm_, ycm_, zcm_), I, 
        SM2Size_.xsizIn_, SM2Size_.ysizIn_, 
        SM2Size_.xsizOut_, SM2Size_.ysizOut_, 
        SM2Size_.zsiz_);
      
      default: break;
    }
  }
  return (Obj3D*) area_;
}
#endif


  
     
