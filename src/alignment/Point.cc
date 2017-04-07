// $Id: Point.cc,v 1.5 2008/02/20 13:30:37 rgazda Exp $
/*!
   \file    Point.cc
   \brief   x,y point object
   \author  Hugo Pereira
   \version $Revision: 1.5 $
   \date    $Date: 2008/02/20 13:30:37 $
*/
#include "Point.h"
#include "TMath.h"
#include "Defs.h"
#include <iostream>


//________________________________________________________________
const Point Point::Orig = Point();
Point Point::Rotate( TMatrix M )
{ return Point( x*M(0,0)+y*M(0,1), x*M(1,0)+y*M(1,1), z ); }

//________________________________________________________________
Point Point::Rotate( double ang_deg )
{ 
  double cosang = cos( ang_deg*PI/180 );
  double sinang = sin( ang_deg*PI/180 );
  return Point( x*cosang-y*sinang, x*sinang+y*cosang, z ); 
}

//________________________________________________________________
//! to dump object to screen
std::ostream &operator << (std::ostream &o,const Point &p)
{
  printf("(%10.4f,%10.4f,%10.4f)", p.x, p.y, p.z );
  return o;
}
