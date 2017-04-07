// $Id: Point.cc,v 1.3 2008/04/10 14:36:39 rgazda Exp $
/*!
   \file    Point.cc
   \brief   x,y point object
   \author  Hugo Pereira
   \version $Revision: 1.3 $
   \date    $Date: 2008/04/10 14:36:39 $
*/
#include "TMath.h"
#include "Point.h"
#include "Defs.h"

#include <iostream>

using namespace std;

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
ostream &operator << (ostream &o,const Point &p)
{
  char* text = new char[200];

  snprintf(text, 200, "(%10.4f,%10.4f,%10.4f)", p.x, p.y, p.z );
  o<<text;
  SafeDelete(text);

  return o;
}
