// $Id: Point.h,v 1.8 2008/03/05 20:36:26 ybedfer Exp $
#ifndef Point_h
#define Point_h
/*!
   \file    Point.h
   \brief   x,y point object
   \author  Hugo Pereira
   \version $Revision: 1.8 $
   \date    $Date: 2008/03/05 20:36:26 $
*/
#include <stdio.h>
#include <math.h>
#include <TROOT.h>
#include <TMatrix.h>

/*!
   \class   Point
   \brief   x,y point object
*/
// class ostream;
class Point {
  public:
  Point( double x=0, double y=0, double z=0 ): x(x), y(y), z(z) {}  //!< constructor
  Point Rotate( TMatrix M );                      //!< rotate point coordinates perp to z, using rotation matrix M
  Point Rotate( double ang_deg );                 //!< rotate point coordinates perp to z, using angle in deg
  Point operator+ (const Point &p) { return Point( x+p.x, y+p.y); }  //!< adds point x, y coordinates
  bool  operator< (const Point &p) { return z < p.z; }               //!< compares points according to z
  double x;  //!< x coordinate
  double y;  //!< y coordinate
  double z;  //!< z coordinate 

  static const Point Orig; //!< (0,0,0)

  //! Used to dump informations needed for the alignment
  friend std::ostream &operator << (std::ostream &o,const Point &p);
};

#endif
