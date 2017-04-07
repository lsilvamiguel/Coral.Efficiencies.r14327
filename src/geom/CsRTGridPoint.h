// $Id: CsRTGridPoint.h,v 1.5 2002/03/09 23:29:17 hpereira Exp $

/*!
   \file    CsRTGridPoint.h
   \brief   Simple Grid Point structure for CsRTRelation object
   \author  Hugo Pereira
   \version $Revision: 1.5 $
   \date    $Date: 2002/03/09 23:29:17 $
*/

#ifndef CsRTGridPoint_h
#define CsRTGridPoint_h

class RTGridPoint {
  public:
  double t;      //!< time of RT Grid Point
  double r;      //!< associated distance to wire
  double res;     //!< error on r

  RTGridPoint() : t(0), r(0), res(0) {};                              //!< default creator
  RTGridPoint( const double t, const double r, const double res ):    //!< usable constructor
    t( t ), r( r ), res( res ) {};

  friend std::istream& operator>> (std::istream& in, RTGridPoint &gp) //!< used for DB reading
  { in >> gp.t; in >> gp.r; in >> gp.res; return in; }
};

#endif
