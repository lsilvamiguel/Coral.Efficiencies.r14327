// $Id: RTInfo.h,v 1.2 2007/06/01 09:13:52 conrad Exp $

/*!
   \file    RTInfo.h
   \brief   to store detector associated RT relation informations
   \author  Hugo Pereira
   \version $Revision: 1.2 $
   \date    $Date: 2007/06/01 09:13:52 $
*/

#ifndef RTInfo_h
#define RTInfo_h

#include <string>
#include <vector>
#include <iostream>

#include <TROOT.h>
#include <TObject.h>

/*!
   \class   RTGridPoint
   \brief   to store r,t points of RT relations in case of grid parametrisation
*/

class RTGridPoint {
  public:
  double t;      //!< time of RT Grid Point
  double r;      //!< associated distance to wire
  double res;     //!< error on r

  RTGridPoint() : t(0), r(0), res(0) {};                              //!< default creator
  RTGridPoint( const double t, const double r, const double res ):    //!< usable constructor
    t( t ), r( r ), res( res ) {};

  friend std::istream& operator>> (std::istream& in, RTGridPoint &gp)           //!< used for DB reading
  { in >> gp.t; in >> gp.r; in >> gp.res; return in; }
};

/*!
   \class   RTInfo
   \brief   to store detector associated RT relation informations
*/

class RTInfo: public TObject {
  public:
  RTInfo( std::string TBName = "" );

  inline bool HasTBName( void ) { return TBName_.size(); }
  inline bool HasT0Par( void )   { return !(t0Par_.empty()); }
  inline bool HasRTPar( void )   { return !(rtPar_.empty()); }
  inline bool HasRTGrid( void )  { return !(rtGrid_.empty()); }
  void DumpT0Par( void ); 
  void DumpRTPar( void ); 
  void DumpRTGrid( void );   
  bool SetValidity( const char* start = "YYYY-MM-DD-hh:mm:ss", const char* stop = "YYYY-MM-DD-hh:mm:ss" );
  bool WriteRTParToDB(  const bool useMyT0 = false, const double myT0 = 0 );
  bool WriteRTGridToDB( const bool useMyT0 = false, const double myT0 = 0 );
  
  bool _CheckValidityFormat( const char* stamp );
  void _AddRTGridPoint( const double t, const double r, const double res  );

  // copied from coral/src/geom/CsRTRelation
  double GetRfromT( const double t, bool &error );  

  std::string TBName_;
  std::vector< double > t0Par_;
  std::vector< double > rtPar_;
  std::vector< RTGridPoint > rtGrid_;
  std::string valid_start_;
  std::string valid_stop_;

  ClassDef(RTInfo,0)

};
#endif
  
