// $Id: RTInfo.h,v 1.10 2002/08/17 20:54:26 hpereira Exp $

#ifndef RTInfo_h
#define RTInfo_h

#include <string>
#include <vector.h>
#include <TROOT.h>
#include <TObject.h>

class RTGridPoint;
class RTInfo: public TObject {
  public:
  RTInfo( string TBName = "" );

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

  //=== copied from coral/src/geom/CsRTRelation
  double GetRfromT( const double t, bool &error );  

  string TBName_;
  vector< double > t0Par_;
  vector< double > rtPar_;
  vector< RTGridPoint > rtGrid_;
  string valid_start_;
  string valid_stop_;

  ClassDef(RTInfo,0)

};
#endif
  
