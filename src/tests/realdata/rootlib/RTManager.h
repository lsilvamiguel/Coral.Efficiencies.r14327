// $Id: RTManager.h,v 1.20 2002/09/06 09:00:59 hpereira Exp $

#ifndef RTManager_h
#define RTManager_h

#include "HManager.h"

#include <string>

#include <TROOT.h>
#include <TObject.h>
#include <vector>

class DetInfo;
class RTGridPoint;
class RTInfo;
class TH2S;

class RTManager: public HManager {
  public:
  RTManager( void );             //!< default constructor
  RTManager( const HManager h ); //!< 'Extend' constructor
  RTManager( const char* file, bool batch = false ); //!< usable constructor
  
  //!< ShortCuts to DetInfo fit methods
  void DumpRTInfo       ( const char* selection = "*" );
  void DrawT0           ( const char* selection = "*", TCut cut = "T_fnd" );
  void FitT0            ( const char* selection, TCut cut, double tMin, double tMax );
  void FitRTStraights   ( const char* selection = "*", TCut cut = "T_fnd", double rMin = 1, double rMax = -1 );
  void FitRTRelations   ( const char* selection = "*", TCut cut = "T_fnd", const char* detFile = "" );
  void FitRTGrids       ( const char* selection = "*", TCut cut = "T_fnd", unsigned int nBin = 0 );
  void CheckResolutions ( const char* selection = "*", TCut cut = "T_fnd", unsigned int nEvt = 1000000000 );

  void DrawSingleRT( const char* selection = "*", TCut cut = "T_fnd" );
  void FitSingleRTStraight( double rMin = 1, double rMax = -1  );
  void FitSingleRTGrid( unsigned int nBins = 0 );

  //!< Shortcuts to DetInfo I/O methods
  void WriteRTParsToDB ( const char* selection = "*", const bool useMyT0 = false, const double myT0 = 0 );
  void WriteRTGridsToDB( const char* selection = "*", const bool useMyT0 = false, const double myT0 = 0 );
  void SetValidity        ( const char* selection = "*", const char* start = "YYYY-MM-DD-hh:mm:ss", const char* stop = "YYYY-MM-DD-hh:mm:ss" );
 
  private:  
  void _WriteT0ToLog( DetInfo* detI );
  TH2S* H_single_RT_;    
  RTInfo *single_rtInfo_;
  
  ClassDef(RTManager,0)
  
};

#endif
