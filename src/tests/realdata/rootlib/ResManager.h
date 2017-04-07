// $Id: ResManager.h,v 1.20 2002/11/06 14:03:40 hpereira Exp $
#ifndef ResManager_h
#define ResManager_h

#include <string>
#include <vector>
#include <TROOT.h>
#include <TObject.h>

#include "HManager.h"

class DetInfo;
class Fit;

//_______________________________________________________________________________
class ResManager: public HManager {

  public:
  ResManager( void );                                          //!< default constructor
  ResManager( const HManager h );                              //!< 'Extend' constructor
  ResManager( const char* fileSelection, bool batch = false ); //!< usable constructor
  
  //!< ShortCuts to DetInfo fit methods
  void FitDU(           const char* selection = "*", TCut cut = "T_fnd" );
  void FitDU_MWPC(      const char* selection = "*", TCut cut = "T_fnd" );
  void FitDUvsU(        const char* selection = "*", TCut cut = "T_fnd" );
  void FitDUvsV(        const char* selection = "*", TCut cut = "T_fnd" );
  void DrawDU_LR(       const char* selection = "*", TCut cut = "T_fnd", 
    double OffsetN = 0, 
    double OffsetP = 0 );
  void DrawWRSProfiles( const char* selection = "*", TCut cut = "" );
  bool FitResolutions( const char* selection, TCut cut = "T_fnd", const char* detFileNamed = 0 );
  
  protected:
  double _ResCalc( unsigned int it, vector< double > res );
  friend class Fit;
   
  ClassDef(ResManager,0)
  
};

#endif
