// $Id: EffManager.h,v 1.18 2002/09/06 09:00:58 hpereira Exp $
#ifndef EffManager_h
#define EffManager_h

#include <string>
#include <vector>
#include <iostream.h>

//=== all Root Objects ===
#include <TROOT.h>
#include <TObject.h>

#include "HManager.h"
class DetInfo;

//_______________________________________________________________________________
class EffManager: public HManager {
  public:
  EffManager( void );                               //!< default constructor
  EffManager( const HManager h );                   //!< 'Extend' constructor
  EffManager( const char* file, bool batch=false ); //!< usable constructor

  //!< ShortCuts to DetInfo fit methods
  void DrawEff( const char* varSel, const char* selection = "*", TCut cut="1", const char* opt="" );
  void DrawEffUProf(  const char* selection = "*", TCut cut="T_inActive", const char* opt="" );
  inline void SetEffCut( TCut cut ) { cut_ = cut; }
  
  TCut cut_;
  ClassDef(EffManager,0)
  
};

#endif
