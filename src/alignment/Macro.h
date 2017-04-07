// $Id: Macro.h,v 1.12 2006/06/16 15:21:47 conrad Exp $
#ifndef Macro_h
#define Macro_h

/*!
   \file    Macro.h
   \brief   some root macro to make predefined plots, fits, etc.
   \author  Hugo Pereira
   \version $Revision: 1.12 $
   \date    $Date: 2006/06/16 15:21:47 $
*/

#include <TROOT.h>
#include <TObject.h>
#include <TCanvas.h>
#include <string>
#include <vector>
#include "Utils.h"

class TH2;
class TH1;
class DetFileManager;

/*!
   \class Macro
   \brief some root macro to make predefined plots, fits, etc.
*/
class Macro: public TObject {
  public:

  //! add U offsets to gems when going from alignment run to physics run to take lorentz angle into accound
  static void AddGemOffset( DetFileManager &df );

  //! add U offsets to micromegas when going from alignment run to physics run to take lorentz angle into accound
  static void AddMMOffset( DetFileManager &df );

  //! put correct value for MB tiSli, t0 and gate
  static void SetMBTiSli( const char* fileselection );

  //! Change GM02 Order: U (5185.592) < V (5185.642) < Y (5209.222) < X (5209.272)
  static void ReorderGM02( const char* fileselection );
  
  //! define Gem Dead zone size as a 5cm disk
  static void SetGemDeadZoneSize( const char* fileselection );
  
  //! define Gem Dead zone size as a 5cm disk
  static void SetDCDeadZoneSize( const char* fileselection );
  
  //! Yann Bedfer macro to draw Chi2 for tracks bridged accross diff zones
  static TH2* Chi2( const char* fileselection ){ return Chi2( Utils::GetFiles( fileselection ) ); }

  //! make bit pattern plots using Yann Bedfer TrafDic histograms
  static TCanvas* BitPattern( const char* fileselection, TCanvas *cv = 0, int zone_mask = 0x1f )
  { return BitPattern( Utils::GetFiles( fileselection ), cv, zone_mask ); }
  
  //! Yann Bedfer macro to have tracks, beam and vertex statistics
  static void Eval( const char* filename, int trigger_mask = 0x7 );
  
  //! retrieve FI01/02/15 residual histograms obtained using halo muon with trafdic
  static TCanvas* HaloFI( const char* fileselection ); 
  
  //! retrieve SI residual histograms obtained using halo muon with trafdic
  static TCanvas* HaloSI( const char* fileselection ); 
  
  #ifndef __CINT__ 
  static TH2* Chi2( std::vector< std::string > files );
  static TCanvas* BitPattern( std::vector<std::string>detfiles, TCanvas *cv = 0, int zone_mask = 0x1f );
  #endif
  
  ClassDef(Macro,0)              //!< Root Macro
};

#endif
