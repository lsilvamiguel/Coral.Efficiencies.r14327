// $Id: Macro.h,v 1.15 2002/10/14 21:02:16 hpereira Exp $
#ifndef Macro_h
#define Macro_h

#include <TROOT.h>
#include <TObject.h>
#include <TH2.h>
#include <TTree.h>
#include <TCanvas.h>
#include <fstream.h>
#include <strstream.h>
#include <vector>
#include <string>

#include "HManager.h"

class DetInfo;
class AlManager;

//___________________________
class Macro: public TObject {
  public:
  Macro( ): TObject() { };
    
  static void PropTime( const char* fileSelection, 
    const char* cut = "abs(T_vLx)<300", 
    const char* outFile = "PropTime.out" );
  
  static void DUFit( const char* fileSelection,
    const char* cut = "",
    const char* tag="",
    bool oneGaus = false,
    bool UseCorrectedVals = false );
  
  static void RTGrid( const char* fileSelection,
    const char* cut = "",
    const char* tag="",
    unsigned int n=0 );
  
  static void EffvsXY( const char* fileSelection,
    const char* cut = "",
    const char* tag="");
  
  static void EffvsUV( const char* fileSelection,
    const char* cut = "",
    const char* tag="");
    
  static void ProfileToRates( TH2* h, 
    unsigned int nEvt = 1, 
    double tMin=0, double tMax=0 );
  
  private:
  static vector< string > _GetFiles( const char* fileSelection );
  
  ClassDef(Macro,0)

};

#endif
