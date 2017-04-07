// $Id: HManager.h,v 1.21 2002/11/06 14:04:54 hpereira Exp $

#ifndef HManager_h
#define HManager_h

#include <string>
#include <iostream.h>

//=== all Root Objects ===
#include <TROOT.h>
#include <TObject.h>
#include <TCut.h>

//=== needed to define list vectors etc ===
#ifdef COMPASS_USE_OSPACE_STD
# include <ospace/std/vector>
#else
# include <vector>
#endif 

class DetInfo;
class TFile;
class TCanvas;

class HManager: public TObject {
  public:
  HManager( void );                                          //!< default constructor
  HManager( const char* fileSelection, bool batch = false ); //!< usable constructor
  void Draw( const char* varSel, const char* selection = "*", TCut cut="1", const char* opt="" ); 
  void Fit( const char* varSel, const char* selection = "*", TCut cut="1", const char* func="gaus", const char* opt="" );
  void DumpTBNames( const char* selection = "*" );
  inline HManager* GetPointer() { return ptr_; }
  static HManager* ptr_;

  inline TCanvas* GetTCanvas( void ) { return tcv_; } //!< To have access to TCanvas throught CINT
  DetInfo* GetDetInfo( const char* TBName = "" );
  void SetDetFile( const char* file );
  string fileSelection_;
  bool UseCorrectedVals_;     //<! To be passed to all DetInfos
  TCanvas* tcv_;              //<! Main canvas
  bool OneGaus_;              //<! To be passed to all DetInfos
  bool Batch_;                //<! No TCanvas drawn
  void _LoadTBNames( TFile* tf );

  vector< DetInfo* > selDetI_; //!< current list of selected detector information objects
  
  //=== protected functions
  
  protected:
  vector< DetInfo* > GetSelection( const char* selection );
  vector< DetInfo* > detInf_;  //!< detector information objects
  vector< string > _GetFiles( const char* fileSelection );
  void _SortDetInfos( void );
  struct sortDetInfos_;
  void _DivideCanvas( const unsigned int nPads ); 
  
  int color_; //! color index for 'same' diagrams 
  
  ClassDef(HManager,0)
  
};


#endif

