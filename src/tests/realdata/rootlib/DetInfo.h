// $Id: DetInfo.h,v 1.38 2002/11/06 14:02:13 hpereira Exp $

#ifndef DetInfo_h
#define DetInfo_h

#include <string>
#include <iostream.h>
#include <strstream.h>
#include <fstream.h>
#include <unistd.h>
#include <TROOT.h>
#include <TObject.h>
#include <TChain.h>
#include <TTree.h>
#include <TCut.h>
#include <TH2.h>
#include <TH1.h>
#include <TMatrix.h>
#include <TProfile.h>
#include <TText.h>
#include <TGraphErrors.h>

#include <vector>

class DetFileInfo;
class RTInfo;
class Utils;
class TH2Fitter;
class DetInfo: public TObject {
  public:
  DetInfo( string TBName="", string detFile="" ); 
  DetInfo( int ID, string detFile );     
  
  void AddToChain( const char* treeName, const char* fileName );
  bool SetDetFile( const char* detFile );


  inline bool HasTBName( void )  { return bool(TBName_.size()); }

  TF1* FitTvsV( TCut cut="T_fnd", double p0 = 0, double p1 = 0 );
  bool FitT0 (        TCut cut, double tMin, double tMax );
  bool FitRTStraight( TCut cut="T_fnd", double rMin = 1, double rMax = -1 );
  bool FitRTRelation( TCut cut="T_fnd" );
  bool FitRTGrid (    TCut cut="T_fnd", unsigned int nBins = 0 );


  //=== Residuals Managment
  bool FitDU(         TCut cut="T_fnd" );
  bool FitDU_MWPC(    TCut cut="T_fnd" );
  bool FitDUvsU(      TCut cut="T_fnd" );
  bool FitDUvsV(      TCut cut="T_fnd" );
  bool DrawDU_LR(     TCut cut="T_fnd", double OffsetN = 0, double OffsetP = 0 );
  bool DrawWRSProfile( TCut cut="" );

  string TBName_;

  //=== TreeForAll
  TChain* T_eff_;   //!< All histograms (except GUIs) are switched off when tree is found.
 
  //=== Geometry informations ===
  bool hasDetFile_;
  string detFile_;
  vector< DetFileInfo* > detFList_;
  DetFileInfo* detFI_;
  
  //=== residual from fit
  double resid_;   //=== sigma value for T_duMin in the detector
  
  //=== RTRelation Informations
  RTInfo* rtInfo_;

  //=== Various Switches
  bool OneGaus_;          //!< use only on gauss to fit resolution if true
  bool UseCorrectedVals_; //!< use propagation time corrected du if true
  bool CorrectFromAngle_; //!< correct du from track angle effects (du_cor=uDet*cos(atan(T_tuLx))*udet-uLx) if true
  
  private:
  void _AddDetFileInfo( DetFileInfo *dfi );
  bool _ReadDetFile( void );  
  bool _GetTBNameFromDetFile( int ID );  
    
  ClassDef(DetInfo,0)
   
};
#endif
