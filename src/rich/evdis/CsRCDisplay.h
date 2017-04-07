// $Id:

/*!
   \file    CsDisplay.h
   \brief   CORAL Event Display Package.
   \version $Revision: 1.2 $
   \author  Take-Aki TOEDA
   \date    $Date: 2010/06/18 10:44:21 $
*/

#ifndef CsRCDisplay_h
#define CsRCDisplay_h

#include <iostream>

#include "coral_config.h"
#include "TObject.h"
#include "TGFrame.h"
#include "TGButton.h"
#include "TLine.h"
#include "TRint.h"
class TPaveText;
class TH2;
class TStyle;
class TFile;
class TCanvas;
//class TButton;
class TControlBar;
class TText;
class TLatex;
class TPaveLabel;
class TView;
class TList;
class TGeometry;
class TGFrame;
class TGLabel;
class TGTextView;
class TGTab;
class TGClient;
//class TGButton;
class TGMenu;
class TRootEmbeddedCanvas;
class TGMsgBox;
class TGTextEntry;
class TBRIK;
class TNode;
class TPostScript;
class TH2F;
class TGTextBuffer;
class TGTextFrame;
class TMaterial;
//#include <TGTextWindow.h>

#include "CsRCEDisplay.h"
#include "CsRCGTrack.h"
#include "CsRCGRingTrack.h"
#include "CsRCGHit.h"
#include "CsRCGCathode.h"
#include "CsRCGRing.h"

enum RCommandsId {
RF_QUIT,
RF_NEXT,
RF_PREVIUS,

RF_NEXT_RING,
RF_PREVIUS_RING,

RF_ZOOM_IN,
RF_ZOOM_OUT,

RF_SCALE_INC,
RF_SCALE_DEC,

RF_EXPLOR_UP,
RF_EXPLOR_DOWN,
RF_EXPLOR_LEFT,
RF_EXPLOR_RIGHT,

RF_RESET_VIEW,
RF_RELOAD_AF,

RF_TOP_VIEW,
RF_SIDE_VIEW,
RF_FRONT_VIEW,

RF_GRAPHIC_OUT,

RF_COLOR_BACK,

RF_DET_ON_OFF,
RF_DIGIT_ON_OFF,
RF_HIT_ON_OFF,
RF_CLUSTER_ON_OFF,
RF_TRACK_ON_OFF,

RF_X3D,

}; 



/*! \class CsRCDisplay
    \brief   CORAL Event Display Package.
    Event display class 
*/
enum RDet {RWith,RWithout};
enum RInit_mode {REmbedded,RStandalone};

class CsRCDisplay : public TGMainFrame
{
 public:
  CsRCDisplay(RInit_mode);
  virtual ~CsRCDisplay();

  virtual Bool_t ProcessMessage(Long_t msg, Long_t parm1, Long_t parm2);

  void View();
  void SetRange(float32,float32,float32,float32,float32,float32);
  void SetView(float,float,float);
  void BuildControlPanel(RInit_mode);
  void BackColor();
  //  void BuildAxis();
  //  void DrawAxis();
  void DrawHits();
  void DrawClusters();
  void DrawCathodes();
  void DrawRings();
  void DrawRingsIpo(Int_t);
  void DrawInfo();
  void ReloadAF();
  void RingOn();
  void RingOff();
  void HitOn();
  void HitOff();
  void ClusterOn();
  void ClusterOff();
  void TrackOn();
  void TrackOff();
  void DrawTracks();

  void TopView();
  void SideView();
  void FrontView();

  bool GetContinueDisplay();
  
  //EButtonState GetObjectState(RDet);

  void NextEvent();
  void PreviusEvent();
  void NextRing();
  void UpdateDisplay();
  void PreviusRing();
  void OpenEvent();
  void GetRichDisplay();
  void PrintRichDisplay();
  void MakeGraphicFile();
  void ClearParticleInfo();
  void ReadEvent();
  void SaveViewDirection();
  void SetHit(CsRCGHit *);
  void SetCluster(CsRCGTrack *);
  void SetTrack(CsRCGTrack *);
  void SetRingTrack(CsRCGRingTrack *);
  void SetCathode(CsRCGCathode *);
  void SetRing(CsRCGRing *);
  void SetRingIpo(CsRCGRing *[3]);

  void SetRunNumber(uint64 runNumber)
    { RunNumber = runNumber;}
  void SetEventNumberInRun(uint64 eventNumberInRun)
    { EventNumberInRun = eventNumberInRun;}
  void SetEventNumberInBurst(uint64 eventNumberInBurst)
    { EventNumberInBurst = eventNumberInBurst;}

  void ClearViewPad();
  void ExploreUp();
  void ExploreDown();
  void ExploreLeft();
  void ExploreRight();
  void Redraw();
  void ResetView();
  void ScaleInc();
  void ScaleDec();
  void X3D();
  void ZoomIn();
  void ZoomOut();
  void Quit();

 private:

  TGMainFrame *fMF;
  TGCompositeFrame *fF0,*fF0h,*fF1,*fF2,*fFPIC;
  TGCompositeFrame *Tf0,*Tf1,*Tf2,*Tf3,*Tf4;
  TGCompositeFrame *TIf0,*TIf1;
  TGCompositeFrame *Tif0,*Tif1;
  TGCompositeFrame *fTab3_1,*fTab3_2,*fTab3_3,*fTab3_4,*fTab3_5;
  TGCompositeFrame *fTab3_6,*fTab3_7,*fTab3_8,*fTab3_9,*fTab3_10;
  TGCompositeFrame *fTab3_11,*fTab3_12,*fTab3_13;
  TGGroupFrame *fT1,*fT2,*fT21,*fT22,*fT23,*fT24,*fT3,*fT31,*fT32,*fT33,*fT4;
  TGGroupFrame *fT1_1,*fT1_2,*fT1_3,*fT1_4;
  TGGroupFrame *fT3_3_1,*fT3_3_2,*fT3_3_3,*fT3_3_4,*fT3_3_5;        ;
  TGGroupFrame *fT3_3_6,*fT3_3_7,*fT3_3_8,*fT3_3_9,*fT3_3_10;
  TGGroupFrame *fT3_3_11,*fT3_3_12,*fT3_3_13;

  static const int NtheHist=4;
  TH1F *theHist[NtheHist];
  TRootEmbeddedCanvas *fCanvasA,*fCanvasB,*fCanvasC,*fCanvasD,*fCanvasE;
  TCanvas *m_canvas,*info_canvas,*phi_canvas,*ring_canvas;
  TPad    *m_pad;
  TPad    *m_buttons;
  TControlBar *controller;
  TView *view,*phiview;
  TGGroupFrame *fButtonFrame,*fExploreFrame;
  TGTextButton *fExitButton,*fNextEventButton,*fPreviusEventButton;
  TGTextButton *fNextRingButton,*fPreviusRingButton;
  TGTextButton *fZoomInButton,*fZoomOutButton;
  TGTextButton *fScaleIncButton,*fScaleDecButton;
  TGTextButton *fResetViewButton,*fDetOnOff;
  TGTextButton *fRegionUpButton,*fRegionDownButton;
  TGTextButton *fRegionLeftButton,*fRegionRightButton;
  TGTextButton *fReloadAF,*fTopView,*fSideView,*fFrontView;
  TGTextButton *fGraphicOutButton,*fBackColor,*fX3D,*fDigitOnOff,*fHitOnOff,*fClusterOnOff;
  TGTextButton *fTrackOnOff;
  TGTextEntry  *fGraFileName;
  TGTextBuffer *fGraFileBuf;
  TGTab *fTab,*fITab,*fiTab,*fTab3;
  TGTextFrame *fTextFrame;
  TGTextView *fTextView;
  TList *fAxesItems;

  CsRCEDisplay *edisplay,*hdisplay,*phidisplay,*rdisplay,*idisplay;

  bool ContinueDisplay;

  //Flags to draw or not
  bool drawHits;
  bool drawClusters;
  bool drawTracks;
  bool drawAxis;
  bool drawDigits;
  bool drawRings;

  bool drawRICH1;
  bool drawRICH2;

  TGPictureButton *fPicBut; 
 
  TMaterial *mat;
  TBRIK *HallShape;

  TList *listTrack;
  TList *listRingTrack;
  TList *listHit;
  TList *listCluster;
  TList *listCathode;
  TList *listRing;
  TList *listRingIpo;

  TText  *text; 
  TLatex *tlatex; 

  float32 Phi,Theta,Psi;
  float32 min[3],max[3];
  float32 Longitude,Latitude;

  static const int32 cPion   = 2;
  static const int32 cKaon   = 3;
  static const int32 cProton = 4;

  uint64 RunNumber,EventNumberInRun,EventNumberInBurst;
  Int_t  nLoopRing;

  uint32 nParticles,nRings,nRingPhot;
  bool nFlag1,nFlag2,nFlag3;
  float ringMom;
  float ringThe,ringTheReco,ringTheRFit;
  double   raPart[3];
  double  thePart[3];
  double  chiPart[3];
  double likePart[3];
  uint32 nMuonPlus,nMuonMinus,nGamma,nProton,nElectron;

  enum RView {RTop,RSide,RFront} ViewDirection;
  enum RColor {RBlack = 1,RWhite = 0} BGcolor;
  
 public:
    ClassDef(CsRCDisplay,0)
};

extern TRint *theAppR;

#endif
