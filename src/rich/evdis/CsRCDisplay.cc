// $Id:

/*!
  \file    CsDisplay.cc
  \brief   CORAL Event Display Package.
  \version $Revision: 1.8 $
  \author  Take-Aki TOEDA
  \date    $Date: 2010/10/06 08:28:45 $
*/
#include "coral_config.h"
#include "CsEvent.h"
// --- rich includes ---
#include "../CsRichOne.h"
#include "../CsRCEventPads.h"
#include "../CsRCPad.h"
#include "../CsRCCathode.h"
#include "../CsRCDetectors.h"
#include "../CsRCEventClusters.h"
#include "../CsRCCluster.h"
#include "../CsRCEventParticles.h"
#include "../CsRCPartPhotons.h"
#include "../CsRCPhoton.h"
#include "../CsRCEventRings.h"
#include "../CsRCRing.h"
#include "../CsRCEventAnalysis.h"
#include "../CsRCExeKeys.h"
#include "../CsRCRecConst.h"
#include "../CsRCHistos.h"
#include "../CsRCUty.h"
// --- end of rich includes --- 
#include "CsGeom.h"
#include "CsField.h"
#include "CsRCDisplay.h"
 
#include "RVersion.h"
#include "TROOT.h"
#include "TRint.h"
#include "TPaveText.h"
#include "TFrame.h"
#include "TH2.h"
#include "TStyle.h"
#include "TFile.h"
#include "TCanvas.h"
#include "TButton.h"
#include "TControlBar.h"
#include "TText.h"
#include "TLatex.h"
#include "TPaveLabel.h"
#include "TView.h"
#include "TList.h"
#include "TGeometry.h"
#include "TGaxis.h"
#include "TGLabel.h"
#include "TGTextView.h"
#include "TGTab.h"
#include "TGClient.h"
#include "TGMenu.h"
#include "TRootEmbeddedCanvas.h"
#include "TGMsgBox.h"
#include "TGTextEntry.h"
#include "TBRIK.h"
#include "TNode.h"
#include "TPostScript.h"
#include "TPolyMarker3D.h"

#include <sstream>
#include <fstream>
#include <iomanip>

using namespace std;
using namespace CLHEP;

extern TRint *theAppR;

ClassImp(CsRCDisplay)

  TNode *hall;
Float_t ZoomFactor = 0.7;
Float_t ExplorFactor = 0.1;
Float_t d_min[3] = {-900,-900,4000},d_max[3]={900,900,9000};
Float_t d_Phi =   0,d_Theta = 0,d_Psi = -90;

int32 cDisplay = 1;
int32 cCanvas  = 1;

CsEvent *event;
RInit_mode Mode;

CsRCDisplay::CsRCDisplay(RInit_mode mode) : TGMainFrame(gClient->GetRoot(),1,1){
  Mode = mode;

  if(mode == RStandalone) event = CsEvent::Instance();

//--- from "CsRCRecConst"
  CsRCRecConst *cons = CsRCRecConst::Instance();

  int mcaScan = cons->mcaScan();
  float xlScan = cons->xlScan();
  float xuScan = cons->xuScan();

  theHist[0] = new TH1F("the_ring0","",mcaScan,xlScan,xuScan);
  theHist[2] = new TH1F("the_ring0all","",mcaScan,xlScan,xuScan);
  theHist[1] = new TH1F("the_ring1","",mcaScan,xlScan,xuScan);
  theHist[3] = new TH1F("the_ring1all","",mcaScan,xlScan,xuScan);

  theHist[0]->SetFillColor(2);
  theHist[1]->SetFillColor(2);
  theHist[2]->SetFillColor(3);
  theHist[3]->SetFillColor(3);

  for(int i=0;i<NtheHist;i++){
    theHist[i]->SetLineColor(4);
    theHist[i]->SetAxisColor(4,"X");
    theHist[i]->SetAxisColor(4,"Y");
    theHist[i]->SetLabelColor(4,"X");
    theHist[i]->SetLabelColor(4,"Y");
    theHist[i]->SetLabelSize(.07,"X");
    theHist[i]->SetLabelSize(.07,"Y");
    theHist[i]->SetLabelFont(82,"X");
    theHist[i]->SetLabelFont(82,"Y");
    theHist[i]->SetStats(kFALSE);
  }

  nLoopRing = 0;
  listTrack         = new TList();
  listRingTrack     = new TList();
  listHit           = new TList();
  listCluster       = new TList();
  listCathode       = new TList();
  listRing          = new TList();
  listRingIpo       = new TList();

  ContinueDisplay   = kTRUE;
  drawHits          = kTRUE;
  drawClusters      = kFALSE;
  drawTracks        = kTRUE;
  drawRings         = kTRUE;
  drawDigits        = kFALSE;
  drawAxis          = kFALSE;

  ViewDirection     = RFront;
  BGcolor           = RBlack;

  //----------------------------------------------------------

  //------- Set default view drection ------------------------
  Phi  = d_Phi;
  Theta = d_Theta;
  Psi = d_Psi;
  for( int ii=0; ii<3; ii++ ) { min[ii]=d_min[ii]; max[ii]=d_max[ii]; }
  //----------------------------------------------------------

  //---- Build Window ----------------------------------------
  fF0 = new TGCompositeFrame(this,900, 600, kHorizontalFrame);
  AddFrame(fF0, new TGLayoutHints(kLHintsLeft | kLHintsTop | kLHintsExpandX ));

  fF1 = new TGCompositeFrame(this, 900, 80, kHorizontalFrame);
  AddFrame(fF1, new TGLayoutHints(kLHintsLeft | kLHintsTop | kLHintsExpandX | kLHintsExpandY));

  BuildControlPanel(mode);
  //--------------------------------------------------------------  
  fITab = new TGTab(fF0,300,600);
  fF0->AddFrame(fITab, new TGLayoutHints(kLHintsLeft | kLHintsTop | kLHintsExpandX));

  //--------------------------------------------------------------  
  TIf0 = fITab->AddTab("Theta"); 
  fCanvasA = new TRootEmbeddedCanvas("canvasA",TIf0,300,580);
  TIf0->AddFrame(fCanvasA, 
		new TGLayoutHints(kLHintsTop | kLHintsLeft | kLHintsExpandX | kLHintsExpandY,5,5,0,0));
  info_canvas = fCanvasA->GetCanvas();
  info_canvas->SetEditable(false);
  info_canvas->SetFillColor(cCanvas);   
  info_canvas->Draw();
  //  info_canvas->Update();
  info_canvas->cd();
  hdisplay = new CsRCEDisplay();
  hdisplay->SetBorderSize(2);
  hdisplay->SetLineColor(4);
  hdisplay->SetFillColor(cDisplay);
  hdisplay->Divide(1,2);
  hdisplay->Draw();
  hdisplay->Modified();
  hdisplay->cd(1);
  gPad->SetGridx();
  gPad->SetGridy();
  gPad->SetFrameLineColor(4);
  hdisplay->cd(2);
  gPad->SetGridx();
  gPad->SetGridy();
  gPad->SetFrameLineColor(4);
  hdisplay->Update();
  hdisplay->cd();
  //--------------------------------------------------------------  
  TIf1 = fITab->AddTab("Phi"); 
  fCanvasD = new TRootEmbeddedCanvas("canvasD",TIf1,300,580);
  TIf1->AddFrame(fCanvasD, 
		new TGLayoutHints(kLHintsTop | kLHintsLeft | kLHintsExpandX | kLHintsExpandY,5,5,0,0));
  phi_canvas = fCanvasD->GetCanvas();
  phi_canvas->SetEditable(false);
  phi_canvas->SetFillColor(cCanvas);   
  phi_canvas->Draw();
  phi_canvas->Update();
  phi_canvas->cd();
  phidisplay = new CsRCEDisplay();
  phidisplay->SetBorderSize(2);
  phidisplay->SetLineColor(4);
  phidisplay->SetFillColor(cDisplay);
  phidisplay->Draw();
  phidisplay->Modified();

  int32 irep1;
  phidisplay->SetPhi(-90.);
  phidisplay->SetTheta(-90.);

  phidisplay->Update();
  phidisplay->cd();
  //--------------------------------------------------------------  

  fCanvasB = new TRootEmbeddedCanvas("canvasB",fF0,600,600);
  fF0->AddFrame(fCanvasB, new TGLayoutHints(kLHintsTop | kLHintsLeft,1,5,5,5));
  m_canvas = fCanvasB->GetCanvas();
  m_canvas->SetEditable(false);
  m_canvas->SetFillColor(cCanvas);   
  m_canvas->cd();
  edisplay = new CsRCEDisplay();
  edisplay->Draw();
  edisplay->Modified();
  edisplay->SetFillColor(cDisplay);
  edisplay->SetBorderSize(2);
  edisplay->SetPhi(-90.-Phi);
  edisplay->SetTheta(90.-Theta);

  m_canvas->cd();
   
  //  BuildAxis();
  m_canvas->Update();
 //--------------------------------------------------------------  
  fiTab = new TGTab(fF0,200,600);
  fF0->AddFrame(fiTab, new TGLayoutHints(kLHintsLeft | kLHintsTop | kLHintsExpandX | kLHintsExpandY));
  
  Tif0 = fiTab->AddTab("Event"); 
  fCanvasC = new TRootEmbeddedCanvas("canvasC",Tif0,200,580);
  Tif0->AddFrame(fCanvasC, 
		new TGLayoutHints(kLHintsTop | kLHintsLeft | kLHintsExpandX | kLHintsExpandY,5,5,0,0));

  info_canvas = fCanvasC->GetCanvas();
  info_canvas->SetEditable(false);
  info_canvas->SetFillColor(cCanvas);   

  idisplay = new CsRCEDisplay();
  idisplay->SetBorderSize(2);
  idisplay->SetFillColor(cDisplay);
  idisplay->Draw();
  idisplay->Modified();
  idisplay->cd();
  idisplay->Update();
  //--------------------------------------------------------------  
  Tif1 = fiTab->AddTab("Ring"); 
  fCanvasE = new TRootEmbeddedCanvas("canvasE",Tif1,200,580);
  Tif1->AddFrame(fCanvasE, 
		new TGLayoutHints(kLHintsTop | kLHintsLeft | kLHintsExpandX | kLHintsExpandY,5,5,0,0));

  ring_canvas = fCanvasE->GetCanvas();
  ring_canvas->SetEditable(false);
  ring_canvas->SetFillColor(cCanvas);   

  rdisplay = new CsRCEDisplay();
  rdisplay->SetBorderSize(2);
  rdisplay->SetFillColor(cDisplay);
  rdisplay->Draw();
  rdisplay->Modified();
  rdisplay->cd();
  rdisplay->Update();
  //--------------------------------------------------------------  
  hdisplay->cd();
  hdisplay->Update();

  MapSubwindows();
  SetWindowName("COMPASS RICH1 Offline Display");
  Resize(GetDefaultSize());  
  MapWindow();
}


void CsRCDisplay::View(){
  edisplay->cd();
  edisplay->SetFillColor(cDisplay);
  edisplay->Clear();
}
/*
void CsRCDisplay::BuildAxis(){

  fAxesItems = new TList();
  Float_t xArray[] = {0,0,0,50,0,0};
  Float_t yArray[] = {0,0,0,0,50,0};
  Float_t zArray[] = {0,0,0,0,0,50};

  TPolyLine3D *xAxis = new TPolyLine3D(2,xArray);
  xAxis->SetLineColor(4); // blue
  fAxesItems->AddLast(xAxis);

  TPolyLine3D *yAxis = new TPolyLine3D(2,yArray);
  yAxis->SetLineColor(5); // yellow
  fAxesItems->AddLast(yAxis);

  TPolyLine3D *zAxis = new TPolyLine3D(2,zArray);
  zAxis->SetLineColor(7); // sky blue
  fAxesItems->AddLast(zAxis);                                                
}

void CsRCDisplay::DrawAxis(){
  TObject *obj;
  TIter nextAxesItem(fAxesItems);
  while ((obj = nextAxesItem())) {
    if (obj->InheritsFrom(TPolyLine3D::Class())) {
      TPolyLine3D *item = (TPolyLine3D*) obj;
      item->Draw("same");
    }
  }  
}
*/

void CsRCDisplay::DrawHits(){
  if(drawHits){
    TIter nHit(listHit);
    CsRCGHit *hit;
    while ((hit = (CsRCGHit*)nHit())){
      hit->Draw();    
    }
  }
}

void CsRCDisplay::DrawClusters(){
  if(drawClusters){
    TIter nClu(listCluster);
    CsRCGTrack *clu;
    while ((clu = (CsRCGTrack*)nClu())){
      clu->Draw();    
    }
  }
}

void CsRCDisplay::DrawCathodes(){

  TIter nCat(listCathode);
  CsRCGCathode *cat;
  while ((cat = (CsRCGCathode*)nCat())){
    cat->Draw();    
  }

}

void CsRCDisplay::DrawTracks(){
  if(drawTracks){
    TIter nTrk(listTrack);
    CsRCGTrack *trk;
    while ((trk = (CsRCGTrack*)nTrk())){
      trk->Draw();    
    }
  }
}

void CsRCDisplay::DrawRings(){

  if(drawRings){
    TIter nRing(listRing);
    CsRCGRing *ring;
    while ((ring = (CsRCGRing*)nRing())){
      ring->Draw();    
    }
    TIter nRTrk(listRingTrack);
    CsRCGRingTrack *rtrk;
    while ((rtrk = (CsRCGRingTrack*)nRTrk())){
      rtrk->Draw();    
    }
  }
}

void CsRCDisplay::DrawRingsIpo(Int_t nring){

  TIter nRing(listRingIpo);
  CsRCGRing *ring;
  Int_t iring;
  while ((ring = (CsRCGRing*)nRing())){
    iring++;
    if(iring>(nring-1)*3 && iring <=(nring*3))	 ring->Draw();    
  }
}

void CsRCDisplay::SetView(float longitude, float latitude, float psi){
  int32 irep;
  view->SetView(longitude, latitude, psi,irep);
}

void CsRCDisplay::UpdateDisplay(){

  hdisplay->SetFillColor(cDisplay);
  hdisplay->cd(1);
  gPad->SetFillColor(cDisplay);
  gPad->Clear();
  hdisplay->cd(2);
  gPad->SetFillColor(cDisplay);
  gPad->Clear();
  hdisplay->Modified();
  hdisplay->Update();

  phidisplay->SetFillColor(cDisplay);
  phidisplay->cd();
  phidisplay->Clear();
  phidisplay->Modified();
  phidisplay->Update();

  idisplay->SetFillColor(cDisplay);
  idisplay->Modified();
  idisplay->Update();
  idisplay->cd();

  rdisplay->SetFillColor(cDisplay);
  rdisplay->Modified();
  rdisplay->Update();

  edisplay->SetFillColor(cDisplay);
  edisplay->Modified();
  edisplay->Update();

}
void CsRCDisplay::NextEvent(){
  SaveViewDirection();
  theAppR->Terminate(0);

  UpdateDisplay();
  rdisplay->cd();

}

void CsRCDisplay::PreviusEvent(){
}

void CsRCDisplay::NextRing(){

  
  UpdateDisplay();
  phidisplay->cd();
  gPad->Clear();
  PrintRichDisplay();
  edisplay->cd();

  CsRCRecConst *cons = CsRCRecConst::Instance();
  static float CFRefInd = cons->CFRefInd();

  CsRCGHit  *hit;
  CsRCGRing *ring0,*ring0ipo[3];
  CsRCGRing *ring1,*ring1ipo[3];

  for(int i=0;i<3;i++){
    raPart[i]=0; thePart[i]=0;
   chiPart[i]=0;likePart[i]=0;
  }

  double raIpo = 0.;
  
  for(int i=0;i<NtheHist;i++){
    theHist[i]->Reset();
  }

  CsRCEventPads* evpads = CsRCEventPads::Instance();
  CsRCDetectors * dets = CsRCDetectors::Instance();

  list<CsRCCathode*> lCathodes = dets->lCathodes();
  int nCathode = lCathodes.size();
  float nPadx = lCathodes.front()->nPadx();
  float nPady = lCathodes.front()->nPady();
  float padx  = lCathodes.front()->padx();
  float pady = lCathodes.front()->pady();
  float hCatx = padx * nPadx / 2. ;
  float hCaty = pady * nPady / 2. ;

  Int_t nhit  = 0;
  Int_t nring = 0;

  //----- reconstructed rings ----------------------------------
  CsRCEventRings* evrings = CsRCEventRings::Instance();
  list<CsRCRing*> lRings = evrings->lRings();
  int nRingEv = lRings.size();
  list<CsRCRing*>::iterator ir;

  int kRing = -1;

  nLoopRing++;

  if(nRingEv > 0) ir=lRings.begin();

  if(nRingEv>=nLoopRing){
    for( int iring=1;iring<nLoopRing;iring++) {
      ir++;
    }
  } else {
    nLoopRing=0;
    NextEvent();
    return;
  }

    nFlag1=(*ir)->flag();
    nFlag2=(*ir)->flagReco();
    nFlag3=(*ir)->flagOverThrs();
 
    ringThe     = (*ir)->the();
    ringTheReco = (*ir)->theReco();
    ringTheRFit = (*ir)->thetaRFit();

    CsRCParticle* part = (*ir)->pPartPhot()->pPart();
    double theRing = (*ir)->the();
//---- part probs --------------------------------------------------
    chiPart[0]    = (*ir)->partProbs()[12];
    chiPart[1]    = (*ir)->partProbs()[13];
    chiPart[2]    = (*ir)->partProbs()[14];
    likePart[0]   = (*ir)->partProbs()[1];
    likePart[1]   = (*ir)->partProbs()[2];
    likePart[2]   = (*ir)->partProbs()[3];
//------------------------------------------------------------------
    list<CsRCPhoton*> lPhotonsAll = (*ir)->pPartPhot()->lPhotons();
    list<CsRCPhoton*> lPhotons = (*ir)->lPhotons();
    list<CsRCPhoton*>::iterator ih;

    float partMom = part->mom();
    double raRing = 0.;

    int nPhoRing = lPhotons.size();
    nRingPhot=nPhoRing;
    float PHPho = 100.;
    //@@----------------------------
    int kPhoUp = 0;
    int kPhoDw = 0;
    double xPa[2];
    double yPa[2];
//----- particle on detectors :
    int kDetPart = (*ir)->pPartPhot()->kDetPart();
    double xPade = (*ir)->pPartPhot()->vPoPaDetW()[kDetPart].x();
    double yPade = (*ir)->pPartPhot()->vPoPaDetW()[kDetPart].y();

    for( ih=lPhotonsAll.begin(); ih!=lPhotonsAll.end(); ih++ ) {

      CsRCCluster* clu = (*ih)->ptToClu();
      int iCat = clu->ic();
      float xClu = clu->xc();
      float yClu = clu->yc();

      float thePAD=(*ih)->the();
      float PH_PAD=(*ih)->PH();

      theHist[2]->Fill(thePAD);
      theHist[3]->Fill(thePAD,PH_PAD);

    }

    Int_t mClu  = nPhoRing;
    Float_t    pClu[mClu*3];
    Float_t  phiClu[mClu*3];
    int iClu=0;

    for( ih=lPhotons.begin(); ih!=lPhotons.end(); ih++ ) {
      
      CsRCCluster* clu = (*ih)->ptToClu();
      int iCat = clu->ic();
      float xClu = clu->xc();
      float yClu = clu->yc();

      float thePAD=(*ih)->the();
      float phiPAD=(*ih)->phi();
      float PH_PAD=(*ih)->PH();

      theHist[0]->Fill(thePAD);
      theHist[1]->Fill(thePAD,PH_PAD);

      double xde = xClu;
      double yde = yClu;

      pClu[3*iClu+0] = xClu;
      pClu[3*iClu+1] = yClu;
      pClu[3*iClu+2] = 6000.;

      phiClu[3*iClu+0] = ringThe-thePAD;
      phiClu[3*iClu+1] = phiPAD;
      phiClu[3*iClu+2] = 0.;

      iClu++;

      if( part->vPade().size() == 0 ) continue;
      double xPaDet;
      double yPaDet;
      if( yde >= 0. )  {
	kPhoUp++;
	xPaDet = part->vPade()[0].x();
	yPaDet = part->vPade()[0].y();
	xPa[0] = xPaDet ;
	yPa[0] = yPaDet ;
       }
      else  {
	kPhoDw++;
	xPaDet = part->vPade()[1].x();
	yPaDet = part->vPade()[1].y();
	xPa[1] = xPaDet ;
	yPa[1] = yPaDet ;
      }
      raRing += (xde-xPaDet)*(xde-xPaDet) +
	(yde-yPaDet)*(yde-yPaDet);
    }   /* end loop on Photons */
    if( nPhoRing > 0 )  { 
      //------------- reconstructed ring radius
      raRing = sqrt( raRing / float ( nPhoRing ));

      int iPaTy=0;
      for( int kPaTy=8; kPaTy <= 14; kPaTy+=3 ) {
	float massIpo = part->mass( kPaTy );
	double betaIpo = partMom /
	  sqrt( partMom*partMom + massIpo*massIpo );
	double cosTheW = 1./ (betaIpo * CFRefInd);
	if( cosTheW <= 1.) {
	  double theIpo = acos( cosTheW ) * 1000.;
	  //--------------- mass hypo calculated ring radius (approximate)
	  if( theRing > 0. ) raIpo = raRing * theIpo / theRing;
	  raPart[iPaTy] = raIpo;
	  thePart[iPaTy] = theIpo;
	}
	iPaTy++;
      }
      
      double PI = 3.141592653;
      Int_t nPoints = 200;
      Double_t point0[nPoints*3],point0i0[nPoints*3],point0i1[nPoints*3],point0i2[nPoints*3];
      Double_t point1[nPoints*3],point1i0[nPoints*3],point1i1[nPoints*3],point1i2[nPoints*3];
      Int_t mPoints = 1;
      Float_t  mPoint0[mPoints*3];
      Float_t  mPoint1[mPoints*3];
      
      if(kPhoUp > 0) {
	for( int i = 0; i < nPoints; i++ ) {
	  float csin = sin(2*PI / (double)nPoints  * (double)i);
	  float ccos = cos(2*PI / (double)nPoints  * (double)i);
	  point0[3*i+0] = xPa[0] + raRing * ccos;
	  point0[3*i+1] = yPa[0] + raRing * csin;
	  point0[3*i+2] = 6000. ;

	  point0i0[3*i+0] = xPa[0] + raPart[0] * ccos;
	  point0i0[3*i+1] = yPa[0] + raPart[0] * csin;
	  point0i0[3*i+2] = 6000. ;
	  point0i1[3*i+0] = xPa[0] + raPart[1] * ccos;
	  point0i1[3*i+1] = yPa[0] + raPart[1] * csin;
	  point0i1[3*i+2] = 6000. ;
	  point0i2[3*i+0] = xPa[0] + raPart[2] * ccos;
	  point0i2[3*i+1] = yPa[0] + raPart[2] * csin;
	  point0i2[3*i+2] = 6000. ;

	}

	ring0 = new CsRCGRing(nPoints,point0,6,1);

	ring0ipo[0] = new CsRCGRing(nPoints,point0i0,cPion,1);
	ring0ipo[1] = new CsRCGRing(nPoints,point0i1,cKaon,1);
	ring0ipo[2] = new CsRCGRing(nPoints,point0i2,cProton,1);

	ring0->Draw();


	for(Int_t ipo=0;ipo<3;ipo++){
	  if(raPart[ipo]>0.)  ring0ipo[ipo]->Draw();
	}


	for (Int_t i=0;i<mPoints;i++){
	  mPoint0[3*i+0] = xPa[0];
	  mPoint0[3*i+1] = yPa[0];
	  mPoint0[3*i+2] = 6000.;
	}
	TPolyMarker3D *xx = new TPolyMarker3D(mPoints,mPoint0,2);
	xx->SetMarkerColor(6);
	xx->SetMarkerSize(2);
	xx->Draw();
	nring++;

      }

      if(kPhoDw > 0) {
	for( int i = 0; i < nPoints; i++ ) {
	  float csin = sin(2*PI / (double)nPoints  * (double)i);
	  float ccos = cos(2*PI / (double)nPoints  * (double)i);
	  point1[3*i+0] = xPa[1] + raRing * csin;
	  point1[3*i+1] = yPa[1] + raRing * ccos;
	  point1[3*i+2] = 6000. ;

	  point1i0[3*i+0] = xPa[1] + raPart[0] * ccos;
	  point1i0[3*i+1] = yPa[1] + raPart[0] * csin;
	  point1i0[3*i+2] = 6000. ;
	  point1i1[3*i+0] = xPa[1] + raPart[1] * ccos;
	  point1i1[3*i+1] = yPa[1] + raPart[1] * csin;
	  point1i1[3*i+2] = 6000. ;
	  point1i2[3*i+0] = xPa[1] + raPart[2] * ccos;
	  point1i2[3*i+1] = yPa[1] + raPart[2] * csin;
	  point1i2[3*i+2] = 6000. ;

	}

	ring1 = new CsRCGRing(nPoints,point1,6,1);

	ring1ipo[0] = new CsRCGRing(nPoints,point1i0,cPion,1);
	ring1ipo[1] = new CsRCGRing(nPoints,point1i1,cKaon,1);
	ring1ipo[2] = new CsRCGRing(nPoints,point1i2,cProton,1);

	ring1->Draw();
	for(Int_t ipo=0;ipo<3;ipo++){
	  if(raPart[ipo]>0.)  ring1ipo[ipo]->Draw();
	}

	for (Int_t i=0;i<mPoints;i++){
	  mPoint1[3*i+0] = xPa[1];
	  mPoint1[3*i+1] = yPa[1];
	  mPoint1[3*i+2] = 6000.;
	}
	TPolyMarker3D *xx = new TPolyMarker3D(mPoints,mPoint1,2);
	xx->SetMarkerColor(6);
	xx->SetMarkerSize(2);
	xx->Draw();
	nring++;
      }

    TPolyMarker3D *yy = new TPolyMarker3D(mClu,pClu,5);
    yy->SetMarkerColor(3);
    yy->SetMarkerSize(0.7);
    yy->Draw();


    }
    
    
    ringMom = partMom;

    double probPion = (*ir)->partProbs()[1];
    double probKaon = (*ir)->partProbs()[2];
    double probProton = (*ir)->partProbs()[3];
    double xPaDetW = 0.;
    double yPaDetW = 0.;
    bool bPade = part->vPade().size() > 0;
  //------------------------------------------------------------
    edisplay->Modified();
    edisplay->Update();
    edisplay->cd();
    hdisplay->cd(1);
    theHist[2]->Draw();
    theHist[0]->Draw("same");
    hdisplay->cd(2);
    theHist[3]->Draw();
    theHist[1]->Draw("same");
    hdisplay->Modified();
    hdisplay->Update();

    hdisplay->cd();
    DrawInfo();

    phidisplay->cd();
    int32 irep1;
#if ROOT_VERSION_CODE >= ROOT_VERSION(5,15,0)
    const Double_t rmin[]={-4,60,-4};
    const Double_t rmax[]={ 4,300, 4};
    phiview = TView::CreateView(1,rmin,rmax);
#else
    const Float_t rmin[]={-4,60,-4};
    const Float_t rmax[]={ 4,300, 4};
    phiview = new TView(rmin,rmax,1);
#endif
    phiview->SetView(Phi, Theta, Psi,irep1);
    TText *t1 = new TText(0,0,"a");

    t1->SetTextFont(42);
    t1->SetTextAlign(22);
    t1->SetTextSize(.07);
    t1->SetTextColor(4);

    TPolyMarker3D *pr = new TPolyMarker3D(mClu,phiClu,2);
    Int_t nphiax=2;
    Double_t phiax[nphiax*3];
    phiax[0]=0;phiax[1]=  0;phiax[2]=0.;
    phiax[3]=0;phiax[4]=360;phiax[5]=0.;
    TPolyLine3D *prl = new TPolyLine3D(nphiax,phiax);
    prl->SetLineColor(4);
    prl->Draw();
    TMarker3DBox *box= new TMarker3DBox(0,180,0,5,180,1,0,0);
    box->SetLineColor(4);
    box->Draw();
 
    t1->DrawText(-.01,-.90,"0");
    t1->SetTextAlign(12);
    t1->DrawText( .77,-.90,"5");
    t1->SetTextAlign(32);
    t1->DrawText(-.77,-.90,"-5");
    t1->DrawText(-.77,+.85,"360");
    t1->DrawText(-.77,+.00,"180");
    t1->DrawText(-.77,-.85,"0");

    phidisplay->cd();
    pr->SetMarkerColor(2);
    pr->SetMarkerSize(2);
    pr->Draw();

    phidisplay->Modified();
    phidisplay->Update();
    phidisplay->cd();
    edisplay->cd();

    edisplay->cd();

}
void CsRCDisplay::PreviusRing(){
  nLoopRing-=2;
  NextRing();
}

void CsRCDisplay::SaveViewDirection(){
  Phi = view->GetLongitude();
  Theta = view->GetLatitude();
  Psi = view->GetPsi();
  view->GetRange(min,max);
}

void CsRCDisplay::SetHit(CsRCGHit *hit){
  listHit->Add(hit);
}

void CsRCDisplay::SetCluster(CsRCGTrack *clu){
  listCluster->Add(clu);
}

void CsRCDisplay::SetTrack(CsRCGTrack *trk){
  listTrack->Add(trk);
}

void CsRCDisplay::SetRingTrack(CsRCGRingTrack *rtrk){
  listRingTrack->Add(rtrk);
}

void CsRCDisplay::SetCathode(CsRCGCathode *cat){
  listCathode->Add(cat);
}

void CsRCDisplay::SetRing(CsRCGRing *ring){
  listRing->Add(ring);
}

void CsRCDisplay::SetRingIpo(CsRCGRing *ring[3]){
  for(int i=0;i<3;i++){
    listRingIpo->Add(ring[i]);
  }
}



void CsRCDisplay::OpenEvent(){

  if(Mode == RStandalone) ReadEvent();

  GetRichDisplay();

  theAppR->Run(kTRUE);    
  listTrack->Clear();
  listRingTrack->Clear();
  listHit->Clear();
  listCluster->Clear();
  listRing->Clear();
  ClearParticleInfo();
  hdisplay->cd(1);
  gPad->Clear();
  hdisplay->cd(2);
  gPad->Clear();
  for(int i=0;i<3;i++){
    raPart[i]=0; thePart[i]=0;
   chiPart[i]=0;likePart[i]=0;
  }
}

void CsRCDisplay::Quit(){
  ContinueDisplay = kFALSE;
  theAppR->Terminate(0);
  cout << endl;
  string str = 
    "RICHONE, CsRCDisplay::Quit() : end of Event Display";
  CsErrLog::Instance()->mes( elFatal, str );
}


void CsRCDisplay::PrintRichDisplay(){

  ClearViewPad();
  SetRange(min[0],min[1],min[2],max[0],max[1],max[2]);    
  SetView(Phi,Theta,Psi);

  DrawCathodes();
  DrawHits();
  DrawClusters();
  DrawTracks();
  DrawRings();

  edisplay->Modified();
  edisplay->Update();
  edisplay->cd();
  DrawInfo();

}

void CsRCDisplay::GetRichDisplay(){

  CsRCRecConst *cons = CsRCRecConst::Instance();
  static float CFRefInd = cons->CFRefInd();
  static int mcaScan = cons->mcaScan();

  CsRCGHit       *hit;
  CsRCGCathode   *cat;
  CsRCGTrack     *trk;
  CsRCGTrack     *clu;
  CsRCGRingTrack *rtrk;

  CsRCGRing *ring0,*ring0ipo[3];
  CsRCGRing *ring1,*ring1ipo[3];

  CsRCEventPads * evpads = CsRCEventPads::Instance();
  CsRCDetectors   * dets = CsRCDetectors::Instance();

  CsRCExeKeys *key = CsRCExeKeys::Instance();

  list<CsRCCathode*> lCathodes = dets->lCathodes();
  int nCathode = lCathodes.size();
  float nPadx  = lCathodes.front()->nPadx();
  float nPady  = lCathodes.front()->nPady();
  float padx   = lCathodes.front()->padx();
  float pady   = lCathodes.front()->pady();
  float hCatx  = padx * nPadx / 2. ;
  float hCaty  = pady * nPady / 2. ;
   
  list<CsRCCathode*>::const_iterator ic;

  Int_t nhit  = 0;
  Int_t ncat  = 0;
  Int_t ntrk  = 0;
  Int_t nring = 0;

  CsRCCathode* catP = NULL;


  static bool firstCall = true;
  if( firstCall ) { 
    firstCall = false;
    key->acknoMethod( "EventROOTDisplay" );
    for( ic=lCathodes.begin(); ic!=lCathodes.end(); ic++ ) {
      Hep3Vector	  vOffCatW_ = ( (*ic)->vOffCatW() ) ;
      float xcat = vOffCatW_.x();
      float ycat = vOffCatW_.y();
      float zcat = vOffCatW_.z();
      cat = new CsRCGCathode(nhit,xcat,ycat,zcat,hCatx,hCaty,1.,4);
      SetCathode(cat);
      ncat++;
    }
  }

  float PHPad = 1.;
  //@@----------------------
  list<CsRCPad*> lPadAll = CsRCEventPads::Instance()->lPads();
  list<CsRCPad*>::iterator ia;

  for( ia=lPadAll.begin(); ia!=lPadAll.end(); ia++ ) {
    int iCat = (*ia)->ic();
    catP = dets->ptrToCat( iCat );
    if( catP->isPMT() ) {
      int nPPadx = catP->nPPadx();
      int nPPady = catP->nPPady();
      int ixPa = (*ia)->ix();
      int iyPa = (*ia)->iy();
      PHPad = (*ia)->PH();
      double xde = (ixPa*catP->ppadx() - catP->hCatx()+catP->ppadx()/2) + dets->vOffCatW( iCat ).x();
      double yde = (iyPa*catP->ppady() - catP->hCaty()+catP->ppady()/2) + dets->vOffCatW( iCat ).y();
      double zde = dets->vOffCatW( iCat ).z();
      hit = new CsRCGHit(nhit,xde,yde,zde,catP->ppadx()/2,catP->ppadx()/2,1.,2);
      SetHit(hit);
      nhit++;
    } else {
      int ixPa = (*ia)->ix();
      int iyPa = (*ia)->iy();
      PHPad = (*ia)->PH();
      double xde = (ixPa*padx - hCatx+padx/2) + dets->vOffCatW( iCat ).x();
      double yde = (iyPa*pady - hCaty+pady/2) + dets->vOffCatW( iCat ).y();
      double zde = dets->vOffCatW( iCat ).z();
      hit = new CsRCGHit(nhit,xde,yde,zde,padx/2,pady/2,1.,2);
      SetHit(hit);
      nhit++;
    }
  }

    list<CsRCCluster*> lClusters =  CsRCEventClusters::Instance()->lClusters();
    list<CsRCCluster*>::iterator cp;

    int iClu=0;
    Float_t pCluster[3];
    for( cp=lClusters.begin(); cp!=lClusters.end(); cp++ ) {
      int   iCat = (*cp)->ic();
      float xClu = (*cp)->xc();
      float yClu = (*cp)->yc();

      pCluster[0] = xClu;
      pCluster[1] = yClu;
      pCluster[2] = 6000.;

      Float_t size=0.7;
      clu = new CsRCGTrack(1,pCluster,5,5,size);
      SetCluster(clu);

      iClu++;
    }

  //------------------------------------------------------------
  int kEvent = CsRichOne::Instance()->kEvent();
  CsRCEventAnalysis* ana = CsRCEventAnalysis::Instance();
  // ---- reconstructed tracks ---------------------------------
  CsRCEventParticles* evparts = CsRCEventParticles::Instance();
  list<CsRCParticle*> lParticles = evparts->lParticles();
  list<CsRCParticle*>::iterator ipart;
  nParticles = lParticles.size();

  for( ipart=lParticles.begin(); ipart!=lParticles.end(); ipart++ ){
    if( (*ipart)->vPade().size() == 0 ) continue;

      Int_t mPoints = 2;
      Float_t  PaDet[3];
      for (Int_t i=0;i<mPoints;i++){
	double zzmir=8700;
	double yout=(*ipart)->ya() + (*ipart)->md() * (zzmir-(*ipart)->za());

	//	cout << " yout " << (*ipart)->yo() << " " << yout << " flag " << (*ipart)->flag() << 
	//	  " "<< (*ipart)->vPade()[i].y() << endl;

	if(yout>0 && i==0){
	  PaDet[0] = (*ipart)->vPade()[i].x();
	  PaDet[1] = (*ipart)->vPade()[i].y();
	  PaDet[2] = 6000.;
	} else if (yout<0 && i==1){
	  PaDet[0] = (*ipart)->vPade()[i].x();
	  PaDet[1] = (*ipart)->vPade()[i].y();
	  PaDet[2] = 6000.;
	}
      }
      trk = new CsRCGTrack(1,PaDet,2,5,3);
      SetTrack(trk);
      ntrk++;
  }
  //------------------------------------------------------------

  //----- reconstructed rings ----------------------------------
  CsRCEventRings* evrings = CsRCEventRings::Instance();
  list<CsRCRing*> lRings = evrings->lRings();
  int nRingEv = lRings.size();
  list<CsRCRing*>::iterator ir;
  
  nRings=nRingEv;
  
  int kRing = -1;
  for( ir=lRings.begin(); ir!=lRings.end(); ir++ ) {

    //if( !(*ir)->flag() ) continue;
    kRing++;                                  //   020906
    CsRCParticle* part = (*ir)->pPartPhot()->pPart();
    double theRing = (*ir)->the();
    list<CsRCPhoton*> lPhotons = (*ir)->lPhotons();
    list<CsRCPhoton*>::iterator ih;
    int nPhoRing = lPhotons.size();

    float partMom = part->mom();
    double raIpo = 0.;
    
    float PHPho = 100.;
    //@@----------------------------
    double raRing = 0.;
    int kPhoUp = 0;
    int kPhoDw = 0;
    double xPa[2];
    double yPa[2];
 //----- particle on detectors :
    int kDetPart = (*ir)->pPartPhot()->kDetPart();
    double xPade = (*ir)->pPartPhot()->vPoPaDetW()[kDetPart].x();
    double yPade = (*ir)->pPartPhot()->vPoPaDetW()[kDetPart].y();

    
    Int_t mClu = nPhoRing;
    Float_t  pClu[mClu*3];
    int iClu=0;

    for( ih=lPhotons.begin(); ih!=lPhotons.end(); ih++ ) {
      
      CsRCCluster* clu = (*ih)->ptToClu();
      int iCat = clu->ic();
      float xClu = clu->xc();
      float yClu = clu->yc();
      double xde = xClu;
      double yde = yClu;

      pClu[3*iClu+0] = xClu;
      pClu[3*iClu+1] = yClu;
      pClu[3*iClu+2] = 6000.;

      iClu++;

      if( part->vPade().size() == 0 ) continue;
      double xPaDet;
      double yPaDet;
      if( yde >= 0. )  {
	kPhoUp++;
	xPaDet = part->vPade()[0].x();
	yPaDet = part->vPade()[0].y();
	xPa[0] = xPaDet ;
	yPa[0] = yPaDet ;
       }
      else  {
	kPhoDw++;
	xPaDet = part->vPade()[1].x();
	yPaDet = part->vPade()[1].y();
	xPa[1] = xPaDet ;
	yPa[1] = yPaDet ;
      }
      raRing += (xde-xPaDet)*(xde-xPaDet) +
	(yde-yPaDet)*(yde-yPaDet);
    }  
    if( nPhoRing > 0 )  { 
      //------------- reconstructed ring radius
      raRing = sqrt( raRing / float ( nPhoRing ));
      int iPaTy=0;
      for( int kPaTy=8; kPaTy <= 14; kPaTy+=3 ) {
	float massIpo = part->mass( kPaTy );
	double betaIpo = partMom /
	  sqrt( partMom*partMom + massIpo*massIpo );
	double cosTheW = 1./ (betaIpo * CFRefInd);
	if( cosTheW <= 1.) {
	  double theIpo = acos( cosTheW ) * 1000.;
	  //--------------- mass hypo calculated ring radius (approximate)
	  if( theRing > 0. ) raIpo = raRing * theIpo / theRing;
	  raPart[iPaTy] = raIpo;
	  thePart[iPaTy] = theIpo;
	}
	iPaTy++;
      }
      
      double PI = 3.141592653;
      Int_t nPoints = 200;
      Double_t point0[nPoints*3],point0i0[nPoints*3],point0i1[nPoints*3],point0i2[nPoints*3];
      Double_t point1[nPoints*3],point1i0[nPoints*3],point1i1[nPoints*3],point1i2[nPoints*3];
      Int_t mPoints = 1;
      Float_t  mPoint0[mPoints*3];
      Float_t  mPoint1[mPoints*3];
      
      if(kPhoUp > 0) {
	for( int i = 0; i < nPoints; i++ ) {
	  float csin = sin(2*PI / (double)nPoints  * (double)i);
	  float ccos = cos(2*PI / (double)nPoints  * (double)i);
	  point0[3*i+0] = xPa[0] + raRing * ccos;
	  point0[3*i+1] = yPa[0] + raRing * csin;
	  point0[3*i+2] = 6000. ;

	  point0i0[3*i+0] = xPa[0] + raPart[0] * ccos;
	  point0i0[3*i+1] = yPa[0] + raPart[0] * csin;
	  point0i0[3*i+2] = 6000. ;
	  point0i1[3*i+0] = xPa[0] + raPart[1] * ccos;
	  point0i1[3*i+1] = yPa[0] + raPart[1] * csin;
	  point0i1[3*i+2] = 6000. ;
	  point0i2[3*i+0] = xPa[0] + raPart[2] * ccos;
	  point0i2[3*i+1] = yPa[0] + raPart[2] * csin;
	  point0i2[3*i+2] = 6000. ;

	}

	ring0ipo[0] = new CsRCGRing(nPoints,point0i0,cPion,1);
	ring0ipo[1] = new CsRCGRing(nPoints,point0i1,cKaon,1);
	ring0ipo[2] = new CsRCGRing(nPoints,point0i2,cProton,1);

	SetRingIpo(ring0ipo);

	if( (*ir)->flag() ) {
	  ring0 = new CsRCGRing(nPoints,point0,3,1);
	  SetRing(ring0);

	  for (Int_t i=0;i<mPoints;i++){
	    mPoint0[3*i+0] = xPa[0];
	    mPoint0[3*i+1] = yPa[0];
	    mPoint0[3*i+2] = 6000.;
	  }
	  Float_t size=2;
	  rtrk = new CsRCGRingTrack(mPoints,mPoint0,3,3,size);
	  SetRingTrack(rtrk);
	  nring++;
	}

      }

      if(kPhoDw > 0) {
	for( int i = 0; i < nPoints; i++ ) {
	  float csin = sin(2*PI / (double)nPoints  * (double)i);
	  float ccos = cos(2*PI / (double)nPoints  * (double)i);
	  point1[3*i+0] = xPa[1] + raRing * csin;
	  point1[3*i+1] = yPa[1] + raRing * ccos;
	  point1[3*i+2] = 6000. ;

	  point1i0[3*i+0] = xPa[1] + raPart[0] * ccos;
	  point1i0[3*i+1] = yPa[1] + raPart[0] * csin;
	  point1i0[3*i+2] = 6000. ;
	  point1i1[3*i+0] = xPa[1] + raPart[1] * ccos;
	  point1i1[3*i+1] = yPa[1] + raPart[1] * csin;
	  point1i1[3*i+2] = 6000. ;
	  point1i2[3*i+0] = xPa[1] + raPart[2] * ccos;
	  point1i2[3*i+1] = yPa[1] + raPart[2] * csin;
	  point1i2[3*i+2] = 6000. ;

	}

	ring1ipo[0] = new CsRCGRing(nPoints,point1i0,cPion,1);
	ring1ipo[1] = new CsRCGRing(nPoints,point1i1,cKaon,1);
	ring1ipo[2] = new CsRCGRing(nPoints,point1i2,cProton,1);

	SetRingIpo(ring1ipo);

	if( (*ir)->flag() ) {
	  ring1 = new CsRCGRing(nPoints,point1,3,1);
	  SetRing(ring1);

	  for (Int_t i=0;i<mPoints;i++){
	    mPoint1[3*i+0] = xPa[1];
	    mPoint1[3*i+1] = yPa[1];
	    mPoint1[3*i+2] = 6000.;
	  }
	  rtrk = new CsRCGRingTrack(mPoints,mPoint1,3,3,2.0);
	  SetRingTrack(rtrk);
	  nring++;
	}

      }

      if( (*ir)->flag() ) {
	rtrk = new CsRCGRingTrack(mClu,pClu,5,3,0.8);
	SetRingTrack(rtrk);
      }
    }
  } 
  
  PrintRichDisplay();
  //------------------------------------------------------------


}

bool CsRCDisplay::GetContinueDisplay(){
  return ContinueDisplay;
}

void CsRCDisplay::ClearViewPad(){
  edisplay->cd();
  edisplay->Clear();
  edisplay->SetFillColor(cDisplay);
#if ROOT_VERSION_CODE >= ROOT_VERSION(5,15,0)
  view = TView::CreateView(1);
#else
  view = new TView(1);    
#endif
  phidisplay->cd();
  phidisplay->Clear();
  phidisplay->SetFillColor(cDisplay);
#if ROOT_VERSION_CODE >= ROOT_VERSION(5,15,0)
  phiview = TView::CreateView(1);
#else
  phiview = new TView(1);    
#endif

  edisplay->cd();
}

CsRCDisplay::~CsRCDisplay(){
}

Bool_t CsRCDisplay::ProcessMessage(Long_t msg, Long_t parm1, Long_t parm2)
{
  uint32 retval; 

  switch(GET_MSG(msg)) {
  case kC_COMMAND:
    switch(GET_SUBMSG(msg)) {
    case kCM_BUTTON: 
      switch(parm1){
      case RF_QUIT: 
	Quit();
	break;
      case RF_NEXT: 
	NextEvent();
	break;
      case RF_PREVIUS: 
	PreviusEvent();
	break;
      case RF_NEXT_RING: 
	NextRing();
	break;
      case RF_PREVIUS_RING: 
	PreviusRing();
	break;
      case RF_ZOOM_IN: 
	ZoomIn();       
	break;
      case RF_ZOOM_OUT: 
	ZoomOut();
	break;
      case RF_EXPLOR_UP: 
	ExploreUp();
	break;
      case RF_EXPLOR_DOWN: 
	ExploreDown();
	break;
      case RF_EXPLOR_RIGHT: 
	ExploreRight();
	break;
      case RF_EXPLOR_LEFT: 
	ExploreLeft();
	break;
      case RF_RESET_VIEW: 
	ResetView();
	break;
      case RF_RELOAD_AF: 
	ReloadAF();
	break;      
      case RF_TOP_VIEW: 
	TopView();
	break;      
      case RF_SIDE_VIEW: 
	SideView();
	break;      
      case RF_FRONT_VIEW: 
	FrontView();
	break;      
      case RF_GRAPHIC_OUT: 
	MakeGraphicFile();
	break;      
      case RF_COLOR_BACK: 
	BackColor();
	break;      
      case RF_X3D: 
	X3D();
	break;      
      case RF_SCALE_INC: 
	ScaleInc();
	break;      
      case RF_SCALE_DEC: 
	ScaleDec();
	break;      
      case RF_DIGIT_ON_OFF: 
	switch(drawRings){
	case kFALSE:
	  RingOn();
	  break;      
	case kTRUE:
	  RingOff();
	  break;      
	default:
	  break;
	}
	break;
      case RF_HIT_ON_OFF: 
	switch(drawHits){
	case kFALSE:
	  HitOn();
	  break;      
	case kTRUE:
	  HitOff();
	  break;      
	default:
	  break;
	}
	break;
      case RF_CLUSTER_ON_OFF: 
	switch(drawClusters){
	case kFALSE:
	  ClusterOn();
	  break;      
	case kTRUE:
	  ClusterOff();
	  break;      
	default:
	  break;
	}
	break;
      case RF_TRACK_ON_OFF: 
	switch(drawTracks){
	case kFALSE:
	  TrackOn();
	  break;      
	case kTRUE:
	  TrackOff();
	  break;      
	default:
	  break;
	}
	break;
      default:
	break;
      }
    case kCM_RADIOBUTTON: 
      switch(parm1){
      default:
	break;    
      }
    }
  default:
    break;
  }
  return kTRUE; 
}

void CsRCDisplay::ZoomIn(){
  Float_t Cmin,Cmax,Center; 
  Phi = view->GetLongitude();
  Theta = view->GetLatitude();
  view->GetRange(min,max);
  Psi = view->GetPsi();
  for(int32 i=0;i<3;i++) {
    Center = (min[i]+max[i])/2.;
    Cmin = ZoomFactor*TMath::Abs(Center - min[i]);
    Cmax = ZoomFactor*TMath::Abs(Center - max[i]);
    min[i] = Center - Cmin;
    max[i] = Center + Cmax;
  } 
  Redraw();
}

void CsRCDisplay::ZoomOut(){
  Float_t Cmin,Cmax,Center; 
  Phi = view->GetLongitude();
  Theta = view->GetLatitude();
  Psi = view->GetPsi();
  view->GetRange(min,max);
  for(int32 i=0;i<3;i++) {
    Center = (min[i]+max[i])/2.;
    Cmin = 1. / ZoomFactor*TMath::Abs(Center - min[i]);
    Cmax = 1. / ZoomFactor*TMath::Abs(Center - max[i]);
    min[i] = Center - Cmin;
    max[i] = Center + Cmax;
  }  
  Redraw();
}

void CsRCDisplay::ExploreUp(){
  Phi = view->GetLongitude();
  Theta = view->GetLatitude();
  Psi = view->GetPsi();
  view->GetRange(min,max);
  if(ViewDirection == RTop){
    float wind=max[0]-min[0];
    min[0]+=wind*ExplorFactor;
    max[0]+=wind*ExplorFactor;
  }else{
    if(ViewDirection == RSide || ViewDirection == RFront){
      float wind=max[1]-min[1];
      min[1]+=wind*ExplorFactor;
      max[1]+=wind*ExplorFactor;
    }      
  }
  Redraw();
}

void CsRCDisplay::ExploreDown(){
  Phi = view->GetLongitude();
  Theta = view->GetLatitude();
  Psi = view->GetPsi();
  view->GetRange(min,max);
  if(ViewDirection == RTop){
    float wind=max[0]-min[0];
    min[0]-=wind*ExplorFactor;
    max[0]-=wind*ExplorFactor;
  }else{
    if(ViewDirection == RSide || ViewDirection == RFront){
      float wind=max[1]-min[1];
      min[1]-=wind*ExplorFactor;
      max[1]-=wind*ExplorFactor;
    }      
  }
  Redraw();
}

void CsRCDisplay::ExploreRight(){
  Phi = view->GetLongitude();
  Theta = view->GetLatitude();
  Psi = view->GetPsi();
  view->GetRange(min,max);
  if(ViewDirection == RFront){
    float wind=max[0]-min[0];
    min[0]+=wind*ExplorFactor;
    max[0]+=wind*ExplorFactor;
  }else{
    if(ViewDirection == RTop || ViewDirection == RSide){
      float wind=max[1]-min[1];
      min[2]-=wind*ExplorFactor;
      max[2]-=wind*ExplorFactor;
    }      
  }
  Redraw();
}

void CsRCDisplay::ExploreLeft(){
  Phi = view->GetLongitude();
  Theta = view->GetLatitude();
  Psi = view->GetPsi();
  view->GetRange(min,max);
  if(ViewDirection == RFront){
    float wind=max[0]-min[0];
    min[0]-=wind*ExplorFactor;
    max[0]-=wind*ExplorFactor;
  }else{
    if(ViewDirection == RTop || ViewDirection == RSide){
      float wind=max[1]-min[1];
      min[2]+=wind*ExplorFactor;
      max[2]+=wind*ExplorFactor;
    }      
  }
  Redraw();
}

void CsRCDisplay::ResetView(){
  ClearViewPad();
  edisplay->SetPhi(-90.-d_Phi);
  edisplay->SetTheta(90.-d_Theta);
  for( int ii=0; ii<3; ii++ ) { min[ii]=d_min[ii]; max[ii]=d_max[ii]; }
  PrintRichDisplay();
}

void CsRCDisplay::SetRange(Float_t xmin,Float_t ymin,Float_t zmin,
			   Float_t xmax,Float_t ymax,Float_t zmax){
  view->SetRange(xmin,ymin,zmin,xmax,ymax,zmax);
}

void CsRCDisplay::Redraw(){
  SetView(Phi,Theta,Psi);
  PrintRichDisplay();

  if(nLoopRing>0) {
    nLoopRing--;
    NextRing();
  }

}


void CsRCDisplay::ReloadAF(){
  SaveViewDirection();
  ClearViewPad();
  SetRange(min[0],min[1],min[2],max[0],max[1],max[2]);
  SetView(Phi,Theta,Psi);
  PrintRichDisplay();
}


void CsRCDisplay::TopView(){
  ViewDirection = RTop;
  SaveViewDirection();
  Phi = -90;
  Theta = -90;
  Psi = -90;
  ClearViewPad();
  edisplay->SetPhi(-90.-Phi);
  edisplay->SetTheta(90.-Theta);
  SetRange(min[0],min[1],min[2],max[0],max[1],max[2]);
  SetView(Phi,Theta,Psi);
  PrintRichDisplay();
}

void CsRCDisplay::SideView(){
  ViewDirection = RSide;
  SaveViewDirection();
  Phi = 0;
  Theta = -90;
  Psi = -90;
  ClearViewPad();
  edisplay->SetPhi(-90.-Phi);
  edisplay->SetTheta(90.-Theta);
  SetRange(min[0],min[1],min[2],max[0],max[1],max[2]);
  SetView(Phi,Theta,Psi);
  PrintRichDisplay();
}

void CsRCDisplay::FrontView(){
  ViewDirection = RFront;
  SaveViewDirection();
  Phi = 0;
  Theta = -180;
  Psi = -90;
  ClearViewPad();
  edisplay->SetPhi(-90.-Phi);
  edisplay->SetTheta(90.-Theta);
  SetRange(min[0],min[1],min[2],max[0],max[1],max[2]);
  SetView(Phi,Theta,Psi);
  PrintRichDisplay();
}

void CsRCDisplay::DrawInfo(){

  idisplay->cd();
  idisplay->Clear();

  int32 cTextParticles   = 7;
  int32 cTextNumParticle = 7;
  int32 cTextPartPlus    = 3;
  int32 cTextPartMinus   = 5;
  int32 cTextPartNeutral = 4;
  Float_t xLabelParticle = 0.1,yLabelParticle = 0.95;
  Float_t xTextParticle  = 0.9,yTextParticle;
  Float_t sTextDY = 0.06;
  Float_t sTextParticles = 0.08;
  Float_t sTextRunNumber = 0.08;

  Char_t note[100];

  sprintf(note,"Run ");
  text = new TText(xLabelParticle,yLabelParticle,note);
  text->SetTextSize(sTextRunNumber);
  text->SetTextColor(0);
  text->Draw();

  sprintf(note,"%Lu",RunNumber);
  yTextParticle = yLabelParticle; 
  text = new TText(xTextParticle,yTextParticle,note);
  text->SetTextAlign(31);
  text->SetTextSize(sTextRunNumber);
  text->SetTextColor(0);
  text->Draw();

  sprintf(note,"Event ");
  yLabelParticle-=sTextDY; 
  text = new TText(xLabelParticle,yLabelParticle,note);
  text->SetTextSize(sTextRunNumber);
  text->SetTextColor(0);
  text->Draw();

  sprintf(note,"%Lu",EventNumberInRun);
  yTextParticle = yLabelParticle; 
  text = new TText(xTextParticle,yTextParticle,note);
  text->SetTextAlign(31);
  text->SetTextSize(sTextRunNumber);
  text->SetTextColor(0);
  text->Draw();

  sprintf(note,"Event Burst");
  yLabelParticle-=sTextDY; 
  text = new TText(xLabelParticle,yLabelParticle,note);
  text->SetTextSize(sTextRunNumber);
  text->SetTextColor(0);
  text->Draw();

  sprintf(note,"%Lu",EventNumberInBurst);
  yTextParticle = yLabelParticle; 
  text = new TText(xTextParticle,yTextParticle,note);
  text->SetTextAlign(31);
  text->SetTextSize(sTextRunNumber);
  text->SetTextColor(0);
  text->Draw();

  sprintf(note,"Tracks");
  yLabelParticle-=sTextDY; 
  text = new TText(xLabelParticle,yLabelParticle,note);
  text->SetTextSize(sTextParticles);
  text->SetTextColor(cTextParticles);
  text->Draw();

  sprintf(note,"%d",nParticles);
  yTextParticle = yLabelParticle; 
  text = new TText(xTextParticle,yTextParticle,note);
  text->SetTextAlign(31);
  text->SetTextSize(sTextParticles);
  text->SetTextColor(cTextNumParticle);
  text->Draw();

  sprintf(note,"Rings");
  yLabelParticle-=sTextDY; 
  text = new TText(xLabelParticle,yLabelParticle,note);
  text->SetTextSize(sTextParticles);
  text->SetTextColor(cTextPartPlus);
  text->Draw();

  sprintf(note,"%d",nRings);
  yTextParticle = yLabelParticle;
  text = new TText(xTextParticle,yTextParticle,note);
  text->SetTextAlign(31);
  text->SetTextSize(sTextParticles);
  text->SetTextColor(cTextNumParticle);
  text->Draw();

  sprintf(note,"Current Ring:");
  yLabelParticle-=sTextDY; 
  text = new TText(xLabelParticle,yLabelParticle,note);
  text->SetTextSize(sTextParticles);
  text->SetTextColor(cTextPartMinus);
  text->Draw();

  sprintf(note,"%d",nLoopRing);
  yTextParticle = yLabelParticle;
  text = new TText(xTextParticle,yTextParticle,note);
  text->SetTextAlign(31);
  text->SetTextSize(sTextParticles);
  text->SetTextColor(cTextNumParticle);
  text->Draw();

  sprintf(note,"Photons:");
  yLabelParticle-=sTextDY; 
  text = new TText(xLabelParticle,yLabelParticle,note);
  text->SetTextSize(sTextParticles);
  text->SetTextColor(cTextPartPlus);
  text->Draw();

  sprintf(note,"%d",nRingPhot);
  yTextParticle = yLabelParticle;
  text = new TText(xTextParticle,yTextParticle,note);
  text->SetTextAlign(31);
  text->SetTextSize(sTextParticles);
  text->SetTextColor(cTextNumParticle);
  text->Draw();


  sprintf(note,"flag:");
  yLabelParticle-=sTextDY; 
  text = new TText(xLabelParticle,yLabelParticle,note);
  text->SetTextSize(sTextParticles);
  text->SetTextColor(cTextPartNeutral);
  text->Draw();

  sprintf(note,"%d",nFlag1);
  yTextParticle = yLabelParticle;
  text = new TText(xTextParticle,yTextParticle,note);
  text->SetTextAlign(31);
  text->SetTextSize(sTextParticles);
  text->SetTextColor(cTextNumParticle);
  text->Draw();

  sprintf(note,"flagReco:");
  yLabelParticle-=sTextDY; 
  text = new TText(xLabelParticle,yLabelParticle,note);
  text->SetTextSize(sTextParticles);
  text->SetTextColor(cTextPartPlus);
  text->Draw();

  sprintf(note,"%d",nFlag2);
  yTextParticle = yLabelParticle;
  text = new TText(xTextParticle,yTextParticle,note);
  text->SetTextAlign(31);
  text->SetTextSize(sTextParticles);
  text->SetTextColor(cTextNumParticle);
  text->Draw();

  sprintf(note,"flagOverThrs");
  yLabelParticle-=sTextDY; 
  text = new TText(xLabelParticle,yLabelParticle,note);
  text->SetTextSize(sTextParticles);
  text->SetTextColor(cTextPartMinus);
  text->Draw();

  sprintf(note,"%d",nFlag3);
  yTextParticle = yLabelParticle;
  text = new TText(xTextParticle,yTextParticle,note);
  text->SetTextAlign(31);
  text->SetTextSize(sTextParticles);
  text->SetTextColor(cTextNumParticle);
  text->Draw();


  sprintf(note,"Momentum");
  yLabelParticle-=sTextDY; 
  text = new TText(xLabelParticle,yLabelParticle,note);
  text->SetTextSize(sTextParticles);
  text->SetTextColor(cTextPartMinus);
  text->Draw();

  sprintf(note,"%.2f",ringMom);
  yTextParticle = yLabelParticle;
  text = new TText(xTextParticle,yTextParticle,note);
  text->SetTextAlign(31);
  text->SetTextSize(sTextParticles);
  text->SetTextColor(cTextNumParticle);
  text->Draw();


  sprintf(note,"Theta   ");
  yLabelParticle-=sTextDY; 
  text = new TText(xLabelParticle,yLabelParticle,note);
  text->SetTextSize(sTextParticles);
  text->SetTextColor(cTextPartMinus);
  text->Draw();

  sprintf(note,"%.2f",ringThe);
  yTextParticle = yLabelParticle;
  text = new TText(xTextParticle,yTextParticle,note);
  text->SetTextAlign(31);
  text->SetTextSize(sTextParticles);
  text->SetTextColor(cTextNumParticle);
  text->Draw();

  sprintf(note,"ThetaRec  ");
  yLabelParticle-=sTextDY; 
  text = new TText(xLabelParticle,yLabelParticle,note);
  text->SetTextSize(sTextParticles);
  text->SetTextColor(cTextPartMinus);
  text->Draw();

  sprintf(note,"%.2f",ringTheReco);
  yTextParticle = yLabelParticle;
  text = new TText(xTextParticle,yTextParticle,note);
  text->SetTextAlign(31);
  text->SetTextSize(sTextParticles);
  text->SetTextColor(cTextNumParticle);
  text->Draw();

  sprintf(note,"Theta Fit ");
  yLabelParticle-=sTextDY; 
  text = new TText(xLabelParticle,yLabelParticle,note);
  text->SetTextSize(sTextParticles);
  text->SetTextColor(cTextPartMinus);
  text->Draw();

  sprintf(note,"%.2f",ringTheRFit);
  yTextParticle = yLabelParticle;
  text = new TText(xTextParticle,yTextParticle,note);
  text->SetTextAlign(31);
  text->SetTextSize(sTextParticles);
  text->SetTextColor(cTextNumParticle);
  text->Draw();

  idisplay->Modified();
  idisplay->Update();
  idisplay->cd();


  rdisplay->cd();
  rdisplay->Clear();

  xLabelParticle = 0.1;
  yLabelParticle = 0.95;
  xTextParticle  = 0.9;

  sprintf(note,"Momentum");
  yLabelParticle-=sTextDY; 
  text = new TText(xLabelParticle,yLabelParticle,note);
  text->SetTextSize(sTextParticles);
  text->SetTextColor(cTextPartMinus);
  text->Draw();

  sprintf(note,"%.2f",ringMom);
  yTextParticle = yLabelParticle;
  text = new TText(xTextParticle,yTextParticle,note);
  text->SetTextAlign(31);
  text->SetTextSize(sTextParticles);
  text->SetTextColor(cTextNumParticle);
  text->Draw();

  sprintf(note,"Theta   ");
  yLabelParticle-=sTextDY; 
  text = new TText(xLabelParticle,yLabelParticle,note);
  text->SetTextSize(sTextParticles);
  text->SetTextColor(cTextPartMinus);
  text->Draw();

  sprintf(note,"%.2f",ringThe);
  yTextParticle = yLabelParticle;
  text = new TText(xTextParticle,yTextParticle,note);
  text->SetTextAlign(31);
  text->SetTextSize(sTextParticles);
  text->SetTextColor(cTextNumParticle);
  text->Draw();

  sprintf(note,"Theta Fit ");
  yLabelParticle-=sTextDY; 
  text = new TText(xLabelParticle,yLabelParticle,note);
  text->SetTextSize(sTextParticles);
  text->SetTextColor(cTextPartMinus);
  text->Draw();

  sprintf(note,"%.2f",ringTheRFit);
  yTextParticle = yLabelParticle;
  text = new TText(xTextParticle,yTextParticle,note);
  text->SetTextAlign(31);
  text->SetTextSize(sTextParticles);
  text->SetTextColor(cTextNumParticle);
  text->Draw();

  sprintf(note,"Theta Pion ");
  yLabelParticle-=sTextDY; 
  text = new TText(xLabelParticle,yLabelParticle,note);
  text->SetTextSize(sTextParticles);
  text->SetTextColor(cPion);
  text->Draw();

  sprintf(note,"%.2f",thePart[0]);
  yTextParticle = yLabelParticle;
  text = new TText(xTextParticle,yTextParticle,note);
  text->SetTextAlign(31);
  text->SetTextSize(sTextParticles);
  text->SetTextColor(cPion);
  text->Draw();

  sprintf(note,"Theta Kaon ");
  yLabelParticle-=sTextDY; 
  text = new TText(xLabelParticle,yLabelParticle,note);
  text->SetTextSize(sTextParticles);
  text->SetTextColor(cKaon);
  text->Draw();

  sprintf(note,"%.2f",thePart[1]);
  yTextParticle = yLabelParticle;
  text = new TText(xTextParticle,yTextParticle,note);
  text->SetTextAlign(31);
  text->SetTextSize(sTextParticles);
  text->SetTextColor(cKaon);
  text->Draw();

  sprintf(note,"Theta Proton ");
  yLabelParticle-=sTextDY; 
  text = new TText(xLabelParticle,yLabelParticle,note);
  text->SetTextSize(sTextParticles);
  text->SetTextColor(cProton);
  text->Draw();

  sprintf(note,"%.2f",thePart[2]);
  yTextParticle = yLabelParticle;
  text = new TText(xTextParticle,yTextParticle,note);
  text->SetTextAlign(31);
  text->SetTextSize(sTextParticles);
  text->SetTextColor(cProton);
  text->Draw();

  sprintf(note,"#chi^{2} (#pi)");
  yLabelParticle-=sTextDY; 
  tlatex = new TLatex(xLabelParticle,yLabelParticle,note);
  tlatex->SetTextSize(sTextParticles);
  tlatex->SetTextColor(cPion);
  tlatex->Draw();

  sprintf(note,"%.2f",chiPart[0]);
  yTextParticle = yLabelParticle;
  text = new TText(xTextParticle,yTextParticle,note);
  text->SetTextAlign(31);
  text->SetTextSize(sTextParticles);
  text->SetTextColor(cPion);
  text->Draw();

  sprintf(note,"#chi^{2} (K)");
  yLabelParticle-=sTextDY; 
  tlatex = new TLatex(xLabelParticle,yLabelParticle,note);
  tlatex->SetTextSize(sTextParticles);
  tlatex->SetTextColor(cKaon);
  tlatex->Draw();

  sprintf(note,"%.2f",chiPart[1]);
  yTextParticle = yLabelParticle;
  text = new TText(xTextParticle,yTextParticle,note);
  text->SetTextAlign(31);
  text->SetTextSize(sTextParticles);
  text->SetTextColor(cKaon);
  text->Draw();

  sprintf(note,"#chi^{2} (p)");
  yLabelParticle-=sTextDY; 
  tlatex = new TLatex(xLabelParticle,yLabelParticle,note);
  tlatex->SetTextSize(sTextParticles);
  tlatex->SetTextColor(cProton);
  tlatex->Draw();

  sprintf(note,"%.2f",chiPart[2]);
  yTextParticle = yLabelParticle;
  text = new TText(xTextParticle,yTextParticle,note);
  text->SetTextAlign(31);
  text->SetTextSize(sTextParticles);
  text->SetTextColor(cProton);
  text->Draw();

  rdisplay->Modified();
  rdisplay->Update();
  rdisplay->cd();

}


void CsRCDisplay::ClearParticleInfo(){
  nParticles = 0;
  nRings     = 0;
  nRingPhot  = 0;
  nFlag1     = 0;
  nFlag2     = 0;
  nFlag3     = 0;
  nLoopRing  = 0;
  ringMom    =0.;
  ringThe    =0.;
}


void CsRCDisplay::MakeGraphicFile(){

  int32 index = 0;
  char s1[10];
  const char *output =  fGraFileBuf->GetString();
  for(int32 i=0;i<(int32)strlen(output);i++){
    if(output[i] == '.') index = i; 
  }

  static Int_t loop_gif=0;
  static Int_t loop_eps=0;
  static Int_t loop_ps =0;

  std::ostringstream psfile0,psfile1,psfile2,psfile3,psfile4;
  
  char s2[index+1];
  strncpy(s2,output,index);
  s2[index]='\0';
  output+=index;
  strcpy(s1,output);

  if(strcmp(s1,".ps") == 0){
    BGcolor = RWhite;
    cDisplay = RWhite;
    cout << "Change Back Ground Color : white " << endl; 

    psfile0 << s2 << "." << loop_ps << ".ps";
    std::string  tid0 = psfile0.str();
    const char * cid0 = tid0.c_str();

    cout << "PS file : " 
	 << cid0
	 << " has been created." <<  endl; 

    TPostScript myps(cid0,111);

    myps.NewPage();
    myps.Range(20,20); 
    edisplay->SetFillColor(cDisplay);
    edisplay->Modified();
    edisplay->Update();
    edisplay->cd();
    myps.NewPage();

    myps.Range(10,20); 
    hdisplay->SetFillColor(cDisplay);
    hdisplay->cd(1);
    gPad->SetFillColor(cDisplay);
    theHist[2]->Draw();
    theHist[0]->Draw("same");
    hdisplay->cd(2);
    gPad->SetFillColor(cDisplay);
    theHist[3]->Draw();
    theHist[1]->Draw("same");
    hdisplay->Modified();
    hdisplay->Update();
    hdisplay->cd();
    myps.NewPage();

    myps.Range(10,20); 
    phidisplay->SetFillColor(cDisplay);
    phidisplay->Modified();
    phidisplay->Update();
    phidisplay->cd();
    myps.NewPage();

    myps.Range(10,20); 
    idisplay->SetFillColor(cDisplay);
    idisplay->Modified();
    idisplay->Update();
    idisplay->cd();
    myps.NewPage();

    myps.Range(10,20); 
    rdisplay->SetFillColor(cDisplay);
    rdisplay->Modified();
    rdisplay->Update();
    rdisplay->cd();
    myps.Close();

    BGcolor = RBlack;
    cDisplay = RBlack;

    edisplay->cd();

    loop_ps++;

  }

  if(strcmp(s1,".eps") == 0){
    BGcolor = RWhite;
    cDisplay = RWhite;
    edisplay->cd();
    edisplay->Modified();
    edisplay->Update();

    psfile0 << s2 << "." << loop_eps << ".eps";
    std::string  tid0 = psfile0.str();
    const char * cid0 = tid0.c_str();

    TPostScript myeps(cid0,113);
    //    TPostScript myeps((Text_t*)fGraFileBuf->GetString(),113);
    myeps.Range(20,20); 
    edisplay->SetFillColor(cDisplay);
    edisplay->Modified();
    edisplay->Update();
    cout << "EPS file : " 
	 << fGraFileBuf->GetString() 
	 << " has been created." <<  endl; 
    myeps.Close();

    BGcolor = RBlack;
    cDisplay = RBlack;
    edisplay->cd();
    edisplay->Update();

    loop_eps++;
  }

  if(strcmp(s1,".gif") == 0){

    psfile0 << s2 << ".inf." << loop_gif << ".gif";
    psfile1 << s2 << ".r-i." << loop_gif << ".gif";
    psfile2 << s2 << ".the." << loop_gif << ".gif";
    psfile3 << s2 << ".phi." << loop_gif << ".gif";
    psfile4 << s2 << ".rng." << loop_gif << ".gif";

    std::string  tid0 = psfile0.str();
    std::string  tid1 = psfile1.str();
    std::string  tid2 = psfile2.str();
    std::string  tid3 = psfile3.str();
    std::string  tid4 = psfile4.str();
    const char * cid0 = tid0.c_str();
    const char * cid1 = tid1.c_str();
    const char * cid2 = tid2.c_str();
    const char * cid3 = tid3.c_str();
    const char * cid4 = tid4.c_str();
    idisplay->Print(cid0,"gif");
    rdisplay->Print(cid1,"gif");
        gPad->Print(cid2,"gif");
  phidisplay->Print(cid3,"gif");
    edisplay->Print(cid4,"gif");

    psfile0.clear();
    psfile1.clear();
    psfile2.clear();
    psfile3.clear();
    psfile4.clear();

    loop_gif++;

  }
}

void CsRCDisplay::BuildControlPanel(RInit_mode mode){

  //------- COMPASS logo -----------------------------------------------
  fFPIC = new TGCompositeFrame(fF1, 100, 100, kVerticalFrame);
  fF1->AddFrame(fFPIC, new TGLayoutHints(kLHintsLeft | kLHintsTop));

  fPicBut = new TGPictureButton(fFPIC, 
				gClient->GetPicture("${CORAL}/src/evdis/compass.xpm"), 3);
  fPicBut->SetCommand("printf(\"Hello world!\\n\");");
  fFPIC->AddFrame(fPicBut,
		  new TGLayoutHints(kLHintsLeft | kLHintsTop,10,5,5,5));
  //---------------------------------------------------------------------

  //------------------  Controll Window Frame ---------------------------
  fF2 = new TGCompositeFrame(fF1, 700, 100, kHorizontalFrame);
  fF1->AddFrame(fF2, new TGLayoutHints(kLHintsLeft | kLHintsTop | kLHintsExpandX | kLHintsExpandY));

  //------------------- TAB ---------------------------------------------
  fTab = new TGTab(fF2,500,150);
  fF2->AddFrame(fTab, new TGLayoutHints(kLHintsLeft | kLHintsTop | kLHintsExpandX | kLHintsExpandY));
  
  Tf0 = fTab->AddTab("General"); 
  fT1 = new TGGroupFrame(Tf0, "General Control Panel", kHorizontalFrame);
  Tf0->AddFrame(fT1,new TGLayoutHints(kLHintsTop | kLHintsLeft | kLHintsExpandX | kLHintsExpandY,5,2,5,0)); 
  //----------------------------------------------------------------------
  fT1_1 = new TGGroupFrame(fT1, "Coral Control", kVerticalFrame);
  fT1->AddFrame(fT1_1,new TGLayoutHints(kLHintsTop | kLHintsLeft | kLHintsExpandX | kLHintsExpandY,5,0,5,0)); 

  fExitButton = new TGTextButton(fT1_1,"    CLOSE   ", RF_QUIT); 
  fExitButton->Associate(this);
  fT1_1->AddFrame(fExitButton,
		  new TGLayoutHints(kLHintsLeft | kLHintsTop | kLHintsExpandX | kLHintsExpandY,10,5,5,5));  
 
  fNextEventButton = new TGTextButton(fT1_1," Next  Event ", RF_NEXT); 
  fNextEventButton->Associate(this);
  fT1_1->AddFrame(fNextEventButton,
		  new TGLayoutHints(kLHintsLeft | kLHintsTop | kLHintsExpandX | kLHintsExpandY,10,5,5,5)); 
  
  if(mode == RStandalone){
    fPreviusEventButton = new TGTextButton(fT1_1,
					   " Previus Event ", RF_PREVIUS); 
    fPreviusEventButton->Associate(this);
    fT1_1->AddFrame(fPreviusEventButton,
		    new TGLayoutHints(kLHintsLeft | kLHintsTop | kLHintsExpandX | kLHintsExpandY,10,5,5,5)); 
  }
  //----------------------------------------------------------------------
  fT1_2 = new TGGroupFrame(fT1, "Ring Control Panel", kVerticalFrame);
  fT1_2->SetLayoutManager(new TGMatrixLayout(fT1_2,2,2,10)); 
  fT1->AddFrame(fT1_2,new TGLayoutHints(kLHintsTop | kLHintsLeft | kLHintsExpandX | kLHintsExpandY,5,0,5,0)); 

  fNextRingButton = new TGTextButton(fT1_2,"Next Ring",RF_NEXT_RING); 
  fNextRingButton->Associate(this);
  fT1_2->AddFrame(fNextRingButton,
		  new TGLayoutHints(kLHintsLeft | kLHintsTop | kLHintsExpandX | kLHintsExpandY,5,5,5,5)); 
  fPreviusRingButton = new TGTextButton(fT1_2,"Previus Ring",RF_PREVIUS_RING); 
  fPreviusRingButton->Associate(this);
  fT1_2->AddFrame(fPreviusRingButton,
		  new TGLayoutHints(kLHintsLeft | kLHintsTop | kLHintsExpandX | kLHintsExpandY,5,5,5,5)); 

  fGraphicOutButton = new TGTextButton(fT1_2,
				       "Graphic Out ",RF_GRAPHIC_OUT); 
  fGraphicOutButton->Associate(this);
  fT1_2->AddFrame(fGraphicOutButton,
		  new TGLayoutHints(kLHintsLeft | kLHintsTop | kLHintsExpandX | kLHintsExpandY,10,5,5,5)); 

  fGraFileBuf = new TGTextBuffer(100);
  fGraFileBuf->AddText(0,"evdis.ps");
  fGraFileName = new TGTextEntry(fT1_2,fGraFileBuf, -1);
  fT1_2->AddFrame(fGraFileName, new TGLayoutHints(kLHintsTop | kLHintsLeft | kLHintsExpandX | kLHintsExpandY, 10, 2, 2, 2));
  fGraFileName->Resize(100, fGraFileName->GetDefaultHeight());

  //----------------------------------------------------------------------

  fT1_3 = new TGGroupFrame(fT1, "Background Color", kVerticalFrame);
  fT1_3->SetLayoutManager(new TGMatrixLayout(fT1_3,1,0,10)); 
  fT1->AddFrame(fT1_3,new TGLayoutHints(kLHintsTop | kLHintsLeft | kLHintsExpandX | kLHintsExpandY,5,0,5,0)); 

  fBackColor = new TGTextButton(fT1_3,"  Black or White   ",RF_COLOR_BACK);
  fBackColor->Associate(this);
  fT1_3->AddFrame(fBackColor,
		  new TGLayoutHints(kLHintsTop | kLHintsLeft | kLHintsExpandX | kLHintsExpandY, 5, 5, 5, 5));

  //----------------------------------------------------------------------

  fT1_4 = new TGGroupFrame(fT1, "Select Infos", kVerticalFrame);
  fT1_4->SetLayoutManager(new TGMatrixLayout(fT1_4,2,2,10)); 
  fT1->AddFrame(fT1_4,new TGLayoutHints(kLHintsTop | kLHintsLeft | kLHintsExpandX | kLHintsExpandY,5,0,5,0)); 

  fTrackOnOff = new TGTextButton(fT1_4,"  TRACK   ",RF_TRACK_ON_OFF);
  fTrackOnOff->Associate(this);
  fT1_4->AddFrame(fTrackOnOff,
		new TGLayoutHints(kLHintsTop | kLHintsLeft | kLHintsExpandX | kLHintsExpandY, 5, 5, 5, 5));

  fDigitOnOff = new TGTextButton(fT1_4,"  RINGS   ",RF_DIGIT_ON_OFF);
  fDigitOnOff->Associate(this);
  fT1_4->AddFrame(fDigitOnOff,
		new TGLayoutHints(kLHintsTop | kLHintsLeft | kLHintsExpandX | kLHintsExpandY, 5, 5, 5, 5));

  fHitOnOff = new TGTextButton(fT1_4,"   HIT    ",RF_HIT_ON_OFF);
  fHitOnOff->Associate(this);
  fT1_4->AddFrame(fHitOnOff,
		new TGLayoutHints(kLHintsTop | kLHintsLeft | kLHintsExpandX | kLHintsExpandY, 5, 5, 5, 5));

  fClusterOnOff = new TGTextButton(fT1_4,"   CLUSTER    ",RF_CLUSTER_ON_OFF);
  fClusterOnOff->Associate(this);
  fT1_4->AddFrame(fClusterOnOff,
		new TGLayoutHints(kLHintsTop | kLHintsLeft | kLHintsExpandX | kLHintsExpandY, 5, 5, 5, 5));

  //-------------------- View tab -------------------------------------
  Tf0 = fTab->AddTab("View"); 

  fT2 = new TGGroupFrame(Tf0, "View Control Panel", kHorizontalFrame);
  Tf0->AddFrame(fT2,new TGLayoutHints(kLHintsTop | kLHintsLeft | kLHintsExpandX | kLHintsExpandY,5,0,5,0)); 

  fT21 = new TGGroupFrame(fT2, "View Direction", kVerticalFrame);
  fT21->SetLayoutManager(new TGMatrixLayout(fT21,3,0,10)); 
  fT2->AddFrame(fT21,new TGLayoutHints(kLHintsTop | kLHintsLeft | kLHintsExpandX | kLHintsExpandY,5,0,5,0)); 

  fTopView = new TGTextButton(fT21,"    TOP  VIEW    ",RF_TOP_VIEW);
  fTopView->Associate(this);
  fT21->AddFrame(fTopView,
		 new TGLayoutHints(kLHintsTop | kLHintsLeft | kLHintsExpandX | kLHintsExpandY, 5, 5, 5, 5)); 

  fSideView = new TGTextButton(fT21,"    SIDE VIEW    ",RF_SIDE_VIEW);
  fSideView->Associate(this);
  fT21->AddFrame(fSideView,
		 new TGLayoutHints(kLHintsTop | kLHintsLeft | kLHintsExpandX | kLHintsExpandY, 5, 5, 5, 5)); 

  fFrontView = new TGTextButton(fT21,"    FRONT VIEW    ",RF_FRONT_VIEW);
  fFrontView->Associate(this);
  fT21->AddFrame(fFrontView,
		 new TGLayoutHints(kLHintsTop | kLHintsLeft | kLHintsExpandX | kLHintsExpandY, 5, 5, 5, 5)); 


  fT22 = new TGGroupFrame(fT2, "Zooming", kVerticalFrame);
  fT22->SetLayoutManager(new TGMatrixLayout(fT22,2,0,10)); 
  fT2->AddFrame(fT22,new TGLayoutHints(kLHintsTop | kLHintsLeft | kLHintsExpandX | kLHintsExpandY,5,0,5,0)); 


  fZoomInButton = new TGTextButton(fT22, "    ZOOM  IN     ", RF_ZOOM_IN);
  fZoomInButton->Associate(this);
  fT22->AddFrame(fZoomInButton, 
		 new TGLayoutHints(kLHintsTop | kLHintsLeft | kLHintsExpandX | kLHintsExpandY, 5, 5, 5, 5));  

  fZoomOutButton = new TGTextButton(fT22, "    ZOOM OUT   ", RF_ZOOM_OUT);
  fZoomOutButton->Associate(this);
  fT22->AddFrame(fZoomOutButton, 
		 new TGLayoutHints(kLHintsTop | kLHintsLeft | kLHintsExpandX | kLHintsExpandY, 5, 5, 5, 5));  


  fT23 = new TGGroupFrame(fT2, "Explore Buttons", kHorizontalFrame);
  fT23->SetLayoutManager(new TGMatrixLayout(fT23,2,2,10,10));
  fT2->AddFrame(fT23, 
		new TGLayoutHints(kLHintsTop | kLHintsLeft | kLHintsExpandX | kLHintsExpandY, 5,5,5,5));
 

  fRegionUpButton = new TGTextButton(fT23, "     UP      ", RF_EXPLOR_UP);
  fRegionUpButton->Associate(this);
  fT23->AddFrame(fRegionUpButton, 
		 new TGLayoutHints(kLHintsTop | kLHintsLeft | kLHintsExpandX | kLHintsExpandY, 5, 5, 5, 5));  

  fRegionDownButton = new TGTextButton(fT23, "   DOWN   ",RF_EXPLOR_DOWN);
  fRegionDownButton->Associate(this);
  fT23->AddFrame(fRegionDownButton, 
		 new TGLayoutHints(kLHintsTop | kLHintsLeft | kLHintsExpandX | kLHintsExpandY, 5, 5, 5, 5));  

  fRegionLeftButton = new TGTextButton(fT23, "    LEFT    ", RF_EXPLOR_LEFT);
  fRegionLeftButton->Associate(this);
  fT23->AddFrame(fRegionLeftButton, 
		 new TGLayoutHints(kLHintsTop | kLHintsLeft | kLHintsExpandX | kLHintsExpandY, 5, 5, 5, 5));  

  fRegionRightButton = new TGTextButton(fT23,"   RIGHT    ",RF_EXPLOR_RIGHT);
  fRegionRightButton->Associate(this);
  fT23->AddFrame(fRegionRightButton, 
		 new TGLayoutHints(kLHintsTop | kLHintsLeft | kLHintsExpandX | kLHintsExpandY, 5, 5, 5, 5));  

  fT24 = new TGGroupFrame(fT2, "Scale of Z",kVerticalFrame);
  fT24->SetLayoutManager(new TGMatrixLayout(fT24,2,0,10,10));
  fT2->AddFrame(fT24, 
		new TGLayoutHints(kLHintsTop | kLHintsLeft | kLHintsExpandX | kLHintsExpandY, 5,5,5,5));
 
  fScaleIncButton = new TGTextButton(fT24, "    INCREASE    ", RF_SCALE_INC);
  fScaleIncButton->Associate(this);
  fT24->AddFrame(fScaleIncButton, 
		 new TGLayoutHints(kLHintsTop | kLHintsLeft | kLHintsExpandX | kLHintsExpandY, 5, 5, 5, 5));  

  fScaleDecButton = new TGTextButton(fT24,"   DECREASE    ",RF_SCALE_DEC);
  fScaleDecButton->Associate(this);
  fT24->AddFrame(fScaleDecButton, 
		 new TGLayoutHints(kLHintsTop | kLHintsLeft | kLHintsExpandX | kLHintsExpandY, 5, 5, 5, 5));  

  // Now not implemented. 
  fX3D = new TGTextButton(fT2, "  X3D  ", RF_X3D);
  fX3D->Associate(this);
  fT2->AddFrame(fX3D, new TGLayoutHints(kLHintsTop | kLHintsLeft | kLHintsExpandX | kLHintsExpandY,
  					5, 5, 5, 5));	

  fResetViewButton = new TGTextButton(fT2, "   RESET VIEW  ", RF_RESET_VIEW);
  fResetViewButton->Associate(this);
  fT2->AddFrame(fResetViewButton, 
		new TGLayoutHints(kLHintsTop | kLHintsLeft | kLHintsExpandX | kLHintsExpandY, 5, 5, 5, 5));  

}

void CsRCDisplay::BackColor(){
  switch(BGcolor){
  case RBlack:
    BGcolor = RWhite;
    cDisplay = BGcolor;
    cout << "Change Back Ground Color : white " << endl; 
    Redraw();
    break;
  case RWhite:
    BGcolor = RBlack;
    cDisplay = BGcolor;
    cout << "Change Back Ground Color : black" << endl; 
    Redraw();
    break;
  defult:
    break;
  }
}

void CsRCDisplay::X3D(){
  edisplay->cd();
  TView *view = edisplay->GetView();
  if (!view) return; edisplay->x3d("OPENGL");
}


void CsRCDisplay::ReadEvent(){

  event->getNextEvent();
  SetRunNumber(event->getRunNumber());  
  SetEventNumberInRun(event->getEventNumberInRun());  
  SetEventNumberInBurst(event->getEventNumberInBurst());  

}


void CsRCDisplay::ScaleInc(){
  Float_t Cmin,Cmax,Center; 
  Phi = view->GetLongitude();
  Theta = view->GetLatitude();
  view->GetRange(min,max);
  Psi = view->GetPsi();
  Center = (min[2]+max[2])/2.;
  Cmin = ZoomFactor*TMath::Abs(Center - min[2]);
  Cmax = ZoomFactor*TMath::Abs(Center - max[2]);
  min[2] = Center - Cmin;
  max[2] = Center + Cmax;
 
  ClearViewPad();
  SetRange(min[0],min[1],min[2],max[0],max[1],max[2]);
  SetView(Phi,Theta,Psi);
  edisplay->Modified();
  edisplay->Update();
  edisplay->cd();
}

void CsRCDisplay::ScaleDec(){
  Float_t Cmin,Cmax,Center;
  Phi = view->GetLongitude();
  Theta = view->GetLatitude();
  Psi = view->GetPsi();
  view->GetRange(min,max);
  Center = (min[2]+max[2])/2.;
  Cmin = 1. / ZoomFactor*TMath::Abs(Center - min[2]);
  Cmax = 1. / ZoomFactor*TMath::Abs(Center - max[2]);
  min[2] = Center - Cmin;
  max[2] = Center + Cmax;
  ClearViewPad();
  SetRange(min[0],min[1],min[2],max[0],max[1],max[2]);
  SetView(Phi,Theta,Psi);
  edisplay->Modified();
  edisplay->Update();
  edisplay->cd();
}


void CsRCDisplay::RingOn(){ 
  drawRings = kTRUE;
  Redraw();
}

void CsRCDisplay::RingOff(){ 
  drawRings = kFALSE;
  Redraw();
}

void CsRCDisplay::HitOn(){ 
  drawHits = kTRUE;
  Redraw();
}

void CsRCDisplay::HitOff(){ 
  drawHits = kFALSE;
  Redraw();
}

void CsRCDisplay::ClusterOn(){ 
  drawClusters = kTRUE;
  Redraw();
}

void CsRCDisplay::ClusterOff(){ 
  drawClusters = kFALSE;
  Redraw();
}

void CsRCDisplay::TrackOn(){ 
  drawTracks = kTRUE;
  Redraw();
}

void CsRCDisplay::TrackOff(){ 
  drawTracks = kFALSE;
  Redraw();
}
