#include "PlaneDriftChamberPanel.h"
#include <TMath.h>
#include <TPad.h>
#include <TF1.h>
#include <TGFileDialog.h>
#include <TGMsgBox.h>

#include <stdio.h>
#include <iostream>

ClassImp(PlaneDriftChamberPanel);
ClassImp(FitData);

//_________________________________________________________
PlaneDriftChamberPanel::PlaneDriftChamberPanel(
  const TGWindow *p,const TGWindow *main,
  UInt_t w, UInt_t h,
  PlaneDriftChamber *plane):
  Plane1VPanel( p, main, w, h, plane),
  fHtvsch(0)
{
   CreateControlTab();
   fFitDatas.clear();
}

//_________________________________________________________
PlaneDriftChamberPanel::PlaneDriftChamberPanel(
  const TGWindow *p,const TGWindow *main,
  UInt_t w, UInt_t h,
  PlaneDriftChamber2V *plane):
  Plane1VPanel( p, main, w, h, plane),
  fHtvsch(0)
{
   CreateControlTab();
   fFitDatas.clear();
}

//_________________________________________________________
void PlaneDriftChamberPanel::CreateControlTab()
{
//  std::cout << "PlaneDriftChamberPanel::CreateControlTab.\n";

  //control tab
  TGCompositeFrame *tf=fTab->AddTab("Time Fit");
  TGCompositeFrame *hFrame=new TGCompositeFrame(tf, 60, 60, kHorizontalFrame);

  //=== Root canvas ===
  TGGroupFrame *DpyFrame=new TGGroupFrame(hFrame, "Display", kVerticalFrame);
  fTimeDpy = new TRootEmbeddedCanvas("chDpy",DpyFrame, 300, 300);
  DpyFrame->AddFrame(fTimeDpy);

  //=== Parameter Group frame ===
  TGGroupFrame *ParFrame=new TGGroupFrame(hFrame, "Parameters", kVerticalFrame);
  TGHorizontalFrame *LbFrame = new TGHorizontalFrame( ParFrame, 60, 60 );
  LbFrame->SetLayoutManager(new TGMatrixLayout(LbFrame,0,2));

  //=== channel step
  TGLabel *chstL = new TGLabel(LbFrame,"Channel step ");
  fChStepTE = new TGTextEntry(LbFrame,"8      ");
  LbFrame->AddFrame(chstL);
  LbFrame->AddFrame(fChStepTE);

  Variable* v = fPlane->GetVariable("_t_on_trigger");
  char* val = new char[16];

  //=== TMin
  TGLabel *tMinL = new TGLabel(LbFrame,"       T Min ");
  sprintf( val, "%7.0f", ( (v) ? v->GetMin():0 ) ); fTMinTE = new TGTextEntry(LbFrame, val);
  LbFrame->AddFrame(tMinL);
  LbFrame->AddFrame(fTMinTE);

  //=== TMax
  TGLabel *tMaxL = new TGLabel(LbFrame,"       T Max ");
  sprintf( val, "%7.0f", ( (v) ? v->GetMax():0 ) ); fTMaxTE = new TGTextEntry(LbFrame, val);
  LbFrame->AddFrame(tMaxL);
  LbFrame->AddFrame(fTMaxTE);
  delete val;

  //=== Buttons
  TGCompositeFrame *BtFrame = new TGCompositeFrame( ParFrame, 60, 60, kVerticalFrame );

  fFitBt = new TGTextButton(BtFrame, "  &Fit  ",0);
  fFitBt->Connect("Clicked()","PlaneDriftChamberPanel",this,"Fit()");
  BtFrame->AddFrame(fFitBt, new TGLayoutHints(  kLHintsExpandX ));

  fWriteBt = new TGTextButton(BtFrame, " &Write to File ",0);
  fWriteBt->Connect("Clicked()","PlaneDriftChamberPanel",this,"WriteDBFile()");
  BtFrame->AddFrame(fWriteBt, new TGLayoutHints(  kLHintsExpandX ) );

  //=== Put everything together

  ParFrame->AddFrame( LbFrame, new TGLayoutHints(kLHintsTop |
						      kLHintsLeft |
                  kLHintsExpandX,
						      2, 2, 2, 2));
  ParFrame->AddFrame( BtFrame, new TGLayoutHints(kLHintsTop |
						      kLHintsLeft |
						      kLHintsExpandX,
						      2, 2, 2, 2));

  hFrame->AddFrame(ParFrame,new TGLayoutHints(kLHintsTop |
						      kLHintsLeft,
						      2, 2, 2, 2));

  hFrame->AddFrame(DpyFrame,new TGLayoutHints(kLHintsTop |
						      kLHintsLeft |
						      kLHintsExpandY |
                  kLHintsExpandX,
						      2, 2, 2, 2));

  tf->AddFrame(hFrame);


  fMain->MapSubwindows();
  fMain->Resize(fMain->GetDefaultSize());
  fMain->Layout();

  fHtvsch = (TH2F*) fPlane->GetHisto("_tVSch");
  if( fHtvsch ){
    fHtvsch->ProjectionY()->Draw();
    TPad::Pad()->Modified();
    TPad::Pad()->Update();
  }

}


//_________________________________________________________
void PlaneDriftChamberPanel::DeleteFitDatas( void )
{
  while( !fFitDatas.empty() ) {
    FitData *fd = fFitDatas.back();
    fFitDatas.pop_back();
    delete fd;
  }

  return;

}

//_________________________________________________________
PlaneDriftChamberPanel::~PlaneDriftChamberPanel()
{
  //=== delete FitDatas
  DeleteFitDatas();

  //=== delete buttons ===
  if( fFitBt )    delete fFitBt;
  if( fWriteBt )  delete fWriteBt;

  //=== delete Text Entries ===
  if( fChStepTE ) delete fChStepTE;
  if( fTMinTE ) delete fTMinTE;
  if( fTMaxTE ) delete fTMaxTE;

  //=== delete ROOT canvas ===
  if( fTimeDpy ) delete fTimeDpy;


}

//_________________________________________________________
void PlaneDriftChamberPanel::Fit( int id )
{

  //=== Check objects
  if( fTMinTE == 0 || fTMaxTE == 0 || fChStepTE == 0) {
    std::cout << "PlaneDriftChamberPanel::Fit - ERROR: cannot read TextEntries.\n";
    return;
  }

  //=== Check histogram
  if( !fHtvsch ) {
    std::cout << "PlaneDriftChamberPanel::Fit - ERROR: cannot find histogram.\n";
    return;
  }


  //=== Get Fit Range
  float tMin; sscanf( fTMinTE->GetText(),   "%f", &tMin );
  float tMax; sscanf( fTMaxTE->GetText(),   "%f", &tMax );
  int chStep; sscanf( fChStepTE->GetText(), "%i", &chStep );

  //=== check time limits
  if( tMin >= tMax ) {
    TGMsgBox( gClient->GetRoot(), fMain, "Error", "ERROR: wrong time limits.", 0 );
    return;
  }

  //=== create function
  TF1 *f = new TF1("TimeFitFunction", TimeFitFunction , tMin, tMax, 4);
  std::cout << "PlaneDriftChamberPanel::TimeFit - INFO:  TMin="<< tMin << "  TMax=" << tMax << " ChStep=" << chStep << std::endl;

  //=== Delete datas from previous fit
  DeleteFitDatas();

  for( int i = 0; i < fHtvsch->GetNbinsX() ; i += chStep ) {

    //=== Get and draw the slice histogram
    TH1* h = fHtvsch->ProjectionY( "_py",i, i+chStep-1 );
    h->Draw();
    if( !h->GetEntries() ) continue;

    //=== Set TimeFitFunction parameters
    double p0 = double( h->GetBinContent( h->GetXaxis()->FindBin( tMin ) ) );
    double p1 = double( h->GetMaximum() - p0 );
    double p2 = 0.5*(tMax+tMin);
    double p3 = double(tMax-tMin);

    f->SetParameter( 0, p0 );
    f->SetParameter( 1, p1 );
    f->SetParameter( 2, p2 );
    f->SetParameter( 3, p3 );
    f->SetLineWidth( 4 );
    f->SetLineColor( 2 );

    //=== Do the fit. Plot the result
    h->Fit(f,"0Q","", tMin, tMax );
    f->Draw("same");

    TPad::Pad()->Modified();
    TPad::Pad()->Update();

    //=== Put the result in fitData list
    FitData * fd = new FitData( i, i+chStep-1 );
    fd->entries = (unsigned int) h->GetEntries();
    for( unsigned int ip = 0; ip< 4; ip++ )
    fd->pars.push_back( (double) f->GetParameter( ip ) );
    fd->chi2 = (double) f->GetChisquare();
    fd->ndf =  (unsigned int) f->GetNDF();
    fFitDatas.push_back( fd );
  }

  TGMsgBox( gClient->GetRoot(), fMain, "Info", "Done.", 0 );
  return;
}

//_________________________________________________________
const char* PlaneDriftChamberPanel::DBFileTypes[] = {
  "DB files",   "*.DB",
  "All files",  "*.*",
  0, 0 };

//_________________________________________________________
void PlaneDriftChamberPanel::WriteDBFile( int id )
{

  //=== check fit Datas exist
  if( fFitDatas.empty() ) {
    TGMsgBox( gClient->GetRoot(), fMain, "Error", "ERROR: no fit data.", 0 );
    return;
  }

  //=== Get Filename from dialog box
  //=== Generate default file Name
  TGFileInfo fi;
  fi.fFileTypes =  CHAR_TGFILEINFO(DBFileTypes);
  new TGFileDialog( gClient->GetRoot(), fMain, kFDSave, &fi);

  if( !fi.fFilename ) return;

//   std::ofstream out( fi.fFilename, std::ios::out );
  FILE* out = fopen( fi.fFilename, "r" );
  if( !out ) {
    std::string msg = std::string("ERROR: cannot write to file \"") + std::string(fi.fFilename)+"\".";
    TGMsgBox( gClient->GetRoot(), fMain, "Error", msg.c_str(), 0 );
    delete fi.fFilename;
    return;
  }

  delete fi.fFilename;

  for( std::list<FitData*>::iterator i = fFitDatas.begin(); i != fFitDatas.end(); i++ ) {
//     out.form( "%4i %4i  %7i %15.8g %15.8g %15.8g %15.8g %15.8g %6i \n",
//       (*i)->iMin, (*i)->iMax,
//       (*i)->entries,
//       (*i)->pars[0], (*i)->pars[1], (*i)->pars[2], (*i)->pars[3],
//       (*i)->chi2,
//       (*i)->ndf);
    fprintf(out, "%4i %4i  %7i %15.8g %15.8g %15.8g %15.8g %15.8g %6i \n",
      (*i)->iMin, (*i)->iMax,
      (*i)->entries,
      (*i)->pars[0], (*i)->pars[1], (*i)->pars[2], (*i)->pars[3],
      (*i)->chi2,
      (*i)->ndf);
  }

//   out.form( "\n%4s %4s  %7s %15s %15s %15s %15s %15s %6s \n",
//     "iMin", "iMax",
//     "entries",
//     "bg", "sig", "t0", "slope",
//     "chi2",
//     "ndf");
  fprintf(out, "\n%4s %4s  %7s %15s %15s %15s %15s %15s %6s \n",
    "iMin", "iMax",
    "entries",
    "bg", "sig", "t0", "slope",
    "chi2",
    "ndf");

//   out.close();
  fclose(out);
  return;
}

//_______________________________________________________________________________
double TimeFitFunction(double *x,double *par)
{
  double t = x[0];
  double N = par[0]+par[1]*TMath::Erf((t-par[2])/par[3]);
  return N;
}
