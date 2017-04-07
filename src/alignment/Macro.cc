// $Id: Macro.cc,v 1.21 2006/06/16 15:21:47 conrad Exp $

/*!
   \file    Macro.cc
   \brief   some root macro to make predefined plots, fits, etc.
   \author  Hugo Pereira
   \version $Revision: 1.21 $
   \date    $Date: 2006/06/16 15:21:47 $
*/

#include "Macro.h"
#include <TH1.h>
#include <TF1.h>
#include <TH2.h>
#include <TFile.h>
#include <TCanvas.h>
#include <TPad.h>
#include <TStyle.h>
#include <TPaveText.h>
#include <TPaveStats.h>
#include <TStyle.h>
#include <vector>
#include <string>

#include "DetFileManager.h"
#include "DetectorInfo.h"
#include "DeadZoneInfo.h"
#include "Fit.h"
#include "Utils.h"

using namespace std;

//________________________________________
ClassImp(Macro)

//________________________________________
void Macro::SetMBTiSli( const char* fileselection )
{
  vector<string> files = Utils::GetFiles( fileselection );
  for( unsigned int i=0; i<files.size(); i++ ) {
    DetFileManager df( files[i].c_str() );
    vector<DetectorInfo*> dets = df.GetDetSelection("MB");
    for( unsigned int j = 0; j < dets.size(); j++ ) 
    for( unsigned int k = 0; k < dets[j]->sub_.size(); k++ ) { 
      dets[j]->sub_[k]->vel_   = 0.0075;
      dets[j]->sub_[k]->t0_    = -1245.0;
      dets[j]->sub_[k]->thRes_ = 900.0;
      dets[j]->sub_[k]->spSli_ = 260.0;
      dets[j]->sub_[k]->tiSli_ = 0.12892;
    }
    df.DumpToFile();
  }
  return;
}
    
  
//________________________________________
void Macro::SetGemDeadZoneSize( const char* fileselection )
{
  vector<string> files = Utils::GetFiles( fileselection );
  for( unsigned int i=0; i<files.size(); i++ ) {
    DetFileManager df( files[i].c_str() );
    vector<DetectorInfo*> dets = df.GetDetSelection("GM");
    for( unsigned int j = 0; j< dets.size(); j++ ) {
      DeadZoneInfo* dead = dets[j]->GetMain()->dead_;
      if( !dead ) continue; 
      dead->xsiz_ = 50;
      dead->ysiz_ = 0.05;
    }
    df.DumpToFile();
  }
  return;
}

//________________________________________
void Macro::SetDCDeadZoneSize( const char* fileselection )
{
  vector<string> files = Utils::GetFiles( fileselection );
  for( unsigned int i=0; i<files.size(); i++ ) {
    DetFileManager df( files[i].c_str() );
    vector<DetectorInfo*> dets = df.GetDetSelection("DC01");
    for( unsigned int j = 0; j< dets.size(); j++ ) {
      DeadZoneInfo* dead = dets[j]->GetMain()->dead_;
      if( dead ) dead->xsiz_ = 300;
    }
    df.DumpToFile();
  }
  return;
}

//________________________________________
void Macro::ReorderGM02( const char* fileselection )
{
  vector<string> files = Utils::GetFiles( fileselection );
  for( unsigned int i=0; i<files.size(); i++ ) {
    DetFileManager df( files[i].c_str() );
    DetectorInfo* det;
    double newZ;
    det = df.GetDetInfo( "GM02U" );
    if( det ) {
      newZ = 5185.592;
      det->zcm_ = newZ;
      if( det->dead_ ) det->dead_->zcm_ = newZ;
    }

    det = df.GetDetInfo( "GM02V" );
    if( det ) {
      newZ = 5185.642;
      det->zcm_ = newZ;
      if( det->dead_ ) det->dead_->zcm_ = newZ;
    }

    det = df.GetDetInfo( "GM02Y" );
    if( det ) {
      newZ = 5209.222;
      det->zcm_ = newZ;
      if( det->dead_ ) det->dead_->zcm_ = newZ;
    }

    det = df.GetDetInfo( "GM02X" );
    if( det ) {
      newZ = 5209.272;
      det->zcm_ = newZ;
      if( det->dead_ ) det->dead_->zcm_ = newZ;
    }
    
    // make all deadZones z match detector z
    vector< DetectorInfo* > dets = df.GetDetSelection( "*" );
    for( unsigned int i=0; i< dets.size(); i++ )
    if( dets[i]->dead_ ) dets[i]->dead_->zcm_ = dets[i]->zcm_;
    
    df.Sort( "GM" );
    df.DumpToFile();
  }
  
  
  return;
}
  
//________________________________________
TCanvas* Macro::HaloFI( const char* fileselection )
{
  TCanvas *cv = new TCanvas("HaloFI","HaloFI", 700, 700 );
  gROOT->SetStyle("Plain");
  gStyle->SetOptFit(1);
  cv->Divide( 2, 3, 0.005, 0.005, 0 );
  const char* h_names[] = {
    "Traffic/RDmonitor/FI01X1  _res1",
    "Traffic/RDmonitor/FI01Y1  _res1",
    "Traffic/RDmonitor/FI15X1  _res1",
    "Traffic/RDmonitor/FI15Y1  _res1",
    "Traffic/RDmonitor/FI02X1  _res1",
    "Traffic/RDmonitor/FI02Y1  _res1",
    0 };
  
  for( unsigned int i=0; h_names[i]; i++ ) {
    if( !h_names[i] ) break;
    TH2* h2 = Utils::Get2D( fileselection, h_names[i] );
    if( !h2 ) continue;
    TH1* h1 = h2->ProjectionY();
    cv->cd( i+1 );
    h1->Draw();
    
    double min = h1->GetXaxis()->GetXmin();
    double max = h1->GetXaxis()->GetXmax();
    TF1* f = new TF1("GP0", Fit::GP0, min, max, 4);
    f->SetLineColor(4);
    f->SetLineWidth(1);
    f->SetParameter( 0, 0.5*h1->GetMaximum() ) ;
    f->SetParameter( 1, h1->GetBinCenter( h1->GetMaximumBin() ) );
    f->SetParameter( 2, 0.5*h1->GetRMS( ) );
    f->SetParameter( 3, 0 );
    h1->Fit( f, "0" ); 
    f->Draw("same");
    
    cv->Update();
  }
  
  return cv;
}
 
//________________________________________
TCanvas* Macro::HaloSI( const char* fileselection )
{
  TCanvas *cv = new TCanvas("HaloSI","HaloSI", 700, 700 );
  gROOT->SetStyle("Plain");
  gStyle->SetOptFit(1);
  cv->Divide( 4, 2, 0.005, 0.005, 0 );
  const char* h_names[] = {
    "Traffic/RDmonitor/SI01U1  _res1",
    "Traffic/RDmonitor/SI01V1  _res1",
    "Traffic/RDmonitor/SI01X1  _res1",
    "Traffic/RDmonitor/SI01Y1  _res1",
    "Traffic/RDmonitor/SI02U1  _res1",
    "Traffic/RDmonitor/SI02V1  _res1",
    "Traffic/RDmonitor/SI02X1  _res1",
    "Traffic/RDmonitor/SI02Y1  _res1",
    0 };
  
  for( unsigned int i=0; h_names[i]; i++ ) {
    if( !h_names[i] ) break;
    TH2* h2 = Utils::Get2D( fileselection, h_names[i] );
    if( !h2 ) continue;
    TH1* h1 = h2->ProjectionY();
    cv->cd( i+1 );
    h1->Draw();
    
    double min = h1->GetXaxis()->GetXmin();
    double max = h1->GetXaxis()->GetXmax();
    TF1* f = new TF1("GP0", Fit::GP0, min, max, 4);
    f->SetLineColor(4);
    f->SetLineWidth(1);
    f->SetParameter( 0, 0.5*h1->GetMaximum() ) ;
    f->SetParameter( 1, h1->GetBinCenter( h1->GetMaximumBin() ) );
    f->SetParameter( 2, 0.5*h1->GetRMS( ) );
    f->SetParameter( 3, 0 );
    h1->Fit( f, "0" ); 
    f->Draw("same");
    
    cv->Update();
  }
  
  return cv;
}

#ifndef __CINT__   
//________________________________________
TH2* Macro::Chi2( vector< string > files )
{
  TH2* histo = Utils::Get2D(files,"Traffic/RDmonitor/tr_11" );
  if( !histo ) return 0;
  
  // SLICES 
  // Set name, title, line color, entries (for entries are not properly
  // handled by ->Projection())
  gROOT->SetStyle("Plain");
  int N1 = histo->GetNbinsX()+1;
  
  TH1D *p0x3 = new TH1D(*(TH1D*)histo->ProjectionX("_0x3",4,4));
  p0x3->SetLineColor(2);
  p0x3->SetName("0x3"); 
  p0x3->SetTitle("0x3");
  p0x3->SetEntries(histo->Integral(0,N1,4,4));
  
  TH1D *p0x7 = new TH1D(*(TH1D*)histo->ProjectionX("_0x7",8,8));
  p0x7->SetLineColor(4);
  p0x7->SetName("0x7"); 
  p0x7->SetTitle("0x7");
  p0x7->SetEntries(histo->Integral(0,N1,8,8));
  
  TH1D *p0xf = new TH1D(*(TH1D*)histo->ProjectionX("_0xf",16,16));
  p0xf->SetLineColor(6);
  p0xf->SetName("0xf"); 
  p0xf->SetTitle("0xf");
  p0xf->SetEntries(histo->Integral(0,N1,16,16));
  
  TH1D *p0x6 = new TH1D(*(TH1D*)histo->ProjectionX("_0x6",7,7));
  p0x6->SetLineColor(2);
  p0x6->SetName("0x6"); 
  p0x6->SetTitle("0x6");
  p0x6->SetEntries(histo->Integral(0,N1,7,7));
  
  TH1D *p0xe = new TH1D(*(TH1D*)histo->ProjectionX("_0xe",15,15));
  p0xe->SetLineColor(3);
  p0xe->SetName("0xe"); 
  p0xe->SetTitle("0xe");
  p0xe->SetEntries(histo->Integral(0,N1,15,15));

  // 0x3 OR 
  TH1D *all0x3 = new TH1D(*p0x3);
  all0x3->Add(p0x7); all0x3->Add(p0xf);
  int len = strlen(histo->GetTitle()) + 1 + strlen(" &0x3");
  char *title = new char[len]; sprintf(title,"%s &0x3",histo->GetTitle());
  all0x3->SetLineColor(1);  // Overwrite color inherited from parent "p0x3"
  all0x3->SetName("&0x3"); 
  all0x3->SetTitle(title);
  
  // 0x6 OR 
  TH1D *all0x6 = new TH1D(*p0x6);
  all0x6->Add(p0x7); 
  all0x6->Add(p0xe); 
  all0x6->Add(p0xf);
  sprintf(title,"%s &0x6",histo->GetTitle());
  all0x6->SetLineColor(1);  // Overwrite color inherited from parent "p0x3"
  all0x6->SetName("&0x6"); 
  all0x6->SetTitle(title);
  
  //all0x3->SetEntries(histo->Integral(0,N1,7,7)+histo->Integral(0,N1,8,8)+histo->Integral(0,N1,15,15)+histo->Integral(0,N1,16,16));
  // 0xe OR 
  TH1D *all0xe = new TH1D(*p0xe);
  all0xe->Add(p0xf);
  sprintf(title,"%s &0xe",histo->GetTitle());
  all0xe->SetLineColor(1);  // Overwrite color inherited from parent "p0x3"
  all0xe->SetName("&0xe"); 
  all0xe->SetTitle(title);

  // TCanvas 
  TCanvas *C1; TPad *c1;
  C1 = new TCanvas("c1","TTrack",  900, 450);
  c1 = (TPad*)C1;
  c1->Clear(); c1->Divide(2,1);

  // WHICH OptStat? 
  int OptStat = gStyle->GetOptStat();
  printf("OptStat %d\n",OptStat);

  int csiz = 0;

  // 0x3 OR 
  c1->cd(1); 
  gStyle->SetOptStat(OptStat);
  all0x3->Draw();
  gPad->Update();
  {              // Edit title
    TPaveText *titre = (TPaveText*)gPad->GetPrimitive("title");
    Coord_t x1 = titre->GetX1NDC(), x2 = titre->GetX2NDC(), dx = x2-x1;
    x2 += .30*dx; titre->SetX2NDC(x2);
    Coord_t y1 = titre->GetY1NDC(), y2 = titre->GetY2NDC(), dy = y2-y1;
    y1 -= .8*dy; titre->SetY1NDC(y1);
    titre->SetTextSize(csiz);
    titre->Draw(); gPad->Update();
  }
  {
     TPaveStats *st = (TPaveStats*)gPad->GetPrimitive("stats");
    Coord_t x1 = st->GetX1NDC(), x2 = st->GetX2NDC(), dx = x2-x1;
    Coord_t X1 = st->GetX1(),    X2 = st->GetX2(),    dX = X2-X1;
    x1 -= .55*dx; x2 += (gPad->GetX2()-X2)*dx/dX;
    st->SetX1NDC(x1); st->SetX2NDC(x2);
    Coord_t y1 = st->GetY1NDC(), y2 = st->GetY2NDC(), dy = y2-y1;
    Coord_t Y1 = st->GetY1(),    Y2 = st->GetY2(),    dY = Y2-Y1;
    y1 -= -.20*dy;
    y2 += (gPad->GetY2()-Y2)*dy/dY;
    st->SetY1NDC(y1); st->SetY2NDC(y2);
    st->SetTextSize(csiz);
    st->Draw();
    TPaveStats *st2 = (TPaveStats*)st->Clone("sAll");
    dx = x2-x1; x1 -= dx; x2 -= dx; st2->SetX1NDC(x1); st2->SetX2NDC(x2);
    st2->Draw();
    gPad->Update();
  }
  p0x3->Draw("sames");
  {
    gPad->Update(); TPaveStats *st = (TPaveStats*)gPad->GetPrimitive("stats");
    TPaveStats *st2 = (TPaveStats*)st->Clone("s0x3");
    st2->SetTextColor(2); st2->Draw();
    gPad->Update();
  }
  p0x7->Draw("sames");
  {
    gPad->Update(); TPaveStats *st = (TPaveStats*)gPad->GetPrimitive("stats");
    TPaveStats *st2 = (TPaveStats*)st->Clone("s0x7");
    Coord_t y1 = st2->GetY1NDC(), y2 = st2->GetY2NDC(), dy = y2-y1;
    y1 -= dy; y2 -= dy; st2->SetY1NDC(y1); st2->SetY2NDC(y2);
    st2->SetTextColor(4); st2->Draw();
    gPad->Update();
  }
  p0xf->Draw("sames");
  {
    gPad->Update(); TPaveStats *st = (TPaveStats*)gPad->GetPrimitive("stats");
    Coord_t y1 = st->GetY1NDC(), y2 = st->GetY2NDC(), dy = y2-y1;
    y1 -= 2*dy; y2 -= 2*dy; st->SetY1NDC(y1); st->SetY2NDC(y2);
    st->SetTextColor(6); st->Draw();
    gPad->Update();
  }

  // 0x6 OR 

  c1->cd(2);
  all0x6->Draw();
  gPad->Update();
  {              // Edit title
    TPaveText *titre = (TPaveText*)gPad->GetPrimitive("title");
    Coord_t x1 = titre->GetX1NDC(), x2 = titre->GetX2NDC(), dx = x2-x1;
    x2 += .30*dx; titre->SetX2NDC(x2);
    Coord_t y1 = titre->GetY1NDC(), y2 = titre->GetY2NDC(), dy = y2-y1;
    y1 -= .8*dy; titre->SetY1NDC(y1);
    titre->SetTextSize(csiz);
    titre->Draw(); gPad->Update();
  }
  {
     TPaveStats *st = (TPaveStats*)gPad->GetPrimitive("stats");
    Coord_t x1 = st->GetX1NDC(), x2 = st->GetX2NDC(), dx = x2-x1;
    Coord_t X1 = st->GetX1(),    X2 = st->GetX2(),    dX = X2-X1;
    x1 -= .55*dx; x2 += (gPad->GetX2()-X2)*dx/dX;
    st->SetX1NDC(x1); st->SetX2NDC(x2);
    Coord_t y1 = st->GetY1NDC(), y2 = st->GetY2NDC(), dy = y2-y1;
    Coord_t Y1 = st->GetY1(),    Y2 = st->GetY2(),    dY = Y2-Y1;
    y1 -= -.20*dy;
    y2 += (gPad->GetY2()-Y2)*dy/dY;
    st->SetY1NDC(y1); st->SetY2NDC(y2);
    st->SetTextSize(csiz);
    st->Draw();
    TPaveStats *st2 = (TPaveStats*)st->Clone("sAll");
    dx = x2-x1; x1 -= dx; x2 -= dx; st2->SetX1NDC(x1); st2->SetX2NDC(x2);
    st2->Draw();
    gPad->Update();
  }
  p0x6->Draw("sames");
  {
    gPad->Update(); TPaveStats *st = (TPaveStats*)gPad->GetPrimitive("stats");
    TPaveStats *st2 = (TPaveStats*)st->Clone("s0x6");
    st2->SetTextColor(2); st2->Draw();
    gPad->Update();
  }
  p0x7->Draw("sames");
  {
    gPad->Update(); TPaveStats *st = (TPaveStats*)gPad->GetPrimitive("stats");
    TPaveStats *st2 = (TPaveStats*)st->Clone("s0x7");
    Coord_t y1 = st2->GetY1NDC(), y2 = st2->GetY2NDC(), dy = y2-y1;
    dy = y2-y1; y1 -= dy; y2 -= dy; st2->SetY1NDC(y1); st2->SetY2NDC(y2);
    st2->SetTextColor(4); st2->Draw();
    gPad->Update();
  }
  p0xe->Draw("sames");
  {
    gPad->Update(); TPaveStats *st = (TPaveStats*)gPad->GetPrimitive("stats");
    TPaveStats *st2 = (TPaveStats*)st->Clone("s0xe");
    Coord_t y1 = st2->GetY1NDC(), y2 = st2->GetY2NDC(), dy = y2-y1;
    dy = y2-y1; y1 -= 2*dy; y2 -= 2*dy; st2->SetY1NDC(y1); st2->SetY2NDC(y2);
    st2->SetTextColor(3); st2->Draw();
    gPad->Update();
  }
  p0xf->Draw("sames");
  {
    gPad->Update(); TPaveStats *st = (TPaveStats*)gPad->GetPrimitive("stats");
    Coord_t y1 = st->GetY1NDC(), y2 = st->GetY2NDC(), dy = y2-y1;
    dy = y2-y1; y1 -= 3*dy; y2 -= 3*dy; st->SetY1NDC(y1); st->SetY2NDC(y2);
    st->SetTextColor(6); st->Draw();
    gPad->Update();
  }
  
  return histo;

}
#endif
   
#ifndef __CINT__    
//________________________________________
TCanvas* Macro::BitPattern( vector<string> files, TCanvas *cv , int zone_mask)
{
  vector< TH1* > exp;
  vector< TH1* > fnd;
  vector< TH1* > dets;
  
  // first value for dets histograms, to have proper scale
  double firstDet[] = {0.1,0.4,0.6,0.8};
  
  // expected bitmaps
  if( zone_mask & (1<<0) ) exp.push_back( Utils::Get(files, "Traffic/TTrackRefine/Expected0" ) );
  if( zone_mask & (1<<1) ) exp.push_back( Utils::Get(files, "Traffic/TTrackRefine/Expected1" ) );
  if( zone_mask & (1<<2) ) exp.push_back( Utils::Get(files, "Traffic/TTrackRefine/Expected2" ) );
  if( zone_mask & (1<<3) ) exp.push_back( Utils::Get(files, "Traffic/TTrackRefine/Expected3" ) );
  
  // expected bitmaps
  if( zone_mask & (1<<0) ) fnd.push_back( Utils::Get(files, "Traffic/TTrackRefine/Found0" ) );
  if( zone_mask & (1<<1) ) fnd.push_back( Utils::Get(files, "Traffic/TTrackRefine/Found1" ) );
  if( zone_mask & (1<<2) ) fnd.push_back( Utils::Get(files, "Traffic/TTrackRefine/Found2" ) );
  if( zone_mask & (1<<3) ) fnd.push_back( Utils::Get(files, "Traffic/TTrackRefine/Found3" ) );

  // detector type
  if( zone_mask & (1<<0) ) dets.push_back( Utils::Get(files, "Traffic/TTrackRefine/DetTyp0" ) );
  if( zone_mask & (1<<1) ) dets.push_back( Utils::Get(files, "Traffic/TTrackRefine/DetTyp1" ) );
  if( zone_mask & (1<<2) ) dets.push_back( Utils::Get(files, "Traffic/TTrackRefine/DetTyp2" ) );
  if( zone_mask & (1<<3) ) dets.push_back( Utils::Get(files, "Traffic/TTrackRefine/DetTyp3" ) );
    
  TCanvas *c1 = (cv) ? cv:new TCanvas("c1","TTrack",  700, 700);
  if( !cv ) {
    gROOT->SetStyle("Plain");
    c1->Clear(); 
  }
  
  Utils::DivideCanvas( c1, exp.size() );
  
  for( unsigned int i=0; i < exp.size(); i++ ) 
  if( !( exp[i] && fnd[i] ) ) { 
    cout << "Macro::BitPattern - nothing found for zone " << i << ".\n";
    continue;
  } else {
    c1->cd(i+1);
    TH1* eff = (TH1*) exp[i]->Clone();
    Utils::HDiv( fnd[i], exp[i], eff );
    eff->Draw("H");
    if( dets[i] ) {
      dets[i]->SetLineColor(2);
      
      // scale histogram so that first bin matches firstDet.
      dets[i]->Scale( firstDet[i]/(double) dets[i]->GetBinContent( 1 ) );       
      dets[i]->Draw("same");
    }
  }
  
  return c1;
}
#endif

//________________________________________________________  
void Macro::Eval( const char* fileselection, int trigger_mask )
{

  vector<string> files = Utils::GetFiles( fileselection );
  TH2* tr_11 = Utils::Get2D( fileselection, "Traffic/RDmonitor/tr_11");
  TH2* tr_13 = Utils::Get2D( fileselection, "Traffic/RDmonitor/tr_13");
  TH2* tr_15 = Utils::Get2D( fileselection, "Traffic/RDmonitor/tr_15");
  TH1* hZTmm = Utils::Get( fileselection, "CsAverPattern/Distribution/hZTmm");
  TH1* hZ    = Utils::Get( fileselection, "CsAverPattern/Secondaries/hZ");
  TH1* hNTrk = Utils::Get( fileselection, "CsKalmanFitting/hNTrk");
  
  int    N1 = (tr_11) ? tr_11->GetNbinsX()+1:0; 
  int    sum1 = 0, sum2 = 0, sum3 = 0, sum_t = 0;
  double chi1 = 0, chi2 = 0, chi3 = 0, chi_t = 0;
  
  if( tr_11 ) {
    // 1st spectro: sum slices 0x3, 0x7, 0xf
    sum1 += int(tr_11->Integral(0,N1,4,4));
    sum1 += int(tr_11->Integral(0,N1,8,8));
    sum1 += int(tr_11->Integral(0,N1,16,16));
    
    chi1 += tr_11->ProjectionX("_0x3",4,4)->GetMean()  *int(tr_11->Integral(0,N1,4,4));
    chi1 += tr_11->ProjectionX("_0x3",8,8)->GetMean()  *int(tr_11->Integral(0,N1,8,8));
    chi1 += tr_11->ProjectionX("_0x3",16,16)->GetMean()*int(tr_11->Integral(0,N1,16,16));
    chi1 /= double( sum1 );
    
    // 2nd spectro: sum slices 0x6, 0x7, 0xe, 0xf
    sum2 += int(tr_11->Integral(0,N1,7,7));
    sum2 += int(tr_11->Integral(0,N1,8,8));
    sum2 += int(tr_11->Integral(0,N1,15,15));
    sum2 += int(tr_11->Integral(0,N1,16,16));
    
    chi2 += tr_11->ProjectionX("_px",7,7)->GetMean()  *int(tr_11->Integral(0,N1,7,7));
    chi2 += tr_11->ProjectionX("_px",8,8)->GetMean()  *int(tr_11->Integral(0,N1,8,8));
    chi2 += tr_11->ProjectionX("_px",15,15)->GetMean()*int(tr_11->Integral(0,N1,15,15));
    chi2 += tr_11->ProjectionX("_px",16,16)->GetMean()*int(tr_11->Integral(0,N1,16,16));
    chi2 /= double( sum2 );
    
    // muons: sum slices 0xe, 0xf
    sum3 += int(tr_11->Integral(0,N1,15,15));
    sum3 += int(tr_11->Integral(0,N1,16,16));
    chi3 += tr_11->ProjectionX("_px",15,15)->GetMean()*int(tr_11->Integral(0,N1,15,15));
    chi3 += tr_11->ProjectionX("_px",16,16)->GetMean()*int(tr_11->Integral(0,N1,16,16));
    chi3 /= double( sum3 );
    
    // Everything
    sum_t+= sum1;
    sum_t+= int(tr_11->Integral(0,N1,7,7));
    sum_t+= int(tr_11->Integral(0,N1,15,15));
    
    chi_t+= chi1*sum1;
    chi_t+= tr_11->ProjectionX("_px",7,7)->GetMean()  *int(tr_11->Integral(0,N1,7,7));
    chi_t+= tr_11->ProjectionX("_px",15,15)->GetMean()*int(tr_11->Integral(0,N1,15,15));
    chi_t/= double( sum_t );
  }
  
  // dump number of tracks
  printf("TTtracks = %d %d %d %d\n",
    sum_t,
    sum1,
    sum2,
    sum3 );
  
  // get number of triggers  
  int n_trigs = 0;
  if( tr_13 ) {
    if( trigger_mask )
    for( int itrig = 0; itrig<32; itrig++ ) 
    if( itrig&trigger_mask ) n_trigs += int(tr_13->Integral(itrig+1,itrig+1));
    
    if( n_trigs==0 ) {
      n_trigs = int(tr_13->GetEntries());
      trigger_mask = 0xffffffff;
    }
  }
  
  // dump number of tracks/trigger
  if( n_trigs ) 
  printf("TTtracks = %.2f %.2f %.2f %.3f\n",
	  (double)sum_t/n_trigs,
    (double)sum1/n_trigs,
    (double)sum2/n_trigs,
	  (double)sum3/n_trigs );
  
  
  // dump average chi2
  printf("chi2     = %.2f %.2f %.2f %.3f\n",
	  chi_t,
    chi1,
    chi2,
	  chi3 );
     
  // vertices
  int nmumup =   (hZTmm) ? int(hZTmm->Integral())  :0;    // number of vertices
  int nmumup_t = (hZTmm) ? int(hZTmm->GetEntries()):0;    // Number of mu/mu'
  int nvsec =    (hZ)    ? int(hZ->GetEntries())   :0;    // Number of Secondaries
  int nvprim =   (hNTrk) ? int(hNTrk->Integral())  :0;    // Number of primary vertices
  float ntracks = (hNTrk) ? hNTrk->GetMean():0;           // Number of tracks in vertex
  printf("Vertices (prim(tracks),mu/mu',sec) = %d(%.3f) %d(%d) %d\n", 
    nvprim,
    ntracks,
    nmumup,
    nmumup_t,
    nvsec );
  if( n_trigs )  
  printf("Vertices (prim(tracks),mu/mu',sec) = %.1f(%.3f) %.1f(%.1f) %.3f\n",
	  100.*nvprim/n_trigs,
    ntracks,
    100.*nmumup/n_trigs,
	  100.*nmumup_t/n_trigs,
	  100.*nvsec/n_trigs  );
  
  // beam
  int n_beam0 = (tr_15) ? int(tr_15->Integral(1,1)):0;
  int n_beam1 = (tr_15) ? int(tr_15->Integral(2,2)):0;
  int n_beams = (tr_15) ? int(tr_15->GetEntries()-n_beam0-n_beam1):0;
  
  // Correct for non processed triggers
  if( tr_13 ) n_beam0 -= int(tr_13->GetEntries()-n_trigs);
  
  // dump
  if( tr_15 ) printf("Beam (0,1,>1) =  %5d %5d %5d\n",n_beam0,n_beam1,n_beams);
  if( tr_15 && n_trigs ) printf("Beam (0,1,>1) =   %.2f  %.2f  %.2f\n",
	   (double)n_beam0/n_trigs,
	   (double)n_beam1/n_trigs,
     (double)n_beams/n_trigs);
}
