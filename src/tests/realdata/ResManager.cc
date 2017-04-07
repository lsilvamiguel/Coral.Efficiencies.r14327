// $Id: ResManager.cc,v 1.6 2009/04/16 09:34:55 aaust Exp $
 
/*!
   \file    ResManager.cc
   \brief   Resolution Managment Interface Class.
   \author  Hugo Pereira
   \version $Revision: 1.6 $
   \date    $Date: 2009/04/16 09:34:55 $
*/

#include "ResManager.h"
#include "DetectorInfo.h"
#include "DetFileManager.h"
#include "Utils.h"
#include "Fit.h"
#include "TH2Fit.h"

#include <TChain.h>
#include <TCut.h>
#include <TH1.h>
#include <TH2.h>
#include <TF1.h>
#include <TCanvas.h>
#include <TMinuit.h>
#include <TText.h>
#include <TVirtualPad.h>

#include <iostream>
#include <fstream>

using namespace std;

ClassImp( ResManager )

//__________________________________________________________________
#define resout "res.out"
void ResManager::FitDU( const char* selection, const char* cut, int nEnt )
{
  if( !_Select( selection ) ) return;
  
  // handle canvas
  cv_->Clear();
  Utils::DivideCanvas( cv_, dets_.size() );
  
  
  for( unsigned int iInf = 0; iInf < dets_.size(); iInf++ ) 
  {
    double resid = 0;
    double du = 0;
    
    DetectorInfo* d = dets_[iInf];
    cv_->cd( iInf+1 );
    TH1 *h = _Get("T_duVect", dets_[iInf], cut, nEnt );

    //TH1 *h2 = _Get("T_dvVect", dets_[iInf], cut, nEnt );
    //    h2->SetLineColor(2);
        
    cout << "ResManager::FitDU - det \"" 
      << d->TBName_.c_str() 
      << "\": " 
      << int((h)?h->GetEntries():-1)
      << " entries - ";
    if( !( h&&h->GetEntries() ) ) {
      cout << endl;
      continue;
    }
    
    h->Draw();

    //    if( d->TBName_.substr(0,2) == "GP" && d->TBName_.substr(4,1) == "P") h2->Draw("same");
    
    if( OneGaus_ ) {

      h->Fit( "gaus", "0Q" );     
      TF1* f = h->GetFunction("gaus");
      if( f ) {
        f->SetLineColor(4);
        f->SetLineWidth(1);
        f->Draw("same");
        resid = fabs(f->GetParameter(2));
        du    = f->GetParameter(1);
      }
    
    } else {
      
      // define function used for the fit and limits 
      const unsigned int nP = 6;
      double xMin = h->GetXaxis()->GetBinCenter( 1 );
      double xMax = h->GetXaxis()->GetBinCenter( h->GetNbinsX() );
      char name[50];
      snprintf(name, 50, "du1DFit%d",iInf);
      TF1* f = new TF1(name, Fit::GG , xMin, xMax, nP);
      f->SetLineColor(4);
      f->SetLineWidth(1);
      
      // define starting parameters  
      // first gaussian
      f->SetParameter( 0, h->GetMaximum() );
      f->SetParameter( 1, h->GetBinCenter( h->GetMaximumBin() ) );
      f->SetParameter( 2, 0.5*h->GetRMS( ) );
    
      // second gaussian
      f->SetParameter( 3, 0.2*h->GetMaximum() );
      f->SetParameter( 4, h->GetBinCenter( h->GetMaximumBin() ) );
      f->SetParameter( 5, 2*h->GetRMS( ) );
      
      // do the fit
      h->Fit( f, "0" ); 
      f->Draw("same");
    
      // Draw The largest gaussian separately
      TF1* fg = new TF1("gaus", (const char*) "gaus(0)", xMin, xMax );
      fg->SetLineColor(2);
      fg->SetLineWidth(1);
      if( f->GetParameter( 0 ) < f->GetParameter( 3 ) ) {
        fg->SetParameter( 0, f->GetParameter( 0 ) );
        fg->SetParameter( 1, f->GetParameter( 1 ) );
        fg->SetParameter( 2, f->GetParameter( 2 ) );
      } else {
        fg->SetParameter( 0, f->GetParameter( 3 ) );
        fg->SetParameter( 1, f->GetParameter( 4 ) );
        fg->SetParameter( 2, f->GetParameter( 5 ) );
      }
      fg->Draw("same");
  
      // store result for alignment
      resid = fabs( ( f->GetParameter(0) > f->GetParameter(3) )? f->GetParameter(2):f->GetParameter(5) );
      du    = ( f->GetParameter(0) > f->GetParameter(3) )? f->GetParameter(1):f->GetParameter(4);
    }
  
    char* text = new char[200];
    sprintf(text,"%6.1fum",1000.0*resid);
    TText* tt = new TText();
    tt->SetTextColor( 4 ); 
    tt->SetTextSize(0.06);
    tt->DrawTextNDC( 0.15, 0.15, text );
    cv_->Update();
    
    cout << " res=" << resid << "mm. ave=" << du << "mm.\n";
    d->resid_ = resid;
    
    // dump to file
    ofstream out( resout, ios::app );
    snprintf(text, 200, "ResManager::FitDU - \"%s\" [%10.3f] [ %10.4f %10.4f ] -> [ %10.4f %10.4f ].\n",
      d->TBName_.c_str(),
      du,
      d->xcm_/10,
      d->ycm_/10,
      (d->xcm_-d->rotM_(0,0)*du)/10,
      (d->ycm_-d->rotM_(1,0)*du)/10 );
    out<<text;
    out.close();
    SafeDelete( text );
  }
  
  return;
}

//_________________P_I_X_E_L_____________________________________________
void ResManager::FitDV( const char* selection, const char* cut, int nEnt )
{
  if( !_Select( selection ) ) return;
  
  // handle canvas
  cv_->Clear();
  Utils::DivideCanvas( cv_, dets_.size() );
  
  
  for( unsigned int iInf = 0; iInf < dets_.size(); iInf++ ) 
  {
    double resid = 0;
    double du = 0;
    
    DetectorInfo* d = dets_[iInf];
    cv_->cd( iInf+1 );
    TH1 *h = _Get("T_dvVect", dets_[iInf], cut, nEnt );

    cout << "ResManager::FitDV - det \"" 
      << d->TBName_.c_str() 
      << "\": " 
      << int((h)?h->GetEntries():-1)
      << " entries - ";
    if( !( h&&h->GetEntries() ) ) {
      cout << endl;
      continue;
    }
    
    h->Draw();

    //    if( d->TBName_.substr(0,2) == "GP" && d->TBName_.substr(4,1) == "P") h2->Draw("same");
    
    if( OneGaus_ ) {

      h->Fit( "gaus", "0Q" );     
      TF1* f = h->GetFunction("gaus");
      if( f ) {
        f->SetLineColor(4);
        f->SetLineWidth(1);
        f->Draw("same");
        resid = fabs(f->GetParameter(2));
        du    = f->GetParameter(1);
      }
    
    } else {
      
      // define function used for the fit and limits 
      const unsigned int nP = 6;
      double xMin = h->GetXaxis()->GetBinCenter( 1 );
      double xMax = h->GetXaxis()->GetBinCenter( h->GetNbinsX() );
      char name[50];
      snprintf(name, 50, "dv1DFit%d",iInf);
      TF1* f = new TF1(name, Fit::GG , xMin, xMax, nP);
      f->SetLineColor(4);
      f->SetLineWidth(1);
      
      // define starting parameters  
      // first gaussian
      f->SetParameter( 0, h->GetMaximum() );
      f->SetParameter( 1, h->GetBinCenter( h->GetMaximumBin() ) );
      f->SetParameter( 2, 0.5*h->GetRMS( ) );
    
      // second gaussian
      f->SetParameter( 3, 0.2*h->GetMaximum() );
      f->SetParameter( 4, h->GetBinCenter( h->GetMaximumBin() ) );
      f->SetParameter( 5, 2*h->GetRMS( ) );
      
      // do the fit
      h->Fit( f, "0" ); 
      f->Draw("same");
    
      // Draw The largest gaussian separately
      TF1* fg = new TF1("gaus", (const char*) "gaus(0)", xMin, xMax );
      fg->SetLineColor(2);
      fg->SetLineWidth(1);
      if( f->GetParameter( 0 ) < f->GetParameter( 3 ) ) {
        fg->SetParameter( 0, f->GetParameter( 0 ) );
        fg->SetParameter( 1, f->GetParameter( 1 ) );
        fg->SetParameter( 2, f->GetParameter( 2 ) );
      } else {
        fg->SetParameter( 0, f->GetParameter( 3 ) );
        fg->SetParameter( 1, f->GetParameter( 4 ) );
        fg->SetParameter( 2, f->GetParameter( 5 ) );
      }
      fg->Draw("same");
  
      // store result for alignment
      resid = fabs( ( f->GetParameter(0) > f->GetParameter(3) )? f->GetParameter(2):f->GetParameter(5) );
      du    = ( f->GetParameter(0) > f->GetParameter(3) )? f->GetParameter(1):f->GetParameter(4);
    }
  
    char* text = new char[200];
    sprintf(text,"%6.1fum",1000.0*resid);
    TText* tt = new TText();
    tt->SetTextColor( 4 ); 
    tt->SetTextSize(0.06);
    tt->DrawTextNDC( 0.15, 0.15, text );
    cv_->Update();
    
    cout << " res=" << resid << "mm. ave=" << du << "mm.\n";
    d->resid_ = resid;
    
    // dump to file
    ofstream out( resout, ios::app );
    snprintf(text, 200, "ResManager::FitDV - \"%s\" [%10.3f] [ %10.4f %10.4f ] -> [ %10.4f %10.4f ].\n",
      d->TBName_.c_str(),
      du,
      d->xcm_/10,
      d->ycm_/10,
      (d->xcm_+d->rotM_(1,0)*du)/10,
      (d->ycm_-d->rotM_(0,0)*du)/10 );
    out<<text;
    out.close();
    SafeDelete( text );
  }
  
  return;
}

//__________________________________________________________________
void ResManager::FitDU_MWPC( const char* selection, const char* cut, int nEnt )
{
  if( !_Select( selection ) ) return;
  
  // handle canvas
  cv_->Clear();
  Utils::DivideCanvas( cv_, dets_.size() );
  
  
  for( unsigned int iInf = 0; iInf < dets_.size(); iInf++ ) 
  {
     
    DetectorInfo* d = dets_[iInf];
    cv_->cd( iInf+1 );
    TH1 *h = _Get("T_duVect", dets_[iInf], cut, nEnt );
        
    cout << "ResManager::FitDU_MWPC - det \"" 
      << d->TBName_.c_str() << "\": " 
      << int((h)?h->GetEntries():-1)
      << " entries - ";
    if( !( h&&h->GetEntries() ) ) {
      cout << endl;
      continue;
    }
    
    h->Draw();
    
    double resid, du;
    
    if( OneGaus_ ) {
      // define function used for the fit and limits 
      const unsigned int nP = 4;
      double xMin = h->GetXaxis()->GetBinCenter( 1 );
      double xMax = h->GetXaxis()->GetBinCenter( h->GetNbinsX() );
      char name[50];
      snprintf(name, 50, "du1DFit_MWPC%d",iInf);
      TF1* f = new TF1(name, Fit::EE , xMin, xMax, nP); // two Erf
      f->SetLineColor(4);
      f->SetLineWidth(1);
      
      // define starting parameters  
      f->SetParameter( 0, h->GetMaximum() );
      f->SetParameter( 1, h->GetBinCenter(h->GetMaximumBin() ) );
      f->SetParameter( 2, 0.5*h->GetRMS( ) );
      f->SetParameter( 3, 1.0 );
     
      h->Fit( f, "0" ); 
      f->Draw("same");
    
      // store result for alignment
      resid = f->GetParameter(2)*2.0/sqrt(12.0);
      du =    f->GetParameter(1);
    } else {
      // define function used for the fit and limits 
      const unsigned int nP = 7;
      double xMin = h->GetXaxis()->GetBinCenter( 1 );
      double xMax = h->GetXaxis()->GetBinCenter( h->GetNbinsX() );
      char name[50];
      snprintf(name, 50, "du1DFit_MWPC%d",iInf);
      TF1* f = new TF1(name, Fit::EEG , xMin, xMax, nP);  // two Erf and a gaussian
      f->SetLineColor(4);
      f->SetLineWidth(1);
      
      // define starting parameters  
      f->SetParameter( 0, h->GetMaximum() );
      f->SetParameter( 1, h->GetBinCenter(h->GetMaximumBin() ) );
      f->SetParameter( 2, 0.5*h->GetRMS( ) );
      f->SetParameter( 3, 1.0 );
    
      f->SetParameter( 4, 0.2*h->GetMaximum() );
      f->SetParameter( 5, h->GetBinCenter(h->GetMaximumBin() ) );
      f->SetParameter( 6, 2*h->GetRMS() );
    
      h->Fit( f, "0" ); 
      f->Draw("same");
    
      // store result for alignment
      resid = f->GetParameter(2)*2.0/sqrt(12.0);
      du =    f->GetParameter(1);
    }
    
    char* text = new char[200];
    sprintf(text,"%6.1fum",1000.0*resid);
    TText* tt = new TText();
    tt->SetTextColor( 4 ); 
    tt->SetTextSize(0.06);
    tt->DrawTextNDC( 0.15, 0.15, text );
    cv_->Update();
  
    cout << " res=" << resid << "mm. ave=" << du << "mm.\n";
    d->resid_ = resid;
    
    // dump to file
    ofstream out( resout, ios::app );
    snprintf(text, 200, "ResManager::FitDU_MWPC - \"%s\" [%10.3f] [ %10.4f %10.4f ] -> [ %10.4f %10.4f ].\n",
      d->TBName_.c_str(),
      du,
      d->xcm_/10,
      d->ycm_/10,
      (d->xcm_+d->rotM_(0,0)*du)/10,
      (d->ycm_+d->rotM_(1,0)*du)/10 );
    out << text;
    out.close();
    SafeDelete( text );
  }
     
  return;
}

//__________________________________________________________________
void ResManager::FitDUvsU( const char* selection, const char* cut, int nEnt )
{
  if( !df_ )    { cout << "ResManager::FitDUvsU - no detectorfile loaded.\n"; return; }
  if( !chain_ ) { cout << "ResManager::FitDUvsU - no tracks loaded.\n"; return; }
  
  // Get Selection
  cout << "ResManager::FitDUvsU - " 
    << _Select( selection, cut, nEnt ) 
    << " detector(s) selected.\n";
  if( !dets_.size() ) return;
  
  // handle canvas
  cv_->Clear();
  Utils::DivideCanvas( cv_, dets_.size() );
  
  
  for( unsigned int iInf = 0; iInf < dets_.size(); iInf++ ) 
  {
     
    DetectorInfo* d = dets_[iInf];
    cv_->cd( iInf+1 );
    TH1 *h = _Get("T_duVect:T_uLx", dets_[iInf], cut, nEnt );
        
    cout << "ResManager::FitDUvsU - det \"" 
      << d->TBName_.c_str() << "\": " 
      << int((h)?h->GetEntries():-1)
      << " entries - ";
    if( !( h&&h->GetEntries() ) ) {
      cout << endl;
      continue;
    }
    
    h->Draw("box");
    
    // define fit function and histo 2D fitter
    double uMin =  h->GetXaxis()->GetBinCenter( 1 );
    double uMax =  h->GetXaxis()->GetBinCenter( h->GetNbinsX() );
    TF1* f = new TF1("P1", Fit::P1, uMin, uMax, 2);
    f->SetParameter( 0, 0 );
    f->SetParameter( 1, 0 );
    f->SetLineColor( 2 );
    f->SetLineWidth( 2 );
    
    TH2Fit fit( f, 2 );
    fit.Fit( (TH2*) h, uMin, uMax );
    
    f->Draw("same");
  
    char* text = new char[100];
    sprintf(text,"dU=%.3g+%.3g*U", f->GetParameter(0), f->GetParameter(1) );
    TText* tt = new TText();
    tt->SetTextColor( 4 ); 
    tt->SetTextSize(0.06);
    tt->DrawTextNDC( 0.15, 0.15, text );
    SafeDelete( text );

    gPad->SetLogz();
    cv_->Update();

    cout << " ave=" << f->GetParameter( 0 ) << "mm - slope=" << f->GetParameter(1) << ".\n";
  }
}

//__________________________________________________________________
void ResManager::FitDUvsV( const char* selection, const char* cut, int nEnt  )
{
  if( !_Select( selection ) ) return;
 
  // handle canvas
  cv_->Clear();
  Utils::DivideCanvas( cv_, dets_.size() );
  
  for( unsigned int iInf = 0; iInf < dets_.size(); iInf++ ) 
  {
     
    DetectorInfo* d = dets_[iInf];
    cv_->cd( iInf+1 );
    TH1 *h = _Get("T_duVect:T_vLx", dets_[iInf], cut, nEnt );
        
    cout << "ResManager::FitDUvsV - det \"" 
      << d->TBName_.c_str() << "\": " 
      << int((h)?h->GetEntries():-1)
      << " entries - ";
    if( !( h&&h->GetEntries() ) ) {
      cout << endl;
      continue;
    }
    
    h->Draw("box");
    
    // define fit function and histo 2D fitter
    double uMin =  h->GetXaxis()->GetBinCenter( 1 );
    double uMax =  h->GetXaxis()->GetBinCenter( h->GetNbinsX() );
    TF1* f = new TF1("P1", Fit::P1, uMin, uMax, 2);
    f->SetParameter( 0, 0 );
    f->SetParameter( 1, 0 );
    f->SetLineColor( 2 );
    f->SetLineWidth( 2 );
    
    TH2Fit fit( f, 2 );
    fit.Fit( (TH2*) h, uMin, uMax );
    
    f->Draw("same");
  
    char* text = new char[100];
    sprintf(text,"dV=%.3g+%.3g*U", f->GetParameter(0), f->GetParameter(1) );
    TText* tt = new TText();
    tt->SetTextColor( 4 ); 
    tt->SetTextSize(0.06);
    tt->DrawTextNDC( 0.15, 0.15, text );
    SafeDelete( text );
    
    gPad->SetLogz();
    cv_->Update();
    
    cout << " ave=" << f->GetParameter( 0 ) << "mm - slope=" << f->GetParameter(1) << ".\n";
    
  }
    
}

//__________________________________________________________________
void ResManager::DrawDU_LR( const char* selection, const char* cut,
    double OffsetN,
    double OffsetP,
    int nEnt )
{
  if( !_Select( selection ) ) return;
  
  // handle canvas
  cv_->Clear();
  Utils::DivideCanvas( cv_, dets_.size() );
  
  for( unsigned int iInf = 0; iInf < dets_.size(); iInf++ ) 
  {
     
    DetectorInfo* d = dets_[iInf];
    cv_->cd( iInf+1 );
    
    string nameN = d->TBName_ + "_N";
    string nameP = d->TBName_ + "_P";
    
    TH1 *h;
    if( (h = (TH1*) gDirectory->Get(nameN.c_str())) ) SafeDelete(h);
    if( (h = (TH1*) gDirectory->Get(nameP.c_str())) ) SafeDelete(h);
   
    // Get r<0 histograms
    char *varN = new char[100];
    if( OffsetN ) sprintf( varN, "T_duVect+%f", OffsetN );
    else sprintf( varN, "T_duVect" );
    
    TCut cutN = TCut( cut ) && TCut("T_rVect<0") ;
    TH1 *hN = _Get( varN, dets_[iInf], (const char*) cutN, nEnt );
    SafeDelete( varN );
    
    hN->SetName(  nameN.c_str() );
    hN->SetTitle( nameN.c_str() );
    
    // Get r>0 histograms
    char *varP = new char[100];
    if( OffsetN ) sprintf( varP, "T_duVect+%f", OffsetP );
    else sprintf( varP, "T_duVect" );
    
    TCut cutP = TCut( cut ) && TCut("T_rVect>0") ;
    TH1 *hP = _Get( varP, dets_[iInf], (const char*) cutP, nEnt );
    SafeDelete( varP );
   
    hP->SetName(  nameP.c_str() );
    hP->SetTitle( nameP.c_str() );

    // check histograms
    if( (!(hN && hN->GetEntries())) && (!(hP && hP->GetEntries())) ) continue;
    
    // Change histogram colors
    if( hN ) hN->SetLineColor( 4 );
    if( hP ) hP->SetLineColor( 2 );
    
    // Get Histo with max value
    TH1* h1, *h2;
    if( hN && ((!hP) || hP->GetMaximum()<hN->GetMaximum() ) ) { h1 = hN; h2 = hP; }
    else{ h1 = hP; h2 = hN; }
    
    //! Dump number of entries
    cout << "ResManager::DrawDU_LR - det \"" 
      << d->TBName_.c_str() 
      << "\": " 
      << int((hN)?hN->GetEntries():-1)
      << "/" 
      << int((hP)?hP->GetEntries():-1)
      << " entries.\n";
    
    //!Draw histograms
    if( h1 ) {
      h1->Draw();
      if( h2 ) h2->Draw("same");
    } else if( h2 ) h2->Draw();
    cv_->Update();
  }

}

//_______________________________________________________________________________
bool ResManager::FitResolutions( const char* selection, const char* cut, int nEnt )
{
  if( !_Select( selection ) ) return false;
 
  // set residuals according to selection using cut
  cv_->Clear();
  FitDU( selection, cut, nEnt );
  
  // dump resolutions to screen
  for( unsigned int i = 0; i < dets_.size(); i++ )
  if( dets_[i]->resid_ < 0 ) { 
    cout << "ResManager::FitResolutions - ERROR: troubles for \"" << dets_[i]->TBName_ << "\" residual.\n";
    return false;
  } else cout << "ResManager::FitResolution - INFO: \"" << dets_[i]->TBName_ << "\" -> " << dets_[i]->resid_ << "mm.\n";
  
  // initialise minuit
  TMinuit *gMinuit_ = new TMinuit((int)dets_.size()); 
  gMinuit_->mninit(5,6,7);     // Logical units
  gMinuit_->SetFCN(Fit::resFcn);

  cout << "ResManager::FitResolution - INFO: TMinuit initialization performed." << endl;
  
  double dSig = .00001;
  int error;
    
  for( unsigned int i=0; i< dets_.size(); i++ ) {
    char pName[16];
    sprintf(pName,"res_%s",dets_[i]->TBName_.c_str() );
    gMinuit_->mnparm(i,pName,dets_[i]->res_,dSig, 0,0,error);
    if (error) {
      cout << "ResManager::FitResolution - INFO: Troubles defining parameter " << i << ".\n";
      return false;
    }
  }
  
  double flag0 = 0, flag1 = 1, flag3 = 3;
  gMinuit_->mnexcm("CALL FCN",  &flag1 ,1,error);
  gMinuit_->mnexcm("SET PRINT", &flag0 ,1,error);
  gMinuit_->mnexcm("SET NOGradient", &flag0 ,1,error);
  gMinuit_->mnexcm("MIGRAD",   &flag0 ,0,error);
  gMinuit_->mnexcm("MINOS",    &flag0 ,0,error);
  gMinuit_->mnexcm("CALL FCN", &flag3 ,1,error);
  
  cout << endl << "--------------------------------------------" << endl;
  cout << "ResManager::FitResolution - INFO: TMinuit minimisation performed." << endl;
  cout << "ResManager::FitResolution - INFO: Asuming Fast_Minimisation." << endl;
  
  // store Minuit results into vectors
  vector< double > resV;    resV.clear();
  vector< double > resErrV; resErrV.clear();
  for( unsigned int i=0; i<dets_.size(); i++ ) {
    double res, resErr;
    gMinuit_->GetParameter(i, res, resErr);
    resV.push_back( res );
    resErrV.push_back( resErr );
  }

  // dump to screen
  ofstream out( resout, ios::app );
  out << "===";
  for( unsigned int i = 0; i < fileSelection_.size(); i++ ) out << " " << fileSelection_[i];
  out << endl;
  out << "===" << string( Utils::GetTimeStamp( ) ) << endl;

  char* text = new char[200];
  
  for( unsigned int i=0; i<dets_.size(); i++ ) {
    snprintf(text, 200, "%8s %10.5f %10.5f %10.5f %10.5f %10.5f\n",
      dets_[i]->TBName_.c_str(), 
      resV[i],  resErrV[i], 
      _ResCalc( i, resV ),
      sqrt( pow( _ResCalc( i, resV ), 2 ) + pow( resV[i], 2 ) ),
      dets_[i]->resid_  
    );
    cout<<text<<flush;
    out<<text;
  }
  
  snprintf(text, 200, "%8s %10s %10s %10s %10s %10s\n\n", "TBName", "res", "error", "tr_res", "residF" ,"residI");
  cout<<text<<flush;
  out<<text;
  out.close();
  SafeDelete( text );
  SafeDelete(gMinuit_);
  return true;
}

//_____________________________________________________________________________
double ResManager::_ResCalc( unsigned int it, vector< double > res )
{
  // check vector sizes
  if( res.size() != dets_.size() ) { 
    cout << "AlManager::_ResCalc - FATAL: resolution and detInfo sizes does not match.\n";
    return -1;
  }
      
  // check <it>
  if( it >= dets_.size() ) {
    cout << "AlManager::_ResCalc - FATAL: wrong vector index " << it << ".\n";
    return -1;
  }
  
  TMatrix iV_Fast(2,2); 
  for( unsigned int i=0; i < 2; i++ )   
  for( unsigned int j=0; j < 2; j++ )   
  iV_Fast(i,j) = 0;
   
  for( unsigned int i = 0; i < dets_.size(); i++ ) {
    if( it == i ) continue;
    TMatrix iR = dets_[i]->irotM_;
    double r2 = pow( res[i], 2 );
    iV_Fast(0,0) += iR(0,0)*iR(0,0)/r2;
    iV_Fast(0,1) += iR(0,0)*iR(0,1)/r2;
    iV_Fast(1,1) += iR(0,1)*iR(0,1)/r2;
    iV_Fast(1,0) += iR(0,0)*iR(0,1)/r2;
  }
    
  TMatrix V_Fast = TMatrix( TMatrix::kInverted, iV_Fast );
  TMatrix iR = dets_[it]->irotM_;
  return sqrt( 
      V_Fast(0,0)*iR(0,0)*iR(0,0)
    + V_Fast(1,1)*iR(0,1)*iR(0,1)
    + 2*V_Fast(0,1)*iR(0,0)*iR(0,1)
  );
  
}    

