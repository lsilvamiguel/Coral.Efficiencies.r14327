#include "Fits.h"
#include "TF1.h"
#include <iostream>

void FitCTime(TH1F* Data,float &sigma, float &mean) {

  //STEP 1: Generates theoretical function
  Double_t params[4] = {0,0,0,0};
  TF1 *theory = new TF1("theory","pol0(0)+gaus(1)",-400,400);
  
  //STEP 3: Estimates background parameters using a gaussian
  
  Data->Fit("pol0","q0","",-400,-250);
  
  //STEP 4: Subtract estimated background to original data
  // Creates a temporary histogram and fit a gaussian
  TH1F *htemp = (TH1F*)Data->Clone();
  htemp->Reset();
  TF1 *eback = Data->GetFunction("pol0");
  for (Int_t bin=1;bin<=Data->GetNbinsX();bin++) {
    Float_t x = Data->GetBinCenter(bin);
    Double_t fval = eback->Eval(x);
    Double_t diff = Data->GetBinContent(bin)-fval;
    if(diff>0)
      htemp->Fill(x,diff);
    else 
      htemp->Fill(x,0);      
  }
  htemp->Fit("gaus","q0");
  TF1 *esig = htemp->GetFunction("gaus");
  
  //STEP 5: Fit background + signal
  eback->GetParameters(&params[0]);
  esig->GetParameters(&params[1]);
  theory->SetParameters(params);
  theory->SetParLimits(0,params[0]-20,params[0]+20);
  Data->Fit("theory","BRq0");
  
  sigma=theory->GetParameter(3);
  mean=theory->GetParameter(2);

}


void FitResiduals(TH1F* Data, float &sigma, float &mean) {
  //STEP 1: Generates theoretical function
  Double_t params[4] = {0,0,0,0};
  TF1 *theory = new TF1("theory","pol0(0)+gaus(1)",
			Data->GetXaxis()->GetXmin(),
			Data->GetXaxis()->GetXmax());
   
  //STEP 3: Estimates background parameters using a gaussian

  Data->Fit("pol0","q0","",0.2,1.5);
  
  //STEP 4: Subtract estimated background to original data
  // Creates a temporary histogram and fit a gaussian
  TH1F *htemp = (TH1F*)Data->Clone();
  htemp->Reset();
  TF1 *eback = Data->GetFunction("pol0");
  for (Int_t bin=1;bin<=Data->GetNbinsX();bin++) {
    Float_t x = Data->GetBinCenter(bin);
    Double_t fval = eback->Eval(x);
    Double_t diff = TMath::Abs(fval - Data->GetBinContent(bin));
    htemp->Fill(x,diff);
  }
  htemp->Fit("gaus","q0");
  TF1 *esig = htemp->GetFunction("gaus");

  //STEP 5: Fit background + signal
  eback->GetParameters(&params[0]);
  esig->GetParameters(&params[1]);
  theory->SetParameters(params);
  Data->Fit("theory","q");
  Data->GetFunction("theory")->SetLineColor(5);

  sigma=theory->GetParameter(3);
  mean=theory->GetParameter(2);
}



