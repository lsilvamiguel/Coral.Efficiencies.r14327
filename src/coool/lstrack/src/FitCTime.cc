#include "Fits.h"
#include "TF1.h"

void FitCTime(TH1F* Data, float &sigma, float &mean, float min,
	      float max,float minbg,float maxbg) {
  //STEP 1: Generates theoretical function
  const Int_t npar = 4;
  Double_t params[npar] = {0,0,0,0};
  TF1 *theory = new TF1("theory","pol0(0)+gaus(1)",min,max);
   
  //STEP 3: Estimates background parameters using a gaussian

  Data->Fit("pol0","q0","",minbg,maxbg);
  
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
  Data->Fit("theory","q0");

  sigma=theory->GetParameter(3);
  mean=theory->GetParameter(2);
}


