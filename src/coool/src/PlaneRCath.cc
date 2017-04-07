#include "PlaneRCath.h"

ClassImp(PlaneRCath);

void PlaneRCath::Init(TTree* tree) {

  Plane2V::Init(tree);

  // Add Non weighted histos to Plane2V list

//    string name = fName + "_2Dhits";
//    fH2dhits = new TH2F(name.c_str(),name.c_str(),
//  		      fVcol->GetNbins(),
//  		      fVcol->GetMin(),
//  		      fVcol->GetMax(),
//  		      fVrow->GetNbins(),
//  		      fVrow->GetMin(),
//  		      fVrow->GetMax());
//    AddHistogram(fH2dhits);

  fHrc1->SetOption("col"); 
  fHrca->SetOption("col"); 
}

void PlaneRCath::EndEvent(const CS::DaqEvent &event) {

  if (thr_flag) TThread::Lock();

  //fHrca->SetOption("scat"); 
  
  for (int i=0; i<fNhits; i++) {
    
    int row=fRow[i];
    int col=fCol[i];
    int amp=fAmp[i];
 
    if(fVrow->Test(row) &&
       fVcol->Test(col) &&
       fVamp->Test(amp)) {    

      fHa->Fill(amp);
      fHrc1->Fill(col,row);
      fHrca->Fill(col,row,amp);

      fVrow->Store(row);
      fVcol->Store(col);
      fVamp->Store(amp);

      fNhitsKept++;
    }
  }
  fHhit->Fill(fNhitsKept);
  
  if (thr_flag) TThread::UnLock();
}






