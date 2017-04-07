#include "PlaneVBOX.h"

ClassImp(PlaneVBOX);

void PlaneVBOX::Init(TTree* tree) {
  
  PlaneHCAL1::Init(tree);
  

  lAmpHistOffset=this->GetHistoList().size();
  std::cout<<"PlaneVBOX: fNrows: "<<fNrows<<", fNcols:"<<fNcols<<", histOfs: "<<lAmpHistOffset<<std::endl;


  for (int col=0; col < fNcols; col++){
    for (int row=0; row < fNrows; row++){
#if __GNUC__ > 3 || __GNUC__ == 3
      std::ostringstream nametmp;
#else
      ostrstream nametmp;
#endif

      nametmp<<fName << "_tamp_" << col << "-" << row << std::ends;
      std::string name = nametmp.str();
      TH1F *fHlocal=new TH1F_Ref(name.c_str(),name.c_str(),
                      fVamp->GetNbins(), fVamp->GetMin(),
                       fVamp->GetMax(), fRateCounter);
      AddHistogram(fHlocal);
    }
  }
  std::cout<<"PlaneVBOX: fNrows: "<<fNrows<<", fNcols:"<<fNcols<<", histoListSize: "<<this->GetHistoList().size()<<std::endl;

}

void PlaneVBOX::EndEvent(const CS::DaqEvent &event) {

  if (thr_flag) TThread::Lock();

  PlaneHCAL1::EndEvent(event);

  for (int i=0; i < fNhits; i++) {
    int row=fRow[i];
    int col=fCol[i];
    int amp=fAmp[i];

    if (fVrow->Test(row) &&
        fVcol->Test(col)
       ) {
//	fVamp->Test(amp)) {
     
//  std::cout<<"PlaneVBOX: Fill "<<col<<","<<row<<";"<<amp<<std::endl;
	   (dynamic_cast<TH1F*> (this->GetHistoList()[lAmpHistOffset + col*fNrows + row]))->Fill(amp); 
    }
  }

  
  if (thr_flag) TThread::UnLock();
}


