#include "Plane2V.h"
#include "TProfile.h"
#include "Plane2VPanel.h"

ClassImp(Plane2V);

const int Plane2V::fMAX_MULT = 8;

Plane2V::Plane2V(const char *detname,int ncol, int nrow, int center, int width)
  : Plane(detname),fNrows(nrow),fNcols(ncol),fNchan(nrow*ncol) {

  fRow = new int[fNchan*fMAX_MULT];
  fCol = new int[fNchan*fMAX_MULT];
  fAmp = new int[fNchan*fMAX_MULT];
  std::string rowname=fName+"_row";
  std::string colname=fName+"_col";
  std::string ampname=fName+"_amp";

  fVrow=AddVariable(rowname.c_str(),fNrows,0,fNrows ,fNchan*fMAX_MULT);
  fVcol=AddVariable(colname.c_str(),fNcols,0,fNcols ,fNchan*fMAX_MULT);
  fVamp=AddVariable(ampname.c_str(),100,center-width,center+width,
		    fNchan*fMAX_MULT);
}

void Plane2V::Init(TTree* tree) {

  // 0 multiplicity within the time cut on the whole detector
  std::string hitsname = fName + "_hits";
  fHhit=new TH1F_Ref(hitsname.c_str(),hitsname.c_str(),500,0,500, fRateCounter);
  ((TH1F_Ref*)fHhit)->SetReference(fReferenceDirectory);
  AddHistogram(fHhit);
  std::string hitsleavlist = hitsname + "/I";

  // 1 histogram : row % col, amp weight
  std::string name01 = fVrow->GetName() + ":" + fVcol->GetName() + "_amp";
  fHrca=new TH2F(name01.c_str(),name01.c_str(),
		 fVcol->GetNbins(),
		 fVcol->GetMin(),
		 fVcol->GetMax(),
		 fVrow->GetNbins(),
		 fVrow->GetMin(),
		 fVrow->GetMax());
  fHrca->SetOption("box");
  AddHistogram(fHrca);

  // histogram : row % col,
  std::string name02 = fVrow->GetName() + ":" + fVcol->GetName();
  fHrc1=new TH2F(name02.c_str(),name02.c_str(),
		 fVcol->GetNbins(),
		 fVcol->GetMin(),
		 fVcol->GetMax(),
		 fVrow->GetNbins(),
		 fVrow->GetMin(),
		 fVrow->GetMax());
  fHrc1->SetOption("box");
  AddHistogram(fHrc1);

  // 2 histogram : amplitude all pads
  std::string ampname = fName + "_amp";
  fHa=new TH1F_Ref(fVamp->GetName().c_str(),fVamp->GetName().c_str(),
	       fVamp->GetNbins(),
	       fVamp->GetMin(),
	       fVamp->GetMax(), fRateCounter);
  ((TH1F_Ref*)fHa)->SetReference(fReferenceDirectory);
  AddHistogram(fHa);

  // 3 row vs column
  std::string avsadr = fName + "_adrVSamp";
  fHavsadr = new TH2F(avsadr.c_str(), avsadr.c_str(),
		fNrows*fNcols, 0, fNrows*fNcols,
		fVamp->GetNbins(), fVamp->GetMin(),
		fVamp->GetMax());
  AddHistogram(fHavsadr);

  std::string leavlist0 =fVrow->GetName()  + "[" + hitsname + "]/F";
  std::string leavlist1 =fVcol->GetName()  + "[" + hitsname + "]/F";
  std::string leavlist2 =fVamp->GetName()  + "[" + hitsname + "]/F";

  if(tree) {
    fIsInTree = true;
    tree->Branch(hitsname.c_str(),&fNhitsKept,
		 hitsleavlist.c_str(),32000);
    tree->Branch(fVrow->GetName().c_str(),fVrow->GetValues(),
		 leavlist0.c_str(),32000);
    tree->Branch(fVcol->GetName().c_str(),fVcol->GetValues(),
		 leavlist1.c_str(),32000);
    tree->Branch(fVamp->GetName().c_str(),fVamp->GetValues(),
		 leavlist2.c_str(),32000);

  }
}

Plane2V::~Plane2V() {
  delete fRow; delete fCol; delete fAmp;
}

void Plane2V::StoreDigit(int col, int row, int amp) {

  if (fNhits < fNchan*fMAX_MULT) {
    fRow[fNhits]=row;
    fCol[fNhits]=col;
    fAmp[fNhits]=amp;
    fNhits++;
  }
}

void Plane2V::StoreDigit(CS::Chip::Digit* digit) {
  std::vector<float> data=digit->GetNtupleData();
  if(data.size()>2)
    this->StoreDigit((int) data[0],(int) data[1],(int) data[2]);
}

void Plane2V::EndEvent(const CS::DaqEvent &event) {

  if (thr_flag) TThread::Lock();

  for (register int i=0; i<fNhits; i++) {

    register int row=fRow[i];
    register int col=fCol[i];
    register int amp=fAmp[i];

    if(fVrow->Test(row) &&
       fVcol->Test(col) &&
       fVamp->Test(amp)) {

      register int adr=fRow[i] + fCol[i]*fNrows;

      fHrca->Fill(col,row,amp);
      fHavsadr->Fill(adr,amp);
      fHrc1->Fill(col,row);
      fHa->Fill(amp);

      fVrow->Store(row);
      fVcol->Store(col);
      fVamp->Store(amp);

      fNhitsKept++;
    }
  }
  fHhit->Fill(fNhitsKept);

  if (thr_flag) TThread::UnLock();
}

void Plane2V::ControlPanel(const TGWindow *p, const TGWindow *main) {

  if (!fControlPanel) fControlPanel = new Plane2VPanel(p, main, 100, 100, this);
}

void Plane2V::AmpSpectrum() {

  std::string hdname = fName + "_1ch";
  if(fCurChan<fNchan) {
    //fHistList.pop_back();
    fHavsadr->ProjectionY(hdname.c_str(),fCurChan+1,fCurChan+1,"")->Draw();
  }
}





