#include "PlaneRICH.h"
#include "TF1.h"

ClassImp(PlaneRICH);

//const int PlaneRICH::fMAX_MULT = 1;

PlaneRICH::PlaneRICH(const char *detname, size_t ncaths, int ncol, int nrow,
		     int center, int width)
  : Plane(detname) {


  std::cout<<"Plane: "<<detname
      <<" #cathodes: "<<ncaths
      <<" #rows: "<<nrow
      <<" #cols: "<<ncol
      <<" center: "<<center
      <<" width: "<<width<<std::endl;

  // creating Photocathodes
  for (register unsigned int i=0; i<ncaths; i++) {
    char ext[4];
    sprintf(ext,"%s%d","_",i);
    std::string catname = fName+ext;
    fCathode.push_back(new PlaneRCath(catname.c_str(),ncol,nrow,center,width));
  }
}


void PlaneRICH::Init(TTree* tree) {

  static int icout=0;
  // intialize cathodes
  for (register size_t i=0; i<fCathode.size(); i++) {
#if USE_DATABASE == 1
    fCathode[i]->setDBpt(fDataBase);
#endif
    fCathode[i]->OpenReference();
    fCathode[i]->Init(tree);
    std::vector<TH1*>& cathists=fCathode[i]->GetHistoList();
    for (register size_t j=0; j<cathists.size(); j++) {
      if (icout == 0) {
	std::cout << "Sizes: " << fCathode.size() << "  " << cathists.size() << std::endl;
	icout = 1;
      }
      AddHistogram(cathists[j]);
    }
      std::vector<Variable*>& catvars=fCathode[i]->GetVariables();
    for(size_t j=0; j<catvars.size(); j++)
      AddVariable(catvars[j]);
  }

  //book histos for rich itself ...
  std::string strtmp;
  strtmp = fName + "_UpPD";
  fH2u=new TH2F(strtmp.c_str(),"Upper Photon Detectors",288,0,288,144,0,144);
  strtmp = fName + "_LoPD";
  fH2l=new TH2F(strtmp.c_str(),"Lower Photon Detectors",288,0,288,144,0,144);
  strtmp = fName + "_Nhit";
  fH1s=new TH1F_Ref(strtmp.c_str(),"nHits Sum on each Cathode",16,0,16,fRateCounter);
  ((TH1F_Ref*)fH1s)->SetReference(fReferenceDirectory);
  strtmp = fName + "_Nhiw";
  fH1h=new TH1F(strtmp.c_str(),"nHits/ev on each Cathode",16,0,16);
  fH2u->SetOption("colz");
  fH2l->SetOption("colz");

  fH2u->SetNdivisions(-4,"X");
  fH2l->SetNdivisions(-4,"X");
  fH2u->SetNdivisions(-2,"Y");
  fH2l->SetNdivisions(-2,"Y");

  AddHistogram(fH2u);
  AddHistogram(fH2l);
  AddHistogram(fH1s);
  AddHistogram(fH1h);
}


void PlaneRICH::Reset() {
  for(register size_t i=0; i<fCathode.size(); i++) {
    fCathode[i]->Reset();
  }
}


void PlaneRICH::ResetHistograms() {
  for(register size_t i=0; i<fCathode.size(); i++) {
    fCathode[i]->ResetHistograms();
  }
}


PlaneRICH::~PlaneRICH() {
}


void PlaneRICH::StoreDigit(CS::Chip::Digit* digit) {

  std::vector<float> data=digit->GetNtupleData();
  if(data.size()>3)
    if(data[0]<fCathode.size()) {
      fCathode[(int) data[0]]->StoreDigit((int) data[1],
				    (int) data[2],
				    (int) data[3]);
    }
}


void PlaneRICH::EndEvent(const CS::DaqEvent &event) {

  static int eventn = 0;

//    size_t ct0=0,cath_h=0;
//    if (fCathode.size() != 0 ) {
//      std::vector<TH1*>& cathists=fCathode[ct0]->GetHistoList();
//      cath_h=fCathode.size()*cathists.size();
//    }

  for (register size_t i=0; i<fCathode.size(); i++) {
    fCathode[i]->UpdateRateCounter();
    fCathode[i]->EndEvent(event);

    std::vector<Variable*>& v0=((Plane*)fCathode[i])->GetVariables(); // Plane* cast is needed by gcc 2.95.2
    int nhits=((Plane*)fCathode[i])->GetNhits();
    float* row=v0[0]->GetValues();
    float* col=v0[1]->GetValues();
    float* amp=v0[2]->GetValues();

//      TH2F *fH2u=dynamic_cast<TH2F*> (fHistList[cath_h+0]);
//      TH2F *fH2l=dynamic_cast<TH2F*> (fHistList[cath_h+1]);
//      TH1F *fH1s=dynamic_cast<TH1F*> (fHistList[cath_h+2]);
//      TH1F *fH1h=dynamic_cast<TH1F*> (fHistList[cath_h+3]);

    fH1s->Fill(i,nhits);

    if((eventn%20)==0) {

      if ( i == 0 ) {
	fH2u->Reset();
	fH2l->Reset();
	fH1h->Reset();
      }

      fH1h->Fill(i,nhits);
      for (register int j=0; j<nhits; j++) {
	int    off_c=(3 - (i%8)/2)*72;
	int    off_r=(1-(i%2))*72;
	if ( i > 7 ) fH2l->Fill(col[j]+off_c,row[j]+off_r,amp[j]);
	if ( i < 8 ) fH2u->Fill(col[j]+off_c,row[j]+off_r,amp[j]);
      }
    }
  }
  eventn++;
}



void PlaneRICH::TextOutput(ostream& out)
/// get of means in hits histograms by Damien, july03
{
  for (register size_t i=0; i<fCathode.size(); i++) {
    TH1F* hhit = fCathode[i]->GetHhit();
    out << "\thitmean "<<fCathode[i]->GetName()<<" "<< i << " "<<hhit->GetMean(1)<<std::endl;
  }
  out <<"\n"<<std::endl;


/// get of slopes in amplitudes histograms by Anton Nicolas summer student, july03   

  for (register size_t i=0; i<fCathode.size(); i++) {
      
    TH1F* hamp = fCathode[i]->GetAmp();
  
    hamp->Fit("expo","V","V",10,40);
    TF1* fit1 = hamp->GetFunction("expo");
    if (fit1 != NULL) {
      Double_t pente = fit1->GetParameter(1);
      out << "\tpente "<<fCathode[i]->GetName()<<" "<< i << " "<<pente<<std::endl;
    } else {
      out << "\tpente "<<fCathode[i]->GetName()<<" *** "<<std::endl;
    }
  } 
}

void PlaneRICH::ControlPanel(const TGWindow *p, const TGWindow *main) {

  if (!fControlPanel) fControlPanel = new PlanePanel(p, main, 100, 100, this);
}








