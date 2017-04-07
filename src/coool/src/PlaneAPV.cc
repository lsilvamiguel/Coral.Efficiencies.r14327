#include "PlaneAPV.h"
#include "PlanePanel.h"

ClassImp(PlaneAPV);

const int PlaneAPV::fMAX_MULT = 1;

const CS::Chip::Maps* PlaneAPV::fRunMaps = 0;

PlaneAPV::PlaneAPV(const char *detname,int nchan, int center, int width, bool pixel)
  : Plane(detname),fNchan(nchan) {

  fPixel = pixel;

  fVch = AddVariable("_ch",fNchan,0,fNchan,fNchan*fMAX_MULT);
  fVa0 = AddVariable("_a0",500,0,500,fNchan*fMAX_MULT);
  fVa1 = AddVariable("_a1",500,0,500,fNchan*fMAX_MULT);
  fVa2 = AddVariable("_a2",500,0,500,fNchan*fMAX_MULT);
  fVt = AddVariable("_t",100,0,30000,fNchan*fMAX_MULT);
}

void PlaneAPV::Init(TTree* tree) {

  fRateCounter = 0;

  // Book Histograms for single strips, create tree leaves

  // Number of hit strips
  std::string hitsname = fName + "_hitMultiplicity";
  fHhits=new TH1F_Ref(hitsname.c_str(),hitsname.c_str(),
                      300,-0.5,299.5,fRateCounter);
  ((TH1F_Ref*)fHhits)->SetReference(fReferenceDirectory);
  AddHistogram(fHhits);
  std::string hitsleavlist = hitsname + "/I";

  // Hit profile
  std::string chname;
  std::string chleavlist;
  
  chname = fName + fVch->GetName();
  fHch=new TH1F_Ref(chname.c_str(),chname.c_str(),
          fVch->GetNbins(),
		  fVch->GetMin(),
		  fVch->GetMax(),fRateCounter);	
  ((TH1F_Ref*)fHch)->SetReference(fReferenceDirectory);
  AddHistogram(fHch);
  chleavlist = chname + "[" + hitsname +"]/F";

  // Amplitude 1st sample
  std::string a0name = fName + fVa0->GetName();
  fHa0=new TH1F_Ref(a0name.c_str(),a0name.c_str(),
		   fVa0->GetNbins(),
		   fVa0->GetMin(),
		   fVa0->GetMax(),fRateCounter);
  ((TH1F_Ref*)fHa0)->SetReference(fReferenceDirectory);
  AddHistogram(fHa0);
  std::string a0leavlist = a0name + "[" + hitsname +"]/F";

  // Amplitude 2nd sample
  std::string a1name = fName + fVa1->GetName();
  fHa1=new TH1F_Ref(a1name.c_str(),a1name.c_str(),
		fVa1->GetNbins(),
		fVa1->GetMin(),
		fVa1->GetMax(),fRateCounter);
  ((TH1F_Ref*)fHa1)->SetReference(fReferenceDirectory);
  AddHistogram(fHa1);
  std::string a1leavlist = a1name + "[" + hitsname +"]/F";

  // Amplitude 3rd sample
  std::string a2name = fName + fVa2->GetName();
  fHa2=new TH1F_Ref(a2name.c_str(),a2name.c_str(),
		   fVa2->GetNbins(),
		   fVa2->GetMin(),
		   fVa2->GetMax(),fRateCounter);
  ((TH1F_Ref*)fHa2)->SetReference(fReferenceDirectory);
  AddHistogram(fHa2);
  std::string a2leavlist = a2name + "[" + hitsname +"]/F";

  // Time from single strip amplitudes
  std::string tname = fName + fVt->GetName();
  fHt=new TH1F(tname.c_str(),tname.c_str(),
		   fVt->GetNbins(),
		   fVt->GetMin(),
		   fVt->GetMax());
  AddHistogram(fHt);
  std::string tleavlist = tname + "[" + hitsname +"]/F";

  if(tree) {
    fIsInTree = true;
    tree->Branch(hitsname.c_str(),&fNhitsKept,
		 hitsleavlist.c_str(),32000);
    tree->Branch(chname.c_str(),fVch->GetValues(),
         chleavlist.c_str(),32000);
    tree->Branch(a0name.c_str(),fVa0->GetValues(),
		 a0leavlist.c_str(),32000);
    tree->Branch(a1name.c_str(),fVa1->GetValues(),
		 a1leavlist.c_str(),32000);
    tree->Branch(a2name.c_str(),fVa2->GetValues(),
		 a2leavlist.c_str(),32000);
    tree->Branch(tname.c_str(),fVt->GetValues(),
		 tleavlist.c_str(),32000);
  }
}


PlaneAPV::~PlaneAPV() {}

void PlaneAPV::ResetHistograms() {
  Plane::ResetHistograms();
}


// for pixels of pixelGEM like detectors
void PlaneAPV::StoreDigit(int x, int y, int channel, int amp0, int amp1, int amp2, int cmc0, int cmc1, int cmc2, int chip) {
  if (fNhits==0) fDigits.clear();
  float a12ratio = (amp2==0) ? 0 : float(amp1)/float(amp2);
  float a02ratio = (amp2==0) ? 0 : float(amp0)/float(amp2);
  std::vector<double> data;
  data.push_back(x);
  data.push_back(y);
  data.push_back(amp0);
  data.push_back(amp1);
  data.push_back(amp2);
  data.push_back(a12ratio);
  data.push_back(a02ratio);
  data.push_back(cmc0);
  data.push_back(cmc1);
  data.push_back(cmc2);
  data.push_back(chip);
  //data.push_back(a02ratio);
  fDigits.insert(make_pair(channel,CDigit1(channel,data)));
  fNhits++;
}


// for strips
void PlaneAPV::StoreDigit(int channel, int amp0, int amp1, int amp2, int cmc0, int cmc1, int cmc2, int wirepos, int chip) {
  if (fNhits==0) fDigits.clear();
  float a12ratio = (amp2==0) ? 0 : float(amp1)/float(amp2);
  float a02ratio = (amp2==0) ? 0 : float(amp0)/float(amp2);
  std::vector<double> data;
  data.push_back(amp0);
  data.push_back(amp1);
  data.push_back(amp2);
  data.push_back(a12ratio);
  data.push_back(a02ratio);
  data.push_back(cmc0);
  data.push_back(cmc1);
  data.push_back(cmc2);
  data.push_back(a02ratio);
  data.push_back(wirepos);
  data.push_back(chip);
  fDigits.insert(make_pair(channel,CDigit1(channel,data)));
  fNhits++;
}


void PlaneAPV::StoreDigit(CS::Chip::Digit* digit) {
  std::vector<float> data=digit->GetNtupleData();
  if (fPixel) {
    if (data.size() >= 8) {
      this->StoreDigit((int)data[0],  (int)data[1],  (int)data[2],
                       (int)data[5],  (int)data[6],  (int)data[7],
                       (int)data[12], (int)data[13], (int)data[14],
                       (int)data[3]);
    }
  } else {
    if(data.size() >= 6) {
      this->StoreDigit((int)data[0],(int)data[3],(int)data[4],(int)data[5],
	               (int)data[10],(int)data[11],(int)data[12],(int)data[13],(int)data[1]); 
    }
  }
}

void PlaneAPV::EndEvent(const CS::DaqEvent &event) {

  if (thr_flag) TThread::Lock();

  for (std::map<int,CDigit1>::const_iterator id = fDigits.begin();
       id!=fDigits.end(); ) {

    const CDigit1& digit = id->second;

    int channel=digit.ch; //fChannel[i];
    int a0=(int)digit.dt[0];   //fAmp1[i];
    int a1=(int)digit.dt[1];   //fAmp2[i];
    int a2=(int)digit.dt[2];   //fAmp3[i];
    int a12ratio =(int)digit.dt[3];
    int a02ratio =(int)digit.dt[4];
    int cmc0 =(int)digit.dt[5];
    int cmc1 =(int)digit.dt[6];
    int cmc2 =(int)digit.dt[7];

    printf("APV data: %5d %5d %5d %5d %5d %5d %5d %5d %5d"
	      "Normally this method is overloaded by detector classes - "
	      "how did you get here?\n",
	      channel, a0, a1, a2, a12ratio, a02ratio, cmc0, cmc1, cmc2);
  }
  if (thr_flag) TThread::UnLock();
}

void PlaneAPV::ControlPanel(const TGWindow* p, const TGWindow* main) {
    if (!fControlPanel) fControlPanel = new PlanePanel(p, main, 100, 100, this);
  }





