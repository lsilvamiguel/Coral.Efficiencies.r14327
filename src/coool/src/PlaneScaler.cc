#include "PlaneScaler.h"
#include "PlanePanel.h"
#include "TThread.h"
#include "Riostream.h"
#include <string>

ClassImp(PlaneScaler);

const int PlaneScaler::fNchan = 32;


PlaneScaler::PlaneScaler(const char *detname, int ebe) 
  : Plane(detname), fIsEbe(ebe),  fBeamhitspspill(0), fSpillnum(0) {

  fChanMult = new int[PlaneScaler::fNchan]; 
  fLastCounts = new int[PlaneScaler::fNchan]; 
  fCounts = new int[PlaneScaler::fNchan]; 

  for(int i=0; i<PlaneScaler::fNchan; i++) {
    fChanMult[i]=0;
    fLastCounts[i]=0;
    fCounts[i]=0;  
  }
  fTime = 0;
  fLastTime = 0;
}

void PlaneScaler::Reset() {
  for(int i=0; i<PlaneScaler::fNchan; i++) {
    fChanMult[i]=0;
    fLastCounts[i]=fCounts[i];
    fCounts[i]=0;
  }
  fLastTime = fTime;
}

void PlaneScaler::NewRun() {
  // it is meaningless to go on filling this histo when starting a new run...
  if (thr_flag) TThread::Lock();
  fHIntensity->Reset();   
  if (thr_flag) TThread::UnLock();
}

void PlaneScaler::End() {
  // fill last spill intensity in fHIntensity
  if(fHIntensity) 
    fHIntensity -> Fill(fSpillnum, fBeamhitspspill);
}

void PlaneScaler::Init(TTree* tree) {
  
  fRateCounter = 0;

  std::string cname = fName + "_counts";
  fHcounts = new TH1F_Ref(cname.c_str(),cname.c_str(),
		      PlaneScaler::fNchan,0,
                      PlaneScaler::fNchan,fRateCounter);
  ((TH1F_Ref*)fHcounts)->SetReference(fReferenceDirectory);
  AddHistogram(fHcounts);

  std::string tname = fName + "_times";
  fHtimes = new TH1F(tname.c_str(),tname.c_str(),
		      60,0,6);
  AddHistogram(fHtimes);


  std::string rname = fName + "_rates";
  fHrates = new TProfile(rname.c_str(),rname.c_str(),
		     PlaneScaler::fNchan,0,PlaneScaler::fNchan,
		     0,1e12);
  AddHistogram(fHrates);

  std::string sname = fName + "_spillprof";
  fHSpillProfile = new TProfile(sname.c_str(),sname.c_str(),
				60,0,6,
				0,1e8);
  AddHistogram(fHSpillProfile);
    
  std::string iname = fName + "_intensity"; 
  fHIntensity = new TH1F(iname.c_str(),iname.c_str(),
			 200,1,201);
  AddHistogram(fHIntensity);
    
  if(tree) {
    std::string mname = fName + "_mults";
    fIsInTree = true;

    std::string mleavlist=mname + "[32]/I";    
    std::string tname=fName + "_t";
    std::string tleavlist=tname + "/I";
    std::string dtname=fName + "_dt";
    std::string dtleavlist=tname + "/F";
    std::string pname=fName + "_pattern";
    std::string pleavlist=pname + "/I";

    //      tree->Branch(mname.c_str(),fChanMult,"c01/F:c02:c03:c04:c05:c06:c07:c08:c09:c10:c11:c12:c13:c14:c15:c16:c17:c18:c19:c20:c21:c22:c23:c24:c25:c26:c27:c28:c29:c30:c31:c32",32000);
    tree->Branch(mname.c_str(),fChanMult,mleavlist.c_str(),32000);
    tree->Branch(tname.c_str(),&fTime,tleavlist.c_str(),32000);
    tree->Branch(dtname.c_str(),&fDeltaTime,dtleavlist.c_str(),32000);
    tree->Branch(pname.c_str(),&fPattern,pleavlist.c_str(),32000);
  }
}


PlaneScaler::~PlaneScaler() {
  delete fChanMult;
  delete fLastCounts;
  delete fCounts;
}


void PlaneScaler::StoreDigit(CS::Chip::Digit* digit) {

  std::vector<float> data=digit->GetNtupleData();

  if(data.size() != 2)
    throw CS::Exception("PlaneScaler::StoreDigit(): Scaler data not recognized.");
  
  int chan = static_cast<int>(data[0]);
  int counts = static_cast<int>(data[1]);

  //cout<<fName<<" "<<chan<<" "<<counts<<endl;

  if(chan <= 31)
    fCounts[chan] = counts;
  else if(chan == 33) // time
    fTime = counts;
  else if(chan == 32) // pattern
     fPattern = counts;
   else 
    cout << "Wrong mapping of Scaler: " << fName.c_str() << "  Channel:  " << chan << "  is not allowed" << endl;
}

void PlaneScaler::EndEvent(const CS::DaqEvent &event) {
// cerr<<"PlaneScaler::EndEvent fRateCounter "<<fRateCounter<<endl;

  if (thr_flag) TThread::Lock();
  fDeltaTime = (fTime-fLastTime)/38.88E6;
  int beamhits=0;
  
  if(!fSpillnum) fSpillnum = event.GetBurstNumber();

  for(int i=0; i<PlaneScaler::fNchan; i++) {
    int ch=i;

    fChanMult[ch]=fCounts[i]-fLastCounts[i];

    if (fDeltaTime>0) {
      fHcounts->Fill(ch,fChanMult[ch]);
      fHrates->Fill(ch,fChanMult[ch]/fDeltaTime);    
      beamhits+=fChanMult[ch];
    }
  }
  if (fHSpillProfile && fDeltaTime>0) {
    fHSpillProfile->Fill(fTime/38.88E6,beamhits/fDeltaTime);
  }
  
  if(fDeltaTime<0) { // first event after a reset (ie new spill) 
    fHIntensity -> Fill(fSpillnum, fBeamhitspspill);
    fBeamhitspspill = 0;
    fSpillnum = event.GetBurstNumber();
  } else { // still in the same spill - go on counting
    fBeamhitspspill += beamhits;
  }

  fHtimes->Fill(fTime/38.88E6);
  if (thr_flag) TThread::UnLock();
}

void PlaneScaler::ControlPanel(const TGWindow* p, const TGWindow* main) {

  if (!fControlPanel) fControlPanel = new PlanePanel(p, main, 100, 100, this);
}











