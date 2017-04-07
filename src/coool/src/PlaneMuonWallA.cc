#include "PlaneMuonWallA.h"
#include "ChipF1.h"
#include "TriggerTime.h"
#include "TF1.h"
#include "TFile.h"
#include "TThread.h"
#include "PlanePanel.h"
#include "Reference.h"

ClassImp(PlaneMuonWallA);

#define TMIN -11000
#define TMAX  -7000
#define DT   2000

const int   PlaneMuonWallA::fMAX_MULT = 8;
const float PlaneMuonWallA::fF1_TICK  = 128.e-9; // is it so??:?
const int   PlaneMuonWallA::fRATE_UPDATE = 100;

const float PlaneMuonWallA::fHIT_TMIN = TMIN - DT;
const float PlaneMuonWallA::fHIT_TMAX = TMAX + DT;
const float PlaneMuonWallA::fMT_T0    = (fHIT_TMIN+fHIT_TMAX)/2;
const float PlaneMuonWallA::fCL_TMIN  =  - DT; // temporary
const float PlaneMuonWallA::fCL_TMAX  =  + DT; // temporary

TFile* PlaneMuonWallA::refFile = NULL;

bool PlaneMuonWallA::OpenRefFile(const char* fileName)
{
  refFile = new TFile(fileName);
  if(!refFile || !refFile->IsOpen())
    {
      std::cout << "Cannot open reference file: " << fileName << std::endl;
      return false;
    }
  return true;
}

void PlaneMuonWallA::CloseRefFile()
{
  if(refFile)
    delete refFile;
}


PlaneMuonWallA::PlaneMuonWallA(const char *detname, int nchan, int center, int width)
  :
  Plane(detname) {
//   fIsHR = false;
  int fNdet=800;
  if(strstr(detname,"X")) {
    fNchan = 456;
    fNdet  = 57;
  }
  else if(strstr(detname,"Y")) {
    fNchan = 392;
    fNdet  = 49;
  }
//    cerr<<"Plane: "<<detname
//        <<" #channels: "<<nchan
//        <<" center: "<<center
//        <<" width: "<<width<<std::endl;
  // fVch used now to store Digits and print hits for efiiciency estimation
  // It is not go to histograms
  fVch    = AddVariable("_ch",fNchan,0,fNchan,fNchan*fMAX_MULT);
  fVch_up = AddVariable("_ch_up",fNchan,0,fNchan,fNchan*fMAX_MULT);
  fVch_down = AddVariable("_ch_down",fNchan,0,fNchan,fNchan*fMAX_MULT);
  fVdet_up = AddVariable("_det_up",    fNdet,0,fNdet,fNdet*fMAX_MULT);
  fVdet_down = AddVariable("_det_down",fNdet,0,fNdet,fNdet*fMAX_MULT);
  fVt = AddVariable("_t",100,center-width,center+width,fNchan*fMAX_MULT);

//   fChannel=new int[fNchan*fMAX_MULT]; // Used to store data for efficiency
//   fTime = new int[fNchan*fMAX_MULT];  // estimation
//   fKey  = new int[fNchan*fMAX_MULT];

}

void PlaneMuonWallA::Init(TTree *tree) {
  fRateCounter = 0;

  float posmax = 20;
  if(fGeom)
    posmax = fGeom->GetHalfSize();

  //  fVcch = AddVariable("_cch",fNchan,-posmax,posmax,fNchan);
  //  fVct = AddVariable("_ct",100,-DT,DT, fNchan);

  //multiplicity within the time cut on the whole detector
  std::string hitsname = fName + "_hits";
  fHhits=new TH1F_Ref(hitsname.c_str(),hitsname.c_str(),10,0,10,fRateCounter);
  AddHistogram(fHhits);
  std::string hitsleavlist = hitsname + "/I";

  //branch : hit map
  std::string chuname = fName + fVch_up->GetName();
  fHch_up=new TH1F_Ref(chuname.c_str(),chuname.c_str(),
		       fVch_up->GetNbins(),
		       fVch_up->GetMin(),
		       fVch_up->GetMax(),
		       fRateCounter);
  AddHistogram(fHch_up);
  std::string chuleavlist = chuname + "[" + hitsname +"]/F";

  std::string chdname = fName + fVch_down->GetName();
  fHch_down=new TH1F_Ref(chdname.c_str(),chdname.c_str(),
		     fVch_down->GetNbins(),
		     fVch_down->GetMin(),
		     fVch_down->GetMax(),
		     fRateCounter);
  AddHistogram(fHch_down);
  std::string chdleavlist = chdname + "[" + hitsname +"]/F";

  // groups by 8 wires
  std::string duname = fName + fVdet_up->GetName();
  fHdet_up=new TH1F_Ref(duname.c_str(),duname.c_str(),
                fVdet_up->GetNbins(),
                fVdet_up->GetMin(),
                fVdet_up->GetMax(),
		fRateCounter);
  AddHistogram(fHdet_up);
  std::string duleavlist = duname + "[" + hitsname +"]/F";

  std::string ddname = fName + fVdet_down->GetName();
  fHdet_down=new TH1F_Ref(ddname.c_str(),ddname.c_str(),
                fVdet_down->GetNbins(),
                fVdet_down->GetMin(),
                fVdet_down->GetMax(),
		fRateCounter);
  AddHistogram(fHdet_down);
  std::string ddleavlist = ddname + "[" + hitsname +"]/F";

  //branch : time
  std::string tname = fName + fVt->GetName();
  fHt=new TH1F_Ref(tname.c_str(),tname.c_str(),
		   fVt->GetNbins(),
		   fVt->GetMin(),
		   fVt->GetMax(),
		   fRateCounter);
  AddHistogram(fHt);
  std::string tleavlist = tname + "[" + hitsname +"]/F";

  //Cluster multiplicity
//    int ccol=8; // color for cluster histos
//    std::string chitsname = fName+"_chits";
//    std::string title = fName+" Cluster Multiplicities";
//    fHchits=new TH1F(chitsname.c_str(),title.c_str(),11,-0.5,10.5);
//    fHchits->SetFillColor(ccol);
//    AddHistogram(fHchits);
//    std::string chitsleavlist = chitsname + "/I";

//    //Cluster profile
//    std::string cchname = fName+"_cch";
//    title = fName+" Cluster Profile";
//    //fHcch=new TH1F(cchname.c_str(),title.c_str(),2*fNchan-1,0,fNchan);
//    fHcch=new TH1F(cchname.c_str(),title.c_str(),100,-posmax,posmax);
//    fHcch->SetFillColor(ccol);
//    AddHistogram(fHcch);
//    std::string cchleavlist = cchname + "[" + chitsname +"]/F";

//    //Cluster time
//    std::string cd1name = fName+"_ct"; title = fName+" Cluster Time";
//    fHct=new TH1F(cd1name.c_str(),title.c_str(),100, -2000., +2000.);
//    fHct->SetFillColor(ccol);
//    AddHistogram(fHct);
//    std::string cd1leavlist = cd1name + "[" + chitsname +"]/F";

  if(tree) {
    fIsInTree = true;
    tree->Branch(hitsname.c_str(),&fNhitsKept,
                 hitsleavlist.c_str(),32000);
    tree->Branch(chuname.c_str(),fVch_up->GetValues(),
                 chuleavlist.c_str(),32000);
    tree->Branch(chdname.c_str(),fVch_down->GetValues(),
                 chdleavlist.c_str(),32000);
    tree->Branch(duname.c_str(),fVdet_up->GetValues(),
                 duleavlist.c_str(),32000);
    tree->Branch(ddname.c_str(),fVdet_down->GetValues(),
                 ddleavlist.c_str(),32000);
    tree->Branch(tname.c_str(),fVt->GetValues(),
                 tleavlist.c_str(),32000);
  }
  ((TH1F_Ref*)fHhits)->SetReference(fReferenceDirectory);
  ((TH1F_Ref*)fHch_down)->SetReference(fReferenceDirectory);
  ((TH1F_Ref*)fHch_up)->SetReference(fReferenceDirectory);
  ((TH1F_Ref*)fHdet_down)->SetReference(fReferenceDirectory);
  ((TH1F_Ref*)fHdet_up)->SetReference(fReferenceDirectory);
  ((TH1F_Ref*)fHt)->SetReference(fReferenceDirectory);
  CloseRefFile();
}


PlaneMuonWallA::~PlaneMuonWallA() {
//   delete fChannel; delete fTime; delete fKey;
}

void PlaneMuonWallA::StoreDigit(int channel, int time) {
  std::cerr<<"Plane1V::StoreDigit(int,...): Error ! this method is obsolete !!!\n";
//   if (fNhits == 0) fDigits.clear();
//   if (fNhits < fNchan*fMAX_MULT) {
//     fTime[fNhits]=time;
//     fChannel[fNhits]=channel;
//     fNhits++;
//   }
}

void PlaneMuonWallA::StoreDigit(int channel, const std::vector<double>& data) {
  std::cerr<<"Plane1V::StoreDigit(int,...) ("<<GetName()<<"): Error ! this method is obsolete !!!\n";
//   if (fNhits == 0) fDigits.clear();

  //  fDigits.insert(make_pair(channel,CDigit1(channel,data)));
//   float time = data[0];
//   if (fNhits < fNchan*fMAX_MULT) {
//     fTime[fNhits]=int(time);
//     fChannel[fNhits]=channel;
//     fKey[fNhits]= int(data[1]);
//    fNhits++;
//   }
}

void PlaneMuonWallA::StoreDigit(CS::Chip::Digit* digit) {
  lDigits.push_back(digit);
  if (fNhits == 0) fDigits.clear();
  fNhits++;
//   std::vector<float> data=digit->GetNtupleData();               // ch, time, pos_key, time_unit
//   std::vector<double> d;
//   d.push_back(data[1]); // time
//   d.push_back(data[2]); // signature for wire: 0/1/-1
//   this->StoreDigit((int)data[0], d);
//   CS::ChipF1::Digit *df1=dynamic_cast<CS::ChipF1::Digit*>(digit);
//   if( df1!=NULL )
//     if( fabs(df1->GetTimeUnit()-CS::ChipF1::GetUnitHigh())<1e-6 )
//       fIsHR = true;
//     else
//     if( fabs(df1->GetTimeUnit()-CS::ChipF1::GetUnitNorm())<1e-6 )
//       fIsHR = false;
}


void PlaneMuonWallA::ResetHistograms() {
  Plane::ResetHistograms();
}


void PlaneMuonWallA::EndEvent(const CS::DaqEvent &event) {
  // const float reftime = event.GetTT().GetTimeNorm(); // may be later TSC
//   const float reftime = fIsHR ? event.GetTT().GetTimeHigh() :
//                                   event.GetTT().GetTimeNorm();

  if (thr_flag) TThread::Lock();

//   for (int i=0; i<fNhits; i++) {
//     int time = int(CS::ChipF1::TimeDifference(fTime[i],reftime));
//     int channel=fChannel[i];
  typedef std::list<CS::Chip::Digit*>::iterator lDIter;
  for (lDIter ii = lDigits.begin(); ii != lDigits.end(); ii++) {
    CS::ChipF1::Digit* iii = dynamic_cast<CS::ChipF1::Digit*> (*ii);
    if (!iii) {
      std::cerr<<"PlaneMuonWallA::EndEvent ("<<GetName()<<"): a digit is not a F1 one, strange...\n";
      continue;
    }
    register int time = (int) (iii->GetTimeDecoded() / iii->GetTimeUnit());
    register int channel = iii->GetChannel();
    register int key = iii->GetChannelPos();

    if( fVch->Test(channel) ) {
      if( fVt->Test(time) ) {
        fVch->Store(channel);
        fVt->Store(time);
        //fHch->Fill(channel);
        fHt->Fill(time);
	fNhitsKept++;
// 	if(fKey[i] >= 0) {
	if(key >= 0) {
	  fHch_up->Fill(channel);
	  fHdet_up->Fill(channel/8);
	}
	else {
	  fHch_down->Fill(channel);
	  fHdet_down->Fill(channel/8);
	}
  	// std::cout<<"TSV: "<<GetName() // Please, don't delete these lines, tsv
	// <<" "<<channel           // print hits for efficiency estimation
	//  <<" "<<time             // by mw1_eff code
	//  <<" "<<fKey[i]
	//  <<'\n';
      }
    }
    float timeT0 = fMT_T0;
    std::vector<double> data;
    data.push_back(time - timeT0);
//     data.push_back(fKey[i]);
    data.push_back(key);
    fDigits.insert(make_pair(channel, CDigit1(channel,data)));
  }

  fHhits->Fill(fNhitsKept);

  if (thr_flag) TThread::UnLock();
}

void PlaneMuonWallA::ControlPanel(const TGWindow *p, const TGWindow *main) {

  if (!fControlPanel) fControlPanel = new PlanePanel(p, main, 100, 100, this);
}

void PlaneMuonWallA::Clusterize(void) {
  if(!fGeom) return;

//    set<CCluster1>& clusters = fGeom->Clusterize(fDigits);
//      if (thr_flag) TThread::Lock();
//    int nclusters=0;

//    typedef set<CCluster1>::iterator IC;
//    for (IC i = clusters.begin(); i!=clusters.end(); i++) {
//      float cpos = i->pos;

//        float time = i->dt[0];
//      //int time = int(CS::ChipF1::TimeDifference(i->dt[0],reftime));
//      if (fVct->Test(time) &&
//          fVcch->Test(cpos)) {
//        fHct->Fill(time);
//        fHcch->Fill(cpos);
//        nclusters++;
//      }
//      else {
//        clusters.erase(i);
//      }
//    }
//    fNCl = nclusters; // Update number of clusters after cut's been applied
//    fHchits->Fill(fNCl);
  //  if (thr_flag) TThread::UnLock();
}

void PlaneMuonWallA::TextOutput(ostream &) {

}



