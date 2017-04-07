#include "PlaneMwpc.h"
#include "Plane1VPanel.h"
#include "ChipF1.h"
#include "TriggerTime.h"
#include "Reference.h"
#include "TProfile.h"
#include "MwpcReconst.h"


ClassImp(PlaneMwpc);

//float f1_reso = 0.123;
float f1_reso = 0.12892;

extern CChamber* chamber[10];

struct MWPC_MAP_Struc
{
  char planeName[16];
  int  nStation;
  int  nPlane;
  int  low_time;
  int  hi_time;
};

static MWPC_MAP_Struc MwpcMap[34] =
{
  { "PS01X",  1, 0, -1550, -1425 },
  { "PS01Y",  1, 1, -1550, -1425 },
  { "PS01U",  1, 2, -1550, -1425 },
  { "PS01V",  1, 3, -1550, -1425 },
  { "PA01X",  2, 0, -1275, -1160 },
  { "PA01U",  2, 2, -1275, -1160 },
  { "PA01V",  2, 3, -1275, -1160 },
  { "PA02X",  3, 0, -1275, -1150 },
  { "PA02U",  3, 2, -1275, -1150 },
  { "PA02V",  3, 3, -1275, -1150 },
  { "PA03X",  4, 0, -1250, -1125 },
  { "PA03U",  4, 2, -1250, -1125 },
  { "PA03V",  4, 3, -1250, -1125 },
  { "PA04X",  5, 0, -1250, -1125 },
  { "PA04U",  5, 2, -1250, -1125 },
  { "PA04V",  5, 3, -1250, -1125 },
  { "PA05X",  6, 0, -1250, -1125 },
  { "PA05U",  6, 2, -1250, -1125 },
  { "PA05V",  6, 3, -1250, -1125 },
  { "PA06X",  7, 0, -1225, -1100 },
  { "PA06U",  7, 2, -1225, -1100 },
  { "PA06V",  7, 3, -1225, -1100 },
  { "PB01X",  8, 0, -1200, -1075 },
  { "PB01U",  8, 2, -1200, -1075 },
  { "PB02V",  8, 3, -1200, -1075 },
  { "PB03X",  9, 0, -1180, -1060 },
  { "PB03U",  9, 2, -1180, -1060 },
  { "PB04V",  9, 3, -1180, -1060 },
  { "PB05X", 10, 0, -1175, -1050 },
  { "PB05U", 10, 2, -1175, -1050 },
  { "PB06V", 10, 3, -1175, -1050 },
  { "PA11X", 11, 0, -1225, -1100 },
  { "PA11U", 11, 2, -1225, -1100 },
  { "PA11V", 11, 3, -1225, -1100 }
};


PlaneMwpc::PlaneMwpc(const char *detname,int nchan, int center, int width):
  Plane1V(detname,nchan,center,width) {
  INeedGeom();

  fNstation = -1;
  fNplane = -1;
  for(int i = 0; i < 34; i++)
    {
      if(strncmp(detname,MwpcMap[i].planeName,5) == 0)
	{
	  fNstation = MwpcMap[i].nStation;
	  fNplane   = MwpcMap[i].nPlane;
	  fLowTime  = MwpcMap[i].low_time;
	  fHiTime   = MwpcMap[i].hi_time;
	  fNwires   = (fNplane == 1)?520:760;
	  break;
	}
    }
}


void PlaneMwpc::Init(TTree *tree) {

  fFillfDigits = false;

  // Plane1V::Init(tree);
  ///**********************************************///
  fRateCounter = 0;

  // hits --------------------------------------------------------

  //multiplicity within the time cut on the whole detector
  std::string hitsname = fName + "_hits";
  fHhits=new TH1F_Ref(hitsname.c_str(),hitsname.c_str(),10,0,10,fRateCounter);
  AddHistogram(fHhits);
  std::string hitsleavlist = hitsname + "/I";

  //branch : hit map
  std::string chname = fName + fVch->GetName();
  //  fHch=new TH1F_Ref(chname.c_str(),chname.c_str(),
  fHch=new TH1F_Ref(chname.c_str(),chname.c_str(),
		   fVch->GetNbins(),
		   fVch->GetMin(),
		   fVch->GetMax(),
		   fRateCounter);
  AddHistogram(fHch);
  std::string chleavlist = chname + "[" + hitsname +"]/F";

  //branch : time
  std::string tname = fName + fVt->GetName();
  fHt=new TH1F_Ref(tname.c_str(),tname.c_str(),
		   fVt->GetNbins(),
		   fVt->GetMin(),
		   fVt->GetMax(),
		   fRateCounter);
  AddHistogram(fHt);
  std::string tleavlist = tname + "[" + hitsname +"]/F";

  // time vs channel histo
  std::string tcname = fName + "_tVSch";
  fHtvsch=new TH2F(tcname.c_str(),tcname.c_str(),
		   fVch->GetNbins(),
		   fVch->GetMin(),
		   fVch->GetMax(),
		   fVt->GetNbins(),
		   fVt->GetMin(),
		   fVt->GetMax());
  fHtvsch->SetOption("col");
  AddHistogram(fHtvsch);

  // On-Trigger Time Distribution
  std::string onttname = fName + "_t_on_trig";
  fHonTT=new TH1F(onttname.c_str(),onttname.c_str(),   //onTT is "off Trigger Time"
		  fOnTrigTVarID->GetNbins(),
		  fOnTrigTVarID->GetMin(),
		  fOnTrigTVarID->GetMax());
  AddHistogram(fHonTT);

  // Off-Trigger Time Distribution
  std::string offttname = fName + "_t_off_trig";
  fHoffTT=new TH1F(offttname.c_str(),offttname.c_str(),   //offTP is "off Trigger Time"
		   fOffTrigTVarID->GetNbins(),
		   fOffTrigTVarID->GetMin(),
		   fOffTrigTVarID->GetMax());
  AddHistogram(fHoffTT);


  // On-Trigger Profile
  std::string offtcname = fName + "_ch_on_trig";
  fHonTP=new TH1F(offtcname.c_str(),offtcname.c_str(),   //onTP is "on Trigger Profile"
		  fVch->GetNbins(),
		  fVch->GetMin(),
		  fVch->GetMax());
  AddHistogram(fHonTP);

  // Off-Trigger Profile
  std::string offtpname = fName + "_ch_off_trig";
  fHoffTP=new TH1F(offtpname.c_str(),offtpname.c_str(),   //offTP is "off Trigger Profile"
		   fVch->GetNbins(),
		   fVch->GetMin(),
		   fVch->GetMax());
  AddHistogram(fHoffTP);

  // Corrected On-Trigger Profile
  std::string ontcorrname = fName + "_ch_on_corr_trig";
  fHoncorrTP=new TH1F(ontcorrname.c_str(),ontcorrname.c_str(),
		      fVch->GetNbins(),
		      fVch->GetMin(),
		      fVch->GetMax());
  AddHistogram(fHoncorrTP);

  // Rates HISTOGRAMS
  std::string rname = fName + "_rates";
  fHrates=new TH1F_Ref(rname.c_str(),rname.c_str(),
		   fVch->GetNbins(),
		   fVch->GetMin(),
		   fVch->GetMax(), fRateCounter, true);
  fHrates->SetYTitle("rates per channel (kHz)");
  fHrates->SetTitleOffset(1.2,"Y");
  ((TH1F_Ref*)fHrates)->SetReference(fReferenceDirectory);
  AddHistogram(fHrates);

  if(tree) {
    fIsInTree = true;
    tree->Branch(hitsname.c_str(),&fNhitsKept,
		 hitsleavlist.c_str(),32000);
    tree->Branch(chname.c_str(),fVch->GetValues(),
		 chleavlist.c_str(),32000);
    tree->Branch(tname.c_str(),fVt->GetValues(),
		 tleavlist.c_str(),32000);
  }
  ///**********************************************///
  // if geometry has been set, it shall be used...
  float posmax = 50;
  if(fGeom)
    posmax = fGeom->GetHalfSize();

  fVcch = AddVariable("_cch",100,-posmax,posmax,fNchan*fMAX_MULT);
  fVct = AddVariable("_ct",100,fCenter-fWidth,fCenter+fWidth,fNchan*fMAX_MULT);
  fVcs = AddVariable("_cs",10,0,10,fNchan*fMAX_MULT);

  // clusters ----------------------------------------------------------

  //multiplicity within the time cut on the whole detector
  std::string chitsname = fName + "_chits";
  fHchits=new TH1F(chitsname.c_str(),chitsname.c_str(),10,0,10);
  AddHistogram(fHchits);
  std::string chitsleavlist = chitsname + "/I";

  //branch : cluster map
  std::string cchname = fName + fVcch->GetName();
  fHcch=new TH1F(cchname.c_str(),cchname.c_str(),
		 fVcch->GetNbins(),
		 fVcch->GetMin(),
		 fVcch->GetMax());
  AddHistogram(fHcch);
  std::string cchleavlist = cchname + "[" + chitsname +"]/F";

  //branch : cluster time
  std::string ctname = fName + fVct->GetName();
  fHct=new TH1F(ctname.c_str(),ctname.c_str(),
		fVct->GetNbins(),
		fVct->GetMin(),
		fVct->GetMax());
  AddHistogram(fHct);
  std::string ctleavlist = ctname + "[" + chitsname +"]/F";

  //Cluster size
  std::string csname = fName+"_cs";
  fHcs=new TH1F(csname.c_str(),csname.c_str(),
		fVcs->GetNbins(),
		fVcs->GetMin(),
		fVcs->GetMax());
  AddHistogram(fHcs);
  std::string csleavlist = csname + "[" + chitsname +"]/F";


  if(tree) {
    tree->Branch(chitsname.c_str(),&fNclustKept,
		 chitsleavlist.c_str(),32000);
    tree->Branch(cchname.c_str(),fVcch->GetValues(),
		 cchleavlist.c_str(),32000);
    tree->Branch(ctname.c_str(),fVct->GetValues(),
		 ctleavlist.c_str(),32000);
    tree->Branch(csname.c_str(),fVcs->GetValues(),
		 csleavlist.c_str(),32000);
  }
  ((TH1F_Ref*)fHhits)->SetReference(fReferenceDirectory);
  ((TH1F_Ref*)fHch)->SetReference(fReferenceDirectory);
  ((TH1F_Ref*)fHt)->SetReference(fReferenceDirectory);

  if(fNstation != -1)
    {
      std::string histname = fName + "_eff";
      fPeff = new TProfile(histname.c_str(),histname.c_str(),
			   fNwires,0,fNwires,0,100);
      AddHistogram(fPeff);
      fPeff->GetXaxis()->SetTitle("Wire number");

      histname = fName + "_noise";
      fHns = new TH1F(histname.c_str(),histname.c_str(),fNwires, 0, fNwires);
      AddHistogram(fHns);

      fHns->GetXaxis()->SetTitle("Wire number");
    }
}


void PlaneMwpc::ResetHistograms() {
  Plane::ResetHistograms();
}


void PlaneMwpc::Reset() {
  Plane1V::Reset();
  fNclustKept=0;
}

void PlaneMwpc::StoreDigit(CS::Chip::Digit* digit) {
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


#ifndef __CINT__
void PlaneMwpc::EndEvent(const CS::DaqEvent &event)
{
//   const float reftime = event.GetTT().GetTimeNorm();
  fDigits.clear();
  //int max_wire = -1;
  std::vector< std::pair<int,int> > hits;
  //int clust_mult = 0;
  //int cluster_size = 0;
  //int last_hit = -2;
  //int last_time = 0;
  if (thr_flag) TThread::Lock();
  //nevent++;
#ifdef MWPC_DEBUG
  std::cout << std::endl << "Plane " << GetName() << " hits:" << std::endl;
  std::cout << "PLane " << GetName() << ": time range set to "
        << fVt->GetMin() << " -> " << fVt->GetMax()
	<< "; t0 = " << fVt->GetMin()+20 << std::endl;
#endif
//   for (int i=0; i<fNhits; i++) {
//     //       int time=fTime[i]-reftime;
//     int time = (int)CS::ChipF1::TimeDifference(fTime[i],reftime);
//     int channel=fChannel[i];
//     //       if(time>fUP_BOUND) time-=fMAX_CLOCK;
//     //       if (time<fLOW_BOUND) time+=fMAX_CLOCK;
// 
//     time = (int)(f1_reso*time);  // Time in ns
  typedef std::list<CS::Chip::Digit*>::iterator lDIter;
  for (lDIter ii = lDigits.begin(); ii != lDigits.end(); ii++) {
    CS::ChipF1::Digit* iii = dynamic_cast<CS::ChipF1::Digit*> (*ii);
    if (!iii) {
      std::cerr<<"PlaneMwpc::EndEvent: a digit is not a F1 one, strange...\n";
      continue;
    }
    register int time = (int) iii->GetTimeDecoded();
    register int channel = iii->GetChannel();

    if( fVch->Test(channel) ) {

      if( fVt->Test(time) ) {
	
	std::vector<double> data;
	data.push_back(time);
	fDigits.insert(make_pair(channel,CDigit1(channel,data)));

	fVch->Store(channel);
	fVt->Store(time);
	fHch->Fill(channel);
	fHt->Fill(time);
	fHtvsch->Fill(channel,time);
	fNhitsKept++;
      }

      if( fOnTrigTVarID->Test(time) )  {
	fOnTrigTVarID->Store(time);
	fHonTT->Fill(time);        // on trigger time distribution
	fHonTP->Fill(channel);  // on trigger profile
	fNhitsKeptOnTrig++;
      }

      if( fOffTrigTVarID->Test(time) )  {
	fOffTrigTVarID->Store(time);
	fHoffTT->Fill(time);        // off trigger time distribution
	fHoffTP->Fill(channel);  // off trigger profile
	fNhitsKeptOffTrig++;
      }
	
    }
    if(fNstation != -1 && chamber[fNstation-1])
      {
	if(time > fLowTime && time < fHiTime)
	  {
	    chamber[fNstation-1]->PutHit(fNplane,channel,time);
	  }
      CPlane* plane = chamber[fNstation-1]->planes[fNplane];
      int nRecTrig = plane->GetNRecEvents();
      if(fRateCounter%1000 == 0 && nRecTrig)
	{
	  float tot_noise = 0.;
	  for(int i = 0; i < fNwires; i++)
	    {
	      CWire& w = plane->wires[i];
	      if(w.GetNTrig())
		{
		  float nHits = w.GetNEff()+w.GetNIneff();
		  if(nHits > .1)
		    {
		      float efficiency = w.GetNEff()/nHits*100.;
		      fPeff->Fill(i,efficiency,w.GetNTrig());
		    }
		  float noise = w.GetNNoise();
		  tot_noise += noise;
		  fHns->Fill(i,noise);
		}
	    }
	  plane->ClearStatistic();
	}
      }
  }
  fHhits->Fill(fNhitsKept);


  if((fRateCounter%Plane1V::fRATE_UPDATE)==0) {

    // corrected on-trigger profile
    if( fOffTrigTVarID->GetMax() != fOffTrigTVarID->GetMin() ) {
      float scale = (float) ( fOnTrigTVarID->GetMax() - fOnTrigTVarID->GetMin() )/( fOffTrigTVarID->GetMax() - fOffTrigTVarID->GetMin() );
      for(register int i=1; i<=fHonTP->GetNbinsX(); i++)
	fHoncorrTP->SetBinContent(i, (float) fHonTP->GetBinContent(i) - scale*fHoffTP->GetBinContent(i) );
    }

    // rates histogram
    // A. Ferrero - 20/05/2004: removed fF1_TICK from rate computation 
    //float timewin = ( fOffTrigTVarID->GetMax() - fOffTrigTVarID->GetMin() )*fF1_TICK*fRateCounter;
    float timewin = ( fOffTrigTVarID->GetMax() - fOffTrigTVarID->GetMin() )*fRateCounter*1e-6;  // times already in ns
    if( timewin > 0 )
      for(register int i=1; i<=fHonTP->GetNbinsX(); i++) {
	fHrates->SetBinContent(i,(float) fHoffTP->GetBinContent(i)/timewin);
	fHrates->SetBinError(i,(float) fHoffTP->GetBinError(i)/timewin);
      }
  }

  if (thr_flag) TThread::UnLock();

}
#endif

void PlaneMwpc::Clusterize() {

  if(!fGeom) return;
  std::set<CCluster1>& clusters = fGeom->Clusterize(fDigits);

  if (thr_flag) TThread::Lock();
  int nclusters=0;
  typedef std::set<CCluster1>::iterator IC;
  for (IC i = clusters.begin(); i!=clusters.end(); i++) {

    float cpos = i->pos;
    float size = i->size;
    float time = i->dt[0];


    if (fVct->Test(time) &&
      	fVcch->Test(cpos) &&
      	fVcs->Test(size)) {


      fHcch->Fill(cpos);
      fHct->Fill(time);
      fHcs->Fill(size);

      fVcch->Store(cpos);
      fVct->Store(time);

      nclusters++;
    }
    else {
      //update set of clusters
//      clusters.erase(i);
    }
  }
  fNclustKept = nclusters; // Update number of clusters after cut's been applied
  fHchits->Fill(nclusters);

  if (thr_flag) TThread::UnLock();
}







