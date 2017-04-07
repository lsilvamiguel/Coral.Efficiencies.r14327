#include <cmath>
#include "PlaneMumega.h"
#include "PlaneTrigHodo.h"
#include "ChipF1.h"
#include "TriggerTime.h"

#include <string>

#include "TProfile.h"
#include "TF1.h"

ClassImp(PlaneMumega);

const int PlaneMumega::fMAX_MULT = 8;
const float PlaneMumega::fF1_TICK = 128.e-9; // in ms
//const int PlaneMumega::fF1_GATE = 3500; // in ms - setup file ~saclay/f1/narrow.f1
//const int PlaneMumega::fF1_GATE = 20000; // in ms - setup file ~saclay/f1/wide.f1
const int PlaneMumega::fRATE_UPDATE = 100;
const float PlaneMumega::fL_WGHT = .6;
const float PlaneMumega::fT_WGHT = .4;
const float PlaneMumega::fHIT_TMIN = -10000-500;
const float PlaneMumega::fHIT_TMAX = -10000+500;
const float PlaneMumega::fMT_T0 = (fHIT_TMIN+fHIT_TMAX)/2;
const float PlaneMumega::fCL_TMIN = -10000-250;
const float PlaneMumega::fCL_TMAX = -10000+300;

PlaneMumega::PlaneMumega(const char *detname,int nchan, int center, int width, int parity)
  : Plane(detname), fNchan(nchan), fNCl(0),
    fPARITY_LEAD(parity) {

  INeedGeom();


  fVch = AddVariable("_ch",fNchan,0,fNchan,fNchan*fMAX_MULT);
  fVt = AddVariable("_t",100,center-width,center+width,fNchan*fMAX_MULT);

  fVhch = AddVariable("_hch",fNchan,0,fNchan,fNchan*fMAX_MULT);
  fVht = AddVariable("_ht",100,center-width,center+width,fNchan*fMAX_MULT);
  fVhlt = AddVariable("_hlt",100,-99999,99999,fNchan*fMAX_MULT);
  fVhtt = AddVariable("_htt",100,-99999,99999,fNchan*fMAX_MULT);
  fVhtot = AddVariable("_htot",100,center-width,center+width,fNchan*fMAX_MULT);
  fVhtcalib = AddVariable("_htcalib",100,-width,width,fNchan*fMAX_MULT);
}


void PlaneMumega::Init(TTree* tree) {

  fRateCounter = 0;


  // leading edge histograms -----------------------------------

  int lcol=7; // leading edge color

  // lead edge multiplicity
  std::string name = fName+"_lhits";
  std::string title = fName+" Lead Edge Mult";
  fHlhits = new TH1F(name.c_str(),title.c_str(),100,0,100);
  fHlhits->SetFillColor(lcol);
  AddHistogram(fHlhits);

  // lead edge profile
  name = fName+"_lch";
  title = fName+" Lead Edge Profile";
  fHlch=new TH1F(name.c_str(),title.c_str(),fNchan,0,fNchan);
  fHlch->SetFillColor(lcol);
  AddHistogram(fHlch);

  // lead edge time
  name = fName+"_lt"; title = fName+" Lead Edge Time";
  fHlt=new TH1F(name.c_str(),title.c_str(),
	     fVt->GetNbins(),
	     fVt->GetMin(),
	     fVt->GetMax());
  fHlt->SetFillColor(lcol);
  AddHistogram(fHlt);

  // ********** TRAILING HISTOGRAMS **********

  // trail edge time
  name = fName+"_tt"; title = fName+" Trail Edge Time";
  fHtt=new TH1F(name.c_str(),title.c_str(),
	     fVt->GetNbins(),
	     fVt->GetMin(),
	     fVt->GetMax());
  fHtt->SetFillColor(lcol);
  AddHistogram(fHtt);


  // hit histograms ---------------------------------------------------

  int hitcol=3; // color for hit histograms

  // Hit multiplicity
  std::string hitsname = fName + "_hhits";
  title = fName + " Hit Multiplicities";
  fHhhits=new TH1F_Ref(hitsname.c_str(),title.c_str(),50,0,50,fRateCounter);
  fHhhits->SetFillColor(hitcol);
  ((TH1F_Ref*)fHhhits)->SetReference(fReferenceDirectory);
  AddHistogram(fHhhits);
  std::string hitsleavlist = hitsname + "/I";

  // Hit profile
  std::string chname = fName + "_hch";
  title = fName + " Hit Profile";
  fHhch=new TH1F_Ref(chname.c_str(),title.c_str(),
	     fVhch->GetNbins(),
	     fVhch->GetMin(),
	     fVhch->GetMax(),fRateCounter);
  fHhch->SetFillColor(hitcol);
  ((TH1F_Ref*)fHhch)->SetReference(fReferenceDirectory);
  AddHistogram(fHhch);
  std::string chleavlist = chname + "[" + hitsname +"]/F";

  // Hit time
  std::string d1name = fName + "_ht"; title = fName + " Hit Time";
  fHht=new TH1F_Ref(d1name.c_str(),title.c_str(),
	            fVht->GetNbins(),
	            fVht->GetMin(),
	            fVht->GetMax(),fRateCounter);
  fHht->SetFillColor(hitcol);
  ((TH1F_Ref*)fHht)->SetReference(fReferenceDirectory);
  AddHistogram(fHht);
  std::string d1leavlist = d1name + "[" + hitsname +"]/F";

  // Hit time vs ch
  std::string htvschname = fName + "_htvsch"; title = fName + " Hit Time VS Chan";
  fHhtvsch=new TH2F(htvschname.c_str(),title.c_str(),
		    fVhch->GetNbins(),
		    fVhch->GetMin(),
		    fVhch->GetMax(),
		    fVht->GetNbins(),
		    fVht->GetMin(),
		    fVht->GetMax());
  fHhtvsch->SetOption("colz");
  AddHistogram(fHhtvsch);

  // calibrated hit time vs ch
  std::string htcalibvschname = fName + "_htcalibvsch"; title = fName + " Hit Time VS Chan, calibrated";
  fHhtcalibvsch = new TH2F(htcalibvschname.c_str(),title.c_str(),
			   fVhch->GetNbins(),
			   fVhch->GetMin(),
			   fVhch->GetMax(),
			   fVhtcalib->GetNbins(),
			   fVhtcalib->GetMin(),
			   fVhtcalib->GetMax() );
  fHhtcalibvsch->SetOption("colz");
  AddHistogram(fHhtcalibvsch);


  // Hit ToT
  std::string d2name = fName + "_htot";
  title = fName + " Hit ToT";
  fHhtot=new TH1F_Ref(d2name.c_str(),title.c_str(),
	              fVhtot->GetNbins(),
	              fVhtot->GetMin(),
	              fVhtot->GetMax(),fRateCounter);
  fHhtot->SetFillColor(hitcol);
  ((TH1F_Ref*)fHhtot)->SetReference(fReferenceDirectory);
  AddHistogram(fHhtot);
  std::string d2leavlist = d2name + "[" + hitsname +"]/F";

  // Hit rates
  std::string rname = fName + "_hrates";
  title = fName + " Hit Rates";
  fHhrates=new TH1F_Ref(rname.c_str(),title.c_str(),
		        fVch->GetNbins(),
		        fVch->GetMin(),
		        fVch->GetMax(),fRateCounter,true);
  fHhrates->SetYTitle("rates per channel (kHz)");
  fHhrates->SetTitleOffset(1.2,"Y");
  fHhrates->SetFillColor(hitcol);
  fHhrates->SetOption("hist");
  ((TH1F_Ref*)fHhrates)->SetReference(fReferenceDirectory);
  AddHistogram(fHhrates);

  std::string ltname = fName + "_hlt";
  std::string ltleavlist = fName + "[" + hitsname + "]";
  std::string ttname = fName + "_htt";
  std::string ttleavlist = fName + "[" + hitsname + "]";

  // Hit time with calibration
  std::string tcalibname = fName + "_htcalib"; title = fName + " Hit Time with Calibration";
  fHhtcalib=new TH1F_Ref(tcalibname.c_str(),title.c_str(),
	                 fVhtcalib->GetNbins(),
	                 fVhtcalib->GetMin(),
	                 fVhtcalib->GetMax(),fRateCounter);
  fHhtcalib->SetFillColor(hitcol);
  ((TH1F_Ref*)fHhtcalib)->SetReference(fReferenceDirectory);
  AddHistogram(fHhtcalib);

  // cluster histos -----------------------------------------------

  int ccol=8; // color for cluster histos

  // if geometry has been set, it shall be used...
  float posmax = 20;
  if(fGeom)
    posmax = fGeom->GetHalfSize();

  fVcch = AddVariable("_cch",fNchan,-posmax,posmax,fNchan);
  fVct = AddVariable("_ct",100,-500,500,fNchan);
  fVctot = AddVariable("_ctot",100,0,5000,fNchan); // *=.13=650 ns
  fVcs = AddVariable("_cs",100,0,fNchan,fNchan);
  fVcamp = AddVariable("_camp",100,0,1000,fNchan);
  fVcres = AddVariable("_cres",100,0,1,fNchan);

  //Cluster multiplicity
  std::string chitsname = fName+"_chits";
  title = fName+" Cluster Multiplicities";
  fHchits=new TH1F_Ref(chitsname.c_str(),title.c_str(),51,-0.5,50.5,fRateCounter);
  fHchits->SetFillColor(ccol);
  ((TH1F_Ref*)fHchits)->SetReference(fReferenceDirectory);
  AddHistogram(fHchits);
  std::string chitsleavlist = chitsname + "/I";

  //Cluster profile
  std::string cchname = fName+"_cch";
  title = fName+" Cluster Profile";
  //fHcch=new TH1F(cchname.c_str(),title.c_str(),2*fNchan-1,0,fNchan);
  fHcch=new TH1F_Ref(cchname.c_str(),title.c_str(),100,-posmax,posmax,fRateCounter);
  fHcch->SetFillColor(ccol);
  ((TH1F_Ref*)fHcch)->SetReference(fReferenceDirectory);
  AddHistogram(fHcch);
  std::string cchleavlist = cchname + "[" + chitsname +"]/F";

  //Cluster time
  std::string cd1name = fName+"_ct"; title = fName+" Cluster Time";
  fHct=new TH1F_Ref(cd1name.c_str(),title.c_str(),100,-500,500,fRateCounter);
  fHct->SetFillColor(ccol);
  ((TH1F_Ref*)fHct)->SetReference(fReferenceDirectory);
  AddHistogram(fHct);
  std::string cd1leavlist = cd1name + "[" + chitsname +"]/F";

  //Cluster time vs  channel
  std::string ctvschname = fName+"_ctVSch"; title = fName+" Cluster Time vs Channel";
  fHctvsch= new TH2F(ctvschname.c_str(),title.c_str(),
		     50,-posmax,posmax,
		     100,-500,500);
  fHctvsch->SetOption("colz");
  AddHistogram(fHctvsch);

  //Cluster ToT vs time
  std::string ctotvstname = fName+"_ctotVSt"; title = fName+" Cluster ToT VS Time";
  fHctotvst= new TH2F(ctotvstname.c_str(),title.c_str(),
		      100,-2000,2000,
		      100,0,4000);
  fHctotvst->SetOption("colz");
  AddHistogram(fHctotvst);

  //Cluster ToT
  std::string cd2name = fName+"_ctot"; title = fName+" Cluster ToT";
  fHctot=new TH1F_Ref(cd2name.c_str(),title.c_str(),100,0,3000,fRateCounter);
  fHctot->SetFillColor(ccol);
  ((TH1F_Ref*)fHctot)->SetReference(fReferenceDirectory);
  AddHistogram(fHctot);
  std::string cd2leavlist = cd2name + "[" + chitsname +"]/F";

  //Cluster size
  std::string cd3name = fName+"_cs"; title = fName+" Cluster Size";
  fHcs=new TH1F_Ref(cd3name.c_str(),title.c_str(),32,-.5,31.5,fRateCounter);
  fHcs->SetFillColor(ccol);
  ((TH1F_Ref*)fHcs)->SetReference(fReferenceDirectory);
  AddHistogram(fHcs);
  std::string cd3leavlist = cd3name + "[" + chitsname +"]/F";

  //Cluster amplitude
  std::string campname = fName+"_camp"; title = fName+" Cluster Amplitude";
  fHcamp=new TH1F_Ref(campname.c_str(),title.c_str(),100,0,500,fRateCounter);
  fHcamp->SetXTitle("Q * 1000 e-");
  fHcamp->SetFillColor(ccol);
  ((TH1F_Ref*)fHcamp)->SetReference(fReferenceDirectory);
  AddHistogram(fHcamp);
  std::string campleavlist = campname + "[" + chitsname +"]/F";

  fHcampsig=new TH1F_Ref((campname+"sig").c_str(),(title + " cs>1").c_str(),100,0,500,fRateCounter);
  fHcampsig->SetXTitle("Q * 1000 e-");
  fHcampsig->SetFillColor(ccol);
  ((TH1F_Ref*)fHcampsig)->SetReference(fReferenceDirectory);
  AddHistogram(fHcampsig);


  //Cluster resolution
  std::string cresname = fName+"_cres"; title = fName+" Cluster Resolution";
  fHcres=new TH1F(cresname.c_str(),title.c_str(),100,0,0.1);
  fHcres->SetXTitle("Q * 1000 e-");
  fHcres->SetFillColor(ccol);
  AddHistogram(fHcres);
  std::string cresleavlist = cresname + "[" + chitsname +"]/F";

  if (fExpertHistos) {
    // Trigger type vs hit profile
    std::string htrvschname = fName+"_htrVSch"; title = fName+" Trigger type vs Hit profile";
    fHhtrigvsch = new TH2F(htrvschname.c_str(),title.c_str(),
		           fVhch->GetNbins(),
		           fVhch->GetMin(),
		           fVhch->GetMax(),
		           fMaxTriggerNumber, 0., fMaxTriggerNumber*1.);
    fHhtrigvsch->SetOption("colz");
    AddHistogram(fHhtrigvsch);

    // Trigger type vs hit rate
    char istring[12];
    for (register int i=0; i < fMaxTriggerNumber; i++) {
      sprintf(istring,"%d",i+1);
      std::string htrvsratesname = fName+"_hrates_trig_"+istring;
      title = fName+" Hit rates for trigger type "+istring;
      fHhtrigvsrates[i] = new TH1F(htrvsratesname.c_str(),title.c_str(),
		                   fVhch->GetNbins(),
		                   fVhch->GetMin(),
		                   fVhch->GetMax());
      fHhtrigvsrates[i]->SetYTitle("rates per channel (kHz)");
      fHhtrigvsrates[i]->SetOption("hist");
      fHhtrigvsrates[i]->SetFillColor(hitcol);
      AddHistogram(fHhtrigvsrates[i]);

      fHhNtrig[i] = 0;
      fNtrig[i] = 0;
    }

  }


  if(tree) {
    fIsInTree = true;
#define uM_CLUSTERS_TREE
#ifndef uM_CLUSTERS_TREE
    tree->Branch(hitsname.c_str(),&fNhitsKept,
		 hitsleavlist.c_str(),32000);
    tree->Branch(chname.c_str(),fVhch->GetValues(),
		 chleavlist.c_str(),32000);
    tree->Branch(d1name.c_str(),fVht->GetValues(),
		 d1leavlist.c_str(),32000);
    tree->Branch(ltname.c_str(),fVhlt->GetValues(),
		 ltleavlist.c_str(),32000);
    tree->Branch(ttname.c_str(),fVhtt->GetValues(),
		 ttleavlist.c_str(),32000);
    tree->Branch(d2name.c_str(),fVhtot->GetValues(),
		 d2leavlist.c_str(),32000);
#else
    tree->Branch(chitsname.c_str(),&fNCl,
		 chitsleavlist.c_str(),32000);
    tree->Branch(cchname.c_str(),fVcch->GetValues(),
		 cchleavlist.c_str(),32000);
    tree->Branch(cd1name.c_str(),fVct->GetValues(),
		 cd1leavlist.c_str(),32000);
    tree->Branch(cd2name.c_str(),fVctot->GetValues(),
		 cd2leavlist.c_str(),32000);
    tree->Branch(cd3name.c_str(),fVcs->GetValues(),
		 cd3leavlist.c_str(),32000);
    tree->Branch(campname.c_str(),fVcamp->GetValues(),
		 campleavlist.c_str(),32000);
    tree->Branch(cresname.c_str(),fVcres->GetValues(),
		 cresleavlist.c_str(),32000);
#endif
  }
}


PlaneMumega::~PlaneMumega() {
}


void PlaneMumega::StoreDigit(int channel, int data) {
  std::cerr<<"PlaneMumega::StoreDigit(int, int): Error ! this method is obsolete !!!\n";
//   if (fNhits < fNchan*fMAX_MULT) {
//     fChannel[fNhits]=channel;
//     fData[fNhits]=data;
//     fNhits++;
//   }
}


void PlaneMumega::StoreDigit(CS::Chip::Digit* digit) {
//   std::vector<float> data=digit->GetNtupleData();
//   this->StoreDigit((int)data[0],(int)data[1]);
  lDigits.push_back(digit);
}


void PlaneMumega::ResetHistograms() {
  Plane::ResetHistograms();
  if (fExpertHistos) for (register int i = 0; i < fMaxTriggerNumber ; i++) {
    fHhNtrig[i] = 0;
    fNtrig[i] = 0;
  }
}


void PlaneMumega::EndEvent(const CS::DaqEvent &event) {

//   const float reftime = event.GetTT().GetTimeNorm();

// increase counts of trigger types if needed
  register unsigned int triggerwd = event.GetTrigger() ;
  if (fExpertHistos) {
    register int m=1;
    for (register int i=0; i<fMaxTriggerNumber; i++) {
      if ((triggerwd & m) >> i) {
        fNtrig[i]++;
      }
      m <<= 1;    
    }
  }

  // Reset the divisions of the X axis. Otherwise labels are overlapping.
  // These same instructions, when executed @ booking time, have no effect (!?)
  if (thr_flag) TThread::Lock();
  fHht->GetXaxis()->SetNdivisions(505);
  fHct->GetXaxis()->SetNdivisions(505);

  int nedges, lastchan, prvchan; static int lasttime;
  fDigits.clear();
  lastchan=prvchan = -1;
  nedges = 0;
//   for (i = 0; i<fNhits; i++) {
  typedef std::list<CS::Chip::Digit*>::iterator lDIter;
  for (lDIter ii = lDigits.begin(); ii != lDigits.end(); ii++) {
    CS::ChipF1::Digit* iii = dynamic_cast<CS::ChipF1::Digit*> (*ii);
    if (!iii) {
      std::cerr<<"PlaneMumega::EndEvent: a digit is not a F1 one, strange...\n";
      continue;
    }

//     int time = static_cast<int>(CS::ChipF1::TimeDifference(fData[i],reftime));
//     int chan = fChannel[i];
    int time = (int) (iii->GetTimeDecoded() / iii->GetTimeUnit());
    int chan = iii->GetChannel();


#define REMOVE_NOISY
#ifdef REMOVE_NOISY
    if(fUseCalib && calib_data[chan].on == 1) continue;
#endif

    if(fVch->Test(chan) && fVt->Test(time)) {

      if (IsLeading(iii->GetTime())) {
	
	fHlch->Fill(chan);
	fHlt->Fill(time);
	
	lasttime = time;
	lastchan = chan; // Open search for trailing
	
	nedges++;
      }
      else {
	
	fHtt->Fill(time);
	
	if (chan==lastchan) { // ***** HIT = LEADING + TRAILING
	
//  	  if ((!fNoisyCh.empty()) && (fNoisyCh.find(chan) != fNoisyCh.end()))
//  	    continue;

	  lastchan = -1; // Close search for trailing	
	  float mean = lasttime*fL_WGHT+time*fT_WGHT;  //MOD
	  int tot=time-lasttime;

	  if(fVht->Test(mean) &&
	     fVhlt->Test(lasttime) &&
	     fVhtt->Test(time) &&
	     fVhch->Test(chan) &&
	     fVhtot->Test(tot)) {

	    fHhch->Fill(chan);
	    fVhch->Store(chan);
	    if (fExpertHistos) {
              register int m=1;
//              for(register size_t i=0; i<sizeof(int)*8; i++) { // jardinage
              for (register int i=0; i<fMaxTriggerNumber; i++) {
                if ((triggerwd & m) >> i) {
 	          fHhtrigvsch->Fill(chan, i);
                  fHhNtrig[i]++;
                }
                m <<= 1;    
              }
            }

	    fVhlt->Store(lasttime);
	    fVhtt->Store(time);

	    fHht->Fill(mean);
	    fVht->Store(mean);
	
	    fHhtvsch->Fill(chan,mean);

	    fHhtot->Fill(tot);
	    fVhtot->Store(tot);

            float timeT0 = fMT_T0;
            if (fUseCalib) {
              timeT0 = calib_data[chan].t0;
	      // std::cout<<fName<<" "<<chan<<" "<<timeT0<<std::endl;

	      float calt = mean - timeT0;
	      fHhtcalib->Fill(calt);
	      fHhtcalibvsch->Fill(chan, calt);
	      fVhtcalib->Store(calt);
            }


	    fNhitsKept++;

//              if ((!fNoisyCh.empty()) && (fNoisyCh.find(chan) != fNoisyCh.end()))
//                continue;

	    // ***** SELECT HITS WITHIN TIME WINDOW - RETAIN THAT CLOSEST TO T0
	    if (chan!=prvchan) {
	      // NEW WIRE => SET PATTERN for CLUSTER SEARCH and ToT and Time
	      prvchan = chan;
	
	      std::vector<double> data;
	      data.push_back(mean - timeT0);
	      data.push_back(tot);
	      fDigits.insert(make_pair(chan,CDigit1(chan,data)));	
	    }
	    else {
	      CDigit1& digit=fDigits[chan];
	      if (fabs(mean-timeT0)<fabs(digit.dt[0])) {
		//std::cout<<"rewriting !"<<std::endl;
		// RETAIN HIT CLOSEST TO T0 => UPDATE tot and Time
		digit.ch = chan;
		digit.dt[0] = mean - timeT0;
		digit.dt[1] = tot;
	      }
	    }
	  }
	}
      }
    }
  }

  // Fill multiplicities
  fHlhits->Fill(nedges);
  fHhhits->Fill(fNhitsKept);

  if((fRateCounter%PlaneMumega::fRATE_UPDATE)==0) {
    float timewin = (fVht->GetMax()-fVht->GetMin())*fF1_TICK*fRateCounter;
    float timewintrig = (fVht->GetMax()-fVht->GetMin())*fF1_TICK;
    if(timewin) {
      for(register int i=1; i<=fHhch->GetNbinsX(); i++) {
	fHhrates->SetBinContent(i,fHhch->GetBinContent(i)/timewin);
	fHhrates->SetBinError(i,fHhch->GetBinError(i)/timewin);
	if (fExpertHistos) {
          for(register int j=1; j<=fHhtrigvsch->GetNbinsY(); j++) {
	    if (fNtrig[j-1] > 0) {
              fHhtrigvsrates[j-1]->SetBinContent(i,
                  fHhtrigvsch->GetBinContent(i,j) / timewintrig / fNtrig[j-1]);
	      fHhtrigvsrates[j-1]->SetBinError(i,
                  fHhtrigvsch->GetBinError(i,j) / timewintrig / fNtrig[j-1]);
            }
          }
        }
      }
      fHhrates->SetEntries(fHhch->GetEntries());
      if (fExpertHistos) for (register int i=0; i < fMaxTriggerNumber; i++) {
          fHhtrigvsrates[i]->SetEntries(fHhNtrig[i]);
      }
    }
  }

  if (thr_flag) TThread::UnLock();
}


void PlaneMumega::Clusterize() {

  //  static int n= 0;

  if(!fGeom) return;

  std::set<CCluster1>& clusters = fGeom->Clusterize(fDigits);

  int nclusters=0;
  if (thr_flag) TThread::Lock();
  typedef std::set<CCluster1>::iterator IC;
  for (IC i = clusters.begin(); i!=clusters.end(); i++) {

    float cpos = i->pos;
    float size = i->size;
    float time = i->dt[0];
    float tot  = i->dt[1];
    float amp  = i->dt[2];
    float res  = i->res;

    if (size < 0 || size > 1024) continue;  // bad cluster...

    if (fVct->Test(time) &&
	fVcch->Test(cpos) &&
	fVcs->Test(size) &&
	fVctot->Test(tot) &&
	fVcamp->Test(amp) &&
	fVcres->Test(res)) {

      TH2F *h2 = dynamic_cast<TH2F*> (fHctvsch);
      if(h2) h2->Fill(cpos,time);
      fHct->Fill(time);
      fHcch->Fill(cpos);
      fHctot->Fill(tot);
      TH2F *h22 = dynamic_cast<TH2F*> (fHctotvst);
      if(h22) h22->Fill(time,tot);
      fHcs->Fill(size);
      if(size>1) // remove noise
	fHcampsig->Fill(amp);
      fHcamp->Fill(amp);
      fHcres->Fill(res);

	// to be used in GroupMumega
      fVcch->Store(cpos);
      fVct->Store(time);
      fVctot->Store(tot);
      fVcs->Store(size);
      fVcamp->Store(amp);
      fVcres->Store(res);

      nclusters++;
    }
    else {
      //update set of clusters
      // clusters.erase(i);
    }
  }
  fNCl = nclusters; // Update number of clusters after cut's been applied
  fHchits->Fill(fNCl);

  if (thr_flag) TThread::UnLock();
}


void PlaneMumega::ControlPanel(const TGWindow *p, const TGWindow *main) {

  if (!fControlPanel) fControlPanel = new PlanePanel(p, main, 100, 100, this);
}


#if USE_DATABASE == 1
void PlaneMumega::ReadCalib(const tm &t)
{
  // read-in corresponding calibration constants
  try{
    ReadFromDataBase(calib_data,t);
//     std::cout<<"PlaneMumega::ReadCalib() ==> "<<this->GetName()<<" calibrations are found !"<<std::endl;

    if(calib_data.size() != (unsigned) fNchan) {
      std::cerr<<GetName()<<": Size of Calibration File is not correct ! Should be : "
	  <<fNchan<<" Is "<<calib_data.size()<<" "
	  <<t.tm_mday<<"."<<t.tm_mon+1<<"."<<t.tm_year+1900<<" "
	  <<t.tm_hour<<":"<<t.tm_min<<":"<<t.tm_sec<<std::endl;
//       for (std::vector<ChannelCalib>::iterator ii = calib_data.begin(); ii != calib_data.end(); ii++) {
//         std::cerr<<"ch "<<(*ii).ch<<" t0 "<<(*ii).t0<<" on "<<(*ii).on<<std::endl;
//       }

    }
    else
      fUseCalib = true;
  }

  catch(CS::Exception& e) {
    std::cerr<<e.what()<<std::endl;
  }
  catch(const std::exception &e) {
    std::cerr<<e.what()<<std::endl;
  }
  catch(...) {
    std::cout<<"PlaneMumega::ReadCalib() ==> "<<GetName()<<" calibrations, valid for ";
    std::cout<<t.tm_mday<<"."<<t.tm_mon+1<<"."<<t.tm_year+1900<<" "
	<<t.tm_hour<<":"<<t.tm_min<<":"<<t.tm_sec
	<<", not found in DB"<<std::endl;
  }

  // uncomment not to use calibrations
  // fUseCalib = false;
}
#endif //USE_DATABASE

void PlaneMumega::TextOutput(ostream& out) {

  if(ChannelProblems()) {
    out<<"\tmissing";
//     out<<fMissingChans<<std::endl;
    for (std::set<int>::iterator it=fMissingChans.begin(); it!=fMissingChans.end(); it++)
      out<<" "<<(*it);
    out<<std::endl;
    out<<"\tnoisy";
//     out<<fNoisyChans<<std::endl;
    for (std::set<int>::iterator it=fNoisyChans.begin(); it!=fNoisyChans.end(); it++)
      out<<" "<<(*it);
    out<<std::endl;
  }
  else {
    out<<"\tmissing ***"<<std::endl;
    out<<"\tnoisy ***"<<std::endl;
  }
  if(float gain = Gain())
    out<<"\tgain "<<gain<<std::endl;
  else
    out<<"\tgain ***"<<std::endl;
}


bool PlaneMumega::ChannelProblems() {
  // fHhch must have at least one entry otherwise fit will fail
  if (fHhch->GetEntries()==0)
    return false;

  const int minaverage = 30;

  // missing channels are found in hit profile histogram, if we have enough stat
  TF1 pol0("pol0","pol0");
  fHhch->Fit(&pol0,"WNQ");
  float average = pol0.GetParameter(0);
  if(average<minaverage) return false;

  fMissingChans.clear();
  fNoisyChans.clear();

  float min = average/20;
  float max = average*5;

  for(int i=1; i<=fHhch->GetNbinsX(); i++) {
    if(fHhch->GetBinContent(i)<min) fMissingChans.insert(i-1);
    if(fHhch->GetBinContent(i)>max) fNoisyChans.insert(i-1);
  }

  return true;
}

float PlaneMumega::Gain() {

  static const int statcut = 1000;
  if(fHcampsig->GetEntries() < statcut) return 0;

  fHcampsig->Fit("landau","Q0");
  return fHcampsig->GetFunction("landau")->GetParameter(1);
}




