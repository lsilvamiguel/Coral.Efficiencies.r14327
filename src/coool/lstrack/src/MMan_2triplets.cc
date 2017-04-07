#include "MMan_2triplets.h"
#include "Fits.h"
#include "MManGeometry.h"
#include <iostream>
#include <stddef.h>
#include <stdio.h>
#include <sys/types.h>
#include <dirent.h>

#include <string>

#include "TH2.h"
#include "TF1.h"
#include "TStyle.h"
#include "TCanvas.h"
#include "TText.h"



ClassImp(MMan);

MMan *ThisMM = 0;

void wrap(int &npar, double *gin, double &f, double *X, int iflag);

bool MMan::TreatEvent() {

  if(!fAllocated) {
    std::cerr<<"System not allocated... Sorry"<<std::endl;
    return false;
  }

  fTracker->Reset();

  if(fHitsInTree) {

    for(int i=0; i<MM01X1___hhits; i++) {
      if(MM01X1___ht[i]>-11000 && MM01X1___ht[i]<-9500)
	fTracker->AddDigit(2,MM01X1___hch[i],MM01X1___ht[i],MM01X1___htot[i]);
    }
    for(int i=0; i<MM01Y1___hhits; i++) {
      if(MM01Y1___ht[i]>-11000 && MM01Y1___ht[i]<-9500)
	fTracker->AddDigit(3,MM01Y1___hch[i],MM01Y1___ht[i],MM01Y1___htot[i]);
    }
    for(int i=0; i<MM01U1___hhits; i++) {
      if(MM01U1___ht[i]>-11000 && MM01U1___ht[i]<-9500)
	fTracker->AddDigit(1,MM01U1___hch[i],MM01U1___ht[i],MM01U1___htot[i]);
    }
    for(int i=0; i<MM01V1___hhits; i++) {
      if(MM01V1___ht[i]>-11000 && MM01V1___ht[i]<-9500)
	fTracker->AddDigit(0,MM01V1___hch[i],MM01V1___ht[i],MM01V1___htot[i]);
    }
    for(int i=0; i<MM03X1___hhits; i++) {
      if(MM03X1___ht[i]>-11000 && MM03X1___ht[i]<-9500)
	fTracker->AddDigit(6,MM03X1___hch[i],MM03X1___ht[i],MM03X1___htot[i]);
    }
    for(int i=0; i<MM03Y1___hhits; i++) {
      if(MM03Y1___ht[i]>-11000 && MM03Y1___ht[i]<-9500)
	fTracker->AddDigit(7,MM03Y1___hch[i],MM03Y1___ht[i],MM03Y1___htot[i]);
    }
    for(int i=0; i<MM03U1___hhits; i++) {
      if(MM03U1___ht[i]>-11000 && MM03U1___ht[i]<-9500)
	fTracker->AddDigit(5,MM03U1___hch[i],MM03U1___ht[i],MM03U1___htot[i]);
    }
    for(int i=0; i<MM03V1___hhits; i++) {
      if(MM03V1___ht[i]>-11000 && MM03V1___ht[i]<-9500)
	fTracker->AddDigit(4,MM03V1___hch[i],MM03V1___ht[i],MM03V1___htot[i]);
    }
  }
  else {
    for(int i=0; i<MM01X1___hhits; i++) {
      fTracker->AddClusterWRS(103,MM01X1___hch[i],MM01X1___ht[i],
			      MM01X1___hs[i], MM01X1___hres[i]);
    }
    for(int i=0; i<MM01Y1___hhits; i++) {
      fTracker->AddClusterWRS(104,MM01Y1___hch[i],MM01Y1___ht[i],
			      MM01Y1___hs[i], MM01Y1___hres[i]);
    }
    for(int i=0; i<MM01U1___hhits; i++) {
      fTracker->AddClusterWRS(102,MM01U1___hch[i],MM01U1___ht[i],
			      MM01U1___hs[i], MM01U1___hres[i]);
    }
    for(int i=0; i<MM01V1___hhits; i++) {
      fTracker->AddClusterWRS(101,MM01V1___hch[i],MM01V1___ht[i],
			      MM01V1___hs[i], MM01V1___hres[i]);
    }
    for(int i=0; i<MM02X1___hhits; i++) {
      fTracker->AddClusterWRS(107,MM02X1___hch[i],MM02X1___ht[i],
			      MM02X1___hs[i], MM02X1___hres[i]);
    }
    for(int i=0; i<MM02Y1___hhits; i++) {
      fTracker->AddClusterWRS(108,MM02Y1___hch[i],MM02Y1___ht[i],
			      MM02Y1___hs[i], MM02Y1___hres[i]);
    }
    for(int i=0; i<MM02U1___hhits; i++) {
      fTracker->AddClusterWRS(106,MM02U1___hch[i],MM02U1___ht[i],
			      MM02U1___hs[i], MM02U1___hres[i]);
    }
    for(int i=0; i<MM02V1___hhits; i++) {
      fTracker->AddClusterWRS(105,MM02V1___hch[i],MM02V1___ht[i],
			      MM02V1___hs[i], MM02V1___hres[i]);
    }
    for(int i=0; i<MM03X1___hhits; i++) {
      fTracker->AddClusterWRS(111,MM03X1___hch[i],MM03X1___ht[i],
			      MM03X1___hs[i], MM03X1___hres[i]);
    }
    for(int i=0; i<MM03Y1___hhits; i++) {
      fTracker->AddClusterWRS(112,MM03Y1___hch[i],MM03Y1___ht[i],
			      MM03Y1___hs[i], MM03Y1___hres[i]);
    }
    for(int i=0; i<MM03U1___hhits; i++) {
      fTracker->AddClusterWRS(110,MM03U1___hch[i],MM03U1___ht[i],
      			      MM03U1___hs[i], MM03U1___hres[i]);
    }
    for(int i=0; i<MM03V1___hhits; i++) {
      fTracker->AddClusterWRS(109,MM03V1___hch[i],MM03V1___ht[i],
			      MM03V1___hs[i], MM03V1___hres[i]);
    }
  }
  fTracker->Clusterize(fMult) ;
  int ncombok = fTracker->Candidates();
    
  return ncombok;
}


void MMan::Test() {
  fTracker->GetTracker(0)->AddCluster(512,0,1,1);
  fTracker->GetTracker(1)->AddCluster(512,0,1,1);
  fTracker->GetTracker(2)->AddCluster(576,0,1,1);
  fTracker->GetTracker(4)->AddCluster(512,0,1,1);
  fTracker->GetTracker(5)->AddCluster(512,0,1,1);
  fTracker->GetTracker(6)->AddCluster(512,0,1,1);
  fTracker->Clusterize(0);
  fTracker->Candidates();
  fTracker->Draw(1,1,1,1);
}

void MMan::PlotTracking(int n) {
  TCanvas *tracking = new TCanvas("tracking","tracking", 800,600);

  TText *filename=new TText(0.01,0.01,fFile->GetName());
  tracking->Divide(1,2);
  fTracker->ResetHistos(); 
  Loop(n,0);
  tracking->cd(1);
  TPad *mpad = dynamic_cast<TPad*> (TPad::Pad());
  mpad->Divide(4,1);
  mpad->cd(1);
  fTracker->PlotProbchi2();
  mpad->cd(2);
  fTracker->PlotProfile();
  mpad->cd(3);
  fTracker->PlotYp();
  mpad->cd(4);
  fTracker->PlotZp();
  tracking->cd(2);
  fTracker->PlotCTimes();
  filename->Draw();
}

TH1F* MMan::PlotResiduals(int id,int n,const char *opt) {
  TPad *mpad = dynamic_cast<TPad*> (TPad::Pad());

  if(!mpad) {
    mpad = new TCanvas();
  }
  mpad->cd();
  if(strcmp("same",opt))  mpad->Clear();

  if(Tracker *tracker=fTracker->GetTracker(id)) {
    tracker->Activate(0);
    Loop(n,0);
    TH1F *h=new TH1F(*(tracker->fHres));
    if(strcmp("",opt)) h->SetLineColor(4);
    h->Draw(opt);
    mpad->Modified();
    mpad->Update();
    tracker->Activate(1);
    return tracker->fHres;
  }
  else {
    std::cerr<<"no tracker with id "<<id<<std::endl;
    return 0;
  }
}

void MMan::Center(int id, int n) {
  TH1F *h = PlotResiduals(id,n);
  if(h) {
    TF1 *func=new TF1("func","pol0+gaus(1)",-0.5,0.5);
    float mean, sigma;
    //    FitCTime(h,sigma,mean,-0.5,0.5,0.5,1);
    h->Fit("gaus");
    mean = h->GetFunction("gaus")->GetParameter(1);
    h->Reset();
    if(Tracker *tracker=fTracker->GetTracker(id)) {
      tracker->MovePerp(-mean);
      PlotResiduals(id,n,"same");
    }
  }
}

void MMan::Efficiency(const char* sel,
		      const char* outfile, int nevent, const char* opt) {
  

  for(unsigned i=0;i<fHresSave.size(); i++)
    delete fHresSave[i];
  fHresSave.clear();

  for(unsigned i=0;i<fHproSave.size(); i++)
    delete fHproSave[i];
  fHproSave.clear();

  std::vector<Tracker*>& planes = fTracker->GetTrackers(sel);
  
  // all planes activated
  for(unsigned i=0; i<planes.size();i++) {
    planes[i]->Activate(1);
  }

  std::string datoutname = outfile;
  datoutname += ".dat";
  ofstream out(datoutname.c_str());

  std::string rootoutname = outfile;
  rootoutname += ".root";
  TDirectory *ld = gDirectory;
  TFile rootout(rootoutname.c_str(),"recreate");
  ld->cd();

  TCanvas *residuals=0;
  if(strcmp(opt,"batch")) {
    residuals=new TCanvas("residuals","residuals",800,600);
    residuals->Divide(4,3);
  }
  TCanvas *profiles=0;
  if(strcmp(opt,"batch")) {
    profiles=new TCanvas("profiles","profiles",800,600);
    profiles->Divide(4,3);
  }


  // now, each plane is deactivated for efficiency calculation
  for(unsigned i=0; i<planes.size();i++) {


    fTracker->ResetHistos();

    planes[i]->Activate(0);

    Loop(nevent,0);    

    out<<planes[i]->GetName()<<"\t"<<planes[i]->GetEfficiency()<<"\t"
       <<planes[i]->GetDEfficiency()<<"\t"
       <<planes[i]->GetNTot()<<"\t"
       <<planes[i]->GetNEff()<<"\t"
       <<planes[i]->GetNBackg()<<std::endl;
    

    if(i<12 && strcmp(opt,"batch")) {
      residuals->cd(i+1);
      TH1F *hres = new TH1F(*(planes[i]->fHres));
      fHresSave.push_back(hres);
      hres->Draw();
      //  planes[i]->DrawResCuts();
      residuals->Modified();
      residuals->Update();


      profiles->cd(i+1);
      TH2F *hpro = new TH2F(*(planes[i]->EffProfile()));
      fHproSave.push_back(hpro);
      hpro->Draw("colz");
      planes[i]->WriteEfficiency();
      profiles->Modified();
      profiles->Update();
      
      TH2F h1(*planes[i]->fHeffcor);
      TH1F h2(*planes[i]->fHres);
      TH1F h3(*planes[i]->fHtres);
      TDirectory *ld = gDirectory;
      rootout.cd();
      h1.Write();
      h2.Write();
      h3.Write();
      ld->cd();
    }
    planes[i]->Activate(1);
  }  
   
  out.close();
  rootout.Close();
}

void MMan::EffVsChi2(const char* outfile,float chi2min, float chi2max, float step) {
  
  ofstream out(outfile);


  float chi2=chi2min;
  while(chi2<chi2max) {
    fTracker->SetPChi2(chi2);

    std::vector<Tracker*>& planes = fTracker->GetTrackers();
    
    // all planes activated
    for(unsigned i=0; i<planes.size();i++) {
      planes[i]->Activate(1);
    }
    
    // now, each plane is deactivated for efficiency calculation
    for(unsigned i=0; i<planes.size();i++) {
      fTracker->ResetHistos();
      planes[i]->Activate(0);
      Loop(20000,0);    
      planes[i]->Activate(1);
      out<<chi2<<"\t"<<planes[i]->GetName()<<"\t"<<planes[i]->GetEfficiency()<<std::endl;
    }  
    chi2 += step;
  }

  out.close();
}

void MMan::TimeCuts(int n) {
  std::vector<Tracker*>& planes = fTracker->GetTrackers();

  for(unsigned i=0; i<planes.size(); i++) {
    std::string histname("MM/");
    histname += planes[i]->GetName();
    histname += "/";
    histname += planes[i]->GetName();    
    histname += "_ct";
    TH1F *h=(TH1F*)fFile->Get(histname.c_str());
    if(h) {
      float sigma; float mean;
      FitCTime(h,sigma, mean);
      //planes[i]->SetCTimeRange(mean-n*sigma,mean+n*sigma);
      planes[i]->SetCTimeRange(0,500);
    }
    else {
      std::cerr<<"MMan::TimeCuts() Cannot find histogram "<<histname<<". Using : tmin = "<<planes[i]->GetTMin()<<", tmax = "<<planes[i]->GetTMax()<<std::endl;
    }
  }
}

void MMan::Loop(Int_t nentries, bool select) {
  if (fChain == 0) return;
  
  if (!nentries) nentries = Int_t(fChain->GetEntries());
  
  Int_t nbytes = 0, nb = 0;
  
  int nsel=0;
  for (Int_t jentry=0; jentry<nentries;jentry++) {
    
    if((jentry%1000)==0) std::cout<<jentry<<std::endl;
    
    LoadTree(jentry); //in case of a TChain, ientry is the entry number in the current file
    nb = fChain->GetEntry(jentry);   nbytes += nb;
    // if (Cut(ientry) < 0) continue;
    
    if(select) {
      fEventMask.push_back(TreatEvent());
      if(fEventMask.back()) {
	nsel++;
	if(fDoTree) fOutTree->Fill();
      }
    }
    else if(!fEventMask.empty() && fEventMask[jentry] || fEventMask.empty()) {
      TreatEvent(); 
      nsel++;
      if(fDoTree) fOutTree->Fill();
    }
  }
  std::cout<<nsel<<" events selected"<<std::endl;
}

void MMan::Allocate(int run) {


  // open config file 
  if(fDetdat == "") { 
    std::string configdir = "Config/";
    char file[20];
    sprintf(file,"%s%d%s","mm",run,".dat");

    std::string configfile = configdir + file;
    ifstream in(configfile.c_str());
    std::cout<<"opening config file : "<<configfile<<std::endl;
  
    fTracker = new TrackMaker("MM",2,fChi2,in,fFile);
  }
  else {
    MManGeometry *geom = 0;
    try {
      geom=new MManGeometry(fDetdat.c_str());
    }  
    catch(const char* errmsg) {
      std::cerr<<"FATAL MMan::Allocate. Geometry problem"<<std::endl;
      std::cerr<<errmsg<<std::endl;
      delete this;
      return;
    } 
    fTracker = new TrackMaker("MM",2,fChi2);
    // MM
    fTracker->AddTracker(geom->GetTracker("MM01V1__"));
    fTracker->AddTracker(geom->GetTracker("MM01U1__"));
    fTracker->AddTracker(geom->GetTracker("MM01X1__"));
    fTracker->AddTracker(geom->GetTracker("MM01Y1__"));
    fTracker->AddTracker(geom->GetTracker("MM02V1__"));
    fTracker->AddTracker(geom->GetTracker("MM02U1__"));
    fTracker->AddTracker(geom->GetTracker("MM02X1__"));
    fTracker->AddTracker(geom->GetTracker("MM02Y1__"));
    fTracker->AddTracker(geom->GetTracker("MM03V1__"));
    fTracker->AddTracker(geom->GetTracker("MM03U1__"));
    fTracker->AddTracker(geom->GetTracker("MM03X1__"));
    fTracker->AddTracker(geom->GetTracker("MM03Y1__"));

    fTracker->Init();
  }

  if(dynamic_cast<TChain*>(fChain)) {
    TTree *tree = fChain->GetTree();
    fOutTree = tree->CloneTree(0);
  }
  else 
    fOutTree = fChain->CloneTree(0);
  
  fTracker->OutputTree(fOutTree);
  fAllocated = true;
  
  //TimeCuts(1);
}

MMan::~MMan()
{
  //  if (!fChain) return;
  // delete fChain->GetCurrentFile();

  //delete gGeometry;
  delete fTracker;
}


void MMan::Entry(Int_t i) {
  
  //fH4hs->Reset();
  GetEntry(i);
  TreatEvent();
}

MMan::MMan(const char* dir,int run, const char* detdat, bool hitsintree, int mult, float chi2) 
  : fHitsInTree(hitsintree),fMult(mult), fAllocated(false), fChi2(chi2), 
  fDoTree(false),fDetdat(detdat),  fMinuit(0), fNpar(0), fNeventOpt(0),
  fFile(0) {

  ThisMM = this;

  std::string dirname = dir;
  dirname += "/";

  DIR *dp;
  struct dirent *ep;

  TChain *chain = new TChain("T");

  dp = opendir (dir);
  if (dp != NULL) {
    std::cout<<dir<<" opened. Using files :"<<std::endl;
    while (ep = readdir (dp)) {
      std::string filename(ep->d_name);
      char runnum[10];
      sprintf(runnum,"%d",run);
      if(filename.rfind(runnum)!=-1) {
	std::string fullname = dirname+filename;
	if(!fFile)
	  fFile = new TFile(fullname.c_str());
	std::cout<<fullname<<std::endl;
	chain->Add(fullname.c_str());
      }
    }
    (void) closedir (dp);
    Init(chain);
  }
  else {
    std::cerr<<"Cannot open directory "<<dir<<std::endl;
    delete this;
  }
  return;  
}

MMan::MMan(TFile *file, TTree* tree,bool hitsintree, int mult, float chi2) 
  : fHitsInTree(hitsintree),fMult(mult), fAllocated(false), fChi2(chi2), 
  fDoTree(false),fMinuit(0),fNpar(0),fNeventOpt(0),fFile(file) {
  
  ThisMM = this;

  if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("../root/mm.root");
      if (!f) {
         f = new TFile("../root/mm.root");
      }
      tree = (TTree*)gDirectory->Get("T");

   }
   Init(tree);
}

void MMan::Reset() {
  fTracker->ResetHistos();
  fOutTree->Reset();
}

Int_t MMan::GetEntry(Int_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}

Int_t MMan::LoadTree(Int_t entry)
{
// Set the environment to read one entry
   if (!fChain) return -5;
   Int_t centry = fChain->LoadTree(entry);
   if (centry < 0) return centry;
   if (fChain->IsA() != TChain::Class()) return centry;
   TChain *chain = (TChain*)fChain;
   if (chain->GetTreeNumber() != fCurrent) {
      fCurrent = chain->GetTreeNumber();
      Notify();
   }
   return centry;
}

void MMan::Init(TTree *tree)
{
//   Set branch addresses
   if (tree == 0) return;
   fChain    = tree;
   fCurrent = -1;

   fChain->SetBranchAddress("runno",&runno);
   fChain->SetBranchAddress("spill",&spill);
   fChain->SetBranchAddress("evtinspl",&evtinspl);
   fChain->SetBranchAddress("timesec",&timesec);
   fChain->SetBranchAddress("timemic",&timemic);
   fChain->SetBranchAddress("TT_hits",&TT_hits);
   fChain->SetBranchAddress("TT_ch",TT_ch);
   fChain->SetBranchAddress("TT_t",TT_t);

   if(fHitsInTree) {

     fChain->SetBranchAddress("MM03Y1___hhits",&MM03Y1___hhits);
     fChain->SetBranchAddress("MM03Y1___hch",MM03Y1___hch);
     fChain->SetBranchAddress("MM03Y1___ht",MM03Y1___ht);
     fChain->SetBranchAddress("MM03Y1___htot",MM03Y1___htot);
     
     fChain->SetBranchAddress("MM03X1___hhits",&MM03X1___hhits);
     fChain->SetBranchAddress("MM03X1___hch",MM03X1___hch);
     fChain->SetBranchAddress("MM03X1___ht",MM03X1___ht);
     fChain->SetBranchAddress("MM03X1___htot",MM03X1___htot);
     
     fChain->SetBranchAddress("MM03U1___hhits",&MM03U1___hhits);
     fChain->SetBranchAddress("MM03U1___hch",MM03U1___hch);
     fChain->SetBranchAddress("MM03U1___ht",MM03U1___ht);
     fChain->SetBranchAddress("MM03U1___htot",MM03U1___htot);
   
     fChain->SetBranchAddress("MM03V1___hhits",&MM03V1___hhits);
     fChain->SetBranchAddress("MM03V1___hch",MM03V1___hch);
     fChain->SetBranchAddress("MM03V1___ht",MM03V1___ht);
     fChain->SetBranchAddress("MM03V1___htot",MM03V1___htot);
     
     fChain->SetBranchAddress("MM01Y1___hhits",&MM01Y1___hhits);
     fChain->SetBranchAddress("MM01Y1___hch",MM01Y1___hch);
     fChain->SetBranchAddress("MM01Y1___ht",MM01Y1___ht);
     fChain->SetBranchAddress("MM01Y1___htot",MM01Y1___htot);

     fChain->SetBranchAddress("MM01X1___hhits",&MM01X1___hhits);
     fChain->SetBranchAddress("MM01X1___hch",MM01X1___hch);
     fChain->SetBranchAddress("MM01X1___ht",MM01X1___ht);
     fChain->SetBranchAddress("MM01X1___htot",MM01X1___htot);
     
     fChain->SetBranchAddress("MM01U1___hhits",&MM01U1___hhits);
     fChain->SetBranchAddress("MM01U1___hch",MM01U1___hch);
     fChain->SetBranchAddress("MM01U1___ht",MM01U1___ht);
     fChain->SetBranchAddress("MM01U1___htot",MM01U1___htot);
     
     fChain->SetBranchAddress("MM01V1___hhits",&MM01V1___hhits);
     fChain->SetBranchAddress("MM01V1___hch",MM01V1___hch);
     fChain->SetBranchAddress("MM01V1___ht",MM01V1___ht);
     fChain->SetBranchAddress("MM01V1___htot",MM01V1___htot);
   }     
   else {
     fChain->SetBranchAddress("MM03Y1___chits",&MM03Y1___hhits);
     fChain->SetBranchAddress("MM03Y1___cch",MM03Y1___hch);
     fChain->SetBranchAddress("MM03Y1___ct",MM03Y1___ht);
     fChain->SetBranchAddress("MM03Y1___ctot",MM03Y1___htot);
     fChain->SetBranchAddress("MM03Y1___cs",MM03Y1___hs);
     fChain->SetBranchAddress("MM03Y1___cres",MM03Y1___hres);
     
     fChain->SetBranchAddress("MM03X1___chits",&MM03X1___hhits);
     fChain->SetBranchAddress("MM03X1___cch",MM03X1___hch);
     fChain->SetBranchAddress("MM03X1___ct",MM03X1___ht);
     fChain->SetBranchAddress("MM03X1___ctot",MM03X1___htot);
     fChain->SetBranchAddress("MM03X1___cs",MM03X1___hs);
     fChain->SetBranchAddress("MM03X1___cres",MM03X1___hres);
     
     fChain->SetBranchAddress("MM03U1___chits",&MM03U1___hhits);
     fChain->SetBranchAddress("MM03U1___cch",MM03U1___hch);
     fChain->SetBranchAddress("MM03U1___ct",MM03U1___ht);
     fChain->SetBranchAddress("MM03U1___ctot",MM03U1___htot);
     fChain->SetBranchAddress("MM03U1___cs",MM03U1___hs);
     fChain->SetBranchAddress("MM03U1___cres",MM03U1___hres);
   
     fChain->SetBranchAddress("MM03V1___chits",&MM03V1___hhits);
     fChain->SetBranchAddress("MM03V1___cch",MM03V1___hch);
     fChain->SetBranchAddress("MM03V1___ct",MM03V1___ht);
     fChain->SetBranchAddress("MM03V1___ctot",MM03V1___htot);
     fChain->SetBranchAddress("MM03V1___cs",MM03V1___hs);
     fChain->SetBranchAddress("MM03V1___cres",MM03V1___hres);
     
     fChain->SetBranchAddress("MM02Y1___chits",&MM02Y1___hhits);
     fChain->SetBranchAddress("MM02Y1___cch",MM02Y1___hch);
     fChain->SetBranchAddress("MM02Y1___ct",MM02Y1___ht);
     fChain->SetBranchAddress("MM02Y1___ctot",MM02Y1___htot);
     fChain->SetBranchAddress("MM02Y1___cs",MM02Y1___hs);
     fChain->SetBranchAddress("MM02Y1___cres",MM02Y1___hres);
     
     fChain->SetBranchAddress("MM02X1___chits",&MM02X1___hhits);
     fChain->SetBranchAddress("MM02X1___cch",MM02X1___hch);
     fChain->SetBranchAddress("MM02X1___ct",MM02X1___ht);
     fChain->SetBranchAddress("MM02X1___ctot",MM02X1___htot);
     fChain->SetBranchAddress("MM02X1___cs",MM02X1___hs);
     fChain->SetBranchAddress("MM02X1___cres",MM02X1___hres);
     
     fChain->SetBranchAddress("MM02U1___chits",&MM02U1___hhits);
     fChain->SetBranchAddress("MM02U1___cch",MM02U1___hch);
     fChain->SetBranchAddress("MM02U1___ct",MM02U1___ht);
     fChain->SetBranchAddress("MM02U1___ctot",MM02U1___htot);
     fChain->SetBranchAddress("MM02U1___cs",MM02U1___hs);
     fChain->SetBranchAddress("MM02U1___cres",MM02U1___hres);
   
     fChain->SetBranchAddress("MM02V1___chits",&MM02V1___hhits);
     fChain->SetBranchAddress("MM02V1___cch",MM02V1___hch);
     fChain->SetBranchAddress("MM02V1___ct",MM02V1___ht);
     fChain->SetBranchAddress("MM02V1___ctot",MM02V1___htot);
     fChain->SetBranchAddress("MM02V1___cs",MM02V1___hs);
     fChain->SetBranchAddress("MM02V1___cres",MM02V1___hres);
     
     fChain->SetBranchAddress("MM01Y1___chits",&MM01Y1___hhits);
     fChain->SetBranchAddress("MM01Y1___cch",MM01Y1___hch);
     fChain->SetBranchAddress("MM01Y1___ct",MM01Y1___ht);
     fChain->SetBranchAddress("MM01Y1___ctot",MM01Y1___htot);
     fChain->SetBranchAddress("MM01Y1___cs",MM01Y1___hs);
     fChain->SetBranchAddress("MM01Y1___cres",MM01Y1___hres);

     fChain->SetBranchAddress("MM01X1___chits",&MM01X1___hhits);
     fChain->SetBranchAddress("MM01X1___cch",MM01X1___hch);
     fChain->SetBranchAddress("MM01X1___ct",MM01X1___ht);
     fChain->SetBranchAddress("MM01X1___ctot",MM01X1___htot);
     fChain->SetBranchAddress("MM01X1___cs",MM01X1___hs);
     fChain->SetBranchAddress("MM01X1___cres",MM01X1___hres);
     
     fChain->SetBranchAddress("MM01U1___chits",&MM01U1___hhits);
     fChain->SetBranchAddress("MM01U1___cch",MM01U1___hch);
     fChain->SetBranchAddress("MM01U1___ct",MM01U1___ht);
     fChain->SetBranchAddress("MM01U1___ctot",MM01U1___htot);
     fChain->SetBranchAddress("MM01U1___cs",MM01U1___hs);
     fChain->SetBranchAddress("MM01U1___cres",MM01U1___hres);
     
     fChain->SetBranchAddress("MM01V1___chits",&MM01V1___hhits);
     fChain->SetBranchAddress("MM01V1___cch",MM01V1___hch);
     fChain->SetBranchAddress("MM01V1___ct",MM01V1___ht);
     fChain->SetBranchAddress("MM01V1___ctot",MM01V1___htot);
     fChain->SetBranchAddress("MM01V1___cs",MM01V1___hs);
     fChain->SetBranchAddress("MM01V1___cres",MM01V1___hres);
   }

   Notify();
}

Bool_t MMan::Notify()
{
//   called when loading a new file
//   get branch pointers

  // removing all detectors except MM
  fChain->SetBranchStatus("*",0);
  fChain->SetBranchStatus("MM*",1);
  
  fChain->SetBranchStatus("runno",1);
  fChain->SetBranchStatus("spill",1);
  fChain->SetBranchStatus("evtinspl",1);
  fChain->SetBranchStatus("timesec",1);
  fChain->SetBranchStatus("timemic",1);

  b_runno = fChain->GetBranch("runno");
  b_spill = fChain->GetBranch("spill");
  b_evtinrun = fChain->GetBranch("evtinrun");
  b_evtinspl = fChain->GetBranch("evtinspl");
  b_timesec = fChain->GetBranch("timesec");
  b_timemic = fChain->GetBranch("timemic");
  b_TT_hits = fChain->GetBranch("TT_hits");
  b_TT_ch = fChain->GetBranch("TT_ch");
  b_TT_t = fChain->GetBranch("TT_t");


  if(fHitsInTree) {
    //hits in input
    b_MM03X1___hhits = fChain->GetBranch("MM03X1___hhits");    
    b_MM03X1___hch = fChain->GetBranch("MM03X1___hch");
    b_MM03X1___ht = fChain->GetBranch("MM03X1___ht");
    b_MM03X1___htot = fChain->GetBranch("MM03X1___htot");
    b_MM03U1___hhits = fChain->GetBranch("MM03U1___hhits");
    b_MM03U1___hch = fChain->GetBranch("MM03U1___hch");
    b_MM03U1___ht = fChain->GetBranch("MM03U1___ht");
    b_MM03U1___htot = fChain->GetBranch("MM03U1___htot");
    b_MM03V1___hhits = fChain->GetBranch("MM03V1___hhits");
    b_MM03V1___hch = fChain->GetBranch("MM03V1___hch");
    b_MM03V1___ht = fChain->GetBranch("MM03V1___ht");
    b_MM03V1___htot = fChain->GetBranch("MM03V1___htot");
    b_MM01Y1___hhits = fChain->GetBranch("MM01Y1___hhits");
    b_MM01Y1___hch = fChain->GetBranch("MM01Y1___hch");
    b_MM01Y1___ht = fChain->GetBranch("MM01Y1___ht");
    b_MM01Y1___htot = fChain->GetBranch("MM01Y1___htot");
    b_MM01X1___hhits = fChain->GetBranch("MM01X1___hhits");
    b_MM01X1___hch = fChain->GetBranch("MM01X1___hch");
    b_MM01X1___ht = fChain->GetBranch("MM01X1___ht");
    b_MM01X1___htot = fChain->GetBranch("MM01X1___htot");
    b_MM01U1___hhits = fChain->GetBranch("MM01U1___hhits");
    b_MM01U1___hch = fChain->GetBranch("MM01U1___hch");
    b_MM01U1___ht = fChain->GetBranch("MM01U1___ht");
    b_MM01U1___htot = fChain->GetBranch("MM01U1___htot");
    b_MM01V1___hhits = fChain->GetBranch("MM01V1___hhits");
    b_MM01V1___hch = fChain->GetBranch("MM01V1___hch");
    b_MM01V1___ht = fChain->GetBranch("MM01V1___ht");
    b_MM01V1___htot = fChain->GetBranch("MM01V1___htot");
  }
  else {
    //clusters in input
    std::cout<<"The tree contains clusters"<<std::endl;
    
    b_MM03X1___hhits = fChain->GetBranch("MM03X1___chits");
    b_MM03X1___hch = fChain->GetBranch("MM03X1___cch");
    b_MM03X1___ht = fChain->GetBranch("MM03X1___ct");
    b_MM03X1___htot = fChain->GetBranch("MM03X1___ctot");
    b_MM03X1___hs = fChain->GetBranch("MM03X1___cs");
    b_MM03X1___hres = fChain->GetBranch("MM03X1___cres");
    
    b_MM03U1___hhits = fChain->GetBranch("MM03U1___chits");
    b_MM03U1___hch = fChain->GetBranch("MM03U1___cch");
    b_MM03U1___ht = fChain->GetBranch("MM03U1___ct");
    b_MM03U1___htot = fChain->GetBranch("MM03U1___ctot");
    b_MM03U1___hs = fChain->GetBranch("MM03U1___cs");
    b_MM03U1___hres = fChain->GetBranch("MM03U1___cres");
    
    b_MM03V1___hhits = fChain->GetBranch("MM03V1___chits");
    b_MM03V1___hch = fChain->GetBranch("MM03V1___cch");
    b_MM03V1___ht = fChain->GetBranch("MM03V1___ct");
    b_MM03V1___htot = fChain->GetBranch("MM03V1___ctot");
    b_MM03V1___hs = fChain->GetBranch("MM03V1___cs");
    b_MM03V1___hres = fChain->GetBranch("MM03V1___cres");
    
    b_MM03V1___hhits = fChain->GetBranch("MM03Y1___chits");
    b_MM03V1___hch = fChain->GetBranch("MM03Y1___cch");
    b_MM03V1___ht = fChain->GetBranch("MM03Y1___ct");
    b_MM03V1___htot = fChain->GetBranch("MM03Y1___ctot");
    b_MM03V1___hs = fChain->GetBranch("MM03Y1___cs");
    b_MM03V1___hres = fChain->GetBranch("MM03Y1___cres");
    
    b_MM01Y1___hhits = fChain->GetBranch("MM01Y1___chits");
    b_MM01Y1___hch = fChain->GetBranch("MM01Y1___cch");
    b_MM01Y1___ht = fChain->GetBranch("MM01Y1___ct");
    b_MM01Y1___htot = fChain->GetBranch("MM01Y1___ctot");
    b_MM01Y1___hs = fChain->GetBranch("MM01Y1___cs");
    b_MM01Y1___hres = fChain->GetBranch("MM01Y1___cres");
    
    b_MM01X1___hhits = fChain->GetBranch("MM01X1___chits");
    b_MM01X1___hch = fChain->GetBranch("MM01X1___cch");
    b_MM01X1___ht = fChain->GetBranch("MM01X1___ct");
    b_MM01X1___htot = fChain->GetBranch("MM01X1___ctot");
    b_MM01X1___hs = fChain->GetBranch("MM01X1___cs");
    b_MM01X1___hres = fChain->GetBranch("MM01X1___cres");
    
    b_MM01U1___hhits = fChain->GetBranch("MM01U1___chits");
    b_MM01U1___hch = fChain->GetBranch("MM01U1___cch");
    b_MM01U1___ht = fChain->GetBranch("MM01U1___ct");
    b_MM01U1___htot = fChain->GetBranch("MM01U1___ctot");
    b_MM01U1___hs = fChain->GetBranch("MM01U1___cs");
    b_MM01U1___hres = fChain->GetBranch("MM01U1___cres");
    
    b_MM01V1___hhits = fChain->GetBranch("MM01V1___chits");
    b_MM01V1___hch = fChain->GetBranch("MM01V1___cch");
    b_MM01V1___ht = fChain->GetBranch("MM01V1___ct");
    b_MM01V1___htot = fChain->GetBranch("MM01V1___ctot");
    b_MM01V1___hs = fChain->GetBranch("MM01V1___cs");
    b_MM01V1___hres = fChain->GetBranch("MM01V1___cres");
    
    b_MM01Y1___hhits = fChain->GetBranch("MM02Y1___chits");
    b_MM02Y1___hch = fChain->GetBranch("MM02Y1___cch");
    b_MM02Y1___ht = fChain->GetBranch("MM02Y1___ct");
    b_MM02Y1___htot = fChain->GetBranch("MM02Y1___ctot");
    b_MM02Y1___hs = fChain->GetBranch("MM02Y1___cs");
    b_MM02Y1___hres = fChain->GetBranch("MM02Y1___cres");
    
    b_MM02X1___hhits = fChain->GetBranch("MM02X1___chits");
    b_MM02X1___hch = fChain->GetBranch("MM02X1___cch");
    b_MM02X1___ht = fChain->GetBranch("MM02X1___ct");
    b_MM02X1___htot = fChain->GetBranch("MM02X1___ctot");
    b_MM02X1___hs = fChain->GetBranch("MM02X1___cs");
    b_MM02X1___hres = fChain->GetBranch("MM02X1___cres");

    b_MM02U1___hhits = fChain->GetBranch("MM02U1___chits");
    b_MM02U1___hch = fChain->GetBranch("MM02U1___cch");
    b_MM02U1___ht = fChain->GetBranch("MM02U1___ct");
    b_MM02U1___htot = fChain->GetBranch("MM02U1___ctot");
    b_MM02U1___hs = fChain->GetBranch("MM02U1___cs");
    b_MM02U1___hres = fChain->GetBranch("MM02U1___cres");
    
    b_MM02V1___hhits = fChain->GetBranch("MM02V1___chits");
    b_MM02V1___hch = fChain->GetBranch("MM02V1___cch");
    b_MM02V1___ht = fChain->GetBranch("MM02V1___ct");
    b_MM02V1___htot = fChain->GetBranch("MM02V1___ctot");
    b_MM02V1___hs = fChain->GetBranch("MM02V1___cs");
    b_MM02V1___hres = fChain->GetBranch("MM02V1___cres");
  }   
  
  if(!fAllocated) {
    GetEntry(0);
    Allocate(runno);
  }
  
  return kTRUE;
}


void MMan::InitMinuit(int nevent,float dpos) {
  // for now, will optimize only the positions ... 

  fNeventOpt = nevent;

  int dummy;
  std::cout<<"**********      Geometry optimization      ************"<<std::endl
      <<std::endl
      <<std::endl
      <<"STEP 1 : event selection"<<std::endl;

  //Loop(nevent,1);
  TCanvas *c1=new TCanvas();
  TH1F *pchi2bck = new TH1F(*(fTracker->fHprobchi2));
  pchi2bck->Draw();
  c1->Modified();
  c1->Update();
  
  fTracker->NoCuts();
  fTracker->SetChi2(100);

  std::cout<<"STEP 2 : minuit initialization"<<std::endl;

  std::vector<Tracker*> &Trackers = fTracker->GetTrackers();
  double p0=0, p1=1, p2=2;
  
  // parameters : (y0,z0) for each plane, <y>, <z>
  int npar = Trackers.size() * 2 + 2;
  TMinuit *fMinuit = new TMinuit(npar); 
  fMinuit->mninit(5,6,7);     // Logical units
  fMinuit->SetFCN(wrap);  

  int ierflg=0;
  double mposy=0;
  double mposz=0;
  fNpar=0;
  for (unsigned i=0; i<Trackers.size(); i++) {
//      float angle = Trackers[i]->GetAngle();
//      if(fabs(fabs(angle)-90.0) > 1.0) {
//        // not a y plane
//        std::cerr<<Trackers[i]->GetName()<<" !y"<<std::endl;
//        mposy += Trackers[i]->GetY();
//        std::string name(Trackers[i]->GetName());
//        fMinuit->mnparm( fNpar++,(name+"Y0").c_str(),Trackers[i]->GetY(),dpos, 0,0,ierflg);    
//      } 

//      if(fabs(angle) > 1.0) {
//        // this is not an x plane
//        std::cerr<<Trackers[i]->GetName()<<" !x"<<std::endl;
//        mposz += Trackers[i]->GetZ();
//        std::string name(Trackers[i]->GetName());
//        fMinuit->mnparm( fNpar++,(name+"Z0").c_str(),Trackers[i]->GetZ(),dpos, 0,0,ierflg);
//      }
    //fMinuit->mnparm( fNpar++,Trackers[i]->GetName(),Trackers[i]->GetOffset(),dpos, 0,0,ierflg);
    fMinuit->mnparm( fNpar++,Trackers[i]->GetName(),0,dpos, 0,0,ierflg);
  }
  

//    mposy = mposy/Trackers.size();
//    mposz = mposz/Trackers.size();  
  fMinuit->FixParameter(0);
  fMinuit->FixParameter(1);
  fMinuit->FixParameter(fNpar-2);
  fMinuit->FixParameter(fNpar-1);
  


//    fMinuit->mnparm( ipar,"<Y0>"   ,mposy,dpos, 0,0,ierflg);
//    double pcur = ipar+1;
//    fMinuit->mnexcm("FIX", &pcur ,1,ierflg);
//    ipar++;
//    fMinuit->mnparm( ipar,"<Z0>"   ,mposz,dpos, 0,0,ierflg);
//    pcur = ipar+1;
//    fMinuit->mnexcm("FIX", &pcur ,1,ierflg);

  if (ierflg) {
    Printf(" UNABLE TO DEFINE PARAMETER NO.");
  }

  // minuit options

  std::cout<<"STEP 3 : fit. go on ?"<<std::endl;

  //fMinuit->SetMaxIterations(10);
  fMinuit->Migrad();
  //  gMinuit->mnexcm("MINOS", &p0 ,0,ierflg);
  // fMinuit->mnexcm("CALL FCN", &p1 ,1,ierflg);

  Double_t amin,edm,errdef;
  Int_t nvpar,nparx,icstat;
  fMinuit->mnstat(amin,edm,errdef,nvpar,nparx,icstat);
  fMinuit->mnprin(3,amin);
  
  std::cout<<"the end"<<std::endl;
  
  fTracker->fHprobchi2->SetLineColor(4);
  fTracker->fHprobchi2->Draw("same");
  c1->Modified();
  c1->Update();
}

void MMan::Fcn(int &npar, double *gin, double &f, double *par, int iflag){
  fTracker->ResetHistos();
  std::vector<Tracker*> &Trackers = fTracker->GetTrackers();

  int ipar=0;
  for(int i=0; i<Trackers.size(); i++) {
    if(ipar==fNpar) {
      std::cerr<<"TrackMaker::Fcn : not enough parameters !"<<std::endl;
      break;
    }
    Trackers[i]->MoveBack();
    Trackers[i]->MovePerp(par[ipar++]);
//      float angle = Trackers[i]->GetAngle();
//      if(fabs(fabs(angle)-90.0) > 1.0) {
//        // not a y plane
//        Trackers[i]->SetY(par[ipar++]);    
//      } 
//      if(fabs(angle) > 1.0) {
//        // this is not an x plane
//        Trackers[i]->SetZ(par[ipar++]);      
//      }
  }
  Loop(fNeventOpt,0);

  f = fTracker->GetMeanChi2();
  std::cout<<f<<std::endl;
}  

void MMan::Show(Int_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}

Int_t MMan::Cut(Int_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}

void wrap(int &npar, double *gin, double &f, double *X, int iflag) {
  ThisMM->Fcn(npar,gin,f,X,iflag);
}

