#include "MMan_2triplets.h"
#include "Fits.h"
#include "MManGeometry.h"

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

void MMan::TreatEvent() {

  if(!fAllocated) {
    cerr<<"System not allocated... Sorry"<<endl;
    return;
  }

  fTracker->Reset();

  if(fHitsInTree) {

//      for(int i=0; i<MM01X1___hhits; i++) {
//        if(MM01X1___ht[i]>-11000 && MM01X1___ht[i]<-9500)
//  	fTracker->AddDigit(2,MM01X1___hch[i],MM01X1___ht[i],MM01X1___htot[i]);
//      }
//      for(int i=0; i<MM01Y1___hhits; i++) {
//        if(MM01Y1___ht[i]>-11000 && MM01Y1___ht[i]<-9500)
//  	fTracker->AddDigit(3,MM01Y1___hch[i],MM01Y1___ht[i],MM01Y1___htot[i]);
//      }
//      for(int i=0; i<MM01U1___hhits; i++) {
//        if(MM01U1___ht[i]>-11000 && MM01U1___ht[i]<-9500)
//  	fTracker->AddDigit(1,MM01U1___hch[i],MM01U1___ht[i],MM01U1___htot[i]);
//      }
//      for(int i=0; i<MM01V1___hhits; i++) {
//        if(MM01V1___ht[i]>-11000 && MM01V1___ht[i]<-9500)
//  	fTracker->AddDigit(0,MM01V1___hch[i],MM01V1___ht[i],MM01V1___htot[i]);
//      }
//      for(int i=0; i<MM03X1___hhits; i++) {
//        if(MM03X1___ht[i]>-11000 && MM03X1___ht[i]<-9500)
//  	fTracker->AddDigit(6,MM03X1___hch[i],MM03X1___ht[i],MM03X1___htot[i]);
//      }
//      for(int i=0; i<MM03Y1___hhits; i++) {
//        if(MM03Y1___ht[i]>-11000 && MM03Y1___ht[i]<-9500)
//  	fTracker->AddDigit(7,MM03Y1___hch[i],MM03Y1___ht[i],MM03Y1___htot[i]);
//      }
//      for(int i=0; i<MM03U1___hhits; i++) {
//        if(MM03U1___ht[i]>-11000 && MM03U1___ht[i]<-9500)
//  	fTracker->AddDigit(5,MM03U1___hch[i],MM03U1___ht[i],MM03U1___htot[i]);
//      }
//      for(int i=0; i<MM03V1___hhits; i++) {
//        if(MM03V1___ht[i]>-11000 && MM03V1___ht[i]<-9500)
//  	fTracker->AddDigit(4,MM03V1___hch[i],MM03V1___ht[i],MM03V1___htot[i]);
//      }
  }
  else {
//      for(int i=0; i<MM01X1___hhits; i++) {
//  	fTracker->AddCluster(103,MM01X1___hch[i],MM01X1___ht[i],
//  			     MM01X1___htot[i], MM01X1___hs[i]);
//      }
//      for(int i=0; i<MM01Y1___hhits; i++) {
//  	fTracker->AddCluster(3,MM01Y1___hch[i],MM01Y1___ht[i],
//  			     MM01Y1___htot[i], MM01Y1___hs[i]);
//      }
//      for(int i=0; i<MM01U1___hhits; i++) {
//  	fTracker->AddCluster(102,MM01U1___hch[i],MM01U1___ht[i],
//  			     MM01U1___htot[i], MM01U1___hs[i]);
//      }
//      for(int i=0; i<MM01V1___hhits; i++) {
//  	fTracker->AddCluster(101,MM01V1___hch[i],MM01V1___ht[i],
//  			     MM01V1___htot[i], MM01V1___hs[i]);
//      }
//      for(int i=0; i<MM03X1___hhits; i++) {
//  	fTracker->AddCluster(111,MM03X1___hch[i],MM03X1___ht[i],
//  			     MM03X1___htot[i], MM03X1___hs[i]);
//      }
//      for(int i=0; i<MM03Y1___hhits; i++) {
//  	fTracker->AddCluster(7,MM03Y1___hch[i],MM03Y1___ht[i],
//  			     MM03Y1___htot[i], MM03Y1___hs[i]);
//      }
//      for(int i=0; i<MM03U1___hhits; i++) {
//  	fTracker->AddCluster(110,MM03U1___hch[i],MM03U1___ht[i],
//  			     MM03U1___htot[i], MM03U1___hs[i]);
//      }
//      for(int i=0; i<MM03V1___hhits; i++) {
//  	fTracker->AddCluster(109,MM03V1___hch[i],MM03V1___ht[i],
//  			     MM03V1___htot[i], MM03V1___hs[i]);
//      }

    for(int i=0; i<PB01U1___chits; i++) {
      fTracker->AddClusterWRS(3012,PB01U1___cch[i],0,0);
    }
    for(int i=0; i<PB01X1___chits; i++) {
      fTracker->AddClusterWRS(3011,PB01X1___cch[i],0,0);
    }
    for(int i=0; i<PB02V1___chits; i++) {
      fTracker->AddClusterWRS(3021,PB02V1___cch[i],0,0);
    }    
    for(int i=0; i<PB03U1___chits; i++) {
      fTracker->AddClusterWRS(3032,PB03U1___cch[i],0,0);
    }
    for(int i=0; i<PB03X1___chits; i++) {
      fTracker->AddClusterWRS(3031,PB03X1___cch[i],0,0);
    }
    for(int i=0; i<PB04V1___chits; i++) {
      fTracker->AddClusterWRS(3041,PB04V1___cch[i],0,0);
    }
    
    for(int i=0; i<PB05U1___chits; i++) {
      fTracker->AddClusterWRS(3052,PB05U1___cch[i],0,0);
    }
    for(int i=0; i<PB05X1___chits; i++) {
      fTracker->AddClusterWRS(3051,PB05X1___cch[i],0,0);
    }
    for(int i=0; i<PB06V1___chits; i++) {
      fTracker->AddClusterWRS(3061,PB06V1___cch[i],0,0);
    }    
  }
  fTracker->Clusterize(fMult);


//    fTracker->GetDoublet(0)->Point();
//    fTracker->GetDoublet(1)->Point();
//    fTracker->GetDoublet(2)->Point();
//    fTracker->GetDoublet(3)->Point();
  //fTracker->GetStation(0)->Coinc3(0);
  //fTracker->GetStation(1)->Coinc3(0);
  fTracker->Candidates();
  //fTracker->Residuals();

  if(fDoTree) 
    fOutTree->Fill();
  //fTracker->CheckCandidates(1000,0);
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
  Loop(n);
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

void MMan::Efficiency(const char* chamber,
		      const char* outfile, const char* opt) {
  
  //static bool firstpass = true;

  string sc(chamber);

  vector<Tracker*>& planes = fTracker->GetTrackers();
  
  // all planes activated
  for(unsigned i=0; i<planes.size();i++) {
    planes[i]->Activate(1);
  }

  ofstream out(outfile);

  TCanvas *efficiencies=0;
  if(strcmp(opt,"batch")) {
    efficiencies=new TCanvas("efficencies","efficencies",800,600);
    if(sc.empty())
      efficiencies->Divide(3,2);
  }

  vector<TH1F*> residuals;

  // now, each plane is deactivated for efficiency calculation
  for(unsigned i=0; i<planes.size();i++) {

    if(!sc.empty() && sc!=planes[i]->GetName()) continue;
    fTracker->ResetHistos();

    planes[i]->Activate(0);
    //    planes[i]->SetRMin(3);
    //planes[i]->SetCheckRange(-10000,10000);
    Loop(0);    
    planes[i]->Activate(1);
    out<<planes[i]->GetName()<<"\t"<<planes[i]->GetEfficiency()<<"\t"
       <<planes[i]->GetDEfficiency()<<"\t"
       <<planes[i]->GetNTot()<<"\t"
       <<planes[i]->GetNEff()<<"\t"
       <<planes[i]->GetNBackg()<<endl;
    
    if(i<6 && strcmp(opt,"batch")) {
      if(sc.empty())
	efficiencies->cd(i+1);

      (TPad::Pad())->SetLogy();
      residuals.push_back(new TH1F(*(planes[i]->fHres)));

      //   if(firstpass) {
//  	cout<<"Determining cuts, restart for an efficiency measurement"<<endl;
//  	float mean; float sigma;
//  	FitResiduals(residuals.back(),sigma,mean);
//  	residuals.back()->Draw();
//  	planes[i]->SetCheckRange(mean-16*sigma,mean+16*sigma);
	// }
	//else 
	//residuals.back()->Draw();

      residuals.back()->Draw();
      planes[i]->DrawResCuts();
      efficiencies->Modified();
      efficiencies->Update();
    }
  }  

//    if(firstpass)
//      firstpass = false;
//    else 
//      firstpass = true;
    
  out.close();
}

void MMan::EffVsChi2(const char* outfile,float chi2min, float chi2max, float step) {
  
  ofstream out(outfile);


  float chi2=chi2min;
  while(chi2<chi2max) {
    fTracker->SetPChi2(chi2);

    vector<Tracker*>& planes = fTracker->GetTrackers();
    
    // all planes activated
    for(unsigned i=0; i<planes.size();i++) {
      planes[i]->Activate(1);
    }
    
    // now, each plane is deactivated for efficiency calculation
    for(unsigned i=0; i<planes.size();i++) {
      fTracker->ResetHistos();
      planes[i]->Activate(0);
      Loop(20000);    
      planes[i]->Activate(1);
      out<<chi2<<"\t"<<planes[i]->GetName()<<"\t"<<planes[i]->GetEfficiency()<<endl;
    }  
    chi2 += step;
  }

  out.close();
}

void MMan::TimeCuts(int n) {
  vector<Tracker*>& planes = fTracker->GetTrackers();

  for(unsigned i=0; i<planes.size(); i++) {
    string histname = planes[i]->GetName(); 
    histname += "___ct";
    TH1F *h=(TH1F*)fFile->Get(histname.c_str());
    if(h) {
      float sigma; float mean;
      FitCTime(h,sigma, mean);
      planes[i]->SetCTimeRange(mean-n*sigma,mean+n*sigma);
    }
    else {
      cerr<<"MMan::TimeCuts() Cannot find histogram "<<histname<<". Using : tmin = "<<planes[i]->GetTMin()<<", tmax = "<<planes[i]->GetTMax()<<endl;
    }
  }
}

void MMan::Loop(Int_t nentries)
{

   if (fChain == 0) return;

   if (!nentries) nentries = Int_t(fChain->GetEntries());

   Int_t nbytes = 0, nb = 0;
   for (Int_t jentry=0; jentry<nentries;jentry++) {

     if((jentry%1000)==0) cout<<jentry<<endl;

      LoadTree(jentry); //in case of a TChain, ientry is the entry number in the current file
      nb = fChain->GetEntry(jentry);   nbytes += nb;
      // if (Cut(ientry) < 0) continue;
      TreatEvent();
   }
}

void MMan::Allocate(int run) {


  // open config file 
  if(fDetdat == "") { 
    string configdir = "Config/";
    char file[20];
    sprintf(file,"%s%d%s","mm",run,".dat");

    string configfile = configdir + file;
    ifstream in(configfile.c_str());
    cout<<"opening config file : "<<configfile<<endl;
  
    fTracker = new TrackMaker(2,fChi2,in,fFile);
  }
  else {
    MManGeometry *geom = 0;
    try {
      geom=new MManGeometry(fDetdat.c_str());
    }  
    catch(...) {
      cerr<<"FATAL MMan::Allocate. Geometry problem"<<endl;
      delete this;
      return;
    } 
    fTracker = new TrackMaker(2,fChi2);
    // MM
//      fTracker->AddTracker(geom->GetTracker("MM01V1__"));
//      fTracker->AddTracker(geom->GetTracker("MM01U1__"));
//      fTracker->AddTracker(geom->GetTracker("MM01X1__"));
//      fTracker->AddTracker(geom->GetTracker("MM03V1__"));
//      fTracker->AddTracker(geom->GetTracker("MM03U1__"));
//      fTracker->AddTracker(geom->GetTracker("MM03X1__"));

    // PB
    fTracker->AddTracker(geom->GetTracker("PB01X1__"));
    fTracker->AddTracker(geom->GetTracker("PB01U1__"));
    fTracker->AddTracker(geom->GetTracker("PB02V1__"));
    fTracker->AddTracker(geom->GetTracker("PB03X1__"));
    fTracker->AddTracker(geom->GetTracker("PB03U1__"));
    fTracker->AddTracker(geom->GetTracker("PB04V1__"));
    fTracker->AddTracker(geom->GetTracker("PB05X1__"));
    fTracker->AddTracker(geom->GetTracker("PB05U1__"));
    fTracker->AddTracker(geom->GetTracker("PB06V1__"));

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
  fDoTree(false),fDetdat(detdat), fFile(0) {

  string dirname = dir;
  dirname += "/";

  DIR *dp;
  struct dirent *ep;

  TChain *chain = new TChain("T");

  dp = opendir (dir);
  if (dp != NULL) {
    cout<<dir<<" opened. Using files :"<<endl;
    while (ep = readdir (dp)) {
      string filename(ep->d_name);
      char runnum[10];
      sprintf(runnum,"%d",run);
      if(filename.rfind(runnum)!=-1) {
	string fullname = dirname+filename;
	if(!fFile)
	  fFile = new TFile(fullname.c_str());
	cout<<fullname<<endl;
	chain->Add(fullname.c_str());
      }
    }
    (void) closedir (dp);
    Init(chain);
  }
  else {
    cerr<<"Cannot open directory "<<dir<<endl;
    delete this;
  }
  return;  
}

MMan::MMan(TFile *file, TTree* tree,bool hitsintree, int mult, float chi2) 
  : fHitsInTree(hitsintree),fMult(mult), fAllocated(false), fChi2(chi2), 
  fDoTree(false),fFile(file) {
  
  //TTree *tree = (TTree*) file->Get("T");

  if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("../root/mm.root");
      if (!f) {
         f = new TFile("../root/mm.root");
      }
      tree = (TTree*)gDirectory->Get("T");

   }
   Init(tree);
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
   fChain->SetBranchAddress("evtinrun",&evtinrun);
   fChain->SetBranchAddress("evtinspl",&evtinspl);
   fChain->SetBranchAddress("timesec",&timesec);
   fChain->SetBranchAddress("timemic",&timemic);
   fChain->SetBranchAddress("TT_hits",&TT_hits);
   fChain->SetBranchAddress("TT_ch",TT_ch);
   fChain->SetBranchAddress("TT_t",TT_t);

   if(fHitsInTree) {

//       fChain->SetBranchAddress("MM03Y1___hhits",&MM03Y1___hhits);
//       fChain->SetBranchAddress("MM03Y1___hch",MM03Y1___hch);
//       fChain->SetBranchAddress("MM03Y1___ht",MM03Y1___ht);
//       fChain->SetBranchAddress("MM03Y1___htot",MM03Y1___htot);
     
//       fChain->SetBranchAddress("MM03X1___hhits",&MM03X1___hhits);
//       fChain->SetBranchAddress("MM03X1___hch",MM03X1___hch);
//       fChain->SetBranchAddress("MM03X1___ht",MM03X1___ht);
//       fChain->SetBranchAddress("MM03X1___htot",MM03X1___htot);
     
//       fChain->SetBranchAddress("MM03U1___hhits",&MM03U1___hhits);
//       fChain->SetBranchAddress("MM03U1___hch",MM03U1___hch);
//       fChain->SetBranchAddress("MM03U1___ht",MM03U1___ht);
//       fChain->SetBranchAddress("MM03U1___htot",MM03U1___htot);
   
//       fChain->SetBranchAddress("MM03V1___hhits",&MM03V1___hhits);
//       fChain->SetBranchAddress("MM03V1___hch",MM03V1___hch);
//       fChain->SetBranchAddress("MM03V1___ht",MM03V1___ht);
//       fChain->SetBranchAddress("MM03V1___htot",MM03V1___htot);
     
//       fChain->SetBranchAddress("MM01Y1___hhits",&MM01Y1___hhits);
//       fChain->SetBranchAddress("MM01Y1___hch",MM01Y1___hch);
//       fChain->SetBranchAddress("MM01Y1___ht",MM01Y1___ht);
//       fChain->SetBranchAddress("MM01Y1___htot",MM01Y1___htot);

//       fChain->SetBranchAddress("MM01X1___hhits",&MM01X1___hhits);
//       fChain->SetBranchAddress("MM01X1___hch",MM01X1___hch);
//       fChain->SetBranchAddress("MM01X1___ht",MM01X1___ht);
//       fChain->SetBranchAddress("MM01X1___htot",MM01X1___htot);
     
//       fChain->SetBranchAddress("MM01U1___hhits",&MM01U1___hhits);
//       fChain->SetBranchAddress("MM01U1___hch",MM01U1___hch);
//       fChain->SetBranchAddress("MM01U1___ht",MM01U1___ht);
//       fChain->SetBranchAddress("MM01U1___htot",MM01U1___htot);
     
//       fChain->SetBranchAddress("MM01V1___hhits",&MM01V1___hhits);
//       fChain->SetBranchAddress("MM01V1___hch",MM01V1___hch);
//       fChain->SetBranchAddress("MM01V1___ht",MM01V1___ht);
//       fChain->SetBranchAddress("MM01V1___htot",MM01V1___htot);
   }     
   else {
//       fChain->SetBranchAddress("MM03Y1___chits",&MM03Y1___hhits);
//       fChain->SetBranchAddress("MM03Y1___cch",MM03Y1___hch);
//       fChain->SetBranchAddress("MM03Y1___ct",MM03Y1___ht);
//       fChain->SetBranchAddress("MM03Y1___ctot",MM03Y1___htot);
//       fChain->SetBranchAddress("MM03Y1___cs",MM03Y1___hs);
     
//       fChain->SetBranchAddress("MM03X1___chits",&MM03X1___hhits);
//       fChain->SetBranchAddress("MM03X1___cch",MM03X1___hch);
//       fChain->SetBranchAddress("MM03X1___ct",MM03X1___ht);
//       fChain->SetBranchAddress("MM03X1___ctot",MM03X1___htot);
//       fChain->SetBranchAddress("MM03X1___cs",MM03X1___hs);
     
//       fChain->SetBranchAddress("MM03U1___chits",&MM03U1___hhits);
//       fChain->SetBranchAddress("MM03U1___cch",MM03U1___hch);
//       fChain->SetBranchAddress("MM03U1___ct",MM03U1___ht);
//       fChain->SetBranchAddress("MM03U1___ctot",MM03U1___htot);
//       fChain->SetBranchAddress("MM03U1___cs",MM03U1___hs);
   
//       fChain->SetBranchAddress("MM03V1___chits",&MM03V1___hhits);
//       fChain->SetBranchAddress("MM03V1___cch",MM03V1___hch);
//       fChain->SetBranchAddress("MM03V1___ct",MM03V1___ht);
//       fChain->SetBranchAddress("MM03V1___ctot",MM03V1___htot);
//       fChain->SetBranchAddress("MM03V1___cs",MM03V1___hs);
     
//       fChain->SetBranchAddress("MM01Y1___chits",&MM01Y1___hhits);
//       fChain->SetBranchAddress("MM01Y1___cch",MM01Y1___hch);
//       fChain->SetBranchAddress("MM01Y1___ct",MM01Y1___ht);
//       fChain->SetBranchAddress("MM01Y1___ctot",MM01Y1___htot);
//       fChain->SetBranchAddress("MM01Y1___cs",MM01Y1___hs);

//       fChain->SetBranchAddress("MM01X1___chits",&MM01X1___hhits);
//       fChain->SetBranchAddress("MM01X1___cch",MM01X1___hch);
//       fChain->SetBranchAddress("MM01X1___ct",MM01X1___ht);
//       fChain->SetBranchAddress("MM01X1___ctot",MM01X1___htot);
//       fChain->SetBranchAddress("MM01X1___cs",MM01X1___hs);
     
//       fChain->SetBranchAddress("MM01U1___chits",&MM01U1___hhits);
//       fChain->SetBranchAddress("MM01U1___cch",MM01U1___hch);
//       fChain->SetBranchAddress("MM01U1___ct",MM01U1___ht);
//       fChain->SetBranchAddress("MM01U1___ctot",MM01U1___htot);
//       fChain->SetBranchAddress("MM01U1___cs",MM01U1___hs);
     
//       fChain->SetBranchAddress("MM01V1___chits",&MM01V1___hhits);
//       fChain->SetBranchAddress("MM01V1___cch",MM01V1___hch);
//       fChain->SetBranchAddress("MM01V1___ct",MM01V1___ht);
//       fChain->SetBranchAddress("MM01V1___ctot",MM01V1___htot);
//       fChain->SetBranchAddress("MM01V1___cs",MM01V1___hs);
     fChain->SetBranchAddress("PB01U1___chits",&PB01U1___chits);
     fChain->SetBranchAddress("PB01U1___cch",PB01U1___cch);
     fChain->SetBranchAddress("PB01U1___ct",PB01U1___ct);
     
     fChain->SetBranchAddress("PB01X1___chits",&PB01X1___chits);
     fChain->SetBranchAddress("PB01X1___cch",PB01X1___cch);
     fChain->SetBranchAddress("PB01X1___ct",PB01X1___ct);
     
     fChain->SetBranchAddress("PB02V1___chits",&PB02V1___chits);
     fChain->SetBranchAddress("PB02V1___cch",PB02V1___cch);
     fChain->SetBranchAddress("PB02V1___ct",PB02V1___ct);
     
     fChain->SetBranchAddress("PB03U1___chits",&PB03U1___chits);
     fChain->SetBranchAddress("PB03U1___cch",PB03U1___cch);
     fChain->SetBranchAddress("PB03U1___ct",PB03U1___ct);
     
     fChain->SetBranchAddress("PB03X1___chits",&PB03X1___chits);
     fChain->SetBranchAddress("PB03X1___cch",PB03X1___cch);
     fChain->SetBranchAddress("PB03X1___ct",PB03X1___ct);
     
     fChain->SetBranchAddress("PB04V1___chits",&PB04V1___chits);
     fChain->SetBranchAddress("PB04V1___cch",PB04V1___cch);
     fChain->SetBranchAddress("PB04V1___ct",PB04V1___ct);

     fChain->SetBranchAddress("PB05U1___chits",&PB05U1___chits);
     fChain->SetBranchAddress("PB05U1___cch",PB05U1___cch);
     fChain->SetBranchAddress("PB05U1___ct",PB05U1___ct);
     
     fChain->SetBranchAddress("PB05X1___chits",&PB05X1___chits);
     fChain->SetBranchAddress("PB05X1___cch",PB05X1___cch);
     fChain->SetBranchAddress("PB05X1___ct",PB05X1___ct);
     
     fChain->SetBranchAddress("PB06V1___chits",&PB06V1___chits);
     fChain->SetBranchAddress("PB06V1___cch",PB06V1___cch);
     fChain->SetBranchAddress("PB06V1___ct",PB06V1___ct);
   }

   Notify();
}

Bool_t MMan::Notify()
{
//   called when loading a new file
//   get branch pointers

  // removing all detectors except PB
  fChain->SetBranchStatus("*",0);
  fChain->SetBranchStatus("PB*",1);
  fChain->SetBranchStatus("runno",1);
  fChain->SetBranchStatus("spill",1);
  fChain->SetBranchStatus("evtinrun",1);
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
//       b_MM03X1___hhits = fChain->GetBranch("MM03X1___hhits");    
//       b_MM03X1___hch = fChain->GetBranch("MM03X1___hch");
//       b_MM03X1___ht = fChain->GetBranch("MM03X1___ht");
//       b_MM03X1___htot = fChain->GetBranch("MM03X1___htot");
//       b_MM03U1___hhits = fChain->GetBranch("MM03U1___hhits");
//       b_MM03U1___hch = fChain->GetBranch("MM03U1___hch");
//       b_MM03U1___ht = fChain->GetBranch("MM03U1___ht");
//       b_MM03U1___htot = fChain->GetBranch("MM03U1___htot");
//       b_MM03V1___hhits = fChain->GetBranch("MM03V1___hhits");
//       b_MM03V1___hch = fChain->GetBranch("MM03V1___hch");
//       b_MM03V1___ht = fChain->GetBranch("MM03V1___ht");
//       b_MM03V1___htot = fChain->GetBranch("MM03V1___htot");
//       b_MM01Y1___hhits = fChain->GetBranch("MM01Y1___hhits");
//       b_MM01Y1___hch = fChain->GetBranch("MM01Y1___hch");
//       b_MM01Y1___ht = fChain->GetBranch("MM01Y1___ht");
//       b_MM01Y1___htot = fChain->GetBranch("MM01Y1___htot");
//       b_MM01X1___hhits = fChain->GetBranch("MM01X1___hhits");
//       b_MM01X1___hch = fChain->GetBranch("MM01X1___hch");
//       b_MM01X1___ht = fChain->GetBranch("MM01X1___ht");
//       b_MM01X1___htot = fChain->GetBranch("MM01X1___htot");
//       b_MM01U1___hhits = fChain->GetBranch("MM01U1___hhits");
//       b_MM01U1___hch = fChain->GetBranch("MM01U1___hch");
//       b_MM01U1___ht = fChain->GetBranch("MM01U1___ht");
//       b_MM01U1___htot = fChain->GetBranch("MM01U1___htot");
//       b_MM01V1___hhits = fChain->GetBranch("MM01V1___hhits");
//       b_MM01V1___hch = fChain->GetBranch("MM01V1___hch");
//       b_MM01V1___ht = fChain->GetBranch("MM01V1___ht");
//       b_MM01V1___htot = fChain->GetBranch("MM01V1___htot");
   }
   else {
     //clusters in input
     cout<<"The tree contains clusters"<<endl;

//       b_MM03X1___hhits = fChain->GetBranch("MM03X1___chits");
//       b_MM03X1___hch = fChain->GetBranch("MM03X1___cch");
//       b_MM03X1___ht = fChain->GetBranch("MM03X1___ct");
//       b_MM03X1___htot = fChain->GetBranch("MM03X1___ctot");
//       b_MM03X1___hs = fChain->GetBranch("MM03X1___cs");

//       b_MM03U1___hhits = fChain->GetBranch("MM03U1___chits");
//       b_MM03U1___hch = fChain->GetBranch("MM03U1___cch");
//       b_MM03U1___ht = fChain->GetBranch("MM03U1___ct");
//       b_MM03U1___htot = fChain->GetBranch("MM03U1___ctot");
//       b_MM03U1___hs = fChain->GetBranch("MM03U1___cs");

//       b_MM03V1___hhits = fChain->GetBranch("MM03V1___chits");
//       b_MM03V1___hch = fChain->GetBranch("MM03V1___cch");
//       b_MM03V1___ht = fChain->GetBranch("MM03V1___ct");
//       b_MM03V1___htot = fChain->GetBranch("MM03V1___ctot");
//       b_MM03V1___hs = fChain->GetBranch("MM03V1___cs");

//       b_MM01Y1___hhits = fChain->GetBranch("MM01Y1___chits");
//       b_MM01Y1___hch = fChain->GetBranch("MM01Y1___cch");
//       b_MM01Y1___ht = fChain->GetBranch("MM01Y1___ct");
//       b_MM01Y1___htot = fChain->GetBranch("MM01Y1___ctot");
//       b_MM01Y1___hs = fChain->GetBranch("MM01Y1___cs");

//       b_MM01X1___hhits = fChain->GetBranch("MM01X1___chits");
//       b_MM01X1___hch = fChain->GetBranch("MM01X1___cch");
//       b_MM01X1___ht = fChain->GetBranch("MM01X1___ct");
//       b_MM01X1___htot = fChain->GetBranch("MM01X1___ctot");
//       b_MM01X1___hs = fChain->GetBranch("MM01X1___cs");

//       b_MM01U1___hhits = fChain->GetBranch("MM01U1___chits");
//       b_MM01U1___hch = fChain->GetBranch("MM01U1___cch");
//       b_MM01U1___ht = fChain->GetBranch("MM01U1___ct");
//       b_MM01U1___htot = fChain->GetBranch("MM01U1___ctot");
//       b_MM01U1___hs = fChain->GetBranch("MM01U1___cs");

//       b_MM01V1___hhits = fChain->GetBranch("MM01V1___chits");
//       b_MM01V1___hch = fChain->GetBranch("MM01V1___cch");
//       b_MM01V1___ht = fChain->GetBranch("MM01V1___ct");
//       b_MM01V1___htot = fChain->GetBranch("MM01V1___ctot");
//       b_MM01V1___hs = fChain->GetBranch("MM01V1___cs");

     b_PB01U1___chits = fChain->GetBranch("PB01U1___chits");
     b_PB01U1___cch = fChain->GetBranch("PB01U1___cch");
     b_PB01U1___ct = fChain->GetBranch("PB01U1___ct");

     b_PB01X1___chits = fChain->GetBranch("PB01X1___chits");
     b_PB01X1___cch = fChain->GetBranch("PB01X1___cch");
     b_PB01X1___ct = fChain->GetBranch("PB01X1___ct");

     b_PB02V1___chits = fChain->GetBranch("PB02V1___chits");
     b_PB02V1___cch = fChain->GetBranch("PB02V1___cch");
     b_PB02V1___ct = fChain->GetBranch("PB02V1___ct");

     b_PB03U1___chits = fChain->GetBranch("PB03U1___chits");
     b_PB03U1___cch = fChain->GetBranch("PB03U1___cch");
     b_PB03U1___ct = fChain->GetBranch("PB03U1___ct");

     b_PB03X1___chits = fChain->GetBranch("PB03X1___chits");
     b_PB03X1___cch = fChain->GetBranch("PB03X1___cch");
     b_PB03X1___ct = fChain->GetBranch("PB03X1___ct");

     b_PB04V1___chits = fChain->GetBranch("PB04V1___chits");
     b_PB04V1___cch = fChain->GetBranch("PB04V1___cch");
     b_PB04V1___ct = fChain->GetBranch("PB04V1___ct");

     b_PB05U1___chits = fChain->GetBranch("PB05U1___chits");
     b_PB05U1___cch = fChain->GetBranch("PB05U1___cch");
     b_PB05U1___ct = fChain->GetBranch("PB05U1___ct");

     b_PB05X1___chits = fChain->GetBranch("PB05X1___chits");
     b_PB05X1___cch = fChain->GetBranch("PB05X1___cch");
     b_PB05X1___ct = fChain->GetBranch("PB05X1___ct");

     b_PB06V1___chits = fChain->GetBranch("PB06V1___chits");
     b_PB06V1___cch = fChain->GetBranch("PB06V1___cch");
     b_PB06V1___ct = fChain->GetBranch("PB06V1___ct");

   }   

   if(!fAllocated) {
     GetEntry(0);
     Allocate(runno);
   }

   return kTRUE;
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


