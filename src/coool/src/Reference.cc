
#include <iostream>

#include <TFile.h>
#include <TWebFile.h>
#include <TTree.h>
#include <TKey.h>
#include "Reference.h"

//ClassImp(RefHist);
// ClassImp(TH1F_Ref);

std::map<std::string,ReferenceDirectory*> ReferenceDirectory::mapRefDir;


// TH1F_Ref::TH1F_Ref(const char* name, const char* title,
// 		   Int_t nbinsx, Axis_t xlow, Axis_t xup, bool isOccup):
//   TH1F(name, title, nbinsx, xlow, xup)
// {
//   char refname[256];
//   strcpy(refname,name);
//   strcat(refname,"_ref");
//   ref = new RefHist(refname, nbinsx, xlow, xup, this, *default_fCounter, isOccup);
//   ref->SetLineColor(2);
//   ref->SetLineStyle(2);
//   //  ref->SetLineWidth(3);
// }


TH1F_Ref::TH1F_Ref(const char* name, const char* title,
		   Int_t nbinsx, Axis_t xlow, Axis_t xup,
		   int& counter, bool isOccup) :
  TH1F(name, title, nbinsx, xlow, xup)
{
  char refname[256];
  strcpy(refname,name);
  strcat(refname,"_ref");
  if (!default_fCounter) default_fCounter=&counter;
  ref = new RefHist(refname, nbinsx, xlow, xup, this, counter, isOccup);
  ref->SetLineColor(2);
  ref->SetLineStyle(2);
  //  ref->SetLineWidth(3);
}

int *TH1F_Ref::default_fCounter=NULL;
bool TH1F_Ref::showRef = true;


TH1F_Ref::TH1F_Ref()
{
  ref = NULL;
}


TH1F_Ref::~TH1F_Ref()
{
  //  delete ref;
}


void TH1F_Ref::Draw(Option_t* option)
{
  TH1F::Draw(option);
  if (ref && showRef) ref->Draw("same");
}


void TH1F_Ref::SetBins(Int_t nx, Axis_t xmin, Axis_t xmax)
{
  TH1F::SetBins(nx,xmin,xmax);
  ref->SetBins(nx,xmin,xmax);
}


void TH1F_Ref::SetReference(TDirectory* dir)
{
  const char *histName;
  histName=GetName();
  SetReference(dir,histName);
}


void TH1F_Ref::SetReference(TDirectory* dir, const char* histName)
{
  if (!dir) return; // no dir, so we can't find anything
  if (!showRef) return; // we don't show reference plots
  TDirectory *old_dir = gDirectory;
  dir->cd();
  TH1F* hits = FindHisto("TT_hits");
  dir->cd();
  TH1F* rhist = FindHisto(histName);
  if(rhist)
    ref->SetReference((int)(hits->GetEntries()), rhist);
  else
    {
      std::cerr<<"TH1F_Ref::SetReference: Cannot find the histogram "<<histName;
      std::cerr<<" into reference file (dir path: "<<dir->GetPath()<<" )"<<std::endl;
    }
  dir->cd();
  old_dir->cd();
}


RefHist::RefHist(const char* name, Int_t nbinsx, Axis_t xlow, Axis_t xup,
		 TH1F* mainHist, int& counter, bool isOccup) :
  TH1F(name, "Reference", nbinsx, xlow, xup)
{
  isOccupHisto = isOccup;
  fCounter = &counter;
  main = mainHist;
  reference = NULL;
  maximum = 0.;
  isNewReference = true;
}


void RefHist::Rescale()
{
  if(maximum != 0.)
    {
      if(!isOccupHisto)
	{
	  register float scale = float(*fCounter)/maximum;
	  for(register int i = 0; i < GetNbinsX(); i++)
	    SetBinContent(i+1,reference[i]*scale);
	}
      else if(isNewReference)
	{
	  isNewReference = false;
	  for(register int i = 0; i < GetNbinsX(); i++)
	    SetBinContent(i+1,reference[i]);	
	}
    }
}


void RefHist::Paint(Option_t* option)
{
  Rescale();
  TH1F::Paint(option);
}


void RefHist::SetReference(int nEvents, TH1F* rhist)
{
  if(!nEvents) return;
  isNewReference = true;
  maximum = (float)nEvents;
  reference = new float[rhist->GetNbinsX()+1];
  for(register int i = 0; i <= rhist->GetNbinsX(); i++)
    reference[i] = rhist->GetBinContent(i+1);
}


Double_t TH1F_Ref::Compare()
{
  Double_t difference=0;
  difference=KolmogorovTest(ref,"UOM");
  return difference;
}



TH1F* TH1F_Ref::FindHisto(const char* keyname)
{
  TH1F* ret = NULL;

  Short_t  cycle;
  char     name[256];

// ROOT sucks, sometimes...
#if ROOT_VERSION_CODE >= ROOT_VERSION(5,34,12)
  gDirectory->DecodeNameCycle(keyname, name, cycle, 256-1);
#else
  gDirectory->DecodeNameCycle(keyname, name, cycle);
#endif

  TIter next(gDirectory->GetListOfKeys());
  TKey *key;
  while ((key = (TKey *)next())) {
    if (!strcmp(name, key->GetName())) {
      if (cycle == 9999)             return (TH1F*)key->ReadObj();
      if (cycle >= key->GetCycle())  return (TH1F*)key->ReadObj();
    }
  }
  //try with subdirectories
  next.Reset();
  while ((key = (TKey *) next())) {
    // starting with ROOT 5.16 TDirectory seems to be an abstract class
    // everything that is called TDirectory with older ROOT versions
    // now seems to be called TDirectoryFile
    if ( strcmp(key->GetClassName(),"TDirectory")==0 || strcmp(key->GetClassName(),"TDirectoryFile")==0 ) {
      gDirectory->cd(key->GetName());
      ret = FindHisto(keyname);
      if(ret) return ret;
      else gDirectory->cd("..");
    }
  }
  return ret;
}


void ReferenceDirectory::LoadReferenceFile(const char *fileName)
{
  if (!referenceFileDoesntOpen && !referenceFile) {
    TDirectory *old_dir=gDirectory;
    //     char *fileName="http://pccoeb03.cern.ch/rootfile.php/runnb=20820";
    //    char *fileName="/afs/cern.ch/compass/detector/monitor/References/standard.root";
    if (strncmp(fileName,"http://",7)==0) {
      referenceFile=new TWebFile(fileName);
    } else {
      referenceFile= new TFile(fileName,"READ");
    }
    if(!referenceFile || !referenceFile->IsOpen()) {
      std::cout << "ReferenceDirectory::LoadReferenceFile: Cannot open reference file: " << fileName << std::endl;
      referenceFileDoesntOpen=true;
      referenceFile=NULL;
    } else {
      std::cout << "Reference file " << fileName << " opened" << std::endl;
    }
    old_dir->cd();
  }
}

ReferenceDirectory* ReferenceDirectory::GiveReferenceDir(const char *fileName)
{
  std::string fname(fileName);

  std::map<std::string,ReferenceDirectory*>::iterator i=mapRefDir.find(fname);
  if(i!=mapRefDir.end()) return i->second;

  ReferenceDirectory* refdir = new ReferenceDirectory(fname.c_str());
  mapRefDir[fname] = refdir;
  return refdir;
}

TProfile_Ref::TProfile_Ref(const char* name, const char* title,
          		   Int_t nbinsx, Axis_t xlow, Axis_t xup) :
  TProfile(name, title, nbinsx, xlow, xup)
{
  char refname[256];
  strcpy(refname,name);
  strcat(refname,"_ref");
  ref = 0;
}

TProfile_Ref::TProfile_Ref()
{
  ref = NULL;
}


TProfile_Ref::~TProfile_Ref()
{
  //  delete ref;
}


void TProfile_Ref::Draw(Option_t* option)
{
  TProfile::Draw(option);
  if (GetShowRef() && ref!=0) ref->Draw("same");
}

void TProfile_Ref::Paint(Option_t* option)
{
  if (GetShowRef() && ref!=0) {
    SetMinimum(); // ROOT can soooooo much do what you do not want it to do
    SetMaximum(); // for example return as maximum and minimum not actual values
                  // but what you set

    double min=std::min(GetMinimum(), ref->GetMinimum());
    double max=std::max(GetMaximum(), ref->GetMaximum());

    double margin=0.05*(max-min);

    SetMinimum(min-margin);
    SetMaximum(max+margin);
  }
  TProfile::Paint(option);
}

void TProfile_Ref::SetBins(Int_t nx, Axis_t xmin, Axis_t xmax)
{
  TProfile::SetBins(nx,xmin,xmax);
  if (ref!=0) ref->SetBins(nx,xmin,xmax);
}


void TProfile_Ref::SetReference(TDirectory* dir)
{
  const char *histName;
  histName=GetName();
  SetReference(dir,histName);
}


void TProfile_Ref::SetReference(TDirectory* dir, const char* histName)
{
  if (!dir) return; // no dir, so we can't find anything
  if (!GetShowRef()) return; // we don't show reference plots
  TDirectory *old_dir = gDirectory;
  dir->cd();
  TProfile* rhist = (TProfile*) FindHisto(histName);
  if(rhist) {
    ref = rhist;
    ref->SetLineColor(2);
    ref->SetLineStyle(2);
    ref->SetLineWidth(2);
  } else {
      std::cerr<<"TProfile_Ref::SetReference: Cannot find the histogram "<<histName;
      std::cerr<<" into reference file (dir path: "<<dir->GetPath()<<" )"<<std::endl;

      ref = 0;
    }
  dir->cd();
  old_dir->cd();
}

Double_t TProfile_Ref::Compare()
{
  Double_t difference=0;
  difference=KolmogorovTest(ref,"UOM");
  return difference;
}

TH1* TProfile_Ref::FindHisto(const char* keyname)
{
  TH1* ret = NULL;

  Short_t  cycle;
  char     name[256];

// ROOT sucks, sometimes...
#if ROOT_VERSION_CODE >= ROOT_VERSION(5,34,12)
  gDirectory->DecodeNameCycle(keyname, name, cycle, 256-1);
#else
  gDirectory->DecodeNameCycle(keyname, name, cycle);
#endif

  TIter next(gDirectory->GetListOfKeys());
  TKey *key;
  while ((key = (TKey *)next())) {
    if (!strcmp(name, key->GetName())) {
      if (cycle == 9999)             return (TH1*)key->ReadObj();
      if (cycle >= key->GetCycle())  return (TH1*)key->ReadObj();
    }
  }
  //try with subdirectories
  next.Reset();
  while ((key = (TKey *) next())) {
    // starting with ROOT 5.16 TDirectory seems to be an abstract class
    // everything that is called TDirectory with older ROOT versions
    // now seems to be called TDirectoryFile
    if ( strcmp(key->GetClassName(),"TDirectory")==0 || strcmp(key->GetClassName(),"TDirectoryFile")==0 ) {
      gDirectory->cd(key->GetName());
      ret = FindHisto(keyname);
      if(ret) return ret;
      else gDirectory->cd("..");
    }
  }
  return ret;
}





