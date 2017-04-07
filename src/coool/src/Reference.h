#ifndef __REFERENCE_H__
#define __REFERENCE_H__

#include <vector>
#include <TH1.h>
#include <TProfile.h>
#include <stdio.h>
#include <string>
#include <map>

class ReferenceDirectory;

class ReferenceDirectory {
private:
  /// file for reference histos
  /// and a flag that complains if the file cant be opened
  TFile *referenceFile;
  std::string referenceFileName;
  bool referenceFileDoesntOpen;
  static std::map<std::string,ReferenceDirectory*> mapRefDir;

  /// private constructor, the new ReferenceDirectory must be taken from GiveReferenceDir
  ReferenceDirectory(const char *fileName) : referenceFile(NULL),referenceFileName(fileName),referenceFileDoesntOpen(false) {
    LoadReferenceFile(referenceFileName.c_str());}

  /// load the reference file
  void LoadReferenceFile(const char *fileName);

public:

  /// give an existing or new ReferenceDirectory corresponding to the fileName
  static ReferenceDirectory* GiveReferenceDir(const char *fileName);

  /// return the reference dir
  TDirectory* GiveReferenceFile()   { return (TDirectory*) referenceFile; }

  /// give the name of the open file
  const char *GetRefFileName() {return referenceFileName.c_str();}
};

class RefHist : public TH1F
{
  float*   reference;       //!
  float    maximum;         //!
  TH1F*    main;            //!
  int*     fCounter;        //!
  bool     isOccupHisto;    //!
  bool     isNewReference;  //!
public:
  RefHist(const char* name, Int_t nbinsx, Axis_t xlow, Axis_t xup,
	  TH1F* mainHist, int& counter, bool isOccup);
  ~RefHist() { delete reference; }
  void Paint(Option_t* option);
  void SetReference(int nEvents, TH1F* rhist);
  float* GetRef(void){ return reference;}
protected:
  void Rescale();
  //ClassDef(RefHist,0)
};

class TH1F_Ref : public TH1F
{
  RefHist* ref;  //!
  static int *default_fCounter;
  static bool showRef;
public:
  TH1F_Ref(const char* name, const char* title, Int_t nbinsx, Axis_t xlow, Axis_t xup,
	   int& counter, bool isOccup = false);
//   TH1F_Ref(const char* name, const char* title, Int_t nbinsx, Axis_t xlow, Axis_t xupr, bool isOccup = false);
  TH1F_Ref();
  ~TH1F_Ref();
  void Draw(Option_t* option = "");
  void SetReference(TDirectory* refDir, const char* histName);
  void SetReference(TDirectory* refDir);
  static void SetShowRef(bool state)  { showRef = state; }
  static bool GetShowRef()            { return showRef; }
  void SetBins(Int_t nx, Axis_t xmin, Axis_t xmax);
  virtual Double_t Compare();		// compares hist with reference, smaller number means worse
  const float* GetRef(void) const {if(ref!=NULL) return ref->GetRef(); else return NULL;}
protected:
  TH1F* FindHisto(const char* name);

// Please leave the following line commented !!! Otherwise TH1F_Ref objects are
// no more seen as TH1F in root files and they become not usable (in particular
// with the logbook)
//
//  ClassDef(TH1F_Ref,0) this must be kept commented !!!
};

class TProfile_Ref : public TProfile
{
  TProfile* ref;  //!
public:
  TProfile_Ref(const char* name, const char* title, Int_t nbinsx, Axis_t xlow, Axis_t xup);
  TProfile_Ref();
  ~TProfile_Ref();
  void Draw(Option_t* option = "");
  void Paint(Option_t* option = "");
  void SetReference(TDirectory* refDir, const char* histName);
  void SetReference(TDirectory* refDir);
  static void SetShowRef(bool state)  { TH1F_Ref::SetShowRef(state); }
  static bool GetShowRef()            { return TH1F_Ref::GetShowRef(); }
  void SetBins(Int_t nx, Axis_t xmin, Axis_t xmax);
  virtual Double_t Compare();		// compares hist with reference, smaller number means worse
  const double* GetRef(void) const { if (ref==0) return 0; else return ref->GetArray(); }
protected:
  TH1* FindHisto(const char* name);

// Please leave the following line commented !!! Otherwise TProfile_Ref objects are
// no more seen as TProfile in root files and they become not usable (in particular
// with the logbook)
//
//  ClassDef(TProfile_Ref,0) this must be kept commented !!!
};


#endif // #ifndef __REFERENCE_H__




