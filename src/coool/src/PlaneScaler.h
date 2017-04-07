#ifndef __PlaneScaler__
#define __PlaneScaler__

#include "TProfile.h"
#include "Plane.h"

/*! \brief class for Scalers
  \author Colin Bernet
*/
class PlaneScaler : public  Plane {
  
 private:
  static const int fDepth;
  static const int fNchan;

  bool fIsEbe;
  
  //          don't delete comment-looking "[fNchan]" (needed by rootcint)
  Int_t*  fChanMult;    //[fNchan]
  Int_t*  fLastCounts;  //[fNchan]
  Int_t*  fCounts;      //[fNchan]

  int  fPattern;
  int  fTime;
  int  fLastTime;
  Float_t fDeltaTime;

  /// times histogram
  TH1F *fHtimes;

  /// counts histogram
  TH1F *fHcounts;
  
  /// rates histogram 
  TProfile *fHrates;

  /// beam profile
  TProfile *fHSpillProfile;

  /// number of muons per spill 
  TH1F *fHIntensity;

  // number of muons per spill
  int  fBeamhitspspill;      

  // spill number
  int  fSpillnum;      


 public:
  
  /*! \brief constructor
    \param detname detector name
    \param ebe 1 if the scaler is read event by event, 0 otherwise

    \todo only event by event mode implemented right now
  */
  PlaneScaler(const char *detname,int ebe);
  virtual ~PlaneScaler();
  
#ifndef __CINT__
  void StoreDigit(CS::Chip::Digit* digit);
  void EndEvent(const CS::DaqEvent &event);
#endif
  void Init(TTree* tree =0);

  void Reset();

  void NewRun();
  void End();
 
  void ControlPanel(const TGWindow* p, const TGWindow* main); 

  int GetTime() const {return fTime;}
  float GetTimeSeconds() const {return fTime/38.88E6;}
  float GetTimeInterval() const {return fDeltaTime;}
  int GetPattern() const {return fPattern;}
  int* GetCounts() const {return fChanMult;}

  ClassDef(PlaneScaler,1)
};

#endif




