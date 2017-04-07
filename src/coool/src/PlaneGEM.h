#ifndef __PlaneGEM__
#define __PlaneGEM__

#include "PlaneAPV.h"
#include "TThread.h"

// forward declaration of ROOT classes
class TProfile;

/// Plane for GEM detector (readout with APV chips).


class PlaneGEM : public PlaneAPV {
  
 protected:
  // TCS phase
  float TCSphase;
  
  // hit histograms
  TH1F*     fHoccup;      // occupancy
  TH1F*     fHsumAmp;
  TH1F*     fHmeanAmp;    // mean hit amplitude per strip
  TH2F*     fHampRatio;   // amplitude ratio
  TH2F*     fHchvsa2;     // channel vs amplitude
  TProfile* fHchipvscmc2; // channel vs common mode correction

  std::vector<TH1D*> vHchipcmc2s; // common mode corrections per chip

  // number of clusters
  int fNclustKept;

  // cluster variables
  Variable *fVcPos, *fVcA2, *fVcSize, *fVcTime;

  // cluster histograms
  TH1F *fHcHits,               *fHcHitsT;     // multiplicity
  TH1F *fHcPos,                *fHcPosT;      // position
  TH1F *fHcA0, *fHcA1, *fHcA2, *fHcA2T;       // amplitudes
  TH1F *fHcSize,               *fHcSizeT;     // size
  TH2F *fHcAmpRatio,           *fHcAmpRatioT; // amplitude ratio
  TH1F *fHcTime,               *fHcTimeT;     // time
  TH1F *fHcSumAmp,             *fHcSumAmpT;
  TH1F *fHcMeanAmp,            *fHcMeanAmpT;  // mean cluster amplitude

  bool fPrintedMissCalib;

 public:
  
  /// calibration
#if USE_DATABASE == 1
  class APVchannel {
   public:
    int ch, flag;
    float ped, sigma, calped, calsigma;
    int ldcId, srcId, adcId, chipId;
    APVchannel() : ch(0), flag(1), ped(700), sigma(4.5), calped(700), 
      calsigma(4.5), ldcId(0), srcId(0), adcId(0), chipId(0) {}
    APVchannel(int f, float p, float s, float cp, float cs,
	       int l, int sr, int ad, int ci) {
      flag=f; ped=p; sigma=s; calped=cp; calsigma=cs; ch=0;
      ldcId=l; srcId=sr, adcId=ad, chipId=ci;
    }
  };
  class APVcalib { 
   public:
    std::vector<APVchannel> channel;
  };
  std::map<unsigned int, PlaneGEM::APVcalib> calib_data;
  std::map<int, float> maxBaselines;

  std::vector<float> calib_time;
#endif //USE_DATABASE
  
  PlaneGEM(const char *detname,int nchan, int center, int width, bool pixel=false);

  virtual void Reset();
  
#ifndef  __CINT__
  virtual void EndEvent(const CS::DaqEvent &event);
#endif 

  virtual void Init(TTree* tree =0);
  
  virtual void Clusterize();
  
#if USE_DATABASE == 1
  /// method to read calibration data
  void ReadCalib(const tm& t);
#endif

  virtual void TextOutput(ostream& out);
  float Gain();

  ClassDef(PlaneGEM,0)
};

#if USE_DATABASE == 1
istream& operator>>(istream& in, std::map<unsigned int, PlaneGEM::APVcalib> &c);
#endif

#endif

