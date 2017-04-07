#ifndef __PlanePMumega__
#define __PlanePMumega__

#ifndef __CINT__
#include "config.h"
#endif

#include <string>
#include <vector>
#include <map>
#include <set>

#include "PlaneAPV.h"
#include "PlanePanel.h"
#include "TThread.h"

#define MAX_STAT_MumegaAPV 6144
#define MAX_APV_NBCH 128
#define MAX_NBAPV_MumegaAPV 10
#define MAX_AMP_DEVIATION 20.

/// Plane for Pixel Micromegas detectors (readout with APV chips)

class PlanePMumega : public PlaneAPV {

 protected:

  typedef struct { double sum, sum2; int nb; short xpx, ypx; } stat_s;
  typedef struct { int flag; float pedestal, sigma; } calib_s;
  typedef struct { double amp[MAX_APV_NBCH]; } ampS_s;
  typedef struct { ampS_s ampSort[3]; int ampFlag[MAX_APV_NBCH]; } ampS3;

  ampS3 ampTot[MAX_NBAPV_MumegaAPV];

  // TCS phase
  float TCSphase;
  double fCMapv[MAX_NBAPV_MumegaAPV][3];
  
  // hit histograms
  TH1F *fHoccup;      // occupancy
  TH1F *fHa1a2Ratio;
  TH1F *fHa1a2RatioCM;
  TH1D *fHavgamp;
  TH1D *fHsigamp;
  TH1D *fHavgampCM;
  TH1D *fHsigampCM;
  TH1D *fHhitVsPix, *fHoccupVsPix;
  TH2F *fHampRatio;   // amplitude ratio
  TH2F *fHchvsa2;     // channel vs amplitude
  TH2F *fHchvscmc2;   // channel vs common mode correction
  TH2F *fHa0s;
  TH2F *fHa1s;
  TH2F *fHa2s;
  TH2F *fHa1a2s;
  TH2F *fHtvsTCSph;
  TH2F *fHtimeRatio;
  TH2F *fHchvsa2ped;
  TH2F *fHchvsa2CM;
  TH2F *fHa1a2CMs;


  stat_s fStat[MAX_STAT_MumegaAPV];
  stat_s fStatCM[MAX_STAT_MumegaAPV];


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

 private:
 
  int pixeltype;
  bool fPixelMM;
  bool fPixelMMsimplified;
  // hit histograms
  TH2F *fHhitMapPix;                // hit map in case of pixel plane
  TH2F *fHoccupPix;                 // occupancy in case of pixel plane
  TH2F *fHhitMapPixNorm;            // normalized hit map in case of pixel plane
  TH2F *fHoccupPixNorm;             // normalized occupancy in case of pixel plane

  // cluster variables
  Variable *fVcPosPixX, *fVcPosPixY; // cluster position in case of pixel plane

  // cluster histograms
  TH2F *fHcPosPix, *fHcPosPixT;     // position in case of pixel plane

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
  std::map<unsigned int, PlanePMumega::APVcalib> calib_data;

  std::vector<float> calib_time;

  // arrays with calibration data
  calib_s *calib_arr;

#endif //USE_DATABASE

  PlanePMumega(const char *detname, int nchan, int center, int width, int pixeltype);

  virtual void Reset();

  void ResetHistograms();

#ifndef  __CINT__
  virtual void EndEvent(const CS::DaqEvent &event);
  
#endif

  virtual void Init(TTree* tree = 0);

  virtual void Clusterize();

  #if USE_DATABASE == 1
  /// method to read calibration data
  void ReadCalib(const tm& t);
  #endif

  virtual void TextOutput(ostream& out);
  float Gain();

  void CalcCommonModeCorrection(double cut);

#ifndef  __CINT__
  void StoreDigit(CS::Chip::Digit* digit);
#endif

#ifndef __CINT__
  void StoreDigit(float x, float y, int channel, int amp0, int amp1, int amp2, int cmc0, int cmc1, int cmc2, int
  chip, int pix, int connb);
# endif

  ClassDef (PlanePMumega,0)
};

#if USE_DATABASE == 1
istream& operator>>(istream& in, std::map<unsigned int, PlanePMumega::APVcalib> &c);
#endif
#endif
