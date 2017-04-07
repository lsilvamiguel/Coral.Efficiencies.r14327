#ifndef __PlaneRiAPV__
#define __PlaneRiAPV__

#include "PlaneAPV.h"

#define MAX_STAT_RiAPV 6144
#define MAX_APV_NBCH 128
#define MAX_NBAPV_RiAPV 48
#define MAX_AMP_DEVIATION 20.

/// Plane for Rich detector APV electronics (readout with APV chips).


class PlaneRiAPV : public PlaneAPV {

 private:

  typedef struct { int x, y; } xy_s ;
  typedef struct { double sum, sum2; int nb; short xpx, ypx; } stat_s ;
  typedef struct { int flag; float pedestal, sigma; } calib_s ;
  typedef struct { double amp[MAX_APV_NBCH]; } ampS_s;
  typedef struct { ampS_s ampSort[3]; int ampFlag[MAX_APV_NBCH]; } ampS3;

  int fMinChan;
  ampS3 ampTot[MAX_NBAPV_RiAPV];

  // Channel histograms
  TH2F *fHampRatio, *fHampSubRatio, *fHtimeRatio;
  TH1F *fHa1a2Ratio;
  TH2F *fHchvsa2, *fHchvsa2ped, *fHchvsa2CM;
  TH1F *fHoccup, *fHampSubRap;
  TH1F *fHa0CM, *fHa1CM, *fHa2CM;
  TH2F *fHa0s, *fHa1s, *fHa2s;
  TH2F *fHa0CMs, *fHa1CMs, *fHa2CMs;
  TH2F *fHch2D, *fHchamp2D;
  TH1D *fHavgamp, *fHsigamp, *fHavgampCM, *fHsigampCM;
  TH2D *fHavgamp2D, *fHsigamp2D, *fHavgampCM2D, *fHsigampCM2D;
  TH1D *fHCMa2;
  TH2D *fHCMa2apv;
  TH2D *fHevtamp2D, *fHevtampcut2D;
  TH2F *fHtvsTCSph;
  TH2F *fHapvvsa2mean;
  TH1D *fHa2tmp;

  // name of the temporary histo used for a2 fit for each APV
  std::string fHa2tmpname;

  stat_s fStat[MAX_STAT_RiAPV];
  stat_s fStatCM[MAX_STAT_RiAPV];

  double tcsphase;

  /// variables
//   Variable *fVa2CM;

  /// common mode correction for each chip
  double fCMapv[MAX_NBAPV_RiAPV][3];

  // Cluster variables
//   double tcsphase;
//   Variable *fVcPos, *fVcA0, *fVcA1, *fVcA2, *fVcSize, *fVcTimeCorr;

  // Cluster histograms
//   TH1F *fHcHits, *fHcHitsT;
//   TH1F *fHcPos, *fHcPosT;
//   TH1F *fHcA0, *fHcA1, *fHcA2, *fHcA2T;
//   TH1F *fHcSize, *fHcSizeT;
//   TH2F *fHcAmpRatio, *fHcAmpRatioT;
//   TH1F *fHcTime1, *fHcTime1T;
//   TH1F *fHcTime2, *fHcTime2T;
//   TH2F *fHcTime1vs2, *fHcTime1vs2T;
//   TH1F *fHcTime, *fHcTimeT;
//   TH2F *fHcTimevsTCS;
//   TH1F *fHcTimeCorr, *fHcTimeCorrT;

//   int fNclustKept;

  /// additional variables
  Variable *fVxpx, *fVypx;

  /// current channel looked at (see PlaneRiAPVPanel)
  int fCurChan;

  /// says if the 20 not used channels are at the end of the APV channels (true) or at the both edges (false)
  bool fOldAPVBoard;

  /// says if the card is of S type (true) or J type (false)
  bool fCardRotated;


 public:

  /// calibration

//   int *flag;
//   float *pedestal;
//   float *sigma;
//
//   std::vector<float> calib_time;

#if USE_DATABASE == 1

  class APVchannel {
  public:
    int flag;
    float ped, sigma, calped, calsigma;
    APVchannel() : flag(1), ped(700), sigma(4.5), calped(700),
      calsigma(4.5) {}
    APVchannel(int f, float p, float s, float cp, float cs) {
      flag=f; ped=p; sigma=s; calped=cp; calsigma=cs; 
    }
  };
  class APVref {
  public:
    int ldcId, srcId, adcId, chipId;
    APVref(int l, int sr, int ad, int ci) : ldcId(l), srcId(sr), adcId(ad), chipId(ci) {}
    APVref(const APVref& ref) { ldcId=ref.ldcId; srcId=ref.srcId; adcId=ref.adcId; chipId=ref.chipId; }
    bool operator< (const APVref& ot) const {
      return (srcId*10000+adcId*100+chipId) < (ot.srcId*10000+ot.adcId*100+ot.chipId);
    }
  };
  class APVcalib {
  public:
    std::vector<APVchannel> channel;
    APVref ref;
    APVcalib() : ref(0,0,0,0) {}
    APVcalib(const APVref& ref) : ref(ref) { channel.clear(); }
    APVcalib(int l, int sr, int ad, int ci) : ref(l, sr, ad, ci) { channel.clear(); }
  };
  friend istream& operator>>(istream& in, std::map<const PlaneRiAPV::APVref, PlaneRiAPV::APVcalib> &c);
  std::map<const APVref,APVcalib> calib_data;

  /// arrays with calibration data
  calib_s *calib_arr;

#endif //USE_DATABASE

  PlaneRiAPV(const char *detname,int nchan, int center, int width);

#ifndef  __CINT__
  void EndEvent(const CS::DaqEvent &event);
#endif

  void Init(TTree* tree =0);

  /// Resets all histograms (and associated counters)
  void ResetHistograms();

//   void Clusterize();

#if USE_DATABASE == 1
  /// method to read calibration data
  void ReadCalib(const tm& t);
#endif

//   void TextOutput(ostream& out);
//   float Gain();

/// Calculates the common mode correction to apply to latch all data for each chip
void calcCommonModeCorrection(double cut);

/// Projects the a2 amplitude vs chan histogram to the a2 axis in a predefined channel
void ChannelA2Spectrum();

/// sets the channel number for the projection of the time vs chan histogram
void SetChannel(int channel) {fCurChan=channel;}

/// create the control panel
virtual void ControlPanel(const TGWindow* p, const TGWindow* main);

/// read the RiAPV Digits
#ifndef  __CINT__
void StoreDigit(CS::Chip::Digit* digit); 
#endif

/// write the text output file (containing result of the slope fit of a2)
void TextOutput(ostream& out);


private: 

/// StoreDigit(channel, amp,...) forbidden with PlaneRiAPV
void StoreDigit(int channel, int amp1, int amp2, int amp3) {  };

 friend class GroupRiAPV;


  ClassDef(PlaneRiAPV,0)
};

#endif










