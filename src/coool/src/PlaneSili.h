#ifndef __PlaneSili__
#define __PlaneSili__

#include "PlaneAPV.h"
#include "TThread.h"
#include "TRandom.h"


#define MAX_STAT 1536
#define MAX_APV_NBCH 128
#define MAX_NBAPV 12
#define MAX_AMP_DEVIATION 20.


/// Plane for detectors with a readout based on an Sili chip.

class PlaneSili : public PlaneAPV {

 private:
  typedef struct { int x, y; } xy_s ;
  typedef struct { double sum, sum2; int nb; } stat_s ;
  typedef struct { int flag; float pedestal, sigma; } calib_s ;

#ifndef __CINT__
  // Silicon time calibration reading
  class ChannelCalib {
  public:
    int ch;
    float t0;
    int flag;
    float ped, sigma, calped, calsigma;
    int ldcId, srcId, adcId, chipId;
    ChannelCalib() : ch(0),t0(0),flag(0) {}
    ChannelCalib(const char *s) {
      if(3 != sscanf(s,"%d%f%d",&ch,&t0,&flag)) {
	throw CS::Exception("PlaneSili::ChannelCalib : bad line \"%s\"",s);
	std::cerr<<"bad line, exception not caught !"<<std::endl;
      }
    }
  /// calibration for every channel
    std::vector<ChannelCalib> calib_data;  //???????
  }; 
    /// calibration for every channel
  std::vector<ChannelCalib> calib_data; //????????????

  friend istream& operator>>(istream& in, PlaneSili::ChannelCalib &c) {
    in>>c.ch;
    in>>c.t0;
    in>>c.flag;
    return in;
  }

  std::vector<float> calib_time;

#endif  // __CINT__


  int fMinChan;
  static xy_s fXY[];
  static bool fXYinit;

  /// histograms
  TH2F *fHchVsTis;
  TH2F *fHampRatio, *fHampSubRatio, *fHtimeRatio;
  TH1F *fHoccup, *fHampSubRap;;
  TH1F *fHchits, *fHcch, *fHca2, *fHcs;
  TH1F *fHhitpos;
  TH2F *fChAmp, *fHtiming;
  TH2F *fHAmpTCStime;
  TH1F *fHa0CM, *fHa1CM, *fHa2CM;
  TH2F *fHa0s, *fHa1s, *fHa2s;
  TH2F *fHa0CMs, *fHa1CMs, *fHa2CMs;
  TH2F *fHch2D, *fHchamp2D;
  TH1D *fHavgamp, *fHsigamp, *fHavgampCM, *fHsigampCM;
  TH2D *fHavgamp2D, *fHsigamp2D, *fHavgampCM2D, *fHsigampCM2D;
  TH2D *fHevtampCM2D, *fHevtampCMcut2D;
  TH2F *fHchvscmc2;



  /// cluster variables
  double tcsphase;

  bool timing_initialised;
  double fPar0[6];
  double fPar1[6];

  Variable *fVcch, *fVca2, *fVcs, *ft_cor;
  //!!
  TH2F *fHctiming;
  TH1F *fHctime1, *fHctime2;
  TH2F *fHctimeCorr;
  TH1F *fHctime;
  TH2F *fHctimeCorrClusterTCS;
  TH2F *fHcRatio0TCS;
  TH2F *fHcRatio1TCS;
  TH2F *fHcTiTi;
  TH1F *fHcRightTime;
  TH1F *fHcCalibTime;


  stat_s fStat[MAX_STAT];
  stat_s fStatCM[MAX_STAT];

  /// common mode correction for each chip
  double cCMC_fCMapv[MAX_NBAPV][3];

  xy_s *ch2xy(int chan); // return x, y
  static void initconndata(void); // initialize connector structure

  int fNclustKept;

 public:

  PlaneSili(const char *detname,int nchan, int center, int width)
    : PlaneAPV(detname,nchan,center,width), fMinChan(0) {INeedGeom();}
 
  /// arrays with calibration data
  calib_s *calib_arr;

#ifndef  __CINT__
  void EndEvent(const CS::DaqEvent &event);
#endif
  void Init(TTree* tree =0);

  void Clusterize();

  double GetTCSphase(){return tcsphase;};

  /// Calculates the common mode correction to apply to latch all data for each chip
  void calcCommonModeCorrection(double cut);
  

#if USE_DATABASE == 1
  /// method to read calibration data
  void ReadCalib(const tm& t);

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
    std::vector<ChannelCalib> channel;
  };


#endif

  ClassDef(PlaneSili,0)
};

#endif


