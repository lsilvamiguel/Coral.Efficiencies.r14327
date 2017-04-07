#ifndef __PlaneECAL2__
#define __PlaneECAL2__

#include "Plane.h"
#include "PlanePanel.h"
#include "TProfile2D.h"
#include "TThread.h"
//=========== ECAL2 parameters =========================
#define nb_ECAL2_SADCcnl  48*64
#define nMAX_ECAL2_SADC_Samples 32
//=======================================================

class PlaneECAL2 : public Plane {
  
  private:
// number of ECAL2 rows/columns
  const int    fNrows, fNcols;

// number of channels:  total/FIADC/SADC
  int fNchan, nSADCchan;

// number of tree SADS samples per channels
  int fNsamples;

//  SADC Ped/LED windows for calib and physical events 
    int pedstart, pedend;
    int pulsstart, pulsend;
    int calpedstart, calpedend;
    int ledstart, ledend;

// current channel looked at in PlaneECAL2Panel
  int fCurChan;

  int fledratecounter;
  /// data for each SADC channel        
  int  fSADCRow[nb_ECAL2_SADCcnl], fSADCCol[nb_ECAL2_SADCcnl];
  int  fSADCev[nb_ECAL2_SADCcnl][nMAX_ECAL2_SADC_Samples];

//   FI/SADC event hit numbers
 int  fNSADChits;

 int SADCampcut, SADCLEDcut;  // SADC  amplitude/LED  cuts
 int FIampcut, FILEDcut;      // FIADC     amplitude/LED  cuts
 double wwwled_perif, wwwled_cent;
 //time offset as read from the calibrations file
 float fT0;
 bool fUseCalib;
 float fTimeWindow;

  /// variables
  Variable *fVrow, *fVcol, *fVamp;
  Variable *fVped, *fVpuls, *fVledped,*fVled;
  Variable *fVscuts, *fVfcuts;

  // markers for histograms
  TH1F *hSADChit,*hSADCsmpl;
  TH2F *fHrca, *fHrc1, *hnhitcut;
  TH2F *fHavsadr, *fHxyac;
  TH1F *fHa, *fHch;
  TH1F *hPed, *hPedRMS;
  TH2F *hPedXY;
  TH2F *hLed, *hLedChnl,*hrefXY;
  TH2F *hoverflow, *fHchac;
  TH1F *hSADCsum, *hSADCcutsum;
  TH1F *hSADCcuthits;
  TH2F *hxynoise, *hampnoise,*hxynoisecut, *hampnoisecut;
  TH1F *hrnd, *hledlen;
  TH1F *hLedTm, *hPhysTm;
  TH1F *hmnled, *hmnledrms, *hmnledamp, *hmnledr;
  TH1F *hmnped, *hmnpedr;
  TH1F *hprofx, *hprofy;
  TH2F *hTimech;
  TProfile2D *hTimechprof;
  TH2F *hLEDTimech;

  
//  std::vector<TH1F*> hTime;
  TProfile *hLedProf,*hPedProf, *href;

  std::vector<float> fTcalib;
  std::vector<Float_t>* ledcalibCFD;
  template <class T>
  void ReadFromDataBase(T& v, const tm& t,const char* typecalib) ;
  class TimeCalib {
    public:
      int x,y;
      float time, sigma,stat; 
      TimeCalib() :x(0),y(0),time(0),sigma(0),stat(0) {}
      TimeCalib(const char *s) {
        if(5 != sscanf(s,"%d%d%f%f%f",&x,&y,&time,&sigma,&stat)) {
          throw CS::Exception("PlaneECAL2::ChannelCalib : bad line \"%s\"",s);
          std::cerr<<"bad line, exception not caught !"<<std::endl;
        }
      }
      void Print() {
        printf("X:  %2d  Y:  %2d  time:  %3.2f  sigma:  %2.1f  stat:  %9.0f\n",x,y,time,sigma,stat);
      }
  };

  friend istream& operator>>(istream& in, PlaneECAL2::TimeCalib &c)
  {
    in>>c.x;
    in>>c.y;
    in>>c.time;
    in>>c.sigma;
    in>>c.stat;
    return in;
  }





  
double CFDtime(int ipbg, int ipend, int *data);
double Hmaxtime(int ipbg, int ipend, int *data);
double SADCped(int ipbg, int  ipend, int *data,double &RMS);
double SADCpuls(int ipbg, int  ipend, int *data,double ped,
                                            double &sum,double &time);
void makeProfiles(TProfile *prof, TH1F *mean, TH1F *RMS, TH1F *integral,
                  TH1F *rmsint, TH1F *rmsint1);
void makeReference(TProfile *prof, TH1F *mean, TProfile *ref);

 public:
 /*! \brief constructor
    \param detname detector name
    \param nchan number of channels
    \param center center of the time band
    \param width width of the time band
  */
  PlaneECAL2(const char *detname,int ncols, int nrows, int center, int width);
 ~PlaneECAL2() {}
  
  /// Passes a digit to the plane
  void StoreDigit(int col, int row,  std::vector<float>& data);
  void StoreDigit(int col, int row,  int amp);
  void TextOutput(ostream& out);

#ifndef __CINT__
   void StoreDigit(CS::Chip::Digit* digit);
   void EndEvent(const CS::DaqEvent &event);  
#endif

  void readledcalibCFD();
  /// book histograms and branch the tree
  void Init(TTree* tree =0);

  /// sets the channel number for the projection of the amp vs chan histogram
  void SetChannel(int channel) {fCurChan=channel;}  

  void AmpSpectrum();

  /// \return number of channels
  int GetNchannels() {return fNrows*fNcols;}
  int GetNrows() {return fNrows;}
  int GetNcols() {return fNcols;}

  /// get hit histograms
  TH1F* GetHhit() {return hSADChit;}

  /// get amplitude
  TH1F* GetAmp() {return fHa;}
  
  void ControlPanel(const TGWindow *p, const TGWindow *main);
  void ReadCalib(const tm&);

  ClassDef(PlaneECAL2,1)
};

#endif



