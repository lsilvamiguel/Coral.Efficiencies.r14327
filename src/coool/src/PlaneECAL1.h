#ifndef __PlaneECAL1__
#define __PlaneECAL1__

#include "Plane.h"
#include "PlanePanel.h"
#include <vector>

#include "TThread.h"
//=========== ECAL1 parameters =========================
//       #define nb_ECAL1_SADCcnl  608
#define nb_ECAL1_SADCcnl  1056           // new GAMS central hole 20X12 -> 208 channels plus
#define nMAX_ECAL1_SADC_Samples 32
//=======================================================

using namespace std;

class PlaneECAL1 : public Plane {

  private:
/// maximum multiplicity per adc channel
  static const int fMAX_MULT;   

// number of ECAL1 rows/columns
  const int    fNrows, fNcols;

// number of channels:  total/FIADC/SADC
  int fNchan, nSADCchan;

// number of tree SADS samples per channels
  int fNsamples;

//  ampl./led cuts
  int SADCampcut, SADCLEDcut;

//  SADC Ped/LED windows for calib and physical events
    int pedstart, pedend;
    int pulsstart, pulsend;
    int calpedstart, calpedend;
    int ledstart, ledend;
    int nCounterX,nCounterY;
    int RefLedMin,RefLedMax;
    double wwwled;

// current channel looked at in PlaneECAL1Panel
  int fCurChan;

  /// data for each SADC channel
  int  fSADCRow[nb_ECAL1_SADCcnl], fSADCCol[nb_ECAL1_SADCcnl];
  int  fSADCev[nb_ECAL1_SADCcnl][nMAX_ECAL1_SADC_Samples];

  /// row for each digit          
  int  *fRow; //[fNchan*fMAX_MULT]    comment needed by rootcint -don't delete
  /// column for each digit          
  int  *fCol; //[fNchan*fMAX_MULT]

  /// pulse amplitude for each cell          
  double  *fPulseAmp; //[fNchan]

  /// amplitude for each digit        
  int  *fAmp; //[fNchan*fMAX_MULT]
  /// amplitude for each digit        
  int  *fAmpsum; //[fNchan*fMAX_MULT]

// variables
  Variable *fVampcut,*fVped,*fVpuls,*fVledped,*fVled;
  Variable *fVrow, *fVcol, *fVamp, *fVsum;
  Variable *fNcounter;

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
  TH2F *hxynoise, *hampnoise;
  TH1F *hrnd, *hledlen;
  TH1F *hLedTm, *hPhysTm;
  TH1F *hmnled, *hmnledrms, *hmnledamp, *hmnledr;
  TH1F *hmnped, *hmnpedr;
  TH1F *hprofx, *hprofy;
  TProfile *hLedProf,*hPedProf, *href;
  vector<TProfile*> hLasPuls;
  vector<TH1F*> hLasAmpl;

double SADCped(int ipbg, int  ipend, int *data,double &RMS);
double SADCpuls(int ipbg, int  ipend, int *data,double ped,
                                            double &sum,double &time);
void makeProfiles(TProfile *prof, TH1F *mean, TH1F *RMS, TH1F *integral,
                  TH1F *rmsint,TH1F *rmsint1);
void makeReference(TProfile *prof, TH1F *mean, TProfile *ref);

 public:
 /*! \brief constructor
    \param detname detector name
    \param nchan number of channels
    \param center center of the time band
    \param width width of the time band
  */
  PlaneECAL1(const char *detname,int ncols, int nrows, int center, int width);
// ~PlaneECAL1() {}
  virtual ~PlaneECAL1();

  /// Passes a digit to the plane
  void StoreDigit(int col, int row,  std::vector<float>& data);
  void TextOutput(ostream& out);

#ifndef __CINT__
   void StoreDigit(CS::Chip::Digit* digit);
   void EndEvent(const CS::DaqEvent &event);
#endif

  /// book histograms and branch the tree
  void Init(TTree* tree =0);

  /// sets the channel number for the projection of the amp vs chan histogram
  void SetChannel(int channel) {fCurChan=channel;}

  void AmpSpectrum();

  /// \return number of channels
  int GetNchannels() const {return fNrows*fNcols;}
  int GetNrows() const {return fNrows;}
  int GetNcols() const {return fNcols;}

  /// get pulse amplitude array
  const double* GetPulseAmp() const { return fPulseAmp; }

  /// get hit histograms
  TH1F* GetHhit() {return hSADChit;}

  /// get amplitude
  TH1F* GetAmp() {return fHa;}

  void ControlPanel(const TGWindow *p, const TGWindow *main);

  ClassDef(PlaneECAL1,1)
};

#endif



