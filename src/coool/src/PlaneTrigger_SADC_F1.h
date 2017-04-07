#ifndef __PlaneTrigger_SADC_F1__
#define __PlaneTrigger_SADC_F1__

#include "Plane.h"
#include "PlanePanel.h"
#include <vector>

#ifndef __CINT__
#include "ChipF1.h"
#include "ChipSADC.h"
using namespace CS;
#endif

using namespace std;

class PlaneTrigger_SADC_F1 : public Plane
{

 public:

  /*! \brief constructor
      \param detname detector name
      \param nchan number of channels
      \param center center of the time band
      \param width width of the time band
   */
  PlaneTrigger_SADC_F1(const char *detname, int nAdcChan, int nTdcChan, int center, int width);
  ~PlaneTrigger_SADC_F1(void);

  virtual void Init(TTree* tree = 0);

#ifndef __CINT__
  virtual void Reset(void);
  virtual void StoreDigit(CS::Chip::Digit* digit);
  virtual void EndEvent(const CS::DaqEvent &event);
  virtual void ControlPanel(const TGWindow* p, const TGWindow* main);
#endif

 private:

#ifndef __CINT__
  struct ChannelInfo {
    CS::uint16 min, max;
    float offs, rms, integral;
    bool pileup;
    const ChipSADC::Digit* digit_sadc;
    //vector<const ChipF1::Digit*> digits_f1;

    ChannelInfo() : min(0), max(0), offs(0), rms(0), integral(0), pileup(false), digit_sadc(0) {}

    void clear()
    {
      //digits_f1.clear();
      //digits_f1.reserve(10);
    }
  };

  void put_in_hist(TH2F* h2,bool subtract_offset, const ChannelInfo& ci);
  void StoreDigitSADC (const CS::ChipSADC::Digit*);
  void StoreDigitF1 (const CS::ChipF1::Digit*);
#endif

#ifndef __CINT__
  map<int, ChannelInfo> channel_infos;
  vector<const ChipSADC::Digit*> digits_sadc;
  vector<const ChipF1::Digit*> digits_f1;
#endif
  
  static const float fF1_TICK;

  /// maximum multiplicity per tdc channel
  static const int fMAX_MULT;

  /// frequency of updating rate histograms
  static const int fRATE_UPDATE;

  int fNAdcChan;         // number of ADC channels
  int fNTdcChan;         // number of TDC channels
  float fCenter, fWidth; // definition of histogram size

  bool inTime;
  int tdc_mult;

  Variable* fVch;
  Variable* fVt;
  Variable* fVtNs;
  Variable* fVtOnTrig;
  Variable* fVtOffTrig;

  float* max_sadc;  //[fNAdcChan]

  //histograms
  vector<TH2F*> h2_ADC_samples;
  vector<TH2F*> h2_ADC_samples_tdc_cor;
  vector<TH2F*> h2_ADC_samples_wo_pileup;
  vector<TH2F*> h2_ADC_samples_wo_pileup_tdc_cor;

  vector<TH1F*> h1_ADC_samples_max;
  vector<TH1F*> h1_ADC_samples_max_tdc_cor;
  vector<TH1F*> h1_ADC_samples_max_wo_pileup;
  vector<TH1F*> h1_ADC_samples_max_wo_pileup_tdc_cor;
  
  vector<TH2F*> h2_ADC_samples_off;
  vector<TH2F*> h2_ADC_samples_off_tdc_cor;
  vector<TH2F*> h2_ADC_samples_off_wo_pileup;
  vector<TH2F*> h2_ADC_samples_off_wo_pileup_tdc_cor;

  vector<TH1F*> h1_ADC_samples_max_off;
  vector<TH1F*> h1_ADC_samples_max_off_tdc_cor;
  vector<TH1F*> h1_ADC_samples_max_off_wo_pileup;
  vector<TH1F*> h1_ADC_samples_max_off_wo_pileup_tdc_cor;

  vector<TH1F*> h1_ADC_samples_integral;
  vector<TH1F*> h1_ADC_samples_integral_tdc_cor;
  vector<TH1F*> h1_ADC_samples_integral_wo_pileup;
  vector<TH1F*> h1_ADC_samples_integral_wo_pileup_tdc_cor;

  vector<TH1F*> h1_ADC_samples_diff;

  TH1F *h1_adc_ch;

  TH1F* h1_tdc_time;
  TH1F* h1_tdc_time_ns;

  TH1F *h1_tdc_ch;

  TH1F *h1_tdc_mult;


  TH1F *h1_Analog_Sum_max;

  TH2F *h2_Analog_Sum;

  TH2F *h2_tdc_t_vs_ch;
  TH1F* h1_tdc_time_on_trig;
  TH1F* h1_tdc_time_off_trig;
  TH1F* h1_tdc_ch_on_trig;
  TH1F* h1_tdc_ch_off_trig;
  TH1F* h1_tdc_cor_ch_on_trig;
  TH1F* h1_tdc_rates;

  ClassDef(PlaneTrigger_SADC_F1,1)

};

#endif
