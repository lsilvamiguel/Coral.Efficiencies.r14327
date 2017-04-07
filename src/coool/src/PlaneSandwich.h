#ifndef __PlaneSandwich__
#define __PlaneSandwich__

#include <vector>
#include "Plane.h"

#ifndef __CINT__
#include "ChipF1.h"
#include "ChipSADC.h"
using namespace CS;
#endif

using namespace std;

class PlaneSandwich : public Plane
{
  public:  

                        PlaneSandwich           (const char *detname);
                       ~PlaneSandwich           (void);

    void                Init                    (TTree* tree = 0);

    #ifndef __CINT__
    void                Reset                   (void);

    /// create the control panel
    virtual void        ControlPanel            (const TGWindow* p, const TGWindow* main);
    int                 GetNchannels            () { return 12; }

    void                StoreDigit              (CS::Chip::Digit* digit);
    void                EndEvent                (const CS::DaqEvent &event);
    #endif

  private :

    #ifndef __CINT__
    struct ChannelInfo {
      CS::uint16 min, max;
      float offs, rms, integral;
      bool pileup;
      const ChipSADC::Digit* digit_sadc;
      vector<const ChipF1::Digit*> digits_f1;

      ChannelInfo() : min(0), max(0), offs(0), rms(0), integral(0),
           pileup(false), digit_sadc(0) {}
    };


    void                put_in_hist             (TH2F* h2, bool subtract_offs,
                                                 const ChannelInfo& ci);
    void                StoreDigitSADC          (const CS::ChipSADC::Digit*);
    void                StoreDigitF1            (const CS::ChipF1::Digit*);
    #endif

    Variable *fT; // time_range field, used for histograms that
                  // have a time axis.

    static const unsigned int           channels_total = 12;

    #ifndef __CINT__
    map<int, ChannelInfo>               channel_infos;

    vector<const ChipSADC::Digit*>      digits_sadc;
    vector<const ChipF1::Digit*>        digits_f1;
    #endif

    // Total ADC entries
    TH1F *                          h1_adc_chan;
    // Number of ADC entries in event
    TH1F *                          h1_adc_multiplicity;
    // Total TDC hits per channesl
    TH1F *                          h1_tdc_chan;
    // Channel has TDC hits in how many events
    TH1F *                          h1_tdc_chan_hit;
    // TDC times
    TH1F *                          h1_tdc_time;
    // Channel vs TDC times
    TH2F *                              h2_tdc_chan_time;
    // Times for events with 12xADC
    TH2F *                              h2_tdc_chan_time_12;
    // Time for first entry
    TH2F *                              h2_tdc_chan_time_first;
    // Time for entries in histograms with echos
    TH2F *                              h2_tdc_chan_time_echo;

    // #TDC vs #ADC 
    TH2F *                              h2_tdc_vs_adc;
    TProfile *                          p_tdc_vs_adc;
    string                              name_tdc_vs_adc;
    
    // Scope picture, amplitudes
    vector<TH1F*>                   vh1_adc_ampl_max;
    vector<TH2F*>                       vh2_adc_samples;
    
    // Same, but with pileup and offset removed.
    vector<TH1F*>                   vh1_adc_ampl_max_offs;
    vector<TH2F*>                       vh2_adc_samples_offs;

    // Integrals of signals
    vector<TH1F*>                   vh1_adc_integrals;

    // ADC spectrum for events with < 12 ADC entries
    vector<TH2F*>                       vh2_adc_samples_offs_not_12;
    // Amplitude vs time for channels with 1 TDC entry & ADC data
    vector<TH2F*>                       vh2_ampl_vs_time;
    // Amplitude vs time for channels with >1 TDC entries & ADC data
    vector<TH2F*>                       vh2_ampl_vs_time_2;

    // ADC spectrum for events with MIP in #7
    vector<TH2F*>                       vh2_adc_mip_company;

    // Channels with hits in both adc and tdc ...
    TH1F *                          h1_adc_and_tdc;
    // ... channels with hits in either ...
    TH1F *                          h1_adc_or_tdc;
    // ... used to calculate this histogram, which is the quotient
    // and which is displayed in coool.
    TH1F *                          h1_adc_tdc_comparison;

    // Time difference between successive TDC hits in the same event
    vector<TH1F *>                  vh1_tdc_time_diff;

    #ifndef __CINT__
    struct earlier
      : public binary_function<const CS::ChipF1::Digit*,
                               const CS::ChipF1::Digit*, bool>
    {
      bool operator()(const CS::ChipF1::Digit* x, const CS::ChipF1::Digit* y)
      {
        return x->GetTimeDecoded() < y->GetTimeDecoded();
      }
    };
    #endif



    ClassDef(PlaneSandwich,1)
};

#endif
