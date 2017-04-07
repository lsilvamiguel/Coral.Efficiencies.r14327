// assimilated from PlaneSandwich
// modified by jasinski@kph.uni-mainz.de 22.06.09

#ifndef __PlaneCEDAR__
#define __PlaneCEDAR__

#include <vector>
#include "Plane.h"

#ifndef __CINT__
#include "ChipF1.h"
#include "ChipSADC.h"
#include "ChipGandalf.h"
#include "PlaneGandalf.h"
using namespace CS;
#endif

using namespace std;

class PlaneCEDAR : public Plane
{
public:

	PlaneCEDAR (const char *detname);
	~PlaneCEDAR (void);

	void Init (TTree* tree = 0);

#ifndef __CINT__
	void Reset (void);

	/// create the control panel
	//virtual void        ControlPanel            (const TGWindow* p, const TGWindow* main);
	int GetNchannels () {return 8;}

	void StoreDigit (CS::Chip::Digit* digit);
	void EndEvent (const CS::DaqEvent &event);
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
	struct g_ChannelInfo {
		CS::uint16 min, max;
		float offs, rms, integral;
		bool pileup;
		vector<const ChipGandalf::DigitGADC*> digits;

		g_ChannelInfo() : min(0), max(0), offs(0), rms(0), integral(0),
		pileup(false) {}
	};
	bool hit_tdc_cut[8];
	bool hit_g_cut[8];
	int  multiplicity_counter[8];
	float adc_amplitudes[8];

	void put_in_hist (TH2F* h2, bool subtract_offs,
			const ChannelInfo& ci);
	void StoreDigitSADC (const CS::ChipSADC::Digit*);
	void StoreDigitF1 (const CS::ChipF1::Digit*);
	void StoreDigitGandalf (const CS::ChipGandalf::DigitGADC*);
#endif

	Variable *fT[2]; // time_range field, used for histograms that
	// have a time axis. We have two CEDARs

	int numb_cedar; // is it CEDAR 1 or 2?

	static const unsigned int channels_total = 8;

#ifndef __CINT__
	map<int, ChannelInfo> channel_infos;
	map<int, g_ChannelInfo> g_channel_infos;
	map<int, PlaneGandalf*> g_planes;
	vector<const ChipSADC::Digit*> digits_sadc;
	vector<const ChipF1::Digit*> digits_f1;
	vector<const ChipGandalf::DigitGADC*> digits_gandalf;
#endif

	// Total ADC entries
	TH1F * h1_adc_chan;
	// Number of ADC entries in event
	TH1F * h1_adc_multiplicity;
	// Total TDC hits per channesl
	TH1F * h1_tdc_chan;
	// Channel has TDC hits in how many events
	TH1F * h1_tdc_chan_hit;
	// TDC times
	TH1F * h1_tdc_time;
	// Channel vs TDC times
	TH2F * h2_tdc_chan_time;
	TH2F * h2_tdc_chan_time_all;
	TH2F * h2_maj6_tdc_alignment; // to show the alignment of the CEDAR
	// Times for events with 12xADC
	//TH2F *                              h2_tdc_chan_time_12;
	// Time for first entry
	//TH2F *                              h2_tdc_chan_time_first;
	// Time for entries in histograms with echos
	//TH2F *                              h2_tdc_chan_time_echo;

	// #TDC vs #ADC
	TH2F * h2_tdc_vs_adc;
	TProfile * p_tdc_vs_adc;
	string name_tdc_vs_adc;

	// Scope picture, amplitudes
	vector<TH1F*> vh1_adc_ampl_max;
	vector<TH2F*> vh2_adc_samples;

	// Same, but with pileup and offset removed.
	vector<TH1F*> vh1_adc_ampl_max_offs;
	vector<TH2F*> vh2_adc_samples_offs;

	// adc max ampl spectrum cut by tdc response
	vector<TH1F*> vh1_adc_ampl_max_offs_when_tdc_cut;
	vector<TH1F*> vh1_adc_ampl_max_offs_when_tdc;

	// recording majorities
	TH2F* h2_tdc_maj_vs_adc_maj;
	TH1F* h1_tdc_maj_cut;
	TH1F* h1_tdc_maj_6_chan;
	//TH1F* h1_tdc_mean_photons_7_over_8_cut;
	//TH1F* h1_tdc_mean_photons_6_over_8_cut;

	// Integrals of signals
	//vector<TH1F*>                   vh1_adc_integrals;

	// ADC spectrum for events with < 12 ADC entries
	//vector<TH2F*>                       vh2_adc_samples_offs_not_12;
	// Amplitude vs time for channels with 1 TDC entry & ADC data
	//vector<TH2F*>                       vh2_ampl_vs_time;
	// Amplitude vs time for channels with >1 TDC entries & ADC data
	//vector<TH2F*>                       vh2_ampl_vs_time_2;

	// ADC spectrum for events with MIP in #7
	//vector<TH2F*> vh2_adc_mip_company;

	// Channels with hits in both adc and tdc ...
	TH1F * h1_adc_and_tdc;
	// ... channels with hits in either ...
	TH1F * h1_adc_or_tdc;
	// ... used to calculate this histogram, which is the quotient
	// and which is displayed in coool.
	TH1F * h1_adc_tdc_comparison;

	// Time difference between successive TDC hits in the same event
	//vector<TH1F *>                  vh1_tdc_time_diff;



	//////////////////////////////////////////////////////////////////////////////////////////////////
	///////////////////////////////////////  Plots for Gandalf readout


	// Total hits per channel
	TH1F * h1_ghits_chan;
	// Total hits per channel
	TH2F * h2_g_mult_chan;

	// times
	TH1F * h1_g_time;
	// Channel vs amplitude
	TH2F * h2_g_chan_ampl;

        // Channel vs integral
	TH2F * h2_g_chan_int;
        
	TH2F** vh2_g_amp_f1time;
	TH2F** vh2_g_time_f1time;

	// Channel vs times
	TH2D * h2_g_chan_time_all;

	// Number of entries in event
	TH1F * h1_g_multiplicity;


	// Channels with hits in both adc and tdc ...
	TH1F * h1_gandalf_and_tdc;
	// ... channels with hits in either ...
	TH1F * h1_gandalf_or_tdc;
	// ... used to calculate this histogram, which is the quotient
	// and which is displayed in coool.
	TH1F * h1_gandalf_tdc_comparison;

	// Channels with hits in both adc and tdc ...
	TH1F * h1_gandalf_at_and_tdc;
	// ... channels with hits in either ...
	TH1F * h1_gandalf_at_or_tdc;
	// ... used to calculate this histogram, which is the quotient
	// and which is displayed in coool.
	TH1F * h1_gandalf_at_tdc_comparison;
	TH2I *h2DataTransmission;

	float g_amplitudes[channels_total];

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
	struct earlierG
	: public binary_function<const CS::ChipGandalf::DigitGADC*,
	const CS::ChipGandalf::DigitGADC*, bool>
	{
		bool operator()(const CS::ChipGandalf::DigitGADC* x, const CS::ChipGandalf::DigitGADC* y)
		{
			if (x->GetTimeDecoded() && y->GetTimeDecoded())
				return x->GetTimeDecoded() < y->GetTimeDecoded();
			//std::cout << "Undecoded digit: " << x->GetDetID() << "  "<< x->getChannel()<< "  " << x->GetDataID() << std::endl;
			return false;
		}
	};

#endif

	ClassDef(PlaneCEDAR,1)
};

#endif
