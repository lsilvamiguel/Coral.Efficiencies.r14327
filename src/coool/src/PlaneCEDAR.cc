// assimilated from PlaneSandwich
// modified by jasinski@kph.uni-mainz.de 22.06.09
// modified by tosello@to.infn.it 2015 July 26

#ifndef ROOT_TROOT
#include "TROOT.h"
#endif

#include "PlaneCEDAR.h"

//#include "PlaneCEDARPanel.h"

#include <algorithm>

ClassImp( PlaneCEDAR);

// unsigned int PlaneCEDAR::channels_total = 12;

int GANDALF_THRESH = 250;


PlaneCEDAR::PlaneCEDAR(const char *detname) :
  Plane(detname) {
  // single precision mode until 28.10.09
  //printf("PlaneCEDAR::PlaneCEDAR(): %s\n",detname);
  //fT[0] = AddVariable("time_range_CEDAR1", 2448 - 2430, -2448, -2430, 2448
  //		- 2430);
  //fT[1] = AddVariable("time_range_CEDAR2", 2466 - 2448, -2466, -2448, 2466
  //		- 2448);
  // double precision mode since 28.10.09 after trigger timeshift

  //fT[0] = AddVariable("time_range_CEDAR1", 1730 - 1680, 1680, 1730, 1730 - 1680);
  //fT[1] = AddVariable("time_range_CEDAR2", 1700 - 1650, 1650, 1700, 1700 - 1650);
  // new trigger delay in 2012
  // fT[0] = AddVariable("time_range_CEDAR1", 1460 - 1440, -1460, -1440, 1460 - 1440);
  // fT[1] = AddVariable("time_range_CEDAR2", 1475 - 1455, -1475, -1455, 1475 - 1455);
  //
  // new trigger delay in 2015 
  fT[0] = AddVariable("time_range_CEDAR1", 2000 - 1400,  1400,  2000, 2000 - 1400);
  fT[1] = AddVariable("time_range_CEDAR2", 2000 - 1400,  1400,  2000, 2000 - 1400);

  for (int chan = 0; chan < GetNchannels(); chan++) {
    adc_amplitudes[chan] = 0.;
    hit_tdc_cut[chan] = false;
  }
}

PlaneCEDAR::~PlaneCEDAR(void) {
  printf("PlaneCEDAR::~PlaneCEDAR()\n");
  delete h1_adc_or_tdc;
  delete h1_adc_and_tdc;

  delete h2_tdc_vs_adc;
}

void PlaneCEDAR::Init(TTree* tree) {
  // 	cout << " initialization: " << endl;
  const char *channel_names[channels_total] = { "PM1", "PM2", "PM3", "PM4",
						"PM5", "PM6", "PM7", "PM8" };
  bool found(false);
  std::string name, title;
  numb_cedar = 0;

  float t_min[2] = { fT[0]->GetMin(), fT[1]->GetMin() };
  float t_max[2] = { fT[0]->GetMax(), fT[1]->GetMax() };
  int t_nbins[2] = { fT[0]->GetNbins(), fT[1]->GetNbins() };

  if (strncmp(fName.c_str(), "CE01", 4) == 0) {
    //		cout << " booking histograms for CEDAR 1 " << endl;
    numb_cedar = 0;
    found = true;
  }
  if (strncmp(fName.c_str(), "CE02", 4) == 0) {
    //else {
    //		cout << " booking histograms for CEDAR 2 " << endl;
    numb_cedar = 1;
    found = true;
  }
  if (!found) {
    std::cout << " error in PlaneCEDAR::Init(): detector name " << fName
	      << " not known !" << endl;
  }
  name = fName + "_adc_chan";
  title = fName + " SADC no of entries";
  h1_adc_chan = new TH1F_Ref(name.c_str(), title.c_str(), channels_total, 1,
			     channels_total + 1, fRateCounter);
  h1_adc_chan->SetMinimum(0);
  //AddHistogram(h1_adc_chan);
  if (fReferenceDirectory)
    ((TH1F_Ref*) h1_adc_chan)->SetReference(fReferenceDirectory);

  name = fName + "_adc_multiplicity";
  title = fName + " SADC multiplicity distribution";
  h1_adc_multiplicity = new TH1F_Ref(name.c_str(), title.c_str(),
				     channels_total, 1, channels_total + 1, fRateCounter);
  h1_adc_multiplicity->SetMinimum(0);
  //AddHistogram(h1_adc_multiplicity);
  if (fReferenceDirectory)
    ((TH1F_Ref*) h1_adc_multiplicity)->SetReference(fReferenceDirectory);


  ////////////////////////////////// Gandalf Plots
  ////////////////////////////////////////////////////////////////

  name = fName + "_g_hits_chan";
  title = fName + " Gandalf no of entries";
  h1_ghits_chan = new TH1F_Ref(name.c_str(), title.c_str(), channels_total, 1,
			       channels_total + 1, fRateCounter);
  h1_ghits_chan->SetMinimum(0);
  AddHistogram(h1_ghits_chan);
  if (fReferenceDirectory)
    ((TH1F_Ref*) h1_ghits_chan)->SetReference(fReferenceDirectory);

  name = fName + "_g_multiplicity";
  title = fName + " Gandalf hit-multiplicity distribution";
  h1_g_multiplicity = new TH1F_Ref(name.c_str(), title.c_str(),
				   50, 0, 50, fRateCounter);
  //			channels_total, 1, channels_total + 1, fRateCounter);
  h1_g_multiplicity->SetMinimum(0);
  AddHistogram(h1_g_multiplicity);
  if (fReferenceDirectory)
    ((TH1F_Ref*) h1_g_multiplicity)->SetReference(fReferenceDirectory);

  name = fName + "_g_multiplicity_vs_channel";
  title = fName + " Gandalf hit-multiplicity per channel";
  h2_g_mult_chan = new TH2F(name.c_str(), title.c_str(), channels_total, 1,
			    channels_total + 1, 10, 0, 10);
  //			channels_total + 1, 5, 0, 5);
  h2_g_mult_chan->SetOption("col");
  AddHistogram(h2_g_mult_chan);

  name = fName + "_g_time";
  title = fName + " Gandalf time distribution";
  h1_g_time = new TH1F_Ref(name.c_str(), title.c_str(),
			   1000,0,0,fRateCounter);
  AddHistogram(h1_g_time);
  if (fReferenceDirectory)
    ((TH1F_Ref*) h1_g_time)->SetReference(fReferenceDirectory);

  name = fName + "_g_chan_ampl";
  title = fName + " Gandalf amplitude distribution ";
  h2_g_chan_ampl = new TH2F(name.c_str(), title.c_str(), channels_total, 1,
			    channels_total + 1, 1000, 0, 4069);
  h2_g_chan_ampl->SetOption("colz");
  AddHistogram(h2_g_chan_ampl);

  name = fName + "_g_chan_integral";
  title = fName + " Gandalf integral distribution ";
  h2_g_chan_int = new TH2F(name.c_str(), title.c_str(), channels_total, 1,
			   channels_total + 1, 1000, 0, 100000);
  h2_g_chan_int->SetOption("colz");
  AddHistogram(h2_g_chan_int);
        
  name = fName + "_g_chan_time_all";
  title = fName + " Gandalf time distribution for channels; ; GetTimeDecoded() [ns]";
  h2_g_chan_time_all = new TH2D(name.c_str(), title.c_str(),
				channels_total, 1, channels_total + 1, 
				10000, -300.0e9, 300.0e9);
  h2_g_chan_time_all->SetOption("colz");
  AddHistogram(h2_g_chan_time_all);

  name = fName + "_SendDataGanF1VsChannel";
  title = fName + " Send Data from Gan/F1  vs channel";
  h2DataTransmission=new TH2I(name.c_str(),title.c_str(),channels_total,0,channels_total, 4, 0, 4);
  h2DataTransmission->SetOption("text");
  h2DataTransmission->GetXaxis()->SetTitle("Channel");
  h2DataTransmission->GetYaxis()->SetTitle("Hit Difference");
  h2DataTransmission->GetYaxis()->SetBinLabel(1, "g&f");
  h2DataTransmission->GetYaxis()->SetBinLabel(2, "g&!f");
  h2DataTransmission->GetYaxis()->SetBinLabel(3, "!g&f");
  h2DataTransmission->GetYaxis()->SetBinLabel(4, "!g&!f");
  AddHistogram(h2DataTransmission);

  name = fName + "_Gandalf_and_TDC";
  title = fName + " Gandalf && TDC";
  h1_gandalf_and_tdc = new TH1F(name.c_str(), title.c_str(), channels_total, 1,
				channels_total + 1);

  name = fName + "_Gandalf_or_TDC";
  title = fName + " Gandalf || TDC";
  h1_gandalf_or_tdc = new TH1F(name.c_str(), title.c_str(), channels_total, 1,
			       channels_total + 1);

  name = fName + "_g_tdc_comparison";
  title = fName + " (Gandalf && tdc) / (Gandalf || tdc)";
  h1_gandalf_tdc_comparison = new TH1F_Ref(name.c_str(), title.c_str(),
					   channels_total, 1, channels_total + 1, fRateCounter, 1);
  h1_gandalf_tdc_comparison->SetMaximum(1);
  h1_gandalf_tdc_comparison->SetMinimum(0);
  AddHistogram(h1_gandalf_tdc_comparison);


  name = fName + "_Gandalf_at_and_TDC";
  title = fName + " Gandalf above Thresh && TDC";
  h1_gandalf_at_and_tdc = new TH1F(name.c_str(), title.c_str(), channels_total, 1,
				   channels_total + 1);

  name = fName + "_Gandalf_at_or_TDC";
  title = fName + " Gandalf above Thresh || TDC";
  h1_gandalf_at_or_tdc = new TH1F(name.c_str(), title.c_str(), channels_total, 1,
				  channels_total + 1);

  name = fName + "_g_at_tdc_comparison";
  title = fName + " (Gandalf above Thresh && tdc) / (Gandalf above Thresh || tdc)";
  h1_gandalf_at_tdc_comparison = new TH1F_Ref(name.c_str(), title.c_str(),
					      channels_total, 1, channels_total + 1, fRateCounter, 1);
  h1_gandalf_at_tdc_comparison->SetMaximum(1);
  h1_gandalf_at_tdc_comparison->SetMinimum(0);
  AddHistogram(h1_gandalf_at_tdc_comparison);

  //////////////////////////////////////////////////////////////////////


  name = fName + "_tdc_chan";
  title = fName + " TDC channel all hits";
  h1_tdc_chan = new TH1F_Ref(name.c_str(), title.c_str(), channels_total, 1,
			     channels_total + 1, fRateCounter);
  h1_tdc_chan->SetMinimum(0);
  AddHistogram(h1_tdc_chan);
  if (fReferenceDirectory)
    ((TH1F_Ref*) h1_tdc_chan)->SetReference(fReferenceDirectory);

  name = fName + "_tdc_chan_cut";
  title = fName + " TDC channel with time cut";
  h1_tdc_chan_hit = new TH1F_Ref(name.c_str(), title.c_str(), channels_total,
				 1, channels_total + 1, fRateCounter);
  h1_tdc_chan_hit->SetMinimum(0);
  AddHistogram(h1_tdc_chan_hit);
  if (fReferenceDirectory)
    ((TH1F_Ref*) h1_tdc_chan_hit)->SetReference(fReferenceDirectory);

  name = fName + "_tdc_time";
  title = fName + " time distribution ; GetTimeDecoded() [ns]";
  h1_tdc_time = new TH1F_Ref(name.c_str(), title.c_str(),
			     t_nbins[numb_cedar], t_min[numb_cedar], t_max[numb_cedar],
			     fRateCounter);
  AddHistogram(h1_tdc_time);
  if (fReferenceDirectory)
    ((TH1F_Ref*) h1_tdc_time)->SetReference(fReferenceDirectory);

  name = fName + "_tdc_chan_time";
  title = fName + " time distribution for channels within time cut";
  h2_tdc_chan_time = new TH2F(name.c_str(), title.c_str(), channels_total, 1,
			      channels_total + 1, t_nbins[numb_cedar], t_min[numb_cedar],
			      t_max[numb_cedar]);
  h2_tdc_chan_time->SetOption("col");
  AddHistogram(h2_tdc_chan_time);

  name = fName + "_tdc_chan_time_all";
  title = fName + " time distribution for channels";
  h2_tdc_chan_time_all = new TH2F(name.c_str(), title.c_str(),
				  channels_total, 1, channels_total + 1, 5000, -3000, 2000);
  h2_tdc_chan_time_all->SetOption("col");
  AddHistogram(h2_tdc_chan_time_all);

  name = fName + "_tdc_maj_cut";
  title = fName + " majority distribution (0...8)";
  h1_tdc_maj_cut = new TH1F(name.c_str(), title.c_str(), 
			    channels_total+1, 0, channels_total + 1);
  // channels_total, 1, channels_total + 1);
  AddHistogram(h1_tdc_maj_cut);
  //if (fReferenceDirectory) ((TH1F_Ref*)h1_tdc_time)->SetReference(fReferenceDirectory);

  name = fName + "_tdc_maj_6_chan";
  title = fName + " channels fired for majority 6";
  h1_tdc_maj_6_chan = new TH1F(name.c_str(), title.c_str(), channels_total, 1,
			       channels_total + 1);
  AddHistogram(h1_tdc_maj_6_chan);

  name = fName + "_tdc_mean_photons_7_over_8_cut";
  title = fName
    + "calculated mean photo electrons from maj 7/8 within time cut";
  //h1_tdc_mean_photons_7_over_8_cut = new TH1F(name.c_str(), title.c_str(),
  //		t_nbins[numb_cedar], channels_total, 1, channels_total + 1);
  //AddHistogram(h1_tdc_mean_photons_7_over_8_cut);

  name = fName + "_tdc_mean_photons_6_over_8_cut";
  title = fName
    + "calculated mean photo electrons from maj 6/8 within time cut";
  //h1_tdc_mean_photons_6_over_8_cut = new TH1F(name.c_str(), title.c_str(),
  //		t_nbins[numb_cedar], channels_total, 1, channels_total + 1);
  //AddHistogram(h1_tdc_mean_photons_6_over_8_cut);

  name = fName + "_tdc_chan_time_first";
  title = fName + " time distribution for first TDC per event";
  //h2_tdc_chan_time_first = new TH2F(name.c_str(),title.c_str(),channels_total,1,channels_total+1,t_nbins[numb_cedar],t_min[numb_cedar],t_max[numb_cedar]);
  //h2_tdc_chan_time_first->SetOption("col");
  //AddHistogram(h2_tdc_chan_time_first);

  if (fExpertHistos) {
    name = fName + "_tdc_chan_time_12";
    title = fName + " time distribution for channels for 12x SADC";
    //h2_tdc_chan_time_12 = new TH2F(name.c_str(),title.c_str(),channels_total,1,channels_total+1,t_nbins[numb_cedar],t_min[numb_cedar],t_max[numb_cedar]);
    //h2_tdc_chan_time_12->SetOption("col");
    //AddHistogram(h2_tdc_chan_time_12);

    name = fName + "_tdc_chan_time_echo";
    title = fName + " time distribution for channels for ev w/ echo";
    //h2_tdc_chan_time_echo = new TH2F(name.c_str(),title.c_str(),channels_total,1,channels_total+1,t_nbins[numb_cedar],t_min[numb_cedar],t_max[numb_cedar]);
    //h2_tdc_chan_time_echo->SetOption("col");
    //AddHistogram(h2_tdc_chan_time_echo);

    name = fName + "_tdc_vs_adc";
    title = fName + " #TDC hits vs. #ADC hits";
    h2_tdc_vs_adc = new TH2F(name.c_str(), title.c_str(), 15, 0, 15, 15, 0,
			     15);
    name_tdc_vs_adc = fName + "_p_tdc_vs_adc";
    //p_tdc_vs_adc = h2_tdc_vs_adc->ProfileY(name_tdc_vs_adc.c_str());
    //AddHistogram((TH1D*) p_tdc_vs_adc);

    name = fName + "_tdc_maj_vs_adc_maj";
    title = fName + " tdc maj vs adc maj";
    h2_tdc_maj_vs_adc_maj = new TH2F(name.c_str(), title.c_str(),
				     channels_total, 1, channels_total + 1, channels_total, 1,
				     channels_total + 1);
    h2_tdc_maj_vs_adc_maj->SetOption("col");
    //AddHistogram(h2_tdc_maj_vs_adc_maj);

    name = fName + "_maj6_tdc_alignment";
    title = fName + " alignment based on maj 6 events of the 8 PMs";
    h2_maj6_tdc_alignment = new TH2F(name.c_str(), title.c_str(), 4, 0.5, 4.5, 4, 0.5, 4.5);
    h2_maj6_tdc_alignment->SetOption("box text");
    AddHistogram(h2_maj6_tdc_alignment);
  }

  assert(vh1_adc_ampl_max.empty());

  vh2_g_amp_f1time = new TH2F*[channels_total];
  vh2_g_time_f1time = new TH2F*[channels_total];

  const int ampl_max = 1024;
  for (unsigned int ch = 1; ch <= channels_total; ch++) {

    char b1[222], b2[222];
    if (fExpertHistos) {
      sprintf(b1, "%s_adc_ch%d_ampl_max", fName.c_str(), ch);
      sprintf(b2, "%s Maximum (in samples) amplitude for channel %d %s",
	      fName.c_str(), ch, channel_names[ch - 1]);
      vh1_adc_ampl_max.push_back(new TH1F(b1, b2, ampl_max, 0, ampl_max));
      //AddHistogram(vh1_adc_ampl_max.back());

      sprintf(b1, "%s_adc_ch%d", fName.c_str(), ch);
      sprintf(b2, "%s ADC samples for channel %d %s", fName.c_str(), ch,
	      channel_names[ch - 1]);
      vh2_adc_samples.push_back(new TH2F(b1, b2, 32, 0, 32, ampl_max, 0,
					 ampl_max));
      //AddHistogram(vh2_adc_samples.back());
    }


    sprintf(b1, "%s_g_Amplitude_vs_f1time_chan_%d", fName.c_str(), ch);
    vh2_g_amp_f1time[ch-1] = (new TH2F(b1, b1, 4096, 0, 4095, t_nbins[numb_cedar], t_min[numb_cedar], t_max[numb_cedar]));
    vh2_g_amp_f1time[ch-1]->SetOption("col");
    AddHistogram(vh2_g_amp_f1time[ch-1]);

    sprintf(b1, "%s_g_time_vs_f1time_chan_%d", fName.c_str(), ch);
    vh2_g_time_f1time[ch-1] = (new TH2F(b1, b1, 1000, -80, 40, t_nbins[numb_cedar], t_min[numb_cedar], t_max[numb_cedar]));
    vh2_g_time_f1time[ch-1]->SetOption("colz");
    AddHistogram(vh2_g_time_f1time[ch-1]);



    sprintf(b1, "%s_adc_ch%d_ampl_max_offs", fName.c_str(), ch);
    sprintf(
	    b2,
	    "%s Maximum (in samples, offset subtracted) amplitude for channel %d %s",
	    fName.c_str(), ch, channel_names[ch - 1]);
    vh1_adc_ampl_max_offs.push_back(new TH1F_Ref(b1, b2, ampl_max, -10,
						 ampl_max, fRateCounter));
    //AddHistogram(vh1_adc_ampl_max_offs.back());

    if (fReferenceDirectory)
      ((TH1F_Ref*) vh1_adc_ampl_max_offs.back())->SetReference(
							       fReferenceDirectory);

    sprintf(b1, "%s_adc_ch%d_ampl_max_offs_when_tdc_cut", fName.c_str(), ch);
    sprintf(
	    b2,
	    "%s Maximum (in samples, offset subtracted) amplitude for channel %d %s when TDC responded",
	    fName.c_str(), ch, channel_names[ch - 1]);
    vh1_adc_ampl_max_offs_when_tdc_cut.push_back(new TH1F_Ref(b1, b2,
							      ampl_max, -10, ampl_max, fRateCounter));
    //AddHistogram(vh1_adc_ampl_max_offs_when_tdc_cut.back());
    if (fReferenceDirectory)
      ((TH1F_Ref*) vh1_adc_ampl_max_offs_when_tdc_cut.back())->SetReference(
									    fReferenceDirectory);

    if (fExpertHistos) {
      sprintf(b1, "%s_adc_offs_ch%d", fName.c_str(), ch);
      sprintf(b2, "%s ADC samples (offset subtracted) for channel %d %s",
	      fName.c_str(), ch, channel_names[ch - 1]);
      vh2_adc_samples_offs.push_back(new TH2F(b1, b2, 32, 0, 32,
					      ampl_max, -10, ampl_max));
      //AddHistogram(vh2_adc_samples_offs.back());

      sprintf(b1, "%s_adc_integrals_ch%d", fName.c_str(), ch);
      sprintf(b2, "%s ADC samples integral for channel %d %s",
	      fName.c_str(), ch, channel_names[ch - 1]);
      //vh1_adc_integrals.push_back( new TH1F(b1,b2,50, 0, 2000) );
      //AddHistogram(vh1_adc_integrals.back());

      sprintf(b1, "%s_adc_offs_not_12_ch%d", fName.c_str(), ch);
      sprintf(
	      b2,
	      "%s ADC samples (offset subtracted, not mult. 12) for channel %d %s",
	      fName.c_str(), ch, channel_names[ch - 1]);
      //vh2_adc_samples_offs_not_12.push_back( new TH2F(b1,b2,32,0,32,ampl_max,0,ampl_max) );
      //AddHistogram(vh2_adc_samples_offs_not_12.back());

      sprintf(b1, "%s_tdc_time_diff_ch%d", fName.c_str(), ch);
      sprintf(b2, "%s time diff between TDC hits for channel %d %s",
	      fName.c_str(), ch, channel_names[ch - 1]);
      //vh1_tdc_time_diff.push_back(new TH1F(b1, b2, 50, 0, 200));
      //AddHistogram(vh1_tdc_time_diff.back());

      sprintf(b1, "%s_ampl_vs_time_%d", fName.c_str(), ch);
      sprintf(b2,
	      "%s amplitude vs time for single entry events, ch %d %s",
	      fName.c_str(), ch, channel_names[ch - 1]);
      //vh2_ampl_vs_time.push_back(new TH2F(b1, b2, ampl_max+10, -10, ampl_max,
      //        t_nbins[numb_cedar], t_min[numb_cedar], t_max[numb_cedar]));
      //AddHistogram(vh2_ampl_vs_time.back());

      sprintf(b1, "%s_ampl_vs_time_2_%d", fName.c_str(), ch);
      sprintf(
	      b2,
	      "%s amplitude vs time for non-single entry events, ch %d %s",
	      fName.c_str(), ch, channel_names[ch - 1]);
      //vh2_ampl_vs_time_2.push_back(new TH2F(b1, b2, ampl_max+10, -10, ampl_max,
      //        t_nbins[numb_cedar], t_min[numb_cedar], t_max[numb_cedar]));
      //AddHistogram(vh2_ampl_vs_time_2.back());

      sprintf(b1, "%s_adc_mip_company_ch%d", fName.c_str(), ch);
      sprintf(b2, "%s ADC samples where MIP in #7 for channel %d %s",
	      fName.c_str(), ch, channel_names[ch - 1]);
      //vh2_adc_mip_company.push_back(new TH2F(b1, b2, 32, 0, 32, ampl_max,
      //		0, ampl_max));
      //AddHistogram(vh2_adc_mip_company.back());

      sprintf(b1, "%s_adc_ch%d_ampl_max_offs_when_tdc", fName.c_str(), ch);
      sprintf(
	      b2,
	      "%s Maximum (in samples, offset subtracted) amplitude for channel %d %s when TDC responded",
	      fName.c_str(), ch, channel_names[ch - 1]);
      //vh1_adc_ampl_max_offs_when_tdc.push_back(new TH1F_Ref(b1, b2,
      //		ampl_max, -10, ampl_max, fRateCounter));
      //AddHistogram(vh1_adc_ampl_max_offs_when_tdc.back());
      //if (fReferenceDirectory) ((TH1F_Ref*)vh1_adc_ampl_max_offs_when_tdc_cut.back())->SetReference(fReferenceDirectory);
    }
  }

  name = fName + "_ADC_and_TDC";
  title = fName + " ADC && TDC";
  h1_adc_and_tdc = new TH1F(name.c_str(), title.c_str(), channels_total, 1,
			    channels_total + 1);

  name = fName + "_ADC_or_TDC";
  title = fName + " ADC || TDC";
  h1_adc_or_tdc = new TH1F(name.c_str(), title.c_str(), channels_total, 1,
			   channels_total + 1);

  name = fName + "_adc_tdc_comparison";
  title = fName + " (adc && tdc) / (adc || tdc)";
  h1_adc_tdc_comparison = new TH1F_Ref(name.c_str(), title.c_str(),
				       channels_total, 1, channels_total + 1, fRateCounter, 1);
  h1_adc_tdc_comparison->SetMaximum(1);
  h1_adc_tdc_comparison->SetMinimum(0);
  //AddHistogram(h1_adc_tdc_comparison);
  if (fReferenceDirectory)
    ((TH1F_Ref*) h1_adc_tdc_comparison)->SetReference(fReferenceDirectory);

  // Set channel names
  // 	cout << " and setting channels names " << endl;
  for (unsigned int i = 1; i <= channels_total; i++) {
    char buf[55];
    sprintf(buf, "%d %s", i, channel_names[i - 1]);
    h1_adc_chan -> GetXaxis() -> SetBinLabel(i, buf);
    h1_tdc_chan -> GetXaxis() -> SetBinLabel(i, buf);
    h1_tdc_chan_hit -> GetXaxis() -> SetBinLabel(i, buf);
    h2_tdc_chan_time -> GetXaxis() -> SetBinLabel(i, buf);
    h2_tdc_chan_time_all -> GetXaxis() -> SetBinLabel(i, buf);
    h1_tdc_maj_6_chan -> GetXaxis() -> SetBinLabel(i, buf);
    //h2_tdc_chan_time_first -> GetXaxis() -> SetBinLabel(i,buf);
    h1_adc_and_tdc -> GetXaxis() -> SetBinLabel(i, buf);
    h1_adc_or_tdc -> GetXaxis() -> SetBinLabel(i, buf);
    h1_adc_tdc_comparison -> GetXaxis() -> SetBinLabel(i, buf);
    if (fExpertHistos) {
      h2_tdc_maj_vs_adc_maj -> GetXaxis() -> SetBinLabel(i, buf);
      h2_tdc_maj_vs_adc_maj -> GetYaxis() -> SetBinLabel(i, buf);
    }

    h1_ghits_chan -> GetXaxis() -> SetBinLabel(i, buf);
    h2_g_chan_time_all -> GetXaxis() -> SetBinLabel(i, buf);
    h2_g_chan_ampl -> GetXaxis() -> SetBinLabel(i, buf);
    h2_g_chan_int -> GetXaxis() -> SetBinLabel(i, buf);
  }

  //////////////////////////////////////////////////////////////
  /////////////////// PlaneGandalf plots for experts

  if (fExpertHistos) {
    for(int i=0;i<8;i++) {
      stringstream name;
      name<< fName << "g_" << i;
      g_planes[i] = new PlaneGandalf(name.str().c_str());
      g_planes[i]->Init(tree);
      for( vector<TH1*>::iterator it = g_planes[i]->GetHistoList().begin(); it != g_planes[i]->GetHistoList().end();it++)
	AddHistogram((*it));
    }
  }
  // 	cout << "done" << endl;
}

void PlaneCEDAR::Reset(void) {
  Plane::Reset();

  channel_infos.clear();
  g_channel_infos.clear();

  digits_sadc.clear();
  digits_f1.clear();
  digits_gandalf.clear();
  for (int chan = 0; chan < GetNchannels(); chan++) {
    multiplicity_counter[chan] = 0.;
  }
  if (fExpertHistos) {
    for (map<int, PlaneGandalf*>::iterator plane=g_planes.begin();plane!=g_planes.end();++plane)
      (*plane).second->Reset();
  }
}

/*
  void PlaneCEDAR::ControlPanel(const TGWindow* p, const TGWindow* main) {
  if (!fControlPanel)
  fControlPanel = new PlaneCEDARPanel(p, main, 100, 100, this);
  }
*/

void PlaneCEDAR::put_in_hist(TH2F* h2, bool subtract_offs,
			     const ChannelInfo& ci) {
  if (!ci.digit_sadc)
    return;
  const vector<CS::uint16> vsamples = ci.digit_sadc->GetSamples();
  float offs = subtract_offs ? ci.offs : 0;
  double nentries = h2->GetEntries();
  for (size_t i = 0; i < vsamples.size(); i++)
    h2->Fill(i, vsamples[i] - offs);
  h2->SetEntries(nentries + 1);
}


void PlaneCEDAR::StoreDigitGandalf(const CS::ChipGandalf::DigitGADC* digit)
{

  // TODO: for now, skip frame digits
  if (digit->getNumSamples()>0) return;

  //add digit to digits list
  digits_gandalf.push_back(digit);

  // time part, mostly taken over from PlaneCEDAR;
  // TODO: what does it have todo with sandwich ??
  int ch = digit->getChannel() + 1; // since Sandwich is counted 1-8 and not 0-7 like CEDAR
  g_channel_infos[ch].digits.push_back(digit);

  h1_ghits_chan->Fill(ch);
  h1_g_time->Fill(digit->GetTimeDecoded());

  h2_g_chan_time_all->Fill(ch, digit->GetTimeDecoded());

  h2_g_chan_ampl->Fill(ch, digit->getMaxAmplitude());
        
  h2_g_chan_int->Fill(ch, digit->getIntegral());

  g_amplitudes[ch - 1] = digit->getMaxAmplitude();

  // TODO: what is a pileup ?
  // only write the sample histos when we have a dbg frame
  //	if (fExpertHistos && !chi.pileup && (digitG->getOpMode()==CS::ChipGandalf::GADC_DEBUG
  //			|| digitG->getOpMode()==CS::ChipGandalf::GADC_DEBUG_IL)) {
  //		vh1_adc_ampl_max[ch - 1]->Fill(digitG->getMaxAmplitude());
  //		//vh1_adc_integrals[ch-1]->Fill(integral);
  //		put_in_hist(vh2_adc_samples[ch - 1], false, digitG);
  //		put_in_hist(vh2_adc_samples_offs[ch - 1], true, digitG);
  //	}

}

void PlaneCEDAR::StoreDigitSADC(const CS::ChipSADC::Digit* d_sadc) {
  digits_sadc.push_back(d_sadc);
  //if( debug )
  //d_sadc->Print();

  if (d_sadc->GetY() != 0)
    printf("PlaneCEDAR::StoreDigit(): y=%d: it is ignored!!\n",
	   d_sadc->GetY());

  int ch = d_sadc->GetX() + 1; // SADC channels 0-7
  if (ch < 1 || ch > (int) channels_total) {
    printf("PlaneCEDAR::StoreDigit(): bad channel number %d: ignored!!\n",
	   ch);
    return;
  }

  h1_adc_chan->Fill(ch);

  const size_t window_width = 8; // First window_width samples are used to determine the offset.
  double sum = 0;
  double sum2 = 0;
  float integral = 0;
  CS::uint16 ampl_max = 0;
  CS::uint16 ampl_min = 1023;
  const vector<CS::uint16>& vsamples = d_sadc->GetSamples();
  for (size_t i = 0; i < vsamples.size(); i++) {
    // Calculate average + fluctuation only for first few samples
    if (i < window_width) {
      sum += vsamples[i];
      sum2 += vsamples[i] * vsamples[i];
    }

    integral += vsamples[i];

    if (ampl_max <= vsamples[i])
      ampl_max = vsamples[i];
    if (ampl_min >= vsamples[i])
      ampl_min = vsamples[i];
  }

  float offset = sum / window_width;
  float rms = sum2 / window_width - offset * offset;

  integral -= vsamples.size() * offset;

  ChannelInfo& chi = channel_infos[ch];
  assert(!chi.digit_sadc);

  chi.min = ampl_min;
  chi.max = ampl_max;
  chi.offs = offset;
  chi.rms = rms;
  chi.integral = integral;
  chi.pileup = rms >= 4;
  chi.digit_sadc = d_sadc;

  vh1_adc_ampl_max_offs[ch - 1]->Fill(ampl_max - offset);
  adc_amplitudes[ch - 1] = ampl_max - offset;
  if (fExpertHistos && !chi.pileup) {
    vh1_adc_ampl_max[ch - 1]->Fill(ampl_max);
    //vh1_adc_integrals[ch-1]->Fill(integral);
    put_in_hist(vh2_adc_samples[ch - 1], false, chi);
    put_in_hist(vh2_adc_samples_offs[ch - 1], true, chi);
  }
}

void PlaneCEDAR::StoreDigitF1(const CS::ChipF1::Digit* d_f1) {
  digits_f1.push_back(d_f1);
  //std::cout << " f1 channel " << d_f1->GetChannel() << std::endl;
  assert(d_f1->GetChannel() >= 0 && d_f1->GetChannel()
	 <= (int) channels_total - 1);// since Sandwich is counted 1-8 and not 0-7 like CEDAR
  //printf("chan: %d %g\n",d_f1->GetChannel(),d_f1->GetTimeDecoded());

  int ch = d_f1->GetChannel() + 1; // since Sandwich is counted 1-8 and not 0-7 like CEDAR
  channel_infos[ch].digits_f1.push_back(d_f1);

  h1_tdc_chan->Fill(ch);
  if (fT[numb_cedar]->Test(d_f1->GetTimeDecoded())) {
    h1_tdc_chan_hit->Fill(ch);
    hit_tdc_cut[ch - 1] = true;
    h2_tdc_chan_time->Fill(ch, d_f1->GetTimeDecoded());
  } // looking at events in time cut
  h1_tdc_time->Fill(d_f1->GetTimeDecoded());
  h2_tdc_chan_time_all->Fill(ch, d_f1->GetTimeDecoded());
}

void PlaneCEDAR::StoreDigit(CS::Chip::Digit* digit) {
  //printf("PlaneCEDAR::StoreDigit(): call!\n");
  //digit->Print();
  const CS::ChipSADC::Digit* d_sadc =
    dynamic_cast<const CS::ChipSADC::Digit*> (digit);
  if (d_sadc != NULL) {
    //std::cout << " found sadc digit " << std::endl;
    StoreDigitSADC(d_sadc);
    return;
  }

  const CS::ChipF1::Digit* d_f1 =
    dynamic_cast<const CS::ChipF1::Digit*> (digit);
  if (d_f1 != NULL) {
    //std::cout << " found tdc digit " << std::endl;
    StoreDigitF1(d_f1);
    return;
  }

  const CS::ChipGandalf::DigitGADC* d_g =
    dynamic_cast<const CS::ChipGandalf::DigitGADC*> (digit);
  if (d_g != NULL) {
    //std::cout << " found tdc digit " << std::endl;
    StoreDigitGandalf(d_g);
    if (fExpertHistos)
      g_planes[d_g->getChannel()]->StoreDigit(digit);
    return;
  }

  cout << " Error in PlaneCEDAR::StoreDigit(): digit could not get casted "
       << endl;
}

void PlaneCEDAR::EndEvent(const CS::DaqEvent &event) {
  //vector<const ChipSADC::Digit*>::const_iterator sadc;
  //vector<const ChipF1::Digit*>::const_iterator f1;

  //printf("PlaneCEDAR::EndEvent(): call!\n");
  //printf("%d %d\n",digits_sadc.size(),digits_f1.size());
  if (fExpertHistos) {
    if (digits_f1.size() > 0 || digits_sadc.size() > 0)
      h2_tdc_vs_adc->Fill(digits_f1.size(), digits_sadc.size());

    h2_tdc_vs_adc->ProfileY(name_tdc_vs_adc.c_str());
  }

  if (digits_sadc.size() > 0) {
    h1_adc_multiplicity->Fill(digits_sadc.size());
    if (digits_sadc.size() > 8)
      puts("More ADC hits in CEDAR than channels: impossible.\n");

    if (fExpertHistos) {
      if (digits_sadc.size() == 8) {
	for (size_t i = 0; i < digits_f1.size(); i++) {
	  //h2_tdc_chan_time_12->Fill(digits_f1[i]->GetChannel()+1,
	  //digits_f1[i]->GetTimeDecoded());
	}
      } else {
	for (map<int, ChannelInfo>::const_iterator i =
	       channel_infos.begin(); i != channel_infos.end(); i++) {
	  // int ch = i->first;
	  const ChannelInfo& chi = i->second;
	  if (chi.digit_sadc && !chi.pileup) {
	    //put_in_hist(vh2_adc_samples_offs_not_12[ch-1], true, chi);
	  }
	}
      }
    }
  }

  // Sort by arrival time.
  sort(digits_f1.begin(), digits_f1.end(), earlier());
  for (size_t i = 1; i <= channels_total; i++) {
    sort(channel_infos[i].digits_f1.begin(),
	 channel_infos[i].digits_f1.end(), earlier());
  }

  if (digits_f1.size() > 0 && digits_sadc.size() < 8) {
    // int ch = digits_f1[0]->GetChannel() + 1;
    double t = digits_f1[0]->GetTimeDecoded();
    //h2_tdc_chan_time_first->Fill(ch, t);

    if (fExpertHistos) {
      // Look at 'echos'.
      bool found = false;
      for (vector<const CS::ChipF1::Digit*>::const_iterator f1 =
	     channel_infos[5].digits_f1.begin(); f1
	     != channel_infos[5].digits_f1.end(); f1++) {
	t = (*f1)->GetTimeDecoded();
	if (t > -888 && t < -874) {
	  //h2_tdc_chan_time_echo->Fill(5, t);
	  found = true;
	  break;
	}
      }

      if (found) {
	vector<const CS::ChipF1::Digit*> digits_in_window;
	for (size_t i = 0; i < digits_f1.size(); i++) {
	  // int ch = digits_f1[i]->GetChannel() + 1;
	  double t = digits_f1[i]->GetTimeDecoded();
	  if (t < -930 || t > -900)
	    continue;
	  //h2_tdc_chan_time_echo->Fill(ch, t);
	  //digits_in_window.push_back(digits_f1[i]);
	}
#if 0
	if (digits_in_window.size() == 1) {
	  int ch = digits_in_window[0]->GetChannel()+1;
	  double t = digits_in_window[0]->GetTimeDecoded();
	  //h2_tdc_chan_time_echo->Fill(ch, t);
	}
#endif
      }
    }
  }

  if (digits_f1.size() > 1) {
    double first[channels_total];
    bool found_first[channels_total];

    memset(found_first, 0, sizeof(found_first));
    for (size_t i = 1; i < digits_f1.size(); i++) {
      int ch = digits_f1[i]->GetChannel() + 1;
      //std::cout << first[ch-1] << " ";
      //if (fT[numb_cedar]->Test(digits_f1[i]->GetTimeDecoded())) {h1_tdc_chan_hit->Fill(ch);} // looking at events in time cut
      if (!found_first[ch - 1]) {
	found_first[ch - 1] = true;
	first[ch - 1] = digits_f1[i]->GetTimeDecoded();
	//if (fT[numb_cedar]->Test(first[ch-1])) {h1_tdc_chan_hit->Fill(ch);}
      } else {
	if (fExpertHistos) {
	  // double current = digits_f1[i]->GetTimeDecoded();
	  //vh1_tdc_time_diff[ch-1]->Fill(current - first[ch-1]);
	}
      }
    }
    //std::cout << endl;
  }

  // Compare number of events where both adc and tdc fired to number of events
  // where either fired.
  double nentries = h1_adc_tdc_comparison->GetEntries();
  for (size_t ch = 1; ch <= channels_total; ch++) {
    bool tdc_hit = channel_infos[ch].digits_f1.size() > 0;
    bool adc_hit = channel_infos[ch].digit_sadc != 0;
    if (tdc_hit || adc_hit) {
      h1_adc_or_tdc->Fill(ch);
      if (tdc_hit && adc_hit)
	h1_adc_and_tdc->Fill(ch);
      h1_adc_tdc_comparison ->SetBinContent(ch,
					    h1_adc_and_tdc->GetBinContent(ch)
					    / h1_adc_or_tdc->GetBinContent(ch));
      h1_adc_tdc_comparison->SetEntries(nentries + 1);
    }
  }

  if (fExpertHistos) {
    // For channels with exactly two TDC hit and SADC data, plot
    // maximum against time
    for (size_t ch = 1; ch <= channels_total; ch++) {
      int tdc_hit = channel_infos[ch].digits_f1.size();
      bool adc_hit = channel_infos[ch].digit_sadc != 0;

      if (!tdc_hit || !adc_hit)
	continue;

      // Find first TDC hit and the ADC hit for this channel.
      const CS::ChipF1::Digit* f1 = channel_infos[ch].digits_f1.front();
      assert(f1);

      if (channel_infos[ch].pileup)
	continue; // Ignore pileup.

      // float offset = channel_infos[ch].offs;
      // CS::uint16 ampl_max = channel_infos[ch].max;
      // CS::uint16 ampl_min = channel_infos[ch].min;

      // float time = f1->GetTimeDecoded();

      if (tdc_hit == 1) {
	//vh2_ampl_vs_time[ch-1]->Fill(ampl_max - offset, time);
	//vh2_ampl_vs_time[ch-1]->Fill(ampl_min - offset, time);
      } else {
	//vh2_ampl_vs_time_2[ch-1]->Fill(ampl_max - offset, time);
	//vh2_ampl_vs_time_2[ch-1]->Fill(ampl_min - offset, time);
      }
    }

    // See if there was a Mip in channel 7 and no pileup, if yes
    // plot scope pictures where available.
    if (!channel_infos[7].pileup && channel_infos[7].max
	- channel_infos[7].offs >= 20 && channel_infos[7].max
	- channel_infos[7].offs <= 35) {
      for (int ch = 1; ch <= 12; ch++) {
	const ChannelInfo& chi = channel_infos[ch];
	if (chi.digit_sadc) {
	  const vector<CS::uint16>& vsamples =
	    chi.digit_sadc->GetSamples();
	  for (size_t j = 0; j < vsamples.size(); j++) {
	    //	vh2_adc_mip_company[ch-1]->Fill(j, vsamples[j] - chi.offs);
	  }
	}
      }
    }

  }

  int tdc_mult(0);
  int adc_mult(0);
  for (int chan = 0; chan < GetNchannels(); chan++) {
    if (hit_tdc_cut[chan]) {
      tdc_mult++;
    }
    if (adc_amplitudes[chan] != 0.) {
      adc_mult++;
    }
  }
  // show which channels fired in the case of 6 fold multiplicities
  if (tdc_mult >= 6)
    for (int chan = 0; chan < GetNchannels(); chan++) {
      if (hit_tdc_cut[chan]) {
	int _pm = chan+1;
	h1_tdc_maj_6_chan->Fill(_pm);
	// show the alignment for the 6-fold case
	if (fExpertHistos) switch (_pm){
	  case 1: h2_maj6_tdc_alignment->Fill(3,4); break; // PM1
	  case 2: h2_maj6_tdc_alignment->Fill(4,3); break; // PM2
	  case 3: h2_maj6_tdc_alignment->Fill(4,2); break;
	  case 4: h2_maj6_tdc_alignment->Fill(3,1); break;
	  case 5: h2_maj6_tdc_alignment->Fill(2,1); break;
	  case 6: h2_maj6_tdc_alignment->Fill(1,2); break;
	  case 7: h2_maj6_tdc_alignment->Fill(1,3); break;
	  case 8: h2_maj6_tdc_alignment->Fill(2,4); break; // PM8
	  }

      }
    }

  if (fExpertHistos) {
    h2_tdc_maj_vs_adc_maj->Fill(tdc_mult, adc_mult);
  }

  if (tdc_mult == 0) h1_tdc_maj_cut->Fill(0.0);
  for (int chan = 0; chan < GetNchannels(); chan++) {
    if (hit_tdc_cut[chan] && adc_amplitudes[chan] != 0.) {
      vh1_adc_ampl_max_offs_when_tdc_cut[chan]->Fill(adc_amplitudes[chan]);
    }
    // using chan index also for multiplicities
    if ((chan+1) <= tdc_mult) {
      h1_tdc_maj_cut->Fill(chan+1);
      // h1_tdc_maj_cut->Fill(chan);
    }
    adc_amplitudes[chan] = 0.;
    hit_tdc_cut[chan] = false;
  }

  //////////////////////////////////////////////////////////////////////////////////////
  ///////////////////////////////////// Gandalf Readout
  if (digits_gandalf.size() > 0) {
    h1_g_multiplicity->Fill(digits_gandalf.size());
  }

  // Sort by arrival time.
  sort(digits_gandalf.begin(), digits_gandalf.end(), earlierG());
  for (size_t i = 1; i <= channels_total; i++) {
    sort(g_channel_infos[i].digits.begin(),
	 g_channel_infos[i].digits.end(), earlierG());
    h2_g_mult_chan->Fill(i,g_channel_infos[i].digits.size());
  }

  nentries = h1_gandalf_tdc_comparison->GetEntries();
  for (size_t ch = 1; ch <= channels_total; ch++) {
    bool tdc_hit = channel_infos[ch].digits_f1.size() > 0;
    bool adc_hit = g_channel_infos[ch].digits.size() > 0;
    if (tdc_hit || adc_hit) {
      h1_gandalf_or_tdc->Fill(ch);
      if (tdc_hit && adc_hit) {
	h1_gandalf_and_tdc->Fill(ch);

	//compare the first gandalf hit with the f1 time
	for (vector<const CS::ChipF1::Digit*>::const_iterator f1 =
	       channel_infos[ch].digits_f1.begin(); f1
	       != channel_infos[ch].digits_f1.end(); ++f1) {
	  vh2_g_amp_f1time[ch-1]->Fill( g_channel_infos[ch].digits[0]->getMaxAmplitude(), (*f1)->GetTimeDecoded());
	  vh2_g_time_f1time[ch-1]->Fill( g_channel_infos[ch].digits[0]->GetTimeDecoded(), (*f1)->GetTimeDecoded());
	}
      }
      h1_gandalf_tdc_comparison ->SetBinContent(ch,
						h1_gandalf_and_tdc->GetBinContent(ch)
						/ h1_gandalf_or_tdc->GetBinContent(ch));
      h1_gandalf_tdc_comparison->SetEntries(nentries + 1);
    }

    if (adc_hit && tdc_hit)
      h2DataTransmission->Fill(ch-1,0);
    if (adc_hit && !tdc_hit)
      h2DataTransmission->Fill(ch-1,1);
    if (!adc_hit && tdc_hit)
      h2DataTransmission->Fill(ch-1,2);
    if (!adc_hit && !tdc_hit)
      h2DataTransmission->Fill(ch-1,3);
  }

  //// DEAD CODE ! ///////////////////////////////////

  if (digits_gandalf.size() > 0) {
    const CS::ChipGandalf::DigitGADC* dig = *digits_gandalf.begin();
    // int ch = dig->getChannel() + 1;
    double t = dig->GetTimeDecoded();
    //h2_tdc_chan_time_first->Fill(ch, t);

    if (fExpertHistos) {
      // Look at 'echos'.
      bool found = false;
      // TODO: change variable names, review if this code makes sense for gandalf
      for (vector<const CS::ChipGandalf::DigitGADC*>::const_iterator f1 =
	     g_channel_infos[5].digits.begin(); f1
	     != g_channel_infos[5].digits.end(); f1++) {

	t = (*f1)->GetTimeDecoded();
	if (t > -888 && t < -874) {
	  //h2_tdc_chan_time_echo->Fill(5, t);
	  found = true;
	  break;
	}
      }

      if (found) {
	vector<const CS::ChipGandalf::DigitGADC*> digits_in_window;
	for (std::vector<const CS::ChipGandalf::DigitGADC*>::iterator i = digits_gandalf.begin(); i != digits_gandalf.end(); ++i) {
	  // int ch = (*i)->getChannel() + 1;
	  double t = (*i)->GetTimeDecoded();
	  if (t < -930 || t > -900)
	    continue;
	  //h2_tdc_chan_time_echo->Fill(ch, t);
	  //digits_in_window.push_back(digits_f1[i]);
	}
      }
    }
  }

  if (digits_gandalf.size() > 1) {
    double first[channels_total];
    bool found_first[channels_total];

    memset(found_first, 0, sizeof(found_first));
    for (std::vector<const CS::ChipGandalf::DigitGADC*>::iterator i = digits_gandalf.begin()++; i != digits_gandalf.end(); ++i) {
      int ch = (*i)->getChannel() + 1;

      //std::cout << first[ch-1] << " ";
      //if (fT[numb_cedar]->Test(digits_f1[i]->GetTimeDecoded())) {h1_tdc_chan_hit->Fill(ch);} // looking at events in time cut
      if (!found_first[ch - 1]) {
	found_first[ch - 1] = true;
	first[ch - 1] = (*i)->GetTimeDecoded();
	//if (fT[numb_cedar]->Test(first[ch-1])) {h1_tdc_chan_hit->Fill(ch);}
      } else {
	if (fExpertHistos) {
	  // double current = (*i)->GetTimeDecoded();
	  //vh1_tdc_time_diff[ch-1]->Fill(current - first[ch-1]);
	}
      }
    }
  }

  if (fExpertHistos) {
    // For channels with exactly two TDC hit and SADC data, plot
    // maximum against time
    for (size_t ch = 1; ch <= channels_total; ch++) {
      if (g_channel_infos[ch].digits.size()==0) continue;

      // Find first hit for this channel.
      const CS::ChipGandalf::DigitGADC* hit = g_channel_infos[ch].digits.front();
      assert(hit);

      // TODO: what is pileup ??
      //			if (channel_infos[ch].pileup)
      //				continue; // Ignore pileup.

      // float offset = hit->getBaseLine();
      // CS::uint16 ampl_max = hit->getMaxAmplitude();
      // TODO: gandalf doesn't calculate MinAmplitude CS::uint16 ampl_min = hit->get;

      // float time = hit->GetTimeDecoded();

    }
    //		// TODO:
    //		// See if there was a Mip in channel 7 and no pileup, if yes
    //		// plot scope pictures where available.
    //		if (!channel_infos[7].pileup && channel_infos[7].max
    //				- channel_infos[7].offs >= 20 && channel_infos[7].max
    //				- channel_infos[7].offs <= 35) {
    //			for (int ch = 1; ch <= 12; ch++) {
    //				const ChannelInfo& chi = channel_infos[ch];
    //				if (chi.digit_sadc) {
    //					const vector<CS::uint16>& vsamples =
    //							chi.digit_sadc->GetSamples();
    //					for (size_t j = 0; j < vsamples.size(); j++) {
    //						//	vh2_adc_mip_company[ch-1]->Fill(j, vsamples[j] - chi.offs);
    //					}
    //				}
    //			}
    //		}

    ////////////////////////////////////////////////


    // fill the perchannel expert histos
    for (int i=0;i<8;i++) {
      g_planes[i]->EndEvent(event);
    }
  }

  Plane::EndEvent(event);

}
