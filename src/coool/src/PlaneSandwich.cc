#include "PlaneSandwich.h"

#include "PlaneSandwichPanel.h"

#include <algorithm>

ClassImp(PlaneSandwich);

// unsigned int PlaneSandwich::channels_total = 12;

PlaneSandwich::PlaneSandwich(const char *detname) :
    Plane(detname)
{
  //printf("PlaneSandwich::PlaneSandwich(): %s\n",detname);
  fT = AddVariable("time_range", 100, -960, -860, 100);
}

PlaneSandwich::~PlaneSandwich(void)
{
    //printf("PlaneSandwich::~PlaneSandwich()\n");
  delete h1_adc_or_tdc;
  delete h1_adc_and_tdc;

  delete h2_tdc_vs_adc;
}

void PlaneSandwich::Init(TTree* tree)
{
    const char *channel_names[channels_total]={"TopJura","TopJC","TopCenter","TopCS","TopSaleve",
                                               "Jura","Saleve",
                                               "BotJura","BotJC","BotCenter","BotCS","BotSaleve"};

    std::string name,title;
    
    float t_min = fT->GetMin();
    float t_max = fT->GetMax();
    int t_nbins = fT->GetNbins();
    
    name = fName + "_adc_chan";
    title = fName + " SADC no of entries";
    h1_adc_chan = new TH1F_Ref(name.c_str(),title.c_str(),channels_total,1,channels_total+1,fRateCounter);
    h1_adc_chan->SetMinimum(0);
    AddHistogram(h1_adc_chan);
    if (fReferenceDirectory) ((TH1F_Ref*)h1_adc_chan)->SetReference(fReferenceDirectory);

    name = fName + "_adc_multiplicity";
    title = fName + " SADC multiplicity distribution";
    h1_adc_multiplicity = new TH1F_Ref(name.c_str(),title.c_str(),channels_total , 1,channels_total+1,fRateCounter);
    h1_adc_multiplicity->SetMinimum(0);
    AddHistogram(h1_adc_multiplicity);
    if (fReferenceDirectory) ((TH1F_Ref*)h1_adc_multiplicity)->SetReference(fReferenceDirectory);

    name = fName + "_tdc_chan";
    title = fName + " TDC channel total no hits";
    h1_tdc_chan = new TH1F_Ref(name.c_str(),title.c_str(),channels_total,1,channels_total+1,fRateCounter);
    h1_tdc_chan->SetMinimum(0);
    AddHistogram(h1_tdc_chan);
    if (fReferenceDirectory) ((TH1F_Ref*)h1_tdc_chan)->SetReference(fReferenceDirectory);

    name = fName + "_tdc_chan_yes_no";
    title = fName + " TDC channel which send data";
    h1_tdc_chan_hit = new TH1F_Ref(name.c_str(),title.c_str(),channels_total,1,channels_total+1,fRateCounter);
    h1_tdc_chan_hit->SetMinimum(0);
    AddHistogram(h1_tdc_chan_hit);
    if (fReferenceDirectory) ((TH1F_Ref*)h1_tdc_chan_hit)->SetReference(fReferenceDirectory);

    name = fName + "_tdc_time";
    title = fName + " time distribution";
    h1_tdc_time = new TH1F_Ref(name.c_str(),title.c_str(),t_nbins,t_min,t_max,fRateCounter);
    AddHistogram(h1_tdc_time);
    if (fReferenceDirectory) ((TH1F_Ref*)h1_tdc_time)->SetReference(fReferenceDirectory);

    name = fName + "_tdc_chan_time";
    title = fName + " time distribution for channels";
    h2_tdc_chan_time = new TH2F(name.c_str(),title.c_str(),channels_total,1,channels_total+1,t_nbins,t_min,t_max);
    h2_tdc_chan_time->SetOption("col");
    AddHistogram(h2_tdc_chan_time);

    name = fName + "_tdc_chan_time_first";
    title = fName + " time distribution for first TDC per event";
    h2_tdc_chan_time_first = new TH2F(name.c_str(),title.c_str(),channels_total,1,channels_total+1,t_nbins,t_min,t_max);
    h2_tdc_chan_time_first->SetOption("col");
    AddHistogram(h2_tdc_chan_time_first);

    if (fExpertHistos) {
      name = fName + "_tdc_chan_time_12";
      title = fName + " time distribution for channels for 12x SADC";
      h2_tdc_chan_time_12 = new TH2F(name.c_str(),title.c_str(),channels_total,1,channels_total+1,t_nbins,t_min,t_max);
      h2_tdc_chan_time_12->SetOption("col");
      AddHistogram(h2_tdc_chan_time_12);

      name = fName + "_tdc_chan_time_echo";
      title = fName + " time distribution for channels for ev w/ echo";
      h2_tdc_chan_time_echo = new TH2F(name.c_str(),title.c_str(),channels_total,1,channels_total+1,t_nbins,t_min,t_max);
      h2_tdc_chan_time_echo->SetOption("col");
      AddHistogram(h2_tdc_chan_time_echo);

      name = fName + "_tdc_vs_adc";
      title = fName + " #TDC hits vs. #ADC hits";
      h2_tdc_vs_adc = new TH2F(name.c_str(),title.c_str(),15,0,15,15,0,15);
      name_tdc_vs_adc = fName + "_p_tdc_vs_adc";
      p_tdc_vs_adc = h2_tdc_vs_adc->ProfileY(name_tdc_vs_adc.c_str());
      AddHistogram((TH1D*)p_tdc_vs_adc);
    }

    assert(vh1_adc_ampl_max.empty());
    
    const int ampl_max=1024;
    for( unsigned int ch=1; ch<=channels_total; ch++ )
    {
        char b1[222],b2[222];
	if (fExpertHistos) {
	  sprintf(b1,"%s_adc_ch%d_ampl_max",fName.c_str(),ch);
	  sprintf(b2,"%s Maximum (in samples) amplitude for channel %d %s",fName.c_str(),ch,channel_names[ch-1]);
	  vh1_adc_ampl_max.push_back(new TH1F(b1,b2,ampl_max,0,ampl_max));
	  AddHistogram(vh1_adc_ampl_max.back());

	  sprintf(b1,"%s_adc_ch%d",fName.c_str(),ch);
	  sprintf(b2,"%s ADC samples for channel %d %s",fName.c_str(),ch,channel_names[ch-1]);
	  vh2_adc_samples.push_back( new TH2F(b1,b2,32,0,32,ampl_max,0,ampl_max) );
	  AddHistogram(vh2_adc_samples.back());
	}
	
        sprintf(b1,"%s_adc_ch%d_ampl_max_offs",fName.c_str(),ch);
        sprintf(b2,"%s Maximum (in samples, offset subtracted) amplitude for channel %d %s",fName.c_str(),ch,channel_names[ch-1]);
        vh1_adc_ampl_max_offs.push_back(new TH1F_Ref(b1,b2,ampl_max,-10,ampl_max,fRateCounter));
        AddHistogram(vh1_adc_ampl_max_offs.back());
        if (fReferenceDirectory) ((TH1F_Ref*)vh1_adc_ampl_max_offs.back())->SetReference(fReferenceDirectory);

	if (fExpertHistos) {
	  sprintf(b1,"%s_adc_offs_ch%d",fName.c_str(),ch);
	  sprintf(b2,"%s ADC samples (offset subtracted) for channel %d %s",fName.c_str(),ch,channel_names[ch-1]);
	  vh2_adc_samples_offs.push_back( new TH2F(b1,b2,32,0,32,ampl_max,-10,ampl_max) );
	  AddHistogram(vh2_adc_samples_offs.back());

	  sprintf(b1,"%s_adc_integrals_ch%d",fName.c_str(),ch);
	  sprintf(b2,"%s ADC samples integral for channel %d %s",fName.c_str(),ch,channel_names[ch-1]);
	  vh1_adc_integrals.push_back( new TH1F(b1,b2,50, 0, 2000) );
	  AddHistogram(vh1_adc_integrals.back());

	  sprintf(b1,"%s_adc_offs_not_12_ch%d",fName.c_str(),ch);
	  sprintf(b2,"%s ADC samples (offset subtracted, not mult. 12) for channel %d %s",fName.c_str(),ch,channel_names[ch-1]);
	  vh2_adc_samples_offs_not_12.push_back( new TH2F(b1,b2,32,0,32,ampl_max,0,ampl_max) );
	  AddHistogram(vh2_adc_samples_offs_not_12.back());

	  sprintf(b1,"%s_tdc_time_diff_ch%d", fName.c_str(), ch);
	  sprintf(b2,"%s time diff between TDC hits for channel %d %s",fName.c_str(), ch, channel_names[ch-1]);
	  vh1_tdc_time_diff.push_back(new TH1F(b1, b2, 50, 0, 200));
	  AddHistogram(vh1_tdc_time_diff.back());

	  sprintf(b1, "%s_ampl_vs_time_%d", fName.c_str(), ch);
	  sprintf(b2, "%s amplitude vs time for single entry events, ch %d %s",
		  fName.c_str(), ch, channel_names[ch-1]);
	  vh2_ampl_vs_time.push_back(new TH2F(b1, b2, ampl_max+10, -10, ampl_max,
					      t_nbins, t_min, t_max));
	  AddHistogram(vh2_ampl_vs_time.back());

	  sprintf(b1, "%s_ampl_vs_time_2_%d", fName.c_str(), ch);
	  sprintf(b2, "%s amplitude vs time for non-single entry events, ch %d %s",
		  fName.c_str(), ch, channel_names[ch-1]);
	  vh2_ampl_vs_time_2.push_back(new TH2F(b1, b2, ampl_max+10, -10, ampl_max,
					      t_nbins, t_min, t_max));
	  AddHistogram(vh2_ampl_vs_time_2.back());

          sprintf(b1,"%s_adc_mip_company_ch%d",fName.c_str(),ch);
	  sprintf(b2,"%s ADC samples where MIP in #7 for channel %d %s",fName.c_str(),ch,channel_names[ch-1]);
	  vh2_adc_mip_company.push_back( new TH2F(b1,b2,32,0,32,ampl_max,0,ampl_max) );
	  AddHistogram(vh2_adc_mip_company.back());

	}
    }

    name = fName + "_ADC_and_TDC";
    title = fName + " ADC && TDC";
    h1_adc_and_tdc = new TH1F(name.c_str(), title.c_str(),
                                  channels_total, 1, channels_total + 1);

    name = fName + "_ADC_or_TDC";
    title = fName + " ADC || TDC";
    h1_adc_or_tdc = new TH1F(name.c_str(), title.c_str(),
                                 channels_total, 1, channels_total + 1);

    name = fName + "_adc_tdc_comparison";
    title = fName + " (adc && tdc) / (adc || tdc)";
    h1_adc_tdc_comparison = new TH1F_Ref(name.c_str(), title.c_str(),
                                    channels_total, 1, channels_total + 1,
                                    fRateCounter, 1);
    h1_adc_tdc_comparison->SetMaximum(1);
    h1_adc_tdc_comparison->SetMinimum(0);
    AddHistogram(h1_adc_tdc_comparison);
    if (fReferenceDirectory) ((TH1F_Ref*)h1_adc_tdc_comparison)->SetReference(fReferenceDirectory);

    // Set channel names
    for( unsigned int i=1; i<=channels_total; i++ )
    {
        char buf[55];
        sprintf(buf,"%d %s",i,channel_names[i-1]);
        h1_adc_chan -> GetXaxis() -> SetBinLabel(i,buf);
        h1_tdc_chan -> GetXaxis() -> SetBinLabel(i,buf);
        h1_tdc_chan_hit -> GetXaxis() -> SetBinLabel(i,buf);
        h2_tdc_chan_time -> GetXaxis() -> SetBinLabel(i,buf);
        h2_tdc_chan_time_first -> GetXaxis() -> SetBinLabel(i,buf);
        h1_adc_and_tdc -> GetXaxis() -> SetBinLabel(i,buf);
        h1_adc_or_tdc -> GetXaxis() -> SetBinLabel(i,buf);
        h1_adc_tdc_comparison -> GetXaxis() -> SetBinLabel(i,buf);
    }
}

void PlaneSandwich::Reset(void)
{
    Plane::Reset();

    channel_infos.clear();

    digits_sadc.clear();
    digits_f1.clear();
}


void PlaneSandwich::ControlPanel(const TGWindow* p, const TGWindow* main) {
  if (!fControlPanel)
    fControlPanel = new PlaneSandwichPanel(p, main, 100, 100, this);
}


void PlaneSandwich::put_in_hist(TH2F* h2, bool subtract_offs,
                                const ChannelInfo& ci)
{
  if (!ci.digit_sadc)
    return;
  const vector<CS::uint16> vsamples = ci.digit_sadc->GetSamples();
  float offs = subtract_offs ? ci.offs : 0;
  double nentries = h2->GetEntries();
  for (size_t i=0; i < vsamples.size(); i++)
    h2->Fill(i, vsamples[i] - offs);
  h2->SetEntries(nentries + 1);
}


void PlaneSandwich::StoreDigitSADC(const CS::ChipSADC::Digit* d_sadc)
{
  digits_sadc.push_back(d_sadc);
  //if( debug )
  //d_sadc->Print();
        
  if( d_sadc->GetY()!=0 )
    printf("PlaneSandwich::StoreDigit(): y=%d: it is ignored!!\n",d_sadc->GetY());

  int ch = d_sadc->GetX();
  if( ch<1 || ch>(int)channels_total )
    {
      printf("PlaneSandwich::StoreDigit(): bad channel number %d: ignored!!\n",ch);
      return;
    }

  h1_adc_chan->Fill(ch);

  const size_t window_width = 8;  // First window_width samples are used to determine the offset.
  double sum=0;
  double sum2 = 0;
  float integral = 0;
  CS::uint16 ampl_max = 0;
  CS::uint16 ampl_min = 1023;
  const vector<CS::uint16>& vsamples = d_sadc->GetSamples();
  for( size_t i=0; i<vsamples.size(); i++ )
    {
      // Calculate average + fluctuation only for first few samples
      if( i < window_width )
        {
          sum += vsamples[i];
          sum2 += vsamples[i]*vsamples[i];
        }

      integral += vsamples[i];

      if( ampl_max <= vsamples[i] )
        ampl_max = vsamples[i];
      if( ampl_min >= vsamples[i] )
        ampl_min = vsamples[i];
    }

  float offset = sum/window_width;
  float rms = sum2/window_width - offset*offset;

  integral -= vsamples.size()*offset;

  ChannelInfo& chi = channel_infos[ch];
  assert (!chi.digit_sadc);

  chi.min = ampl_min;
  chi.max = ampl_max;
  chi.offs = offset;
  chi.rms = rms;
  chi.integral = integral;
  chi.pileup = rms >= 4;
  chi.digit_sadc = d_sadc;

  vh1_adc_ampl_max_offs[ch-1]->Fill(ampl_max - offset);
  if (fExpertHistos && !chi.pileup) {
    vh1_adc_ampl_max[ch-1]->Fill(ampl_max);
    vh1_adc_integrals[ch-1]->Fill(integral);
    put_in_hist(vh2_adc_samples[ch-1], false, chi);
    put_in_hist(vh2_adc_samples_offs[ch-1], true, chi);
  }
}


void PlaneSandwich::StoreDigitF1(const CS::ChipF1::Digit* d_f1)
{
  digits_f1.push_back(d_f1);
  assert(d_f1->GetChannel()>=1 && d_f1->GetChannel()<=(int)channels_total);
  //printf("chan: %d %g\n",d_f1->GetChannel(),d_f1->GetTimeDecoded());

  int ch = d_f1->GetChannel();
  channel_infos[ch].digits_f1.push_back(d_f1);

  h1_tdc_chan->Fill(d_f1->GetChannel());
  h1_tdc_time->Fill(d_f1->GetTimeDecoded());
  h2_tdc_chan_time->Fill(d_f1->GetChannel(),d_f1->GetTimeDecoded());
}


void PlaneSandwich::StoreDigit(CS::Chip::Digit* digit)
{
    //printf("PlaneSandwich::StoreDigit(): call!\n");
    //digit->Print();
    
    const CS::ChipSADC::Digit* d_sadc = dynamic_cast<const CS::ChipSADC::Digit*>(digit);
    if( d_sadc!=NULL )
    {
      StoreDigitSADC(d_sadc);
      return;
    }

    
    const CS::ChipF1::Digit* d_f1 = dynamic_cast<const CS::ChipF1::Digit*>(digit);
    if( d_f1!=NULL )
    {
      StoreDigitF1(d_f1);
      return;
    }
}


void PlaneSandwich::EndEvent(const CS::DaqEvent &event)
{
  //vector<const ChipSADC::Digit*>::const_iterator sadc;
  //vector<const ChipF1::Digit*>::const_iterator f1;
  
  //printf("PlaneSandwich::EndEvent(): call!\n");
  //printf("%d %d\n",digits_sadc.size(),digits_f1.size());
  if (fExpertHistos) {
    if (digits_f1.size() > 0 || digits_sadc.size() > 0)
      h2_tdc_vs_adc->Fill(digits_f1.size(), digits_sadc.size());

    h2_tdc_vs_adc->ProfileY(name_tdc_vs_adc.c_str());
  }

  if (digits_sadc.size() > 0)
    {
      h1_adc_multiplicity->Fill(digits_sadc.size());
      if (digits_sadc.size() >= 13)
        puts("More ADC hits in sandwich than channels: impossible.\n");


      if (fExpertHistos) {
	if (digits_sadc.size() == 12) {
	  for (size_t i = 0; i < digits_f1.size(); i++) {
	    h2_tdc_chan_time_12->Fill(digits_f1[i]->GetChannel(),
				      digits_f1[i]->GetTimeDecoded());
	  }
	} else {
          for (map<int, ChannelInfo>::const_iterator i = channel_infos.begin();
               i != channel_infos.end(); i++) {
            int ch = i->first;
            const ChannelInfo& chi = i->second;
            if(chi.digit_sadc && !chi.pileup) {
              put_in_hist(vh2_adc_samples_offs_not_12[ch-1], true, chi);
            }
	  }
	}
      }
    }

  // Sort by arrival time.
  sort(digits_f1.begin(), digits_f1.end(), earlier());
  for(size_t i=1; i <= channels_total; i++) {
    sort(channel_infos[i].digits_f1.begin(),
         channel_infos[i].digits_f1.end(), earlier());
  }

  if (digits_f1.size() > 0 && digits_sadc.size() < 12) {
    int ch = digits_f1[0]->GetChannel();
    double t = digits_f1[0]->GetTimeDecoded();
    h2_tdc_chan_time_first->Fill(ch, t);

    if (fExpertHistos) {
      // Look at 'echos'.
      bool found = false;
      for(vector<const CS::ChipF1::Digit*>::const_iterator f1
            = channel_infos[5].digits_f1.begin();
          f1 != channel_infos[5].digits_f1.end(); f1++) {
	t = (*f1)->GetTimeDecoded();
	if (t > -888 && t < -874) {
	  h2_tdc_chan_time_echo->Fill(5, t);
	  found = true;
	  break;
	}
      }

      if (found) {
	vector<const CS::ChipF1::Digit*> digits_in_window;
	for (size_t i = 0; i < digits_f1.size(); i++) {
	  int ch = digits_f1[i]->GetChannel();
	  double t = digits_f1[i]->GetTimeDecoded();
	  if (t < -930 || t > -900)
	    continue;
	  h2_tdc_chan_time_echo->Fill(ch, t);
	  //digits_in_window.push_back(digits_f1[i]);
	}
#if 0
	if (digits_in_window.size() == 1) {
	  int ch = digits_in_window[0]->GetChannel();
	  double t = digits_in_window[0]->GetTimeDecoded();
	  h2_tdc_chan_time_echo->Fill(ch, t);
	}
#endif
      }
    }
  }

  if (digits_f1.size() > 1) {
    double first[channels_total];
    bool found_first[channels_total];

    memset(found_first, 0, sizeof(found_first));
    for(size_t i = 1; i < digits_f1.size(); i++) {
      int ch = digits_f1[i]->GetChannel();
      if (!found_first[ch-1]) {
        found_first[ch-1] = true;
        first[ch-1] = digits_f1[i]->GetTimeDecoded();
        h1_tdc_chan_hit->Fill(ch);
      } else {
	if (fExpertHistos) {
	  double current = digits_f1[i]->GetTimeDecoded();
	  vh1_tdc_time_diff[ch-1]->Fill(current - first[ch-1]);
	}
      }
    }
  }

  // Compare number of events where both adc and tdc fired to number of events
  // where either fired.
  double nentries = h1_adc_tdc_comparison->GetEntries();
  for(size_t ch = 1; ch <= channels_total; ch++) {
    bool tdc_hit = channel_infos[ch].digits_f1.size() > 0;
    bool adc_hit = channel_infos[ch].digit_sadc != 0;
    if (tdc_hit || adc_hit) {
      h1_adc_or_tdc->Fill(ch);
      if (tdc_hit && adc_hit)
        h1_adc_and_tdc->Fill(ch);
      h1_adc_tdc_comparison
        ->SetBinContent(ch,
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
	continue;  // Ignore pileup.

      float offset = channel_infos[ch].offs;
      CS::uint16 ampl_max = channel_infos[ch].max;
      CS::uint16 ampl_min = channel_infos[ch].min;

      float time = f1->GetTimeDecoded();

      if (tdc_hit == 1) {
	vh2_ampl_vs_time[ch-1]->Fill(ampl_max - offset, time);
	vh2_ampl_vs_time[ch-1]->Fill(ampl_min - offset, time);
      } else {
	vh2_ampl_vs_time_2[ch-1]->Fill(ampl_max - offset, time);
	vh2_ampl_vs_time_2[ch-1]->Fill(ampl_min - offset, time);
      }
    }

    // See if there was a Mip in channel 7 and no pileup, if yes
    // plot scope pictures where available.
    if (!channel_infos[7].pileup
        && channel_infos[7].max - channel_infos[7].offs >= 20
        && channel_infos[7].max - channel_infos[7].offs <= 35) {
      for (int ch = 1; ch <= 12; ch++) {
        const ChannelInfo& chi = channel_infos[ch];
        if (chi.digit_sadc) {
          const vector<CS::uint16>& vsamples
            = chi.digit_sadc->GetSamples();
          for (size_t j = 0; j < vsamples.size(); j++) {
            vh2_adc_mip_company[ch-1]
              ->Fill(j, vsamples[j] - chi.offs);
          }
        }
      }
    }

  }

  Plane::EndEvent(event);
}
