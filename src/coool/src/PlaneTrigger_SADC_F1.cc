#include "PlaneTrigger_SADC_F1.h"

ClassImp(PlaneTrigger_SADC_F1)

const float PlaneTrigger_SADC_F1::fF1_TICK   = 128.e-9; // in ms
const int PlaneTrigger_SADC_F1::fMAX_MULT    =   8;
const int PlaneTrigger_SADC_F1::fRATE_UPDATE = 100;

PlaneTrigger_SADC_F1::PlaneTrigger_SADC_F1(const char *detname, int nAdcChan, int nTdcChan, int center, int width) : 
    Plane(detname), fNAdcChan(nAdcChan), fNTdcChan(nTdcChan), fCenter(center), fWidth(width), inTime(false) {
    max_sadc = new float[fNAdcChan];
}

PlaneTrigger_SADC_F1::~PlaneTrigger_SADC_F1(void) {
    delete [] max_sadc;
}

void PlaneTrigger_SADC_F1::Init(TTree* tree) {
    //booking histograms

    // first some TDC stuff

    // TDC multiplicity
    std::string name = fName + "_hits";
    std::string title = fName + " TDC Multiplicity";
    std::string tdc_mult_name     = name;
    std::string tdc_mult_leavlist = tdc_mult_name + "/I";
    h1_tdc_mult = new TH1F_Ref(name.c_str(), title.c_str(), 21, -0.5, 20.5, fRateCounter);
    if (fNTdcChan > 0)
        AddHistogram(h1_tdc_mult);

    // TDC channels
    fVch = AddVariable("_ch",fNTdcChan,0,fNTdcChan,fNTdcChan*fMAX_MULT);
    name = fName + fVch->GetName();
    title = fName + " TDC channels";
    std::string tdc_ch_name     = name;
    std::string tdc_ch_leavlist = tdc_ch_name + "[" + tdc_mult_name + "]/F";
    h1_tdc_ch = new TH1F_Ref(name.c_str(), title.c_str(), fVch->GetNbins(), fVch->GetMin(), fVch->GetMax(),fRateCounter);
    if (fNTdcChan > 0)
        AddHistogram(h1_tdc_ch);

    // TDC timing in ns
    fVtNs = AddVariable("_t_ns",100,fCenter-fWidth,fCenter+fWidth,fNTdcChan*fMAX_MULT);
    name = fName + fVtNs->GetName();
    title = fName + " TDC Time Distribution (in ns)";
    std::string tdc_time_ns_name = name;
    std::string tdc_time_ns_leavlist = tdc_time_ns_name + "[" + tdc_mult_name + "]/F";
    h1_tdc_time_ns = new TH1F_Ref(name.c_str(), title.c_str(), fVtNs->GetNbins(), fVtNs->GetMin(), fVtNs->GetMax(), fRateCounter);
    if (fNTdcChan > 0)
        AddHistogram(h1_tdc_time_ns);

    // only SADC stuff to follow
    name = fName + "sadc_ch";
    title = fName + "SADC nbr of entries";
    std::string sadc_ch =name;
    h1_adc_ch = new TH1F_Ref(name.c_str(), title.c_str(),fNAdcChan,0.5,fNAdcChan+0.5,fRateCounter);
    h1_adc_ch->SetMinimum(0);
    if (fNAdcChan > 0)
        AddHistogram(h1_adc_ch);

    if( fExpertHistos && ( strncmp("HM01P2", fName.c_str(), 6) == 0 ) ) {

        name = fName + "AnalogSumMax";
        title = fName + "Software based Analog Sum Pulse Height";
        h1_Analog_Sum_max = new TH1F_Ref(name.c_str(), title.c_str(), 2048, 0, 2048, fRateCounter);
        if (fNAdcChan > 0)
            AddHistogram(h1_Analog_Sum_max);

        name = fName + "AnalogSum";
        title = fName + "Software based Analog Sum";
        h2_Analog_Sum = new TH2F(name.c_str(), title.c_str(), 32, 0, 32, 1024, 0, 1024);
        if (fNAdcChan > 0)
            AddHistogram(h2_Analog_Sum);


    }

    //single histograms per channel
    for(int ch=1; ch<=fNAdcChan; ch++)
    {
        char buff1[255], buff2[255];

        sprintf(buff1,"%s_sadc_ch%d",fName.c_str(),ch);
        sprintf(buff2,"%s_SADC Samples for channel %d",fName.c_str(),ch);
        h2_ADC_samples.push_back(new TH2F(buff1, buff2,32,0,32,1024,0,1024));
        if (fNAdcChan > 0)
            AddHistogram(h2_ADC_samples.back());
        sprintf(buff1,"%s_sadc_ch%d_tdc_cor",fName.c_str(),ch);
        sprintf(buff2,"%s_SADC Samples for channel %d, time corrected",fName.c_str(),ch);

        h2_ADC_samples_tdc_cor.push_back(new TH2F(buff1, buff2,32,0,32,1024,0,1024));
        if (fNAdcChan > 0)
            AddHistogram(h2_ADC_samples_tdc_cor.back());

        //expert histos
        if(fExpertHistos){
            sprintf(buff1,"%s_sadc_ch%d_wo_pu",fName.c_str(),ch);
            sprintf(buff2,"%s_SADC Samples for channel %d pileup corrected",fName.c_str(),ch);
            h2_ADC_samples_wo_pileup.push_back(new TH2F(buff1, buff2,32,0,32,1024,0,1024));
            if (fNAdcChan > 0)
                AddHistogram(h2_ADC_samples_wo_pileup.back());
            //sprintf(buff1,"%s_sadc_ch%d_wo_pu_tdc_cor",fName.c_str(),ch);
            //sprintf(buff2,"%s_SADC Samples for channel %d pileup and time corrected",fName.c_str(),ch);
            //h2_ADC_samples_wo_pileup_tdc_cor.push_back(new TH2F(buff1, buff2,32,0,32,1024,0,1024));
            //if (fNAdcChan > 0)
            //    AddHistogram(h2_ADC_samples_wo_pileup_tdc_cor.back());


            sprintf(buff1,"%s_sadc_ch%d_off",fName.c_str(),ch);
            sprintf(buff2,"%s_SADC Samples for channel %d, offset corrected",fName.c_str(),ch);
            h2_ADC_samples_off.push_back(new TH2F(buff1, buff2,32,0,32,1024,0,1024));
            if (fNAdcChan > 0)
                AddHistogram(h2_ADC_samples_off.back());
            sprintf(buff1,"%s_sadc_ch%d_off_tdc_cor",fName.c_str(),ch);
            sprintf(buff2,"%s_SADC Samples for channel %d, offset and time corrected",fName.c_str(),ch);
            h2_ADC_samples_off_tdc_cor.push_back(new TH2F(buff1, buff2,32,0,32,1024,0,1024));
            if (fNAdcChan > 0)
                AddHistogram(h2_ADC_samples_off_tdc_cor.back());
            sprintf(buff1,"%s_sadc_ch%d_off_wo_pu",fName.c_str(),ch);
            sprintf(buff2,"%s_SADC Samples for channel %d, offset and pileup corrected",fName.c_str(),ch);
            h2_ADC_samples_off_wo_pileup.push_back(new TH2F(buff1, buff2,32,0,32,1024,0,1024));
            if (fNAdcChan > 0)
                AddHistogram(h2_ADC_samples_off_wo_pileup.back());
            //sprintf(buff1,"%s_sadc_ch%d_off_wo_pu_tdc_cor",fName.c_str(),ch);
            //sprintf(buff2,"%s_SADC Samples for channel %d, offset, pileup and time corrected",fName.c_str(),ch);
            //h2_ADC_samples_off_wo_pileup_tdc_cor.push_back(new TH2F(buff1, buff2,32,0,32,1024,0,1024));
            //if (fNAdcChan > 0)
            //    AddHistogram(h2_ADC_samples_off_wo_pileup_tdc_cor.back());
        }

        sprintf(buff1,"%s_sadc_ch%d_max",fName.c_str(),ch);
        sprintf(buff2,"%s_SADC Maximum amplitude for channel %d",fName.c_str(),ch);
        h1_ADC_samples_max.push_back(new TH1F(buff1, buff2,1024,0,1024));
        if (fNAdcChan > 0)
            AddHistogram(h1_ADC_samples_max.back());
        sprintf(buff1,"%s_sadc_ch%d_max_tdc_cor",fName.c_str(),ch);
        sprintf(buff2,"%s_SADC Maximum amplitude for channel %d_time corrected",fName.c_str(),ch);
        h1_ADC_samples_max_tdc_cor.push_back(new TH1F(buff1, buff2,1024,0,1024));
        if (fNAdcChan > 0)
            AddHistogram(h1_ADC_samples_max_tdc_cor.back());

        //expert histos
        if(fExpertHistos){
            sprintf(buff1,"%s_sadc_ch%d_max_wo_pu",fName.c_str(),ch);
            sprintf(buff2,"%s_SADC Maximum amplitude for channel, pileup corrected %d",fName.c_str(),ch);
            h1_ADC_samples_max_wo_pileup.push_back(new TH1F(buff1, buff2,1024,0,1024));
            if (fNAdcChan > 0)
                AddHistogram(h1_ADC_samples_max_wo_pileup.back());
            //sprintf(buff1,"%s_sadc_ch%d_max_wo_pu_tdc_cor",fName.c_str(),ch);
            //sprintf(buff2,"%s_SADC Maximum amplitude for channel, pileup and time corrected %d",fName.c_str(),ch);
            //h1_ADC_samples_max_wo_pileup_tdc_cor.push_back(new TH1F(buff1, buff2,1024,0,1024));
            //if (fNAdcChan > 0)
            //    AddHistogram(h1_ADC_samples_max_wo_pileup_tdc_cor.back());

            sprintf(buff1,"%s_sadc_ch%d_max_off",fName.c_str(),ch);
            sprintf(buff2,"%s_SADC Maximum amplitude for channel %d, offset corrected",fName.c_str(),ch);
            h1_ADC_samples_max_off.push_back(new TH1F(buff1, buff2,1024,0,1024));
            if (fNAdcChan > 0)
                AddHistogram(h1_ADC_samples_max_off.back());
            sprintf(buff1,"%s_sadc_ch%d_max_off_tdc_cor",fName.c_str(),ch);
            sprintf(buff2,"%s_SADC Maximum amplitude for channel %d, offset and time corrected",fName.c_str(),ch);
            h1_ADC_samples_max_off_tdc_cor.push_back(new TH1F(buff1, buff2,1024,0,1024));
            if (fNAdcChan > 0)
                AddHistogram(h1_ADC_samples_max_off_tdc_cor.back());
            sprintf(buff1,"%s_sadc_ch%d_max_off_wo_pu",fName.c_str(),ch);
            sprintf(buff2,"%s_SADC Maximum amplitude for channel %d, offset and pileup corrected",fName.c_str(),ch);
            h1_ADC_samples_max_off_wo_pileup.push_back(new TH1F(buff1, buff2,1024,0,1024));
            if (fNAdcChan > 0)
                AddHistogram(h1_ADC_samples_max_off_wo_pileup.back());
            //sprintf(buff1,"%s_sadc_ch%d_max_off_wo_pu_tdc_cor",fName.c_str(),ch);
            //sprintf(buff2,"%s_SADC Maximum amplitude for channel %d, offset, pileup and time corrected",fName.c_str(),ch);
            //h1_ADC_samples_max_off_wo_pileup_tdc_cor.push_back(new TH1F(buff1, buff2,1024,0,1024));
            //if (fNAdcChan > 0)
            //    AddHistogram(h1_ADC_samples_max_off_wo_pileup_tdc_cor.back());
        }

        sprintf(buff1,"%s_sadc_ch%d_integral",fName.c_str(),ch);
        sprintf(buff2,"%s_SADC Amplitude integral for channel %d",fName.c_str(),ch);
        h1_ADC_samples_integral.push_back(new TH1F(buff1, buff2,50,0,2000));
        if (fNAdcChan > 0)
            AddHistogram(h1_ADC_samples_integral.back());

        //expert histos
        if(fExpertHistos){
            sprintf(buff1,"%s_sadc_ch%d_integral_wo_pu",fName.c_str(),ch);
            sprintf(buff2,"%s_SADC Amplitude integral for channel %d, pileup corrected",fName.c_str(),ch);
            h1_ADC_samples_integral_wo_pileup.push_back(new TH1F(buff1, buff2,50,0,2000));
            if (fNAdcChan > 0)
                AddHistogram(h1_ADC_samples_integral_wo_pileup.back());

            sprintf(buff1, "%s_sadc_ch%d_diff_max_min", fName.c_str(), ch);
            sprintf(buff2, "%s_SADC Difference between maximal and minimal sample in one event", fName.c_str());
            h1_ADC_samples_diff.push_back(new TH1F(buff1, buff2, 200, 0., 200.));
            if (fNAdcChan > 0)
                AddHistogram(h1_ADC_samples_diff.back());
        }
    }

    // rest of the TDC stuff

    // TDC timing
    fVt = AddVariable("_t",100,fCenter-fWidth,fCenter+fWidth,fNTdcChan*fMAX_MULT);
    name = fName + fVt->GetName();
    title = fName + " TDC Time Distribution";
    std::string tdc_time_name = name;
    std::string tdc_time_leavlist = tdc_time_name + "[" + tdc_mult_name + "]/F";
    h1_tdc_time = new TH1F_Ref(name.c_str(), title.c_str(), fVt->GetNbins(), fVt->GetMin(), fVt->GetMax(), fRateCounter);
    if (fNTdcChan > 0)
        AddHistogram(h1_tdc_time);

    // TDC timing versus channel
    name = fName + "_tVSch";
    title = fName + " TDC Time versus channel";
    h2_tdc_t_vs_ch = new TH2F(name.c_str(), title.c_str(), fVch->GetNbins(), fVch->GetMin(), fVch->GetMax(),
            fVt->GetNbins(),  fVt->GetMin(),  fVt->GetMax());
    h2_tdc_t_vs_ch->SetOption("COL");
    if (fNTdcChan > 0)
        AddHistogram(h2_tdc_t_vs_ch);

    // TDC on-trigger timing
    fVtOnTrig = AddVariable("_t_on_trigger",100,fCenter-fWidth,fCenter+fWidth,fNTdcChan*fMAX_MULT);
    name = fName + "_t_on_trig";
    title = fName + " TDC Time on trigger";
    h1_tdc_time_on_trig = new TH1F_Ref(name.c_str(), title.c_str(), fVtOnTrig->GetNbins(), fVtOnTrig->GetMin(), fVtOnTrig->GetMax(), fRateCounter);
    if (fNTdcChan > 0)
        AddHistogram(h1_tdc_time_on_trig);

    // TDC off-trigger timing
    fVtOffTrig = AddVariable("_t_off_trigger",100,fCenter-fWidth,fCenter+fWidth,fNTdcChan*fMAX_MULT);
    name = fName + "_t_off_trig";
    title = fName + " TDC Time off trigger";
    h1_tdc_time_off_trig = new TH1F(name.c_str(), title.c_str(), fVtOffTrig->GetNbins(), fVtOffTrig->GetMin(), fVtOffTrig->GetMax());
    if (fNTdcChan > 0)
        AddHistogram(h1_tdc_time_off_trig);

    // TDC on trigger channels
    name = fName + "_ch_on_trig";
    title = fName + " TDC Channels on trigger";
    h1_tdc_ch_on_trig = new TH1F_Ref(name.c_str(), title.c_str(), fVch->GetNbins(), fVch->GetMin(), fVch->GetMax(), fRateCounter);
    if (fNTdcChan > 0)
        AddHistogram(h1_tdc_ch_on_trig);

    // TDC off trigger channels
    name = fName + "_ch_off_trig";
    title = fName + " TDC Channels off trigger";
    h1_tdc_ch_off_trig = new TH1F(name.c_str(), title.c_str(), fVch->GetNbins(), fVch->GetMin(), fVch->GetMax());
    if (fNTdcChan > 0)
        AddHistogram(h1_tdc_ch_off_trig);

    // TDC corrected on trigger channels
    name = fName + "_ch_on_corr_trig";
    title = fName + " TDC Corrected on trigger channels";
    h1_tdc_cor_ch_on_trig = new TH1F(name.c_str(),title.c_str(), fVch->GetNbins(), fVch->GetMin(), fVch->GetMax());
    if (fNTdcChan > 0)
        AddHistogram(h1_tdc_cor_ch_on_trig);

    // TDC rates
    name = fName + "_rates";
    title = fName + " TDC Rates";
    h1_tdc_rates = new TH1F_Ref(name.c_str(), title.c_str(), fVch->GetNbins(), fVch->GetMin(), fVch->GetMax(), fRateCounter);
    h1_tdc_rates->SetYTitle("rates per channel (kHz)");
    h1_tdc_rates->SetTitleOffset(1.2,"Y");
    h1_tdc_rates->SetOption("hist");
    if (fNTdcChan > 0)
        AddHistogram(h1_tdc_rates);

    if(tree) {
        fIsInTree = true;

        if (fNTdcChan > 0) {
            tree->Branch(tdc_mult_name.c_str(), &tdc_mult, 
                         tdc_mult_leavlist.c_str(), 32000);
            tree->Branch(tdc_time_name.c_str(), fVt->GetValues(),
                         tdc_time_leavlist.c_str(), 32000);
            tree->Branch(tdc_time_ns_name.c_str(), fVtNs->GetValues(),
                         tdc_time_ns_leavlist.c_str(), 32000);
            tree->Branch(tdc_ch_name.c_str(), fVch->GetValues(),
                         tdc_ch_leavlist.c_str(), 32000);
        }

        for (int ch = 1; ch <= fNAdcChan; ch++) {

            char buff1[255], buff2[255];
            sprintf(buff1,"%s_sadc_ch%d_max",fName.c_str(),ch);
            sprintf(buff2,"%s_sadc_ch%d_max/F",fName.c_str(),ch);

            tree->Branch(buff1,&(max_sadc[ch-1]),buff2,32000);

        }

    }

}

void PlaneTrigger_SADC_F1::put_in_hist(TH2F *h2,bool subtract_offset ,const ChannelInfo& ci)
{
    if (!ci.digit_sadc)
        return;
    const vector<CS::uint16> vsamples = ci.digit_sadc->GetSamples();

    float offs = subtract_offset ? ci.offs : 0;

    double nb_entries = h2->GetEntries();
    for(size_t i = 0; i < vsamples.size(); i++) {
        if ( vsamples.size() != 32 ) continue;
        h2->Fill(i, vsamples[i] - offs);
    }
    h2->SetEntries(nb_entries +1);
}


void PlaneTrigger_SADC_F1::StoreDigitSADC(const CS::ChipSADC::Digit* d_sadc)
{
    digits_sadc.push_back(d_sadc);

    if(d_sadc->GetY() != 0)
        printf("Not in use y=%d", d_sadc->GetY());

    int ch = d_sadc->GetX();
    if(ch < 1 || ch > fNAdcChan)
    {
        printf("Bad Number %d", ch);
        return;
    }

    const size_t window_width = 5;  // First window_width samples are used to determine the offset.
    double sum=0;
    double sum2 = 0;
    float integral = 0;
    CS::uint16 ampl_max = 0;
    CS::uint16 ampl_min = 1023;
    const vector<CS::uint16>& vsamples = d_sadc->GetSamples();
    for( size_t i=0; i<vsamples.size(); i++ )
    {
        if (vsamples.size() != 32) continue;
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


    ChannelInfo& chi = channel_infos[ch-1];
    assert (!chi.digit_sadc);

    chi.min = ampl_min;
    chi.max = ampl_max;
    chi.offs = offset;
    chi.rms = rms;
    chi.integral = integral;
    chi.pileup = rms >= 6;
    chi.digit_sadc = d_sadc;

}


void PlaneTrigger_SADC_F1::StoreDigitF1(const CS::ChipF1::Digit* d_f1)
{

    assert(d_f1->GetChannel()>=0 && d_f1->GetChannel()<=fNTdcChan);

    digits_f1.push_back(d_f1);

}

void PlaneTrigger_SADC_F1::StoreDigit(CS::Chip::Digit* digit)
{
    const CS::ChipSADC::Digit* d_sadc = dynamic_cast<const CS::ChipSADC::Digit*>(digit);
    if (d_sadc != NULL)
    {
        StoreDigitSADC(d_sadc);
        return;
    }

    const CS::ChipF1::Digit* d_f1 = dynamic_cast<const CS::ChipF1::Digit*>(digit);
    if (d_f1 != NULL)
    {
        StoreDigitF1(d_f1);
        return;
    }
}

void PlaneTrigger_SADC_F1::EndEvent(const CS::DaqEvent &event)
{
    if (thr_flag) TThread::Lock();



    //Get time information

    double tdc_channel = 0;
    double time = 0;

    //tdc multiplicity
    tdc_mult = digits_f1.size();
    h1_tdc_mult->Fill(tdc_mult);

    for(unsigned int i = 0; i < digits_f1.size(); i++)
    {
        time = digits_f1[i]->GetTimeDecoded();
        double timeRaw = digits_f1[i]->GetTimeDecoded() / digits_f1[i]->GetTimeUnit();

        tdc_channel = digits_f1[i]->GetChannel();
        if (-960 < time && time < -945) inTime = true;

        //Fill TDC histos
        h1_tdc_time->Fill(timeRaw);
        h1_tdc_time_ns->Fill(time);
        h1_tdc_ch->Fill(tdc_channel);
        h2_tdc_t_vs_ch->Fill(tdc_channel, timeRaw);

        if( fVtOnTrig->Test(timeRaw) )  {
            h1_tdc_time_on_trig->Fill(timeRaw);       // on trigger time distribution
            h1_tdc_ch_on_trig  ->Fill(tdc_channel);   // on trigger profile
        }

        if( fVtOffTrig->Test(timeRaw) )  {
            h1_tdc_time_off_trig->Fill(timeRaw);      // off trigger time distribution
            h1_tdc_ch_off_trig  ->Fill(tdc_channel);  // off trigger profile
        }

        fVt->Store(timeRaw);
        fVtNs->Store(time);
        fVch->Store(tdc_channel);
    }

    if((fRateCounter%fRATE_UPDATE)==0) {
        // corrected on-trigger profile
        if( fVtOffTrig->GetMax() != fVtOffTrig->GetMin() ) {
            float scale = (float) ( fVtOnTrig->GetMax() - fVtOnTrig->GetMin() ) / ( fVtOffTrig->GetMax() - fVtOffTrig->GetMin() );
            for(register int i=1; i<=h1_tdc_ch_on_trig->GetNbinsX(); i++)
                h1_tdc_cor_ch_on_trig->SetBinContent(i, (float) h1_tdc_ch_on_trig->GetBinContent(i) - scale*h1_tdc_ch_off_trig->GetBinContent(i) );
        }

        // rates histogram
        float timewin = ( fVtOffTrig->GetMax() - fVtOffTrig->GetMin() )*fF1_TICK*fRateCounter;
        if( timewin > 0 )
            for(register int i=1; i<=h1_tdc_ch_on_trig->GetNbinsX(); i++) {
                h1_tdc_rates->SetBinContent(i,(float) h1_tdc_ch_off_trig->GetBinContent(i)/timewin);
                h1_tdc_rates->SetBinError(i,(float) h1_tdc_ch_off_trig->GetBinError(i)/timewin);
            }
    }

    //Get SADC informations

    for (int ch = 1; ch <= fNAdcChan; ch++)
    {
        const ChannelInfo& chi = channel_infos[ch-1];

        if (!chi.digit_sadc)
            continue;

        max_sadc[ch-1] = chi.max;

        //----------Drawing histos----------------------

        h1_adc_ch->Fill(ch);

        h1_ADC_samples_max[ch-1]->Fill(chi.max);
        h1_ADC_samples_integral[ch-1]->Fill(chi.integral);

        put_in_hist(h2_ADC_samples[ch-1], false,chi);


        if(fExpertHistos) {

            h1_ADC_samples_max_off[ch-1]->Fill(chi.max - chi.offs);
            put_in_hist(h2_ADC_samples_off[ch-1], true,chi);

            h1_ADC_samples_diff[ch-1]->Fill(chi.max - chi.min);
        }

        if(fExpertHistos && !chi.pileup) {

            h1_ADC_samples_max_wo_pileup[ch-1]->Fill(chi.max);
            h1_ADC_samples_max_off_wo_pileup[ch-1]->Fill(chi.max - chi.offs);
            h1_ADC_samples_integral_wo_pileup[ch-1]->Fill(chi.integral);

            put_in_hist(h2_ADC_samples_wo_pileup[ch-1], false,chi);
            put_in_hist(h2_ADC_samples_off_wo_pileup[ch-1], true,chi);

        }

        if(inTime) {

            h1_ADC_samples_max_tdc_cor[ch-1]->Fill(chi.max);
            //h1_ADC_samples_integral_wo_pileup[ch-1]->Fill(chi.integral);

            put_in_hist(h2_ADC_samples_tdc_cor[ch-1], false,chi);
        }

        if(inTime && fExpertHistos) {

            h1_ADC_samples_max_off_tdc_cor[ch-1]->Fill(chi.max - chi.offs);

        }
    }



    //------------Offline Analog Sum-------------------------------
    if( fExpertHistos && ( strncmp("HM01P2", fName.c_str(), 6) == 0 ) ) {
        if (channel_infos[1].digit_sadc && channel_infos[2].digit_sadc) {
            const vector<CS::uint16> vsample_1 = channel_infos[1].digit_sadc->GetSamples();
            const vector<CS::uint16> vsample_2 = channel_infos[2].digit_sadc->GetSamples();

            vector<CS::uint16> analog_sum;

            CS::uint16 sum_max = 0;

            if (vsample_1.size() == vsample_2.size()) {
                for(size_t i = 0; i< vsample_1.size(); i++) {
                    analog_sum.push_back(vsample_1[i] + vsample_2[i]);

                    if( sum_max <= analog_sum[i] )
                        sum_max = analog_sum[i];
                }

                for(size_t i = 0; i< analog_sum .size(); i++) {
                    if(inTime) 
                        h2_Analog_Sum->Fill(i,analog_sum[i]);
                }

                if(inTime) h1_Analog_Sum_max->Fill(sum_max);

                //max_sadc[0] = sum_max;

            }
        }
    }
    //-------------------------------------


    inTime = false;


    if (thr_flag) TThread::UnLock();
}

void PlaneTrigger_SADC_F1::Reset(void)
{
    Plane::Reset();

    for(int i = 0; i < fNAdcChan; i++) {
        max_sadc[i] = 0;
    }

    for (map<int, ChannelInfo>::iterator it = channel_infos.begin();
            it != channel_infos.end(); it++)
    {
        it->second.clear();
    }

    channel_infos.clear();

    digits_sadc.clear();
    digits_f1.clear();
}


void PlaneTrigger_SADC_F1::ControlPanel(const TGWindow* p, const TGWindow* main) {
    if (!fControlPanel) fControlPanel = new PlanePanel(p, main, 100, 100, this);
}



