#include "PlaneRPD_SADC.h"
#include "PlaneRPD_SADC_Panel.h"
#include "TProfile.h"
#include "ChipSADC.h"
#include <cmath>

ClassImp(PlaneRPD_SADC);
int Modulo(int value, int modulo);

PlaneRPD_SADC::PlaneRPD_SADC(const char *detname) : Plane(detname)
{
  for(unsigned short int iChannel=0; iChannel<N_CHANNEL_SADC; ++iChannel)
    {
      fDataSent[iChannel]=false;
      fCalculatedPedestal[iChannel]=false;
      fPedestal[iChannel]=0;
      fPedestalRelativeError[iChannel]=0;
      fAmplitudeMaximum[iChannel]=0;
      fAdc[iChannel]=0;
      for(unsigned short int iSample=0; iSample<N_SADC_SAMPLE; ++iSample)
	{
	  fAmplitude[iChannel][iSample]=0;
	}
    }
}


PlaneRPD_SADC::~PlaneRPD_SADC()
{
}


void PlaneRPD_SADC::Init(TTree* tree)
{
  TString histname = fName + "_NbOfData";
  TString histtitle = fName + " Number of channels which sent data";
  h1NumberOfSentData=new TH1F_Ref(histname,histtitle,N_CHANNEL_SADC+1,0,N_CHANNEL_SADC+1,fRateCounter);
  h1NumberOfSentData->GetXaxis()->SetTitle("Number of sent data");
  AddHistogram(h1NumberOfSentData);

  histname = fName + "_PedesVsChannel";
  histtitle = fName + " Pedestals spectrum vs channels";
  h2PedestalsVsChannel=new TH2F(histname,histtitle,N_CHANNEL_SADC,0,N_CHANNEL_SADC, 100, 0, 200);
  h2PedestalsVsChannel->SetOption("col");
  h2PedestalsVsChannel->GetXaxis()->SetTitle("Channel");
  h2PedestalsVsChannel->GetYaxis()->SetTitle("Calculated pedestals");
  AddHistogram(h2PedestalsVsChannel);

  histname = fName + "_MaxAmplVsChannel";
  histtitle = fName + " Maximum amplitude spectrum vs channel";
  h2AmplitudeMaximumVsChannel=new TH2F(histname,histtitle,N_CHANNEL_SADC,0,N_CHANNEL_SADC, 100, 0, 1000);
  h2AmplitudeMaximumVsChannel->SetOption("col");
  h2AmplitudeMaximumVsChannel->GetXaxis()->SetTitle("Channel");
  h2AmplitudeMaximumVsChannel->GetYaxis()->SetTitle("Maximum of amplitude");
  AddHistogram(h2AmplitudeMaximumVsChannel);

  histname = fName + "_AdcVsChannel";
  histtitle = fName + " Adc spectrum vs channel";
  h2AdcVsChannel=new TH2F(histname,histtitle,N_CHANNEL_SADC,0,N_CHANNEL_SADC, 100, 0, 4000);
  h2AdcVsChannel->SetOption("col");
  h2AdcVsChannel->GetXaxis()->SetTitle("Channel");
  h2AdcVsChannel->GetYaxis()->SetTitle("Adc");
  AddHistogram(h2AdcVsChannel);

  histname = fName + "_StartVsChannel";
  histtitle = fName + " Beginning Of Signal vs channel";
  h2BeginningOfSignalVsChannel=new TH2F(histname,histtitle,N_CHANNEL_SADC,0,N_CHANNEL_SADC, N_SADC_SAMPLE, 0, N_SADC_SAMPLE);
  h2BeginningOfSignalVsChannel->SetOption("col");
  h2BeginningOfSignalVsChannel->GetXaxis()->SetTitle("Channel");
  h2BeginningOfSignalVsChannel->GetYaxis()->SetTitle("Beginning of signal");
  AddHistogram(h2BeginningOfSignalVsChannel);
  
  histname = fName + "_Channel";
  histtitle = fName + " Channel";
  h1Channel=new TH1F(histname,histtitle,N_CHANNEL_SADC,0,N_CHANNEL_SADC);
  h1Channel->GetXaxis()->SetTitle("Channel");
  AddHistogram(h1Channel);
  
  histname = fName + "_EnergyVsChannel";
  histtitle = fName + " Deposited energy vs channel";
  h2EnergyVsChannel=new TH2F(histname,histtitle,N_CHANNEL_SADC/2,0,N_CHANNEL_SADC/2, 100,0,4000);
  h2EnergyVsChannel->SetOption("col");
  h2EnergyVsChannel->GetXaxis()->SetTitle("Channel");
  h2EnergyVsChannel->GetYaxis()->SetTitle("Deposited energy");
  AddHistogram(h2EnergyVsChannel);
  
  histname = fName + "_WidthVsChannel";
  histtitle = fName + " Number of sample higher than Threshold vs channel";
  h2NumberOfSampleAboveThreshold=new TH2F(histname,histtitle,N_CHANNEL_SADC+1,0,N_CHANNEL_SADC+1, N_SADC_SAMPLE,0,N_SADC_SAMPLE);
  h2NumberOfSampleAboveThreshold->SetOption("col");
  h2NumberOfSampleAboveThreshold->GetXaxis()->SetTitle("Channel");
  h2NumberOfSampleAboveThreshold->GetYaxis()->SetTitle("Number of sample higher than Threshold");
  AddHistogram(h2NumberOfSampleAboveThreshold);
  
//   h2AdcUpVsAdcDown = new TH2F*[N_CHANNEL_SADC/2];
//   for(int iChannel=0; iChannel<N_CHANNEL_SADC_A/2; ++iChannel)
//     {
//       histname = fName + "_aupado"; histname+=iChannel;
//       histtitle = fName + " corr Adc up vs Adc down for A"; histtitle+=iChannel;
//       h2AdcUpVsAdcDown[iChannel]=new TH2F(histname,histtitle,400,0,4000,400,0,4000);
//       //       h2AdcUpVsAdcDown[iChannel]->SetOption("col");
//       h2AdcUpVsAdcDown[iChannel]->GetXaxis()->SetTitle("Adc up");
//       h2AdcUpVsAdcDown[iChannel]->GetYaxis()->SetTitle("Adc down");
//       AddHistogram(h2AdcUpVsAdcDown[iChannel]);
//     }
//   for(int iChannel=N_CHANNEL_SADC_A/2; iChannel<N_CHANNEL_SADC/2; ++iChannel)
//     {
//       histname = fName + "_aupado"; histname+=iChannel;
//       histtitle = fName + " corr Adc up vs Adc down for B"; histtitle+=iChannel-N_CHANNEL_SADC_A/2;
//       h2AdcUpVsAdcDown[iChannel]=new TH2F(histname,histtitle,400,0,4000,400,0,4000);
//       //       h2AdcUpVsAdcDown[iChannel]->SetOption("col");
//       h2AdcUpVsAdcDown[iChannel]->GetXaxis()->SetTitle("Adc up");
//       h2AdcUpVsAdcDown[iChannel]->GetYaxis()->SetTitle("Adc down");
//       AddHistogram(h2AdcUpVsAdcDown[iChannel]);
//     }

  //   h2SignalShape = new TH2F*[N_CHANNEL_SADC];
  //   for(unsigned int iChannel=0; iChannel<N_CHANNEL_SADC_A/2; ++iChannel)
  //     {
  //       histname = fName + "_ssau";
  //       histname += iChannel;
  //       histtitle = fName + " Signal shape for A up ";
  //       histtitle += iChannel;
  //       h2SignalShape[iChannel] = new TH2F(histname,histtitle, N_SADC_SAMPLE, 0, N_SADC_SAMPLE, 100, -100, 800);
  //       h2SignalShape[iChannel]->SetOption("col");
  //       AddHistogram(h2SignalShape[iChannel]);
  //     }
  //   for(unsigned int iChannel=0; iChannel<N_CHANNEL_SADC_A/2; ++iChannel)
  //     {
  //       histname = fName + "_ssad";
  //       histname += Modulo(iChannel,N_CHANNEL_SADC_A/2);
  //       histtitle = fName + " Signal shape for A down ";
  //       histtitle += Modulo(iChannel,N_CHANNEL_SADC_A/2);
  //       h2SignalShape[N_CHANNEL_SADC_A/2+iChannel] = new TH2F(histname,histtitle, N_SADC_SAMPLE, 0, N_SADC_SAMPLE, 100, -100, 800);
  //       h2SignalShape[N_CHANNEL_SADC_A/2+iChannel]->SetOption("col");
  //       AddHistogram(h2SignalShape[N_CHANNEL_SADC_A/2+iChannel]);
  //     }
  //   for(unsigned int iChannel=0; iChannel<N_CHANNEL_SADC_B/2; ++iChannel)
  //     {
  //       histname = fName + "_ssbu";
  //       histname += Modulo(iChannel,N_CHANNEL_SADC_B/2);
  //       histtitle = fName + " Signal shape for B up ";
  //       histtitle += Modulo(iChannel,N_CHANNEL_SADC_B/2);
  //       h2SignalShape[N_CHANNEL_SADC_A+iChannel] = new TH2F(histname,histtitle, N_SADC_SAMPLE, 0, N_SADC_SAMPLE, 100, -100, 800);
  //       h2SignalShape[N_CHANNEL_SADC_A+iChannel]->SetOption("col");
  //       AddHistogram(h2SignalShape[N_CHANNEL_SADC_A+iChannel]);
  //     }
  //   for(unsigned int iChannel=0; iChannel<N_CHANNEL_SADC_B/2; ++iChannel)
  //     {
  //       histname = fName + "_ssbd";
  //       histname += Modulo(iChannel,N_CHANNEL_SADC_B/2);
  //       histtitle = fName + " Signal shape for B down ";
  //       histtitle += Modulo(iChannel,N_CHANNEL_SADC_B/2);
  //       h2SignalShape[N_CHANNEL_SADC_A+N_CHANNEL_SADC_B/2+iChannel] = new TH2F(histname,histtitle, N_SADC_SAMPLE, 0, N_SADC_SAMPLE, 100, -100, 800);
  //       h2SignalShape[N_CHANNEL_SADC_A+N_CHANNEL_SADC_B/2+iChannel]->SetOption("col");
  //       AddHistogram(h2SignalShape[N_CHANNEL_SADC_A+N_CHANNEL_SADC_B/2+iChannel]);
  //     }

  h1SignalShape = new TH1F*[N_CHANNEL_SADC];
  for(unsigned int iChannel=0; iChannel<N_CHANNEL_SADC; ++iChannel)
    {
      histname = fName + "_ss";
      histname += iChannel;
      histtitle = fName + " Signal shape for channel ";
      histtitle += iChannel;
      h1SignalShape[iChannel] = new TH1F(histname,histtitle, N_SADC_SAMPLE, 0, N_SADC_SAMPLE);
      AddHistogram(h1SignalShape[iChannel]);
    }

  if(tree)
    {
      fIsInTree = true;
      TString branchName;
      TString planeName = (TString)fName;
      branchName="fAmplitudeMaximum[";branchName+=N_CHANNEL_SADC;branchName+="]/D";
      tree->Branch("rawAmplitudeMaximum"+planeName,fAmplitudeMaximum,branchName);
      branchName="fAmplitude[";branchName+=N_CHANNEL_SADC;branchName+="]";
      branchName+="[";branchName+=N_SADC_SAMPLE;branchName+="]/D";
      tree->Branch("rawAmplitude"+planeName,fAmplitude,branchName);
      branchName="fDataSent[";branchName+=N_CHANNEL_SADC;branchName+="]/B";
      tree->Branch("rawDataSent"+planeName,fDataSent,branchName);
      branchName="fCalculatedPedestal[";branchName+=N_CHANNEL_SADC;branchName+="]/B";
      tree->Branch("rawCalculatedPedestal"+planeName,fCalculatedPedestal,branchName);
      branchName="fPedestal[";branchName+=N_CHANNEL_SADC;branchName+="]/D";
      tree->Branch("rawPedestal"+planeName,fPedestal,branchName);
      branchName="fBeginningOfSignal[";branchName+=N_CHANNEL_SADC;branchName+="]/I";
      tree->Branch("rawBeginningOfSignal"+planeName,fBeginningOfSignal,branchName);

      tree->Branch("rawEventNumberInRun"+planeName,&fEventNumberInRun,"fEventNumberInRun/I");

    }

  return;
}


void PlaneRPD_SADC::StoreDigit(CS::Chip::Digit* digit)
{
  std::vector<float> data=digit->GetNtupleData();
  if(data.size()<3) return;
  const CS::ChipSADC::Digit* sadcdig = dynamic_cast<const CS::ChipSADC::Digit*>(digit);
  if(sadcdig != NULL )
    {
//JM 15/10 : uncomment the line below to get amplitudes
//      sadcdig->Print(cerr);
      uint16 channel = sadcdig->GetX();
      if(channel>N_CHANNEL_SADC) {cerr<<"channel = sadcdig->GetX()>N_CHANNEL_SADC, end of PlaneRPD_SADC::StoreDigit(CS::Chip::Digit* digit)"; return;}
      fDataSent[channel]=true;
//      uint16 chip = sadcdig->GetChip();
      uint16 chipChannel = sadcdig->GetChannel();
      for(unsigned short int iSample=0; iSample<data.size(); ++iSample)
	{
// eb 	  fAmplitude[channel][iSample]=data[iSample+2];
//	  printf("%f\n",);
	  fAmplitude[channel][iSample]=data[iSample+2];
	}
    }
  return;
}


void PlaneRPD_SADC::EndEvent(const CS::DaqEvent &event)
{
  LookForPedestals();
  CorrectPedestal();
  CalculateAdc();
  LookForMaximumAmplitude();
  CalculateNumberOfSampleAboveThreshold();
  fEventTimeSecond = (double)event.GetTime().first;
  fEventTimeMicroSecond = (double)event.GetTime().second;
  fEventRunNumber = (int)event.GetRunNumber();
  fEventNumberInRun = (int)event.GetEventNumberInRun();
  
  unsigned short int numberOfSentData=0;
  for(unsigned short int iChannel=0; iChannel<N_CHANNEL_SADC; ++iChannel)
    {
//       if(fDataSent[iChannel] && fNumberOfSampleAboveThreshold[iChannel]<5) fDataSent[iChannel]=false; // this line can be use to remove background but reject <1% of events
      if(fDataSent[iChannel])++numberOfSentData;

      if(fDataSent[iChannel])
	{
	  h1Channel->Fill(iChannel);
	  h2AmplitudeMaximumVsChannel->Fill(iChannel, fAmplitudeMaximum[iChannel]);
	  h2AdcVsChannel->Fill(iChannel, fAdc[iChannel]);
	  h2BeginningOfSignalVsChannel->Fill(iChannel, fBeginningOfSignal[iChannel]);
	  h2NumberOfSampleAboveThreshold->Fill(iChannel,fNumberOfSampleAboveThreshold[iChannel]);
	  if(fCalculatedPedestal[iChannel])
	    {
	      h2PedestalsVsChannel->Fill(iChannel, fPedestal[iChannel]);
	    }
	}
    }

  h1NumberOfSentData->Fill(numberOfSentData);

  for(unsigned short int iChannel=0; iChannel<N_CHANNEL_SADC_A/2; ++iChannel)
    {
      if(fDataSent[iChannel] && fDataSent[iChannel+N_CHANNEL_SADC_A/2])
	{
// 	  h2AdcUpVsAdcDown[iChannel]->Fill(fAdc[iChannel],fAdc[iChannel+N_CHANNEL_SADC_A/2]);
	  h2EnergyVsChannel->Fill(iChannel,sqrt(fAdc[iChannel]*fAdc[iChannel+N_CHANNEL_SADC_A/2]));
	}
    }
  for(unsigned short int iChannel=N_CHANNEL_SADC_A/2; iChannel<N_CHANNEL_SADC_A/2+N_CHANNEL_SADC_B/2; ++iChannel)
    {
// modified by E.B. on June 10,2008
// modified by GJ on 03_07_08
      int iUp = iChannel+N_CHANNEL_SADC_A/2;
      int iDown = iChannel+N_CHANNEL_SADC_A/2+N_CHANNEL_SADC_B/2;
      if(fDataSent[iUp]&&fDataSent[iDown])
	{
// 	  h2AdcUpVsAdcDown[iChannel]->Fill(fAdc[iUp],fAdc[iDown]);
	  h2EnergyVsChannel->Fill(iChannel,sqrt(fAdc[iUp]*fAdc[iDown]));
	}
    }
  for(unsigned short int iChannel=0; iChannel<N_CHANNEL_SADC_A/2; ++iChannel)
    {
      // rewritten completely by E.B. on June 13,2008
      int iAu=iChannel;
      int iAd=iChannel+N_CHANNEL_SADC_A/2;
      int iBu=2*iChannel+N_CHANNEL_SADC_A;
      int iBd=2*iChannel+N_CHANNEL_SADC_A+N_CHANNEL_SADC_B/2;
      int iBuPrevious=iBu-1;
      if (iChannel==0) iBuPrevious+=N_CHANNEL_SADC_B/2;
      int iBdPrevious=iBd-1;
      if (iChannel==0) iBdPrevious+=N_CHANNEL_SADC_B/2;
      int iBuNext=iBu+1;
      if (iChannel==11) iBuNext-=N_CHANNEL_SADC_B/2;
      int iBdNext=iBd+1;
      if (iChannel==11) iBdNext-=N_CHANNEL_SADC_B/2;
//       if(fDataSent[iAu] && fDataSent[iAd])
// 	{
// 	  if( fDataSent[iBu] && fDataSent[iBd] )
// 	    {
// 	      h2EnergyVsChannelAVsEnergyB[3*iChannel+1]->Fill(sqrt(fAdc[iAu]*fAdc[iAd]),sqrt(fAdc[iBu]*fAdc[iBd]));
// 	    }
// 	  if(fDataSent[iBuPrevious] && fDataSent[iBdPrevious])
// 	    {
// 	      h2EnergyVsChannelAVsEnergyB[3*iChannel]->Fill(sqrt(fAdc[iAu]*fAdc[iAd]),sqrt(fAdc[iBuPrevious]*fAdc[iBdPrevious]));
// 	    }
// 	  if(fDataSent[iBuNext] && fDataSent[iBdNext])
// 	    {
// 	      h2EnergyVsChannelAVsEnergyB[3*iChannel+2]->Fill(sqrt(fAdc[iAu]*fAdc[iAd]),sqrt(fAdc[iBuNext]*fAdc[iBdNext]));
// 	    }
// 	}
    }

  for(unsigned int iChannel=0; iChannel<N_CHANNEL_SADC; ++iChannel)
    {
      if(fDataSent[iChannel])
	for(unsigned int iSample=2; iSample<N_SADC_SAMPLE; ++iSample)
// 	  h2SignalShape[iChannel]->Fill(iSample, fAmplitude[iChannel][iSample]);
	  h1SignalShape[iChannel]->SetBinContent(iSample, fAmplitude[iChannel][iSample]);
      else
	for(unsigned int iSample=2; iSample<N_SADC_SAMPLE; ++iSample)
	  h1SignalShape[iChannel]->SetBinContent(iSample,0);
    }  
  
  return;
}


// void PlaneRPD_SADC::LookForPedestals()
// {
//   for(unsigned short int iChannel=0; iChannel<N_CHANNEL_SADC; ++iChannel)
//     {
//       if(!fDataSent[iChannel])
// 	{
// 	  continue;
// 	}
//       LookForBeginningOfSignal(iChannel);
// //       if(fBeginningOfSignal[iChannel]!=0)
//       if(fBeginningOfSignal[iChannel]>2)
// 	{
// 	  fCalculatedPedestal[iChannel]=true;
// 	  fPedestal[iChannel]=0;
// 	  for(unsigned short int iSample=2; iSample<fBeginningOfSignal[iChannel]; ++iSample)
// 	    {fPedestal[iChannel]+=fAmplitude[iChannel][iSample];}
// 	  fPedestal[iChannel]=fPedestal[iChannel]/fBeginningOfSignal[iChannel];
// 	  fPedestalRelativeError[iChannel]=1./sqrt(fBeginningOfSignal[iChannel]);
// 	}
//     }
//   return;
// }

void PlaneRPD_SADC::LookForPedestals()
{
  for(unsigned short int iChannel=0; iChannel<N_CHANNEL_SADC; ++iChannel)
    {
      if(!fDataSent[iChannel])
	{
	  continue;
	}
      LookForBeginningOfSignal(iChannel);
      if(1)
	{
	  fPedestal[iChannel]=0;
	  //      printf("pedestal of channel %i = %f\n",iChannel, fPedestal[iChannel]);
	  for(unsigned short int iSample=2; iSample<6; ++iSample)
	    {
	      //        printf("increment of pedestal %f\n", fAmplitude[iChannel][iSample]);
	      fPedestal[iChannel]+=fAmplitude[iChannel][iSample];
	    }
	  fPedestal[iChannel]=fPedestal[iChannel]/4;
	  fPedestalRelativeError[iChannel]=1./sqrt(4);
	  //      printf("pedestal of channel %i = %f\n",iChannel, fPedestal[iChannel]);
	  fCalculatedPedestal[iChannel]=true;
	}
    }
  return;
}

void PlaneRPD_SADC::CorrectPedestal()
{
  for(unsigned int iChannel=0; iChannel<N_CHANNEL_SADC; ++iChannel)
    {
      if(fCalculatedPedestal[iChannel])
	for(unsigned int iSample=0; iSample<N_SADC_SAMPLE; ++iSample)
	  fAmplitude[iChannel][iSample]-=fPedestal[iChannel];
    }
  return;
}


// void PlaneRPD_SADC::LookForBeginningOfSignal(unsigned short int channel)
// {
//   fBeginningOfSignal[channel]=0;
//   for(unsigned short int iSample=LookForSampleOfMaximumAmplitude(channel); iSample>0; --iSample)
//     if(fAmplitude[channel][iSample-1]==fAmplitude[channel][iSample]) {fBeginningOfSignal[channel] = iSample;return;}
// //   LookForSampleOfMaximumAmplitude(channel);
// //   fBeginningOfSignal[channel]=7;
//   return;
// }

void PlaneRPD_SADC::LookForBeginningOfSignal(unsigned short int channel)
{
  fBeginningOfSignal[channel]=0;
  for(unsigned short int iSample=LookForSampleOfMaximumAmplitude(channel); iSample>0; --iSample)
    if(fAmplitude[channel][iSample]>fAmplitude[channel][iSample-1]+3) {
      fBeginningOfSignal[channel] = iSample;
      return;
    }
//   LookForSampleOfMaximumAmplitude(channel);
//   fBeginningOfSignal[channel]=7;
  return;
}

unsigned short int PlaneRPD_SADC::LookForSampleOfMaximumAmplitude(unsigned short int channel)
{
  unsigned short int iMax=0;
  for(unsigned short int iSample=1; iSample<N_SADC_SAMPLE; ++iSample)
    if(fAmplitude[channel][iSample]>fAmplitude[channel][iMax]) iMax=iSample;
//   fAmplitudeMaximum[channel]=fAmplitude[channel][iMax];
  return iMax;
}

void PlaneRPD_SADC::LookForMaximumAmplitude()
{
  for(unsigned int iChannel=0; iChannel<N_CHANNEL_SADC; ++iChannel)
    if(fDataSent[iChannel])
      fAmplitudeMaximum[iChannel]=fAmplitude[iChannel][LookForSampleOfMaximumAmplitude(iChannel)];
  return;
}


void PlaneRPD_SADC::ControlPanel(const TGWindow* p, const TGWindow* main)
{
  if (!fControlPanel) fControlPanel = new PlaneRPD_SADC_Panel(p, main, 100, 100, this);
}

void PlaneRPD_SADC::MaxAmplitudeSpectrum()
{
  std::string hdname = fName + "_1ch";
  if(fCurChan_maxAmplitude<N_CHANNEL_SADC)
    {
      h2AmplitudeMaximumVsChannel->ProjectionY(hdname.c_str(),fCurChan_maxAmplitude+1,fCurChan_maxAmplitude+1,"")->DrawCopy();
    }
}

void PlaneRPD_SADC::AdcSpectrum()
{
  std::string hdname = fName + "_2ch";
  if(fCurChan_adc<N_CHANNEL_SADC)
    {
      h2AdcVsChannel->ProjectionY(hdname.c_str(),fCurChan_adc+1,fCurChan_adc+1,"")->DrawCopy();
    }
}

void PlaneRPD_SADC::PedestalsSpectrum()
{
  std::string hdname = fName + "_3ch";
  if(fCurChan_pedestals<N_CHANNEL_SADC)
    {
      h2PedestalsVsChannel->ProjectionY(hdname.c_str(),fCurChan_pedestals+1,fCurChan_pedestals+1,"")->DrawCopy();
    }
}

void PlaneRPD_SADC::BeginningOfSignalSpectrum()
{
  std::string hdname = fName + "_4ch";
  if(fCurChan_beginningOfSignal<N_CHANNEL_SADC)
    {
      h2BeginningOfSignalVsChannel->ProjectionY(hdname.c_str(),fCurChan_beginningOfSignal+1,fCurChan_beginningOfSignal+1,"")->DrawCopy();
    }
}

void PlaneRPD_SADC::CalculateAdc()
{
  for(unsigned short int iChannel=0; iChannel<N_CHANNEL_SADC; ++iChannel)
    {
      if(fDataSent[iChannel])
//eb	for(unsigned short int iSample=1; iSample<N_SADC_SAMPLE; ++iSample)
// discard first 2 channels because they do not make sense
	for(unsigned short int iSample=2; iSample<N_SADC_SAMPLE; ++iSample)
	  {
	    fAdc[iChannel]+=fAmplitude[iChannel][iSample]
// 	      -fPedestal[iChannel]*fCalculatedPedestal[iChannel]
	      ;
	  }
    }
}

void PlaneRPD_SADC::Reset( void ) {
  Plane::Reset();
  
  for(unsigned short int iChannel=0; iChannel<N_CHANNEL_SADC; ++iChannel)
    {
      fDataSent[iChannel]=false;
      fCalculatedPedestal[iChannel]=false;
      fPedestal[iChannel]=0;
      fPedestalRelativeError[iChannel]=0;
      fAmplitudeMaximum[iChannel]=0;
      fAdc[iChannel]=0;
      fNumberOfSampleAboveThreshold[iChannel]=0;
      for(unsigned short int iSample=0; iSample<N_SADC_SAMPLE; ++iSample)
	{
	  fAmplitude[iChannel][iSample]=0;
	}
    }
  return;
}

void PlaneRPD_SADC::CalculateNumberOfSampleAboveThreshold()
{
  for(unsigned short int iChannel=0; iChannel<N_CHANNEL_SADC; ++iChannel)
    {
      if(fDataSent[iChannel])
	{
	  for(unsigned short int iSample=2; iSample<N_SADC_SAMPLE; ++iSample)
	    {
	      if(fAmplitude[iChannel][iSample]>3)++fNumberOfSampleAboveThreshold[iChannel];
	    }
// 	  if(numberOfSampleAboveThreshold<4)fDataSent[iChannel]=false;
	}
    }
}

int Modulo(int value, int modulo)
{
  while(value < 0)
    value+=modulo;
  while(value > modulo)
    value-=modulo;
  return value;
}
