#ifndef __PlaneRPD_SADC__
#define __PlaneRPD_SADC__

#include "Plane.h"
#include "PlanePanel.h"


#define N_CHANNEL_SADC_A 24
#define N_CHANNEL_SADC_B 48
#define N_CHANNEL_SADC 72
#if (N_CHANNEL_SADC_A+N_CHANNEL_SADC_B) != N_CHANNEL_SADC
#error "In PlaneRPD_SADC :: N_CHANNEL_SADC_A+N_CHANNEL_SADC_B != N_CHANNEL_SADC"
#endif
#define N_SADC_SAMPLE 32

using namespace std;

class PlaneRPD_SADC_Panel;

/*! \brief Describe a plane for SADC cards on the RPD
  \author Guillaume Jegou, Etienne Burtin
*/
class PlaneRPD_SADC : public Plane
{
 public:  

  PlaneRPD_SADC(const char *detname);
  ~PlaneRPD_SADC();

  //! For each channel, contains the N_SADC_SAMPLE amplitudes of the signal. (make sence if fDataSent[iChannel]==true)
  double fAmplitude[N_CHANNEL_SADC][N_SADC_SAMPLE];
  //! For each channel, contains the amplitude maximum. (make sence if fDataSent[iChannel]==true)
  double fAmplitudeMaximum[N_CHANNEL_SADC];
  //! For each channel, contains the adc (ie : the integrated amplitude). (make sence if fDataSent[iChannel]==true)
  double fAdc[N_CHANNEL_SADC];
  //! For each channel, is true if the sadc card sent some data
  bool fDataSent[N_CHANNEL_SADC];
  //! For each channel, is true the a pedestal have been calculated
  bool fCalculatedPedestal[N_CHANNEL_SADC];
  //! For each channel, contains the calculated pedestal. (make sence if fCalculatedPedestal[iChannel]==true)
  double fPedestal[N_CHANNEL_SADC];
  //! For each channel, contains the error on the pedestal which comes from the sample which was used to calculate the pedestal. (make sence if fCalculatedPedestal[iChannel]==true)
  double fPedestalRelativeError[N_CHANNEL_SADC];
  //! Contains the bin corresponding to the beginning of the signal . (make sence if fDataSent[iChannel]==true)
  int fBeginningOfSignal[N_CHANNEL_SADC];
  //! Contains for each channel, the number of sample higher the a given value, after pedestal substraction. (make sence etc etc....)
  int fNumberOfSampleAboveThreshold[N_CHANNEL_SADC];
  //! //! the time of the event (sec)
  double fEventTimeSecond;
  //! the time of the event (Musec)
  double fEventTimeMicroSecond;
  //! the run number
  int fEventRunNumber;
  //! the event number in the run fEventRunNumber
  int fEventNumberInRun;
    
  //! Initialisation of the Plane. Mainly make the instantiation of histograms
  void Init(TTree* tree = 0);
  //! Create a new PlaneRPD_SADC_Panel if it do not exist.
  void ControlPanel(const TGWindow* p, const TGWindow* main) ;
  int GetNchannels() const {return N_CHANNEL_SADC;}

  //! Use to make projection in the controlPanel (contains the current position of the maxAmplitude slider)
  int fCurChan_maxAmplitude;
  //! Use to make projection in the controlPanel (contains the current position of the adc slider)
  int fCurChan_adc;
  //! Use to make projection in the controlPanel (contains the current position of the pedestal slider)
  int fCurChan_pedestals;
  //! Use to make projection in the controlPanel (contains the current position of the beginningOfSignal slider)
  int fCurChan_beginningOfSignal;
  //! Use to make projection in the controlPanel (set the current position of the maxAmplitude slider)
  void SetChannel(int channel) {fCurChan_maxAmplitude=channel;}
  //! Use to make projection in the controlPanel (set the current position of the adc slider)
  void SetChannel_2(int channel) {fCurChan_adc=channel;}
  //! Use to make projection in the controlPanel (set the current position of the pedestal slider)
  void SetChannel_3(int channel) {fCurChan_pedestals=channel;}
  //! Use to make projection in the controlPanel (set the current position of the beginningOfSignal slider)
  void SetChannel_4(int channel) {fCurChan_beginningOfSignal=channel;}
  //! Draw the projection of h2AmplitudeMaximumVsChannel for fCurChan_maxAmplitude bin
  void MaxAmplitudeSpectrum();
  //! Draw the projection of h2AdcVsChannel for fCurChan_adc bin
  void AdcSpectrum();
  //! Draw the projection of h2PedestalsVsChannel for fCurChan_pedestals bin
  void PedestalsSpectrum();
  //! Draw the projection of h2BeginningOfSignalVsChannel for fCurChan_beginningOfSignal bin
  void BeginningOfSignalSpectrum();

#ifndef __CINT__
   //! This function read the data and store it in the lDigits member
  void StoreDigit(CS::Chip::Digit* digit);
  //! called at the end of each event, fill members from lDigits
  void EndEvent(const CS::DaqEvent &event);
#endif
  void Reset();

  //! The histo name should tell you what it contains
  TH1F *h1NumberOfSentData;
  //! The histo name should tell you what it contains
  TH2F *h2PedestalsVsChannel;
  //! The histo name should tell you what it contains
  TH2F *h2AmplitudeMaximumVsChannel;
  //! The histo name should tell you what it contains
  TH2F *h2AdcVsChannel;
  //! The histo name should tell you what it contains
  TH2F *h2BeginningOfSignalVsChannel;
  //! The histo name should tell you what it contains
  TH1F *h1Channel;
  //! The histo name should tell you what it contains
  TH2F *h2EnergyVsChannel;
  //! The histo name should tell you what it contains
/*   TH2F **h2AdcUpVsAdcDown; */
/*   //! The histo name should tell you what it contains */
  TH1F **h1SignalShape;
/*   TH2F **h2SignalShape; */
  //! The histo name should tell you what it contains
  TH2F *h2NumberOfSampleAboveThreshold;
  
  
 private :
  //! calculate the pedestal
  void LookForPedestals();
  //! Apply calculated pedestal
  void CorrectPedestal();
  //! calculate the beginning of the signal
  void LookForBeginningOfSignal(unsigned short int channel);
  //! calculate sample corresponding to the maximum amplitude in the signal
  unsigned short int LookForSampleOfMaximumAmplitude(unsigned short int channel);
  //! calculate the maximum amplitude in the signal
  void LookForMaximumAmplitude();
  //! calculate the adc (integrate the signal)
  void CalculateAdc();
  //! calculated the number of bins above threshold
  void CalculateNumberOfSampleAboveThreshold();

  ClassDef(PlaneRPD_SADC,1)
};

#endif
