#ifndef __PlaneRPD_F1__
#define __PlaneRPD_F1__

#include "Plane1V.h"
#include "PlanePanel.h"

#define N_CHANNEL_PLANE_RPD_F1_A 12
#define N_CHANNEL_PLANE_RPD_F1_B 64
#define N_CHANNEL_PLANE_RPD_F1 64//(N_CHANNEL_PLANE_RPD_F1_A+N_CHANNEL_PLANE_RPD_F1_B)
#define F1_MULTIPLICITY 5

#include <sstream>
#include <string>

#include <TTree.h>
#include <TH1F.h>
#include <TH2F.h>

class PlaneRPD_F1_Panel;

/*! \brief Describe a plane for F1 timming cards on the RPD
  \author Guillaume Jegou, Etienne Burtin
*/

class PlaneRPD_F1 : public Plane1V
{
 public:
  PlaneRPD_F1(const char *detname,int nchan, int center, int width);
  ~PlaneRPD_F1();

  //! For each channel, and for the five first hits of its, fDataSent[iChannel][iHit] is true if a time have been measured
  bool fDataSent[N_CHANNEL_PLANE_RPD_F1][F1_MULTIPLICITY];
  //! For each channel, and for the five first hits of its, fTdcRpd[iChannel][iHit] contains the measured time. Make sence only if fDataSent[iChannel][iHit] is true
  double fTdcRpd[N_CHANNEL_PLANE_RPD_F1][F1_MULTIPLICITY];
  //! Difference between fTdcRpd[iChannel][0] - fTdcRpd[iChannel+N_CHANNEL_PLANE_RPD_F1/2][0]; not calculated for all mutiplicities... to do
  double fTimeRpd[N_CHANNEL_PLANE_RPD_F1][F1_MULTIPLICITY];
  //! the run number
  int fEventRunNumber;
  //! the event number in the run fEventRunNumber
  int fEventNumberInRun;
  //! the time of the event (sec)
  double fEventTimeSecond;
  //! the time of the event (Musec)
  double fEventTimeMicroSecond;
  
  //! Initialisation of the Plane. Mainly make the instantiation of histograms
  void Init(TTree* tree = 0);
  
  //! Create a new PlaneRPD_F1_Panel if it do not exist.
  void ControlPanel(const TGWindow* p, const TGWindow* main) ;
  //! returns the number of channels in this plane
  int GetNchannels() const {return N_CHANNEL_PLANE_RPD_F1;}

  //! Use to make projection in the controlPanel (contains the current position of the slider)
  int fCurChan_tdc;
/*   int fCurChan_multiplicity; */
  //! Use to make projection in the controlPanel (set the current position of the slider)
  void SetChannel(int channel) {fCurChan_tdc=channel;}
/*   void SetChannel_2(int channel) {fCurChan_multiplicity=channel;} */
  //! Draw the projection of h2TdcVsChannel for fCurChan_tdc bin
  void TdcSpectrum();
/*   void MultiplicitySpectrum(); */

#ifndef __CINT__
  //! This function read the data and store it in the lDigits member
  void StoreDigit(CS::Chip::Digit* digit);
  //! called at the end of each event, fill members from lDigits
  void EndEvent(const CS::DaqEvent &event);
#endif

  //! Called between each event to re-initialise members
  void Reset();

  //! Spectrum : numbers of channels which sent data
  TH1F *h1NumberOfSentData;
  //! Tdc spectrum for each channel
  TH2F *h2TdcVsChannel;
  //! Multiplicity spectrum for each channel
  TH2F *h2MultiplicityVsChannel;
/*   TH1F *h1Channel; */

  ClassDef(PlaneRPD_F1,2)
};

#endif





