#ifndef __GroupRPD__
#define __GroupRPD__

#include "Group.h"
#include "PlaneRPD_SADC.h"
#include "PlaneRPD_F1.h"

namespace CS { class DaqEvent; }

class TH1F;
class TH2F;

class PlaneRPD_F1;
class PlaneRPD_SADC;

class GroupRPD : public Group
{

 protected:

  const PlaneRPD_F1 *pf1_A;
  const PlaneRPD_F1 *pf1_Bl;
  const PlaneRPD_F1 *pf1_Bh;
  const PlaneRPD_SADC *psadcl;
  const PlaneRPD_SADC *psadch;
  bool fOkPlanes;


 public:
  GroupRPD(const char* name);
  ~GroupRPD(){};
  void Init(void);

#ifndef __CINT__
  //! Fills the histograms
  void EndEvent(const CS::DaqEvent &event);
#endif

#ifndef __CINT__
  //! Histogram : beta on the X axis and deposited energy on the Y axis
  TH2F **h2EnergyLostInBVsBeta;
  //! Histogram : on the X axis and on the Y axis
  TH2F **h2EnergyLostInAVsBeta;
  //! Histogram : tof spectrum
  TH1F **h1Tof;
  //! Histogram : position for A
  TH1F **h1Za;
  //! Histogram : position for B
  TH1F **h1Zb;
  //! Histogram : time difference in A
  TH1F **h1TupMinusTdownInA;
  //! Histogram : time differnce in B
  TH1F **h1TupMinusTdownInB;
   //! Histogram : cf histo name
   TH2F **h2AdcUpVsAdcDown;
  //! Histogram : extrapolated z vertex
  TH1F *h1ZVertex;
  //! Histogram : cf name
  TH2F **h2AdcVsPosition;   
  //! Histogram : A multiplicity
  TH1I *h1MultiplicityOfTriggerAB;
 //! Histogram : List of cuts
  TH1I *h1EfficiencyAfterCuts;   
  //! Histogram : which AB pair was analysable
  TH1I *h1IdentifiedSectorForTriggerABWithCuts;    
  //! The histo name should tell you what it contains
  TH2F **h2EnergyAVsEnergyB;
  //! Later
  TH1F *h1tdcA0Down;
#endif
    
  //! Returns the software channel corresponding the argument (the hardware channel)
  int Pm2Channel(char scint, unsigned int n, char side, TString type);
    
  //! Returns the offset of time difference between up and down Pm for the scintillator in argument
  float GetTimeDifferenceOffset(char scintillator, unsigned int n);
    
  //! Returns the light speed in the scintillator in argument
  float GetLightSpeed(char scintillator, unsigned int n);

  //! Returns time difference between up and down Pm for the scintillator in argument (for the first hit in Up and Down F1 card)
  float GetTimeDifference(char scintillator, int n, const PlaneRPD_F1 *pf1, const PlaneRPD_SADC *psadc);

  //! Returns position for the scintillator in argument  (for the first hit in Up and Down F1 card)
  float GetPosition(char scintillator, int n, const PlaneRPD_F1 *pf1, const PlaneRPD_SADC *psadc);

  //! Returns time sum of up and down Pm for the scintillator in argument (for the first hit in Up and Down F1 card)
  float GetTimeSum(char scintillator, int iScintillator, const PlaneRPD_F1 *pf1, const PlaneRPD_SADC *psadc);    

  //! Return the time of flight between the two scintillators in arguments (for the first hit in Up and Down F1 card of A and B scintillator)
  float GetTof(int iScintillatorA,int iScintillatorB, const PlaneRPD_F1 *pf1_A, const PlaneRPD_F1 *pf1_B, const PlaneRPD_SADC *psadc);

  //! Returns the offset of the time of flight between the two scintillators in arguments (for the first hit in Up and Down F1 card of A and B scintillator)
  float GetTofOffset(int iScintillatorA, int iScintillatorB);

  //! Returns the calculated beta for the two scintillators in arguments (for the first hit in Up and Down F1 card of A and B scintillator)
  float GetBeta(int iScintillatorA,int iScintillatorB, const PlaneRPD_F1 *pf1_A, const PlaneRPD_F1 *pf1_B, const PlaneRPD_SADC *psadc);

  //! Calculate and return the extrapolated Z position in the target by assuming R=0cm for the two scintillators in arguments (for the first hit in Up and Down F1 card of A and B scintillator)
  float GetZPositionInTarget(int iScintillatorA,int iScintillatorB, const PlaneRPD_F1 *pf1_A, const PlaneRPD_F1 *pf1_B, const PlaneRPD_SADC *psadc);

  //! Returns calculated walk correction for the scintillator in argument  (for the first hit in Up and Down F1 card)
  float GetWalkCorrection(char scintillator, int n, char side , const PlaneRPD_SADC *psadc);
    
  //! Return true if the scintillator Ai with one of the facing B is is able to trig
  bool CanAiTrig(unsigned int iA, const PlaneRPD_F1 *pf1_A, const PlaneRPD_F1 *pf1_B, const PlaneRPD_SADC *psadc, const unsigned int iB=100);


  ClassDef (GroupRPD,1)
};

#endif
