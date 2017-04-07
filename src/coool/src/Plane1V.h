#ifndef __Plane1V__
#define __Plane1V__

#include "Plane.h"

#include <sstream>
#include <string>

#include <TTree.h>


class Plane1VPanel;

/*! \brief Base class for detectors having a (ch,data) digit type
  \author Colin Bernet
*/

class Plane1V : public Plane {

 protected:

#ifndef __CINT__
  std::map<int,CDigit1> fDigits;
  // Default calibration reading: 1 calibration variable per channel
  class ChannelCalib {
  public:
    int ch;
    float t0;
    ChannelCalib() : ch(0),t0(0) {}
    ChannelCalib(const char *s) {
      if(2 != sscanf(s,"%d%f",&ch,&t0)) {
	throw CS::Exception("Plane1V::ChannelCalib : bad line \"%s\"",s);
	std::cerr<<"bad line, exception not caught !"<<std::endl;
      }
    }
  };

  /// fill fDigits structure in StoreDigit method
  bool fFillfDigits;

  /// calibration for every channel
  std::vector<ChannelCalib> calib_data;

  friend istream& operator>>(istream& in, Plane1V::ChannelCalib &c) {
    in>>c.ch;
    in>>c.t0;
    return in;
  }
#endif  // __CINT__

  /// maximum multiplicity per tdc channel
  static const int fMAX_MULT;

  public:

  /// F1 TDC tick in ms
  static const float fF1_TICK;
  
  protected:

  // Period in event numbers of rate histo update
  static const int fRATE_UPDATE;

  /// number of channels
  int    fNchan;

  /// channel for each digit
//   int  *fChannel;

  /// time for each digit
//   int  *fTime;

  /// current channel looked at (see Plane1VPanel)
  int fCurChan;

  /// Variables
  Variable *fVch, *fVt;

  /// cut on data;
  float fCenter, fWidth;

  /// number of hits per event
  TH1F *fHhits;

  /// hit profile
  TH1F *fHch;

  /// time distribution
  TH1F *fHt;

  /// time vs channel
  TH2F *fHtvsch;

  /// On trigger time
  Variable *fOnTrigTVarID;

  /// Off trigger time
  Variable *fOffTrigTVarID;

  /// on-Trigger time distribution
  TH1F *fHonTT;

  /// on-Trigger profile
  TH1F *fHonTP;

  ///off-Trigger time distribution
  TH1F *fHoffTT;

  ///off-Trigger profile
  TH1F *fHoffTP;

  ///corrected on-trigger profile
  TH1F *fHoncorrTP;

  ///rates
  TH1F *fHrates;

  /// number of hits kept in the on-trigger time window
  int fNhitsKeptOnTrig;

  /// number of hits kept in the off-trigger time window
  int fNhitsKeptOffTrig;

  /// missing channels
  std::set<int> fMissingChans;

  /// noisy channels
  std::set<int> fNoisyChans;

  /// missing and noisy channels
  virtual bool ChannelProblems();

 public:
  /*! \brief constructor
    \param detname detector name
    \param nchan number of channels
    \param center center of the time band
    \param width width of the time band
  */
  Plane1V(const char *detname,int nchan, int center, int width);
  virtual ~Plane1V();

  /// Passes a digit to the plane
  virtual void StoreDigit(int channel, int time);

#ifndef __CINT__
  virtual void StoreDigit(CS::Chip::Digit* digit);

  /// Applies the cuts, fills histograms and tree. This function must be called at the end of the event
  virtual void EndEvent(const CS::DaqEvent &event);
#endif
  /// book histograms and branch the tree
  virtual void Init(TTree* tree =0);

  /// Resets plane (called at each event)
  virtual void Reset();

  /// Resets all histograms (and associated counters)
  void ResetHistograms();

  /// create the control panel
  virtual void ControlPanel(const TGWindow* p, const TGWindow* main);

    /// text output
  void TextOutput(ostream& out);

  /// \return number of channels
  int GetNchannels() const {return fNchan;}

  /// \return variables
  float* GetTimeValues() const {return fVt->GetValues();}
  float* GetChanValues() {return fVch->GetValues();}

  /// \return variables
  const Variable& GetTimeVariable() const {return *fVt;}
  const Variable& GetChanVariable() const {return *fVch;}

  /// sets the channel number for the projection of the time vs chan histogram
  void SetChannel(int channel) {fCurChan=channel;}

  /// Projects the time vs chan histogram to the t axis in a predefined channel
  void TimeSpectrum();

  const Variable & GetOffTrigTVarID(void) const {return *fOffTrigTVarID;}


  ClassDef(Plane1V,1)
};



#endif





