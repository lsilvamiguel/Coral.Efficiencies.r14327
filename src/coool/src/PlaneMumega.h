#ifndef __PlaneMumega__
#define __PlaneMumega__

#ifndef __CINT__
#include "config.h"
#endif

#include <string>
#include <vector>
#include <map>
#include <set>

#include "Plane.h"
#include "PlanePanel.h"

#include "GeomPlane.h"

#include "TThread.h"

static const int fMaxTriggerNumber=12;

///class for Micromegas - (ch,data1,data2) digit type
class PlaneMumega : public Plane {

 public:

  /// temporary
  std::set<int> noisy;

  /// class needed for calibration purposes
  class ChannelCalib {
  public:
    int ch;
    float t0;
    int on;
    ChannelCalib() : ch(0),t0(0),on(0) {}
    ChannelCalib(const char *s) {
      int nbvar = sscanf(s,"%d%f%d",&ch,&t0,&on);
      if ((nbvar != 2) || (nbvar != 3))  {
	throw CS::Exception("PlaneMumega::ChannelCalib : bad line \"%s\"",s);
	std::cerr<<"bad line, exception not caught !"<<std::endl;
      }
      if (nbvar == 2) on = 0;
    }
  };

  friend istream& operator>>(istream& in, PlaneMumega::ChannelCalib &c) {
    in>>c.ch;
    in>>c.t0;
    in>>c.on;
    return in;
  }

  /*! \brief constructor
    \param detname detector name
    \param nchan number of channels
    \param center center of the time band
    \param width width of the time band
  */
  PlaneMumega(const char *detname,int nchan, int center, int width, int parity=0);
  virtual ~PlaneMumega();

  /// Passes a digit to the plane
  void StoreDigit(int channel, int data);

#ifndef __CINT__
  void StoreDigit(CS::Chip::Digit* digit);

  /*! Applies the cuts, fills histograms and tree.
    This function must be called at the end of the event
  */
  void EndEvent(const CS::DaqEvent &event);
#endif

  /// do the clustering
  void Clusterize();

  /// book histograms and branch the tree
  void Init(TTree* tree =0);

  /// \return number of channels
  int GetNchannels() {return fNchan;}

  /// \return clusters
  //const std::set<CCluster1>& GetClusters() const {return fClusters;}

  /// Resets all histograms (and associated counters)
  void ResetHistograms();

#if USE_DATABASE == 1
  /// method to read calibration data
  void ReadCalib(const tm& t);
#endif

  void ControlPanel(const TGWindow *p, const TGWindow *main);

  /// text output
  void TextOutput(ostream& out);

 private:

  /// maximum multiplicity per tdc channel
  static const int fMAX_MULT;

  /// F1 TDC tick in ms
  static const float fF1_TICK;

  /// F1 time window in F1 clock unit
  static const int fF1_GATE;

  /// Period in event numbers of rate histo update
  static const int fRATE_UPDATE;

  /// Weighting Leading and Trailing
  static const float fL_WGHT, fT_WGHT;

  /// Hit time window
  static const float fHIT_TMIN, fHIT_TMAX;

  /// Coinc time
  static const float fMT_T0;

  /// Cluster time window
  static const float fCL_TMIN, fCL_TMAX;

  /// number of channels
  int    fNchan;

  /// map of digits. indexed by channel number (ie ordered !). Replace with a set ??
  std::map<int,CDigit1> fDigits;

  /// Noisy channels
  //std::set<int> fNoisyCh;

  /// number of clusters passing the cuts.
  int fNCl;

  /// parity of the last bit of tdc data for leading edge.
  int fPARITY_LEAD;

  /// 1D histograms
  TH1F *fHlhits, *fHlch, *fHlt,*fHtt, *fHhhits, *fHhch, *fHht, *fHhtot, *fHhrates, *fHhtcalib, *fHchits, *fHcch, *fHct, *fHctot, *fHcs, *fHcamp, *fHcampsig, *fHcres;
  TH1F *fHhtrigvsrates[fMaxTriggerNumber];

  /// number of hch hits for each trigger type, used by fHhtrigvsrates histos
  int fHhNtrig[fMaxTriggerNumber];

  /// number of each trigger type, used by fHhtrigvsrates histos
  int fNtrig[fMaxTriggerNumber];

  /// 2D histograms
  TH2F *fHhtvsch,*fHhtcalibvsch, *fHctvsch, *fHctotvst, *fHhtrigvsch;

  /// variables
  Variable *fVch, *fVt, *fVhch, *fVht, *fVhlt, *fVhtt, *fVhtcalib, *fVhtot, *fVcch, *fVct, *fVctot, *fVcs, *fVcamp, *fVcres;

  /// calibration for every channel
  std::vector<ChannelCalib> calib_data;

  /// missing channels
  std::set<int> fMissingChans;

  /// noisy channels
  std::set<int> fNoisyChans;

  /// \return true is this data word is leading time
  bool IsLeading(int time) {return ( (time&0x1) == fPARITY_LEAD);}

  /// missing and noisy channels
  bool ChannelProblems();

  /// fit fHcampsig to estimate gain
  float Gain();

  ClassDef(PlaneMumega,0)
};

#endif







