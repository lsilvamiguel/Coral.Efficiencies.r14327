#ifndef __PlaneScifiG__
#define __PlaneScifiG__

#include "Plane.h"
#include "PlanePanel.h"

#include "TThread.h"


/// class for German Scintillating Fibers -  (ch,data1,data2) digit type
class PlaneScifiG : public Plane {

 private:

#ifndef __CINT__
  std::map<int,CDigit1> fDigits;
#endif  // __CINT__

  /// maximum multiplicity per tdc channel
  static const int fMAX_MULT;

  /// number of channels
  int    fNchan;

  /// channel for each digit
//   int  *fChannel;

  /// data for each digit
//   int  *fData;

  /// pos for each digit
//   int  *fPos;

  /// channel variable. 
  Variable *fVch;

  /// time low threshold variable. 
  Variable *fVtl;

  /// time high threshold variable. 
  Variable *fVth;

  /// status varibale
  Variable *fVstat;

  /// time diff variable
  Variable *fVtimediff;
 
 /// number of hits histogram
  TH1F *fHhits;

  /// hit map high threshold  histogram
  TH1F *fHchh;

  /// hit map low threshold  histogram
  TH1F *fHchl;

  /// time spectrum low threshold  histogram
  TH1F *fHtl;

  /// time spectrum high threshold  histogram
  TH1F *fHth;

  /// th-tl vs channel  histogram
  TH2F *fHdtvsch;

  /// tl vs channel  histogram
  TH2F *fHtlvsch;

  /// th vs channel  histogram
  TH2F *fHthvsch;

  /// time diff histogram
  TH1F *fHtimediff;

  /// timediff vs channel histogram
  TH2F *fHtdvsch;

  /// status histogram
  TH1F *fHstat;

  //corrected time
  Variable *fVtc;
  TH1F *fHtc;

  /// hit map high threshold  histogram
  TH2F *fHchhVsTis;

  /// hit map low threshold  histogram
  TH2F *fHchlVsTis;

  /// check if this digit corresponds to the low threshold
  bool IsLowThreshold(int pos) {
    if(pos == -1) return true;
    else return false;
  }

  /// true if this uses high resolution time
//   bool fIsHR;

  /// cluster histograms
  TH1F *fHchits, *fHcch, *fHcs;
  /// cluster variables
  Variable *fVcch, *fVcs;
  int fNclustKept;

#ifndef __CINT__
  // SciFiG calibration reading: 1 calibration variable per channel
  class ChannelCalib {
  public:
    int ch;
    float t0;
    int flag;
    ChannelCalib() : ch(0),t0(0),flag(0) {}
    ChannelCalib(const char *s) {
      if(3 != sscanf(s,"%d%f%d",&ch,&t0,&flag)) {
	throw CS::Exception("PlaneSciFiG::ChannelCalib : bad line \"%s\"",s);
	std::cerr<<"bad line, exception not caught !"<<std::endl;
      }
    }
  };
  
  /// calibration for every channel
  std::vector<ChannelCalib> calib_data;

  friend istream& operator>>(istream& in, PlaneScifiG::ChannelCalib &c) {
    in>>c.ch;
    in>>c.t0;
    in>>c.flag;
    return in;
  }
#endif  // __CINT__



 public:
  /*! \brief constructor
    \param detname detector name
    \param nchan number of channels
    \param center center of the time band
    \param width width of the time band
  */
  PlaneScifiG(const char *detname,int nchan, int center, int width);
  virtual ~PlaneScifiG();

  /// Passes a digit to the plane
  void StoreDigit(int channel, int data, int pos);

#ifndef __CINT__
  void StoreDigit(CS::Chip::Digit* digit);

  /*! Applies the cuts, fills histograms and tree.
    This function must be called at the end of the event
  */
  void EndEvent(const CS::DaqEvent &event);
#endif

  /// Set/Unset high resolution
//   void SetHRTime(bool ishr) {fIsHR=ishr;}

  /// book histograms and branch the tree
  void Init(TTree* tree =0);
  
  void Clusterize();

  /// \return number of channels
  int GetNchannels() {return fNchan;}

  void ControlPanel(const TGWindow *p, const TGWindow *main);

#if USE_DATABASE == 1
  /// method to read calibration data
  void ReadCalib(const tm& t);
#endif


  ClassDef(PlaneScifiG,1)
};

#endif




