#ifndef __PlaneScifiJ__
#define __PlaneScifiJ__

#include "Plane1V.h"
#include "TThread.h"

#include "TTree.h"

/// Plane for japanese fiber hodoscopes

class PlaneScifiJ : public  Plane1V {

private:

#ifndef __CINT__
  // SciFiJ calibration reading: 1 calibration variable per channel
  class ChannelCalib {
  public:
    int ch;
    float t0;
    int flag;
    ChannelCalib() : ch(0),t0(0),flag(0) {}
    ChannelCalib(const char *s) {
      if(3 != sscanf(s,"%d%f%d",&ch,&t0,&flag)) {
	throw CS::Exception("PlaneSciFiJ::ChannelCalib : bad line \"%s\"",s);
	std::cerr<<"bad line, exception not caught !"<<std::endl;
      }
    }
  };
  
  /// calibration for every channel
  std::vector<ChannelCalib> calib_data;

  friend istream& operator>>(istream& in, PlaneScifiJ::ChannelCalib &c) {
    in>>c.ch;
    in>>c.t0;
    in>>c.flag;
    return in;
  }
#endif  // __CINT__
  
  // private:
  /// variables
  Variable  *fVtc, *fVhch, *fVht, *fVhtcalib;
  TH1F *fHtc;
  TH2F *fHtcVSch;
  /// cluster histograms
  TH1F *fHchits, *fHcch, *fHcs;
  /// cluster variables
  Variable *fVcch, *fVcs;
  int fNclustKept;
  // channel versus time in spill
  TH2F *fHchVsTis;

 public:
  
  PlaneScifiJ( const char *detname,int nchan, int center, int width):
    Plane1V(detname,nchan,center,width) {INeedGeom();}
  
  ~PlaneScifiJ() {}
  
  void Init(TTree* tree =0);
  
  void Clusterize();

#ifndef  __CINT__
  void EndEvent(const CS::DaqEvent &event);
#endif

#if USE_DATABASE == 1
  /// method to read calibration data
  void ReadCalib(const tm& t);
#endif

  ClassDef(PlaneScifiJ,1)
};

#endif




