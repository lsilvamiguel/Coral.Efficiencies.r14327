#ifndef __PlaneAPV__
#define __PlaneAPV__

#include "Plane.h"
#include "TThread.h"

/// Plane for detectors with a readout based on an APV chip.

class PlaneAPV : public  Plane {

  //private:


 protected:

  // run maps for all apv detectors
#ifndef __CINT__
  static const CS::Chip::Maps *fRunMaps;
#endif

  std::map<int,CDigit1> fDigits;

  /// maximum multiplicity per channel
  static const int fMAX_MULT;

  int fNchan;

  bool fPixel;

  /// variables
  Variable *fVch, *fVa0, *fVa1, *fVa2, *fVt;

  /// histograms
  TH1F *fHhits, *fHch, *fHa0, *fHa1, *fHa2, *fHt;
  TH2F *fHhitsp;

 public:

  PlaneAPV(const char *detname,int nchan, int center, int width, bool pixel=false);
  virtual ~PlaneAPV();

  virtual void StoreDigit(int channel, int amp1, int amp2, int amp3, int cmc0, int cmc1, int cmc2, int wirepos, int chip);

  virtual void StoreDigit(int x, int y, int channel, int amp1, int amp2, int amp3, int cmc0, int cmc1, int cmc2, int chip);

#ifndef __CINT__
  virtual void StoreDigit(CS::Chip::Digit* digit);

  virtual void EndEvent(const CS::DaqEvent &event);
#endif

  /// Resets all histograms (and associated counters)
  void ResetHistograms();

  /// Sets decoding run maps for calibration reading
#ifndef __CINT__
  static void SetRunMaps(const CS::Chip::Maps* maps) {
    if(!fRunMaps) fRunMaps = maps;
  }
#endif

  virtual void Init(TTree* tree =0);

  virtual void ControlPanel(const TGWindow* p, const TGWindow* main);

  /// \return number of channels
  int GetNchannels() const {return fNchan;}

  ClassDef(PlaneAPV,0)
};

#endif











