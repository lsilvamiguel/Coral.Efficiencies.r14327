#ifndef __Plane2V__
#define __Plane2V__

#include "Plane.h"
#include "PlanePanel.h"

#include "TThread.h"


/// Base class for detectors having a (col,row,data) digit type
class Plane2V : public Plane {
  
 protected:

  /// maximum multiplicity per adc channel
  static const int fMAX_MULT;   

  /// number of rows
  int    fNrows;

  /// number of columns
  int    fNcols;

  /// number of channels
  int fNchan;

  /// current channel looked at in Plane2VPanel
  int fCurChan;

  /// row for each digit          
  int  *fRow; //[fNchan*fMAX_MULT]    comment needed by rootcint -don't delete
  /// column for each digit          
  int  *fCol; //[fNchan*fMAX_MULT]

  /// amplitude for each digit        
  int  *fAmp; //[fNchan*fMAX_MULT]

  /// variables
  Variable *fVrow, *fVcol, *fVamp;

  // markers for histograms
  TH1F *fHhit;
  TH2F *fHrca, *fHrc1;
  TH2F *fHavsadr;
  TH1F *fHa;

 public:
  /*! \brief constructor
    \param detname detector name
    \param nchan number of channels
    \param center center of the time band
    \param width width of the time band
  */
  Plane2V(const char *detname,int ncols, int nrows, int center, int width);
  virtual ~Plane2V();
  
  /// Passes a digit to the plane
  void StoreDigit(int col, int row, int amp);

#ifndef __CINT__
  virtual void StoreDigit(CS::Chip::Digit* digit);

  virtual void EndEvent(const CS::DaqEvent &event);  
#endif

  /// book histograms and branch the tree
  virtual void Init(TTree* tree =0);

  /// sets the channel number for the projection of the amp vs chan histogram
  void SetChannel(int channel) {fCurChan=channel;}  

  void AmpSpectrum();

  /// \return number of channels
  int GetNchannels() {return fNrows*fNcols;}
  int GetNrows() {return fNrows;}
  int GetNcols() {return fNcols;}

  /// get hit histograms
  TH1F* GetHhit() {return fHhit;}



  /// get amplitude
  TH1F* GetAmp() {return fHa;}

  
  
  
  virtual void ControlPanel(const TGWindow *p, const TGWindow *main);

  ClassDef(Plane2V,1)
};

#endif



