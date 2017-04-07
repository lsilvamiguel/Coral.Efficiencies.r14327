#ifndef __PlaneDriftChamberPanel__
#define __PlaneDriftChamberPanel__


#include "Plane1VPanel.h"
#include "PlaneDriftChamber.h"
#include "PlaneDriftChamber2V.h"

double TimeFitFunction(double *x,double *par);

// I wanted the following to be  a private class of  PlaneDriftChamberPanel, but 
// did not manage to do so without root error in Dict.cc and Dict.h
// HUGO
  
class FitData: public TObject {
  public:
  FitData( ):TObject() {};
  FitData( int iMin, int iMax ): TObject(),iMin(iMin), iMax(iMax) {};
  int iMin;
  int iMax;
  std::vector<double> pars;
  double chi2;   
  unsigned int ndf;
  unsigned int entries;   
  ClassDef(FitData,1)
};

class PlaneDriftChamberPanel : public Plane1VPanel {
  public:
  PlaneDriftChamberPanel(const TGWindow *p,const TGWindow *main, UInt_t w, UInt_t h, PlaneDriftChamber *plane);
  PlaneDriftChamberPanel(const TGWindow *p,const TGWindow *main, UInt_t w, UInt_t h, PlaneDriftChamber2V *plane);
  virtual ~PlaneDriftChamberPanel();
  
  /// Creates custom tab for DC drift Time fit
  void CreateControlTab( void );
  void Fit(         int id=-1 );
  void WriteDBFile( int id=-1 );

  private:
  
  static const char* DBFileTypes[];
  TGTextButton *fFitBt, *fWriteBt;
  TRootEmbeddedCanvas *fTimeDpy;
  TGTextEntry *fChStepTE, *fTMinTE, *fTMaxTE;
  TH2F *fHtvsch;
  
  void DeleteFitDatas( void );
  std::list< FitData* > fFitDatas;
  
  ClassDef(PlaneDriftChamberPanel,1)

};

#endif
