#ifndef __PlaneRICH_MAPMT__
#define __PlaneRICH_MAPMT__

#include "Plane2V.h"
#include "PlanePanel.h"
#include "TStyle.h"

/* ---------------------------------------------------------------------
   Maintainer & author: Frank Nerling 
   New Class for monitoring RICH_MAPMT: (x, y, time), created April 2006
   Updated Version: 30.08.2006 
   --------------------------------------------------------------------- */

class PlaneRICH_MAPMT : public Plane2V {
  
 private:

  TH2F *fHtime_vs_chan, *fHrc1b; 
//  TH2F *fHrca_hits;
  TH1F *fHtime_vs_chan_1dim;
  TH1F *fHtc;

 public:
  /*! \brief constructor
    \param detname quadrant name 
    \param ncols number of columns per MAPMT-quadrant
    \param nrows number of rows per MAPMT-quadrant
    \param center center of the time band
    \param width width of the time band
  */
  PlaneRICH_MAPMT(const char *detname, int ncol, int nrow, int center, int width);
  virtual ~PlaneRICH_MAPMT();

#ifndef __CINT__
  void StoreDigit(CS::Chip::Digit* digit);
  void EndEvent(const CS::DaqEvent &event);  
#endif
  
  void Init(TTree* tree =0);

 friend class GroupRICH_MAPMT;

  ClassDef(PlaneRICH_MAPMT,1)
};

#endif


