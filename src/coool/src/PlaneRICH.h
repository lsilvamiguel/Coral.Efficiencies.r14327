#ifndef __PlaneRICH__
#define __PlaneRICH__

#include "Plane.h"
#include "PlanePanel.h"
#include "PlaneRCath.h"


#include "TThread.h"

 

/*! \brief class for RICH : (cathode, x, y, amp) digit type

  \todo test with data
  \todo fill up Variable vector with PlaneRCath variables
  \todo Control Panel -> cathodes control panels

  \author Colin Bernet
*/

class PlaneRICH : public Plane {
  
 private:

  /// vector of photocathodes
  std::vector<PlaneRCath*> fCathode;

  TH2F *fH2u;
  TH2F *fH2l;
  TH1F *fH1s;
  TH1F *fH1h;

 public:
  /*! \brief constructor
    \param detname detector name
    \param ncath number of cathodes   
    \param ncols number of columns per cathode
    \param nrows number of rows per cathode
    \param center center of the amplitude band
    \param width width of the amplitude band
  */
  PlaneRICH(const char *detname,size_t ncath, int ncols, int nrows, int center, int width);
  virtual ~PlaneRICH();

  void Reset();
  void ResetHistograms();

  /// Passes a digit to the plane
  //  void StoreDigit(int cath, int col, int row, int amp);

#ifndef __CINT__
  void StoreDigit(CS::Chip::Digit* digit);
  void EndEvent(const CS::DaqEvent &event);  
#endif

  void TextOutput(ostream& out);
  
  void Init(TTree* tree =0);
  void ControlPanel(const TGWindow *p, const TGWindow *main);

  ClassDef(PlaneRICH,1)
};

#endif


