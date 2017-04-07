#ifndef __PlaneRichWall__
#define __PlaneRichWall__

#include "Plane1V.h"

/// Plane for drift chambers
class PlaneRichWall : public Plane1V 
{
 
  int   fPlaneNumber;
private:

  int  Nevents;
  TH1F_Ref  *fHrates2;
  TH1F_Ref  *fHch2;
  TH1F  *fHnoise_percent;
  TH1F_Ref  *fHtime_rev;

public:
  
  PlaneRichWall( const char *detname,int nchan, int center, int width );

  void Init( TTree* tree =0 );
  
#ifndef __CINT__
  void EndEvent(const CS::DaqEvent &event );
#endif

  void Reset( void );


 void ResetHistograms();


  ClassDef(PlaneRichWall,1)
};

#endif
