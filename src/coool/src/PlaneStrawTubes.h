#ifndef __PlaneStrawTubes__
#define __PlaneStrawTubes__

#include "Plane1V.h"
#include "TTree.h"

/*! \brief Plane for Straw Tubes detectors
    \todo solve the problem in Monitor
    
    \author Alexander Zvyagin
*/

class PlaneStrawTubes : public  Plane1V 
{
  public:
  
    virtual            ~PlaneStrawTubes         (void);
  
                        PlaneStrawTubes         (const char *detname,int center, int width);
    void                Init                    (TTree* tree = 0);

  private:
  

    #ifndef __CINT__
    void                EndEvent                (const CS::DaqEvent &event);
    #endif

  ClassDef(PlaneStrawTubes,2)
};

#endif
