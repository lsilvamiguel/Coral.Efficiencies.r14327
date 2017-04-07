#ifndef __PlaneTcsPhase__
#define __PlaneTcsPhase__

#include "Plane1V.h"

#include <vector>

class PlaneTcsPhase : public  Plane1V {
    public:

        PlaneTcsPhase(const char *detname);
       ~PlaneTcsPhase() {};
    
        void Init(TTree *tree);
        void Reset();

#ifndef __CINT__
        void EndEvent(const CS::DaqEvent &event);
#endif // __CINT__

    private:
#ifndef __CINT__
        /// get a TCS phase for each digit in the TCS detector
        std::vector<float> GetTCSPhase(const CS::DaqEvent &event) const;
#endif // __CINT__

        TH1F* fTCSPhase;
        TH1F* fTCSPhaseDDD;

        float fTCSPhaseVal;

    ClassDef(PlaneTcsPhase,1)
};

#endif


















