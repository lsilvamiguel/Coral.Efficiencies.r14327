#ifndef __GroupStraws__
#define __GroupStraws__

#include "Group.h"

namespace CS { class DaqEvent; }

class TH1F;
class TH2F;

class GroupStraws : public Group
{
  public:
                        GroupStraws             (const char* name);
    void                Init                    (void);

    #ifndef __CINT__
    /// Fills the histograms
    void                EndEvent                (const CS::DaqEvent &event);
    #endif

    TH1F*               fChannels;              // hits profile
    TH1F*               fRates;                 // rates profile
    TH1F*               fHits;                  // hits multiplicity
    TH1F*               fTime;                  // 
    TH1F*               fTime6mm;               // 
    TH1F*               fTime10mm;              // 
    
  ClassDef              (GroupStraws,1)
};

#endif
