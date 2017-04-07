#ifndef CS_Trigger__include
#define CS_Trigger__include

#include <string>
#include "ObjectXML.h"

namespace CS {

/*! \brief The trigger settings.
*/
class Trigger
{
  public:

                        Trigger                 (const ObjectXML &xml=ObjectXML("Trigger"));

  public:

    void                Print                   (const std::string &prefix="") const;

    const std::string&  GetName                 (void) const {return name;}
    void                SetName                 (const std::string &s) {name=s;}

    unsigned            GetBit                  (void) const {return bit;}
    void                SetBit                  (unsigned b) {bit=b;}
   
    float               GetCorrection           (void) const {return mt_correction;}
    void                SetCorrection           (float f) {mt_correction=f;}

    float               GetOffset               (void) const {return mt_offset;}

    float               GetPrecision            (void) const {return precision;}

    bool                IsVerbose               (void) const {return verbose;}

  private:

    unsigned            bit;

    std::string         name;

    // precision in [ns]
    float               precision;

    // Correction to be added to the master time, if this trigger was fired. [ns]
    float               mt_correction;

    // Offset trigger_time-master_time [counts]
    float               mt_offset;

    bool                verbose;
};

////////////////////////////////////////////////////////////////

} // namespace

#endif // CS_Trigger__include
