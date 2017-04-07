/*
   $Source: /afs/cern.ch/user/g/gobbo/w0/coral/src/calorim/Reco/src/ReconstructionLednev.h,v $
   $Date: 2010/06/18 10:44:20 $
   $Revision: 1.4 $
*/

#ifndef ReconstructionLednev_Reco________include
#define ReconstructionLednev_Reco________include

#include <map>

#include "Reco_config.h"

class TH1F;
class TH2F;

#include "Reconstruction.h"

namespace Reco {

class ReconstructionLednev : public Reconstruction {
    // =========================================================================
    // Constructors and destructor
    // =========================================================================

    public:

        /// Default constructor
                                                        ReconstructionLednev   (const Calorimeter* c);

        /// Destructor
        virtual                                        ~ReconstructionLednev   (void) {}

    // ==========================================
    // Methods
    // ==========================================

    public:

        virtual const std::vector<CalorimeterParticle>& DoReconstruction       (const std::vector<CellDataRaw>& s);
        virtual const std::vector<CalorimeterParticle>& DoRepeatReconstruction (const std::vector<CellDataRaw>& s);

        virtual void                                    BookHistograms         ();

        virtual void                                    ReadCalibrations       ();

    // ==========================================
    //  Attributes, data
    // ==========================================

    private:

        /// pointers to histograms
        std::map<std::string, TH1F*>                    fHist1D;
        std::map<std::string, TH2F*>                    fHist2D;
};

} // namespace Reco

#endif // ReconstructionLednev_Reco________include
