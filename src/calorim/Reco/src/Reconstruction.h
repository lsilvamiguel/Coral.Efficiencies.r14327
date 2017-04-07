/*
   $Source: /afs/cern.ch/user/g/gobbo/w0/coral/src/calorim/Reco/src/Reconstruction.h,v $
   $Date: 2011/02/01 22:05:52 $
   $Revision: 1.4 $
*/

#ifndef Reconstruction_Reco________include
#define Reconstruction_Reco________include

#include <vector>

#include "CalorimeterParticle.h"

namespace Reco {

class Calorimeter;
class CellDataRaw;

class Reconstruction {
    // =========================================================================
    // Constructors and destructor
    // =========================================================================

    public:

        /// Default constructor
                                                        Reconstruction         (const Calorimeter* c);

        /// Destructor
        virtual                                        ~Reconstruction         (void) {}

    // ==========================================
    // Methods
    // ==========================================

    public:

        virtual const std::vector<CalorimeterParticle>& DoReconstruction       (const std::vector<CellDataRaw>& s) = 0;
        virtual const std::vector<CalorimeterParticle>& DoRepeatReconstruction (const std::vector<CellDataRaw>& s) = 0;

        virtual void                                    BookHistograms         (void) {}

	/// Read some variables from Calorimeter object into Reconstruction
	/// instance.  Only used in ReconstructionLednev for the moment.
        virtual void                                    ReadCalibrations       (void) = 0;

        const Calorimeter*                              GetCalorimeter         (void) const { return calorimeter_; }

    // ==========================================
    //  Attributes, data
    // ==========================================

    protected:

        std::vector<CalorimeterParticle>    reconstructedParticles;

    private:

        const Calorimeter*                  calorimeter_;
};

} // namespace Reco

#endif // Reconstruction_Reco________include

