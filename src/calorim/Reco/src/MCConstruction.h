/*
   $Source: /afs/cern.ch/user/g/gobbo/w0/coral/src/calorim/Reco/src/MCConstruction.h,v $
   $Date: 2010/03/31 14:50:08 $
   $Revision: 1.1 $
*/

#ifndef MCConstruction_Reco________include
#define MCConstruction_Reco________include

#include <vector>

namespace Reco {

class Calorimeter;
class CalorimeterParticle;
class CellDataRaw;
class OneParticleResponse;

class MCConstruction {
    // =========================================================================
    // Constructors and destructor
    // =========================================================================

    public:

        /// Default constructor
                                                    MCConstruction              (const Calorimeter* c);

        /// Destructor
        virtual                                    ~MCConstruction              (void) {}

    // ==========================================
    // Methods
    // ==========================================

    public:

        virtual const std::vector<CellDataRaw>&     DoDigitization              (const std::vector<CalorimeterParticle>& particles);

        virtual const std::vector<CellDataRaw>&     DoConstruction              (const std::vector<CellDataRaw>& mc_input);

        void                                        ClearExpected               ();

        double                                      AddParticle                 (const CalorimeterParticle& particle);

        const Calorimeter*                          GetCalorimeter              (void) const { return calorimeter_; }

    private:

        void                                        AddToExpected               (const OneParticleResponse& response,
                                                                                 const double& time);

        void                                        ExtractExpected             (const double& extraction_energy_threshold);

        void                                        ImportExpected              (const std::vector<CellDataRaw>& mc_input);

        void                                        FilterOutBadCells           ();

        void                                        AddConstantFluctuations     ();
        void                                        AddStochasticFluctuations   ();
        void                                        AddReadoutFluctuations      ();

    // ==========================================
    //  Attributes, data
    // ==========================================

    private:

        const Calorimeter*                          calorimeter_;

        std::vector<CellDataRaw>                    signals_;

        /// Energy sum per cell expected from particles
        std::vector<std::pair<double, double> >     expected_amplitudes;
        std::vector<double>                         expected_times;

};

} // namespace Reco

#endif // MCConstruction_Reco________include

