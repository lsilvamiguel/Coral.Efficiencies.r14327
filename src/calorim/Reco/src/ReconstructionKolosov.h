/*
   $Source: /afs/cern.ch/user/g/gobbo/w0/coral/src/calorim/Reco/src/ReconstructionKolosov.h,v $
   $Date: 2010/04/15 14:20:02 $
   $Revision: 1.3 $
*/

#ifndef ReconstructionKolosov_Reco________include
#define ReconstructionKolosov_Reco________include

#include "Reconstruction.h"

#include "Calorimeter.h"
#include "OneParticleResponse.h"

namespace Reco {

class ReconstructionKolosov : public Reconstruction {
    // ==========================================
    //  subclasses
    // ==========================================

    // =========================================================================
    // Constructors and destructor
    // =========================================================================

    public:

        /// Default constructor
                                                        ReconstructionKolosov   (const Calorimeter* c);

        /// Destructor
        virtual                                        ~ReconstructionKolosov   (void) {}

    // ==========================================
    // Methods
    // ==========================================

    public:

        virtual const std::vector<CalorimeterParticle>& DoReconstruction        (const std::vector<CellDataRaw>& s);
        virtual const std::vector<CalorimeterParticle>& DoRepeatReconstruction  (const std::vector<CellDataRaw>& s);

        virtual void                                    ReadCalibrations        (void);

    private:

        virtual void                                    ClearArray              (void);
        virtual void                                    ClearExpectedAmplitudes (void);
        virtual void                                    ClearRealAmplitudes     (void);

        virtual void                                    InsertData              (const std::vector<CellDataRaw>& s);

        virtual void                                    InsertBadCellsData      ();

        virtual void                                    ManyParticlesResponse   (void);
        virtual void                                    ManyParticlesShareRealResponse
                                                                                (void);
        virtual void                                    RecoverEnergyDelivery   (const std::vector<size_t>& mp,
                                                                                 const std::vector<size_t>& vcells);

        virtual void                                    AddToExpectedAmplitudes (const OneParticleResponse& partResponse);
        virtual void                                    AddToRealAmplitudes     (const OneParticleResponse& partResponse);

        virtual double                                  GetExpectedAmplitudesSum
                                                                                (void) const;
        virtual double                                  GetRealAmplitudesSum    (void) const;

        virtual bool                                    CompareResponse         (void);

        virtual void                                    ReturnResult            (void);

        virtual void                                    AddCellToCluster        (const size_t& cell_number,
                                                                                 Cluster& cluster,
                                                                                 std::vector<size_t>& in,
                                                                                 const double& cell_porog,
                                                                                 const double& Nstand) const;
        virtual void                                    FindClusters            (std::vector<Cluster> &clusters,
                                                                                 const double& cell_porog = 0.02,
                                                                                 const double& Nstand     = 2.) const;

        virtual void                                    FindGamma0              (const double& GammaThreshold);
        virtual bool                                    IsNicePeak0             (const size_t& cell_number,
                                                                                 const double& Nstand,
                                                                                 const double& CellThreshold) const;


        virtual void                                    ApplyGammaEnergyCuts    (void);
        virtual bool                                    GammaEnergyCut          (const size_t& incell,
                                                                                 const double& x,
                                                                                 const double& y,
                                                                                 const double& e) const;

        virtual void                                    AddTimeToRecoParticles  (const std::vector<CellDataRaw>& s);

        virtual bool                                    CheckPointInCell        (const double& x, const double& y, const size_t& cell) const;

        virtual void                                    Fit                     (const Calorimeter::FitMethod& fitting_type,
                                                                                 OneParticleResponse& partResponse);
        virtual void                                    Fit                     (OneParticleResponse& partResponse);
        virtual void                                    NoFit                   (OneParticleResponse& partResponse);
        virtual void                                    SimpleFit               (OneParticleResponse& partResponse);

    // ==========================================
    //  Attributes, data
    // ==========================================

    private:

        /// Reconstructed particle response
        std::vector<OneParticleResponse>                 manyparticles;

        /// Energy deposit in calorimeter cells
        std::vector<std::pair<double, double> >          real_amplitudes;

        /// Energy sum from all hypothesis
        std::vector<std::pair<double, double> >          expected_amplitudes;

        std::vector<std::pair<double, double> >          array_tmp_;

};

} // namespace Reco

#endif // ReconstructionKolosov_Reco________include

