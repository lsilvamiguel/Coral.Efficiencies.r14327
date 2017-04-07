/*
   $Source: /afs/cern.ch/user/g/gobbo/w0/coral/src/calorim/Reco/src/ReconstructionCombined.h,v $
   $Date: 2011/02/04 17:59:18 $
   $Revision: 1.28 $
*/

#ifndef ReconstructionCombined_Reco________include
#define ReconstructionCombined_Reco________include

#include <RVersion.h>
#if ROOT_VERSION_CODE >= ROOT_VERSION(5,22,00)

#include <limits>
#include <map>

#include <Math/IFunction.h>

#include "Reconstruction.h"

#include "CellDataRaw.h"
#include "CellType.h"

class TH1;
class TF3;

namespace Reco {

class Cluster;

class ReconstructionCombined : public Reconstruction {

    // ==========================================
    // Internal classes
    // ==========================================

    public:

        static const unsigned int nrFitParameters = 4;
        class FitShower {
            public:
                FitShower(const double& e,    const double& x,    const double& y,    const double& t,
                          const double& eErr, const double& xErr, const double& yErr, const double&tErr) :
                    fE(e), fEerr(eErr),
                    fX(x), fXerr(xErr),
                    fY(y), fYerr(yErr),
                    fT(t), fTerr(tErr),
                    fMainCell(std::numeric_limits<size_t>::max()) {}

                double  fE;
                double  fEerr;
                double  fX;
                double  fXerr;
                double  fY;
                double  fYerr;
                double  fT;
                double  fTerr;

                size_t  fMainCell;
        };

        class FitInfo {
            public:
                                        FitInfo(const Calorimeter* c, const CellDataRaw& s);
                const CellDataRaw&      GetData() const { return fData; }
                const CellType&         GetCellType() const;
                const ShowerProfile*    GetShowerProfile() const { return GetCellType().GetShowerProfile(); }
                void                    SetEnergy(double e);
                void                    SetEnergyErr(double err) { fData.SetEnergyErr(err); }
                void                    SetTimeErr(double err) { fData.SetTimeErr(err); }
                double                  GetFittedEnergy(const std::vector<FitShower>& shower) const;
                void                    SetArtificial(const bool& artificial);
                bool                    GetArtificial() const;
                void                    SetAround(const bool& around);
                bool                    GetAround() const;
            private:
                const Calorimeter*      fCalorimeter;
                CellDataRaw             fData;
                bool                    fArtificial;
                bool                    fAround;
        };

        class FitInput : public ROOT::Math::IMultiGradFunction {
            public:
                FitInput(const Calorimeter* c, const std::vector<FitInfo>& cells, const unsigned int& dimensions);

                double                  GetLogLh(const std::vector<FitShower>& shower)                  const;
                void                    GetTotalEnergy(double& eTot, double& eTotErr)                   const;
                void                    Print()                                                         const;

                // implementation of methods from ROOT::Math::IMultiGradFunction
                virtual FitInput*       Clone()                                                         const;
                virtual unsigned int    NDim()                                                          const;
                virtual double          DoEval(const double* par)                                       const;
                virtual double          DoDerivative(const double* par, unsigned int icoord)            const;
                virtual void            Gradient(const double* par, double* grad)                       const;
                virtual void            FdF(const double* par, double& f, double* grad)                 const;

            private:
                double                  CalcEnergyCorr(const double energy, const CellType& ct)         const;
                double                  CalcEnergyCorrDeriv(const double energy, const CellType& ct)    const;

                const Calorimeter*      fCalorimeter;
                std::vector<FitInfo>    fCells;
                unsigned int            fNDim;
        };

        class CorrectionPositionBin {
            public:
                                CorrectionPositionBin() {};
                virtual        ~CorrectionPositionBin() {};

                virtual double  GetCorrection(const double& dist) = 0;

                virtual void    Print(std::ostream& out) = 0;
        };

        class CorrectionPositionBinPol3 : public CorrectionPositionBin {
            public:
                                CorrectionPositionBinPol3(const double& a, const double& b, const double& c);

                double          GetCorrection(const double& dist);

                void            Print(std::ostream& out);

            private:
                double fA;
                double fB;
                double fC;
        };

        class CorrectionPositionBinSin : public CorrectionPositionBin {
            public:
                                CorrectionPositionBinSin(const double& a, const double& b);

                double          GetCorrection(const double& dist);

                void            Print(std::ostream& out);

            private:
                double fA;
                double fB;
        };

        class CorrectionPosition {
            public:
                                CorrectionPosition();
                virtual        ~CorrectionPosition();

                virtual double  GetCorrection(const double& energy, const double& dist);

                virtual void    AddCorrection(const double& minEnergy, CorrectionPositionBin* newBin);

                virtual void    Print(std::ostream& out, const std::string& title="");

                virtual void    SetSmoothLinear();

            private:
                std::map<double, CorrectionPositionBin*> fBins;
                bool                                     fSmoothLinear;
        };

        class CorrectionPositionSinAtan : public CorrectionPosition {
            public:
                                CorrectionPositionSinAtan(const double& a, const double& b, const double& c,
                                                          const double& d, const double& e, const double& f,
                                                          const double& g, const double& h, const double& i,
                                                          const double& j, const double& k);

                double          GetCorrection(const double& energy, const double& dist);

                void            AddCorrection(const double& minEnergy, CorrectionPositionBin* newBin);

                void            Print(std::ostream& out, const std::string& title="");

                void            SetSmoothLinear();

            private:
                double fA;
                double fB;
                double fC;
                double fD;
                double fE;
                double fF;
                double fG;
                double fH;
                double fI;
                double fJ;
                double fK;
        };

        class CorrectionPositionSinSin3Atan : public CorrectionPosition {
            public:
                                CorrectionPositionSinSin3Atan(const double& a, const double& b, const double& c,
                                                              const double& d, const double& e, const double& f,
                                                              const double& g,
                                                              const double& h, const double& i, const double& j,
                                                              const double& k, const double& l, const double& m,
                                                              const double& n, const double& o);

                double          GetCorrection(const double& energy, const double& dist);

                void            AddCorrection(const double& minEnergy, CorrectionPositionBin* newBin);

                void            Print(std::ostream& out, const std::string& title="");

                void            SetSmoothLinear();

            private:
                double fA;
                double fB;
                double fC;
                double fD;
                double fE;
                double fF;
                double fG;
                double fH;
                double fI;
                double fJ;
                double fK;
                double fL;
                double fM;
                double fN;
                double fO;
        };

        class CorrectionPosDepEnergyBin {
            public:
                                CorrectionPosDepEnergyBin() {};
                virtual        ~CorrectionPosDepEnergyBin() {};

                virtual double  GetCorrection(const double& energy,const double& distX, const double& distY) = 0;

                virtual void    Print(std::ostream& out) = 0;
        };

        class CorrectionPosDepEnergyBinQuarter : public CorrectionPosDepEnergyBin {
            public:
                                CorrectionPosDepEnergyBinQuarter(const int& size, const double& stepX, const double& stepY);

                virtual double  GetCorrection(const double& energy,const double& distX, const double& distY);

                void            Read(std::istream& in);
                virtual void    Print(std::ostream& out);

            private:
                int                     fSize;
                double                  fStepX;
                double                  fStepY;
                std::vector<double>     fCalib;
        };

        class CorrectionPosDepEnergyFunctionQuarter : public CorrectionPosDepEnergyBin {
            public:
                                CorrectionPosDepEnergyFunctionQuarter();

                virtual double  GetCorrection(const double& energy,const double& distX, const double& distY);

                void            Read(std::istream& in);
                virtual void    Print(std::ostream& out);

            private:
                TF3*                    fCalibFunction;
        };


        class CorrectionPosDepEnergy {
            public:
                                CorrectionPosDepEnergy();
                               ~CorrectionPosDepEnergy();

                double          GetCorrection(const double& energy, const double& distX, const double& distY);

                void            AddCorrection(const double& minEnergy, CorrectionPosDepEnergyBin* newBin);

                void            Print(std::ostream& out, const std::string& title="");

                void            SetSmoothLinear();

            private:
                std::map<double, CorrectionPosDepEnergyBin*>  fBins;
                bool                                          fSmoothLinear;
        };

    private:

        static const size_t nrCellTypesForHistograms = 10;
        enum HistogramIdentifier {
            h1D_cellEnergy,
            h1D_cellLowEnergy,
            h1D_cellEnergyThreshold,
            h1D_cellEnergyThresholdAdc,
            h1D_cellTime,
            h2D_cellEnergyVsTime,
            h2D_cellLowEnergyVsTime,

            h1D_clusterSize,
            h1D_clusterR4,
            h1D_clusterR9,

            h1D_fitLogLhDiff,
            h1D_fitLogLhDiffRel,

            h1D_fitDiffEnergyFirst,
            h2D_fitDiffEnergyFirstVsEnergy,
            h2D_fitDiffEnergyFirstVsEnergyLow,
            h1D_fitDiffPosXFirst,
            h1D_fitDiffPosYFirst,
            h1D_fitDiffTimeFirst,

            h1D_showerCount,
            h1D_showerDistance,
            h1D_showerEnergy,
            h1D_showerEnergyDesc,
            h1D_showerEnergyDescCells,
            h1D_showerLogLh,
            h1D_showerLogLhOverNdF,
            h1D_showerLogLhOverNdFLarge,
            h1D_showerNdF,
            h1D_showerProb,
            h1D_showerResiduumX,
            h1D_showerResiduumY,
            h1D_showerTime,

            h2D_showerCountVsClusterSize,
            h2D_showerDistance,
            h2D_showerDistanceVsEnergyRatio,
            h2D_showerEnergyDescVsEnergy,
            h2D_showerEnergyDescVsEnergyLow,
            h2D_showerEnergyDescCellsVsEnergy,
            h2D_showerEnergyDescCellsVsEnergyLow,
            h2D_showerLogLhOverNdFVsClusterSize,
            h2D_showerLogLhOverNdFVsShowerCount,
            h2D_showerPosition,
            h2D_showerResiduumXVsEnergy,
            h2D_showerResiduumYVsEnergy,

            h2D_showerEnergyVsShowerCount,
            h2D_showerEnergyVsClusterSize,
            h2D_showerEnergyLowVsShowerCount,
            h2D_showerEnergyLowVsClusterSize,

            h2D_showerEnergyVsTime,
            h2D_showerEnergyLowVsTime,

            h1D_showerTimeDiff,
            h2D_showerTimeDiffVsTime,

            h2D_fitPullsCellTypeFirst,
            h2D_fitPullsCellTypeLast = h2D_fitPullsCellTypeFirst + nrCellTypesForHistograms - 1,
            h2D_fitPullsLowCellTypeFirst,
            h2D_fitPullsLowCellTypeLast = h2D_fitPullsLowCellTypeFirst + nrCellTypesForHistograms - 1,
            h2D_fitNormPullsCellTypeFirst,
            h2D_fitNormPullsCellTypeLast = h2D_fitNormPullsCellTypeFirst + nrCellTypesForHistograms - 1,
            h2D_fitNormPullsLowCellTypeFirst,
            h2D_fitNormPullsLowCellTypeLast = h2D_fitNormPullsLowCellTypeFirst + nrCellTypesForHistograms - 1,

            NUMBER_OF_HISTOGRAMS                    // must be the last entry!
        };

    // ==========================================
    // Constructors and destructor
    // ==========================================

    public:

        /// Default constructor
                                                        ReconstructionCombined (const Calorimeter* c);

        /// Destructor
        virtual                                        ~ReconstructionCombined (void);

    // ==========================================
    // Methods
    // ==========================================

    public:

        /* \brief Reconstruct one calorimeter
         *
         * \param s individual information of hit cells
         *
         * \return Reconstructed particles
         */
        virtual const std::vector<CalorimeterParticle>& DoReconstruction       (const std::vector<CellDataRaw>& s);
        /* \brief NOT IMPLEMENTED
         */
        virtual const std::vector<CalorimeterParticle>& DoRepeatReconstruction (const std::vector<CellDataRaw>& s);

        /*! \brief Initialize the monitoring histograms
        */
        virtual void                                    BookHistograms         ();

        /*! \brief (Re-)read the calibrations
         */
        virtual void                                    ReadCalibrations       (void);

    private:

        /*! \brief Build clusters from the cell information
         *
         * \param s individual information of hit cells
         *
         * \return vector of clusters
         *
         * Build clusters from neighboring cells, for being neighboring, one of the
         * conditions has to be fulfilled:
         * - cells are really neighboring (trivial case)
         * - cells are next neighbors to the same dead cell
         * As the list of neighboring cells created in the Calorimeter class contains
         * also cells of different sizes, clusters may also contain cells of different
         * sizes.
         */
        virtual std::vector<Cluster>                    DoClustering           (const std::vector<CellDataRaw>& s);

        /*! \brief Fit shower parameters to data
         *
         * \param cells Measurement per cell
         * \param shower Parameters of the showers
         *
         * \return log-likelihood of the fit
         *
         * This method takes care of fitting the parameters of the current set of
         * showers to the data. The number of showers is not changed in this method!
         *
         * The information from shower provided when the method is called are used to
         * initialize the fit. In the end of this method the vector is cleared, and the
         * results from the fit are put into the vector.
         */
        virtual double                                  DoFitting              (const std::vector<FitInfo>& cells,
                                                                                std::vector<FitShower>& shower);

        virtual FitShower                               ApplyCorrections       (const FitShower& shower);

        double                                          CalcEnergyError        (const double& energy);
        double                                          CalcTimeError          (const double& energy);

    // ==========================================
    //  Attributes, data
    // ==========================================

    private:

        std::vector<unsigned int>                                               fAllowedNrShowers;
        std::vector<double>                                                     fParamEnergyErr;
        std::vector<double>                                                     fParamTimeErr;

        CorrectionPosition*                                                     fPosCorrX;
        CorrectionPosition*                                                     fPosCorrY;
        std::map<CellType::HwType, CorrectionPosDepEnergy*>                     fPosDepEnergyCorr;

        /// Histograms to monitor the reconstruction
        std::vector<TH1*>                                                       fHist;
        std::map<CellType::HwType, size_t>                                      fMapCellType2HistIndex;

};

} // namespace Reco

#else // ROOT version < 5.22

namespace Reco {

class ReconstructionCombined : public Reconstruction {
    public:
                                                        ReconstructionCombined(const Calorimeter* c) : Reconstruction(c) {
            std::cout << "Reco::ReconstructionCombined requires ROOT-Version >= 5.22/00!" << std::endl;
            exit(1);
        }
        virtual const std::vector<CalorimeterParticle>& DoReconstruction(const std::vector<CellDataRaw>&) {
            std::cout << "Reco::ReconstructionCombined requires ROOT-Version >= 5.22/00!" << std::endl;
            exit(1);
        }
        virtual const std::vector<CalorimeterParticle>& DoRepeatReconstruction(const std::vector<CellDataRaw>&) {
            std::cout << "Reco::ReconstructionCombined requires ROOT-Version >= 5.22/00!" << std::endl;
            exit(1);
        }
        virtual void                                    ReadCalibrations() {
            std::cout << "Reco::ReconstructionCombined requires ROOT-Version >= 5.22/00!" << std::endl;
            exit(1);
        }

};

} // namespace Reco

#endif // ROOT version >= 5.22

#endif // ReconstructionCombined_Reco________include

