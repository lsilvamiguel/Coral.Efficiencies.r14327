/*
   $Source: /afs/cern.ch/user/g/gobbo/w0/coral/src/calorim/Reco/src/ReconstructionCombinedFitInput.cc,v $
   $Date: 2011/02/05 15:46:56 $
   $Revision: 1.5 $
*/

#include <RVersion.h>
#if ROOT_VERSION_CODE >= ROOT_VERSION(5,22,00)

// --- Internal files ---
#include "ReconstructionCombined.h"

#include "Calorimeter.h"
#include "CellDataRaw.h"
#include "ShowerProfile.h"

namespace Reco {

/////////////////////////////////////////////////////////////////////////////

ReconstructionCombined::FitInput::FitInput(const Calorimeter* c, const std::vector<FitInfo>& cells, const unsigned int& dimensions) :
    fCalorimeter(c), fCells(cells), fNDim(dimensions) {
}

/////////////////////////////////////////////////////////////////////////////

unsigned int ReconstructionCombined::FitInput::NDim() const {
    return fNDim;
}

/////////////////////////////////////////////////////////////////////////////

ReconstructionCombined::FitInput* ReconstructionCombined::FitInput::Clone() const {
    return new ReconstructionCombined::FitInput(*this);
}

/////////////////////////////////////////////////////////////////////////////

/*! \brief Calculate the contribution of one shower to the log-likelihood
 */
double ReconstructionCombined::FitInput::DoEval(const double* par) const {
    double logLh(0.);
    for (std::vector<FitInfo>::const_iterator it=fCells.begin(); it!=fCells.end(); it++) {
        double eCalc(0.);
        double tCalc(0.);
        for (unsigned int i=0; i<fNDim/nrFitParameters; i++) {
            double cellE = par[nrFitParameters*i] * it->GetShowerProfile()->GetEnergyInCell(fCalorimeter->GetCells()[it->GetData().GetCellIdx()].GetX()-par[nrFitParameters*i+1],
                                                                                            fCalorimeter->GetCells()[it->GetData().GetCellIdx()].GetY()-par[nrFitParameters*i+2]);

            if (fCalorimeter->GetOptions().correct_longitudinal_leakage)
                cellE *= CalcEnergyCorr(par[nrFitParameters*i], fCalorimeter->GetCells()[it->GetData().GetCellIdx()].GetCellType());

            double cellT = par[nrFitParameters*i+3] * cellE;

            eCalc += cellE;
            tCalc += cellT;
        }
        tCalc /= eCalc;

        double deltaE=it->GetData().GetEnergy()-eCalc;
        double deltaT=it->GetData().GetTime()  -tCalc;

        double nE = it->GetData().GetEnergyErr()*it->GetData().GetEnergyErr();
        double nT = it->GetData().GetTimeErr()  *it->GetData().GetTimeErr();

        if (it->GetAround() && eCalc<it->GetData().GetEnergy())
                deltaE = 0.;

        logLh += 0.5 * ( (deltaE*deltaE)/nE + (deltaT*deltaT)/nT );
    }

    return logLh;
}

/////////////////////////////////////////////////////////////////////////////

void ReconstructionCombined::FitInput::Gradient(const double* par, double* grad) const {
    for (unsigned int i=0; i<fNDim; i++)
        grad[i] = 0.;

    for (std::vector<FitInfo>::const_iterator it=fCells.begin(); it!=fCells.end(); it++) {
        double eCalc(0.);
        double tCalc(0.);
        for (unsigned int i=0; i<fNDim/nrFitParameters; i++) {
            double cellE = par[nrFitParameters*i] * it->GetShowerProfile()->GetEnergyInCell(fCalorimeter->GetCells()[it->GetData().GetCellIdx()].GetX()-par[nrFitParameters*i+1],
                                                                                            fCalorimeter->GetCells()[it->GetData().GetCellIdx()].GetY()-par[nrFitParameters*i+2]);

            if (fCalorimeter->GetOptions().correct_longitudinal_leakage)
                cellE *= CalcEnergyCorr(par[nrFitParameters*i], fCalorimeter->GetCells()[it->GetData().GetCellIdx()].GetCellType());

            double cellT = par[nrFitParameters*i+3] * cellE;

            eCalc += cellE;
            tCalc += cellT;
        }
        tCalc /= eCalc;

        double deltaE=it->GetData().GetEnergy()-eCalc;
        double deltaT=it->GetData().GetTime()  -tCalc;

        double nE = it->GetData().GetEnergyErr()*it->GetData().GetEnergyErr();
        double nT = it->GetData().GetTimeErr()  *it->GetData().GetTimeErr();

        if (it->GetAround() && eCalc<it->GetData().GetEnergy())
            deltaE = 0.;

        for (unsigned int i=0; i<fNDim/nrFitParameters; i++) {
            double tSE = -it->GetShowerProfile()->GetEnergyInCell(fCalorimeter->GetCells()[it->GetData().GetCellIdx()].GetX()-par[nrFitParameters*i+1],
                                                                  fCalorimeter->GetCells()[it->GetData().GetCellIdx()].GetY()-par[nrFitParameters*i+2]);
            double tSX = par[nrFitParameters*i] * it->GetShowerProfile()->GetDeDxInCell(fCalorimeter->GetCells()[it->GetData().GetCellIdx()].GetX()-par[nrFitParameters*i+1],
                                                                                        fCalorimeter->GetCells()[it->GetData().GetCellIdx()].GetY()-par[nrFitParameters*i+2]);
            double tSY = par[nrFitParameters*i] * it->GetShowerProfile()->GetDeDyInCell(fCalorimeter->GetCells()[it->GetData().GetCellIdx()].GetX()-par[nrFitParameters*i+1],
                                                                                        fCalorimeter->GetCells()[it->GetData().GetCellIdx()].GetY()-par[nrFitParameters*i+2]);
            double tEC = par[nrFitParameters*i] * it->GetShowerProfile()->GetEnergyInCell(fCalorimeter->GetCells()[it->GetData().GetCellIdx()].GetX()-par[nrFitParameters*i+1],
                                                                                          fCalorimeter->GetCells()[it->GetData().GetCellIdx()].GetY()-par[nrFitParameters*i+2]);

            if (fCalorimeter->GetOptions().correct_longitudinal_leakage) {
                tSE *= CalcEnergyCorr(par[nrFitParameters*i], fCalorimeter->GetCells()[it->GetData().GetCellIdx()].GetCellType())
                       + par[nrFitParameters*i] * CalcEnergyCorrDeriv(par[nrFitParameters*i], fCalorimeter->GetCells()[it->GetData().GetCellIdx()].GetCellType());
                tSX *= CalcEnergyCorr(par[nrFitParameters*i], fCalorimeter->GetCells()[it->GetData().GetCellIdx()].GetCellType());
                tSY *= CalcEnergyCorr(par[nrFitParameters*i], fCalorimeter->GetCells()[it->GetData().GetCellIdx()].GetCellType());
                tEC *= CalcEnergyCorr(par[nrFitParameters*i], fCalorimeter->GetCells()[it->GetData().GetCellIdx()].GetCellType());
            }

            grad[nrFitParameters*i  ] += deltaE * tSE / nE + deltaT * tSE*(par[nrFitParameters*i+3]-tCalc)/eCalc / nT;
            grad[nrFitParameters*i+1] += deltaE * tSX / nE + deltaT * tSX*(par[nrFitParameters*i+3]-tCalc)/eCalc / nT;
            grad[nrFitParameters*i+2] += deltaE * tSY / nE + deltaT * tSY*(par[nrFitParameters*i+3]-tCalc)/eCalc / nT;
            grad[nrFitParameters*i+3] +=                   - deltaT *                   tEC               /eCalc / nT;
        }
    }
}

/////////////////////////////////////////////////////////////////////////////

double ReconstructionCombined::FitInput::DoDerivative(const double* par, unsigned int coord) const {
    double grad(0.);

    for (std::vector<FitInfo>::const_iterator it=fCells.begin(); it!=fCells.end(); it++) {
        double eCalc(0.);
        double tCalc(0.);
        for (unsigned int i=0; i<fNDim/nrFitParameters; i++) {
            double cellE = par[nrFitParameters*i] * it->GetShowerProfile()->GetEnergyInCell(fCalorimeter->GetCells()[it->GetData().GetCellIdx()].GetX()-par[nrFitParameters*i+1],
                                                                                            fCalorimeter->GetCells()[it->GetData().GetCellIdx()].GetY()-par[nrFitParameters*i+2]);

            if (fCalorimeter->GetOptions().correct_longitudinal_leakage)
                cellE *= CalcEnergyCorr(par[nrFitParameters*i], fCalorimeter->GetCells()[it->GetData().GetCellIdx()].GetCellType());

            double cellT = par[nrFitParameters*i+3] * cellE;

            eCalc += cellE;
            tCalc += cellT;
        }
        tCalc /= eCalc;

        double deltaE=it->GetData().GetEnergy()-eCalc;
        double deltaT=it->GetData().GetTime()  -tCalc;

        double nE = it->GetData().GetEnergyErr()*it->GetData().GetEnergyErr();
        double nT = it->GetData().GetTimeErr()  *it->GetData().GetTimeErr();

        if (it->GetAround() && eCalc<it->GetData().GetEnergy())
            deltaE = 0.;

        for (unsigned int i=0; i<fNDim/nrFitParameters; i++) {
            if (coord%nrFitParameters == 0) {
                double tSE = -it->GetShowerProfile()->GetEnergyInCell(fCalorimeter->GetCells()[it->GetData().GetCellIdx()].GetX()-par[nrFitParameters*i+1],
                                                                      fCalorimeter->GetCells()[it->GetData().GetCellIdx()].GetY()-par[nrFitParameters*i+2]);

                if (fCalorimeter->GetOptions().correct_longitudinal_leakage)
                    tSE *= CalcEnergyCorr(par[nrFitParameters*i], fCalorimeter->GetCells()[it->GetData().GetCellIdx()].GetCellType())
                           + par[nrFitParameters*i] * CalcEnergyCorrDeriv(par[nrFitParameters*i], fCalorimeter->GetCells()[it->GetData().GetCellIdx()].GetCellType());

                grad += deltaE * tSE / nE + deltaT * tSE*(par[nrFitParameters*i+3]-tCalc)/eCalc / nT;
            } else if (coord%nrFitParameters == 1) {
                double tSX = par[nrFitParameters*i] * it->GetShowerProfile()->GetDeDxInCell(fCalorimeter->GetCells()[it->GetData().GetCellIdx()].GetX()-par[nrFitParameters*i+1],
                                                                                            fCalorimeter->GetCells()[it->GetData().GetCellIdx()].GetY()-par[nrFitParameters*i+2]);

                if (fCalorimeter->GetOptions().correct_longitudinal_leakage)
                    tSX *= CalcEnergyCorr(par[nrFitParameters*i], fCalorimeter->GetCells()[it->GetData().GetCellIdx()].GetCellType());

                grad += deltaE * tSX / nE + deltaT * tSX*(par[nrFitParameters*i+3]-tCalc)/eCalc / nT;
            } else if (coord%nrFitParameters == 2) {
                double tSY = par[nrFitParameters*i] * it->GetShowerProfile()->GetDeDyInCell(fCalorimeter->GetCells()[it->GetData().GetCellIdx()].GetX()-par[nrFitParameters*i+1],
                                                                                            fCalorimeter->GetCells()[it->GetData().GetCellIdx()].GetY()-par[nrFitParameters*i+2]);

                if (fCalorimeter->GetOptions().correct_longitudinal_leakage)
                    tSY *= CalcEnergyCorr(par[nrFitParameters*i], fCalorimeter->GetCells()[it->GetData().GetCellIdx()].GetCellType());

                grad += deltaE * tSY / nE + deltaT * tSY*(par[nrFitParameters*i+3]-tCalc)/eCalc / nT;
            } else if (coord%nrFitParameters == 3) {
                double tSE = -it->GetShowerProfile()->GetEnergyInCell(fCalorimeter->GetCells()[it->GetData().GetCellIdx()].GetX()-par[nrFitParameters*i+1],
                                                                      fCalorimeter->GetCells()[it->GetData().GetCellIdx()].GetY()-par[nrFitParameters*i+2]);
                double tEC = par[nrFitParameters*i] * it->GetShowerProfile()->GetEnergyInCell(fCalorimeter->GetCells()[it->GetData().GetCellIdx()].GetX()-par[nrFitParameters*i+1],
                                                                                              fCalorimeter->GetCells()[it->GetData().GetCellIdx()].GetY()-par[nrFitParameters*i+2]);

                if (fCalorimeter->GetOptions().correct_longitudinal_leakage) {
                    tSE *= CalcEnergyCorr(par[nrFitParameters*i], fCalorimeter->GetCells()[it->GetData().GetCellIdx()].GetCellType())
                           + par[nrFitParameters*i] * CalcEnergyCorrDeriv(par[nrFitParameters*i], fCalorimeter->GetCells()[it->GetData().GetCellIdx()].GetCellType());
                    tEC *= CalcEnergyCorr(par[nrFitParameters*i], fCalorimeter->GetCells()[it->GetData().GetCellIdx()].GetCellType());
                }

                grad +=                   - deltaT *                   tEC               /eCalc / nT;
            }
        }
    }

    return grad;
}

/////////////////////////////////////////////////////////////////////////////

void ReconstructionCombined::FitInput::FdF(const double* par, double& f, double* grad) const {
    f=0.;
    for (unsigned int i=0; i<fNDim; i++)
        grad[i] = 0.;

    for (std::vector<FitInfo>::const_iterator it=fCells.begin(); it!=fCells.end(); it++) {
        double eCalc(0.);
        double tCalc(0.);
        for (unsigned int i=0; i<fNDim/nrFitParameters; i++) {
            double cellE = par[nrFitParameters*i] * it->GetShowerProfile()->GetEnergyInCell(fCalorimeter->GetCells()[it->GetData().GetCellIdx()].GetX()-par[nrFitParameters*i+1],
                                                                                            fCalorimeter->GetCells()[it->GetData().GetCellIdx()].GetY()-par[nrFitParameters*i+2]);

            if (fCalorimeter->GetOptions().correct_longitudinal_leakage)
                cellE *= CalcEnergyCorr(par[nrFitParameters*i], fCalorimeter->GetCells()[it->GetData().GetCellIdx()].GetCellType());

            double cellT = par[nrFitParameters*i+3] * cellE;

            eCalc += cellE;
            tCalc += cellT;
        }
        tCalc /= eCalc;

        double deltaE=it->GetData().GetEnergy()-eCalc;
        double deltaT=it->GetData().GetTime()  -tCalc;

        double nE = it->GetData().GetEnergyErr()*it->GetData().GetEnergyErr();
        double nT = it->GetData().GetTimeErr()  *it->GetData().GetTimeErr();

        if (it->GetAround() && eCalc<it->GetData().GetEnergy())
            deltaE = 0.;

        f += 0.5 * ( (deltaE*deltaE)/nE + (deltaT*deltaT)/nT );

        for (unsigned int i=0; i<fNDim/nrFitParameters; i++) {
            double tSE = -it->GetShowerProfile()->GetEnergyInCell(fCalorimeter->GetCells()[it->GetData().GetCellIdx()].GetX()-par[nrFitParameters*i+1],
                                                                  fCalorimeter->GetCells()[it->GetData().GetCellIdx()].GetY()-par[nrFitParameters*i+2]);
            double tSX = par[nrFitParameters*i] * it->GetShowerProfile()->GetDeDxInCell(fCalorimeter->GetCells()[it->GetData().GetCellIdx()].GetX()-par[nrFitParameters*i+1],
                                                                                        fCalorimeter->GetCells()[it->GetData().GetCellIdx()].GetY()-par[nrFitParameters*i+2]);
            double tSY = par[nrFitParameters*i] * it->GetShowerProfile()->GetDeDyInCell(fCalorimeter->GetCells()[it->GetData().GetCellIdx()].GetX()-par[nrFitParameters*i+1],
                                                                                        fCalorimeter->GetCells()[it->GetData().GetCellIdx()].GetY()-par[nrFitParameters*i+2]);
            double tEC = par[nrFitParameters*i] * it->GetShowerProfile()->GetEnergyInCell(fCalorimeter->GetCells()[it->GetData().GetCellIdx()].GetX()-par[nrFitParameters*i+1],
                                                                                          fCalorimeter->GetCells()[it->GetData().GetCellIdx()].GetY()-par[nrFitParameters*i+2]);

            if (fCalorimeter->GetOptions().correct_longitudinal_leakage) {
                tSE *= CalcEnergyCorr(par[nrFitParameters*i], fCalorimeter->GetCells()[it->GetData().GetCellIdx()].GetCellType())
                       + par[nrFitParameters*i] * CalcEnergyCorrDeriv(par[nrFitParameters*i], fCalorimeter->GetCells()[it->GetData().GetCellIdx()].GetCellType());
                tSX *= CalcEnergyCorr(par[nrFitParameters*i], fCalorimeter->GetCells()[it->GetData().GetCellIdx()].GetCellType());
                tSY *= CalcEnergyCorr(par[nrFitParameters*i], fCalorimeter->GetCells()[it->GetData().GetCellIdx()].GetCellType());
                tEC *= CalcEnergyCorr(par[nrFitParameters*i], fCalorimeter->GetCells()[it->GetData().GetCellIdx()].GetCellType());
            }

            grad[nrFitParameters*i  ] += deltaE * tSE / nE + deltaT * tSE*(par[nrFitParameters*i+3]-tCalc)/eCalc / nT;
            grad[nrFitParameters*i+1] += deltaE * tSX / nE + deltaT * tSX*(par[nrFitParameters*i+3]-tCalc)/eCalc / nT;
            grad[nrFitParameters*i+2] += deltaE * tSY / nE + deltaT * tSY*(par[nrFitParameters*i+3]-tCalc)/eCalc / nT;
            grad[nrFitParameters*i+3] +=                   - deltaT *                   tEC               /eCalc / nT;
        }
    }
}

/////////////////////////////////////////////////////////////////////////////

double ReconstructionCombined::FitInput::GetLogLh(const std::vector<FitShower>& shower) const {
    assert(fNDim==nrFitParameters*shower.size());
    double* par = new double[fNDim];
    for (size_t i=0; i<shower.size(); i++) {
        par[nrFitParameters*i]   = shower[i].fE;
        par[nrFitParameters*i+1] = shower[i].fX;
        par[nrFitParameters*i+2] = shower[i].fY;
        par[nrFitParameters*i+3] = shower[i].fT;
    }

    double logLh(DoEval(par));

    delete [] par;

    return logLh;
}

/////////////////////////////////////////////////////////////////////////////

void ReconstructionCombined::FitInput::GetTotalEnergy(double& eTot, double& eTotErr) const {
    eTot    = 0.;
    eTotErr = 0.;
    for (std::vector<FitInfo>::const_iterator it=fCells.begin(); it!=fCells.end(); it++) {
        eTot    += it->GetData().GetEnergy();
        eTotErr += it->GetData().GetEnergyErr() * it->GetData().GetEnergyErr();
    }
    eTotErr = sqrt(eTotErr);
}

/////////////////////////////////////////////////////////////////////////////

void ReconstructionCombined::FitInput::Print() const {
    std::cout << std::endl << "cluster" << std::endl;
    std::cout.width(10);
    std::cout.precision(10);
    for (std::vector<FitInfo>::const_iterator it=fCells.begin(); it!=fCells.end(); it++) {
        std::cout << it->GetData().GetCellIdx() << " " << it->GetData().GetEnergy() << " " << it->GetData().GetEnergyErr() << " "
                  << fCalorimeter->GetCells()[it->GetData().GetCellIdx()].GetX() << " " << fCalorimeter->GetCells()[it->GetData().GetCellIdx()].GetY()
                  << std::endl;
    }
}

/////////////////////////////////////////////////////////////////////////////

double ReconstructionCombined::FitInput::CalcEnergyCorr(const double energy, const CellType& ct) const {
    // no correction if energy is below the critical energy, for
    // 1. this should obviously have no effect, if we cannot even contain
    //    15 MeV photons and electrons/positrons in the calorimeter we are
    //    doomed
    // 2. technical reasons, the log of a number smaller than 1 tends to
    //    quickly approach minus infinity
    if (energy <= ct.GetCriticalEnergy())
        return 1.;

    // using the integral of Formula 27.33 of the PDG 2010 (page 294)
    const double b =  0.5;
    const double t = ct.GetTrueSizeZ() / ct.GetRadiationLength();
    const double y = energy / ct.GetCriticalEnergy();

    double       C =  0.;
    if      ( fCalorimeter->GetOptions().particle_default == CalorimeterParticle::GAMMA )
        C = +0.5;
    else if ( fCalorimeter->GetOptions().particle_default == CalorimeterParticle::POSITRON ||
              fCalorimeter->GetOptions().particle_default == CalorimeterParticle::ELECTRON )
        C = -0.5;
    else {
        std::cerr << "ReconstructionCombined::FitInput::CalcEnergyCorr: Particle with partid " << fCalorimeter->GetOptions().particle_default << " not implemented." << std::endl;
        exit(1);
        // TODO should be fatal exception
        throw Exception("ReconstructionCombined::FitInput::CalcEnergyCorr: Particle with partid %i not implemented.", fCalorimeter->GetOptions().particle_default);
    }

    return TMath::Gamma(b * (log(y) + C) + 1., b * t);
}

/////////////////////////////////////////////////////////////////////////////

double ReconstructionCombined::FitInput::CalcEnergyCorrDeriv(const double energy, const CellType& ct) const {
    // using Richardson's extrapolation
    // ROOT uses 0.001*range as eps, for the range of 0 to 200 GeV, that would
    // be 0.2, I use a value slightly smaller
    const double eps     = 0.05;
    const double epsHalf = eps * 0.5;

    const double DHalf = (CalcEnergyCorr(energy+epsHalf, ct) - CalcEnergyCorr(energy-epsHalf, ct)) / (2.*epsHalf);
    const double D     = (CalcEnergyCorr(energy+eps,     ct) - CalcEnergyCorr(energy-eps,     ct)) / (2.*eps);

    return (4. * DHalf - D) / 3.;
}

/////////////////////////////////////////////////////////////////////////////

} // namespace Reco

#endif // ROOT version >= 5.22

