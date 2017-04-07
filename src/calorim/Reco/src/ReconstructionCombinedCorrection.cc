/*
   $Source: /afs/cern.ch/user/g/gobbo/w0/coral/src/calorim/Reco/src/ReconstructionCombinedCorrection.cc,v $
   $Date: 2011/02/04 18:15:57 $
   $Revision: 1.2 $
*/

#include <RVersion.h>
#if ROOT_VERSION_CODE >= ROOT_VERSION(5,22,00)
#include "TF3.h"

#include <cmath>

// --- Internal files ---
#include "ReconstructionCombined.h"

#include "Calorimeter.h"

namespace Reco {

/////////////////////////////////////////////////////////////////////////////

ReconstructionCombined::CorrectionPosition::CorrectionPosition() : fSmoothLinear(false) {
}

/////////////////////////////////////////////////////////////////////////////

ReconstructionCombined::CorrectionPosition::~CorrectionPosition() {
    for (std::map<double, CorrectionPositionBin*>::const_iterator it=fBins.begin(); it!=fBins.end(); it++)
        delete it->second;
}

/////////////////////////////////////////////////////////////////////////////

double ReconstructionCombined::CorrectionPosition::GetCorrection(const double& energy, const double& dist) {
    if (fBins.size()==0)
        return 0.;

    if (fSmoothLinear) {
        std::map<double, CorrectionPositionBin*>::const_iterator it;
        if ( (it=fBins.lower_bound(energy))!=fBins.begin() ) {
            std::map<double, CorrectionPositionBin*>::const_iterator itLower=it; itLower--;

            if (it==fBins.end())
                return itLower->second->GetCorrection(dist);

            const double x1 = itLower->first;
            const double y1 = itLower->second->GetCorrection(dist);
            const double x2 = it     ->first;
            const double y2 = it     ->second->GetCorrection(dist);
            return y1 + (energy-x1) * (y2-y1)/(x2-x1);
        } else
            return it->second->GetCorrection(dist);
    } else {
        std::map<double, CorrectionPositionBin*>::const_iterator it;
        if ( (it=fBins.lower_bound(energy))!=fBins.begin() ) {
            it--;
            return it->second->GetCorrection(dist);
        }
    }

    return 0.;
}

/////////////////////////////////////////////////////////////////////////////

void ReconstructionCombined::CorrectionPosition::AddCorrection(const double& minEnergy, CorrectionPositionBin* newBin) {
    fBins[minEnergy] = newBin;
}

/////////////////////////////////////////////////////////////////////////////

void ReconstructionCombined::CorrectionPosition::Print(std::ostream& out, const std::string& title) {
    out << title;
    if (fBins.size()==0)
        out << "[none]" << std::endl;
    else {
        bool first(true);
        if (fSmoothLinear) {
            out << "with linear smoothing" << std::endl;
            first = false;
        }
        double lastEnergy(0);
        CorrectionPositionBin* lastPosCorr(0);
        for (std::map<double, CorrectionPositionBin*>::const_iterator it=fBins.begin(); it!=fBins.end(); it++) {
            if (lastPosCorr!=0) {
                if (!first)
                    for (size_t i=0; i<title.length(); i++)
                        out << " ";
                first = false;

                if (fSmoothLinear)
                    out << "at E=" << lastEnergy << ": ";
                else
                    out << lastEnergy << "<=E<" << it->first << ": ";
                lastPosCorr->Print(out);
                out << std::endl;
            }

            lastEnergy  = it->first;
            lastPosCorr = it->second;
        }
        if (!first)
            for (size_t i=0; i<title.length(); i++)
                out << " ";
        if (fSmoothLinear)
            out << "at E=" << lastEnergy << ": ";
        else
            out << "E>=" << lastEnergy << ": ";
        lastPosCorr->Print(out);
        out << std::endl;
    }
}

/////////////////////////////////////////////////////////////////////////////

void ReconstructionCombined::CorrectionPosition::SetSmoothLinear() {
    fSmoothLinear = true;
}

/////////////////////////////////////////////////////////////////////////////

ReconstructionCombined::CorrectionPositionSinAtan::CorrectionPositionSinAtan(const double& a, const double& b, const double& c,
                                                                             const double& d, const double& e, const double& f,
                                                                             const double& g, const double& h, const double& i,
                                                                             const double& j, const double& k) :
    fA(a), fB(b), fC(c), fD(d), fE(e), fF(f), fG(g), fH(h), fI(i), fJ(j), fK(k) {
}

/////////////////////////////////////////////////////////////////////////////

double ReconstructionCombined::CorrectionPositionSinAtan::GetCorrection(const double& energy, const double& dist) {
    return (  fA*atan((energy-fB)/fC)
            + fD*atan((energy-fE)/fF)
            + fG*atan((energy-fH)/fI) + fJ) * sin(fK*dist);

}

/////////////////////////////////////////////////////////////////////////////

void ReconstructionCombined::CorrectionPositionSinAtan::AddCorrection(const double& minEnergy, CorrectionPositionBin* newBin) {
    std::cerr << "ReconstructionCombined::CorrectionPositionSinAtan::AddCorrection: Position correction of this kind does not use individual energy bins!" << std::endl;
    exit(1);
    // TODO: should be fatal exception, but is caught somewhere
    throw Exception("ReconstructionCombined::CorrectionPositionSinAtan::AddCorrection: Position correction of this kind does not use individual energy bins!");
}

/////////////////////////////////////////////////////////////////////////////

void ReconstructionCombined::CorrectionPositionSinAtan::Print(std::ostream& out, const std::string& title) {
    out << title;
    out << "SinAtan: a=" << fA
               << ", b=" << fB
               << ", c=" << fC
               << ", d=" << fD
               << ", e=" << fE
               << ", f=" << fF
               << ", g=" << fG
               << ", h=" << fH
               << ", i=" << fI
               << ", j=" << fJ
               << ", k=" << fK << std::endl;
}

/////////////////////////////////////////////////////////////////////////////

void ReconstructionCombined::CorrectionPositionSinAtan::SetSmoothLinear() {
    std::cerr << "ReconstructionCombined::CorrectionPositionSinAtan::SetSmoothLinear: SinAtan does not need this!" << std::endl;
    exit(1);
    // TODO: should be fatal exception, but is caught somewhere
    throw Exception("ReconstructionCombined::CorrectionPositionSinAtan::SetSmoothLinear: SinAtan does not need this!");
}

/////////////////////////////////////////////////////////////////////////////
//
ReconstructionCombined::CorrectionPositionSinSin3Atan::CorrectionPositionSinSin3Atan(const double& a, const double& b, const double& c,
                                                                                     const double& d, const double& e, const double& f,
                                                                                     const double& g,
                                                                                     const double& h, const double& i, const double& j,
                                                                                     const double& k, const double& l, const double& m,
                                                                                     const double& n, const double& o) :
    fA(a), fB(b), fC(c), fD(d), fE(e), fF(f), fG(g), fH(h), fI(i), fJ(j), fK(k), fL(l), fM(m), fN(n), fO(o) {
}

/////////////////////////////////////////////////////////////////////////////

double ReconstructionCombined::CorrectionPositionSinSin3Atan::GetCorrection(const double& energy, const double& dist) {
    return ( fA*atan((energy-fB)/fC) + fD*atan((energy-fE)/fF) + fG) * sin(fO*dist)
         + ( fH*atan((energy-fI)/fJ) + fK*atan((energy-fL)/fM) + fN) * sin(fO*dist)*sin(fO*dist)*sin(fO*dist);
}

/////////////////////////////////////////////////////////////////////////////

void ReconstructionCombined::CorrectionPositionSinSin3Atan::AddCorrection(const double& minEnergy, CorrectionPositionBin* newBin) {
    std::cerr << "ReconstructionCombined::CorrectionPositionSinSin3Atan::AddCorrection: Position correction of this kind does not use individual energy bins!" << std::endl;
    exit(1);
    // TODO: should be fatal exception, but is caught somewhere
    throw Exception("ReconstructionCombined::CorrectionPositionSinSin3Atan::AddCorrection: Position correction of this kind does not use individual energy bins!");
}

/////////////////////////////////////////////////////////////////////////////

void ReconstructionCombined::CorrectionPositionSinSin3Atan::Print(std::ostream& out, const std::string& title) {
    out << title;
    out << "SinSin3Atan: a=" << fA
                   << ", b=" << fB
                   << ", c=" << fC
                   << ", d=" << fD
                   << ", e=" << fE
                   << ", f=" << fF
                   << ", g=" << fG
                   << ", h=" << fH
                   << ", i=" << fI
                   << ", j=" << fJ
                   << ", k=" << fK
                   << ", l=" << fL
                   << ", m=" << fM
                   << ", n=" << fN
                   << ", o=" << fO << std::endl;
}

/////////////////////////////////////////////////////////////////////////////

void ReconstructionCombined::CorrectionPositionSinSin3Atan::SetSmoothLinear() {
    std::cerr << "ReconstructionCombined::CorrectionPositionSinSin3Atan::SetSmoothLinear: SinSin3Atan does not need this!" << std::endl;
    exit(1);
    // TODO: should be fatal exception, but is caught somewhere
    throw Exception("ReconstructionCombined::CorrectionPositionSinSin3Atan::SetSmoothLinear: SinSin3Atan does not need this!");
}

/////////////////////////////////////////////////////////////////////////////

ReconstructionCombined::CorrectionPositionBinPol3::CorrectionPositionBinPol3(const double& a, const double& b, const double& c) : fA(a), fB(b), fC(c) {
}

/////////////////////////////////////////////////////////////////////////////

double ReconstructionCombined::CorrectionPositionBinPol3::GetCorrection(const double& dist) {
    return fA*dist + fB*dist*dist + fC*dist*dist*dist;
}

/////////////////////////////////////////////////////////////////////////////

void ReconstructionCombined::CorrectionPositionBinPol3::Print(std::ostream& out) {
    out << "Pol3: a=" << fA << ", b=" << fB << ", c=" << fC;
}

/////////////////////////////////////////////////////////////////////////////

ReconstructionCombined::CorrectionPositionBinSin::CorrectionPositionBinSin(const double& a, const double& b) : fA(a), fB(b) {
}

/////////////////////////////////////////////////////////////////////////////

double ReconstructionCombined::CorrectionPositionBinSin::GetCorrection(const double& dist) {
    return fA*sin(fB*dist);
}

/////////////////////////////////////////////////////////////////////////////

void ReconstructionCombined::CorrectionPositionBinSin::Print(std::ostream& out) {
    out << "Sin: a=" << fA << ", b=" << fB;
}

/////////////////////////////////////////////////////////////////////////////

ReconstructionCombined::CorrectionPosDepEnergy::CorrectionPosDepEnergy() : fSmoothLinear(false) {
}

/////////////////////////////////////////////////////////////////////////////

ReconstructionCombined::CorrectionPosDepEnergy::~CorrectionPosDepEnergy() {
    for (std::map<double, CorrectionPosDepEnergyBin*>::const_iterator it=fBins.begin(); it!=fBins.end(); it++)
        delete it->second;
}

/////////////////////////////////////////////////////////////////////////////

double ReconstructionCombined::CorrectionPosDepEnergy::GetCorrection(const double& energy, const double& distX, const double& distY) {
    if (fBins.size()==0)
        return energy;

    if (fSmoothLinear) {
        std::map<double, CorrectionPosDepEnergyBin*>::const_iterator it;
        if ( (it=fBins.lower_bound(energy))!=fBins.begin() ) {
            std::map<double, CorrectionPosDepEnergyBin*>::const_iterator itLower=it; itLower--;

            if (it==fBins.end())
                return itLower->second->GetCorrection(energy, distX, distY);

            const double x1 = itLower->first;
            const double y1 = itLower->second->GetCorrection(energy, distX, distY);
            const double x2 = it     ->first;
            const double y2 = it     ->second->GetCorrection(energy, distX, distY);
            return y1 + (energy-x1) * (y2-y1)/(x2-x1);
        } else
            return it->second->GetCorrection(energy, distX, distY);
    } else {
        std::map<double, CorrectionPosDepEnergyBin*>::const_iterator it;
        if ( (it=fBins.lower_bound(energy))!=fBins.begin() ) {
            it--;
            return it->second->GetCorrection(energy, distX, distY);
        }
    }

    return energy;
}

/////////////////////////////////////////////////////////////////////////////

void ReconstructionCombined::CorrectionPosDepEnergy::AddCorrection(const double& minEnergy, CorrectionPosDepEnergyBin* newBin) {
    fBins[minEnergy] = newBin;
}

/////////////////////////////////////////////////////////////////////////////

void ReconstructionCombined::CorrectionPosDepEnergy::Print(std::ostream& out, const std::string& title) {
    out << title;
    bool first(true);
    if (fSmoothLinear) {
        out << "with linear smoothing" << std::endl;
        first = false;
    }
    double lastEnergy(0);
    CorrectionPosDepEnergyBin* lastCorr(0);
    for (std::map<double, CorrectionPosDepEnergyBin*>::const_iterator it=fBins.begin(); it!=fBins.end(); it++) {
        if (lastCorr!=0) {
            if (!first)
                for (size_t i=0; i<title.length(); i++)
                    out << " ";
            first = false;

            if (fSmoothLinear)
                out << "at E=" << lastEnergy << ": ";
            else
                out << lastEnergy << "<=E<" << it->first << ": ";
            lastCorr->Print(out);
            out << std::endl;
        }

        lastEnergy  = it->first;
        lastCorr    = it->second;
    }
    if (!first)
        for (size_t i=0; i<title.length(); i++)
            out << " ";
    if (fSmoothLinear)
        out << "at E=" << lastEnergy << ": ";
    else
        out << "E>=" << lastEnergy << ": ";
    lastCorr->Print(out);
    out << std::endl;
}

/////////////////////////////////////////////////////////////////////////////

void ReconstructionCombined::CorrectionPosDepEnergy::SetSmoothLinear() {
    fSmoothLinear = true;
}

/////////////////////////////////////////////////////////////////////////////

ReconstructionCombined::CorrectionPosDepEnergyBinQuarter::CorrectionPosDepEnergyBinQuarter(const int& size, const double& stepX, const double& stepY) : fSize(size), fStepX(stepX), fStepY(stepY) {
    fCalib.resize(size*size, 0.);
}

/////////////////////////////////////////////////////////////////////////////

double ReconstructionCombined::CorrectionPosDepEnergyBinQuarter::GetCorrection(const double& energy, const double& distX, const double& distY) {
    const double dX(fabs(distX));
    const double dY(fabs(distY));

    if (dX>(fStepX/2.) || dY>(fStepY/2.))
        return energy;

    const double tX(dX/(fStepX/2.) * fSize);
    const double tY(dY/(fStepY/2.) * fSize);

    int cX((int)tX); if (cX<0) cX=0; if (cX>=fSize) cX=fSize-1;
    int cY((int)tY); if (cY<0) cY=0; if (cY>=fSize) cY=fSize-1;

    return energy * fCalib[cX*fSize + cY];
}

/////////////////////////////////////////////////////////////////////////////

void ReconstructionCombined::CorrectionPosDepEnergyBinQuarter::Read(std::istream& in) {
    for (int i=0; i<fSize; i++)
        for (int j=0; j<fSize; j++) {
            double value;
            in >> value;
            fCalib[i*fSize + j] = value;
        }
}

/////////////////////////////////////////////////////////////////////////////

void ReconstructionCombined::CorrectionPosDepEnergyBinQuarter::Print(std::ostream& out) {
    out << "corrections initialized with " << fSize << "x" << fSize << " values, step sizes of blocks is (" << fStepX << ", " << fStepY << ")";
}

/////////////////////////////////////////////////////////////////////////////

ReconstructionCombined::CorrectionPosDepEnergyFunctionQuarter::CorrectionPosDepEnergyFunctionQuarter() {
}

/////////////////////////////////////////////////////////////////////////////

double ReconstructionCombined::CorrectionPosDepEnergyFunctionQuarter::GetCorrection(const double& energy, const double& distX, const double& distY) {
    const double dX(std::fabs(distX));
    const double dY(std::fabs(distY));
    double retval= fCalibFunction->Eval(dX/10., dY/10., energy);
    if(retval<0)
      retval=1e-3;
    return retval;
}

/////////////////////////////////////////////////////////////////////////////

void ReconstructionCombined::CorrectionPosDepEnergyFunctionQuarter::Read(std::istream& in) {
  std::string value;
  in >> value;
  fCalibFunction = new TF3("fCalibFunctionFunctionQuarter", value.c_str());
}

/////////////////////////////////////////////////////////////////////////////

void ReconstructionCombined::CorrectionPosDepEnergyFunctionQuarter::Print(std::ostream& out) {
    out << "corrections initialized with " << fCalibFunction-> GetExpFormula();
}
/////////////////////////////////////////////////////////////////////////////

} // namespace Reco

#endif // ROOT version >= 5.22

