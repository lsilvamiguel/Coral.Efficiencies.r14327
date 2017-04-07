/*
   $Source: /afs/cern.ch/user/g/gobbo/w0/coral/src/calorim/Reco/src/ReconstructionCombined.cc,v $
   $Date: 2011/02/04 18:15:57 $
   $Revision: 1.55 $
*/

#include <RVersion.h>
#if ROOT_VERSION_CODE >= ROOT_VERSION(5,22,00)

// --- Standard C/C++ library ---
#include <cassert>
#include <cctype>
#include <iostream>
#include <sstream>
#include <set>

#include <Math/Factory.h>
#include <Math/Minimizer.h>

// --- Internal files ---
#include "ReconstructionCombined.h"

#include "Calorimeter.h"
#include "CalorimeterParticle.h"
#include "CellDataRaw.h"
#include "Cluster.h"
#include "Exception.h"
#include "Shower.h"
#include "ShowerProfile.h"

namespace Reco {

/////////////////////////////////////////////////////////////////////////////

ReconstructionCombined::ReconstructionCombined(const Calorimeter* c) : Reconstruction(c), fPosCorrX(0), fPosCorrY(0) {
    fHist.resize(NUMBER_OF_HISTOGRAMS, 0);

    std::cout << "Initializing ReconstructionCombined for " << GetCalorimeter()->GetName() << " with parameters:" << std::endl;
    std::cout << "    cell energy threshold (GeV)                   : " << GetCalorimeter()->GetOptions().cell_energy_threshold << std::endl;
    std::cout << "    cell energy threshold (ADC channels)          : " << GetCalorimeter()->GetOptions().cell_adc_threshold << std::endl;
    std::cout << "    particle energy threshold (GeV)               : " << GetCalorimeter()->GetOptions().particle_energy_threshold << std::endl;
    std::cout << "    minimal distance of two showers (block sizes) : " << GetCalorimeter()->GetOptions().min_shower_distance << std::endl;
    std::cout << "    longitudinal leakage correction               : " << (GetCalorimeter()->GetOptions().correct_longitudinal_leakage ? "activated" : "deactivated") << std::endl;
    std::cout << "    use calibration file for further refinement   : " << (GetCalorimeter()->GetOptions().use_combined_calibration ? "activated" : "deactivated") << std::endl;
    std::cout << "    histogramming                                 : " << (GetCalorimeter()->GetOptions().fill_histos ? "activated" : "deactivated") << std::endl;
    std::cout << "    histogramming level                           : ";
    if      (GetCalorimeter()->GetOptions().histos_level==Calorimeter::Options::HISTLVL_NORMAL)  std::cout << "normal" << std::endl;
    else if (GetCalorimeter()->GetOptions().histos_level==Calorimeter::Options::HISTLVL_VERBOSE) std::cout << "verbose" << std::endl;
    else if (GetCalorimeter()->GetOptions().histos_level==Calorimeter::Options::HISTLVL_DEBUG)   std::cout << "debug" << std::endl;
    else                                                                                         std::cout << "unknown (indicating a bug!)" << std::endl;
}

/////////////////////////////////////////////////////////////////////////////

ReconstructionCombined::~ReconstructionCombined() {
    if (fPosCorrX != 0)
        delete fPosCorrX;
    if (fPosCorrY != 0)
        delete fPosCorrY;
    for (std::map<CellType::HwType, CorrectionPosDepEnergy*>::const_iterator it=fPosDepEnergyCorr.begin(); it!=fPosDepEnergyCorr.end(); it++)
        delete it->second;
}

/////////////////////////////////////////////////////////////////////////////

void ReconstructionCombined::ReadCalibrations() {
    std::cout << "Reading calibrations in ReconstructionCombined for " << GetCalorimeter()->GetName() << ":" << std::endl;

    // ========================================
    // parameterization of energy error
    // ========================================

    std::cout << "    parameterization of energy error              : ";

    if (GetCalorimeter()->GetOptions().param_energy_error.size() > 0)
        fParamEnergyErr = GetCalorimeter()->GetOptions().param_energy_error;
    else {
        fParamEnergyErr.push_back(0.15);  fParamEnergyErr.push_back(0.015); fParamEnergyErr.push_back(0.05);
        std::cout << "[default values] ";
    }

    for (std::vector<double>::const_iterator it=fParamEnergyErr.begin(); it!=fParamEnergyErr.end(); it++)
        std::cout << *it << " ";

    std::cout << std::endl;
    assert(fParamEnergyErr.size() == 3);

    // ========================================
    // parameterization of time error
    // ========================================

    std::cout << "    parameterization of time error                : ";

    if (GetCalorimeter()->GetOptions().param_time_error.size() > 0)
        fParamTimeErr = GetCalorimeter()->GetOptions().param_time_error;
    else {
        fParamTimeErr.push_back(0.9055);  fParamTimeErr.push_back(1.3077);  fParamTimeErr.push_back(0.5831);
        std::cout << "[default values] ";
    }

    for (std::vector<double>::const_iterator it=fParamTimeErr.begin(); it!=fParamTimeErr.end(); it++)
        std::cout << *it << " ";

    std::cout << std::endl;
    assert(fParamEnergyErr.size() == 3);

    // ========================================
    // number of allowed showers per cluster
    // ========================================

    std::cout << "    number of allowed showers for cluster size    : ";

    if (GetCalorimeter()->GetOptions().combined_allowed_showers.size() > 0) {
        fAllowedNrShowers.clear();

        unsigned int curStep = GetCalorimeter()->GetOptions().combined_allowed_showers[0];
        for (unsigned int i=0; i<curStep; i++)
            fAllowedNrShowers.push_back(0);

        unsigned int size(1);
        for (size_t j=1; j<GetCalorimeter()->GetOptions().combined_allowed_showers.size(); j++) {
            unsigned int nextStep = GetCalorimeter()->GetOptions().combined_allowed_showers[j];
            assert(nextStep > curStep);
            for (unsigned int i=curStep; i<nextStep; i++)
                fAllowedNrShowers.push_back(size);
            curStep = nextStep;
            size++;
        }
        fAllowedNrShowers.push_back(size);
    } else {
        // fill the array of allowed number of showers only for cluster sizes
        // up to 9, larger clusters will take the last value
        fAllowedNrShowers.clear();
        for (unsigned int i=0; i<10; i++) {
            unsigned int showers(1);
            while ( ( (showers+1)==2 && i>1 ) ||
                    ( (showers+1)==3 && i>2 ) ||
                    ( (showers+1)==4 && i>3 ) ||
                    ( (showers+1)==5 && i>5 ) ||
                    ( (showers+1)==6 && i>7 ) )
                showers++;

            assert(fAllowedNrShowers.size()==i);
            fAllowedNrShowers.push_back(showers);
        }

        std::cout << "[default behaviour] ";
    }

    for (size_t i=0; i<fAllowedNrShowers.size()-1; i++)
        std::cout << i << "->" << fAllowedNrShowers[i] << " ";
    std::cout << fAllowedNrShowers.size()-1 << ",...->" << fAllowedNrShowers.back() << std::endl;

    // ========================================
    // print the shower profiles in use
    // ========================================

    std::cout << "    shower profiles                               : " << std::endl;
    std::set<CellType::HwType> seen;
    for (std::list<CellType>::const_iterator it=GetCalorimeter()->GetCellsTypeList().begin(); it!=GetCalorimeter()->GetCellsTypeList().end(); it++) {
        const CellType& cellType           = *it;

        if (seen.count(cellType.GetHwType())==0) {
            seen.insert(cellType.GetHwType());

            std::string title;
            for (size_t i=0; i<44-CellType::GetHwTypeString(cellType.GetHwType()).length(); i++)
                title += " ";
            title += "* for " + CellType::GetHwTypeString(cellType.GetHwType()) + ": ";

            const ShowerProfile* showerProfile = cellType.GetShowerProfile();
            showerProfile->Print(std::cout, title);
        }
    }

    // ========================================
    // read a few remaining things from the
    // calibration file
    // ========================================

    // check if the calibrations are to be used
    if (!GetCalorimeter()->GetOptions().use_combined_calibration)
        return;

    std::istringstream calib(GetCalorimeter()->GetOptions().combined_calibration);
    std::string line;
    while (getline(calib, line)) {
        // remove the comments from a string
        // comments start with //
        if (line.find("//") < line.size())
            line.erase(line.find("//"));

        // reduce the nr of white spaces, replace tabs with spaces, replace
        // two spaces by one space, remove leading and trailing whitespaces
        while (line.find("\t") < line.size())
            line.replace(line.find("\t"), 1, " ");
        while (line.find("  ") < line.size())
            line.replace(line.find("  "), 2, " ");
        while (line.find(" ") == 0 && !line.empty())
            line.erase(0, 1);
        while (line.rfind(" ") == line.size()-1 && !line.empty())
            line.erase(line.size()-1);

        if (line.empty() || line==" ")
            continue;

        std::string token;
        std::istringstream curLine(line);
        curLine >> token;
        std::transform(token.begin(), token.end(), token.begin(), ::toupper);
        if (token == "TAG") {
            std::cout << "    tag of calibration file                       : ";
            std::string temp;
            while (curLine >> temp) std::cout << temp << " ";
            std::cout << std::endl;
        } else if (token == "CORRECTIONPOSITIONX") {
            double minEnergy;
            if ( !(curLine >> minEnergy) ) {
                curLine.clear();

                std::string flags;
                if ( !(curLine >> flags) ) {
                    std::cerr << "ReconstructionCombined::ReadCalibrations: Bad format of input line while reading from calibration file:" << std::endl << line << std::endl;
                    exit(1);
                    // TODO: should be fatal exception, but is caught somewhere
                    throw Exception("ReconstructionCombined::ReadCalibrations: Bad format of input line while reading from calibration file:\n%s", line.c_str());
                }
                std::transform(flags.begin(), flags.end(), flags.begin(), ::toupper);
                if (flags == "SMOOTHLINEAR") {
                    if (fPosCorrX == 0)
                        fPosCorrX = new CorrectionPosition();
                    fPosCorrX->SetSmoothLinear();
                } else {
                    std::cerr << "ReconstructionCombined::ReadCalibrations: Unknown flag '" << flags << "' while reading from calibration file:" << std::endl << line << std::endl;
                    exit(1);
                    // TODO: should be fatal exception, but is caught somewhere
                    throw Exception("ReconstructionCombined::ReadCalibrations: Unkown flag '%s' while reading from calibration file:\n%s", flags.c_str(), line.c_str());
                }
            } else {
                double a, b;
                if ( !(curLine >> a >> b) ) {
                    std::cerr << "ReconstructionCombined::ReadCalibrations: Bad format of input line while reading from calibration file:" << std::endl << line << std::endl;
                    exit(1);
                    // TODO: should be fatal exception, but is caught somewhere
                    throw Exception("ReconstructionCombined::ReadCalibrations: Bad format of input line while reading from calibration file:\n%s", line.c_str());
                }

                double c;
                if ( !(curLine >> c) ) {
                    if (fPosCorrX == 0)
                        fPosCorrX = new CorrectionPosition();
                    fPosCorrX->AddCorrection(minEnergy, new CorrectionPositionBinSin(a, b));
                } else {
                    double d, e, f, g, h, i, j, k;
                    if ( !(curLine >> d >> e >> f >> g >> h >> i >> j >> k) ) {
                        if (fPosCorrX == 0)
                            fPosCorrX = new CorrectionPosition();
                        fPosCorrX->AddCorrection(minEnergy, new CorrectionPositionBinPol3(a, b, c));
                    } else {
                        double l, m, n, o;
                        if ( !(curLine >> l >> m >> n >> o) ) {
                            if (fPosCorrX == 0)
                                fPosCorrX = new CorrectionPositionSinAtan(a, b, c, d, e, f, g, h, i, j, k);
                            else {
                                std::cerr << "ReconstructionCombined::ReadCalibrations: Position correction object for X already exists, SinAtan must be the only position correction in X!" << std::endl;
                                exit(1);
                                // TODO: should be fatal exception, but is caught somewhere
                                throw Exception("ReconstructionCombined::ReadCalibrations: Position correction object for X already exists, SinAtan must be the only position correction in X!");
                            }
                        } else {
                            if (fPosCorrX == 0)
                                fPosCorrX = new CorrectionPositionSinSin3Atan(a, b, c, d, e, f, g, h, i, j, k, l, m, n, o);
                            else {
                                std::cerr << "ReconstructionCombined::ReadCalibrations: Position correction object for X already exists, SinSin3Atan must be the only position correction in X!" << std::endl;
                                exit(1);
                                // TODO: should be fatal exception, but is caught somewhere
                                throw Exception("ReconstructionCombined::ReadCalibrations: Position correction object for X already exists, SinSin3Atan must be the only position correction in X!");
                            }
                        }
                    }
                }
            }
        } else if (token == "CORRECTIONPOSITIONY") {
            double minEnergy;
            if ( !(curLine >> minEnergy) ) {
                curLine.clear();

                std::string flags;
                if ( !(curLine >> flags) ) {
                    std::cerr << "ReconstructionCombined::ReadCalibrations: Bad format of input line while reading from calibration file:" << std::endl << line << std::endl;
                    exit(1);
                    // TODO: should be fatal exception, but is caught somewhere
                    throw Exception("ReconstructionCombined::ReadCalibrations: Bad format of input line while reading from calibration file:\n%s", line.c_str());
                }
                std::transform(flags.begin(), flags.end(), flags.begin(), ::toupper);
                if (flags == "SMOOTHLINEAR") {
                    if (fPosCorrY == 0)
                        fPosCorrY = new CorrectionPosition();
                    fPosCorrY->SetSmoothLinear();
                } else {
                    std::cerr << "ReconstructionCombined::ReadCalibrations: Unknown flag '" << flags << "' while reading from calibration file:" << std::endl << line << std::endl;
                    exit(1);
                    // TODO: should be fatal exception, but is caught somewhere
                    throw Exception("ReconstructionCombined::ReadCalibrations: Unkown flag '%s' while reading from calibration file:\n%s", flags.c_str(), line.c_str());
                }
            } else {
                double a, b;
                if ( !(curLine >> a >> b) ) {
                    std::cerr << "ReconstructionCombined::ReadCalibrations: Bad format of input line while reading from calibration file:" << std::endl << line << std::endl;
                    exit(1);
                    // TODO: should be fatal exception, but is caught somewhere
                    throw Exception("ReconstructionCombined::ReadCalibrations: Bad format of input line while reading from calibration file:\n%s", line.c_str());
                }

                double c;
                if ( !(curLine >> c) ) {
                    if (fPosCorrY == 0)
                        fPosCorrY = new CorrectionPosition();
                    fPosCorrY->AddCorrection(minEnergy, new CorrectionPositionBinSin(a, b));
                } else {
                    double d, e, f, g, h, i, j, k;
                    if ( !(curLine >> d >> e >> f >> g >> h >> i >> j >> k) ) {
                        if (fPosCorrY == 0)
                            fPosCorrY = new CorrectionPosition();
                        fPosCorrY->AddCorrection(minEnergy, new CorrectionPositionBinPol3(a, b, c));
                    } else {
                        double l, m, n, o;
                        if ( !(curLine >> l >> m >> n >> o) ) {
                            if (fPosCorrY == 0)
                                fPosCorrY = new CorrectionPositionSinAtan(a, b, c, d, e, f, g, h, i, j, k);
                            else {
                                std::cerr << "ReconstructionCombined::ReadCalibrations: Position correction object for Y already exists, SinAtan must be the only position correction in Y!" << std::endl;
                                exit(1);
                                // TODO: should be fatal exception, but is caught somewhere
                                throw Exception("ReconstructionCombined::ReadCalibrations: Position correction object for Y already exists, SinAtan must be the only position correction in Y!");
                            }
                        } else {
                            if (fPosCorrY == 0)
                                fPosCorrY = new CorrectionPositionSinSin3Atan(a, b, c, d, e, f, g, h, i, j, k, l, m, n, o);
                            else {
                                std::cerr << "ReconstructionCombined::ReadCalibrations: Position correction object for Y already exists, SinSin3Atan must be the only position correction in Y!" << std::endl;
                                exit(1);
                                // TODO: should be fatal exception, but is caught somewhere
                                throw Exception("ReconstructionCombined::ReadCalibrations: Position correction object for Y already exists, SinSin3Atan must be the only position correction in Y!");
                            }
                        }
                    }
                }
            }
        } else if (token == "CORRECTIONENERGYPOSDEP") {
            std::string cellType;
            if ( !(curLine >> cellType) ) {
                std::cerr << "ReconstructionCombined::ReadCalibrations: Bad format of input line while reading from calibration file:" << std::endl << line << std::endl;
                exit(1);
                // TODO: should be fatal exception, but is caught somewhere
                throw Exception("ReconstructionCombined::ReadCalibrations: Bad format of input line while reading from calibration file:\n%s", line.c_str());
            }

            CellType::HwType type = CellType::hwt_unknown;
            try {
                type = CellType::GetHwTypeFromString(cellType);
            } catch (Exception& e) {
                std::cerr << "ReconstructionCombined::ReadCalibrations: Exception catched while trying to get HwType from type name:" << std::endl << e.what() << std::endl;
                exit(1);
                // TODO: should be fatal exception, but is caught somewhere
                throw Exception("ReconstructionCombined::ReadCalibrations: Bad format of input line while reading from calibration file:\n%s", e.what());
            }

            double minEnergy;
            std::string shape;
            if ( !(curLine >> minEnergy >> shape) ) {
                curLine.clear();

                std::string flags;
                if ( !(curLine >> flags) ) {
                    std::cerr << "ReconstructionCombined::ReadCalibrations: Bad format of input line while reading from calibration file:" << std::endl << line << std::endl;
                    exit(1);
                    // TODO: should be fatal exception, but is caught somewhere
                    throw Exception("ReconstructionCombined::ReadCalibrations: Bad format of input line while reading from calibration file:\n%s", line.c_str());
                }
                std::transform(flags.begin(), flags.end(), flags.begin(), ::toupper);
                if (flags == "SMOOTHLINEAR") {
                    if (fPosDepEnergyCorr.count(type)==0)
                        fPosDepEnergyCorr[type] = new CorrectionPosDepEnergy();
                    fPosDepEnergyCorr[type]->SetSmoothLinear();
                } else {
                    std::cerr << "ReconstructionCombined::ReadCalibrations: Unknown flag '" << flags << "' while reading from calibration file:" << std::endl << line << std::endl;
                    exit(1);
                    // TODO: should be fatal exception, but is caught somewhere
                    throw Exception("ReconstructionCombined::ReadCalibrations: Unkown flag '%s' while reading from calibration file:\n%s", flags.c_str(), line.c_str());
                }
            } else {
                transform(shape.begin(), shape.end(), shape.begin(), ::toupper);
                if (shape == "QUARTER") {
                    int size(0);
                    if ( !(curLine >> size) ) {
                        std::cerr << "ReconstructionCombined::ReadCalibrations: Bad format of input line while reading from calibration file:" << std::endl << line << std::endl;
                        exit(1);
                        // TODO: should be fatal exception, but is caught somewhere
                        throw Exception("ReconstructionCombined::ReadCalibrations: Bad format of input line while reading from calibration file:\n%s", line.c_str());
                    }
                    assert(size>0);

                    // search for cell type to get step sizes
                    const std::list<CellType>& allTypes(GetCalorimeter()->GetCellsTypeList());
                    std::list<CellType>::const_iterator thisType(allTypes.end());
                    for (std::list<CellType>::const_iterator it=allTypes.begin(); it!=allTypes.end(); it++)
                        if (it->GetHwType()==type) {
                            thisType=it;
                            break;
                        }
                    if (thisType==allTypes.end()) {
                        std::cerr << "ReconstructionCombined::ReadCalibrations: Could not find the cell type " << cellType << " in the used blocks of this calorimeter, bad it is in the calibration file." << std::endl;
                        exit(1);
                        // TODO: should be fatal exception, but is caught somewhere
                        throw Exception("ReconstructionCombined::ReadCalibrations: Could not find the cell type %s in the used blocks of this calorimeter, bad it is in the calibration file.", cellType.c_str());
                    }

                    CorrectionPosDepEnergyBinQuarter* corr = new CorrectionPosDepEnergyBinQuarter(size, thisType->GetStepX(), thisType->GetStepY());
                    corr->Read(calib);

                    if (fPosDepEnergyCorr.count(type)==0)
                        fPosDepEnergyCorr[type] = new CorrectionPosDepEnergy();

                    fPosDepEnergyCorr[type]->AddCorrection(minEnergy, corr);
                }
                else if (shape == "QUARTERFUNCTION") {
                  const std::list<CellType>& allTypes(GetCalorimeter()->GetCellsTypeList());
                  std::list<CellType>::const_iterator thisType(allTypes.end());
                  for (std::list<CellType>::const_iterator it=allTypes.begin(); it!=allTypes.end(); it++)
                    if (it->GetHwType()==type) {
                      thisType=it;
                      break;
                    }
                  if (thisType==allTypes.end()) {
                    std::cerr << "ReconstructionCombined::ReadCalibrations: Could not find the cell type " << cellType << " in the used blocks of this calorimeter, bad it is in the calibration file." << std::endl;
                    exit(1);
                    // TODO: should be fatal exception, but is caught somewhere
                    throw Exception("ReconstructionCombined::ReadCalibrations: Could not find the cell type %s in the used blocks of this calorimeter, bad it is in the calibration file.", cellType.c_str());
                  }

                  CorrectionPosDepEnergyFunctionQuarter* corr = new CorrectionPosDepEnergyFunctionQuarter();
                  corr->Read(calib);

                  if (fPosDepEnergyCorr.count(type)==0)
                    fPosDepEnergyCorr[type] = new CorrectionPosDepEnergy();

                  fPosDepEnergyCorr[type]->AddCorrection(minEnergy, corr);
                }

            }
        } else {
            std::cerr << "Unknown token '" << token << "'." << std::endl;
        }
    }

    if (fPosCorrX != 0)
        fPosCorrX->Print(std::cout, "    position corrections in X                     : ");
    else
        std::cout << "    position corrections in X                     : [none]" << std::endl;
    if (fPosCorrY != 0)
        fPosCorrY->Print(std::cout, "    position corrections in Y                     : ");
    else
        std::cout << "    position corrections in Y                     : [none]" << std::endl;

    std::cout << "    position dependent energy correction          : ";
    if (fPosDepEnergyCorr.size()==0)
        std::cout << "[none]" << std::endl;
    else {
        std::cout << std::endl;
        for (std::map<CellType::HwType, CorrectionPosDepEnergy*>::const_iterator it=fPosDepEnergyCorr.begin(); it!=fPosDepEnergyCorr.end(); it++) {
            std::string title;
            for (size_t i=0; i<44-CellType::GetHwTypeString(it->first).length(); i++)
                title += " ";
            title += "* for " + CellType::GetHwTypeString(it->first) + ": ";
            it->second->Print(std::cout, title);
        }
    }
}

/////////////////////////////////////////////////////////////////////////////

const std::vector<CalorimeterParticle>& ReconstructionCombined::DoReconstruction(const std::vector<CellDataRaw>& signals) {
    reconstructedParticles.clear();

    std::vector<Cluster> clusters = DoClustering(signals);

    for (std::vector<Cluster>::const_iterator it=clusters.begin(); it!=clusters.end(); it++) {
        // a map for the cell number to the index in the signal vector
        std::map<size_t, size_t> mapCellSignal;
        for (size_t i=0; i<it->GetCluster().size(); i++) {
            if (mapCellSignal.count(it->GetCluster()[i].GetCellIdx())>0)
                throw Exception("ReconstructionCombined::DoReconstruction: Multiple raw data information for the same cell.");

            mapCellSignal.insert(std::pair<size_t, size_t>(it->GetCluster()[i].GetCellIdx(), i));
        }

        std::map<size_t, size_t> mapCellsIndices;
        std::vector<FitInfo> cells;
        std::set<size_t> deadCells;
        // loop over the cells in the cluster to find if we have to add dead
        // cells or boundary cells (zero energy, large error)
        for (std::vector<CellDataRaw>::const_iterator it2=it->GetCluster().begin();
             it2!=it->GetCluster().end(); it2++) {
            // TODO: filter out dead cells

            const std::vector<size_t>& neighbors = GetCalorimeter()->GetCells()[it2->GetCellIdx()].GetNeighbors();
            for (std::vector<size_t>::const_iterator it3=neighbors.begin()+1; it3!=neighbors.end(); it3++)
                if (GetCalorimeter()->CellIsBad(*it3, Calorimeter::OLD))
                    deadCells.insert(*it3);

            mapCellsIndices.insert(std::pair<size_t, size_t>(it2->GetCellIdx(), cells.size()));
            cells.push_back(FitInfo(GetCalorimeter(), *it2));
            cells.back().SetEnergyErr(CalcEnergyError(it2->GetEnergy()));
            cells.back().SetTimeErr(CalcTimeError(it2->GetEnergy()));
        }

        // TODO: make sure no dead cell is amongst the cells in cells

        for (std::vector<FitInfo>::const_iterator it2=cells.begin(); it2!=cells.end(); it2++)
            if (deadCells.count(it2->GetData().GetCellIdx())>0)
                deadCells.erase(it2->GetData().GetCellIdx());

        for (std::set<size_t>::const_iterator it2=deadCells.begin(); it2!=deadCells.end(); it2++) {
            // find the energy in the cells around the dead cell to initialise
            // the energy error
            // also initialize the time error taking into account the timing
            // spread in neighboring cells
            double eAround(0.);
            double tAroundMax(- std::numeric_limits<double>::max());
            double tAroundMin(  std::numeric_limits<double>::max());
            const std::vector<size_t>& neighbors = GetCalorimeter()->GetCells()[*it2].GetNeighbors();
            for (std::vector<size_t>::const_iterator it3=neighbors.begin()+1; it3!=neighbors.end(); it3++) {
                std::map<size_t, size_t>::const_iterator curMap;
                if ( (curMap=mapCellSignal.find(*it3))!=mapCellSignal.end() ) {
                    if (it->GetCluster()[curMap->second].GetEnergy()>eAround)
                        eAround = it->GetCluster()[curMap->second].GetEnergy();
                    tAroundMax = std::max(tAroundMax, it->GetCluster()[curMap->second].GetTime());
                    tAroundMin = std::min(tAroundMin, it->GetCluster()[curMap->second].GetTime());
                }
            }

            double tAroundDiff(tAroundMax - tAroundMin);
            if (tAroundDiff < CalcTimeError(eAround))
                tAroundDiff = CalcTimeError(eAround);

            mapCellsIndices.insert(std::pair<size_t, size_t>(*it2, cells.size()));
            cells.push_back(FitInfo(GetCalorimeter(), CellDataRaw(*it2, 0., 0.)));
            cells.back().SetEnergyErr(eAround);
            cells.back().SetTimeErr(3.*tAroundDiff);
            cells.back().SetArtificial(true);
        }

        // fit parameters
        std::vector<FitShower> shower;
        double logLh(-1.);
        double ndf(-1.);
        // previous fit
        std::vector<FitShower> oldShower;
        double oldLogLh(-1.);
        double oldNdf(-1.);

        // fitting status
        bool doFit(true);

        // in case of only one cell do not fit
        const size_t length=cells.size();
        // TODO: length-deadcells.size() ?
        if (length==1) {
            doFit = false;

            // energy is cell energy corrected for "leakage" into neighboring cells
            double e = cells.begin()->GetData().GetEnergy() / cells.begin()->GetShowerProfile()->GetEnergyInCell(0., 0.);
            double x = GetCalorimeter()->GetCells()[cells.begin()->GetData().GetCellIdx()].GetX();
            double y = GetCalorimeter()->GetCells()[cells.begin()->GetData().GetCellIdx()].GetY();
            double t = cells.begin()->GetData().GetTime();

            // initial error of energy
            double eErr = cells.begin()->GetData().GetEnergyErr();

            // initial error in x and y:
            // cell-size / sqrt(12)
            double xErr=GetCalorimeter()->GetCells()[cells.begin()->GetData().GetCellIdx()].GetCellType().GetTrueSizeX() / sqrt(12.);
            double yErr=GetCalorimeter()->GetCells()[cells.begin()->GetData().GetCellIdx()].GetCellType().GetTrueSizeY() / sqrt(12.);

            // initial error of time
            double tErr = cells.begin()->GetData().GetTimeErr();

            if (e >= GetCalorimeter()->GetOptions().particle_energy_threshold) {
                shower.push_back(FitShower(e, x, y, t, eErr, xErr, yErr, tErr));
                shower.back().fMainCell = cells.begin()->GetData().GetCellIdx();
            }
        } else {
            // surround the cluster with two "rows" of cells containing
            // 1. the cut-off energy and some error
            // 2. energy 0 and some error

            for (size_t i=0; i<length; i++) {
                const std::vector<size_t>& neighbors = GetCalorimeter()->GetCells()[cells[i].GetData().GetCellIdx()].GetNeighbors();
                for (std::vector<size_t>::const_iterator it3=neighbors.begin()+1; it3!=neighbors.end(); it3++) {
                    if (mapCellsIndices.count(*it3)==0) {
                        // cut-off energy is either the cell_energy_threshold
                        // applied in the first step of this reconstruction,
                        // the clustering, ...
                        double thrEnergy(GetCalorimeter()->GetOptions().cell_energy_threshold);
                        // ..., or the minimum energy introduced by the ADC
                        double adcEnergy(GetCalorimeter()->SignalToEnergy(*it3, GetCalorimeter()->GetOptions().cell_adc_threshold, GetCalorimeter()->GetEventID().GetTimeInSpill()));
                        // in any case, the larger of the two defines the
                        // energy
                        if (adcEnergy > thrEnergy)
                            thrEnergy = adcEnergy;

                        CellDataRaw newCell(*it3, thrEnergy, 0.);
                        // TODO realistic energy error, probably related to the ADC cut?
                        newCell.SetEnergyErr(CalcEnergyError(thrEnergy));
                        newCell.SetTimeErr(1e12);
                        mapCellsIndices.insert(std::pair<size_t, size_t>(*it3, cells.size()));
                        cells.push_back(FitInfo(GetCalorimeter(), newCell));
                        cells.back().SetArtificial(true);
                        cells.back().SetAround(true);
                    }
                }
            }
        }

        unsigned int maxShowerCount(fAllowedNrShowers.back());
        if (length<fAllowedNrShowers.size())
            maxShowerCount = fAllowedNrShowers[length];

        while (doFit) {
            // another round of fitting, try to add one more shower
            // in case we do not find enough unfitted energy, also quit the fit
            doFit=false;

            // check whether we have to rescale the shower energies:
            // if the energy described by the showers is larger than the cell energy, then
            // rescale the shower energies (just for the decision whether to add a new shower
            // or not, before the next fit, the original values are put back)
            double scale(1.);
            if (!shower.empty()) {
                // loop over main cells of showers
                for (std::vector<FitShower>::iterator it2=shower.begin(); it2!=shower.end(); it2++) {
                    const FitInfo& cell = cells[mapCellsIndices[it2->fMainCell]];

                    double fittedEnergy(cell.GetFittedEnergy(shower));
                    if (fittedEnergy > cell.GetData().GetEnergy() && (fittedEnergy / cell.GetData().GetEnergy()) > scale)
                        scale = fittedEnergy / cell.GetData().GetEnergy();
                }
                for (std::vector<FitShower>::iterator it2=shower.begin(); it2!=shower.end(); it2++)
                    it2->fE /= scale;
            }

            double maxCompare(0.);
            double e(0.);    double x(0.);    double y(0.);    double t(0.);
            double eErr(0.); double xErr(0.); double yErr(0.); double tErr(0.);
            // initialize the fit with the highest unfitted energy
            for (std::vector<FitInfo>::const_iterator it2=cells.begin(); it2!=cells.end(); it2++) {
                // do not create a new start point in one of the cells that were artificially added
                if (it2->GetArtificial())
                    continue;

                // make sure there is a certain distance to the next shower
                double minDist=10000.;
                for (std::vector<FitShower>::const_iterator it3=shower.begin(); it3!=shower.end(); it3++) {
                    double distX=GetCalorimeter()->GetCells()[it2->GetData().GetCellIdx()].GetX() - it3->fX;
                    double distY=GetCalorimeter()->GetCells()[it2->GetData().GetCellIdx()].GetY() - it3->fY;
                    double dist=sqrt(distX*distX + distY*distY);
                    if (dist<minDist)
                        minDist = dist;
                }
                if (minDist<GetCalorimeter()->GetOptions().min_shower_distance)
                    continue;

                // but is this also the cell with the highest unfitted energy
                // in case of no shower, we search for the highest energy directly
                // otherwise we search for the ratio
                double compare((it2->GetData().GetEnergy()-it2->GetFittedEnergy(shower))/it2->GetData().GetEnergy());
                if (shower.empty())
                    compare = it2->GetData().GetEnergy();

                if (compare>maxCompare) {
                    maxCompare=compare;

                    e=0.; x=0.; y=0.; t=0.;
                    const std::vector<size_t>& neighbors = GetCalorimeter()->GetCells()[it2->GetData().GetCellIdx()].GetNeighbors();
                    for (std::vector<size_t>::const_iterator it3=neighbors.begin(); it3!=neighbors.end(); it3++) {
                        std::map<size_t, size_t>::const_iterator idx=mapCellsIndices.find(*it3);
                        if (idx==mapCellsIndices.end())
                            continue;
                        const FitInfo& info = cells[idx->second];
                        if (info.GetData().GetEnergy()-info.GetFittedEnergy(shower) > 0.) {
                            e+=(info.GetData().GetEnergy()-info.GetFittedEnergy(shower));
                            x+=(info.GetData().GetEnergy()-info.GetFittedEnergy(shower)) * GetCalorimeter()->GetCells()[info.GetData().GetCellIdx()].GetX();
                            y+=(info.GetData().GetEnergy()-info.GetFittedEnergy(shower)) * GetCalorimeter()->GetCells()[info.GetData().GetCellIdx()].GetY();
                            t+=(info.GetData().GetEnergy()-info.GetFittedEnergy(shower)) * info.GetData().GetTime();
                        }
                    }

                    x /= e;
                    y /= e;
                    t /= e;

                    // initial error in energy:
                    eErr = it2->GetData().GetEnergyErr();
                    if (eErr<e/10.) eErr = e/10.;

                    // initial error in x and y:
                    // cell-size / sqrt(12)
                    xErr=GetCalorimeter()->GetCells()[it2->GetData().GetCellIdx()].GetCellType().GetTrueSizeX() / sqrt(12.);
                    yErr=GetCalorimeter()->GetCells()[it2->GetData().GetCellIdx()].GetCellType().GetTrueSizeY() / sqrt(12.);

                    // initialize error in time
                    // time resolution based on the energy deposite in this cell
                    tErr=it2->GetData().GetTimeErr();

                    doFit=true;
                }
            }

            if ( (shower.size()+1) > maxShowerCount ) {
                shower = oldShower;
                logLh  = oldLogLh;
                ndf    = oldNdf;
                break;
            }

            if (!doFit) {
                shower = oldShower;
                logLh  = oldLogLh;
                ndf    = oldNdf;
                break;
            }

            shower.push_back(FitShower(e, x, y, t, eErr, xErr, yErr, tErr));

            for (std::vector<FitShower>::iterator it2=shower.begin(); it2!=shower.end(); it2++) {
                // except for the shower just added
                if (it2->fMainCell!=std::numeric_limits<size_t>::max()) {
                    // initial error in energy:
                    if (it->GetCluster()[mapCellSignal[it2->fMainCell]].HasEnergyErr())
                        it2->fEerr = it->GetCluster()[mapCellSignal[it2->fMainCell]].GetEnergyErr();
                    else {
                        const double e=it->GetCluster()[mapCellSignal[it2->fMainCell]].GetEnergy();
                        it2->fEerr = CalcEnergyError(e);
                    }
                    if (it2->fEerr<it2->fE/10.) it2->fEerr = it2->fE/10.;

                    // initial error in x and y:
                    // cell-size / sqrt(12)
                    it2->fXerr=GetCalorimeter()->GetCells()[it2->fMainCell].GetCellType().GetTrueSizeX() / sqrt(12.);
                    it2->fYerr=GetCalorimeter()->GetCells()[it2->fMainCell].GetCellType().GetTrueSizeY() / sqrt(12.);

                    // initialize error in time
                    // time resolution based on the energy deposite in this cell
                    const double e=it->GetCluster()[mapCellSignal[it2->fMainCell]].GetEnergy();
                    it2->fTerr=CalcTimeError(e);
                }
            }

            logLh = DoFitting(cells, shower);
            // TODO replace cells.size() with length (the cells with >0 energy)
            // or the other way round
            ndf   = length - ((double)nrFitParameters)*shower.size();
// DEBUG
//            std::cout << cells.size() << " " << it->GetCluster().size() << " " << minuit->GetNumPars() << " " << ndf << " " << shower.size() << std::endl;

            if (shower.size()==1 && GetCalorimeter()->GetOptions().fill_histos && GetCalorimeter()->GetOptions().histos_level>=Calorimeter::Options::HISTLVL_DEBUG) {
                assert(fHist[h1D_fitDiffEnergyFirst]!=0);
                fHist[h1D_fitDiffEnergyFirst]->Fill(e - shower[0].fE);
                assert(fHist[h2D_fitDiffEnergyFirstVsEnergy]!=0);
                fHist[h2D_fitDiffEnergyFirstVsEnergy]->Fill(e - shower[0].fE, shower[0].fE);
                assert(fHist[h2D_fitDiffEnergyFirstVsEnergyLow]!=0);
                fHist[h2D_fitDiffEnergyFirstVsEnergyLow]->Fill(e - shower[0].fE, shower[0].fE);
                assert(fHist[h1D_fitDiffPosXFirst]!=0);
                fHist[h1D_fitDiffPosXFirst]->Fill(x - shower[0].fX);
                assert(fHist[h1D_fitDiffPosYFirst]!=0);
                fHist[h1D_fitDiffPosYFirst]->Fill(y - shower[0].fY);
                assert(fHist[h1D_fitDiffTimeFirst]!=0);
                fHist[h1D_fitDiffTimeFirst]->Fill(t - shower[0].fT);
            }

            // revert to the previous fit result
            bool revert(false);

            // oldLogLh is only set to a value larger than 0 after the first fit,
            // so the first time this is true is with two showers, however, when
            // already cleaning up, then also do not check the log-likelihood
            // increase
            if (oldLogLh>0.) {
                if (GetCalorimeter()->GetOptions().fill_histos && GetCalorimeter()->GetOptions().histos_level>=Calorimeter::Options::HISTLVL_DEBUG) {
                    assert(fHist[h1D_fitLogLhDiff]!=0);
                    fHist[h1D_fitLogLhDiff]->Fill(oldLogLh-logLh);
                    assert(fHist[h1D_fitLogLhDiffRel]!=0);
                    fHist[h1D_fitLogLhDiffRel]->Fill((oldLogLh-logLh)/logLh);
                }

                // if the log-likelihood increase is too small, do not continue adding showers
                if ( ((oldLogLh-logLh)/logLh) < 0.02 )
                    revert=true;
            }
            if (revert) {
                shower = oldShower;
                logLh  = oldLogLh;
                ndf    = oldNdf;
                break;
            }

            // test that all showers are within the cells of the cluster
            for (std::vector<FitShower>::iterator it2=shower.begin(); it2!=shower.end(); it2++) {
                bool found(false);
                for (std::vector<FitInfo>::const_iterator it3=cells.begin(); it3!=cells.end(); it3++)
                    if (GetCalorimeter()->GetCells()[it3->GetData().GetCellIdx()].InActiveArea(it2->fX, it2->fY, GetCalorimeter()->GetCells()[it3->GetData().GetCellIdx()].GetZ())) {
                        found = true;
                        break;
                    }

                if (!found) {
                    revert=true;
                    break;
                }
            }
            if (revert) {
                shower = oldShower;
                logLh  = oldLogLh;
                ndf    = oldNdf;
                break;
            }

            // remove showers with an energy lower than the threshold
            // in case such a shower is found, also stop the fitting
            for (std::vector<FitShower>::iterator it2=shower.begin(); it2!=shower.end(); it2++)
                if (it2->fE < GetCalorimeter()->GetOptions().particle_energy_threshold) {
                    revert=true;
                    break;
                }
            if (revert) {
                shower = oldShower;
                logLh  = oldLogLh;
                ndf    = oldNdf;
                break;
            }

            // test if two showers are very close to each other, in which case
            // the previous fit result is restored and the fitting is stopped
            // close is defined as the sum of half the cell size of the cells
            // containing the shower centers
            for (std::vector<FitShower>::iterator it2=shower.begin(); it2!=shower.end(); it2++) {
                const double halfStepIt2(std::max(GetCalorimeter()->GetCells()[it2->fMainCell].GetCellType().GetStepX(),
                                                  GetCalorimeter()->GetCells()[it2->fMainCell].GetCellType().GetStepY()) / 2.);
                std::vector<FitShower>::iterator it3=it2; it3++;
                for (; it3!=shower.end(); it3++) {
                    const double distX = it2->fX - it3->fX;
                    const double distY = it2->fY - it3->fY;
                    double dist  = sqrt(distX*distX + distY*distY);
                    const double halfStepIt3(std::max(GetCalorimeter()->GetCells()[it3->fMainCell].GetCellType().GetStepX(),
                                                      GetCalorimeter()->GetCells()[it3->fMainCell].GetCellType().GetStepY()) / 2.);
                    dist /= halfStepIt2+halfStepIt3;
                    if (dist<GetCalorimeter()->GetOptions().min_shower_distance) {
                        revert=true;
                        break;
                    }
                }

                if (revert)
                    break;
            }
            if (revert) {
                shower = oldShower;
                logLh  = oldLogLh;
                ndf    = oldNdf;
                break;
            }

            // test if two showers are close to each other, and one has only a
            // small ratio of the energy of the more energetic one
            for (std::vector<FitShower>::iterator it2=shower.begin(); it2!=shower.end(); it2++) {
                const double halfStepIt2(std::max(GetCalorimeter()->GetCells()[it2->fMainCell].GetCellType().GetStepX(),
                                                  GetCalorimeter()->GetCells()[it2->fMainCell].GetCellType().GetStepY()) / 2.);
                std::vector<FitShower>::iterator it3=it2; it3++;
                for (; it3!=shower.end(); it3++) {
                    const double distX = it2->fX - it3->fX;
                    const double distY = it2->fY - it3->fY;
                    double dist  = sqrt(distX*distX + distY*distY);
                    const double halfStepIt3(std::max(GetCalorimeter()->GetCells()[it3->fMainCell].GetCellType().GetStepX(),
                                                      GetCalorimeter()->GetCells()[it3->fMainCell].GetCellType().GetStepY()) / 2.);
                    dist /= (halfStepIt2+halfStepIt3) * (it2->fE+it3->fE);
                    dist *= 2*sqrt(it2->fE*it3->fE);
                    if (dist<GetCalorimeter()->GetOptions().min_shower_distance) {
                        revert=true;
                        break;
                    }
                }

                if (revert)
                    break;
            }
            if (revert) {
                shower = oldShower;
                logLh  = oldLogLh;
                ndf    = oldNdf;
                break;
            }

            for (std::vector<FitShower>::iterator it2=shower.begin(); it2!=shower.end(); it2++) {
                std::vector<FitShower>::iterator it3=it2; it3++;
                for (; it3!=shower.end(); it3++) {
                    if (it2->fE > it3->fE) {
                        // energy of shower 2 in cell of shower 3
                        double e2(it2->fE*GetCalorimeter()->GetCells()[it3->fMainCell].GetCellType().GetShowerProfile()->GetEnergyInCell(GetCalorimeter()->GetCells()[it3->fMainCell].GetX()-it2->fX,
                                                                                                                                         GetCalorimeter()->GetCells()[it3->fMainCell].GetY()-it2->fY));
                        // energy of shower 3 in cell of shower 3
                        double e3(it3->fE*GetCalorimeter()->GetCells()[it3->fMainCell].GetCellType().GetShowerProfile()->GetEnergyInCell(GetCalorimeter()->GetCells()[it3->fMainCell].GetX()-it3->fX,
                                                                                                                                         GetCalorimeter()->GetCells()[it3->fMainCell].GetY()-it3->fY));

                        if (e3 < 5.*e2) {
                            revert = true;
                            break;
                        }
                    } else {
                        // energy of shower 2 in cell of shower 2
                        double e2(it2->fE*GetCalorimeter()->GetCells()[it2->fMainCell].GetCellType().GetShowerProfile()->GetEnergyInCell(GetCalorimeter()->GetCells()[it2->fMainCell].GetX()-it2->fX,
                                                                                                                                         GetCalorimeter()->GetCells()[it2->fMainCell].GetY()-it2->fY));
                        // energy of shower 3 in cell of shower 2
                        double e3(it3->fE*GetCalorimeter()->GetCells()[it2->fMainCell].GetCellType().GetShowerProfile()->GetEnergyInCell(GetCalorimeter()->GetCells()[it2->fMainCell].GetX()-it3->fX,
                                                                                                                                         GetCalorimeter()->GetCells()[it2->fMainCell].GetY()-it3->fY));

                        if (e2 < 5.*e3) {
                            revert = true;
                            break;
                        }
                    }
                }

                if (revert)
                    break;
            }
            if (revert) {
                shower = oldShower;
                logLh  = oldLogLh;
                ndf    = oldNdf;
                break;
            }

            // not reverted to previous fit, so prepare for the next fit:
            // - save the current fit status
            oldShower = shower;
            oldLogLh  = logLh;
            oldNdf    = ndf;
        }

        // copy cell indices and original energies to be exported to a CalorimeterParticle
        std::vector<std::pair<size_t, double> > clusterData;
        for (std::vector<CellDataRaw>::const_iterator it2=it->GetCluster().begin(); it2!=it->GetCluster().end(); it2++) {
            clusterData.push_back(std::pair<size_t, double>(it2->GetCellIdx(), it2->GetEnergy()));
        }

        // calculate R4
        // ratio of the highest energetic 2x2 cells containing the most-
        // energetic cell
        //
        //  -+-----+-----+-----+-
        //   |ul   |ul ur|   ur|
        //  -+-----+-----+-----+-
        //   |ul ll|  *  |ur lr|
        //  -+-----+-----+-----+-
        //   |ll   |ll lr|   lr|
        //  -+-----+-----+-----+-
        //
        //  the asterisk marks the highest energtic cell
        //
        size_t maxCell = it->GetMaxCell();
        double r4(-1.);
        double r9(-1.);

        // allow for some imprecision of the positions calculated for the cell positions
        // also allow for this imprecision to check whether all cells in a cluster have
        // the same size
        const double epsilon(0.01);

        if (GetCalorimeter()->XYRegularGrid() || it->IsXYRegular(epsilon)) {
            assert(mapCellsIndices.find(maxCell)!=mapCellsIndices.end());
            double r4_ul(0.);
            double r4_ur(0.);
            double r4_ll(0.);
            double r4_lr(0.);

            r9 = 0;

            const Cell& maxCellC = GetCalorimeter()->GetCells()[maxCell];
            const std::vector<size_t>& neighbors = maxCellC.GetNeighbors();
            for (std::vector<size_t>::const_iterator it2=neighbors.begin(); it2!=neighbors.end(); it2++) {
                if (mapCellsIndices.find(*it2)!=mapCellsIndices.end()) {
                    if (cells[mapCellsIndices[*it2]].GetArtificial())
                        continue;

                    const Cell& cell = GetCalorimeter()->GetCells()[*it2];

                    if (cell.GetX() <= maxCellC.GetX()+epsilon && cell.GetY() <= maxCellC.GetY()+epsilon)
                        r4_ll += cells[mapCellsIndices[*it2]].GetData().GetEnergy();
                    if (cell.GetX() <= maxCellC.GetX()+epsilon && cell.GetY() >= maxCellC.GetY()-epsilon)
                        r4_ul += cells[mapCellsIndices[*it2]].GetData().GetEnergy();
                    if (cell.GetX() >= maxCellC.GetX()-epsilon && cell.GetY() <= maxCellC.GetY()+epsilon)
                        r4_lr += cells[mapCellsIndices[*it2]].GetData().GetEnergy();
                    if (cell.GetX() >= maxCellC.GetX()-epsilon && cell.GetY() >= maxCellC.GetY()-epsilon)
                        r4_ur += cells[mapCellsIndices[*it2]].GetData().GetEnergy();

                    r9 += cells[mapCellsIndices[*it2]].GetData().GetEnergy();
                }
            }

            r4  = std::max(std::max(r4_ul, r4_ur), std::max(r4_ll, r4_lr));
            r4 /= it->AmplitudeTotal();

            r9 /= it->AmplitudeTotal();
        }

        if (shower.empty()) {
            // no fit done or fit too bad to be used
            // use center of gravity instead
            double eTot(0.);
            double eVar(0.);
            double time(0.);
            double timeErr(0.);
            for (std::vector<FitInfo>::const_iterator it2=cells.begin(); it2!=cells.end(); it2++) {
                if (it2->GetArtificial())
                    continue;

                eTot+=it2->GetData().GetEnergy();
                eVar+=it2->GetData().GetEnergyErr()*it2->GetData().GetEnergyErr();

                time   +=it2->GetData().GetTime() / (it2->GetData().GetTimeErr()*it2->GetData().GetTimeErr());
                timeErr+=            1.           / (it2->GetData().GetTimeErr()*it2->GetData().GetTimeErr());
            }

            // only save this cluster if it has an energy larger than given in the options,
            // otherwise go to the next cluster
            if (eTot < GetCalorimeter()->GetOptions().particle_energy_threshold)
                continue;

            eVar   =sqrt(eVar);
            time   =time/timeErr;
            timeErr=1./sqrt(timeErr);

            FitShower shower(eTot, it->MeanX(), it->MeanY(), time,
                             eVar, it->VarX(),  it->VarY(),  timeErr);
            shower.fMainCell = it->GetMaxCell();
            FitShower corrShower(ApplyCorrections(shower));

            // taken analoguesly to ReconstructionLednev.cc
            // TODO: z position only correct for electromagnetic showers
            const double zSurface = GetCalorimeter()->GetCells()[corrShower.fMainCell].GetZ() - GetCalorimeter()->GetCells()[corrShower.fMainCell].GetCellType().GetSizeZ() / 2.;
            const double zShower  = ZmidShowerEM(corrShower.fE) * GetCalorimeter()->GetCells()[corrShower.fMainCell].GetCellType().GetRadiationLength();
            const double z        = zSurface + zShower;
            const double zErr     = 50.;

            reconstructedParticles.push_back(CalorimeterParticle(GetCalorimeter()->GetOptions().particle_default, 1,
                                                                 corrShower.fE,    corrShower.fX,    corrShower.fY,    z,
                                                                 GetCalorimeter(),
                                                                 corrShower.fEerr, corrShower.fXerr, corrShower.fYerr, zErr,
                                                                 0, 0));
            reconstructedParticles.back().SetHitedCell(corrShower.fMainCell);
            reconstructedParticles.back().SetMiscInfo(CalorimeterParticle::VALUE_R4, r4);
            reconstructedParticles.back().SetMiscInfo(CalorimeterParticle::VALUE_R9, r9);
            reconstructedParticles.back().SetMiscInfo(CalorimeterParticle::CLUSTER_SHOWERS, 1);
            reconstructedParticles.back().SetMiscInfo(CalorimeterParticle::UNCORR_X, it->MeanX());
            reconstructedParticles.back().SetMiscInfo(CalorimeterParticle::UNCORR_Y, it->MeanY());
            reconstructedParticles.back().SetMiscInfo(CalorimeterParticle::UNCORR_E, eTot);
            reconstructedParticles.back().SetClusterData(clusterData);
            reconstructedParticles.back().SetTime(corrShower.fT, corrShower.fTerr);

            if (GetCalorimeter()->GetOptions().fill_histos) {
                assert(fHist[h1D_showerEnergy]!=0);
                fHist[h1D_showerEnergy]->Fill(eTot);
                assert(fHist[h1D_showerTime]!=0);
                fHist[h1D_showerTime]->Fill(time);
                assert(fHist[h2D_showerPosition]!=0);
                fHist[h2D_showerPosition]->Fill(it->MeanX(), it->MeanY());

                if (GetCalorimeter()->GetOptions().histos_level>=Calorimeter::Options::HISTLVL_VERBOSE) {
                    assert(fHist[h1D_clusterSize]!=0);
                    fHist[h1D_clusterSize]->Fill(length);
                    assert(fHist[h1D_clusterR4]!=0);
                    fHist[h1D_clusterR4]->Fill(r4);
                    assert(fHist[h1D_clusterR9]!=0);
                    fHist[h1D_clusterR9]->Fill(r9);
                }
            }
        } else {
            double eTotShower(0.);
            for (std::vector<FitShower>::const_iterator it2=shower.begin(); it2!=shower.end(); it2++) {
                eTotShower+=it2->fE;

                FitShower corrShower(ApplyCorrections(*it2));

                // taken analoguesly to ReconstructionLednev.cc
                // TODO: z position only correct for electromagnetic showers
                const double zSurface = GetCalorimeter()->GetCells()[corrShower.fMainCell].GetZ() - GetCalorimeter()->GetCells()[corrShower.fMainCell].GetCellType().GetSizeZ() / 2.;
                const double zShower  = ZmidShowerEM(corrShower.fE) * GetCalorimeter()->GetCells()[corrShower.fMainCell].GetCellType().GetRadiationLength();
                const double z        = zSurface + zShower;
                const double zErr     = 50.;

                reconstructedParticles.push_back(CalorimeterParticle(GetCalorimeter()->GetOptions().particle_default, 1,
                                                                     corrShower.fE,    corrShower.fX,    corrShower.fY,    z,
                                                                     GetCalorimeter(),
                                                                     corrShower.fEerr, corrShower.fXerr, corrShower.fYerr, zErr,
                                                                     0, 0));
                reconstructedParticles.back().SetHitedCell(corrShower.fMainCell);
                reconstructedParticles.back().SetMiscInfo(CalorimeterParticle::VALUE_R4, r4);
                reconstructedParticles.back().SetMiscInfo(CalorimeterParticle::VALUE_R9, r9);
                reconstructedParticles.back().SetMiscInfo(CalorimeterParticle::CLUSTER_SHOWERS, shower.size());
                reconstructedParticles.back().SetMiscInfo(CalorimeterParticle::CHI_GAM, logLh);
                reconstructedParticles.back().SetMiscInfo(CalorimeterParticle::NDF_GAM, ndf);
                reconstructedParticles.back().SetMiscInfo(CalorimeterParticle::UNCORR_X, it2->fX);
                reconstructedParticles.back().SetMiscInfo(CalorimeterParticle::UNCORR_Y, it2->fY);
                reconstructedParticles.back().SetMiscInfo(CalorimeterParticle::UNCORR_E, it2->fE);
                reconstructedParticles.back().SetClusterData(clusterData);
                reconstructedParticles.back().SetTime(corrShower.fT, corrShower.fTerr);

                if (GetCalorimeter()->GetOptions().fill_histos) {
                    assert(fHist[h1D_showerEnergy]!=0);
                    fHist[h1D_showerEnergy]->Fill(it2->fE);
                    assert(fHist[h1D_showerTime]!=0);
                    fHist[h1D_showerTime]->Fill(it2->fT);
                    assert(fHist[h2D_showerEnergyVsTime]!=0);
                    fHist[h2D_showerEnergyVsTime]->Fill(it2->fE, it2->fT);
                    assert(fHist[h2D_showerPosition]!=0);
                    fHist[h2D_showerPosition]->Fill(it2->fX, it2->fY);

                    // fill distance histogram
                    std::vector<FitShower>::const_iterator it3=it2; it3++;
                    for (; it3!=shower.end(); it3++) {
                        double distX = it2->fX-it3->fX;
                        double distY = it2->fY-it3->fY;
                        double dist  = sqrt(distX*distX + distY*distY);
                        assert(fHist[h1D_showerDistance]!=0);
                        fHist[h1D_showerDistance]->Fill(dist);
                        assert(fHist[h2D_showerDistance]!=0);
                        fHist[h2D_showerDistance]->Fill(fabs(distX), fabs(distY));

                        double ratio(-1);
                        if (it2->fE > it3->fE)
                            ratio = it3->fE / it2->fE;
                        else
                            ratio = it2->fE / it3->fE;
                        assert(fHist[h2D_showerDistanceVsEnergyRatio]!=0);
                        fHist[h2D_showerDistanceVsEnergyRatio]->Fill(dist, ratio);

                        if (GetCalorimeter()->GetOptions().histos_level>=Calorimeter::Options::HISTLVL_VERBOSE) {
                            assert(fHist[h1D_showerTimeDiff]!=0);
                            fHist[h1D_showerTimeDiff]->Fill(it2->fT - it3->fT);
                            assert(fHist[h2D_showerTimeDiffVsTime]!=0);
                            fHist[h2D_showerTimeDiffVsTime]->Fill(it2->fT, it2->fT - it3->fT);
                            fHist[h2D_showerTimeDiffVsTime]->Fill(it3->fT, it3->fT - it2->fT);
                        }
                    }

                    if (GetCalorimeter()->GetOptions().histos_level>=Calorimeter::Options::HISTLVL_VERBOSE) {
                        assert(fHist[h1D_clusterSize]!=0);
                        fHist[h1D_clusterSize]->Fill(length);
                        assert(fHist[h1D_clusterR4]!=0);
                        fHist[h1D_clusterR4]->Fill(r4);
                        assert(fHist[h1D_clusterR9]!=0);
                        fHist[h1D_clusterR9]->Fill(r9);

                        assert(fHist[h2D_showerEnergyLowVsTime]!=0);
                        fHist[h2D_showerEnergyLowVsTime]->Fill(it2->fE, it2->fT);

                        assert(fHist[h1D_showerResiduumX]!=0);
                        fHist[h1D_showerResiduumX]->Fill(it2->fX - it->MeanX());
                        assert(fHist[h2D_showerResiduumXVsEnergy]!=0);
                        fHist[h2D_showerResiduumXVsEnergy]->Fill(it2->fX - it->MeanX(), it2->fE);
                        assert(fHist[h1D_showerResiduumY]!=0);
                        fHist[h1D_showerResiduumY]->Fill(it2->fY - it->MeanY());
                        assert(fHist[h2D_showerResiduumYVsEnergy]!=0);
                        fHist[h2D_showerResiduumYVsEnergy]->Fill(it2->fY - it->MeanY(), it2->fE);

                        assert(fHist[h2D_showerEnergyVsShowerCount]!=0);
                        fHist[h2D_showerEnergyVsShowerCount]->Fill(it2->fE, shower.size());
                        assert(fHist[h2D_showerEnergyVsClusterSize]!=0);
                        fHist[h2D_showerEnergyVsClusterSize]->Fill(it2->fE, length);
                        assert(fHist[h2D_showerEnergyLowVsShowerCount]!=0);
                        fHist[h2D_showerEnergyLowVsShowerCount]->Fill(it2->fE, shower.size());
                        assert(fHist[h2D_showerEnergyLowVsClusterSize])    ;
                        fHist[h2D_showerEnergyLowVsClusterSize]->Fill(it2->fE, length);
                    }
                }
            }

            // do not fill the next histograms in cases where there has been no
            // fit, but still showers were created: the shower parameters have
            // then been calculated(!, not fitted) from the cell information.
            // i.e. one cell clusters, ...
            if (logLh >= 0. && GetCalorimeter()->GetOptions().fill_histos) {
                assert(fHist[h1D_showerLogLh]!=0);
                fHist[h1D_showerLogLh]->Fill(logLh);
                assert(fHist[h1D_showerNdF]!=0);
                fHist[h1D_showerNdF]->Fill(ndf);
                assert(fHist[h1D_showerLogLhOverNdF]!=0);
                fHist[h1D_showerLogLhOverNdF]->Fill(logLh/ndf);
                if (GetCalorimeter()->GetOptions().histos_level>=Calorimeter::Options::HISTLVL_VERBOSE) {
                    assert(fHist[h1D_showerLogLhOverNdFLarge]!=0);
                    fHist[h1D_showerLogLhOverNdFLarge]->Fill(logLh/ndf);
                    assert(fHist[h1D_showerProb]!=0);
                    fHist[h1D_showerProb]->Fill(TMath::Prob(logLh, ndf));
                    assert(fHist[h2D_showerLogLhOverNdFVsShowerCount]!=0);
                    fHist[h2D_showerLogLhOverNdFVsShowerCount]->Fill(logLh/ndf, shower.size());
                    assert(fHist[h2D_showerLogLhOverNdFVsClusterSize]!=0);
                    fHist[h2D_showerLogLhOverNdFVsClusterSize]->Fill(logLh/ndf, length);
                }

                assert(fHist[h1D_showerCount]!=0);
                fHist[h1D_showerCount]->Fill(shower.size());

                if (GetCalorimeter()->GetOptions().histos_level>=Calorimeter::Options::HISTLVL_VERBOSE) {
                    assert(fHist[h2D_showerCountVsClusterSize]!=0);
                    fHist[h2D_showerCountVsClusterSize]->Fill(length, shower.size());

                    assert(fHist[h1D_showerEnergyDesc]!=0);
                    fHist[h1D_showerEnergyDesc]->Fill(eTotShower - it->AmplitudeTotal());
                    assert(fHist[h2D_showerEnergyDescVsEnergy]!=0);
                    fHist[h2D_showerEnergyDescVsEnergy]->Fill(eTotShower - it->AmplitudeTotal(), it->AmplitudeTotal());
                    assert(fHist[h2D_showerEnergyDescVsEnergyLow]!=0);
                    fHist[h2D_showerEnergyDescVsEnergyLow]->Fill(eTotShower - it->AmplitudeTotal(), it->AmplitudeTotal());

                    double eTotShowerInCluster(0.);
                    for (std::vector<FitInfo>::const_iterator it2=cells.begin(); it2!=cells.end(); it2++) {
                        if (it2->GetArtificial())
                            continue;

                        eTotShowerInCluster += it2->GetFittedEnergy(shower);
                    }
                    assert(fHist[h1D_showerEnergyDescCells]!=0);
                    fHist[h1D_showerEnergyDescCells]->Fill(eTotShowerInCluster - it->AmplitudeTotal());
                    assert(fHist[h2D_showerEnergyDescCellsVsEnergy]!=0);
                    fHist[h2D_showerEnergyDescCellsVsEnergy]->Fill(eTotShowerInCluster - it->AmplitudeTotal(), it->AmplitudeTotal());
                    assert(fHist[h2D_showerEnergyDescCellsVsEnergyLow]!=0);
                    fHist[h2D_showerEnergyDescCellsVsEnergyLow]->Fill(eTotShowerInCluster - it->AmplitudeTotal(), it->AmplitudeTotal());

                    if (GetCalorimeter()->GetOptions().histos_level>=Calorimeter::Options::HISTLVL_DEBUG) {
                        // fill pulls histograms
                        for (std::vector<FitInfo>::const_iterator it2=cells.begin(); it2!=cells.end(); it2++) {
                            if (it2->GetArtificial())
                                continue;

                            assert(fHist[h2D_fitPullsCellTypeFirst+fMapCellType2HistIndex[it2->GetCellType().GetHwType()]]!=0);
                            fHist[h2D_fitPullsCellTypeFirst+fMapCellType2HistIndex[it2->GetCellType().GetHwType()]]
                                    ->Fill(it2->GetData().GetEnergy(),
                                           (it2->GetFittedEnergy(shower)-it2->GetData().GetEnergy()) / it2->GetData().GetEnergy());
                            assert(fHist[h2D_fitPullsLowCellTypeFirst+fMapCellType2HistIndex[it2->GetCellType().GetHwType()]]!=0);
                            fHist[h2D_fitPullsLowCellTypeFirst+fMapCellType2HistIndex[it2->GetCellType().GetHwType()]]
                                    ->Fill(it2->GetData().GetEnergy(),
                                           (it2->GetFittedEnergy(shower)-it2->GetData().GetEnergy()) / it2->GetData().GetEnergy());
                            assert(fHist[h2D_fitNormPullsCellTypeFirst+fMapCellType2HistIndex[it2->GetCellType().GetHwType()]]!=0);
                            fHist[h2D_fitNormPullsCellTypeFirst+fMapCellType2HistIndex[it2->GetCellType().GetHwType()]]
                                    ->Fill(it2->GetData().GetEnergy(),
                                           (it2->GetFittedEnergy(shower)-it2->GetData().GetEnergy()) / it2->GetData().GetEnergyErr());
                            assert(fHist[h2D_fitNormPullsLowCellTypeFirst+fMapCellType2HistIndex[it2->GetCellType().GetHwType()]]!=0);
                            fHist[h2D_fitNormPullsLowCellTypeFirst+fMapCellType2HistIndex[it2->GetCellType().GetHwType()]]
                                    ->Fill(it2->GetData().GetEnergy(),
                                           (it2->GetFittedEnergy(shower)-it2->GetData().GetEnergy()) / it2->GetData().GetEnergyErr());
                        }
                    }
                }
            } // end of logLh >= 0.
        }
    }

    return reconstructedParticles;
}

/////////////////////////////////////////////////////////////////////////////

const std::vector<CalorimeterParticle>& ReconstructionCombined::DoRepeatReconstruction(const std::vector<CellDataRaw>& signals) {
    throw Exception("Using non implemented method!");
}

/////////////////////////////////////////////////////////////////////////////

void ReconstructionCombined::BookHistograms() {
    // prepare mapping of cell types to entries in histogram array
    size_t count(0);
    fMapCellType2HistIndex.clear();
    for (std::vector<Cell>::const_iterator it=GetCalorimeter()->GetCells().begin(); it!=GetCalorimeter()->GetCells().end(); it++)
        if (fMapCellType2HistIndex.count(it->GetCellType().GetHwType())==0) {
            if (count==nrCellTypesForHistograms) {
                std::cerr << "ReconstructionCombined::BookHistograms: " << GetCalorimeter()->GetName()
                          << ": Maximal number of different cell types (" << nrCellTypesForHistograms << ") exceeded, increase nrCellTypesForHistograms in the source code!" << std::endl;
                exit(1);
                // TODO: should be fatal exception, but is caught somewhere
                throw Exception("ReconstructionCombined::BookHistograms: %s: Maximal number of different cell types (%zu) exceeded, increase nrCellTypesForHistograms in the source code!",
                                GetCalorimeter()->GetName().c_str(), nrCellTypesForHistograms);
            }

            fMapCellType2HistIndex.insert(std::pair<CellType::HwType, size_t>(it->GetCellType().GetHwType(), count));
            count++;
        }

    if (GetCalorimeter()->GetOptions().fill_histos) {
        TDirectory* base     = GetCalorimeter()->GetHistogramsBaseDir();
        TDirectory* recoBase = base->mkdir((GetCalorimeter()->GetName() + "_reconstruction").c_str(),
                (std::string("Histograms for reconstruction of ") + GetCalorimeter()->GetName()).c_str());

        recoBase->cd();

        // individual cell statistics
        fHist[h1D_cellEnergy] = new TH1F((GetCalorimeter()->GetName() + "_cellEnergy").c_str(),
                "Energy for individual cell",
                300, 0., 300.);
        fHist[h1D_cellTime] = new TH1F((GetCalorimeter()->GetName() + "_cellTime").c_str(),
                "Time for individual cell",
                200, -50., 50.);
        fHist[h2D_cellEnergyVsTime] = new TH2F((GetCalorimeter()->GetName() + "_cellEnergyVsTime").c_str(),
                "Energy versus time for individual cell",
                300, 0., 300., 200, -50., 50.);
        fHist[h2D_cellEnergyVsTime]->SetOption("COLZ");

        if (GetCalorimeter()->GetOptions().histos_level>=Calorimeter::Options::HISTLVL_VERBOSE) {
            fHist[h1D_cellLowEnergy] = new TH1F((GetCalorimeter()->GetName() + "_cellLowEnergy").c_str(),
                    "Energy for individual cell",
                    200, 0., 10.);
            fHist[h2D_cellLowEnergyVsTime] = new TH2F((GetCalorimeter()->GetName() + "_cellLowEnergyVsTime").c_str(),
                    "Energy versus time for individual cell",
                    200, 0., 10., 200, -50., 50.);
            fHist[h2D_cellLowEnergyVsTime]->SetOption("COLZ");

            if (GetCalorimeter()->GetOptions().histos_level>=Calorimeter::Options::HISTLVL_DEBUG) {
                fHist[h1D_cellEnergyThreshold] = new TH1F((GetCalorimeter()->GetName() + "_cellEnergyThreshold").c_str(),
                        "Energy threshold for each cell (in GeV)",
                        GetCalorimeter()->GetCells().size(), -0.5, GetCalorimeter()->GetCells().size()-0.5);
                fHist[h1D_cellEnergyThresholdAdc] = new TH1F((GetCalorimeter()->GetName() + "_cellEnergyThresholdAdc").c_str(),
                        "Energy threshold for each cell (in ADC channels)",
                        GetCalorimeter()->GetCells().size(), -0.5, GetCalorimeter()->GetCells().size()-0.5);
            }
        }

        // cluster statistics
        if (GetCalorimeter()->GetOptions().histos_level>=Calorimeter::Options::HISTLVL_VERBOSE) {
            fHist[h1D_clusterSize] = new TH1F((GetCalorimeter()->GetName() + "_clusterSize").c_str(),
                    "Cluster size",
                    50, -0.5, 49.5);
            fHist[h1D_clusterR4] = new TH1F((GetCalorimeter()->GetName() + "_clusterR4").c_str(),
                    "Cluster R4 value (overflow: R4==1, underflow: R4==-1)",
                    50, 0., 1.);
            fHist[h1D_clusterR9] = new TH1F((GetCalorimeter()->GetName() + "_clusterR9").c_str(),
                    "Cluster R9 value (overflow: R9==1, underflow: R9==-1)",
                    50, 0., 1.);
        }

        // shower statistics
        fHist[h1D_showerCount] = new TH1F((GetCalorimeter()->GetName() + "_showerCount").c_str(),
                "Number of reconstructed showers in one cluster",
                10, -0.5, 9.5);
        fHist[h1D_showerEnergy] = new TH1F((GetCalorimeter()->GetName() + "_showerEnergy").c_str(),
                "Energy in reconstructed shower",
                300, 0., 300.);
        fHist[h1D_showerTime] = new TH1F((GetCalorimeter()->GetName() + "_showerTime").c_str(),
                "Time for reconstructed shower",
                200, -50., 50.);
        fHist[h2D_showerPosition] = new TH2F((GetCalorimeter()->GetName() + "_showerPosition").c_str(),
                "Position of showers",
                200, -2000., 2000., 150, -1500., 1500.);
        fHist[h2D_showerPosition]->SetOption("COLZ");

        fHist[h2D_showerEnergyVsTime] = new TH2F((GetCalorimeter()->GetName() + "_showerEnergyVsTime").c_str(),
                "shower energy vs time",
                200, 0., 200., 200, -50., 50.);
        fHist[h2D_showerEnergyVsTime]->SetOption("COLZ");

        fHist[h1D_showerLogLh] = new TH1F((GetCalorimeter()->GetName() + "_showerLogLh").c_str(),
                "log-likelihood of reconstructed shower (one entry per cluster)",
                200, 0., 200.);
        fHist[h1D_showerNdF] = new TH1F((GetCalorimeter()->GetName() + "_showerNdF").c_str(),
                "NdF of reconstructed shower (one entry per cluster)",
                201, -0.5, 200.5);
        fHist[h1D_showerLogLhOverNdF] = new TH1F((GetCalorimeter()->GetName() + "_showerLogLhOverNdF").c_str(),
                "log-likelihood/NdF of reconstructed shower (one entry per cluster)",
                200, 0., 10.);

        fHist[h1D_showerDistance] = new TH1F((GetCalorimeter()->GetName() + "_showerDistance").c_str(),
                "Distance between a pair of showers in the same cluster (in mm)",
                200, 0., 200.);
        fHist[h2D_showerDistance] = new TH2F((GetCalorimeter()->GetName() + "_showerDistance2D").c_str(),
                "Distance between a pair of showers in the same cluster (in mm)",
                200, 0., 200., 200, 0., 200.);
        fHist[h2D_showerDistance]->SetOption("COLZ");
        fHist[h2D_showerDistanceVsEnergyRatio] = new TH2F((GetCalorimeter()->GetName() + "_showerDistanceVsEnergyRatio").c_str(),
                "Distance between a pair of showers in the same cluster vs. ratio of energies",
                200, 0., 200., 200, 0., 1.);
        fHist[h2D_showerDistanceVsEnergyRatio]->SetOption("COLZ");

        if (GetCalorimeter()->GetOptions().histos_level>=Calorimeter::Options::HISTLVL_VERBOSE) {
            fHist[h1D_showerLogLhOverNdFLarge] = new TH1F((GetCalorimeter()->GetName() + "_showerLogLhOverNdFLarge").c_str(),
                    "log-likelihood/NdF of reconstructed shower (one entry per cluster)",
                    200, 0., 100.);
            fHist[h1D_showerProb] = new TH1F((GetCalorimeter()->GetName() + "_showerProb").c_str(),
                    "probability of fit for reconstructed showers (one entry per cluster)",
                    200, 0., 1.);

            fHist[h1D_showerEnergyDesc] = new TH1F((GetCalorimeter()->GetName() + "_showerEnergyDesc").c_str(),
                    "Difference in energy between all showers and cluster (showers - cells)",
                    200, -20., 20.);
            fHist[h1D_showerEnergyDesc]->GetXaxis()->SetTitle("sum of shower energies - cluster energy");
            fHist[h2D_showerEnergyDescVsEnergy] = new TH2F((GetCalorimeter()->GetName() + "_showerEnergyDescVsEnergy").c_str(),
                    "Difference in energy between all showers and cluster (showers - cells) vs. energy",
                    200, -20., 20., 200, 0., 200.);
            fHist[h2D_showerEnergyDescVsEnergy]->GetXaxis()->SetTitle("sum of shower energies - cluster energy");
            fHist[h2D_showerEnergyDescVsEnergy]->GetYaxis()->SetTitle("cluster energy");
            fHist[h2D_showerEnergyDescVsEnergy]->SetOption("COLZ");
            fHist[h2D_showerEnergyDescVsEnergyLow] = new TH2F((GetCalorimeter()->GetName() + "_showerEnergyDescVsEnergyLow").c_str(),
                    "Difference in energy between all showers and cluster (showers - cells) vs. energy",
                    200, -5., 5., 200, 0., 10.);
            fHist[h2D_showerEnergyDescVsEnergyLow]->GetXaxis()->SetTitle("sum of shower energies - cluster energy");
            fHist[h2D_showerEnergyDescVsEnergyLow]->GetYaxis()->SetTitle("cluster energy");
            fHist[h2D_showerEnergyDescVsEnergyLow]->SetOption("COLZ");

            fHist[h1D_showerEnergyDescCells] = new TH1F((GetCalorimeter()->GetName() + "_showerEnergyDescCells").c_str(),
                    "Difference in energy between shower energy in cluster cells and cluster energy",
                    200, -20., 20.);
            fHist[h1D_showerEnergyDescCells]->GetXaxis()->SetTitle("sum of shower energies in cluster cells - cluster energy");
            fHist[h2D_showerEnergyDescCellsVsEnergy] = new TH2F((GetCalorimeter()->GetName() + "_showerEnergyDescCellsVsEnergy").c_str(),
                    "Difference in energy between shower energy in cluster cells and cluster energy vs. energy",
                    200, -20., 20., 200, 0., 200.);
            fHist[h2D_showerEnergyDescCellsVsEnergy]->GetXaxis()->SetTitle("sum of shower energies in cluster cells - cluster energy");
            fHist[h2D_showerEnergyDescCellsVsEnergy]->GetYaxis()->SetTitle("cluster energy");
            fHist[h2D_showerEnergyDescCellsVsEnergy]->SetOption("COLZ");
            fHist[h2D_showerEnergyDescCellsVsEnergyLow] = new TH2F((GetCalorimeter()->GetName() + "_showerEnergyDescCellsVsEnergyLow").c_str(),
                    "Difference in energy between shower energy in cluster cells and cluster energy vs. energy",
                    200, -20., 20., 200, 0., 10.);
            fHist[h2D_showerEnergyDescCellsVsEnergyLow]->GetXaxis()->SetTitle("sum of shower energies in cluster cells - cluster energy");
            fHist[h2D_showerEnergyDescCellsVsEnergyLow]->GetYaxis()->SetTitle("cluster energy");
            fHist[h2D_showerEnergyDescCellsVsEnergyLow]->SetOption("COLZ");

            fHist[h1D_showerResiduumX] = new TH1F((GetCalorimeter()->GetName() + "_showerResiduumX").c_str(),
                    "Difference in shower position and weighted mean in X (shower - mean)",
                    200, -10., 10.);
            fHist[h2D_showerResiduumXVsEnergy] = new TH2F((GetCalorimeter()->GetName() + "_showerResiduumXVsEnergy").c_str(),
                    "Difference in shower position and weighted mean in X versus shower energy",
                    200, -10., 10., 200, 0., 200.);
            fHist[h2D_showerResiduumXVsEnergy]->SetOption("COLZ");
            fHist[h1D_showerResiduumY] = new TH1F((GetCalorimeter()->GetName() + "_showerResiduumY").c_str(),
                    "Difference in shower position and weighted mean in Y (shower - mean)",
                    200, -10., 10.);
            fHist[h2D_showerResiduumYVsEnergy] = new TH2F((GetCalorimeter()->GetName() + "_showerResiduumYVsEnergy").c_str(),
                    "Difference in shower position and weighted mean in Y versus shower energy",
                    200, -10., 10., 200, 0., 200.);
            fHist[h2D_showerResiduumYVsEnergy]->SetOption("COLZ");

            fHist[h2D_showerLogLhOverNdFVsShowerCount] = new TH2F((GetCalorimeter()->GetName() + "_showerLogLhOverNdFVsShowerCount").c_str(),
                    "log-likelihood/NdF of reconstructed shower vs number of showers (one entry per cluster)",
                    200, 0., 100., 10, -0.5, 9.5);
            fHist[h2D_showerLogLhOverNdFVsShowerCount]->SetOption("COLZ");

            fHist[h2D_showerLogLhOverNdFVsClusterSize] = new TH2F((GetCalorimeter()->GetName() + "_showerLogLhOverNdFVsClusterSize").c_str(),
                    "log-likelihood/NdF of reconstructed shower vs cluster size (one entry per cluster)",
                    200, 0., 100., 50, -0.5, 49.5);
            fHist[h2D_showerLogLhOverNdFVsClusterSize]->SetOption("COLZ");

            fHist[h2D_showerCountVsClusterSize] = new TH2F((GetCalorimeter()->GetName() + "_showerCountVsClusterSize").c_str(),
                    "number of showers vs cluster size (one entry per cluster)",
                    50, -0.5, 49.5, 10, -0.5, 9.5);
            fHist[h2D_showerCountVsClusterSize]->SetOption("COLZ");

            fHist[h2D_showerEnergyVsShowerCount] = new TH2F((GetCalorimeter()->GetName() + "_showerEnergyVsShowerCount").c_str(),
                    "shower energy vs number of showers",
                    200, 0., 200., 10, -0.5, 9.5);
            fHist[h2D_showerEnergyVsShowerCount]->SetOption("COLZ");
            fHist[h2D_showerEnergyVsClusterSize] = new TH2F((GetCalorimeter()->GetName() + "_showerEnergyVsClusterSize").c_str(),
                    "shower energy vs cluster size",
                    200, 0., 200., 50, -0.5, 49.5);
            fHist[h2D_showerEnergyVsClusterSize]->SetOption("COLZ");
            fHist[h2D_showerEnergyLowVsShowerCount] = new TH2F((GetCalorimeter()->GetName() + "_showerEnergyLowVsShowerCount").c_str(),
                    "shower energy vs number of showers",
                    200, 0., 10., 10, -0.5, 9.5);
            fHist[h2D_showerEnergyLowVsShowerCount]->SetOption("COLZ");
            fHist[h2D_showerEnergyLowVsClusterSize] = new TH2F((GetCalorimeter()->GetName() + "_showerEnergyLowVsClusterSize").c_str(),
                    "shower energy vs cluster size",
                    200, 0., 10., 50, -0.5, 49.5);
            fHist[h2D_showerEnergyLowVsClusterSize]->SetOption("COLZ");

            fHist[h2D_showerEnergyLowVsTime] = new TH2F((GetCalorimeter()->GetName() + "_showerLowEnergyVsTime").c_str(),
                    "shower energy vs time",
                    200, 0., 10., 200, -50., 50.);
            fHist[h2D_showerEnergyLowVsTime]->SetOption("COLZ");

            fHist[h1D_showerTimeDiff] = new TH1F((GetCalorimeter()->GetName() + "_showerTimeDiff").c_str(),
                    "difference of shower times in the same cluster",
                    400, -100., 100.);
            fHist[h2D_showerTimeDiffVsTime] = new TH2F((GetCalorimeter()->GetName() + "_showerTimeDiffVsTime").c_str(),
                    "difference of shower times in the same cluster vs shower time",
                    200, -50., 50., 400, -100., 100.);
            fHist[h2D_showerTimeDiffVsTime]->SetOption("COLZ");


            // fit internals
            if (GetCalorimeter()->GetOptions().histos_level>=Calorimeter::Options::HISTLVL_DEBUG) {
                fHist[h1D_fitLogLhDiff] = new TH1F((GetCalorimeter()->GetName() + "_fitLogLhDiff").c_str(),
                        "difference in log-likelihood by adding an addition shower (old-new)",
                        200, -100, 100);
                fHist[h1D_fitLogLhDiffRel] = new TH1F((GetCalorimeter()->GetName() + "_fitLogLhDiffRel").c_str(),
                        "relative difference in log-likelihood by adding an addition shower (old-new)/new",
                        300, -1., 2.);

                fHist[h1D_fitDiffEnergyFirst] = new TH1F((GetCalorimeter()->GetName() + "_fitDiffEnergyFirst").c_str(),
                        "difference between start values and fit result for energy of first shower",
                        200, -10., 10.);
                fHist[h1D_fitDiffEnergyFirst]->GetXaxis()->SetTitle("start value - fit result");
                fHist[h2D_fitDiffEnergyFirstVsEnergy] = new TH2F((GetCalorimeter()->GetName() + "_fitDiffEnergyFirstVsEnergy").c_str(),
                        "difference between start values and fit result for energy of first shower vs. shower energy after fit",
                        200, -10., 10., 200, 0., 200.);
                fHist[h2D_fitDiffEnergyFirstVsEnergy]->GetXaxis()->SetTitle("start value - fit result");
                fHist[h2D_fitDiffEnergyFirstVsEnergy]->SetOption("COLZ");
                fHist[h2D_fitDiffEnergyFirstVsEnergyLow] = new TH2F((GetCalorimeter()->GetName() + "_fitDiffEnergyFirstVsEnergyLow").c_str(),
                        "difference between start values and fit result for energy of first shower vs. shower energy after fit",
                        200, -10., 10., 200, 0., 10.);
                fHist[h2D_fitDiffEnergyFirstVsEnergyLow]->GetXaxis()->SetTitle("start value - fit result");
                fHist[h2D_fitDiffEnergyFirstVsEnergyLow]->SetOption("COLZ");
                fHist[h1D_fitDiffPosXFirst] = new TH1F((GetCalorimeter()->GetName() + "_fitDiffPosXFirst").c_str(),
                        "difference between start values and fit result for X-position of first shower",
                        200, -10., 10.);
                fHist[h1D_fitDiffPosYFirst] = new TH1F((GetCalorimeter()->GetName() + "_fitDiffPosYFirst").c_str(),
                        "difference between start values and fit result for Y-position of first shower",
                        200, -10., 10.);
                fHist[h1D_fitDiffTimeFirst] = new TH1F((GetCalorimeter()->GetName() + "_fitDiffTimeFirst").c_str(),
                        "difference between start values and fit result for time of first shower",
                        200, -10., 10.);

                for (std::map<CellType::HwType, size_t>::const_iterator it=fMapCellType2HistIndex.begin(); it!=fMapCellType2HistIndex.end(); it++) {
                    fHist[h2D_fitPullsCellTypeFirst+it->second] = new TH2F((GetCalorimeter()->GetName() + "_fitPulls_" + CellType::GetHwTypeString(it->first)).c_str(),
                            ("energy dependent pulls for cell type " + CellType::GetHwTypeString(it->first)).c_str(),
                            200, 0., 200., 200, -1., 1.);
                    fHist[h2D_fitPullsCellTypeFirst+it->second]->GetXaxis()->SetTitle("measured energy");
                    fHist[h2D_fitPullsCellTypeFirst+it->second]->GetYaxis()->SetTitle("(fitted - measured energy) / measured energy");
                    fHist[h2D_fitPullsCellTypeFirst+it->second]->SetOption("COLZ");
                    fHist[h2D_fitPullsLowCellTypeFirst+it->second] = new TH2F((GetCalorimeter()->GetName() + "_fitPullsLowEnergy_" + CellType::GetHwTypeString(it->first)).c_str(),
                            ("energy dependent pulls for cell type " + CellType::GetHwTypeString(it->first)).c_str(),
                            200, 0., 10., 200, -1., 1.);
                    fHist[h2D_fitPullsLowCellTypeFirst+it->second]->GetXaxis()->SetTitle("measured energy");
                    fHist[h2D_fitPullsLowCellTypeFirst+it->second]->GetYaxis()->SetTitle("(fitted - measured energy) / measured energy");
                    fHist[h2D_fitPullsLowCellTypeFirst+it->second]->SetOption("COLZ");
                    fHist[h2D_fitNormPullsCellTypeFirst+it->second] = new TH2F((GetCalorimeter()->GetName() + "_fitNormPulls_" + CellType::GetHwTypeString(it->first)).c_str(),
                            ("energy dependent pulls for cell type " + CellType::GetHwTypeString(it->first)).c_str(),
                            200, 0., 200., 200, -5., 5.);
                    fHist[h2D_fitNormPullsCellTypeFirst+it->second]->GetXaxis()->SetTitle("measured energy");
                    fHist[h2D_fitNormPullsCellTypeFirst+it->second]->GetYaxis()->SetTitle("(fitted - measured energy) / measurement error");
                    fHist[h2D_fitNormPullsCellTypeFirst+it->second]->SetOption("COLZ");
                    fHist[h2D_fitNormPullsLowCellTypeFirst+it->second] = new TH2F((GetCalorimeter()->GetName() + "_fitNormPullsLowEnergy_" + CellType::GetHwTypeString(it->first)).c_str(),
                            ("energy dependent pulls for cell type " + CellType::GetHwTypeString(it->first)).c_str(),
                            200, 0., 10., 200, -5., 5.);
                    fHist[h2D_fitNormPullsLowCellTypeFirst+it->second]->GetXaxis()->SetTitle("measured energy");
                    fHist[h2D_fitNormPullsLowCellTypeFirst+it->second]->GetYaxis()->SetTitle("(fitted - measured energy) / measurement error");
                    fHist[h2D_fitNormPullsLowCellTypeFirst+it->second]->SetOption("COLZ");
                }
            }
        }
    }
}

/////////////////////////////////////////////////////////////////////////////

std::vector<Cluster> ReconstructionCombined::DoClustering(const std::vector<CellDataRaw>& signals) {
    // vector of clusters (return value)
    std::vector<Cluster> clusters;

    // a map for the cell number to the index in the signal vector
    std::map<size_t, size_t> mapCellSignal;
    for (size_t i=0; i<signals.size(); i++) {
        if (mapCellSignal.count(signals[i].GetCellIdx())>0)
            throw Exception("ReconstructionCombined::DoClustering: Multiple raw data information for the same cell.");

        mapCellSignal.insert(std::pair<size_t, size_t>(signals[i].GetCellIdx(), i));

        if (GetCalorimeter()->GetOptions().fill_histos) {
            assert(fHist[h1D_cellEnergy]!=0);
            fHist[h1D_cellEnergy]->Fill(signals[i].GetEnergy());
            assert(fHist[h1D_cellTime]);
            fHist[h1D_cellTime]->Fill(signals[i].GetTime());
            assert(fHist[h2D_cellEnergyVsTime]!=0);
            fHist[h2D_cellEnergyVsTime]->Fill(signals[i].GetEnergy(), signals[i].GetTime());

            if (GetCalorimeter()->GetOptions().histos_level>=Calorimeter::Options::HISTLVL_VERBOSE) {
                assert(fHist[h1D_cellLowEnergy]!=0);
                fHist[h1D_cellLowEnergy]->Fill(signals[i].GetEnergy());
                assert(fHist[h2D_cellLowEnergyVsTime]!=0);
                fHist[h2D_cellLowEnergyVsTime]->Fill(signals[i].GetEnergy(), signals[i].GetTime());

                if (GetCalorimeter()->GetOptions().histos_level>=Calorimeter::Options::HISTLVL_DEBUG) {
                    assert(fHist[h1D_cellEnergyThreshold]!=0);
                    if (fHist[h1D_cellEnergyThreshold]->GetBinContent(fHist[h1D_cellEnergyThreshold]->FindBin(signals[i].GetCellIdx())) == 0.
                            || fHist[h1D_cellEnergyThreshold]->GetBinContent(fHist[h1D_cellEnergyThreshold]->FindBin(signals[i].GetCellIdx())) > signals[i].GetEnergy())
                        fHist[h1D_cellEnergyThreshold]->SetBinContent(fHist[h1D_cellEnergyThreshold]->FindBin(signals[i].GetCellIdx()), signals[i].GetEnergy());
                    assert(fHist[h1D_cellEnergyThresholdAdc]!=0);
                    if (fHist[h1D_cellEnergyThresholdAdc]->GetBinContent(fHist[h1D_cellEnergyThresholdAdc]->FindBin(signals[i].GetCellIdx())) == 0.
                            || fHist[h1D_cellEnergyThresholdAdc]->GetBinContent(fHist[h1D_cellEnergyThresholdAdc]->FindBin(signals[i].GetCellIdx())) > signals[i].GetAmplitude())
                        fHist[h1D_cellEnergyThresholdAdc]->SetBinContent(fHist[h1D_cellEnergyThresholdAdc]->FindBin(signals[i].GetCellIdx()), signals[i].GetAmplitude());
                }
            }
        }
    }

    std::set<size_t> usedCells; // keep track of the already used cells
    for (size_t i=0; i<signals.size(); i++) {
        // check that this cell was not yet put to any cluster
        if (usedCells.count(signals[i].GetCellIdx())>0)
            continue;

        // TODO: filter dead cells

        // so this is the first cell of a new cluster
        Cluster cl(GetCalorimeter());

        std::set<size_t> checkedCells;
        std::vector<size_t> cellsToCheck;
        cellsToCheck.push_back(signals[i].GetCellIdx());

        while (cellsToCheck.size()>0) {
            size_t curIdx = cellsToCheck.back();
            cellsToCheck.pop_back();

            if (checkedCells.count(curIdx)>0)
                continue;
            checkedCells.insert(curIdx);

            // if a neighboring cell is bad, still search for signals in its neighbors
            if (GetCalorimeter()->CellIsBad(curIdx, Calorimeter::OLD)) {
                const std::vector<size_t>& neighbors = GetCalorimeter()->GetCells()[curIdx].GetNeighbors();
                for (std::vector<size_t>::const_iterator it=neighbors.begin(); it!=neighbors.end(); it++)
                    cellsToCheck.push_back(*it);
            }
            // TODO: filter dead cells

            // if there is a signal at the current cell then add this signal to the cluster
            // and also search in the adjacent cells
            std::map<size_t, size_t>::const_iterator curMap;
            if ( (curMap=mapCellSignal.find(curIdx))!=mapCellSignal.end() ) {
                if (signals[curMap->second].GetEnergy()<GetCalorimeter()->GetOptions().cell_energy_threshold)
                    continue;

                cl.AddCell(signals[curMap->second]);
                usedCells.insert(curIdx);

                const std::vector<size_t>& neighbors = GetCalorimeter()->GetCells()[curIdx].GetNeighbors();
                for (std::vector<size_t>::const_iterator it=neighbors.begin(); it!=neighbors.end(); it++)
                    cellsToCheck.push_back(*it);
            }
        }

        // size 0 might happen, if all remaining cells have a signal below threshold
        if (cl.Size()>0)
            clusters.push_back(cl);
    }

    return clusters;
}

/////////////////////////////////////////////////////////////////////////////

double ReconstructionCombined::DoFitting(const std::vector<FitInfo>& cells, std::vector<FitShower>& shower) {
    // object to hold data to fit to
    FitInput fitInput(GetCalorimeter(), cells, nrFitParameters*shower.size());
// DEBUG
//        fitInput.Print();

    // create fitter
    ROOT::Math::Minimizer* minimizer = ROOT::Math::Factory::CreateMinimizer();
    minimizer->SetFunction(fitInput);
    minimizer->SetErrorDef(0.5);
    //minimizer->SetPrintLevel(3);

    // set start parameters of fitter
    for (size_t i=0; i<shower.size(); i++) {
        char desc[15];

        sprintf(desc, "Energy%02zu", i+1);
        minimizer->SetVariable(nrFitParameters*i+0, desc, shower[i].fE, shower[i].fEerr);
        sprintf(desc, "xPosition%02zu", i+1);
        minimizer->SetVariable(nrFitParameters*i+1, desc, shower[i].fX, shower[i].fXerr);
        sprintf(desc, "yPosition%02zu", i+1);
        minimizer->SetVariable(nrFitParameters*i+2, desc, shower[i].fY, shower[i].fYerr);
        sprintf(desc, "Time%02zu", i+1);
        minimizer->SetVariable(nrFitParameters*i+3, desc, shower[i].fT, shower[i].fTerr);
    }

    // do not print internal error messages from the fitter
    Int_t oldIgnoreLevel = gErrorIgnoreLevel;
    gErrorIgnoreLevel = kError+1;

    // do the fitting
    minimizer->Minimize();

    // set back old ignore level
    gErrorIgnoreLevel = oldIgnoreLevel;

    // extract the shower parameters from the fitter
    shower.clear();
    for (unsigned int i=0; i<fitInput.NDim()/nrFitParameters; i++) {
        double val, valErr;

        FitShower s(0.,0.,0.,0.,0.,0.,0.,0.);

        val = minimizer->X()[nrFitParameters*i+0];
        valErr = minimizer->Errors()[nrFitParameters*i+0];
        if (std::isnan(valErr)) valErr=0.1;
        s.fE=val; s.fEerr=valErr;

        val = minimizer->X()[nrFitParameters*i+1];
        valErr = minimizer->Errors()[nrFitParameters*i+1];
        if (std::isnan(valErr)) valErr=0.1;
        s.fX=val; s.fXerr=valErr;

        val = minimizer->X()[nrFitParameters*i+2];
        valErr = minimizer->Errors()[nrFitParameters*i+2];
        if (std::isnan(valErr)) valErr=0.1;
        s.fY=val; s.fYerr=valErr;

        val = minimizer->X()[nrFitParameters*i+3];
        valErr = minimizer->Errors()[nrFitParameters*i+3];
        if (std::isnan(valErr)) valErr=0.1;
        s.fT=val; s.fTerr=valErr;

        // now find the main cell of the shower,
        // the cell with signal closest to the shower center
        double min(std::numeric_limits<double>::max());
        size_t cell(std::numeric_limits<size_t>::max());
        for (std::vector<FitInfo>::const_iterator it=cells.begin();
             it!=cells.end(); it++) {
            // skip dead cells and other cells artificially added (for example
            // the ones along the border)
            if (it->GetArtificial())
                continue;

            const double distX = GetCalorimeter()->GetCells()[it->GetData().GetCellIdx()].GetX() - s.fX;
            const double distY = GetCalorimeter()->GetCells()[it->GetData().GetCellIdx()].GetY() - s.fY;
            const double dist  = sqrt(distX*distX + distY*distY);

            if (dist < min) {
                min  = dist;
                cell = it->GetData().GetCellIdx();
            }
        }
        s.fMainCell = cell;

// DEBUG
//                std::cout << i << ": " << s.fE << " " << s.fEerr << ", " << s.fX << " " << s.fXerr << ", " << s.fY << " " << s.fYerr << std::endl;

        shower.push_back(s);
    }

    // clean up the fitter
    delete minimizer;

    return fitInput.GetLogLh(shower);
}

/////////////////////////////////////////////////////////////////////////////

ReconstructionCombined::FitShower ReconstructionCombined::ApplyCorrections(const FitShower& shower) {
    FitShower out(shower);

    const size_t    cellIdx (shower.fMainCell);
    const Cell&     cell    (GetCalorimeter()->GetCells()[cellIdx]);
    const CellType& cellType(cell.GetCellType());

    // calculate distance of shower from main cell center
    const double distX(shower.fX - cell.GetX());
    const double distY(shower.fY - cell.GetY());

    if (fPosCorrX != 0)
        out.fX -= fPosCorrX->GetCorrection(shower.fE, distX);

    if (fPosCorrY != 0)
        out.fY -= fPosCorrY->GetCorrection(shower.fE, distY);

    std::map<CellType::HwType, CorrectionPosDepEnergy*>::const_iterator it2;
    if ( (it2=fPosDepEnergyCorr.find(cellType.GetHwType()))!=fPosDepEnergyCorr.end() )
        out.fE = it2->second->GetCorrection(shower.fE, distX, distY);

    return out;
}

/////////////////////////////////////////////////////////////////////////////

double ReconstructionCombined::CalcEnergyError(const double& energy) {
    assert(fParamEnergyErr.size()==3);

    return sqrt(  (fParamEnergyErr[0]*fParamEnergyErr[0]) * (energy)
                + (fParamEnergyErr[1]*fParamEnergyErr[1]) * (energy*energy)
                + (fParamEnergyErr[2]*fParamEnergyErr[2]) );
}

/////////////////////////////////////////////////////////////////////////////

double ReconstructionCombined::CalcTimeError(const double& energy) {
    assert(fParamTimeErr.size()==3);

    return sqrt(  (fParamTimeErr[0]*fParamTimeErr[0]) / (energy)
                + (fParamTimeErr[1]*fParamTimeErr[1]) / (energy*energy)
                + (fParamTimeErr[2]*fParamTimeErr[2]) );
}

/////////////////////////////////////////////////////////////////////////////

ReconstructionCombined::FitInfo::FitInfo(const Calorimeter* c, const CellDataRaw& s) :
    fCalorimeter(c), fData(s), fArtificial(false), fAround(false) {
}

/////////////////////////////////////////////////////////////////////////////

const CellType& ReconstructionCombined::FitInfo::GetCellType() const {
    return fCalorimeter->GetCells()[GetData().GetCellIdx()].GetCellType();
}

/////////////////////////////////////////////////////////////////////////////

void ReconstructionCombined::FitInfo::SetEnergy(double e) {
    if (!GetArtificial())
        std::cerr << "ReconstructionCombined::FitInfo::SetEnergy (" << __FILE__ << ":" << __LINE__ << "): re-setting the energy for a cell with measured energy!" << std::endl;

    fData.SetEnergy(e);
}

/////////////////////////////////////////////////////////////////////////////

void ReconstructionCombined::FitInfo::SetArtificial(const bool& artificial) {
    fArtificial = artificial;
}

/////////////////////////////////////////////////////////////////////////////

bool ReconstructionCombined::FitInfo::GetArtificial() const {
    return fArtificial;
}

/////////////////////////////////////////////////////////////////////////////

void ReconstructionCombined::FitInfo::SetAround(const bool& around) {
    fAround = around;
}

/////////////////////////////////////////////////////////////////////////////

bool ReconstructionCombined::FitInfo::GetAround() const {
    return fAround;
}

/////////////////////////////////////////////////////////////////////////////

double ReconstructionCombined::FitInfo::GetFittedEnergy(const std::vector<FitShower>& shower) const {
    double totE(0.);
    for (std::vector<FitShower>::const_iterator it=shower.begin(); it!=shower.end(); it++) {
        double ecalc=it->fE*fCalorimeter->GetCells()[fData.GetCellIdx()].GetCellType().GetShowerProfile()->GetEnergyInCell(fCalorimeter->GetCells()[fData.GetCellIdx()].GetX()-it->fX,
                                                                                                                           fCalorimeter->GetCells()[fData.GetCellIdx()].GetY()-it->fY);
        totE+=ecalc;
    }

    return totE;
}

/////////////////////////////////////////////////////////////////////////////

} // namespace Reco

#endif // ROOT version >= 5.22

