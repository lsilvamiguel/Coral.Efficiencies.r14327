// test the position dependent energy correction in the Shashlyks used in
// ReconstructionCombined

#include <algorithm>
#include <cassert>
#include <fstream>
#include <iostream>
#include <list>
#include <sstream>

#include <TFile.h>
#include <TGraph.h>
#include <TH1D.h>
#include <TRandom.h>

#include <Calorimeter.h>
#include <Exception.h>
#include <ReconstructionCombined.h>

typedef Reco::CellType                                                 CellType;
typedef Reco::Exception                                                Exception;
typedef Reco::ReconstructionCombined::CorrectionPosDepEnergy           CorrectionPosDepEnergy;
typedef Reco::ReconstructionCombined::CorrectionPosDepEnergyBinQuarter CorrectionPosDepEnergyBinQuarter;

Reco::Calorimeter* calo(0);
Reco::Calorimeter* GetCalorimeter() {
    if (calo==0) {
        calo = new Reco::Calorimeter("TOYCALO", "");

        calo->InsertCellType(CellType("TC001",
                                      38.2, 38.2, 395.5,
                                      38.3, 38.3,
                                      17.65, 381.5,
                                      8.6,
                                      0.065, 0.02, 0.01,
                                      0.1471, 0.,
                                      "SHASHLIK"));
    }

    return calo;
}

void usage(const char* exe) {
    std::cout << "expecting the file with the parameters for the corrections, and an output histogram file as arguments:" << std::endl;
    std::cout << "        " << exe << " [file, e.g. test/testRCecorr.dat] [histogram file, e.g. test/testRCecorr.root]" << std::endl;
}

int main(int argc, char* argv[]) {
    if (argc!=3) {
        usage(argv[0]);
        return 1;
    }
    
    const double stepSize(38.3);
    std::vector<double> borders;
    std::map<Reco::CellType::HwType, Reco::ReconstructionCombined::CorrectionPosDepEnergy*> fPosDepEnergyCorr;

    std::ifstream calib(argv[1]);
    if (!calib.good()) {
        usage(argv[0]);
        return 1;
    }

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
        if (token == "CORRECTIONENERGYPOSDEP") {
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
                    borders.push_back(minEnergy);
                }
            }
        }
    }

    calib.close();

    std::cout << "    position dependent energy correction          : ";
    if (fPosDepEnergyCorr.size()==0)
        std::cout << "[none]" << std::endl;
    else
        std::cout << std::endl;
        for (std::map<Reco::CellType::HwType, Reco::ReconstructionCombined::CorrectionPosDepEnergy*>::const_iterator it=fPosDepEnergyCorr.begin(); it!=fPosDepEnergyCorr.end(); it++) {
            std::string title;
            for (size_t i=0; i<44-Reco::CellType::GetHwTypeString(it->first).length(); i++)
                title += " ";
            title += "* for " + Reco::CellType::GetHwTypeString(it->first) + ": ";
            it->second->Print(std::cout, title);
        }
    std::cout << std::endl;

    TFile* resultFile=new TFile(argv[2], "RECREATE");

    for (std::map<Reco::CellType::HwType, Reco::ReconstructionCombined::CorrectionPosDepEnergy*>::const_iterator it=fPosDepEnergyCorr.begin(); it!=fPosDepEnergyCorr.end(); it++) {
        std::cout << "starting tests of " << Reco::CellType::GetHwTypeString(it->first) << ":" << std::endl;

        TDirectory* thisDir = resultFile->mkdir(Reco::CellType::GetHwTypeString(it->first).c_str());
        thisDir->cd();

        TGraph* g=new TGraph;

        for (unsigned int i=0; i<10000; i++) {
            const double x(0.);
            const double y(0.);
            const double e(200.*i/10000.);
            const double corr = it->second->GetCorrection(e, x, y);

            g->SetPoint(i, e, corr);
        }

        g->Write("graph");
        delete g;

        unsigned int count(0);
        for (std::vector<double>::const_iterator it2=borders.begin(); it2!=borders.end(); it2++) {
            TH1D* hFinalEnergy=new TH1D(TString::Format("hFinalEnergy%02u", count++), "energy after correction", 1000, *it2-10., *it2+10.);

            for (unsigned int i=0; i<10000000; i++) {
                const double x(0.);
                const double y(0.);
                const double e(gRandom->Uniform(std::max(0., *it2-20.), *it2+20.));
                const double corr = it->second->GetCorrection(e, x, y);

                hFinalEnergy->Fill(e*corr);
            }

            hFinalEnergy->Write();
            delete hFinalEnergy;
        }

        std::cout << std::endl;
    }

    std::cout << "finishing..." << std::endl;

    resultFile->Write();
    resultFile->Close();
    delete resultFile;

    for (std::map<Reco::CellType::HwType, Reco::ReconstructionCombined::CorrectionPosDepEnergy*>::const_iterator it=fPosDepEnergyCorr.begin(); it!=fPosDepEnergyCorr.end(); it++) 
        delete it->second;

    return 0;
}

