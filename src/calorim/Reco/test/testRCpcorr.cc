// test the position dependent energy correction in the Shashlyks used in
// ReconstructionCombined

#include <algorithm>
#include <cassert>
#include <fstream>
#include <iostream>
#include <sstream>

#include <TFile.h>
#include <TGraph.h>
#include <TH2D.h>
#include <TRandom.h>

#include <Exception.h>
#include <ReconstructionCombined.h>

typedef Reco::Exception                                             Exception;
typedef Reco::ReconstructionCombined::CorrectionPositionBinSin      CorrectionPositionBinSin;
typedef Reco::ReconstructionCombined::CorrectionPositionBinPol3     CorrectionPositionBinPol3;
typedef Reco::ReconstructionCombined::CorrectionPosition            CorrectionPosition;
typedef Reco::ReconstructionCombined::CorrectionPositionSinAtan     CorrectionPositionSinAtan;
typedef Reco::ReconstructionCombined::CorrectionPositionSinSin3Atan CorrectionPositionSinSin3Atan;

void usage(const char* exe) {
    std::cout << "expecting the file with the parameters for the corrections, and an output histogram file as arguments:" << std::endl;
    std::cout << "        " << exe << " [file, e.g. test/testRCecorr.dat] [histogram file, e.g. test/testRCpcorr.root]" << std::endl;
}

int main(int argc, char* argv[]) {
    if (argc!=3) {
        usage(argv[0]);
        return 1;
    }

    CorrectionPosition* fPosCorrX(0);
    CorrectionPosition* fPosCorrY(0);

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
        if (token == "CORRECTIONPOSITIONX") {
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
        }
    }

    calib.close();

    fPosCorrX->Print(std::cout, "    position corrections in X                     : ");
    fPosCorrY->Print(std::cout, "    position corrections in Y                     : ");
    std::cout << std::endl;

    TFile* resultFile=new TFile(argv[2], "RECREATE");

    TH2D* hPosCorrX = new TH2D("hPosCorrX", "Position correction along X", 200, 0., 200., 200, -20., 20.);
    for (int eb=1; eb<=hPosCorrX->GetNbinsX(); eb++) {
        const double energy = hPosCorrX->GetXaxis()->GetBinCenter(eb);

        for (int db=1; db<=hPosCorrX->GetNbinsY(); db++) {
            const double dist = hPosCorrX->GetYaxis()->GetBinCenter(db);
            const double corr = fPosCorrX->GetCorrection(energy, dist);

            hPosCorrX->Fill(energy, dist, corr);
        }
    }
    hPosCorrX->Write();
    delete hPosCorrX;

    TH2D* hPosCorrY = new TH2D("hPosCorrY", "Position correction along Y", 200, 0., 200., 200, -20., 20.);
    for (int eb=1; eb<=hPosCorrY->GetNbinsX(); eb++) {
        const double energy = hPosCorrY->GetXaxis()->GetBinCenter(eb);

        for (int db=1; db<=hPosCorrY->GetNbinsY(); db++) {
            const double dist = hPosCorrY->GetYaxis()->GetBinCenter(db);
            const double corr = fPosCorrY->GetCorrection(energy, dist);

            hPosCorrY->Fill(energy, dist, corr);
        }
    }
    hPosCorrY->Write();
    delete hPosCorrY;

    std::cout << "finishing..." << std::endl;

    resultFile->Write();
    resultFile->Close();
    delete resultFile;

    delete fPosCorrX;
    delete fPosCorrY;

    return 0;
}

