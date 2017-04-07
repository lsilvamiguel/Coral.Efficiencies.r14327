// test the pi0 calibration

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

#include "ToyCalorimeter.h"

#include <Calorimeter.h>
#include <Exception.h>
#include <ReconstructionCombined.h>

typedef Reco::CellType                                                 CellType;
typedef Reco::Exception                                                Exception;
typedef Reco::ReconstructionCombined::CorrectionPosDepEnergy           CorrectionPosDepEnergy;
typedef Reco::ReconstructionCombined::CorrectionPosDepEnergyBinQuarter CorrectionPosDepEnergyBinQuarter;

void usage(const char* exe) {
    std::cout << "expecting the file with the parameters for the corrections, and an output histogram file as arguments:" << std::endl;
    std::cout << "        " << exe << " [file, e.g. test/testEdepCorr.dat] [histogram file, e.g. test/testEdepCorr.root]" << std::endl;
}

int main(int argc, char* argv[]) {
    if (argc!=3) {
        usage(argv[0]);
        return 1;
    }
    
    std::ifstream calib(argv[1]);
    if (!calib.good()) {
        usage(argv[0]);
        return 1;
    }

    std::ostringstream calibS;
    calibS << calib.rdbuf();    

    calib.close();

    ToyCalorimeter* calo = new ToyCalorimeter;
    calo->InputEdepCorr(calibS.str());

    const unsigned int cellsCount = 12;
    const unsigned int cellsId[cellsCount] = {  0,  1,  2,  3,  4,
                                               15, 16, 17, 18, 19, 20,
                                               30 };

    TGraph* g[cellsCount];
    for (unsigned int cell=0; cell<cellsCount; ++cell)
        g[cell] = new TGraph;

    const double minE = 0.;
    const double maxE = 3.;
    const unsigned int steps = 100;
    for (unsigned int i=1; i<steps; ++i) {
        const double energy = minE + i*(maxE-minE)/steps;

        for (unsigned int cell=0; cell<cellsCount; ++cell)
            g[cell]->SetPoint(g[cell]->GetN(), energy, calo->SignalToEnergy(cellsId[cell], energy, 0.) / energy);
    }

    TFile* file = TFile::Open(argv[2], "RECREATE");
    for (unsigned int cell=0; cell<cellsCount; ++cell) {
        std::ostringstream nameS;
        nameS << "g" << cell;
        g[cell]->Write(nameS.str().c_str());
    }
    file->Close();

    return 0;
}

