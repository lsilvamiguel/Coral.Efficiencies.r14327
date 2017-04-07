#include "ToyCalorimeter.h"
#include <ShowerProfileLednev.h>

#include <sstream>

ToyCalorimeter::ToyCalorimeter(double cellThreshold) : Reco::Calorimeter("toyCalorimeter") {
    // build a calorimeter just out of a simple 15x15 matrix of cells
    // parameters of the cells are similar to GAMS cells
    // matrix is centered on (0,0)
    InsertCellType(Reco::CellType("toyCellType", 37.3, 37.3, 45.,
                                                 38.3, 38.3,
                                                 2.59, 42., 0.014725,
                                                 0.065, 0.02, 0.01,
                                                 0.8686, 0.,
                                                 "GAMS"));
    std::vector<double> params; params.resize(4);
    params[0] = 1.191; params[1] = -0.191; params[2] = 6.664; params[3] = 43.777;
    cells_type.back().SetShowerProfile(new Reco::ShowerProfileLednev(this, cells_type.back(), params));

    for (int y=-7; y<=7; y++)
             for (int x=-7; x<=7; x++) {
                 const double pos_x = x*cells_type.back().GetStepX();
                 const double pos_y = y*cells_type.back().GetStepY();
                 const double pos_z = 0;

                 cells.push_back( Reco::Cell(cells_type.back(), true, pos_x, pos_y, pos_z) );
             }

    InitMemory();
    InitOptions();
    options.SetRecoOptions("RECONSTRUCTION=COMBINED");
    options.SetRecoOptions("HIST=YES");

    std::ostringstream strCellThreshold;
    strCellThreshold << "CELL_ENERGY_THRESHOLD=" << cellThreshold;
    options.SetRecoOptions(strCellThreshold.str());

    options.SetRecoOptions("PARTICLE_ENERGY_THRESHOLD=0.2");
    options.SetRecoOptions("SCALE_CALIBRATION=1.0");
    options.SetRecoOptions("USE_PROB_NOISE=YES");
    options.SetRecoOptions("PROB_CELL_NOISE=0.005");
    options.SetRecoOptions("MIN_SHOWER_DISTANCE=0.7");
    options.SetRecoOptions("USE_EDEPCORR=YES");
    options.SetRecoOptions("USE_TISDEPCORR=YES");
    options.SetOK();
    Init();

    InitXY();
}

ToyCalorimeter::~ToyCalorimeter() {
}

