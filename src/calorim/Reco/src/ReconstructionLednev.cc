/*
   $Source: /afs/cern.ch/user/g/gobbo/w0/coral/src/calorim/Reco/src/ReconstructionLednev.cc,v $
   $Date: 2010/12/05 22:18:58 $
   $Revision: 1.10 $
*/

// --- Standard C/C++ library ---
#include <cassert>
#include <cmath>
#include <iostream>

// --- Internal files ---
#include "ReconstructionLednev.h"

#include "Calorimeter.h"
#include "CellDataRaw.h"
#include "CellType.h"
#include "Shower.h"

//////////////////////////////////////////////////////////////////////////////////
//
// FORTRAN STUFF
//
extern "C" void separa_();

extern struct runpar_cb {
  int nrun_;
  int numbev_;
  int nburst_;
  int numbinburst_;
  int ntrig_;
} runpar_;

extern struct cells_cb {
  int   ncell_;
  float ecell_[3072];
  int   iaddr_[3072];
  float xcell_[3072];
  float ycell_[3072];
  float errcell_[3072];
  float time_[3072];
} cells_;

extern struct gamma_cb {
  int   ngam_;
  float egam_[200];
  float xgam_[200];
  float ygam_[200];
  float chigam_[200];
  int   ndf_[200];
  float gtime_[200];
} gamma_;

extern struct misc_cb {
  float Eerr_[200];
  float Xerr_[200];
  float Yerr_[200];
  int   lnclust_[200];
  int   indmc_[200];
} misc_;

extern struct cluster_cb {
  int ncl_;
  int lencl_[500];
  int listcells_[3072];
  int idx_;
} cluster_;

extern struct iflag_cb {
  int ibad_;
  int icalib_;
  int ifast_;
  int ihist_;
  int ivarthr_;
} iflags_;

extern struct threshold_cb {
   float thrcommon_;
   float cellthr_[3072];
} thresholds_;

extern struct ccut_cb {
  float egcut_;
  float time0_;
  float timew_;
} ccuts_;

extern struct showerprof_cb {
  float a_[4];
  float b_[4];
} showerprof_;

extern struct rawcells_cb {
  float acell_[3072];
} rawcells_;

namespace Reco {

/////////////////////////////////////////////////////////////////////////////

ReconstructionLednev::ReconstructionLednev(const Calorimeter* c) : Reconstruction(c) {
    if (c->GetName() != "EC02P1__") {
        std::cerr << " Sorry! But FortranALReconstruction is valid only for EC02P1__  but not for " << GetCalorimeter()->GetName() << std::endl;
        std::cerr << " Change options for " << GetCalorimeter()->GetName() <<"  please " << std::endl;
        exit(1);
    }
}

/////////////////////////////////////////////////////////////////////////////

const std::vector<CalorimeterParticle>& ReconstructionLednev::DoReconstruction(const std::vector<CellDataRaw>& signals) {
//    const bool debug = true;
    const bool debug = false;

    if (debug)
        std::cout << " This is Calorimeter  shower separation for " << GetCalorimeter()->GetName() << std::endl;

    runpar_.nrun_        = GetCalorimeter()->GetEventID().GetRunNumber();
    runpar_.numbev_      = GetCalorimeter()->GetEventID().GetEventNumberInRun();
    runpar_.nburst_      = GetCalorimeter()->GetEventID().GetBurstNumber();
    runpar_.numbinburst_ = GetCalorimeter()->GetEventID().GetEventNumberInBurst();
    runpar_.ntrig_       = GetCalorimeter()->GetEventID().GetTriggerMask();
    double adcThreshold  = GetCalorimeter()->GetOptions().cell_adc_threshold;
//    std::cout << "adcThreshold ="<<adcThreshold<<std::endl;

    double dcell = 38.3; int ix; int iy;
    int ncell = 0;
    for (std::vector<CellDataRaw>::const_iterator it=signals.begin(); it!=signals.end(); it++) {
        size_t icell = it->GetCellIdx();
        if ( icell>=GetCalorimeter()->NCells() ) {
            std::cerr << "ReconstructionLednev::DoReconstruction: Bad cell address " << icell << " is out of calorimeter size " << GetCalorimeter()->NCells() << "." << std::endl;
            exit(1);
            // TODO: should be fatal exception, but is caught somewhere
            throw Exception("ReconstructionLednev::DoReconstruction: Bad cell address %zu is out of calorimeter size %zu.",
                            icell, GetCalorimeter()->NCells());
        }

        rawcells_.acell_[ncell] = it->GetAmplitude();
        cells_.ecell_[ncell] = it->GetEnergy();
        cells_.xcell_[ncell] = GetCalorimeter()->GetCells()[icell].GetX();
        cells_.ycell_[ncell] = GetCalorimeter()->GetCells()[icell].GetY();
        ix = int((cells_.xcell_[ncell]+32.*dcell)/dcell);
        iy = int((cells_.ycell_[ncell]+24.*dcell)/dcell);
//        std::cout << "ix,xcell ="<<ix<<" "<<cells_.xcell_[ncell]<<std::endl;
        if(ix<0 || ix>63) std::cout << "ix,xcell ="<<ix<<" "<<cells_.xcell_[ncell]<<std::endl;
        if(iy<0 || iy>47) std::cout << "iy,ycell ="<<iy<<" "<<cells_.ycell_[ncell]<<std::endl;
        thresholds_.cellthr_[ix*48+iy]=adcThreshold*cells_.ecell_[ncell]/rawcells_.acell_[ncell];

//        std::cout << "ix,iy,thr ="<<ix<<" "<<iy<<" "<<thresholds_.cellthr_[ix*48+iy]<<std::endl;

        cells_.time_[ncell] = it->GetTime();
        if( !GetCalorimeter()->GetOptions().use_time_in_reconstruction )
            cells_.time_[ncell] = 0.;

        if ( debug )
            std::cout << " Set cell " << icell << " " << GetCalorimeter()->GetCellName(icell).c_str()
                      << " e= " << it->GetEnergy() <<" time = " <<cells_.time_[ncell] << std::endl;
        ncell++;

        if (fHist1D[GetCalorimeter()->GetName() + "_cellEnergy"])   fHist1D[GetCalorimeter()->GetName() + "_cellEnergy"]->Fill(it->GetEnergy());
        if (fHist2D[GetCalorimeter()->GetName() + "_cellPosition"]) fHist2D[GetCalorimeter()->GetName() + "_cellPosition"]->Fill(GetCalorimeter()->GetColumnOfCell(icell), GetCalorimeter()->GetRowOfCell(icell));
    }
//    std::cout <<"ncell = "<<ncell<<std::endl;
    cells_.ncell_ = ncell;

    separa_();

    reconstructedParticles.clear();
//    std::cout << "Output of separa: NumbinBurst = "<<runpar_.numbinburst_<<" Ngamma = "<< gamma_.ngam_<< std::endl;
    for (int i = 0; i < gamma_.ngam_; i++) {
        const double e    = gamma_.egam_[i];
        const double x    = gamma_.xgam_[i];
        const double y    = gamma_.ygam_[i];

        if (e <= 0)
            continue;

        const int iix = misc_.indmc_[i] / 48;
        const int iiy = misc_.indmc_[i] - 48*iix;
        int m   = GetCalorimeter()->GetCellOfColumnRow( iix, iiy);
        if( m < 0 ) {
            std::cerr << " Calorimeter error : cell not found at x="<< x << ", y=" << y
                      << " (main cell of shower suggests it is column=" << iix << ", row=" << iiy << ")." << std::endl;
            exit(1);
            // TODO: should be fatal exception, but is caught somewhere
            throw Exception("ReconstructionLednev::DoReconstruction: main cell of shower at x=%f, y=%f not found, should be at column=%d, row=%d.",
                            x, y, iix, iiy);
        }
//        std::cout <<" igam "<<i<<" E="<<e<<" X="<<x<<" Y="<<y<<" ix="<<iix<<" iy="<<iiy<<" m="<<m<<" Main cell="<<misc_.indmc_[i]<< std::endl;

        const double RadLeng = GetCalorimeter()->GetCells()[m].GetCellType().GetRadiationLength();
        const double SizeZ   = GetCalorimeter()->GetCells()[m].GetCellType().GetSizeZ() / 2.;
        const double z0      = GetCalorimeter()->GetCells()[m].GetZ() - SizeZ;
        const double zw      = ZmidShowerEM(e)*RadLeng;

        const double chi  = gamma_.chigam_[i];
        const double Eerr = misc_.Eerr_[i];
        const double Xerr = misc_.Xerr_[i];
        const double Yerr = misc_.Yerr_[i];
        const double Zerr = 50.;

        reconstructedParticles.push_back (CalorimeterParticle (GetCalorimeter()->GetOptions().particle_default, 1,
                                                               e,    x,    y,    zw+z0,
                                                               GetCalorimeter(),
                                                               Eerr, Xerr, Yerr, Zerr, 0, 0 ));
        reconstructedParticles.back().SetHitedCell((unsigned)m);
        reconstructedParticles.back().SetProb(chi);
        reconstructedParticles.back().SetMiscInfo( CalorimeterParticle::CLUSTER_SIZE, (double)misc_.lnclust_[i]);
        reconstructedParticles.back().SetMiscInfo( CalorimeterParticle::CHI_GAM, (double)gamma_.chigam_[i]);
        reconstructedParticles.back().SetMiscInfo( CalorimeterParticle::NDF_GAM, (double)gamma_.ndf_[i]);

        std::vector<std::pair<size_t, double> > clusterData;

        // fill the cluster data
        // first find the cluster which contains the main cell
        int clLength(-1);
        int clStart(0), clTmpStart(0);
        for (int j=0; j<cluster_.ncl_; j++) {
            for (int k=clTmpStart; k<clTmpStart+cluster_.lencl_[j]; k++) {
                if (cells_.iaddr_[cluster_.listcells_[k]-1]==misc_.indmc_[i]) {
                    clLength=cluster_.lencl_[j];
                    clStart =clTmpStart;
                }
            }

            clTmpStart += cluster_.lencl_[j];
        }

        if (clLength==-1) {
            std::cerr << "ReconstructionLednev::DoReconstruction: Did not find main cell " << misc_.indmc_[i]
                      << " for shower at x="<< x << ", y=" << y << "." << std::endl;
            exit(1);
            // TODO: should be fatal exception, but is caught somewhere
            throw Exception("ReconstructionLednev::DoReconstruction: Did not find main cell %d associated to cluster in %s.",
                            misc_.indmc_[i], GetCalorimeter()->GetName().c_str());
        }

        // then loop over this cluster and get cell indices and cell energies
        for (int k=clStart; k<clStart+clLength; k++) {
            int cellX   = cells_.iaddr_[cluster_.listcells_[k]-1] / 48;
            int cellY   = cells_.iaddr_[cluster_.listcells_[k]-1] - 48*cellX;

            double energy = cells_.ecell_[cluster_.listcells_[k]-1];
            // this was probably a cell added during the reconstruction in either
            // the hole or a dead cell
            if (energy==0.)
                continue;

            int cellIdx = GetCalorimeter()->GetCellOfColumnRow(cellX, cellY);
//            std::cout << "Energy "<<energy<<" cellX "<<cellX<<" cellY "<<cellY<<" cellIdx "<<cellIdx <<std::endl;

            if(cellIdx < 0 ) {
                std::cerr << "ReconstructionLednev::DoReconstruction: Calorimeter cell not found with"
                          << " column=" << cellX
                          << ", row=" << cellY << "." << std::endl;
                // TODO: should be fatal exception, but is caught somewhere
                exit(1);
                throw Exception("ReconstructionLednev::DoReconstruction: Calorimeter cell not found with column=%dn, row=%d.",
                                cellX, cellY);
            }

            size_t icell(cellIdx);

            clusterData.push_back(std::pair<size_t, double>(icell, energy));
        }

        reconstructedParticles.back().SetClusterData(clusterData);

        if (fHist1D[GetCalorimeter()->GetName() + "_showerEnergy"])   fHist1D[GetCalorimeter()->GetName() + "_showerEnergy"]->Fill(e);
        if (fHist2D[GetCalorimeter()->GetName() + "_showerPosition"]) fHist2D[GetCalorimeter()->GetName() + "_showerPosition"]->Fill(iix, iiy);
    }

    if ( debug ) {
        std::cout << "ReconstructionLednev::DoReconstruction: Ngamma = " << reconstructedParticles.size() << std::endl;
        for( unsigned i=0; i < reconstructedParticles.size(); i++) {
            double e = reconstructedParticles[i].GetE();
            double x = reconstructedParticles[i].GetX();
            double y = reconstructedParticles[i].GetY();
            double z = reconstructedParticles[i].GetZ();
            int    m = reconstructedParticles[i].GetMainCells()[0];
            std::cout << " Particle " << i << " E="<< e << " X=" << x << " Y=" << y << " Z=" << z << " Main cell=" << m << std::endl;
        }
    }

    if ( debug )
        std::cout << " shower separation for " << GetCalorimeter()->GetName() << " Ok " << std::endl;

    return reconstructedParticles;
}

////////////////////////////////////////////////////////////////////////////////

const std::vector<CalorimeterParticle>& ReconstructionLednev::DoRepeatReconstruction(const std::vector<CellDataRaw>& signals) {
    if ( GetCalorimeter()->GetOptions().repeat_reconstruction >= 0 )
        throw Exception("ReconstructionLednev::DoRepeatReconstruction: No repetition of reconstruction available for %s.", GetCalorimeter()->GetName().c_str());

    return reconstructedParticles;
}

////////////////////////////////////////////////////////////////////////////////

void ReconstructionLednev::BookHistograms() {
    if (GetCalorimeter()->GetOptions().fill_histos) {
        TDirectory* base     = GetCalorimeter()->GetHistogramsBaseDir();
        TDirectory* recoBase = base->mkdir((GetCalorimeter()->GetName() + "_reconstruction").c_str(),
                                           (std::string("Histograms for reconstruction of ") + GetCalorimeter()->GetName()).c_str());

        recoBase->cd();

        fHist1D[GetCalorimeter()->GetName() + "_cellEnergy"] = new TH1F((GetCalorimeter()->GetName() + "_cellEnergy").c_str(),
                                                                        "Energy deposited in one cell",
                                                                        200, 0., 200.);
        fHist1D[GetCalorimeter()->GetName() + "_showerEnergy"] = new TH1F((GetCalorimeter()->GetName() + "_showerEnergy").c_str(),
                                                                        "Energy in reconstructed shower",
                                                                        200, 0., 200.);

        fHist2D[GetCalorimeter()->GetName() + "_cellPosition"] = new TH2F((GetCalorimeter()->GetName() + "_cellPosition").c_str(),
                                                                          "Position of cells with signal",
                                                                          64, -0.5, 63.5, 48, -0.5, 47.5);
        fHist2D[GetCalorimeter()->GetName() + "_showerPosition"] = new TH2F((GetCalorimeter()->GetName() + "_showerPosition").c_str(),
                                                                             "Position of showers",
                                                                             64, -0.5, 63.5, 48, -0.5, 47.5);

    }
}

////////////////////////////////////////////////////////////////////////////////

void ReconstructionLednev::ReadCalibrations() {
    thresholds_.thrcommon_ = GetCalorimeter()->GetOptions().cell_energy_threshold;

    ccuts_.egcut_ = GetCalorimeter()->GetOptions().particle_energy_threshold;
    ccuts_.time0_ = GetCalorimeter()->GetOptions().ecal2_time_0;
    ccuts_.timew_ = GetCalorimeter()->GetOptions().ecal2_time_width;

    iflags_.ibad_    = GetCalorimeter()->GetOptions().fortran_bad_cell_check;
    iflags_.ifast_   = GetCalorimeter()->GetOptions().fortran_fast_reconstruction;
    iflags_.ihist_   = GetCalorimeter()->GetOptions().fortran_hbook;
    iflags_.ivarthr_ = GetCalorimeter()->GetOptions().fortran_variable_thresholds;


    for (size_t i=0; i<GetCalorimeter()->NCells(); i++) {
        int x = GetCalorimeter()->GetColumnOfCell(i);
        int y = GetCalorimeter()->GetRowOfCell(i);
        int icell = 48*x+y;
        if ( GetCalorimeter()->GetOptions().update_ecut_using_prob_cell_noise ) {
            // the energy cut should only be zero when no calibrations have
            // been read at all, for example in MC case
            if ( GetCalorimeter()->GetEnergyCutInCell(i)==0. )
                thresholds_.cellthr_[icell] = GetCalorimeter()->GetOptions().cell_energy_threshold;
//AL            else
//AL               thresholds_.cellthr_[icell] = GetCalorimeter()->GetEnergyCutInCell(i);
        } //ALelse
//AL            thresholds_.cellthr_[icell] = 0.;
    }

    // set parameters of shower profiles
    if (GetCalorimeter()->GetOptions().fortran_showerprof_a.size() != 0 ||
        GetCalorimeter()->GetOptions().fortran_showerprof_b.size() != 0) {
        if (GetCalorimeter()->GetOptions().fortran_showerprof_a.size() != 4 ||
            GetCalorimeter()->GetOptions().fortran_showerprof_b.size() != 4) {
            std::cout << "ReconstructionLednev::ReadCalibration: Size of one shower profile array is not equal to 4! (a: "
                      << GetCalorimeter()->GetOptions().fortran_showerprof_a.size()
                      << ", b: "
                      << GetCalorimeter()->GetOptions().fortran_showerprof_b.size()
                      << ")" << std::endl;
            exit(1);
            // TODO: should be fatal exception, but is caught somewhere
            throw Exception("ReconstructionLednev::ReadCalibration: Size of one shower profile array is not equal to 4! (a: %zu, b: %zu)",
                            GetCalorimeter()->GetOptions().fortran_showerprof_a.size(),
                            GetCalorimeter()->GetOptions().fortran_showerprof_b.size());
        }

        for (size_t i=0; i<4; i++) {
            showerprof_.a_[i] = GetCalorimeter()->GetOptions().fortran_showerprof_a[i];
            showerprof_.b_[i] = GetCalorimeter()->GetOptions().fortran_showerprof_b[i];
        }
    }
}

////////////////////////////////////////////////////////////////////////////////

} // namespace Reco
