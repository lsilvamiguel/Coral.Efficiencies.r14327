/*
   $Source: /afs/cern.ch/user/g/gobbo/w0/coral/src/calorim/Reco/src/Calibration.cc,v $
   $Date: 2011/01/31 20:35:45 $
   $Revision: 1.25 $
   -------------------------------------------------------------------------

   This file is part of cellular calorimeter reconstruction program.

   Authors:
     Vladimir  Kolosov   ( Kolosov@mx.ihep.su )
     Alexander Zvyagin   ( Alexander.Zviagine@cern.ch, Zvyagin@mx.ihep.su )

   Copyright(C): 1999-2000  V.Kolosov,A.Zvyagin

     This library is free software; you can redistribute it and/or
     modify it under the terms of the GNU Library General Public
     License as published by the Free Software Foundation; either
     version 2 of the License, or (at your option) any later version.

     This library is distributed in the hope that it will be useful,
     but WITHOUT ANY WARRANTY; without even the implied warranty of
     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
     Library General Public License for more details.

     You should have received a copy of the GNU Library General Public
     License along with this library; if not, write to the Free
     Software Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
*/

// --- Standard C/C++ library ---
#include <cassert>
#include <iostream>
#include <cstdlib>
#include <cmath>

#include "TROOT.h"
#include "TH1.h"
#include "TH2.h"
#include "TDirectory.h"
#include "TCanvas.h"
#include "TF1.h"

#include "myROOT_utils.h"

// --- Internal files ----
#include "Calorimeter.h"

#include "CalorimeterHist.h"
#include "CalorimeterParticle.h"
#include "CellDataRaw.h"
#include "CellType.h"
#include "Exception.h"
#include "OneParticleResponse.h"

using namespace std;

namespace Reco {

////////////////////////////////////////////////////////////////////////////////

bool Calorimeter::CalibrateSummEnergy(double esumm,const vector<OneParticleResponse> &hop)
{
//    bool debug=true;
   bool debug=false;
   if( !options.calibration_flag ) return false;
   double summ_measured = 0.;
   for( size_t i=0; i<hop.size(); i++ )
   {
      summ_measured += hop[i].GetParticle().GetE();
   }
   if(debug) cout << " Ebeam = " << esumm << " But measured " << summ_measured << endl;

    // Very far from desired value
   if( fabs(summ_measured-esumm)/esumm < 0.01 || fabs(summ_measured-esumm)/esumm > 100.) return false;

   double factor =  summ_measured/esumm;

   for( size_t i=0; i<hop.size(); i++ )
   {
      const vector<size_t> &c = hop[i].GetCells();
//       const vector<pair<double,double> > &e = hop[i].GetCellsEnergy();
      for( size_t j=0; j != c.size(); j++ )
      {
        // find this cell in calib_data_store
        for (vector<CellDataRaw>::const_iterator it=signals_.begin();
             it!=signals_.end(); it++) {
          if ( it->GetCellIdx()==c[j] ) {
            double weight = it->GetEnergy() / summ_measured;
            if( weight>0.1 )
            {
              cells_info[CALIB][NEW][c[j]].Add( factor , weight );
              cells_info[CALIB][MONITOR][c[j]].Add( factor , weight ); // new monitoring
              if(cells_stat_info[CALIB][MONITOR] != NULL)
                cells_stat_info[CALIB][MONITOR][c[j]].Add( factor , weight );
            }
          }
        }
      }
   }
   return true;
}

////////////////////////////////////////////////////////////////////////////////

void Calorimeter::CalibrateTime(void)
{
  if(!options.enable_time_calibration) return;
  for( size_t it=0; it != reco_particles.size(); it++ )
  {
    if(reco_particles[it].HasTime())
    {
      double time = reco_particles[it].GetTime();
      const vector<size_t> &hitedcells =reco_particles[it].GetMainCells();
      for(size_t i=0; i< hitedcells.size(); i++)
      {
        size_t incell=hitedcells[i];
        if( fabs(time) < options.calibration_time_gate && incell < NCells() )
          cells_info[TIME][NEW][incell].Add(time,1.);
      }
    }
  }
}

////////////////////////////////////////////////////////////////////////////////

void Calorimeter::CalibrateTime(const std::vector<CalorimeterParticle> &particles )
{
  if(!options.enable_time_calibration) return;
  for( size_t it=0; it != particles.size(); it++ )
  {
    if(particles[it].HasTime())
    {
      double time = particles[it].GetTime();
      const vector<size_t> &hitedcells =particles[it].GetHitedCells();
      for(size_t i=0; i< hitedcells.size(); i++)
      {
        size_t incell=hitedcells[i];
        if( fabs(time) < options.calibration_time_gate )
          cells_info[TIME][NEW][incell].Add(time,1.);
      }
    }
  }
}

////////////////////////////////////////////////////////////////////////////////

void Calorimeter::CalibrateAtMesonMass(size_t ipart,size_t jpart,double mass,double mass_sigma,double mass_tab,
                                                const vector<CalorimeterParticle> &particles)
{
//   if(!options.enable_mass_calibration) return;
  if( !options.calibration_flag ) return;
  assert(ipart<particles.size());
  assert(jpart<particles.size());
  double factor =  mass/mass_tab;
  const vector<size_t> main_cells = particles[ipart].GetMainCells();
  if(main_cells.size() ==0 ) return; // information about main cells not provided
  double e = particles[ipart].GetE();
  for( size_t it=0; it != main_cells.size(); it++ )
  {
    size_t mcell = main_cells[it];
    double weight = e;
//     cout << " cell  " << mcell <<" factor " << factor << " wight " << weight << endl;
    if(fabs(mass_tab-mass)/mass_tab < 0.3 )
    {
      cells_info[CALIB][NEW][mcell].Add( factor , weight );
      cells_info[CALIB][MONITOR][mcell].Add( factor , weight ); // new monitoring
      if(cells_stat_info[CALIB][MONITOR] != NULL)
        cells_stat_info[CALIB][MONITOR][mcell].Add( factor , weight );
    }
  }
}

////////////////////////////////////////////////////////////////////////////////

void Calorimeter::CalibrateParticle(size_t ipart,double energy,
                                    const vector<CalorimeterParticle> &particles,
                                             bool fill_default_calib_histo)
{
  if( !options.calibration_flag ) return;
  if(ipart >= particles.size()) assert(false);
  if( energy <= 0 ) return;
  const vector<size_t> main_cells = particles[ipart].GetMainCells();
  if(main_cells.size() ==0 ) return; // information about main cells not provided
  double e = particles[ipart].GetE();
  double factor = e/energy;
  if( factor < options.calibration_min || factor >= options.calibration_max ) return;
  for( size_t it=0; it != main_cells.size(); it++ )
  {
    size_t mcell = main_cells[it];
//    double weight = e;
    double weight = 1.;
    cells_info[CALIB][NEW][mcell].Add( factor , weight );
    cells_info[CALIB][MONITOR][mcell].Add( factor , weight ); // new monitoring
    if(cells_stat_info[CALIB][MONITOR] != NULL)
      cells_stat_info[CALIB][MONITOR][mcell].Add( factor , weight );
    if( fill_default_calib_histo )
    {
      FillCalibrationHisto( mcell, factor, weight);
    }
  }
}

////////////////////////////////////////////////////////////////////////////////

bool Calorimeter::CalibrateCell(size_t mcell,double factor,double weight,bool fill_default_calib_histo)
{
//  bool debug = true;
  bool debug = false;
  if( !options.calibration_flag ) return false;
  if( factor < options.calibration_min || factor >= options.calibration_max ) return false;
  if( mcell > NCells() )
  {
    cerr << "  ERROR!! Calorimeter::CalibrateCell for " << GetName() << " Bad cell index " << mcell << endl;
    return false;
  }
  if( debug ) cout << "  Calorimeter::CalibrateCell for " << GetName() << " cell " << mcell << " factor " <<
                       factor  << " weight " << weight << " make histo=" <<  fill_default_calib_histo << endl;
  cells_info[CALIB][NEW][mcell].Add( factor , weight );
  cells_info[CALIB][MONITOR][mcell].Add( factor , weight );  // new monitoring
  if(cells_stat_info[CALIB][MONITOR] != NULL)
      cells_stat_info[CALIB][MONITOR][mcell].Add( factor , weight );
  if( fill_default_calib_histo )
  {
    FillCalibrationHisto( mcell, factor, weight);
  }
  return true;
}

////////////////////////////////////////////////////////////////////////////////

bool Calorimeter::CalibrateParticle(size_t ipart, const vector<CalorimeterParticle> &particles)
{

//  bool debug = true;
  bool debug = false;
//   if( GetName() == "GDA" ) debug = true;
  if( debug ) cout << " CalibrateParticle " << GetName() <<" option " << options.calibration_flag << endl;
  if( !options.calibration_flag ) return false;
  if( debug ) cout << " CalibrateParticle " << GetName() << endl;
  if(ipart >= particles.size()) assert(false);
  double energy = options.calibration_energy;
  if( energy <= 0 ) return false;
  if( energy < options.calibration_energy_min ) return false;  // check options?? we need this only once!!
  if( energy >= options.calibration_energy_max ) return false; // check options ?? we need this only once!!

  const vector<size_t> main_cells = particles[ipart].GetMainCells();
  if(main_cells.size() ==0 ) return false; // information about main cells not provided
  double e = particles[ipart].GetE();
  if( debug ) cout << " Particle e=" << e <<
        " check for calibration Emin = " << options.calibration_energy_min <<
          " check for calibration Emax = " << options.calibration_energy_max << endl;
  if( e < options.calibration_energy_min ) return false;
  if( e >= options.calibration_energy_max ) return false;
  if( debug ) cout << " Particle e=" << e << " accepted for calibration " << endl;


  double factor = e/energy;
  for( size_t it=0; it != main_cells.size(); it++ )
  {
    size_t mcell = main_cells[it];
//    double weight = e;
    double weight = 1.;
    cells_info[CALIB][NEW][mcell].Add( factor , weight );
    cells_info[CALIB][MONITOR][mcell].Add( factor , weight ); // new monitoring
    if(cells_stat_info[CALIB][MONITOR] != NULL)
      cells_stat_info[CALIB][MONITOR][mcell].Add( factor , weight );
  }
  return true;
}

////////////////////////////////////////////////////////////////////////////////

int Calorimeter::InspectCalibration(void) const
{
//   if(!options.enable_mass_calibration) return;
//   if(!options.enable_time_calibration) return;
//   if(!options.calibration_flag) return;

  int printout_level = 1;

  int istat_max=5;
  int idelta_max=5;
  int iplus_max=2;
  int istat;
  int iplus;
  int idelta;

  int statistic[5][5][2];

  for( int istat=0; istat < istat_max; istat++ )
    for( int idelta=0; idelta < idelta_max; idelta++ )
      for( int iplus=0; iplus < iplus_max; iplus++ )
        statistic[istat][idelta][iplus]=0;


  for( size_t icell=0; icell < NCells(); icell++ )
  {
    istat =0;
    iplus =0;
    idelta=0;
    double stat_around = 0.;
    double dx = cells[icell].GetCellType().GetSizeX();
    double dy = cells[icell].GetCellType().GetSizeY();
    double s_main=dx*dy;
    for( vector<size_t>::const_iterator it=cells[icell].GetNeighbors().begin()+1;
                                        it!=cells[icell].GetNeighbors().end(); it++ )
    {
       double ddx = cells[*it].GetCellType().GetSizeX();
       double ddy = cells[*it].GetCellType().GetSizeY();
       stat_around += cells_info[CALIB][NEW][*it].GetEntries()*s_main/(ddx*ddy);
    }
    stat_around /= (cells[icell].GetNeighbors().size()-1);

    if(cells_info[CALIB][NEW][icell].GetMean()>0)
    {
      double factor  = cells_info[CALIB][NEW][icell].GetMean();
      double stat  = cells_info[CALIB][NEW][icell].GetEntries();
      double cfnew = cells_info[CALIB][OLD][icell].GetMean()/cells_info[CALIB][NEW][icell].GetMean();
//      double cfold = cells_info[CALIB][OLD][icell].GetMean();
      double delta = fabs(factor-1);
      if     (stat > 1000)  istat=4;
      else if(stat >  200)  istat=3;
      else if(stat >   50)  istat=2;
      else if(stat >   10)  istat=1;

      if     (delta < 0.01)  idelta=4;
      else if(delta < 0.05)  idelta=3;
      else if(delta < 0.10)  idelta=2;
      else if(delta < 0.50)  idelta=1;

      if(factor >1)  iplus = 1;

      statistic[istat][idelta][iplus]++;
      if( printout_level >= 0 )
      {
        if(cfnew < 0) cout << " BRED in " << icell << "CELL CALIBRATION " << GetCellName(icell) << endl;
      }
      if( printout_level > 1 )
      {
        if( stat_around > 2*cells_info[CALIB][NEW][icell].GetEntries() )
          cout << " RELATIVELY LOW STATISTIC " <<  cells_info[CALIB][NEW][icell].GetEntries() <<
                      "   STAT AROUND " << stat_around  << " in CELL " << icell <<
                                     " name: " << GetCellName(icell) << endl;
      }
    }
    else
    {
      if( printout_level > 1 )
      {
        cout << " EMPTY CALIBRATION in " << icell << " CELL " << GetCellName(icell) << endl;
      }
      if(cells_info[CALIB][OLD][icell].GetMean()>0)
      {
//         cout << " use Old calibration " << endl;
      }
      else
      {
//         cout << " use Default calibration " << endl;
      }
    }
  }

  int total=0;
  for( int istat=0; istat < istat_max; istat++ )
    for( int idelta=0; idelta < idelta_max; idelta++ )
      for( int iplus=0; iplus < iplus_max; iplus++ )
        total +=statistic[istat][idelta][iplus];

  if( printout_level >= 0 )
  {
    cout << " **** CALORIMETER " << GetName() << " InspectCalibration STATISTIC Total=" << total << endl;
    if( printout_level > 0 && total > 0 )
    {
      cout << "---------------------------------------------------------------------------------------------------------" << endl;
      cout << " Deviation           I        >50%            10-50%          1-10%             <5%             <1%      " << endl;
      cout << "                     I      -       +         -      +       -       +       -      +        -      +    " << endl;
      cout << "---------------------------------------------------------------------------------------------------------" << endl;
      for( int istat=0; istat < istat_max; istat++ )
      {
        if(istat==0) cout << " Statistic <   10    I ";
        if(istat==1) cout << " Statistic <   50    I ";
        if(istat==2) cout << " Statistic <  200    I ";
        if(istat==3) cout << " Statistic < 1000    I ";
        if(istat==4) cout << " Statistic > 1000    I ";

        for( int idelta=0; idelta < idelta_max; idelta++ )
          for( int iplus=0; iplus < iplus_max; iplus++ )
          {
            printf(" %6d ",statistic[istat][idelta][iplus]);
          }
        cout << endl;
      }
      cout << "-----------------------------------------------------------------------------------------" << endl;
    }
    cout << endl;
  }
  return 0;

}

////////////////////////////////////////////////////////////////////////////////

void Calorimeter::MonitorCalibration(void)
{
  bool print_cells_info = false;
  bool print_summary = true;
  if(cells_stat_info[CALIB][MONITOR]==NULL) return;

  int istat_max=5;
  int idelta_max=4;
  int iplus_max=2;
  int istat;
  int iplus;
  int idelta;

  int statistic_empty=0;
  int statistic_relatively_low=0;

  int statistic[5][4][2];

  for( int istat=0; istat < istat_max; istat++ )
    for( int idelta=0; idelta < idelta_max; idelta++ )
      for( int iplus=0; iplus < iplus_max; iplus++ )
        statistic[istat][idelta][iplus]=0;


  for( size_t icell=0; icell < NCells(); icell++ )
  {
    istat =0;
    iplus =0;
    idelta=0;
    double stat_around = 0.;
    double dx = cells[icell].GetCellType().GetSizeX();
    double dy = cells[icell].GetCellType().GetSizeY();
    double s_main=dx*dy;
    for( vector<size_t>::const_iterator it=cells[icell].GetNeighbors().begin()+1;
                                        it!=cells[icell].GetNeighbors().end(); it++ )
    {
       double ddx = cells[*it].GetCellType().GetSizeX();
       double ddy = cells[*it].GetCellType().GetSizeY();
       stat_around += cells_stat_info[CALIB][MONITOR][*it].GetEntries()*s_main/(ddx*ddy);
    }
    stat_around /= (cells[icell].GetNeighbors().size()-1);

    if(cells_stat_info[CALIB][MONITOR][icell].GetMean()>0)
    {
      double factor  = cells_stat_info[CALIB][MONITOR][icell].GetMean();
      double stat  = cells_stat_info[CALIB][MONITOR][icell].GetEntries();
      double cfnew = cells_info[CALIB][OLD][icell].GetMean()/cells_stat_info[CALIB][MONITOR][icell].GetMean();
//      double cfold = cells_info[CALIB][OLD][icell].GetMean();
      double delta = fabs(factor-1);
      if     (stat > 1000)  istat=4;
      else if(stat >  200)  istat=3;
      else if(stat >   50)  istat=2;
      else if(stat >   10)  istat=1;

      if     (delta < 0.01)  idelta=3;
      else if(delta < 0.10)  idelta=2;
      else if(delta < 0.50)  idelta=1;

      if(factor >1)  iplus = 1;

      statistic[istat][idelta][iplus]++;
      if(cfnew < 0) cout << " BRED in " << icell << "CELL CALIBRATION " << GetCellName(icell) << endl;
      if( stat_around > 2*cells_stat_info[CALIB][MONITOR][icell].GetEntries() )
      {
        statistic_relatively_low++;
        if(print_cells_info) cout << " RELATIVELY LOW STATISTIC " <<
                       cells_stat_info[CALIB][MONITOR][icell].GetEntries() <<
                      "   STAT AROUND " << stat_around  << " in CELL " << icell <<
                                     " name: " << GetCellName(icell) << endl;
      }
    }
    else
    {
      statistic_empty++;
      if(print_cells_info) cout << " EMPTY CALIBRATION in " << icell << " CELL " << GetCellName(icell) << endl;
      if(cells_info[CALIB][OLD][icell].GetMean()>0)
      {
//         cout << " use Old calibration " << endl;
      }
      else
      {
//         cout << " use Default calibration " << endl;
      }
    }
  }

  int total=0;
  for( int istat=0; istat < istat_max; istat++ )
    for( int idelta=0; idelta < idelta_max; idelta++ )
      for( int iplus=0; iplus < iplus_max; iplus++ )
        total +=statistic[istat][idelta][iplus];

  if(print_summary)
  {
    cout << " **** CALORIMETER " << GetName() << " MonitorCalibration STATISTIC Total=" << total << endl;
    cout << " EMPTY STATISTIC " << statistic_empty << endl;
    cout << " RELATIVELY LOW STATISTIC " << statistic_relatively_low << endl;
    cout << " Calibration     MORE   50% " << "      10-50%    " << "     1-10%   " << "      LESS 1%   " << endl;

    for( int istat=0; istat < istat_max; istat++ )
    {
      if(istat==0) cout << " Statistic <   10    I " << endl;
      if(istat==1) cout << " Statistic <   50    I " << endl;
      if(istat==2) cout << " Statistic <  200    I " << endl;
      if(istat==3) cout << " Statistic < 1000    I " << endl;
      if(istat==4) cout << " Statistic > 1000    I " << endl;

      cout << "                 I ";
      for( int idelta=0; idelta < idelta_max; idelta++ )
        for( int iplus=0; iplus < iplus_max; iplus++ )
          printf(" %6d ",statistic[istat][idelta][iplus]);
      cout << endl;
    }
  }

}

////////////////////////////////////////////////////////////////////////////////

void Calorimeter::ClearStatInfo(size_t what,size_t when)
{
  if( cells_stat_info[what][when] == NULL )
  {
    cerr << " Calorimeter " << GetName() << " ClearStatInfo: cells_stat_info[" <<
                                    what << "][" << when <<"] not initialized " << endl;
    return;
  }
  for( size_t icell=0; icell < NCells(); icell++ )
  {
    cells_stat_info[what][when][icell].Clear();
  }
}

////////////////////////////////////////////////////////////////////////////////

void Calorimeter::FillCalibrationHisto(size_t icell, double factor, double weight=1.)
{
//   bool debug=true;
  bool debug=false;
  if(debug) cout << " Calorimeter " << GetName() << " FillCalibrationHisto " << endl;
  if( !options.calibration_flag ) return;
  bool ok=false;

  TDirectory *dir_save=gDirectory; // Save current ROOT directory.
  if( calo_hist==NULL ) calo_hist = new CalorimeterHist(this);
  CalorimeterHist &h= *calo_hist;
//  CalorimeterCalibHist &hc=  *calib_hist;  // create a short name for calib_hist
  if( h.root_dir==NULL )
  {
    char dir_name[111];
    ok=(gDirectory->cd("/"));
    assert(ok);
    sprintf(dir_name,"Calorimeter_%s",GetName().c_str());
    h.root_dir = myROOT_utils::TDirectory_checked(dir_name);
  }
  ok=(h.root_dir->cd());
  assert(ok);
  if( h.cells_calib_dir==NULL )
  {
    char dir_name[111];
    sprintf(dir_name,"CellsCalibration_%s",GetName().c_str());
    h.cells_calib_dir = myROOT_utils::TDirectory_checked(dir_name);
    ok=(h.cells_calib_dir->cd());
    assert(ok);

    char hist_name[132];

    options.calib_histo_units_in_gev=false;
    for( size_t i=0; i< NCells(); i++ )
    {
      char name[132];
      sprintf(name,"Calib%s",GetCellName(i).c_str());
      sprintf(hist_name," Default calibration in cell %s ",GetCellName(i).c_str());
      h.h1_CalibCell[i] = myROOT_utils::TH1D_checked(name,hist_name, options.calib_histo_size,
                                                                    0.,options.calib_histo_range);
    }
    if(debug) cout << " Calorimeter " << GetName() << " FillCalibrationHisto Booking OK" << endl;
  }
  ok=(h.cells_calib_dir->cd());
  assert(ok);

  if ( icell >= NCells() ) assert( false);

  if(h.h1_CalibCell[icell]!=NULL) h.h1_CalibCell[icell]->Fill(factor,weight);
  else cout << " h.h1_CalibCell["<< icell <<"]=NULL"<<endl;

  if( dir_save != NULL )dir_save->cd();
  if(debug) cout << " Calorimeter " << GetName() << " FillCalibrationMesonHisto OK" << endl;
}

////////////////////////////////////////////////////////////////////////////////

void Calorimeter::FillCalibrationMesonHisto(size_t ipart,size_t jpart,double mass,double mas_tab,
                                                const vector<CalorimeterParticle> &particles)
{
//   bool debug=true;
  bool debug=false;
  if(debug) cout << " Calorimeter " << GetName() << " FillCalibrationMesonHisto " << endl;
  if( !options.calibration_flag ) return;
  bool ok=false;

  TDirectory *dir_save=gDirectory; // Save current ROOT directory.
  if( calo_hist==NULL ) calo_hist = new CalorimeterHist(this);
  CalorimeterHist &h= *calo_hist;
//  CalorimeterCalibHist &hc=  *calib_hist;  // create a short name for calib_hist
  if( h.root_dir==NULL )
  {
    char dir_name[111];
    ok=(gDirectory->cd("/"));
    assert(ok);
    sprintf(dir_name,"Calorimeter_%s",GetName().c_str());
    h.root_dir = myROOT_utils::TDirectory_checked(dir_name);
  }
  ok=(h.root_dir->cd());
  assert(ok);
  if( h.cells_calib_dir==NULL )
  {
    char dir_name[111];
    sprintf(dir_name,"CellsCalibration_%s",GetName().c_str());
    h.cells_calib_dir = myROOT_utils::TDirectory_checked(dir_name);
    ok=(h.cells_calib_dir->cd());
    assert(ok);

    char hist_name[132];

    options.calib_histo_units_in_gev=true;
//    double emin = 0, emax = options.hist_energy_max;
    for( size_t i=0; i< NCells(); i++ )
    {
      char name[132];
      sprintf(name,"Calib%s",GetCellName(i).c_str());
      sprintf(hist_name," Mass in cell %s ",GetCellName(i).c_str());
      h.h1_CalibCell[i] = myROOT_utils::TH1D_checked(name,hist_name,options.calib_histo_size,0.,1.);
    }
    if(debug) cout << " Calorimeter " << GetName() << " FillCalibrationMesonHisto Booking OK" << endl;
  }
  ok=(h.cells_calib_dir->cd());
  assert(ok);

  if(debug) cout << " Calorimeter " << GetName() << " FillCalibrationMesonHisto Nparticles " <<
                                       particles.size() << " i=" << ipart << " j=" << jpart  << endl;
  assert(ipart < particles.size());
  assert(jpart < particles.size());

  if(debug) cout << " Calorimeter " << GetName() << " FillCalibrationMesonHisto Particle " <<
                             ipart <<" cells " <<  particles[ipart].GetMainCells().size() << endl;
  for( size_t it=0; it< particles[ipart].GetMainCells().size(); it++)
  {
    size_t cell = particles[ipart].GetMainCells()[it];
    assert( cell < NCells() );
    if(h.h1_CalibCell[cell]!=NULL) h.h1_CalibCell[cell]->Fill(mass);
    else cout << " h.h1_CalibCell["<< cell <<"]=NULL"<<endl;
  }

  if(debug) cout << " Calorimeter " << GetName() << " FillCalibrationMesonHisto Particle " <<
                             jpart <<" cells " <<  particles[ipart].GetMainCells().size() << endl;
  for( size_t it=0; it< particles[jpart].GetMainCells().size(); it++)
  {
    size_t cell = particles[jpart].GetMainCells()[it];
    assert( cell < NCells() );
    if(debug) cout << " Calibrate cell " << cell  << " with mass " << mass << endl;
    if(h.h1_CalibCell[cell]!=NULL) h.h1_CalibCell[cell]->Fill(mass);
    else cout << " h.h1_CalibCell["<< cell <<"]=NULL"<<endl;
  }
  if( dir_save != NULL )dir_save->cd();
  if(debug) cout << " Calorimeter " << GetName() << " FillCalibrationMesonHisto OK" << endl;
}

////////////////////////////////////////////////////////////////////////////////

void Calorimeter::FillCalibrationHisto(size_t ipart,const vector<CalorimeterParticle> &particles)
{
  if( !options.calibration_flag ) return;
  if(ipart >= particles.size()) assert(false);
  double energy = options.calibration_energy;
  if( energy <= 0 ) return;
  if( energy < options.calibration_energy_min ) return;
  if( energy >= options.calibration_energy_max ) return;
  double e = particles[ipart].GetE();
  if( e < options.calibration_energy_min ) return;
  if( e >= options.calibration_energy_max ) return;
  FillCalibrationHisto(ipart, energy, particles);
}

////////////////////////////////////////////////////////////////////////////////

void Calorimeter::FillCalibrationHisto(size_t ipart,double energy,
                                 const vector<CalorimeterParticle> &particles)
{
//   bool debug=true;
  bool debug=false;
  if(debug) cout << " Calorimeter " << GetName() << " FillCalibrationMesonHisto " << endl;
  if( !options.calibration_flag ) return;
  bool ok=false;

  TDirectory *dir_save=gDirectory; // Save current ROOT directory.
  if( calo_hist==NULL ) calo_hist = new CalorimeterHist(this);
  CalorimeterHist &h= *calo_hist;
//  CalorimeterCalibHist &hc=  *calib_hist;  // create a short name for calib_hist
  if( h.root_dir==NULL )
  {
    char dir_name[111];
    ok=(gDirectory->cd("/"));
    assert(ok);
    sprintf(dir_name,"Calorimeter_%s",GetName().c_str());
    h.root_dir = myROOT_utils::TDirectory_checked(dir_name);
  }
  ok=(h.root_dir->cd());
  assert(ok);
  if( h.cells_calib_dir==NULL )
  {
    char dir_name[111];
    sprintf(dir_name,"CellsCalibration_%s",GetName().c_str());
    h.cells_calib_dir = myROOT_utils::TDirectory_checked(dir_name);
    ok=(h.cells_calib_dir->cd());
    assert(ok);

    char hist_name[132];

    options.calib_histo_units_in_gev=true;
//    double emin = 0, emax = options.hist_energy_max;
    for( size_t i=0; i< NCells(); i++ )
    {
      char name[132];
      sprintf(name,"Calib%s",GetCellName(i).c_str());
      sprintf(hist_name," Calibration in cell %s ",GetCellName(i).c_str());
      h.h1_CalibCell[i] = myROOT_utils::TH1D_checked(name,hist_name,options.calib_histo_size,0.,options.hist_energy_max);
      if( options.fill_calib_spill_correlations_histo )
      {
        sprintf(name,"CalibInSpill%s",GetCellName(i).c_str());
        sprintf(hist_name," EvInSpill vs Calibration in cell %s ",GetCellName(i).c_str());
	double evmax_in_spill = 10000.;
        h.h2_CalibInSpill[i] = myROOT_utils::TH2D_checked(name,hist_name,options.calib_histo_size,0.,options.hist_energy_max,50, 0., evmax_in_spill);
      }
    }
    if(debug) cout << " Calorimeter " << GetName() << " FillCalibrationHisto Booking OK" << endl;
  }
  ok=(h.cells_calib_dir->cd());
  assert(ok);

  if(debug) cout << " Calorimeter " << GetName() << " FillCalibrationHisto Nparticles " <<
                                       particles.size() << " i=" << ipart <<  endl;
  assert(ipart<particles.size());
  if( energy <= 0 ) return;
  const vector<size_t> main_cells = particles[ipart].GetMainCells();
  if(main_cells.size() ==0 ) return; // information about main cells not provided
  double e = particles[ipart].GetE();
  for( size_t it=0; it != main_cells.size(); it++ )
  {
    size_t mcell = main_cells[it];
    if(h.h1_CalibCell[mcell]!=NULL) h.h1_CalibCell[mcell]->Fill(e);
    else cout << " h.h1_CalibCell["<< mcell <<"]=NULL"<<endl;
    if( options.fill_calib_spill_correlations_histo )
    {
      int evnum = 0;
      if( GetEventID().Initialized() ) evnum = GetEventID().GetEventNumberInBurst();
      if( h.h2_CalibInSpill[mcell] != NULL ) h.h2_CalibInSpill[mcell]->Fill( e, evnum);
    }
  }

  if( dir_save != NULL )dir_save->cd();
  if(debug) cout << " Calorimeter " << GetName() << " FillCalibrationHisto OK" << endl;
}

////////////////////////////////////////////////////////////////////////////////

StatInfoStore *Calorimeter::GetStatInfoStore(size_t what,size_t when,size_t cell)
{
  if(cells_store[what][when] == NULL) return NULL;
  else return &cells_store[what][when][cell];
}

////////////////////////////////////////////////////////////////////////////////

void Calorimeter::InitMonitoring ( int run_start, int run_finish )
{
  for( size_t icell=0; icell!=NCells(); icell++ )
  {
    cells_store[LED][OLD][icell].Init( run_start, run_finish );
    cells_store[CALIB][MONITOR][icell].Init( run_start, run_finish );
  }

  LEDTotal.Init    ( run_start, run_finish );
  CALIBTotal.Init  ( run_start, run_finish );
}

////////////////////////////////////////////////////////////////////////////////

void Calorimeter::FillMonitorHisto ( void )
{
//   bool debug=true;
  bool debug=false;
//  bool find_intervals = true;
  bool find_intervals = false;
  int max_number_of_intervals = 6;
  double hi_break_led = 25;
  double hi_penalty_led = 10;
// Set minimal range, minimal statistic, minimal sigma (additive), minimal sigma (fraction) for LEDs
  int min_range_led = 10;
  int min_stat_led =  10;
  double sigma_min_add_led = 2;
  double sigma_min_frac_led = 0.005;
// Set minimal range, minimal statistic, minimal sigma (additive), minimal sigma (fraction) for CALIBs
  int min_range_calib = 10;
  int min_stat_calib =   5;
  double sigma_min_add_calib = 0.01;
  double sigma_min_frac_calib = 0.005;
  string led_histo_path = "MonitorCellsLED_" + GetName();

  if( debug ) cout << " ********************************* Calorimeter::FillMonitorHisto  " <<
                  GetName() << " debug NCells = " << NCells() << endl;

  if( debug ) cout << " Path for leds " << led_histo_path << endl;
  bool skip_leds = false;
  int monitor_led_histo_level = 1;

  if(cells_store[LED][OLD] != NULL && !skip_leds )
  {
    for( size_t icell=0; icell != NCells(); icell++ )
    {
      if( debug ) cout << " cells_store[LED][OLD][ " << icell <<" ]" << endl;
      cells_store[LED][OLD][icell].SetTogether( min_range_led, min_stat_led, sigma_min_add_led, sigma_min_frac_led );

      cells_store[LED][OLD][icell].BookHisto( "LED_CELL_"+GetCellName(icell), led_histo_path, monitor_led_histo_level );
      if( find_intervals ) cells_store[LED][OLD][icell].FindIntervals(max_number_of_intervals, hi_break_led, hi_penalty_led );
      cells_store[LED][OLD][icell].FillHisto();
      if( find_intervals && options.monitor_histo_show_fit) cells_store[LED][OLD][icell].AddFitToHisto();
    }
  }

  string calib_histo_path = "MonitorCellsCALIB_" + GetName();
//  bool debug_ecal2 = false;

  if( debug ) cout << " Path for calib " << calib_histo_path << endl;
  int monitor_calib_histo_level = 1;

  if(cells_store[CALIB][MONITOR] != NULL )
  {
    for( size_t icell=0; icell != NCells(); icell++ )
    {
      if( debug ) cout << " cells_store[CALIB][MONITOR][ " << icell <<" ]" << endl;
      cells_store[CALIB][MONITOR][icell].SetTogether( min_range_calib, min_stat_calib, sigma_min_add_calib, sigma_min_frac_calib );

      cells_store[CALIB][MONITOR][icell].BookHisto( "CALIB_CELL_"+GetCellName(icell), calib_histo_path, monitor_calib_histo_level );
      if( find_intervals ) cells_store[CALIB][MONITOR][icell].SetIntervals( cells_store[LED][OLD][icell].GetIntervals() );
      if( find_intervals ) cells_store[CALIB][MONITOR][icell].FitInIntervals();
      cells_store[CALIB][MONITOR][icell].FillHisto();
//       if( debug) cout << " cell " << icell << " cells_store[CALIB][MONITOR][icell].GetEntries(0) entries " <<
//                                              cells_store[CALIB][MONITOR][icell].GetStatInfo()[0].GetEntries() << endl;

      if( find_intervals && options.monitor_histo_show_fit) cells_store[CALIB][MONITOR][icell].AddFitToHisto();
    }
  }

  LEDTotal.SetTogether( min_range_led, min_stat_led, sigma_min_add_led, sigma_min_frac_led );

  CALIBTotal.SetTogether( min_range_calib, min_stat_calib, sigma_min_add_calib, sigma_min_frac_calib );

  string histo_path = "MonitorGeneral_" + GetName();
  LEDTotal.BookHisto    ( "LEDTotal", histo_path, 0 );

  CALIBTotal.BookHisto    ( "CALIBTotal", histo_path, 0 );

  if(cells_store[LED][OLD] == NULL )
  {
    cerr << " Error Calorimeter::FillMonitorHisto  cells_store[LED][OLD] == NULL " << endl;
    return;
  }

  if(cells_store[CALIB][MONITOR] == NULL )
  {
    cerr << " Error Calorimeter::FillMonitorHisto  cells_store[CALIB][MONITOR] == NULL " << endl;
    return;
  }

  for( size_t icell=0; icell!=NCells(); icell++ )
  {
      LEDTotal += cells_store[LED][OLD][icell];
      CALIBTotal += cells_store[CALIB][MONITOR][icell];
  }

  if( find_intervals ) LEDTotal.FindIntervals(max_number_of_intervals, hi_break_led, hi_penalty_led );
  LEDTotal.FillHisto();
  if( find_intervals && options.monitor_histo_show_fit) LEDTotal.AddFitToHisto();

  if( find_intervals ) CALIBTotal.SetIntervals( LEDTotal.GetIntervals() );
  if( find_intervals ) CALIBTotal.FitInIntervals();
  CALIBTotal.FillHisto();
  if( find_intervals && options.monitor_histo_show_fit) CALIBTotal.AddFitToHisto();

}

////////////////////////////////////////////////////////////////////////////////

void Calorimeter::StoreMonitorInfo ( int run )
{
  if(cells_store[CALIB][MONITOR] == NULL )
  {
    cerr << " ERROR Calorimeter::StoreMonitorInfo cells_store[CALIB][MONITOR] = NULL " << endl;
    assert(false);
  }
  if(cells_store[LED][OLD] == NULL )
  {
    cerr << " ERROR Calorimeter::StoreMonitorInfo cells_store[LED][OLD] = NULL " << endl;
    assert(false);
  }
//   for( size_t icell=0; icell!=NCells(); icell++ )
//   {
//     if(cells_store[CALIB][MONITOR] != NULL)
//       cells_store[CALIB][MONITOR][icell].SetInfo( run, cells_info[CALIB][MONITOR][icell] );
//     if(cells_store[LED][OLD] != NULL)
//       cells_store[LED][OLD][icell].SetInfo( run, cells_info[LED][OLD][icell] );
//   }
  for( size_t icell=0; icell!=NCells(); icell++ )
  {
    cells_store[CALIB][MONITOR][icell].SetInfo( run, cells_info[CALIB][MONITOR][icell] );
    cells_store[LED][OLD][icell].SetInfo( run, cells_info[LED][OLD][icell] );
  }
}

////////////////////////////////////////////////////////////////////////////////

void Calorimeter::CalculateMonitorCalibration( int run )
{
//  double stat_min = 50;
//   cout << " Calculate calibrations for calorimeter " << GetName() << " for RUN= " << run << endl;
  if(cells_store[CALIB][MONITOR] == NULL) return;
  for( size_t icell=0; icell < NCells(); icell++ )
  {
//     cout << " cell " << icell << " average  in interval " <<
//            cells_store[CALIB][MONITOR][icell].CalculateAverageInInterval( run ) << endl;
//     cout << " cell " << icell <<
//            " Bin value " << cells_store[CALIB][MONITOR][icell].GetBinValue( run ) <<
//            " Fit value " << cells_store[CALIB][MONITOR][icell].CalculateFitValue( run ) << endl;
    StatInfo average = cells_store[CALIB][MONITOR][icell].CalculateAverageInInterval( run );
    double stat = average.GetEntries();
// #warning " TODO Check FitValue is valid "
    double correction = cells_store[CALIB][MONITOR][icell].CalculateFitValue( run );
    double sigma = average.GetSigma();
    cells_info[CALIB][NEW][icell] = StatInfo( stat, correction, sigma);
  }
}

////////////////////////////////////////////////////////////////////////////////

} /// using namespace Reco
