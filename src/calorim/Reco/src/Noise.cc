/*
   $Source: /afs/cern.ch/user/g/gobbo/w0/coral/src/calorim/Reco/src/Noise.cc,v $
   $Date: 2011/01/31 20:35:45 $
   $Revision: 1.22 $
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
#include <set>

// --- Internal files ----
#include "Calorimeter.h"

#include "myROOT_utils.h"
#include "CalorimeterHist.h"
#include "CellDataRaw.h"
#include "CellType.h"

using namespace std;

namespace Reco {

////////////////////////////////////////////////////////////////////////////////

void Calorimeter::ClearCellsNoise(int icell)
{
  if (!options.store_noise_info) return; //  Clear noise info
  if(icell >= 0)
  {
    assert( icell < (int)NCells() );
    cells_info[NOISE][NEW][icell].Clear();
    cells_info[CLUSTER][NEW][icell].Clear();
  }
  else
  {
    cells_noise_statistic=0;
    cells_gamma_noise_statistic=0;
    for( size_t i=0; i<NCells(); i++ )
    {
      cells_info[NOISE][NEW][i].Clear();
      cells_info[CLUSTER][NEW][i].Clear();
    }
  }

  if (!options.store_noise_histo) return; // Reset noise histo
  assert( calo_hist != NULL );
  CalorimeterHist &h=*calo_hist;  // create a short name for calo_hist
  if(icell >= 0)
  {
    if(h.h1_NOISECell[icell]    !=NULL) h.h1_NOISECell[icell]->Reset();
    if(h.h1_NOISEGamCell[icell] !=NULL) h.h1_NOISEGamCell[icell]->Reset();
  }
  else
  {
    for( size_t i=0; i<NCells(); i++ )
    {
      if(h.h1_NOISECell[i]    !=NULL) h.h1_NOISECell[i]->Reset();
      if(h.h1_NOISEGamCell[i] !=NULL) h.h1_NOISEGamCell[i]->Reset();
    }
  }
}

////////////////////////////////////////////////////////////////////////////////

int Calorimeter::InspectCellsGamNoise(void)
{
  bool debug = false;
//   bool debug = true;
  if( debug ) cout << "Calorimeter::InspectCellsGamNoise " << GetName() <<" option " << options.store_gamma_noise_info << endl;
  if( !options.store_gamma_noise_info ) return 0;
  int printout_level = 1;

  double statistic=cells_gamma_noise_statistic;
  if( debug ) cout << "Calorimeter::" << GetName() <<
                      "::InspectCellsGamNoise() statistic=" << statistic << endl;
  double stat_min = 100.;
  assert( calo_hist != NULL );
  CalorimeterHist &h=*calo_hist;  // create a short name for calo_hist
  int empty_statistic = 0;
  int low_statistic = 0;
  int need_up_statistic = 0;
  int need_down_now_statistic = 0;
  int need_down_later_statistic = 0;
  int normal_statistic = 0;

  for( size_t i=0; i<NCells(); i++ )
  {
//     double stat_ngam_in_cell = h.h1_NOISEGamCell[i]->GetEntries();
     double stat_ngam_in_cell = cells_info[CLUSTER][NEW][i].GetEntries();
     if( stat_ngam_in_cell <= 0 ) empty_statistic++;

     double stat_ngam_around = 0.;
     double dx = cells[i].GetCellType().GetSizeX();
     double dy = cells[i].GetCellType().GetSizeY();
     double s_main=dx*dy;
     for( vector<size_t>::const_iterator it=cells[i].GetNeighbors().begin()+1;
                                         it!=cells[i].GetNeighbors().end(); it++ )
     {
       double ddx = cells[*it].GetCellType().GetSizeX();
       double ddy = cells[*it].GetCellType().GetSizeY();
//       stat_ngam_around += h.h1_NOISEGamCell[*it]->GetEntries();
       stat_ngam_around += cells_info[CLUSTER][NEW][*it].GetEntries()*s_main/(ddx*ddy);
     }
     double sigma_stat_ngam_in_cell = sqrt(stat_ngam_in_cell);
     stat_ngam_around /= cells[i].GetNeighbors().size()-1;
     double sigma_stat_ngam_around = sqrt(stat_ngam_around);
     double delta = stat_ngam_in_cell-stat_ngam_around;
     if( stat_ngam_in_cell <= stat_min )
     {
       if( stat_ngam_in_cell > 0 ) low_statistic++;
       if( debug) cout << " Low gamma statistic in cell " << GetCellName(i) << " N=" << stat_ngam_in_cell << endl;
     }
     else
     if(   (delta > 4*(sigma_stat_ngam_in_cell+sigma_stat_ngam_around))&&( delta > 0.1*stat_ngam_in_cell))
     {
// # warning " Calorimeter::InspectCellsNoise:: This cut is only for GAMS valid !!! "
//       if( fabs(cells[i].GetX()) < 38 && fabs(cells[i].GetY()) < 38 ) continue; // Very central cells in GAMS is a very special case

       if( debug )
       {
         cout << " INSPECT GAMMMA NOISE in Calorimeter " << GetName() << " cell  " << GetCellName(i) << " Stat Over NOISY? " << endl;
         cout << " Stat " << stat_ngam_in_cell << " Stat around   " << stat_ngam_around << " Relax Gate " <<
                                                     4*(sigma_stat_ngam_in_cell+sigma_stat_ngam_around) << endl;
       }
       TAxis *axis = h.h1_NOISEGamCell[i]->GetXaxis();
       Int_t ncx   = axis->GetNbins();
//       Double_t xmin = axis->GetXmin();
//       Double_t xmax = axis->GetXmax();
       Double_t sum = 0;
       Double_t sum_before = 0;
       Int_t bin;
       for (bin=ncx+1;bin>=1;bin--)
       {
         sum += h.h1_NOISEGamCell[i]->GetBinContent(bin);
         if( sum > stat_ngam_around ) break;
         sum_before=sum;
       }
       double prob_old=stat_ngam_in_cell/statistic;
       double prob_new=sum/statistic;
       double prob_before=sum_before/statistic;
       double prob_around=stat_ngam_around/statistic;
       if(prob_new > 1.5*prob_around && bin < ncx+1)
       {
         prob_new=prob_before;
         bin += 1;
       }

       Axis_t ecut = h.h1_NOISEGamCell[i]->GetBinCenter(bin);
       if( debug )
       {
         cout << " INSPECT GAMMMA NOISE Recommended cut " << ecut << " Prob with old cut "  <<  prob_old <<
               " Prob above new cut "  <<  prob_new << " Prob around " << prob_around << endl;
       }
       energy_gamma_cut_bad_cells_new[i]=ecut;
       if(!cell_is_bad_[NEW][i])
       {
         if( debug ) cout << " Set cell " << GetCellName(i) << " as bad new ! " << endl;
         bad_cells_[NEW].push_back(i);
         cell_is_bad_[NEW][i]=true;
       }
       if(cell_is_bad_[OLD][i])
       {
          if(energy_gamma_cut_bad_cells_new[i] >  energy_gamma_cut_bad_cells_old[i]) cout << " situation with GAMMA NOISE in cell " << i << " getting worse " << endl;
          if(energy_gamma_cut_bad_cells_new[i] == energy_gamma_cut_bad_cells_old[i]) cout << " situation with GAMMA NOISE in cell " << i << " is stable " << endl;
          if(energy_gamma_cut_bad_cells_new[i] <  energy_gamma_cut_bad_cells_old[i]) cout << " situation with GAMMA NOISE in cell " << i << " slightly getting better " << endl;
       }
       need_up_statistic++;
     }
     else
     if( cell_is_bad_[OLD][i] && !cell_is_bad_[NEW][i] && ( fabs(delta) <= 0.1*stat_ngam_in_cell)) // The situation with bad cells is more or less stable
     {
       cout << " Cell was bad But now ????? " << endl;
       bad_cells_[NEW].push_back(i);
       cell_is_bad_[NEW][i]=true;
       energy_gamma_cut_bad_cells_new[i]=energy_gamma_cut_bad_cells_old[i];
       if( debug )
         cout << " Cell " << i << " stay with gamma ecut " << energy_cut_bad_cells_new[i] << endl;
     }
     else
     if( ( delta < -4*(sigma_stat_ngam_in_cell+sigma_stat_ngam_around))&&(fabs(delta)>0.1*stat_ngam_in_cell))
     {

//        cout << " INSPECT GAMMA NOISE in Calorimeter DEAD CELL? " << GetName() << " cell  " << i << " x= " << cells[i].GetX() << " y= " << cells[i].GetY() <<endl;
//        cout << " Stat " << stat_ngam_in_cell << " Stat around   " << stat_ngam_around << " Stat check   " << cells_info[NOISE][NEW][i].GetEntries() << endl;
//        cout << " Gate " << 4*(sigma_stat_ngam_in_cell+sigma_stat_ngam_around) << endl;
       if(cell_is_bad_[OLD][i])
       {
         need_down_now_statistic++;
       }
       else
       {
         need_down_later_statistic++;
       }

       if( debug )
       {
         cout << " INSPECT GAMMA NOISE in Calorimeter " << GetName() << " cell  " << GetCellName(i) << "LOW NOISE ? "<<endl;
         cout << " Stat " << stat_ngam_in_cell << " Stat around   " << stat_ngam_around << " Stat check   " <<
               cells_info[CLUSTER][NEW][i].GetEntries() << " Gate " << 4*(sigma_stat_ngam_in_cell+sigma_stat_ngam_around) << endl;
         if(cell_is_bad_[OLD][i]) cout << " TOO HIGH threshold declared? seems need to take actions" << endl;
         else cout << " Big noise around? Need to wait a little bit " << endl;
       }
     }
     else
     {
       normal_statistic++;
       if( debug )
       {
         cout << " INSPECT GAMMA NOISE in Calorimeter  " << GetName() << " cell  " << GetCellName(i) << " NORMAL " << endl;
         cout << " Stat " << stat_ngam_in_cell << " Stat around   " << stat_ngam_around << " Stat check   " <<
               cells_info[CLUSTER][NEW][i].GetEntries() << " Gate " << 4*(sigma_stat_ngam_in_cell+sigma_stat_ngam_around) << endl;

       }
     }
  }
  if( printout_level >= 0 )
  {
    cout << " **** CALORIMETER " << GetName() << " InspectCellsGamNoise STATISTIC Total=" <<
                                                               cells_gamma_noise_statistic << endl;
  }
  if( printout_level > 0 && cells_gamma_noise_statistic > 0 )
  {
    cout << "-----------------------------------------------------------------------------------------" << endl;
//     cout << "                     I        >50%            10-50%          1-10%             <1%      " << endl;
    cout << " Normal statistic = " << normal_statistic << " ( " <<
                   100.*(double)(normal_statistic)/(double)(NCells()) << " % of Cells )" << endl;

    cout << " Empty statistic = " << empty_statistic << " ( " <<
                   100.*(double)(empty_statistic)/(double)(NCells()) << " % of Cells )" << endl;

    cout << " Low statistic = " << low_statistic << " ( " <<
                   100.*(double)(low_statistic)/(double)(NCells()) << " % of Cells )" << endl;

    cout << " Up  statistic = " << need_up_statistic << " ( " <<
                   100.*(double)(need_up_statistic)/(double)(NCells()) << " % of Cells )" << endl;

    cout << " Down now  statistic = " << need_down_now_statistic << " ( " <<
                   100.*(double)(need_down_now_statistic)/(double)(NCells()) << " % of Cells )" << endl;

    cout << " Down later  statistic = " << need_down_later_statistic << " ( " <<
                   100.*(double)(need_down_later_statistic)/(double)(NCells()) << " % of Cells )" << endl;

    cout << "-----------------------------------------------------------------------------------------" << endl;
    cout << endl;
  }
  return 0;
}

// ////////////////////////////////////////////////////////////////////////////////

int Calorimeter::InspectPedestal(void) const
{
 bool debug = false;
// bool debug = true;
 if( debug  )
   cout << " Calorimeter::InspectPedestals in " << GetName() << " store_noise_info was disabled " << endl;

 if( debug ) cout << " probability array size " << prob_cut_.size() << endl;
 if( debug ) cout << " ebin array size " << ebin_stat_.size() << endl;

 if( prob_cut_.size() == 0 )
 {
   cerr << " probability array not initialized in " << GetName() << endl;
   exit(1);
 }

  if( calo_hist == NULL )
  {
    cerr << " ERROR Calorimeter::InspectPedestals calo_hist == NULL in " << GetName() << endl;
    return -1;
  }
  CalorimeterHist &h=*calo_hist;  // create a short name for calo_hist

  if(  h.h1_PEDCell.size() == 0 )
  {
    cerr << " ERROR Calorimeter::InspectPedestals in " << GetName();
    cerr << " Storage is empty!!! Filling and booking have never been called !!! " << endl;
    return -2;
  }
  if(  h.h1_PEDCell.size() != NCells() )
  {
    cerr << " ERROR Calorimeter::InspectPedestals in " << GetName();
    cerr << " INTERNAL ERROR in  Storage size !!! " << endl;
    return -2;
  }

  if( debug ) cout << " Calorimeter::InspectPedestals in " << GetName() <<
                                                                 " goto cycle over cells " << endl;
  const unsigned npmax = prob_cut_.size();
  if( npmax <= 0 )
  {
    cerr << " Prob array not initialized " << endl;
    exit(1);
  }

  double max_stat = 0.;
  for( size_t i=0; i<NCells(); i++ )
  {
    for( size_t ip=0; ip<= npmax; ip++ )
    {
      energy_cut_[ip][i]=0;
    }
  }

  for( size_t i=0; i<NCells(); i++ )
  {
     if( debug ) cout << " Cell " << i  << "  " << GetCellName(i) << endl;
     double cf = cells_info[CALIB][OLD][i].GetMean();
     if(  h.h1_PEDCell[i] == NULL )
     {
       cerr << " Calorimeter::InspectPedestals in " << GetName() << endl;
       cerr << " h.h1_PEDCell[" << i <<"] == NULL which is a mess!!! " << endl;
       return -2;
     }
     double statistic_in_cell = h.h1_PEDCell[i]->GetEntries();
     double underflow = h.h1_PEDCell[i]->GetBinContent(0);
     double stat_noise_in_cell = statistic_in_cell - underflow;
         TAxis *axis = h.h1_PEDCell[i]->GetXaxis();
         Int_t ncx   = axis->GetNbins();
//         Double_t xmin = axis->GetXmin();
     Double_t maxbin = 0.;
     Int_t binmax=-1;
     for (Int_t bin=1;bin <= ncx+1;bin++)
     {
       if( h.h1_PEDCell[i]->GetBinContent(bin) > maxbin )
       {
     	 binmax = bin;
         maxbin = h.h1_PEDCell[i]->GetBinContent(bin);
       }
     }
     Axis_t ecut0 = h.h1_PEDCell[i]->GetBinCenter(binmax);

     if( statistic_in_cell > max_stat ) max_stat = statistic_in_cell;
       if( debug )
       {
         cout << " Statistic in cell " << GetCellName(i) << " Entries= " << statistic_in_cell <<
                                             " Max bin " << binmax  <<  " Content " << maxbin << endl;
       }
//      double energy_cut[npmax+1];
     for( size_t ip=0; ip< npmax; ip++ )
     {
//        if(  h.h1_PEDCell[i] == NULL )
//        {
//          cerr << " Calorimeter::InspectPedestals in " << GetName() << endl;
//          cerr << " h.h1_PEDCell[" << i <<"] == NULL which is a mess!!! " << endl;
//          return -2;
//        }
//        double statistic_in_cell = h.h1_PEDCell[i]->GetEntries();
//        double underflow = h.h1_PEDCell[i]->GetBinContent(0);
//        double stat_noise_in_cell = statistic_in_cell - underflow;
//        Double_t maxbin = 0.;
//        Int_t binmax=-1;
//        Int_t bin;
//        for (bin=1;bin <= ncx+1;bin++)
//        {
//          if( h.h1_PEDCell[i]->GetBinContent(bin) > maxbin )
// 	 {
//            binmax = bin;
// 	   maxbin = h.h1_PEDCell[i]->GetBinContent(bin);
// 	 }
//        }
//
//        if( statistic_in_cell > max_stat ) max_stat = statistic_in_cell;

       double stat_cut = prob_cut_[ip]*statistic_in_cell;
       if( debug )
       {
           cout << " Above min " <<  stat_noise_in_cell <<  " For prob " << prob_cut_[ip] <<
                                                         " Cut on statistic " << stat_cut << endl;
       }

       if( stat_noise_in_cell == 0 )
       {
         if(ip==0 ) energy_cut_[0][i]=-3.;
         energy_cut_[ip+1][i]=-3.;
       }
       else if( stat_noise_in_cell > 0 )
       {

//          TAxis *axis = h.h1_PEDCell[i]->GetXaxis();
//          Int_t ncx   = axis->GetNbins();
//          Double_t xmin = axis->GetXmin();
//
         if(ip==0 )
         {
           Double_t sum = 0.;
           Int_t bin;
           for (bin=1;bin <= ncx+1;bin++)
           {
             sum += h.h1_PEDCell[i]->GetBinContent(bin);
             if( sum > 0 ) break;
           }

           Axis_t ecut = 0.;

           ecut = h.h1_PEDCell[i]->GetBinLowEdge(bin);
	   if( debug ) cout << " cell " << i << " summ " << sum <<" bin " << bin << " ecut " << ecut << " ecut+ " << ecut+ h.h1_PEDCell[i]->GetBinWidth(bin) << endl;

           if( stat_noise_in_cell > 10 )
             energy_cut_[0][i]=(ecut-ecut0)*cf;
           else
             energy_cut_[0][i]=-1.;

         }

         if( stat_cut > 5 )
         {
           if( stat_noise_in_cell > stat_cut )
           {
             Double_t sum = 0;
             Double_t sum_before = 0;
             Int_t bin;
             for (bin=ncx+1;bin>=1;bin--)
             {
               sum += h.h1_PEDCell[i]->GetBinContent(bin);
               if( sum > stat_cut ) break;
               sum_before=sum;
             }
//             double prob_new=sum/statistic_in_cell;
//             double prob_before=sum_before/statistic_in_cell;
             if(bin < ncx+1) bin += 1;
             Axis_t ecut = h.h1_PEDCell[i]->GetBinCenter(bin);
             energy_cut_[ip+1][i]=(ecut-ecut0)*cf;
           }
           else  // Predefined threshold is higher, so set at default minimum
           {
             energy_cut_[ip+1][i] = -2.;
           }
         }
         else  // Statistic is not enough to define threshold
         {
           energy_cut_[ip+1][i]=-1.;
         }
       }
       else  // no statistic in cell
       {
         if(ip==0 ) energy_cut_[0][i]=-1.;
         energy_cut_[ip+1][i]=-1.;
       }

       if( debug )
       {
         if( ip == 0 )
         {
           if( debug )
           {
             cout << " Ecut Predefined " << energy_cut_[0][i] << endl;
             cout << " Prob " << prob_cut_[ip] << " Ecut " << energy_cut_[ip+1][i] << endl;
           }
         }
         else
         {
           if( debug )
           {
             cout << " Prob " << prob_cut_[ip] << " Ecut " << energy_cut_[ip+1][i] << endl;
           }
         }
       }

       if( statistic_in_cell > max_stat ) max_stat = statistic_in_cell;
     }
  }

  if( debug ) cout << " Prepare statistic printout " << endl;
  cout << " ebin_stat_.size() = ";
  cout << ebin_stat_.size();
  cout << endl;


  if( ebin_stat_.size() == 0 )
  {
    cerr << " statistic bins not initialized " << endl;
    exit(1);
  }

  int nemax = ebin_stat_.size()+3;
  if( debug ) cout << " nemax= " << nemax << " npmax " << npmax << endl;
  vector < vector < int > >stat_all( npmax+1, vector< int >(nemax, 0) );

  for( size_t i=0; i<NCells(); i++ )
  {
    for( size_t ip=0; ip< npmax+1; ip++ )
    {
      if( debug ) cout << " energy_cut_[" << ip <<"][" << i <<"] " << energy_cut_[ip][i] << endl;
      int iep = -100;
      if( fabs(energy_cut_[ip][i] + 3.) < 0.001 )
      {
        iep = nemax-1;
      }
      else if( fabs(energy_cut_[ip][i]  +1.) < 0.001 )
      {
        iep = nemax-2;
      }
      else
      {
        int ie;
        for ( ie = 0; ie <= (int)ebin_stat_.size(); ie++ )
        {
          if( ie == 0 )
	  {
            if( energy_cut_[ip][i] < ebin_stat_[ie] ) break;
	  }
	  else if ( ie == (int)ebin_stat_.size() )
	  {
            if( energy_cut_[ip][i] > ebin_stat_[ie] ) break;
	  }
	  else
	  {
            if( energy_cut_[ip][i] >= ebin_stat_[ie-1] && energy_cut_[ip][i] < ebin_stat_[ie] ) break;
	  }
        }
        iep = ie;
      }
      if( debug ) cout << " iep = " << iep << endl;

      if( iep >= 0 && iep < nemax )
      {
        stat_all[ip][iep]++;
      }
      else
      {
        cerr << " internal error iep = " << iep << endl;
	exit(1);
      }
    }
  }

  int printout_level = 1;
  if( printout_level > 0  )
  {
    cout << "------------------------------------------------------------------------------------------------------" << endl;
    cout << "  **** CALORIMETER " << GetName() << " InspectPedestals  Max stat " << max_stat << endl;
    cout << "------------------------------------------------------------------------------------------------------" << endl;
    cout << "                      I             -1      -2    -2      -3    -3      -4    -4      -5    -5     -6  " << endl;
    cout << "                      I   Emin    10    5*10    10    5*10    10    5*10    10    5*10    10    5*10   " << endl;
    cout << "------------------------------------------------------------------------------------------------------" << endl;
    cout << " Silent cells ??      I";
    printf (" %6d %6d %6d %6d %6d %6d %6d %6d %6d %6d %6d \n",
                           stat_all[0][nemax-1],stat_all[1][nemax-1],stat_all[2][nemax-1],stat_all[3][nemax-1],stat_all[4][nemax-1],
                           stat_all[5][nemax-1],stat_all[6][nemax-1],stat_all[7][nemax-1],stat_all[8][nemax-1],stat_all[9][nemax-1],stat_all[10][nemax-1] );
    cout << " Not enough statistic I";
    printf (" %6d %6d %6d %6d %6d %6d %6d %6d %6d %6d %6d \n",
                           stat_all[0][nemax-2],stat_all[1][nemax-2],stat_all[2][nemax-2],stat_all[3][nemax-2],stat_all[4][nemax-2],
                           stat_all[5][nemax-2],stat_all[6][nemax-2],stat_all[7][nemax-2],stat_all[8][nemax-2],stat_all[9][nemax-2],stat_all[10][nemax-2] );
    for( int ne=0; ne< nemax-2; ne++)
    {
      if( ne == 0 )
      {
	printf (" Eth < %5.2f          I",ebin_stat_[ne]);
        printf (" %6d %6d %6d %6d %6d %6d %6d %6d %6d %6d %6d \n",
                           stat_all[0][ne],stat_all[1][ne],stat_all[2][ne],stat_all[3][ne],stat_all[4][ne],
                           stat_all[5][ne],stat_all[6][ne],stat_all[7][ne],stat_all[8][ne],stat_all[9][ne],stat_all[10][ne] );
      }
      else if(  ne == nemax-3 )
      {
	printf (" Eth > %5.2f          I",ebin_stat_[ne-1]);
        printf (" %6d %6d %6d %6d %6d %6d %6d %6d %6d %6d %6d \n",
                           stat_all[0][ne],stat_all[1][ne],stat_all[2][ne],stat_all[3][ne],stat_all[4][ne],
                           stat_all[5][ne],stat_all[6][ne],stat_all[7][ne],stat_all[8][ne],stat_all[9][ne],stat_all[10][ne] );
      }
      else
      {
	printf (" Eth %5.2f --  %5.2f  I",ebin_stat_[ne-1],ebin_stat_[ne]);
        printf (" %6d %6d %6d %6d %6d %6d %6d %6d %6d %6d %6d \n",
                           stat_all[0][ne],stat_all[1][ne],stat_all[2][ne],stat_all[3][ne],stat_all[4][ne],
                           stat_all[5][ne],stat_all[6][ne],stat_all[7][ne],stat_all[8][ne],stat_all[9][ne],stat_all[10][ne] );
      }
    }
    cout << "------------------------------------------------------------------------------------------------------" << endl;
    cout << endl;
  }
  return 0;
}

////////////////////////////////////////////////////////////////////////////////

int Calorimeter::InspectCellsNoise4RandomTrigger(void) const
{
 bool debug = false;
// bool debug = true;
 if( debug && !options.store_noise_info )
   cout << " Calorimeter::InspectCellsNoise4RandomTrigger in " << GetName() << " store_noise_info was disabled " << endl;
 if (!options.store_noise_info) return 0; // Store information to investigate noise channels
 if (!options.store_rndm_noise_info) return 0; // Store information to investigate noise channels
 if (!options.store_noise_histo)
 {
   cerr << " Please set MANDATORY(for a while) Option store_noise_histo to store and calculate thresholds by " <<
           " Calorimeter::InspectCellsNoise4RandomTrigger in " <<  GetName() << endl;
   return 0;
 }

 if( debug ) cout << " probability array size " << prob_cut_.size() << endl;
 if( debug ) cout << " ebin array size " << ebin_stat_.size() << endl;

 if( prob_cut_.size() == 0 )
 {
   cerr << " probability array not initialized in " << GetName() << endl;
   exit(1);
 }

  if( calo_hist == NULL )
  {
    cerr << " ERROR Calorimeter::InspectCellsNoise4RandomTrigger calo_hist == NULL in " << GetName() << endl;
    return -1;
  }
  CalorimeterHist &h=*calo_hist;  // create a short name for calo_hist

  if(  h.h1_NOISECell.size() == 0 )
  {
    cerr << " ERROR Calorimeter::InspectCellsNoise4RandomTrigger in " << GetName();
    cerr << " Storage is empty!!! Filling and booking have never been called !!! " << endl;
    return -2;
  }
  if(  h.h1_NOISECell.size() != NCells() )
  {
    cerr << " ERROR Calorimeter::InspectCellsNoise4RandomTrigger in " << GetName();
    cerr << " INTERNAL ERROR in  Storage size !!! " << endl;
    return -2;
  }

  if( debug ) cout << " Calorimeter::InspectCellsNoise4RandomTrigger in " << GetName() <<
                                                                 " goto cycle over cells " << endl;
  const unsigned npmax = prob_cut_.size();
  if( npmax <= 0 )
  {
    cerr << " Prob array not initialized " << endl;
    exit(1);
  }

  double max_stat = 0.;

  for( size_t i=0; i<NCells(); i++ )
  {
     if( debug ) cout << " Cell " << i  << "  " << GetCellName(i) << endl;
//      double energy_cut[npmax+1];
     for( size_t ip=0; ip< npmax; ip++ )
     {
       if(  h.h1_NOISECell[i] == NULL )
       {
         cerr << " Calorimeter::InspectCellsNoise4RandomTrigger in " << GetName() << endl;
         cerr << " h.h1_NOISECell[" << i <<"] == NULL which is a mess!!! " << endl;
         return -2;
       }
       double statistic_in_cell = h.h1_NOISECell[i]->GetEntries();
       double underflow = h.h1_NOISECell[i]->GetBinContent(0);
       double stat_noise_in_cell = statistic_in_cell - underflow;

       if( statistic_in_cell > max_stat ) max_stat = statistic_in_cell;

       double stat_cut = prob_cut_[ip]*statistic_in_cell;
       if( debug )
       {
         cout << " Statistic in cell " << GetCellName(i) << " Entries= " << statistic_in_cell <<
               " Above min " <<  stat_noise_in_cell <<  " For prob " << prob_cut_[ip] <<
                                              " Cut on statistic " << stat_cut << endl;
       }

       if( stat_noise_in_cell == 0 )
       {
         if(ip==0 ) energy_cut_[0][i]=-3.;
         energy_cut_[ip+1][i]=-3.;
       }
       else if( stat_noise_in_cell > 0 )
       {

         TAxis *axis = h.h1_NOISECell[i]->GetXaxis();
         Int_t ncx   = axis->GetNbins();
//         Double_t xmin = axis->GetXmin();

         if(ip==0 )
         {
           Double_t sum = 0.;
           Int_t bin;
           for (bin=1;bin <= ncx+1;bin++)
           {
             sum += h.h1_NOISECell[i]->GetBinContent(bin);
             if( sum > 0 ) break;
           }

           Axis_t ecut = 0.;
//            if(bin > 1 )
//            {
//              bin -= 1;
// //           Axis_t ecut = h.h1_NOISECell[i]->GetBinCenter(bin);
//              ecut = h.h1_NOISECell[i]->GetBinLowEdge(bin) + h.h1_NOISECell[i]->GetBinWidth(bin);
//            }
//            else
//            {
//              ecut = h.h1_NOISECell[i]->GetBinLowEdge(bin);
//            }

           ecut = h.h1_NOISECell[i]->GetBinLowEdge(bin);
	   if( debug ) cout << " cell " << i << " summ " << sum <<" bin " << bin << " ecut " << ecut << " ecut+ " << ecut+ h.h1_NOISECell[i]->GetBinWidth(bin) << endl;

           if( stat_noise_in_cell > 10 )
             energy_cut_[0][i]=ecut;
           else
             energy_cut_[0][i]=-1.;

         }

         if( stat_cut > 5 )
         {
           if( stat_noise_in_cell > stat_cut )
           {
             Double_t sum = 0;
             Double_t sum_before = 0;
             Int_t bin;
             for (bin=ncx+1;bin>=1;bin--)
             {
               sum += h.h1_NOISECell[i]->GetBinContent(bin);
               if( sum > stat_cut ) break;
               sum_before=sum;
             }
//             double prob_new=sum/statistic_in_cell;
//             double prob_before=sum_before/statistic_in_cell;
             if(bin < ncx+1) bin += 1;
             Axis_t ecut = h.h1_NOISECell[i]->GetBinCenter(bin);
             energy_cut_[ip+1][i]=ecut;
           }
           else  // Predefined threshold is higher, so set at default minimum
           {
             energy_cut_[ip+1][i] = -2.;
           }
         }
         else  // Statistic is not enough to define threshold
         {
           energy_cut_[ip+1][i]=-1.;
         }
       }
       else  // no statistic in cell
       {
         if(ip==0 ) energy_cut_[0][i]=-1.;
         energy_cut_[ip+1][i]=-1.;
       }

       if( debug )
       {
         if( ip == 0 )
         {
           if( debug )
           {
             cout << " Ecut Predefined " << energy_cut_[0][i] << endl;
             cout << " Prob " << prob_cut_[ip] << " Ecut " << energy_cut_[ip+1][i] << endl;
           }
         }
         else
         {
           if( debug )
           {
             cout << " Prob " << prob_cut_[ip] << " Ecut " << energy_cut_[ip+1][i] << endl;
           }
         }
       }

       if( statistic_in_cell > max_stat ) max_stat = statistic_in_cell;
     }
  }

  if( debug ) cout << " Prepare statistic printout " << endl;
  cout << " ebin_stat_.size() = ";
  cout << ebin_stat_.size();
  cout << endl;


//   const double ebin_stat[] = { 0.1, 0.4, 1., 2.5, 10. };
  if( ebin_stat_.size() == 0 )
  {
    cerr << " statistic bins not initialized " << endl;
    exit(1);
  }

  int nemax = ebin_stat_.size()+3;
  if( debug ) cout << " nemax= " << nemax << " npmax " << npmax << endl;
  vector < vector < int > > stat_all( npmax+1, vector< int >(nemax, 0) );

  for( size_t i=0; i<NCells(); i++ )
  {
    for( size_t ip=0; ip< npmax+1; ip++ )
    {
      if( debug ) cout << " energy_cut_[ip][i] " << energy_cut_[ip][i] << endl;
      int iep = -100;
      if( fabs(energy_cut_[ip][i] + 3.) < 0.001 )
      {
        iep = nemax-1;
      }
      else if( fabs(energy_cut_[ip][i]  +1.) < 0.001 )
      {
        iep = nemax-2;
      }
      else
      {
        int ie;
        for ( ie = 0; ie <= (int)ebin_stat_.size(); ie++ )
        {
          if( ie == 0 )
	  {
            if( energy_cut_[ip][i] < ebin_stat_[ie] ) break;
	  }
	  else if ( ie == (int)ebin_stat_.size() )
	  {
            if( energy_cut_[ip][i] > ebin_stat_[ie] ) break;
	  }
	  else
	  {
            if( energy_cut_[ip][i] >= ebin_stat_[ie-1] && energy_cut_[ip][i] < ebin_stat_[ie] ) break;
	  }
        }
        iep = ie;
      }
      if( debug ) cout << " iep = " << iep << endl;

      if( iep >= 0 && iep < nemax )
      {
        stat_all[ip][iep]++;
      }
      else
      {
        cerr << " internal error iep = " << iep << endl;
	exit(1);
      }
    }
  }

  int printout_level = 1;
  if( printout_level > 0  )
  {
    cout << "------------------------------------------------------------------------------------------------------" << endl;
    cout << "  **** CALORIMETER " << GetName() << " InspectCellsNoise4RandomTrigger  Max stat " << max_stat << endl;
    cout << "------------------------------------------------------------------------------------------------------" << endl;
    cout << "                      I             -1      -2    -2      -3    -3      -4    -4      -5    -5     -6  " << endl;
    cout << "                      I   Emin    10    5*10    10    5*10    10    5*10    10    5*10    10    5*10   " << endl;
    cout << "------------------------------------------------------------------------------------------------------" << endl;
    cout << " Silent cells ??      I";
    printf (" %6d %6d %6d %6d %6d %6d %6d %6d %6d %6d %6d \n",
                           stat_all[0][nemax-1],stat_all[1][nemax-1],stat_all[2][nemax-1],stat_all[3][nemax-1],stat_all[4][nemax-1],
                           stat_all[5][nemax-1],stat_all[6][nemax-1],stat_all[7][nemax-1],stat_all[8][nemax-1],stat_all[9][nemax-1],stat_all[10][nemax-1] );
    cout << " Not enough statistic I";
    printf (" %6d %6d %6d %6d %6d %6d %6d %6d %6d %6d %6d \n",
                           stat_all[0][nemax-2],stat_all[1][nemax-2],stat_all[2][nemax-2],stat_all[3][nemax-2],stat_all[4][nemax-2],
                           stat_all[5][nemax-2],stat_all[6][nemax-2],stat_all[7][nemax-2],stat_all[8][nemax-2],stat_all[9][nemax-2],stat_all[10][nemax-2] );
    for( int ne=0; ne< nemax-2; ne++)
    {
      if( ne == 0 )
      {
	printf (" Eth < %5.2f          I",ebin_stat_[ne]);
        printf (" %6d %6d %6d %6d %6d %6d %6d %6d %6d %6d %6d \n",
                           stat_all[0][ne],stat_all[1][ne],stat_all[2][ne],stat_all[3][ne],stat_all[4][ne],
                           stat_all[5][ne],stat_all[6][ne],stat_all[7][ne],stat_all[8][ne],stat_all[9][ne],stat_all[10][ne] );
      }
      else if(  ne == nemax-3 )
      {
	printf (" Eth > %5.2f          I",ebin_stat_[ne-1]);
        printf (" %6d %6d %6d %6d %6d %6d %6d %6d %6d %6d %6d \n",
                           stat_all[0][ne],stat_all[1][ne],stat_all[2][ne],stat_all[3][ne],stat_all[4][ne],
                           stat_all[5][ne],stat_all[6][ne],stat_all[7][ne],stat_all[8][ne],stat_all[9][ne],stat_all[10][ne] );
      }
      else
      {
	printf (" Eth %5.2f --  %5.2f  I",ebin_stat_[ne-1],ebin_stat_[ne]);
        printf (" %6d %6d %6d %6d %6d %6d %6d %6d %6d %6d %6d \n",
                           stat_all[0][ne],stat_all[1][ne],stat_all[2][ne],stat_all[3][ne],stat_all[4][ne],
                           stat_all[5][ne],stat_all[6][ne],stat_all[7][ne],stat_all[8][ne],stat_all[9][ne],stat_all[10][ne] );
      }
    }
//     cout << " Eth < 0.1            I";
//     printf (" %6d %6d %6d %6d %6d %6d %6d %6d %6d %6d %6d \n",
//                            stat_all[0][0],stat_all[1][0],stat_all[2][0],stat_all[3][0],stat_all[4][0],
//                            stat_all[5][0],stat_all[6][0],stat_all[7][0],stat_all[8][0],stat_all[9][0],stat_all[10][0] );
//     cout << " Eth   0.1 -- 0.4     I";
//     printf (" %6d %6d %6d %6d %6d %6d %6d %6d %6d %6d %6d \n",
//                            stat_all[0][1],stat_all[1][1],stat_all[2][1],stat_all[3][1],stat_all[4][1],
//                            stat_all[5][1],stat_all[6][1],stat_all[7][1],stat_all[8][1],stat_all[9][1],stat_all[10][1] );
//     cout << " Eth   0.4 -- 1.0     I";
//     printf (" %6d %6d %6d %6d %6d %6d %6d %6d %6d %6d %6d \n",
//                            stat_all[0][2],stat_all[1][2],stat_all[2][2],stat_all[3][2],stat_all[4][2],
//                            stat_all[5][2],stat_all[6][2],stat_all[7][2],stat_all[8][2],stat_all[9][2],stat_all[10][2] );
//     cout << " Eth   1.0 -- 2.5     I";
//     printf (" %6d %6d %6d %6d %6d %6d %6d %6d %6d %6d %6d \n",
//                            stat_all[0][3],stat_all[1][3],stat_all[2][3],stat_all[3][3],stat_all[4][3],
//                            stat_all[5][3],stat_all[6][3],stat_all[7][3],stat_all[8][3],stat_all[9][3],stat_all[10][3] );
//     cout << " Eth   2.5 --10.0     I";
//     printf (" %6d %6d %6d %6d %6d %6d %6d %6d %6d %6d %6d \n",
//                            stat_all[0][4],stat_all[1][4],stat_all[2][4],stat_all[3][4],stat_all[4][4],
//                            stat_all[5][4],stat_all[6][4],stat_all[7][4],stat_all[8][4],stat_all[9][4],stat_all[10][4] );
//     cout << " Eth >10.0            I";
//     printf (" %6d %6d %6d %6d %6d %6d %6d %6d %6d %6d %6d \n",
//                            stat_all[0][5],stat_all[1][5],stat_all[2][5],stat_all[3][5],stat_all[4][5],
//                            stat_all[5][5],stat_all[6][5],stat_all[7][5],stat_all[8][5],stat_all[9][5],stat_all[10][5] );
    cout << "------------------------------------------------------------------------------------------------------" << endl;
    cout << endl;
  }
  return 0;
}

////////////////////////////////////////////////////////////////////////////////

void Calorimeter::DebugChecksCellNoise( int icell ) const
{
  if( icell < 0 )
  {
    cerr << " ERROR in Calorimeter::DebugChecksCellNoise icell = " << icell << " in " << GetName() << endl;
    exit(1);
  }

  if( icell >= (int)NCells() )
  {
    cerr << " ERROR in Calorimeter::DebugChecksCellNoise icell = " << icell << " in " << GetName()
                                                             << " but  NCells() = " << NCells() << endl;
    exit(1);
  }

  if( calo_hist == NULL )
  {
    cerr << " Calorimeter::DebugChecksCellNoise in " << GetName() << endl;
    cerr << " ERROR CalorimeterHist calo_hist == NULL  " << endl;
    exit(1);
  }

  CalorimeterHist &h=*calo_hist;  // create a short name for calo_hist
  if(  h.h1_NOISECell[icell] == NULL )
  {
    cerr << " Calorimeter::DebugChecksCellNoise in " << GetName() << endl;
    cerr << " h.h1_NOISECell[" << icell <<"] == NULL which is a mess!!! " << endl;
    exit(1);
  }

}

////////////////////////////////////////////////////////////////////////////////

double Calorimeter::GetCellNoiseStatEnties( int icell ) const
{
  DebugChecksCellNoise( icell );
  CalorimeterHist &h=*calo_hist;  // create a short name for calo_hist
  return h.h1_NOISECell[icell]->GetEntries();
}

////////////////////////////////////////////////////////////////////////////////

double Calorimeter::GetCellNoiseStatUnderflow( int icell ) const
{
  DebugChecksCellNoise( icell );
  CalorimeterHist &h=*calo_hist;  // create a short name for calo_hist
  return h.h1_NOISECell[icell]->GetBinContent(0);
}

////////////////////////////////////////////////////////////////////////////////

double Calorimeter::GetCellNoiseStatOverflow( int icell ) const
{
  DebugChecksCellNoise( icell );
  CalorimeterHist &h=*calo_hist;  // create a short name for calo_hist
  TAxis *axis = h.h1_NOISECell[icell]->GetXaxis();
  Int_t ncx   = axis->GetNbins();
  return h.h1_NOISECell[icell]->GetBinContent(ncx+1);
}

////////////////////////////////////////////////////////////////////////////////

double Calorimeter::GetCellNoiseStatTotal( int icell ) const
{
  DebugChecksCellNoise( icell );
//  CalorimeterHist &h=*calo_hist;  // create a short name for calo_hist
  return ( GetCellNoiseStatEnties(icell) - GetCellNoiseStatUnderflow(icell) );
}

////////////////////////////////////////////////////////////////////////////////

double Calorimeter::GetCellNoiseStat( int icell, double emin, double emax) const
{
  DebugChecksCellNoise( icell );
  CalorimeterHist &h=*calo_hist;  // create a short name for calo_hist

  TAxis *axis = h.h1_NOISECell[icell]->GetXaxis();
  Int_t ncx   = axis->GetNbins();
//  Double_t xmin = axis->GetXmin();
  Double_t sum = 0.;
  Int_t bin;
  for (bin=1;bin <= ncx+1;bin++)
  {
    Double_t eminbin = h.h1_NOISECell[icell]->GetBinLowEdge(bin);
    Double_t ebinwidth = h.h1_NOISECell[icell]->GetBinWidth(bin);
    Double_t emaxbin = eminbin + ebinwidth;
    if( ebinwidth > emax - emin )
    {
      cerr << " WARNING!!  requested energy gate less than bin width " << endl;
    }
    if( emin > eminbin ) continue;
    if( emax <= emaxbin ) continue;
    sum += h.h1_NOISECell[icell]->GetBinContent(bin);
  }
  return  sum;
}

////////////////////////////////////////////////////////////////////////////////

vector<double>  Calorimeter::GetCellNoiseStatAround( int icell, double emin, double emax) const
{

  vector < double > otvet;
  double entries = GetCellNoiseStatEnties(icell);
  if( entries <= 0 )
  {
    otvet.push_back( 0.);
    otvet.push_back(-1.);
    otvet.push_back( 0.);
    otvet.push_back(-1.);
    otvet.push_back( 0.);
    otvet.push_back(-1.);
    otvet.push_back( 0.);
    otvet.push_back(-1.);
    otvet.push_back(-1.);
    otvet.push_back(0.);
    return otvet;
  }
//  double stat = GetCellNoiseStatTotal(icell);
//  double statunderflow = GetCellNoiseStatUnderflow(icell);
//  double statoverflow = GetCellNoiseStatOverflow(icell);

  double dx = cells[icell].GetCellType().GetSizeX();
  double dy = cells[icell].GetCellType().GetSizeY();
  double s_main=dx*dy/100.;			    // goto cm**2
  double stat_base = 0.;
  stat_base = GetCellNoiseStat( icell, emin, emax);
  double sigma_stat_base = -1.;
  if( stat_base > 0. )
  {
    sigma_stat_base = sqrt( fabs(stat_base) );
  }
  otvet.push_back( stat_base );
  otvet.push_back( sigma_stat_base );

  double stat_around = 0.;
  double statsq_around = 0.;
  double n_neighbors = 0.;

  for( vector<size_t>::const_iterator it=cells[icell].GetNeighbors().begin()+1;
  				      it!=cells[icell].GetNeighbors().end(); it++ )
  {
    if( cell_is_bad_[NEW][*it] ) continue;
    if( emin < energy_cut_bad_cells_new[*it] ) continue;

    double ddx = cells[*it].GetCellType().GetSizeX();
    double ddy = cells[*it].GetCellType().GetSizeY();
    double s_neighbor=ddx*ddy/100.;		  // goto cm**2
    double stat_neighbor =  GetCellNoiseStat( *it, emin, emax);

    stat_around += stat_neighbor;

    statsq_around += stat_neighbor * stat_neighbor;

    n_neighbors++;
  }

  if( n_neighbors == 0 )
  {
    otvet.push_back( 0.);
    otvet.push_back(-1.);
    otvet.push_back( 0.);
    otvet.push_back(-1.);
    otvet.push_back(-1.);
    otvet.push_back(0.);
    return otvet;
  }

  stat_around /= n_neighbors;
  double disp_neighbors = -1.;

  if( n_neighbors > 1 )
  {
    statsq_around /= n_neighbors;
    disp_neighbors = (n_neighbors)/(n_neighbors-1)*(statsq_around - stat_around*stat_around);
  }

  double sigma_neighbors = -1;
  if( disp_neighbors > 0 ) sigma_neighbors = sqrt( disp_neighbors );

  otvet.push_back(stat_around);
  otvet.push_back(sigma_neighbors);
  otvet.push_back(n_neighbors);
  return otvet;

}

////////////////////////////////////////////////////////////////////////////////

class CalorimeterStatNoiseInCellAndNeighbors
{
  public:
        ~CalorimeterStatNoiseInCellAndNeighbors            (void) {}

        /// Default constructor
         CalorimeterStatNoiseInCellAndNeighbors  (int ic, const std::vector < double > &v ) : ic_(ic),v_(v)
	                                         { if( v_.size() != 10 )
						   {
						     std::cerr << " CalorimeterStatNoiseInCellAndNeighbors init error " << std::endl;
						     exit(1);
						   }
						 }

    int    GetCellMain            ( void ) { return ic_;}
    double GetStatMain            ( void ) { return v_[0];}
    double GetSigmaStatMain       ( void ) {
                                             if(v_[1] > 0. ) return v_[1];
					     else	     return 1.;
                                           }
    double GetStatMainsN          ( void ) { return v_[2];}
    double GetSigmaStatMainN      ( void ) {
                                             if(v_[3] > 0. ) return v_[3];
					     else	     return 1.;
                                           }
    double GetStatAround          ( void ) { return v_[4];}
    double GetSigmaNeighbors      ( void ) {
                                             if(v_[5] > 0. ) return v_[5];
					     else	     return 1.;
                                           }
    double GetStatAroundN         ( void ) { return v_[6];}
    double GetSigmaNeighborsN     ( void ) {
                                             if(v_[7] > 0. ) return v_[7];
					     else	     return 1.;
                                           }
    double GetSigmaStatNeighborsN ( void ) {
                                             if(v_[8] > 0. ) return v_[8];
					     else	     return 1.;
                                           }
    double GetGoodNeighbors       ( void ) { return v_[9];}
  private:
   int ic_;
   std::vector < double > v_;

};

int Calorimeter::InspectCellsNoiseContinue(void)
{
//   bool debug = true;
  bool debug = false;
  if( GetName() == "WAD" )
  {
     cout << " We skip WAD in InspectCellsNoiseContinue for a while !! " << endl;
     return 0;
  }
  if( debug ) cout << " Debug Calorimeter::InspectCellsNoiseContinue for " << GetName() << endl;
//  CalorimeterHist &h=*calo_hist;  // create a short name for calo_hist
  double ebase = options.base_energy_noise;

  bad_cells_[NEW].clear();
  for( size_t i=0; i<NCells(); i++ )
  {
    cell_is_bad_[NEW][i]=false;
    energy_cut_bad_cells_new[i]=0.;
  }

  bool local_print = false;

  double erange = ebase;

  for( size_t i=0; i<NCells(); i++ )
  {
//     if( i == 2315 || i == 2308 || i == 2299 || i == 1245 || i == 1387 || i == 1458|| i == 642 )
//       local_print = true;
//     else
//       local_print = false;

    if( debug ) cout << " Cell " << i  << "  " << GetCellName(i) << endl;
    if( local_print ) cout << " Cell " << i  << "  " << GetCellName(i) << endl;
    double entries = GetCellNoiseStatEnties(i);
    if( entries <= 0 ) continue;
//    double stat = GetCellNoiseStatTotal(i);
//    double statunderflow = GetCellNoiseStatUnderflow(i);
//    double statoverflow = GetCellNoiseStatOverflow(i);

    vector< double > stat_noise_cell_and_neighbors = GetCellNoiseStatAround( i, erange, 500. );
    CalorimeterStatNoiseInCellAndNeighbors cellstat( i, stat_noise_cell_and_neighbors);

    double stat_base = cellstat.GetStatMain();
//    double sigma_stat_base = cellstat.GetSigmaStatMain();
    double stat_base_n = cellstat.GetStatMainsN();
    double sigma_stat_base_n = cellstat.GetSigmaStatMainN();
    double stat_around = cellstat.GetStatAround();
    double sigma_neighbors = cellstat.GetSigmaNeighbors();
    double stat_around_n = cellstat.GetStatAroundN();
    double sigma_neighbors_n = cellstat.GetSigmaNeighborsN();
    double sigma_stat_neighbors_n = cellstat.GetSigmaStatNeighborsN();

    if( debug ) cout << " stat_base " << stat_base << " stat_around " << stat_around << " sigma " << sigma_neighbors << endl;
    if( debug ) cout << " stat_base_n " << stat_base_n << " sigma_stat_base_n " << sigma_stat_base_n <<
                 " stat_around_n " << stat_around_n << " sigma_n " << sigma_neighbors_n <<  " sigma_stat_n " << sigma_stat_neighbors_n << endl;
    double sc = 1000000./entries;
    if( debug ) cout << " prob_base_n " << stat_base_n*sc << " sigma_prob_base_n " << sigma_stat_base_n*sc <<
                          " prob_around_n " << stat_around_n*sc <<
		            " prob_sigma_n " << sigma_neighbors_n*sc << endl;

//     if(  stat_around > 10. &&
//          stat_base < stat_around - 5.*sigma_neighbors  &&
//          stat_base_n < stat_around_n - 10.*sigma_neighbors_n  &&
// 	 stat_base_n < stat_around_n - 10.*sigma_stat_base_n )
//     if(  stat_base < stat_around/10. && stat_around > 10. &&
//          stat_base_n < stat_around_n - 3.*sigma_neighbors_n  &&
// 	 stat_base_n < stat_around_n - 5.*sigma_stat_neighbors_n )
//     if(  stat_base < stat_around/10. && stat_around > 10. &&
// 	 stat_base_n < stat_around_n - 6.*sigma_stat_neighbors_n )

    if(  stat_base < stat_around/10. && stat_around > 10. )
    {
//       cout << " Dead cell candidate " << i  << "  " << GetCellName(i) << " LED " << cells_info[LED][NEW][i] << endl;
//       cout << " stat_base " << stat_base << " sigma_stat_base " << sigma_stat_base <<
//               " stat_around " << stat_around << " sigma_neighbors " << sigma_neighbors << endl;
//       cout << " stat_base_n " << stat_base_n << " sigma_stat_base_n " << sigma_stat_base_n <<
//              " stat_around_n " << stat_around_n << " sigma_neighbors_n " << sigma_neighbors_n <<
// 	      " sigma_stat_neighbors_n " << sigma_stat_neighbors_n << endl;
      if( cells_info[LED][NEW][i].GetEntries() == 0 ||
         (cells_info[LED][NEW][i].GetEntries() > 1 && cells_info[LED][NEW][i].GetMean() < 100) )
      {
        cell_is_bad_[NEW][i]=true;
        bad_cells_[NEW].push_back(i);
      }
    }
    else
    {
//       if(  i == 2309 )
//       {
//         cout << " NEPOPAL?? V Dead cell candidate " << i  << "  " << GetCellName(i) << " LED " << cells_info[LED][NEW][i] << endl;
//         cout << " stat_base " << stat_base << " sigma_stat_base " << sigma_stat_base <<
//               " stat_around " << stat_around << " sigma_neighbors " << sigma_neighbors << endl;
//         cout << " stat_base_n " << stat_base_n << " sigma_stat_base_n " << sigma_stat_base_n <<
//              " stat_around_n " << stat_around_n << " sigma_neighbors_n " << sigma_neighbors_n <<
// 	      " sigma_stat_neighbors_n " << sigma_stat_neighbors_n << endl;
//         exit(1);
//       }
    }
  }

  int nestep = 100;
//   double nsigmabig = 3.;
//   double nsigmasmall = 5.;
   double nsigmabig = 6.;
   double nsigmasmall = 10.;
//   double nsigmabig = 1.;
//   double nsigmasmall = 1.;
  double estep = 2*ebase/nestep;
  double stat_around_min = 5.;
  if( debug ) cout << " POROGI ************************************************* Ebase " << ebase << endl;
  double eminmax = 10.;   // Maximum threshold we can imagine
  for( int ie=0; ie<nestep; ie++ )
  {
    double emin = estep*(nestep-ie);
    if( debug ) cout << " Go stedily from top to bottom " << ie << " Emin =" << emin << endl;
    if( debug ) cout << " emin " << emin << endl;
    for( size_t i=0; i<NCells(); i++ )
    {
//     if( i == 2315 || i == 2308 || i == 2299 || i == 1245 || i == 1387 || i == 1458|| i == 642 )
//       local_print = true;
//     else
//       local_print = false;
//     if( local_print ) cout << " Cell " << i  << "  " << GetCellName(i) << endl;
      if( energy_cut_bad_cells_new[i] > 0. ) continue;
      if( cell_is_bad_[NEW][i] ) continue;

      vector< double > stat_noise_cell_and_neighbors = GetCellNoiseStatAround( i, emin, 500. );
      CalorimeterStatNoiseInCellAndNeighbors cellstat( i, stat_noise_cell_and_neighbors);
      if( cellstat.GetGoodNeighbors() < 1. ) continue;
      double stat_base = cellstat.GetStatMain();
      double sigma_stat_base = cellstat.GetSigmaStatMain();
      double stat_base_n = cellstat.GetStatMainsN();
      double sigma_stat_base_n = cellstat.GetSigmaStatMainN();
      double stat_around = cellstat.GetStatAround();
      double sigma_neighbors = cellstat.GetSigmaNeighbors();
      double stat_around_n = cellstat.GetStatAroundN();
      double sigma_neighbors_n = cellstat.GetSigmaNeighborsN();
      if( sigma_neighbors_n <= 0 ) sigma_neighbors_n = 1.;
      double sigma_stat_neighbors_n = cellstat.GetSigmaStatNeighborsN();
      if( sigma_stat_neighbors_n <= 0 ) sigma_stat_neighbors_n = 1.;

      double diff0_n = stat_base_n - stat_around_n;
      double diff0_n_sig1 = diff0_n/sigma_neighbors_n;
      double diff0_n_sig2 = diff0_n/sigma_stat_neighbors_n;


      if(  stat_base > 10. && diff0_n_sig1 > nsigmabig && diff0_n_sig2  > nsigmasmall )
      {
        if( local_print )
	{
        cout << " Noisy cell candidate " << i  << "  " << GetCellName(i) << " Emin = " << emin << endl;
	cout << " diff0_n_sig1 " << diff0_n_sig1 << " diff0_n_sig2 " <<  diff0_n_sig2 << endl;
        cout << " stat_base " << stat_base << " sigma_stat_base " << sigma_stat_base <<
              " stat_around " << stat_around << " sigma_neighbors " << sigma_neighbors << endl;
        cout << " stat_base_n " << stat_base_n << " sigma_stat_base_n " << sigma_stat_base_n <<
             " stat_around_n " << stat_around_n << " sigma_neighbors_n " << sigma_neighbors_n <<
	      " sigma_stat_neighbors_n " << sigma_stat_neighbors_n << endl;
        }
        for( int ie1=ie; ie1 >= 0; ie1-- )
        {
          double emin1 = estep*(nestep-ie1);
          vector< double > stat_noise_cell_and_neighbors1 = GetCellNoiseStatAround( i, emin1, 500. );
          CalorimeterStatNoiseInCellAndNeighbors cellstat1( i, stat_noise_cell_and_neighbors1);
          if( local_print ) cout << " emin1 " << emin1 << endl;
          if( emin1 > eminmax )
	  {
	    if( local_print ) cout << " novyi porog " << emin1 << endl;
            energy_cut_bad_cells_new[i] = emin1;
	    break;
          }
          if( ie1 == 0 )
	  {
	    if( local_print ) cout << " novyi porog " << emin1 << endl;
            energy_cut_bad_cells_new[i] = emin1;
	    break;
          }
          double stat_base1 = cellstat1.GetStatMain();
          double stat_base_n1 = cellstat1.GetStatMainsN();
          double stat_around1 = cellstat1.GetStatAround();
          double stat_around_n1 = cellstat1.GetStatAroundN();
          if( stat_base1 <= 1 || stat_around1 <= stat_around_min )
	  {
	    if( local_print ) cout << " low statistic novyi porog po-neobxodimosti " << emin1 << endl;
            energy_cut_bad_cells_new[i] = emin1;
	    break;
          }
          double sigma_stat_base1 = cellstat1.GetSigmaStatMain();
          double sigma_stat_base_n1 = cellstat1.GetSigmaStatMainN();
          double sigma_neighbors1 = cellstat1.GetSigmaNeighbors();
          double sigma_neighbors_n1 = cellstat1.GetSigmaNeighborsN();
          double sigma_stat_neighbors_n1 = cellstat1.GetSigmaStatNeighborsN();

	  double diff_n = stat_base_n1 - stat_around_n1;
	  double diff_n_sig1 = diff_n/sigma_neighbors_n1;
	  double diff_n_sig2 = diff_n/sigma_stat_neighbors_n1;

	  if( local_print ) cout << " stat_around1 " << stat_around1 << " diff_n_sig1 " << diff_n_sig1 << " diff_n_sig2 " << diff_n_sig2 << endl;
          if(  stat_around1 > stat_around_min &&   diff_n_sig1 < 2.)
          {
	    if( local_print )
	    {
	    cout << " novyi porog " << emin1 << endl;
            cout << " stat_base1 " << stat_base1 << " sigma_stat_base1 " << sigma_stat_base1 <<
                  " stat_around1 " << stat_around1 << " sigma_neighbors1 " << sigma_neighbors1 << endl;
            cout << " stat_base_n1 " << stat_base_n1 << " sigma_stat_base_n1 " << sigma_stat_base_n1 <<
                 " stat_around_n1 " << stat_around_n1 << " sigma_neighbors_n1 " << sigma_neighbors_n1 <<
	          " sigma_stat_neighbors_n1 " << sigma_stat_neighbors_n1 << endl;
            }
            energy_cut_bad_cells_new[i] = emin1;
	    break;
          }
        }

      }
    }
  }

  nestep = 100;
  nsigmabig = 3.;
  nsigmasmall = 5.;

  estep = 2*ebase/nestep;
  stat_around_min = 5.;

  if( debug ) cout << " POROGI ************************************************* Ebase " << ebase << endl;
  eminmax = 8.;   // Maximum threshold we can imagine

  for( int ie=0; ie<nestep; ie++ )
  {
    double emin = estep*(nestep-ie);
    if( debug ) cout << " Go stedily from top to bottom " << ie << " Emin =" << emin << endl;
    if( debug ) cout << " emin " << emin << endl;
    for( size_t i=0; i<NCells(); i++ )
    {
//     if( i == 2315 || i == 2308 || i == 2299 || i == 1245 || i == 1387 || i == 1458|| i == 642 )
//       local_print = true;
//     else
//       local_print = false;
    if( local_print ) cout << " Cell " << i  << "  " << GetCellName(i) << endl;
      if( energy_cut_bad_cells_new[i] > 0. ) continue;
      if( cell_is_bad_[NEW][i] ) continue;

      vector< double > stat_noise_cell_and_neighbors = GetCellNoiseStatAround( i, emin, 500. );
      CalorimeterStatNoiseInCellAndNeighbors cellstat( i, stat_noise_cell_and_neighbors);
      if( cellstat.GetGoodNeighbors() < 1. ) continue;
      double stat_base = cellstat.GetStatMain();
      double sigma_stat_base = cellstat.GetSigmaStatMain();
      double stat_base_n = cellstat.GetStatMainsN();
      double sigma_stat_base_n = cellstat.GetSigmaStatMainN();
      double stat_around = cellstat.GetStatAround();
      double sigma_neighbors = cellstat.GetSigmaNeighbors();
      double stat_around_n = cellstat.GetStatAroundN();
      double sigma_neighbors_n = cellstat.GetSigmaNeighborsN();
      if( sigma_neighbors_n <= 0 ) sigma_neighbors_n = 1.;
      double sigma_stat_neighbors_n = cellstat.GetSigmaStatNeighborsN();
      if( sigma_stat_neighbors_n <= 0 ) sigma_stat_neighbors_n = 1.;

      double diff0_n = stat_base_n - stat_around_n;
      double diff0_n_sig1 = diff0_n/sigma_neighbors_n;
      double diff0_n_sig2 = diff0_n/sigma_stat_neighbors_n;


      if(  stat_base > 10. && diff0_n_sig1 > nsigmabig && diff0_n_sig2  > nsigmasmall )
      {
        if( local_print )
	{
        cout << " Noisy cell candidate " << i  << "  " << GetCellName(i) << " Emin = " << emin << endl;
	cout << " diff0_n_sig1 " << diff0_n_sig1 << " diff0_n_sig2 " <<  diff0_n_sig2 << endl;
        cout << " stat_base " << stat_base << " sigma_stat_base " << sigma_stat_base <<
              " stat_around " << stat_around << " sigma_neighbors " << sigma_neighbors << endl;
        cout << " stat_base_n " << stat_base_n << " sigma_stat_base_n " << sigma_stat_base_n <<
             " stat_around_n " << stat_around_n << " sigma_neighbors_n " << sigma_neighbors_n <<
	      " sigma_stat_neighbors_n " << sigma_stat_neighbors_n << endl;
        }
        for( int ie1=ie; ie1 >= 0; ie1-- )
        {
          double emin1 = estep*(nestep-ie1);
          vector< double > stat_noise_cell_and_neighbors1 = GetCellNoiseStatAround( i, emin1, 500. );
          CalorimeterStatNoiseInCellAndNeighbors cellstat1( i, stat_noise_cell_and_neighbors1);
          if( local_print ) cout << " emin1 " << emin1 << endl;
          if( emin1 > eminmax )
	  {
	    if( local_print ) cout << " novyi porog " << emin1 << endl;
            energy_cut_bad_cells_new[i] = emin1;
	    break;
          }
          if( ie1 == 0 )
	  {
	    if( local_print ) cout << " novyi porog " << emin1 << endl;
            energy_cut_bad_cells_new[i] = emin1;
	    break;
          }
          double stat_base1 = cellstat1.GetStatMain();
          double stat_base_n1 = cellstat1.GetStatMainsN();
          double stat_around1 = cellstat1.GetStatAround();
          double stat_around_n1 = cellstat1.GetStatAroundN();
          if( stat_base1 <= 1 || stat_around1 <= stat_around_min )
	  {
	    if( local_print ) cout << " low statistic novyi porog po-neobxodimosti " << emin1 << endl;
            energy_cut_bad_cells_new[i] = emin1;
	    break;
          }
          double sigma_stat_base1 = cellstat1.GetSigmaStatMain();
          double sigma_stat_base_n1 = cellstat1.GetSigmaStatMainN();
          double sigma_neighbors1 = cellstat1.GetSigmaNeighbors();
          double sigma_neighbors_n1 = cellstat1.GetSigmaNeighborsN();
          double sigma_stat_neighbors_n1 = cellstat1.GetSigmaStatNeighborsN();

	  double diff_n = stat_base_n1 - stat_around_n1;
	  double diff_n_sig1 = diff_n/sigma_neighbors_n1;
	  double diff_n_sig2 = diff_n/sigma_stat_neighbors_n1;

	  if( local_print ) cout << " stat_around1 " << stat_around1 << " diff_n_sig1 " << diff_n_sig1 << " diff_n_sig2 " << diff_n_sig2 << endl;
          if(  stat_around1 > stat_around_min &&   diff_n_sig1 < 2. )
          {
	    if( local_print )
	    {
	    cout << " novyi porog " << emin1 << endl;
            cout << " stat_base1 " << stat_base1 << " sigma_stat_base1 " << sigma_stat_base1 <<
                  " stat_around1 " << stat_around1 << " sigma_neighbors1 " << sigma_neighbors1 << endl;
            cout << " stat_base_n1 " << stat_base_n1 << " sigma_stat_base_n1 " << sigma_stat_base_n1 <<
                 " stat_around_n1 " << stat_around_n1 << " sigma_neighbors_n1 " << sigma_neighbors_n1 <<
	          " sigma_stat_neighbors_n1 " << sigma_stat_neighbors_n1 << endl;
	    }
            energy_cut_bad_cells_new[i] = emin1;
	    break;
          }
        }

      }
    }
  }

  double ecut_too_high = 1.5;
  double ecut_min = 0.15;
  for( size_t i=0; i<NCells(); i++ )
  {
//     if( i == 2315 || i == 2308 || i == 2299 || i == 1245 || i == 1387 || i == 1458|| i == 642 )
//       local_print = true;
//     else
//       local_print = false;
    if( local_print ) cout << " Cell " << i  << "  " << GetCellName(i) << endl;

    if( cell_is_bad_[NEW][i] ) continue;
    if( energy_cut_bad_cells_new[i] > ecut_too_high )
    {
      double emin = ecut_min;
      for( vector<size_t>::const_iterator it=cells[i].GetNeighbors().begin()+1;
  				          it!=cells[i].GetNeighbors().end(); it++ )
      {
        if( cell_is_bad_[NEW][*it] ) continue;
        if( energy_cut_bad_cells_new[*it] > ecut_too_high ) continue;
        if( emin < energy_cut_bad_cells_new[*it] ) emin =  energy_cut_bad_cells_new[*it];
      }
      if( emin < ecut_min ) emin = 0.15;
      if( emin > ecut_too_high ) emin = ecut_too_high;

      if( local_print ) cout << " emin " << emin << endl;

      CalorimeterStatNoiseInCellAndNeighbors cellstat( i,GetCellNoiseStatAround( i, emin, 500. ) );
//      double stat_base_n = cellstat.GetStatMainsN();
      double stat_around_n = cellstat.GetStatAroundN();
      double em = emin;
      double dem = 0.02;
      while( em < 8. )
      {
        em += dem;
        CalorimeterStatNoiseInCellAndNeighbors cellstatnew( i,GetCellNoiseStatAround( i, em, 500. ) );
	double stat_base_n_new = cellstatnew.GetStatMainsN();
	if( local_print ) cout << " Povyshaem em=" << em << " stat " << stat_base_n_new << " around " << stat_around_n << endl;
	if( stat_base_n_new < stat_around_n ) break;
      }
      if( local_print ) cout << "  YEEEEE!!!   !!!!  Staryi porog byl " << energy_cut_bad_cells_new[i] << " My podumali i on stal Enew " << em << " YEEEE!!" << endl;
      energy_cut_bad_cells_new[i] = em;
    }
  }
  return 0;

}

////////////////////////////////////////////////////////////////////////////////

int Calorimeter::InspectCellsNoise(void)
{
 bool debug = false;
//  bool debug = true;
 if( debug && !options.store_noise_info )
   cout << " Calorimeter::InspectCellsNoise in " << GetName() << " store_noise_info was disabled " << endl;
 if (!options.store_noise_info) return 0; // Store information to investigate noise channels
 if (options.store_rndm_noise_info) return 0; // Store information to investigate noise channels
 if (!options.store_noise_histo)
 {
   cerr << " Please set MANDATORY(for a while) Option store_noise_histo to store and calculate thresholds by " <<
           " Calorimeter::InspectCellsNoise in " <<  GetName() << endl;
   return 0;
 }

 if( debug ) cout << " probability array size " << prob_cut_.size() << endl;
 if( debug ) cout << " ebin array size " << ebin_stat_.size() << endl;

 if( prob_cut_.size() == 0 )
 {
   cerr << " probability array not initialized in " << GetName() << endl;
   exit(1);
 }

  if( calo_hist == NULL )
  {
    cerr << " ERROR Calorimeter::InspectCellsNoise calo_hist == NULL in " << GetName() << endl;
    return -1;
  }
  CalorimeterHist &h=*calo_hist;  // create a short name for calo_hist

  if(  h.h1_NOISECell.size() == 0 )
  {
    cerr << " ERROR Calorimeter::InspectCellsNoise in " << GetName();
    cerr << " Storage is empty!!! Filling and booking have never been called !!! " << endl;
    return -2;
  }
  if(  h.h1_NOISECell.size() != NCells() )
  {
    cerr << " ERROR Calorimeter::InspectCellsNoise in " << GetName();
    cerr << " INTERNAL ERROR in  Storage size !!! " << endl;
    return -2;
  }

  if( debug ) cout << " Calorimeter::InspectCellsNoise in " << GetName() <<
                                                                 " goto cycle over cells " << endl;
  const unsigned npmax = prob_cut_.size();
  if( npmax <= 0 )
  {
    cerr << " Prob array not initialized " << endl;
    exit(1);
  }

  double max_stat = 0.;

  for( size_t i=0; i<NCells(); i++ )
  {
     if( debug ) cout << " Cell " << i  << "  " << GetCellName(i) << endl;
//      double energy_cut[npmax+1];
     for( size_t ip=0; ip< npmax; ip++ )
     {
       if(  h.h1_NOISECell[i] == NULL )
       {
         cerr << " Calorimeter::InspectCellsNoise in " << GetName() << endl;
         cerr << " h.h1_NOISECell[" << i <<"] == NULL which is a mess!!! " << endl;
         return -2;
       }
       double statistic_in_cell = h.h1_NOISECell[i]->GetEntries();
       double underflow = h.h1_NOISECell[i]->GetBinContent(0);
       double stat_noise_in_cell = statistic_in_cell - underflow;

       if( statistic_in_cell > max_stat ) max_stat = statistic_in_cell;

       double stat_cut = prob_cut_[ip]*statistic_in_cell;
       if( debug )
       {
         cout << " Statistic in cell " << GetCellName(i) << " Entries= " << statistic_in_cell <<
               " Above min " <<  stat_noise_in_cell <<  " For prob " << prob_cut_[ip] <<
                                              " Cut on statistic " << stat_cut << endl;
       }

       if( stat_noise_in_cell == 0 )
       {
         if(ip==0 ) energy_cut_[0][i]=-3.;
         energy_cut_[ip+1][i]=-3.;
       }
       else if( stat_noise_in_cell > 0 )
       {

         TAxis *axis = h.h1_NOISECell[i]->GetXaxis();
         Int_t ncx   = axis->GetNbins();
//         Double_t xmin = axis->GetXmin();

         if(ip==0 )
         {
           Double_t sum = 0.;
           Int_t bin;
           for (bin=1;bin <= ncx+1;bin++)
           {
             sum += h.h1_NOISECell[i]->GetBinContent(bin);
             if( sum > 0 ) break;
           }

           Axis_t ecut = 0.;
//            if(bin > 1 )
//            {
//              bin -= 1;
// //           Axis_t ecut = h.h1_NOISECell[i]->GetBinCenter(bin);
//              ecut = h.h1_NOISECell[i]->GetBinLowEdge(bin) + h.h1_NOISECell[i]->GetBinWidth(bin);
//            }
//            else
//            {
//              ecut = h.h1_NOISECell[i]->GetBinLowEdge(bin);
//            }

           ecut = h.h1_NOISECell[i]->GetBinLowEdge(bin);
	   if( debug ) cout << " cell " << i << " summ " << sum <<" bin " << bin << " ecut " << ecut << " ecut+ " << ecut+ h.h1_NOISECell[i]->GetBinWidth(bin) << endl;

           if( stat_noise_in_cell > 10 )
             energy_cut_[0][i]=ecut;
           else
             energy_cut_[0][i]=-1.;

         }

         if( stat_cut > 5 )
         {
           if( stat_noise_in_cell > stat_cut )
           {
             Double_t sum = 0;
             Double_t sum_before = 0;
             Int_t bin;
             for (bin=ncx+1;bin>=1;bin--)
             {
               sum += h.h1_NOISECell[i]->GetBinContent(bin);
               if( sum > stat_cut ) break;
               sum_before=sum;
             }
//             double prob_new=sum/statistic_in_cell;
//             double prob_before=sum_before/statistic_in_cell;
             if(bin < ncx+1) bin += 1;
             Axis_t ecut = h.h1_NOISECell[i]->GetBinCenter(bin);
             energy_cut_[ip+1][i]=ecut;
           }
           else  // Predefined threshold is higher, so set at default minimum
           {
             energy_cut_[ip+1][i] = -2.;
           }
         }
         else  // Statistic is not enough to define threshold
         {
           energy_cut_[ip+1][i]=-1.;
         }
       }
       else  // no statistic in cell
       {
         if(ip==0 ) energy_cut_[0][i]=-1.;
         energy_cut_[ip+1][i]=-1.;
       }

       if( debug )
       {
         if( ip == 0 )
         {
           if( debug )
           {
             cout << " Ecut Predefined " << energy_cut_[0][i] << endl;
             cout << " Prob " << prob_cut_[ip] << " Ecut " << energy_cut_[ip+1][i] << endl;
           }
         }
         else
         {
           if( debug )
           {
             cout << " Prob " << prob_cut_[ip] << " Ecut " << energy_cut_[ip+1][i] << endl;
           }
         }
       }

       if( statistic_in_cell > max_stat ) max_stat = statistic_in_cell;
     }
  }

  if( debug ) cout << " Prepare statistic printout " << endl;
  cout << " ebin_stat_.size() = ";
  cout << ebin_stat_.size();
  cout << endl;


//   const double ebin_stat[] = { 0.1, 0.4, 1., 2.5, 10. };
  if( ebin_stat_.size() == 0 )
  {
    cerr << " statistic bins not initialized " << endl;
    exit(1);
  }

  int nemax = ebin_stat_.size()+3;
  if( debug ) cout << " nemax= " << nemax << " npmax " << npmax << endl;
  vector < vector < int > > stat_all( npmax+1, vector < int >(nemax, 0) );

  for( size_t i=0; i<NCells(); i++ )
  {
    for( size_t ip=0; ip< npmax+1; ip++ )
    {
      if( debug ) cout << " energy_cut_[ip][i] " << energy_cut_[ip][i] << endl;
      int iep = -100;
      if( fabs(energy_cut_[ip][i] + 3.) < 0.001 )
      {
        iep = nemax-1;
      }
      else if( fabs(energy_cut_[ip][i]  +1.) < 0.001 )
      {
        iep = nemax-2;
      }
      else
      {
        int ie;
        for ( ie = 0; ie <= (int)ebin_stat_.size(); ie++ )
        {
          if( ie == 0 )
	  {
            if( energy_cut_[ip][i] < ebin_stat_[ie] ) break;
	  }
	  else if ( ie == (int)ebin_stat_.size() )
	  {
            if( energy_cut_[ip][i] > ebin_stat_[ie] ) break;
	  }
	  else
	  {
            if( energy_cut_[ip][i] >= ebin_stat_[ie-1] && energy_cut_[ip][i] < ebin_stat_[ie] ) break;
	  }
        }
        iep = ie;
      }
      if( debug ) cout << " iep = " << iep << endl;

      if( iep >= 0 && iep < nemax )
      {
        stat_all[ip][iep]++;
      }
      else
      {
        cerr << " internal error iep = " << iep << endl;
	exit(1);
      }
    }
  }

  int printout_level = 1;
  if( printout_level > 0  )
  {
    cout << "------------------------------------------------------------------------------------------------------" << endl;
    cout << "  **** CALORIMETER " << GetName() << " InspectCellsNoise  Max stat " << max_stat << endl;
    cout << "------------------------------------------------------------------------------------------------------" << endl;
    cout << "                      I             -1      -2    -2      -3    -3      -4    -4      -5    -5     -6  " << endl;
    cout << "                      I   Emin    10    5*10    10    5*10    10    5*10    10    5*10    10    5*10   " << endl;
    cout << "------------------------------------------------------------------------------------------------------" << endl;
    cout << " Silent cells ??      I";
    printf (" %6d %6d %6d %6d %6d %6d %6d %6d %6d %6d %6d \n",
                           stat_all[0][nemax-1],stat_all[1][nemax-1],stat_all[2][nemax-1],stat_all[3][nemax-1],stat_all[4][nemax-1],
                           stat_all[5][nemax-1],stat_all[6][nemax-1],stat_all[7][nemax-1],stat_all[8][nemax-1],stat_all[9][nemax-1],stat_all[10][nemax-1] );
    cout << " Not enough statistic I";
    printf (" %6d %6d %6d %6d %6d %6d %6d %6d %6d %6d %6d \n",
                           stat_all[0][nemax-2],stat_all[1][nemax-2],stat_all[2][nemax-2],stat_all[3][nemax-2],stat_all[4][nemax-2],
                           stat_all[5][nemax-2],stat_all[6][nemax-2],stat_all[7][nemax-2],stat_all[8][nemax-2],stat_all[9][nemax-2],stat_all[10][nemax-2] );
    for( int ne=0; ne< nemax-2; ne++)
    {
      if( ne == 0 )
      {
	printf (" Eth < %5.2f          I",ebin_stat_[ne]);
        printf (" %6d %6d %6d %6d %6d %6d %6d %6d %6d %6d %6d \n",
                           stat_all[0][ne],stat_all[1][ne],stat_all[2][ne],stat_all[3][ne],stat_all[4][ne],
                           stat_all[5][ne],stat_all[6][ne],stat_all[7][ne],stat_all[8][ne],stat_all[9][ne],stat_all[10][ne] );
      }
      else if(  ne == nemax-3 )
      {
	printf (" Eth > %5.2f          I",ebin_stat_[ne-1]);
        printf (" %6d %6d %6d %6d %6d %6d %6d %6d %6d %6d %6d \n",
                           stat_all[0][ne],stat_all[1][ne],stat_all[2][ne],stat_all[3][ne],stat_all[4][ne],
                           stat_all[5][ne],stat_all[6][ne],stat_all[7][ne],stat_all[8][ne],stat_all[9][ne],stat_all[10][ne] );
      }
      else
      {
	printf (" Eth %5.2f --  %5.2f  I",ebin_stat_[ne-1],ebin_stat_[ne]);
        printf (" %6d %6d %6d %6d %6d %6d %6d %6d %6d %6d %6d \n",
                           stat_all[0][ne],stat_all[1][ne],stat_all[2][ne],stat_all[3][ne],stat_all[4][ne],
                           stat_all[5][ne],stat_all[6][ne],stat_all[7][ne],stat_all[8][ne],stat_all[9][ne],stat_all[10][ne] );
      }
    }
//     cout << " Eth < 0.1            I";
//     printf (" %6d %6d %6d %6d %6d %6d %6d %6d %6d %6d %6d \n",
//                            stat_all[0][0],stat_all[1][0],stat_all[2][0],stat_all[3][0],stat_all[4][0],
//                            stat_all[5][0],stat_all[6][0],stat_all[7][0],stat_all[8][0],stat_all[9][0],stat_all[10][0] );
//     cout << " Eth   0.1 -- 0.4     I";
//     printf (" %6d %6d %6d %6d %6d %6d %6d %6d %6d %6d %6d \n",
//                            stat_all[0][1],stat_all[1][1],stat_all[2][1],stat_all[3][1],stat_all[4][1],
//                            stat_all[5][1],stat_all[6][1],stat_all[7][1],stat_all[8][1],stat_all[9][1],stat_all[10][1] );
//     cout << " Eth   0.4 -- 1.0     I";
//     printf (" %6d %6d %6d %6d %6d %6d %6d %6d %6d %6d %6d \n",
//                            stat_all[0][2],stat_all[1][2],stat_all[2][2],stat_all[3][2],stat_all[4][2],
//                            stat_all[5][2],stat_all[6][2],stat_all[7][2],stat_all[8][2],stat_all[9][2],stat_all[10][2] );
//     cout << " Eth   1.0 -- 2.5     I";
//     printf (" %6d %6d %6d %6d %6d %6d %6d %6d %6d %6d %6d \n",
//                            stat_all[0][3],stat_all[1][3],stat_all[2][3],stat_all[3][3],stat_all[4][3],
//                            stat_all[5][3],stat_all[6][3],stat_all[7][3],stat_all[8][3],stat_all[9][3],stat_all[10][3] );
//     cout << " Eth   2.5 --10.0     I";
//     printf (" %6d %6d %6d %6d %6d %6d %6d %6d %6d %6d %6d \n",
//                            stat_all[0][4],stat_all[1][4],stat_all[2][4],stat_all[3][4],stat_all[4][4],
//                            stat_all[5][4],stat_all[6][4],stat_all[7][4],stat_all[8][4],stat_all[9][4],stat_all[10][4] );
//     cout << " Eth >10.0            I";
//     printf (" %6d %6d %6d %6d %6d %6d %6d %6d %6d %6d %6d \n",
//                            stat_all[0][5],stat_all[1][5],stat_all[2][5],stat_all[3][5],stat_all[4][5],
//                            stat_all[5][5],stat_all[6][5],stat_all[7][5],stat_all[8][5],stat_all[9][5],stat_all[10][5] );
    cout << "------------------------------------------------------------------------------------------------------" << endl;
    cout << endl;
  }

  InspectCellsNoiseContinue();

  return 0;
}

////////////////////////////////////////////////////////////////////////////////

void Calorimeter::StoreCellsNoise(void)
{
//   bool debug = true;
  bool debug = false;
  if( debug ) cout << " Calorimeter::StoreCellsNoise " << GetName() <<
                    " options.store_noise_info=" << options.store_noise_info << endl;
  if (!options.store_noise_info) return; // Store information to investigate noise channels
  if (options.store_rndm_noise_info) return; //  Nu nu
  if( cells_info[NOISE][NEW].size() != NCells() )
  {
    cerr << " ERROR Calorimeter::StoreCellsNoise " << GetName() <<
                    " cells_info[NOISE][NEW] is not properely initialized " << endl;
    return;
  }
  cells_noise_statistic++;
  for (vector<CellDataRaw>::const_iterator it=signals_.begin(); it!=signals_.end(); it++) {
    double e = it->GetEnergy();
    int icell = it->GetCellIdx();
    if( e <= 0.05 ) continue;
    if( e < energy_cut_bad_cells_old[icell] ) continue;
    cells_info[NOISE][NEW][icell].Add(e);
  }
  if( debug ) cout << " Calorimeter::StoreCellsNoise in " <<  GetName() <<
                                          " options.store_noise_histo=" << options.store_noise_histo << endl;
  if (!options.store_noise_histo) return; // Store more information to investigate noise channels
  BookNoiseHisto();
  if( calo_hist == NULL ) assert(false);
  CalorimeterHist &h=*calo_hist;  // create a short name for calo_hist
  if( h.h1_NOISECell.size() != NCells() )
  {
    cerr << " ERROR Calorimeter::StoreCellsNoise " << GetName() <<
                    " h.h1_NOISECell is not properely initialized " << endl;
    return;
  }

  std::set<size_t> cellsFilled;

  // loop over calib_data_store to get hit cells
  for (vector<CellDataRaw>::const_iterator it=signals_.begin(); it!=signals_.end(); it++) {
    h.h1_NOISECell[it->GetCellIdx()]->Fill(it->GetEnergy());
    cellsFilled.insert(it->GetCellIdx());
  }

  // also fill all the other cells (with 0. in that case)
  for (size_t i=0; i<NCells(); i++)
    if ( cellsFilled.count(i)==0 )
      h.h1_NOISECell[i]->Fill(0.);

  // also add hit information to raw data statistics
  StoreRawDataStatistic();
}

////////////////////////////////////////////////////////////////////////////////

void Calorimeter::StoreCellsNoise4RandomTrigger(void)
{
//   bool debug = true;
  bool debug = false;
  if( debug ) cout << " Calorimeter::StoreCellsNoise4RandomTrigger " << GetName() <<
                    " options.store_noise_info=" << options.store_noise_info << endl;
  if (!options.store_noise_info) return; // Store information to investigate noise channels
  if (!options.store_rndm_noise_info) return; //  Nu nu
  if( cells_info[NOISE][NEW].size() != NCells() )
  {
    cerr << " ERROR Calorimeter::StoreCellsNoise4RandomTrigger " << GetName() <<
                    " cells_info[NOISE][NEW] is not properely initialized " << endl;
    return;
  }
  cells_noise_statistic++;
  for (vector<CellDataRaw>::const_iterator it=signals_.begin(); it!=signals_.end(); it++ ) {
    double e=it->GetEnergy();
    int icell = it->GetCellIdx();
    if( e <= 0.05 ) continue;
    if( e < energy_cut_bad_cells_old[icell] ) continue;
    cells_info[NOISE][NEW][icell].Add(e);
  }
  if( debug ) cout << " Calorimeter::StoreCellsNoise4RandomTrigger in " <<  GetName() <<
                                          " options.store_noise_histo=" << options.store_noise_histo << endl;
  if (!options.store_noise_histo) return; // Store more information to investigate noise channels
  BookNoiseHisto();
  if( calo_hist == NULL ) assert(false);
  CalorimeterHist &h=*calo_hist;  // create a short name for calo_hist
  if( h.h1_NOISECell.size() != NCells() )
  {
    cerr << " ERROR Calorimeter::StoreCellsNoise4RandomTrigger " << GetName() <<
                    " h.h1_NOISECell is not properely initialized " << endl;
    return;
  }

  std::set<size_t> cellsFilled;

  // loop over calib_data_store to get hit cells
  for (vector<CellDataRaw>::const_iterator it=signals_.begin(); it!=signals_.end(); it++) {
    h.h1_NOISECell[it->GetCellIdx()]->Fill(it->GetEnergy());
    cellsFilled.insert(it->GetCellIdx());
  }

  // also fill all the other cells (with 0. in that case)
  for (size_t i=0; i<NCells(); i++)
    if ( cellsFilled.count(i)==0 )
      h.h1_NOISECell[i]->Fill(0.);

  // also add hit information to raw data statistics
  StoreRawDataStatistic();
}

////////////////////////////////////////////////////////////////////////////////

void Calorimeter::StoreGammaNoise(void)
{
//   bool debug = true;
  bool debug = false;
  if( !options.store_gamma_noise_info ) return;
  cells_gamma_noise_statistic++;

  if( reco_particles.size() == 0 ) // No data anyway
    return;

  for( size_t it=0; it<reco_particles.size(); it++ )
  {
    const vector<size_t> &hitedcells =reco_particles[it].GetMainCells();
    for(size_t i=0; i< hitedcells.size(); i++)
    {
      size_t incell=hitedcells[i];
      if( reco_particles[it].GetE() < energy_gamma_cut_bad_cells_old[i] ) continue;
      if( debug ) cout << " Calorimeter::StoreGammaNoise " << GetName() << " new entry at cell " <<
                                      GetCellName(incell) << " E=" << reco_particles[it].GetE() << endl;
      cells_info[CLUSTER][NEW][incell].Add(reco_particles[it].GetE());
    }
  }

  if (!options.store_noise_histo) return; // Store more information to investigate noise channels
  BookNoiseHisto();
  if( calo_hist == NULL ) assert(false);
  CalorimeterHist &h=*calo_hist;  // create a short name for calo_hist
  if( h.h1_NOISEGamCell.size() != NCells() )
  {
    cerr << " ERROR Calorimeter::StoreGammaNoise " << GetName() <<
                    " h.h1_NOISECell is not properely initialized " << endl;
    return;
  }
  for( size_t it=0; it<reco_particles.size(); it++ )
  {
    const vector<size_t> &hitedcells =reco_particles[it].GetMainCells();
    for(size_t i=0; i< hitedcells.size(); i++)
    {
      size_t incell=hitedcells[i];
      if( h.h1_NOISEGamCell[incell] == NULL) assert(false);
      h.h1_NOISEGamCell[incell]->Fill(reco_particles[it].GetE());
    }
  }
}

////////////////////////////////////////////////////////////////////////////////

void Calorimeter::BookNoiseHisto(void)
{
  bool ok=false;
  if( gDirectory==NULL )
  {
    cerr << "Calorimeter::BookNoiseHisto():  ROOT file was not opend.\n"
         << "                                    I do not want to work.\n";
    return;
  }
  TDirectory *dir_save=gDirectory; // Save current ROOT directory.
// create test histos directory sructure

  if( calo_hist==NULL ) calo_hist = new CalorimeterHist(this);
  CalorimeterHist &h=*calo_hist;  // create a short name for calo_hist

  if(!noise_histo_booked)
  {

    char dir_name[111];
    char hist_name[111];
    if( h.root_dir==NULL )
    {
      ok=(gDirectory->cd("/"));
      assert(ok);
      sprintf(dir_name,"Calorimeter_%s",GetName().c_str());
      h.root_dir = myROOT_utils::TDirectory_checked(dir_name);
    }
    ok=(h.root_dir->cd());
    assert(ok);
// Create subdirectory  Noise
    sprintf(dir_name,"%s_CellsNoise",GetName().c_str());
    h.noise_dir = myROOT_utils::TDirectory_checked(dir_name);
    ok=(h.noise_dir->cd());
    assert(ok);
  if( h.h1_NOISEGamCell.size() != NCells() || h.h1_NOISECell.size() != NCells())
  {
    cerr << " ERROR Calorimeter::BookNoiseHisto " << GetName() <<
                    " h.h1_NOISECell is not properely initialized " << endl;
    dir_save->cd();
    return;
  }

    for( size_t i=0; i< NCells(); i++ )
    {
      char name[132];
      sprintf(name,"NOISE%s",GetCellName(i).c_str());
      sprintf(hist_name," NOISE Energy in cell %s ",GetCellName(i).c_str());
      h.h1_NOISECell[i] = myROOT_utils::TH1D_checked(name,hist_name,1000,0.005,10.005);
      char name1[132];
      sprintf(name1,"NOISEGam%s",GetCellName(i).c_str());
      sprintf(hist_name," NOISE gamma Energy in main cell %s ",GetCellName(i).c_str());
      h.h1_NOISEGamCell[i] = myROOT_utils::TH1D_checked(name1,hist_name,200,0,10.);
    }

    noise_histo_booked=true;
    dir_save->cd();
  }
}

////////////////////////////////////////////////////////////////////////////////

const StatInfo &Calorimeter::GetCellsNoiseInfo(size_t i) const
{
  assert( i<NCells() );
  return cells_info[NOISE][NEW][i];
}

////////////////////////////////////////////////////////////////////////////////

void Calorimeter::PrintCellsNoiseInfo(void) const
{
  cout << "  Cells noise information for calorimeter " << GetName() << "\n";
  for( size_t i=0; i<NCells(); i++ )
  {
    const StatInfo &c = GetCellsNoiseInfo(i);
    printf("    cell %6zu:       entries=%8d",i, (int)c.GetEntries());
    if( c.GetEntries()>0 )
    {
      printf("    mean=%8.3e   sigma=%8.3e",c.GetMean(),c.GetSigma() );
      if( c.GetSigma()!=0 && c.GetMean()!=0 )
        printf("     sigma/mean = %7.3f",c.GetSigma()/c.GetMean() );
    }
    cout << endl;
  }
  cout << endl;
}

////////////////////////////////////////////////////////////////////////////////

} // using namespace Reco
