/*
   $Source: /afs/cern.ch/user/g/gobbo/w0/coral/src/calorim/Reco/src/ReconstructionKolosov.cc,v $
   $Date: 2010/06/29 15:47:03 $
   $Revision: 1.4 $
   -------------------------------------------------------------------------

   This file is part of cellular calorimeter reconstruction program.

   Authors:
     Vladimir  Kolosov   ( Vladimir.Kolosov@cern.ch, Kolosov@mx.ihep.su )
     Alexander Zvyagin   ( Alexander.Zviagine@cern.ch, Zvyagin@mx.ihep.su )
     Denis     Murashev  ( Denis.Mourachev@cern.ch, Murashev@sirius.ihep.su )

   Copyright(C): 1999-2001  V.Kolosov, A.Zvyagin, D.Murashev

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
#include <cmath>
#include <iostream>

// --- Internal files ---
#include "ReconstructionKolosov.h"

#include "CalorimeterParticle.h"
#include "CellDataRaw.h"
#include "CellType.h"
#include "Cluster.h"
#include "Exception.h"
#include "ParticleSort.h"

using namespace std;

namespace Reco {

////////////////////////////////////////////////////////////////////////////////

ReconstructionKolosov::ReconstructionKolosov(const Calorimeter* c) : Reconstruction(c) {
    for (size_t i=0; i<c->NCells(); i++) {
        real_amplitudes    .push_back(pair<double,double>(0,0));
        expected_amplitudes.push_back(pair<double,double>(0,0));
        array_tmp_         .push_back(pair<double,double>(0,0));
    }
}

////////////////////////////////////////////////////////////////////////////////

const std::vector<CalorimeterParticle>& ReconstructionKolosov::DoReconstruction(const std::vector<CellDataRaw>& signals) {
    bool debug = GetCalorimeter()->GetOptions().debug_reconstruction;

    if( debug )  {
        cout << GetCalorimeter()->GetName() << " Start Reconstruction " << GetCalorimeter()->GetOptions().do_reconstruction << endl;
        cout << " GetCalorimeter()->GetOptions().reco_use_hint_particles_info " << GetCalorimeter()->GetOptions().reco_use_hint_particles_info << endl;
        cout << " GetCalorimeter()->GetOptions().reco_cluster_search " << GetCalorimeter()->GetOptions().reco_cluster_search << endl;
        cout << " GetCalorimeter()->GetOptions().reco_cluster_search_only " << GetCalorimeter()->GetOptions().reco_cluster_search_only << endl;
    }

    // clear list of reconstructed showers
    reconstructedParticles.clear();

    // Return if reconstruction should not be done
    if( !GetCalorimeter()->GetOptions().do_reconstruction )
        return reconstructedParticles;

    // standard insertion and checks
    InsertData(signals);

    // Do not use existing hints for particle clusters
    if (!GetCalorimeter()->GetOptions().reco_use_hint_particles_info) { // No hint particles
        // Do not search for clusters
        if ( !GetCalorimeter()->GetOptions().reco_cluster_search && !GetCalorimeter()->GetOptions().reco_cluster_search_only) {
            // No clusters at all

            // Clear list of particles
            manyparticles.clear();

            if ( debug )
                cout << " Before FindGamma0 Ngamma = " << manyparticles.size() << endl;

            // Search for gammas
            FindGamma0(GetCalorimeter()->GetOptions().cell_energy_threshold);

            if ( debug )
                cout << " After FindGamma0 Ngamma = " << manyparticles.size() << endl;

            int add_new = 0;
more_particles:
            // Calculate real response of particles (execute fit, cuts)
            ManyParticlesResponse(); // N Particles may changed

            if ( debug )
                cout << " After ManyParticlesResponse Ngamma = " << manyparticles.size() << " add_new " << add_new <<endl;

            if( add_new == 0 )
            {
                if( debug ) cout << " Try CompareResponse Ngamma = " << manyparticles.size() << endl;
                // CompareResponse seems to be able to split cluster
                // ?????? The idea behind it is unclear to me
                if( !CompareResponse() )
                {
                    // This can just run once
                    add_new++;
                    if( debug ) cout << " After CompareResponse Ngamma = " << manyparticles.size() << " add_new " << add_new << endl;
                    // !!!!!!!!!! Use other construction
                    goto more_particles;
                }
            }

            if( debug ) cout << " After CompareResponse Ngamma = " << manyparticles.size() << endl;

            // Copy particles from manyparticles to reconstructedParticles (after clear reconstructedParticles)
            ReturnResult();

            if( debug ) cout << " After ReturnResult Ngamma = " << manyparticles.size() << endl;

        }

        // Search for clusters and reconstruct particles
        else if (GetCalorimeter()->GetOptions().reco_cluster_search)
        { // Search for clusters first

            vector<OneParticleResponse> manyparticles_all;
            vector<Cluster> clusters;

            // Search for clusters
            FindClusters(clusters,
                    GetCalorimeter()->GetOptions().cluster_search_cell_energy_threshold,
                    GetCalorimeter()->GetOptions().cluster_search_cell_ampl_deviation);

            // Loop over clusters
            for( vector<Cluster>::iterator it=clusters.begin(); it!=clusters.end(); it++ )
            {
                //  Reconstruction cycle for each cluster

                // Clear particle list
                manyparticles.clear();

                // Add cluster
                InsertData(it->GetCluster());

                // Search for gamma using high threshold
                if(manyparticles.size()==0) FindGamma0(0.4); //  High threshold
                // If not succesfull use lower threshold
                if(manyparticles.size()==0) FindGamma0(GetCalorimeter()->GetOptions().cell_energy_threshold); //  Go to low threshold

                int add_new = 0;
more_more_particles:

                // Calculate real response of particles (execute fit, cuts)
                ManyParticlesResponse(); // N Particles may changed

                if( add_new == 0 )
                {
                    // CompareResponse seems to be able to split cluster
                    // ?????? The idea behind it is unclear to me
                    if( !CompareResponse() )
                    {
                        // This can just run once
                        add_new++;
                        // !!!!!!!!!! Use other construction
                        goto more_more_particles;
                    }
                }

                // Store reconstruction results in temporary manyparticles_all container
                // Append to local vector
                for( vector<OneParticleResponse>::iterator it1=manyparticles.begin();
                        it1!=manyparticles.end(); it1++ )
                    manyparticles_all.push_back(*it1);

            } // End loop over clusters

            // Local to global
            manyparticles = manyparticles_all;

            // Copy particles from manyparticles to reconstructedParticles (after clear reconstructedParticles)
            ReturnResult();
        }
        else if( GetCalorimeter()->GetOptions().reco_cluster_search_only )
        { // One cluster is one particle

            vector<Cluster> clusters;
            // Serach for clusters
            FindClusters(clusters,GetCalorimeter()->GetOptions().cluster_search_cell_energy_threshold,
                    GetCalorimeter()->GetOptions().cluster_search_cell_ampl_deviation  );
            // Clear lists of particles
            reconstructedParticles.clear();
            manyparticles.clear ();
            // Loop over clusters
            for( vector<Cluster>::iterator it=clusters.begin(); it!=clusters.end(); it++ )
            {
                // If cluster amplitude below threshold continue
                if (it->AmplitudeTotal() < GetCalorimeter()->GetOptions().particle_energy_threshold)
                    continue;

                // Add cluster to reco_particle
                reconstructedParticles.push_back(CalorimeterParticle(GetCalorimeter()->GetOptions().particle_default,1,it->AmplitudeTotal(),
                            it->MeanX(),it->MeanY(),it->MeanZ(),
                            GetCalorimeter(),
                            0,it->VarX(),it->VarY(),it->VarZ(),0,0));
                reconstructedParticles.back().SetHitedCell(it->GetMaxCell());

                vector< pair< size_t,double> > data;
                for (vector<CellDataRaw>::const_iterator it2=it->GetCluster().begin(); it2!=it->GetCluster().end(); it2++) {
                  if( it2->GetEnergy() > 0. )
                    data.push_back( pair<size_t,double>(it2->GetCellIdx(), it2->GetEnergy()) );
                }
                reconstructedParticles.back().SetClusterData(data);

                // set manyparticles vector just for completeness
                manyparticles.push_back(OneParticleResponse(reconstructedParticles.back(),GetCalorimeter(),it->GetCluster()));
            }
        }

    }
    else // Using of hint particles
    {
        // Vector particles already contain some answer

        manyparticles.clear ();

        // Use hint particles to calculate response and add to manyparticles
        for (vector<CalorimeterParticle>::const_iterator it= GetCalorimeter()->GetHintParticles().begin();
             it!= GetCalorimeter()->GetHintParticles().end(); it++ )
            manyparticles.push_back(OneParticleResponse((*it),GetCalorimeter()));

        // Recalculate expected amplitudes
        ClearExpectedAmplitudes();
        for( vector<OneParticleResponse>::iterator it=manyparticles.begin();
                it!=manyparticles.end(); it++ )
        {
            it->CalculateExpectedResponse();
            AddToExpectedAmplitudes(*it);
        }
        // Copy particles from manyparticles to reconstructedParticles (after clear reconstructedParticles)
        ReturnResult();
    }

    // Complete particle information with time (if desired)
    if( GetCalorimeter()->GetOptions().add_time_to_reco_particles )
        AddTimeToRecoParticles(signals);

    // Sort particle for black box check
    if( GetCalorimeter()->GetOptions().sort_particles )
        SortPart( reconstructedParticles );

    return reconstructedParticles;
}

////////////////////////////////////////////////////////////////////////////////

const std::vector<CalorimeterParticle>& ReconstructionKolosov::DoRepeatReconstruction(const std::vector<CellDataRaw>& signals) {
//    bool debug = true;
    bool debug = false;
    if ( debug )
        cout << " RepeatReconstruction for " << GetCalorimeter()->GetName() << endl;

    const int repeat_reconstruction_level = GetCalorimeter()->GetOptions().repeat_reconstruction; // 0 Base level => Repeat from scratch

    if ( repeat_reconstruction_level < 0 )
        return reconstructedParticles;

    if ( debug ) {
        cout << " Try to RepeatReconstruction for " << GetCalorimeter()->GetName() << endl;
        cout << " reconstructedParticles.size() = " << reconstructedParticles.size() << endl;
        cout << " manyparticles.size() = " << manyparticles.size() << endl;
        cout << " repeat_reconstruction_level = " << repeat_reconstruction_level << endl;
	for( size_t ip=0; ip < reconstructedParticles.size(); ip++ )
	{
          cout <<" Old ip " << ip <<" E = " << reconstructedParticles[ip].GetE() <<
	                                " X = " << reconstructedParticles[ip].GetX() <<
	                                  " Y = " << reconstructedParticles[ip].GetY() << endl;
	}
    }

    if ( repeat_reconstruction_level == 0 ) {
        // Is done already?  ClearExpectedAmplitudes();
        ClearRealAmplitudes();

        manyparticles.clear();
        reconstructedParticles.clear();

        reconstructedParticles = DoReconstruction(signals);
	if ( debug ) {
	  for( size_t ip=0; ip < reconstructedParticles.size(); ip++ )
	  {
            cout <<" New ip " << ip <<" E = " << reconstructedParticles[ip].GetE() <<
	                                " X = " << reconstructedParticles[ip].GetX() <<
	                                  " Y = " << reconstructedParticles[ip].GetY() << endl;
	  }
	}
    } else if ( repeat_reconstruction_level == 2 ) {
        cout << GetCalorimeter()->GetName() <<" RepeatReconstruction  repeat_reconstruction_level == 2 seems doesnot work " << endl;
        // TODO: copy Calorimeter::reco_particles to reconstructedParticlse
        //       and fill manyparticles????????
        if( debug ) {
            cout << " Try to apply new thresholds for " << GetCalorimeter()->GetName() << endl;
            cout << " Bylo " << manyparticles.size() << endl;
            for (size_t i=0; i<reconstructedParticles.size(); i++ )
                reconstructedParticles[i].Print();
        }

        ClearRealAmplitudes();
        for( vector<OneParticleResponse>::iterator it=manyparticles.begin();
                it!=manyparticles.end(); it++ )
            AddToRealAmplitudes(*it);

        ManyParticlesResponse(); // N Particles may changed
        ReturnResult();

        if( debug ) {
            cout << " Stalo " << manyparticles.size() << endl;
            for(size_t i=0; i<reconstructedParticles.size(); i++ )
                reconstructedParticles[i].Print();
        }
    } else {
        cerr << " Repeat Reconstruction at level = " << repeat_reconstruction_level
             << " not implemented for " << GetCalorimeter()->GetName() << " Calorimeter " << endl;
    }

    if ( debug )
        cout << " After Reconstruction reconstructedParticles.size() = " << reconstructedParticles.size() << endl;

    return reconstructedParticles;
}

////////////////////////////////////////////////////////////////////////////////

void ReconstructionKolosov::ReadCalibrations() {
}

////////////////////////////////////////////////////////////////////////////////

 /*!
    \brief Calculate the response of all particles
    \callgraph
    \callergraph
    \remarks {
    Used GetCalorimeter()->GetOptions():
    fitting_type

    Hardcoded:
    bool debug : Enable/disable debuging output,
    int max_iter : Maximum number of iterations
    }
 */

void ReconstructionKolosov::ManyParticlesResponse(void) {
//    bool debug = true;
    bool debug = false;

    if ( debug )
        cout << " ManyParticlesResponse in " << GetCalorimeter()->GetName() << endl;

    const unsigned int max_iter = 3;

    // iterations loop
    for (unsigned int iter=0; iter!=max_iter; iter++) {

        // Clear vector of expected amplitudes
        ClearExpectedAmplitudes();

        // Loop over particles
        for (vector<OneParticleResponse>::iterator it=manyparticles.begin();
             it!=manyparticles.end(); it++) {
            // Calculate expected amplitudes
            it->CalculateExpectedResponse();
            // Add expected amplitudes to vector
            AddToExpectedAmplitudes(*it);
        }

        // Calculate real responses
        ManyParticlesShareRealResponse();

        // Loop over particles
        for (vector<OneParticleResponse>::iterator it=manyparticles.begin();
             it!=manyparticles.end(); it++) {
            if ( debug )
                cout << " Before fit E " << it->GetParticle().GetE()
                     << " X= " << it->GetParticle().GetX()
                     << " Y= " << it->GetParticle().GetY() << endl;

            // Fit particle
            Fit(GetCalorimeter()->GetOptions().fitting_type, *it);

            if ( debug )
                cout << " After fit E " << it->GetParticle().GetE()
                     << " X= " << it->GetParticle().GetX()
                     << " Y= " << it->GetParticle().GetY() << endl;

            // Check if particle position in main cell
            if ( !CheckPointInCell(it->GetParticle().GetX(), it->GetParticle().GetY(), it->GetParticle().GetMainCells()[0]) ) {
                if ( debug )
                    cout << " Break CheckPoint !!! Waaayyy !!! Blin !! " << endl;

                // todo : Esli fit ushel iz maina nado pustit' ves' cikl po novoi
                // Change main cell
                it->ChangeMain();

                // Refit
                Fit(GetCalorimeter()->GetOptions().fitting_type, *it);

                // Check if particle position in main cell
                if( !CheckPointInCell(it->GetParticle().GetX(), it->GetParticle().GetY(), it->GetParticle().GetMainCells()[0]) ) {
                    cerr << " Warning!! Calorimeter::ManyParticlesResponse:: Double fail in CheckPointInCell Calorim "
                         << GetCalorimeter()->GetName() << endl;
                    it->GetParticle().Print();

                    // If again not in main cell change main cell but skip fitting (!!!! there is something wrong ????? Mark)
                    it->ChangeMain();

                    it->GetParticle().Print();
                }
            } // End loop over particles
        }

        if ( debug )
            cout << " ManyParticlesResponse iter = " << iter << " OK " << endl;
    } // end iterations loop

    if ( debug ) {
        cout << " Before ApplyGammaEnergyCuts ngamma = " << manyparticles.size() <<  endl;
        int index = -1;
        for (vector<OneParticleResponse>::iterator it=manyparticles.begin();
             it!=manyparticles.end(); it++) {
            index++;
            const CalorimeterParticle &p=it->GetParticle();
            double e = p.GetE();
            cout << " X-check particle " << index << " E = " << e << endl;
        }
    }

    // Apply energy cuts to manyparticles
    ApplyGammaEnergyCuts(); // Number of gammas might be changed

    if ( debug )
        cout << " Before ClearExpectedAmplitudes ngamma = " << manyparticles.size() <<  endl;

    // Rerun this functions for changed list of particles
    // !!!!! improve handling in ApplyGammaEnergyCuts() to avoid this
    ClearExpectedAmplitudes();
    for (vector<OneParticleResponse>::iterator it =manyparticles.begin();
         it!=manyparticles.end(); it++) {
        it->CalculateExpectedResponse();
        AddToExpectedAmplitudes(*it);
    }
    ManyParticlesShareRealResponse();

    if ( debug ) {
        cout << " At the end  ngamma = " << manyparticles.size() <<  endl;
        cout << " ManyParticlesResponse  OK !!! " << endl;
    }
}

////////////////////////////////////////////////////////////////////////////////

// Calculate the real response of all found particles
/*!
  \brief Determine real response of cell
  \callgraph
  \callergraph
  \remarks {
  Hardcoded:
  bool debug : debug output
*/
void ReconstructionKolosov::ManyParticlesShareRealResponse(void) {
    bool debug = false;

    if ( debug )
        cout << " ManyParticlesShareRealResponse " << GetCalorimeter()->GetName() << "**************************************************" << endl;

    vector<size_t> problems;
    vector<size_t> particles_with_problems;

    bool ok_all = true;

    // Loop over particles
    for( size_t it=0; it< manyparticles.size(); it++) {
        if( !(manyparticles[it].CalculateRealResponse(problems, real_amplitudes, expected_amplitudes)) ) {
            // Mark particle to have problem
            particles_with_problems.push_back(it);
            if( debug ) cout << " Particle " << it << " NEED TO RECOVER FOR " << problems.size() << " CELLS " << endl;
            ok_all = false;
        } else {
            if ( debug )
                cout << " Particle " << it << " OK " << endl;
        }
    }

    // Check if there were problems encountered
    if( !ok_all )
        // Try to recover for problems (Call CalculateRealResponse again ??????)
        RecoverEnergyDelivery(particles_with_problems, problems);
}

///////////////////////////////////////////////////////////////////////////////

/// Apply gamma energy cut to manyparticles
// TODO: Improve handling of manyparticle list, i.e. just remove bad particles from vector
void ReconstructionKolosov::ApplyGammaEnergyCuts(void) {
    bool debug = false;

    if ( debug )
        cout << GetCalorimeter()->GetName() << " Start ApplyGammaEnergyCuts ngamma " << manyparticles.size() << endl;

    vector<int> good_particles;

    // Loop over particles
    int index = 0;
    for (vector<OneParticleResponse>::iterator it=manyparticles.begin();
         it!=manyparticles.end(); it++, index++) {
        // Get particle
        const CalorimeterParticle &p=it->GetParticle();

        // Get particle energy and position
        double e = p.GetE();
        double x = p.GetX();
        double y = p.GetY();

        // Get hited cells
        const vector<size_t> &hitedcells =p.GetMainCells(); // Dont forget!! This is Main cell really !!
        // Assert single cell cluster ??????
        assert( hitedcells.size() == 1 );

        // Get global cell index
        size_t incell = hitedcells[0];

        if ( debug )
            cout << " Calorimeter::ApplyGammaEnergyCuts Particle= e" << e << " x=" << x << " y=" << y << endl;

        // Add cells which pass cut to good cells
        if( !GammaEnergyCut( incell, x, y, e) )
            good_particles.push_back(index);
    }

    // Return directly if all particles pass cut
    if( manyparticles.size() == good_particles.size() )
        return;

    // Reset manyparticles to the list with good particles
    vector<OneParticleResponse> manyparticles_above_treshold;
    for (size_t it=0; it<good_particles.size(); it++)
        manyparticles_above_treshold.push_back(manyparticles[good_particles[it]]);
    manyparticles=manyparticles_above_treshold;  // Number of gammas might be changed

    if( debug )
        cout << GetCalorimeter()->GetName() << " Finish ApplyGammaEnergyCuts ngamma " << manyparticles.size() << endl;
}

////////////////////////////////////////////////////////////////////////////////

bool ReconstructionKolosov::CheckPointInCell(const double& x, const double& y, const size_t& cell ) const {
  if( !GetCalorimeter()->CellIsValid(cell) )
      return false;
  if( fabs(x-GetCalorimeter()->GetCells()[cell].GetX()) + GetCalorimeter()->GetOptions().tolerance_for_nearby_cells
      > 2.*GetCalorimeter()->GetCells()[cell].GetCellType().GetSizeX())
      return false;
  if( fabs(y-GetCalorimeter()->GetCells()[cell].GetY()) + GetCalorimeter()->GetOptions().tolerance_for_nearby_cells
      > 2.*GetCalorimeter()->GetCells()[cell].GetCellType().GetSizeY())
      return false;

  return true;
}

////////////////////////////////////////////////////////////////////////////////

void ReconstructionKolosov::ReturnResult(void) {
    reconstructedParticles.clear();
    for(size_t ip = 0; ip < manyparticles.size(); ip++) {
        reconstructedParticles.push_back(manyparticles[ip].GetParticle());

        vector< pair<size_t, double> > data;
        const vector<size_t>& id_cells = manyparticles[ip].GetCells();
        const vector<double>& e_cells  = manyparticles[ip].GetCellsEnergy();
        for (size_t it=0; it<id_cells.size(); it++)
            if ( e_cells[it] > 0. )
                data.push_back( pair<size_t, double>(id_cells[it],e_cells[it]) );
        reconstructedParticles.back().SetClusterData(data);
    }
}

////////////////////////////////////////////////////////////////////////////////

bool ReconstructionKolosov::CompareResponse(void) {
    bool debug = false;

    if ( debug )
        cout << " CompareResponse:: in " << GetCalorimeter()->GetName() << " GetCalorimeter()->GetOptions().add_search = " << GetCalorimeter()->GetOptions().add_search << endl;

    if ( GetCalorimeter()->GetOptions().add_search ) {
        double devmax_verybad =0.;
        int cell_devmax_verybad =-1;
        int bad_gamma = 0;

        // Loop over particles
        for(size_t it=0; it!=manyparticles.size(); it++) { // cycle over reconstructed particles
            // Get cells in cluster
            const std::vector<size_t> &hited_cells = manyparticles[it].GetCells();
            // Get main cell index
            int main_cell = hited_cells[0];
            // Get neightboring bad cells
            std::vector<size_t> bad_cells = GetCalorimeter()->BadCellsAround(main_cell);

            // If there are no bad cells continue
            if( bad_cells.size() > 0 ) continue;  // Skip Chisq test in the vicinity of a bad cell

            double devmax = 0.;
            int cell_devmax = -1;
            double hi2 = 0;
            int ndf = 0;

            // Loop over cells in cluster
            for (size_t i =0; i!=hited_cells.size(); i++) {
                // Get global cell index
                int jcell = hited_cells[i];

                // Get real energy and energydeviation
                double e_real = manyparticles[it].GetCellsEnergy()[i];
                double de_real = manyparticles[it].GetCellsDispEnergy()[i];

                // Get expected energy and its deviation
                double e_expected = manyparticles[it].GetExpectedCellsEnergy()[i];
                double de_expected = manyparticles[it].GetExpectedCellsDispEnergy()[i];

                // Calculation deviation from real energy
                double de = e_real-e_expected;
                // Calclulate weight
                double we = 1./(de_real+de_expected);
                // Calculate weighted deviation
                double deviation = de*de*we;

                // Check whether this amplitude was actually measured by electronics at all!
                if( real_amplitudes[jcell].first < GetCalorimeter()->GetEnergyCutSparseMode()[jcell] )
                {
                    //  Real amplitude is below registration threshold, but waht do we expect to see here?
                    if( expected_amplitudes[jcell].first < 2.*GetCalorimeter()->GetEnergyCutSparseMode()[jcell] )
                    {
                        //   So if expected amplitude is only twice above registration threshold, probably we don't care about such discrepancy
                        continue;
                    }
                    else
                    {
                        //   The expected amplitude is more then twice above registration threshold, there is something to care about !?
                        cerr <<" Warning in " << GetCalorimeter()->GetName() << " we expect some amplitude " << expected_amplitudes[jcell].first <<
                            " in the cell " << GetCalorimeter()->GetCellName(jcell) << " But nothing was measured!! Badcell??? " << endl;
                    }
                }

                // If deviation larger than maximal occured deviation
                // Remeber this cell
                if( deviation > devmax )
                {
                    devmax = deviation;
                    cell_devmax = jcell;
                }

                hi2 += deviation;
                ndf++;

            } // End loop over cells

            // varinat 18% (2->1) 0.3% (1->2)       if( devmax > 75. || hi2-devmax > 30. ) bad_gamma++;
            double hicut_totm = 30./50.*GetCalorimeter()->GetOptions().hicut_add_search;
            if ( // Energy over threshold
                    ( GetCalorimeter()->GetOptions().ecut_add_search < manyparticles[it].GetParticle().GetE() ) &&
                    // Deviation over threshold
                    (devmax >  GetCalorimeter()->GetOptions().hicut_add_search*(1. -  1./hicut_totm*(hi2-devmax)) || hi2-devmax > hicut_totm) )
                // Mark bad
                bad_gamma++;

            if( cell_devmax >= 0 )
            {
                if( devmax > devmax_verybad ) cell_devmax_verybad = cell_devmax;
            }

            if( debug ) cout << " kg = " << it << " hi2 = " << hi2 << endl;

        } // End loop over particles

        if( debug )
        {
            cout << " devmax = " << devmax_verybad << " cell_devmax = " << cell_devmax_verybad << endl;
            cout << " Ereal " << real_amplitudes[cell_devmax_verybad].first << " Eexpected " << expected_amplitudes[cell_devmax_verybad].first << endl;
            cout << " SigEreal " << real_amplitudes[cell_devmax_verybad].second << " SigEexpected " << expected_amplitudes[cell_devmax_verybad].second << endl;
            cout << " ErealN " << real_amplitudes[cell_devmax_verybad].first << " EexpectedN " << expected_amplitudes[cell_devmax_verybad].first << endl;
            cout << " SigErealN " << real_amplitudes[cell_devmax_verybad].second << " SigEexpectedN " << expected_amplitudes[cell_devmax_verybad].second << endl;
        }


        // return true if there is a bad gamma
        if( bad_gamma == 0 )
        {
            return true;
        }
        else // no bad gamma
        {
            // Add gamma to the maximum energy cells
            if( manyparticles.size() == 1 )
            {

                int jmax1 = -1;
                double emax1 = 0.;
                int jmax2 = -1;
                double emax2 = 0.;
                // Get cells in cluster
                const std::vector<size_t> &hited_cells = manyparticles[0].GetCells();
                // Loop over cells
                for( int i =0; i !=(int)hited_cells.size(); i++ )
                {
                    // Get global cell index
                    int jcell = hited_cells[i];
                    // Get real amplitude
                    double e = real_amplitudes[jcell].first;

                    // Get the two cells with maxmimum amplitude by updating values
                    if( e > emax1 )
                    {
                        emax2 = emax1;
                        jmax2 = jmax1;
                        emax1 = e;
                        jmax1 = jcell;
                    }
                    else if( e > emax2 )
                    {
                        emax2 = e;
                        jmax2 = jcell;
                    }
                } // End loop cells

                // Found two cells with maximal amplitude
                if( jmax1 >= 0 && jmax2 >= 0 )
                {

                    // Clear list of particles
                    ClearExpectedAmplitudes();
                    manyparticles.clear(); //!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

                    // Index of cell with maximum amplitude
                    int cell_devmax = jmax1;
                    // Get position accordingly
                    double xxg = GetCalorimeter()->GetCells()[cell_devmax].GetX();
                    double yyg = GetCalorimeter()->GetCells()[cell_devmax].GetY();
                    // Get amplitude again
                    double eg = real_amplitudes[cell_devmax].first;
                    // Set amplitude to minimum value (if smaler)
                    if( eg < 0.05 ) eg = 0.05;

                    // Add particle to list
                    manyparticles.push_back(OneParticleResponse( CalorimeterParticle(GetCalorimeter()->GetOptions().particle_default,1,
                                    eg,
                                    xxg,yyg,GetCalorimeter()->GetCells()[cell_devmax].GetZ(),
                                    GetCalorimeter(),
                                    0,20.,20.,50.,   // errors: E,X,Y,Z
                                    0.,0.)   // angles: X,Y
                                ,GetCalorimeter(),GetCalorimeter()->GetCells()[cell_devmax], real_amplitudes));
                    // Index of cell with amplitude second to max
                    cell_devmax = jmax2;
                    // Get position
                    xxg = GetCalorimeter()->GetCells()[cell_devmax].GetX();
                    yyg = GetCalorimeter()->GetCells()[cell_devmax].GetY();
                    // Get energy
                    eg = real_amplitudes[cell_devmax].first;
                    if( eg < 0.05 ) eg = 0.05;
                    // Add particle to list
                    manyparticles.push_back(OneParticleResponse( CalorimeterParticle(GetCalorimeter()->GetOptions().particle_default,1,
                                    eg,
                                    xxg,yyg,GetCalorimeter()->GetCells()[cell_devmax].GetZ(),
                                    GetCalorimeter(),
                                    0,20.,20.,50.,   // errors: E,X,Y,Z
                                    0.,0.)    // angles: X,Y
                                ,GetCalorimeter(),GetCalorimeter()->GetCells()[cell_devmax], real_amplitudes));
                    return false;
                }
                else
                {
                    // Add particle to list
                    cerr << " Polnyi BARDAK!! INTERNAL ERROR !!!! " << endl;
                    int cell_devmax = cell_devmax_verybad;
                    double xxg = GetCalorimeter()->GetCells()[cell_devmax].GetX();
                    double yyg = GetCalorimeter()->GetCells()[cell_devmax].GetY();
                    double eg = real_amplitudes[cell_devmax].first;
                    if( eg < 0.05 ) eg = 0.05;
                    manyparticles.push_back(OneParticleResponse( CalorimeterParticle(GetCalorimeter()->GetOptions().particle_default,1,
                                    eg,
                                    xxg,yyg,GetCalorimeter()->GetCells()[cell_devmax].GetZ(),
                                    GetCalorimeter(),
                                    0,20.,20.,50.,   // errors: E,X,Y,Z
                                    0.,0.)    // angles: X,Y
                                ,GetCalorimeter(),GetCalorimeter()->GetCells()[cell_devmax], real_amplitudes));
                    return false;
                }
            }
            else  // Add gamma to the worst cell
            {
                // Add particle to list
                int cell_devmax = cell_devmax_verybad;
                double xxg = GetCalorimeter()->GetCells()[cell_devmax].GetX();
                double yyg = GetCalorimeter()->GetCells()[cell_devmax].GetY();
                double eg = real_amplitudes[cell_devmax].first;
                if( eg < 0.05 ) eg = 0.05;
                manyparticles.push_back(OneParticleResponse( CalorimeterParticle(GetCalorimeter()->GetOptions().particle_default,1,
                                eg,
                                xxg,yyg,GetCalorimeter()->GetCells()[cell_devmax].GetZ(),
                                GetCalorimeter(),
                                0,20.,20.,50.,   // errors: E,X,Y,Z
                                0.,0.)    // angles: X,Y
                            ,GetCalorimeter(),GetCalorimeter()->GetCells()[cell_devmax], real_amplitudes));
                return false;
            }
        }
    }
    else // !GetCalorimeter()->GetOptions().add_search
    {
        return true;
    }
}

////////////////////////////////////////////////////////////////////////////////

bool ReconstructionKolosov::GammaEnergyCut(const size_t& incell, const double& x, const double& y, const double& e) const {
    bool debug = false;

    if ( debug )
        cout << " e " << e
             << " GetCalorimeter()->GetOptions().particle_energy_threshold " << GetCalorimeter()->GetOptions().particle_energy_threshold
             << " UserEnergyCut(x,y) " << GetCalorimeter()->UserEnergyCut(x,y) << endl;

    if ( !(e >= GetCalorimeter()->GetOptions().particle_energy_threshold && e >= GetCalorimeter()->UserEnergyCut(x,y)) )
        return true;

    if ( debug )
        cout << " energy_cut_bad_cells_old[incell] " << GetCalorimeter()->GetEnergyCutInCell(incell) << endl;

    if ( !(real_amplitudes[incell].first >= GetCalorimeter()->GetEnergyCutInCell(incell)) )
        return true;

    if ( !(e >= GetCalorimeter()->GetEnergyGammaCutInCell(incell)) )
        return true;

    if ( GetCalorimeter()->CellIsBad(incell, Calorimeter::OLD) && GetCalorimeter()->GetBadCellStatus(incell, Calorimeter::OLD)>0 ) {
        if ( debug )
            cout << " cell is bad !!!!! " << incell  <<" status " << GetCalorimeter()->GetBadCellStatus(incell, Calorimeter::OLD) << endl;

        return true;
    }

    return false;
}

////////////////////////////////////////////////////////////////////////////////

/*!
  \brief Try to solve problems while calculating real amplitudes
  \param mp {Vector of particles with problems}
  \param vcells { Vector of cells with problems}
  \remarks {
  Hardcoded:
  bool debug : Enable/disable debug output
  }
*/
void ReconstructionKolosov::RecoverEnergyDelivery(const vector<size_t>& mp, const vector<size_t>& vcells) {
    bool debug = true;

    if( debug )
        cout << " RecoverEnergyDelivery " << GetCalorimeter()->GetName() << " N part " <<  mp.size() << " cells " << vcells.size() << endl;

    // Get cells of particles to recover (loop over particles)
    for ( size_t im=0; im < mp.size(); im ++ )
    {
        // Get particle
        OneParticleResponse &pr = manyparticles[mp[im]];
        // Get cells of particle
        const std::vector<size_t> &gcells = pr.GetCells();
        // Loop over cells
        for( size_t il=0; il < gcells.size(); il++ )
        {
            // Get global index
            size_t jcell = gcells[il];
            // Reset array tmp
            array_tmp_[jcell] = pair < double, double> (0.,0.);
        }
    }

    // Loop oveer cells to recover
    for ( size_t ic=0; ic < vcells.size(); ic ++ )
    {
        // Reset array tmp
        array_tmp_[vcells[ic]] = pair < double, double> (0.,-1.);
    }

    // Loop over particles and reset energy in cells (????? combine with first loop in function)
    for ( size_t im=0; im < mp.size(); im ++ )
    {
        OneParticleResponse &pr = manyparticles[mp[im]];
        const std::vector<size_t> &gcells = pr.GetCells();
        const std::vector<double >   &egcells = pr.GetCellsEnergy();
        // Loop over cells
        for( size_t il=0; il < gcells.size(); il++ )
        {
            // Set energy
            size_t jcell = gcells[il];
            array_tmp_[jcell].first += egcells[il];
        }
    }

    // Some info output of particles (?????? combine with upper loop)
    for ( size_t im=0; im < mp.size(); im ++ )
    {
        OneParticleResponse &pr = manyparticles[mp[im]];
        const std::vector<size_t> &gcells = pr.GetCells();
        for( size_t il=0; il < gcells.size(); il++ )
        {
            size_t jcell = gcells[il];
            if( array_tmp_[jcell].second < -0.5 )
            {
                cout << " cell " << jcell << " er " << real_amplitudes[jcell].first <<
                    " check er " << array_tmp_[jcell].first  << endl;
            }
        }
    }


}

////////////////////////////////////////////////////////////////////////////////

void ReconstructionKolosov::ClearArray() {
    for (size_t it=0; it<GetCalorimeter()->NCells(); it++) {
        double sigmaReadout = GetCalorimeter()->GetCells()[it].GetCellType().GetReadoutTerm();
        double cf = GetCalorimeter()->GetCellInfo(Calorimeter::CALIB, Calorimeter::OLD, it).GetMean();
        real_amplitudes[it].second = sigmaReadout*sigmaReadout;

        if ( GetCalorimeter()->GetOptions().readout_sparsified_ ) {
            double sdelta =  GetCalorimeter()->GetEnergyCutSparseMode()[it] / 3.;
            double sdig = cf/3.;
            real_amplitudes[it].second += sdelta*sdelta + sdig*sdig;
        }

        real_amplitudes    [it].first = 0.;
        expected_amplitudes[it]       = pair<double,double>(0,0);
    }
}

////////////////////////////////////////////////////////////////////////////////

void ReconstructionKolosov::ClearExpectedAmplitudes() {
    for (size_t it=0; it<GetCalorimeter()->NCells(); it++)
        expected_amplitudes[it] = pair<double,double>(0., 0.);
}

////////////////////////////////////////////////////////////////////////////////

void ReconstructionKolosov::ClearRealAmplitudes() {
    for (size_t it=0; it<GetCalorimeter()->NCells(); it++)
        real_amplitudes[it] = pair<double,double>(0., 0.);
}

////////////////////////////////////////////////////////////////////////////////

void ReconstructionKolosov::InsertBadCellsData() {
    // Loop over bad cells
    for(list<size_t>::const_iterator itbad=GetCalorimeter()->GetBadCells(Calorimeter::OLD).begin();
        itbad!=GetCalorimeter()->GetBadCells(Calorimeter::OLD).end(); itbad++) {
        // Get bad_cell's global index
        size_t icell = (*itbad);

        if (real_amplitudes[icell].first!=0. && GetCalorimeter()->GetOptions().print_bad_cells_info)
            cerr << "ReconstructionKolosov::InsertBadCellsData: Not zero BadCellData in " << GetCalorimeter()->GetName()
                 << ", bad cell amplitude should be 0., but e=" << real_amplitudes[icell].first
                 << ", cell=" << icell << ", name=" << GetCalorimeter()->GetCellName(icell) << endl;

        // Get maximum amplitude in neightboring cells
        double emax = 0.;
        for (vector<size_t>::const_iterator it=GetCalorimeter()->GetCells()[icell].GetNeighbors().begin()+1;
             it!=GetCalorimeter()->GetCells()[icell].GetNeighbors().end(); it++)
            if( real_amplitudes[*it].first > emax )
                emax = real_amplitudes[*it].first;

        // Set amplitude to zero and amplitude error to sq of highest neightbor
        real_amplitudes[icell].first   = 0.;
        real_amplitudes[icell].second += emax*emax;
    }
}

////////////////////////////////////////////////////////////////////////////////

void ReconstructionKolosov::InsertData(const std::vector<CellDataRaw>& signals) {
    // Clear arrays
    ClearArray();

    // Loop over cells
    for (std::vector<CellDataRaw>::const_iterator it=signals.begin(); it!=signals.end(); it++) {
        // Get cell address
        size_t icell = it->GetCellIdx();

        // Get cell amplitude
        double ecell = it->GetEnergy();

        // Check if cell index in range
        if ( icell>=GetCalorimeter()->NCells() )
            throw Exception("ReconstructionKolosov::InsertData: Cell Address %d is out of calorimeter %s with size %d",
                            icell, GetCalorimeter()->GetName().c_str(), GetCalorimeter()->NCells() );

        // Set real amplitude for cell
        real_amplitudes[icell].first = ecell;

        // Set amplitude error
        real_amplitudes[icell].second = GetCalorimeter()->EnergyDispersionInCellReadout(ecell,icell);
    }

    // Set data for bad cells (ampl = 0, error = sq(max_neightbor)
    InsertBadCellsData();
}

////////////////////////////////////////////////////////////////////////////////

void ReconstructionKolosov::AddTimeToRecoParticles(const std::vector<CellDataRaw>& signals) {
    for (vector<CalorimeterParticle>::iterator ip=reconstructedParticles.begin();
         ip!=reconstructedParticles.end(); ip++) {
        const vector<size_t> &hitedcells = ip->GetMainCells();
        if (!hitedcells.empty()) {
            size_t icell = hitedcells[0]; // main cell

            for (vector<CellDataRaw>::const_iterator it=signals.begin(); it!=signals.end(); it++)
                if (it->GetCellIdx()==icell)
                    ip->SetTime(it->GetTime(), it->GetTimeErr());
        }
    }
}

////////////////////////////////////////////////////////////////////////////////

double ReconstructionKolosov::GetRealAmplitudesSum(void) const {
    double e=0.;
    for (size_t i=0; i<GetCalorimeter()->NCells(); i++)
        e += real_amplitudes[i].first;

    return e;
}

////////////////////////////////////////////////////////////////////////////////

double ReconstructionKolosov::GetExpectedAmplitudesSum(void) const {
    double e=0.;
    for (size_t i=0; i<GetCalorimeter()->NCells(); i++)
         e += expected_amplitudes[i].first;

    return e;
}

////////////////////////////////////////////////////////////////////////////////

void ReconstructionKolosov::AddToRealAmplitudes(const OneParticleResponse& partResponse) {
    size_t i =0;
    for (vector<size_t>::const_iterator it=partResponse.GetCells().begin();
         it!=partResponse.GetCells().end(); it++ ) {
        real_amplitudes[*it].first  += partResponse.GetCellsEnergy()[i];
        real_amplitudes[*it].second += partResponse.GetCellsDispEnergy()[i];
        i++;
    }
}

////////////////////////////////////////////////////////////////////////////////

/*!
  \brief Translate local information about expected response to global memory
  \callgraph
  \callergraph
 */
void ReconstructionKolosov::AddToExpectedAmplitudes(const OneParticleResponse& partResponse) {
    size_t i =0;
    for (vector<size_t>::const_iterator it=partResponse.GetCells().begin();
         it!=partResponse.GetCells().end(); it++ ) {
        expected_amplitudes[*it].first  += partResponse.GetExpectedCellsEnergy()[i];
        expected_amplitudes[*it].second += partResponse.GetExpectedCellsDispEnergy()[i];
        i++;
    }
}

#if 0

////////////////////////////////////////////////////////////////////////////////

void ReconstructionKolosov::ExtractData(vector<CellDataRaw> &data, Calorimeter::DataType what,double extraction_energy_threshold) const {
    const vector<pair <double,double> > &d = what==DataReal ? real_amplitudes : expected_amplitudes;
    for (size_t it=0; it<GetCalorimeter()->NCells(); it++)
        if ( GetCalorimeter()->GetCells()[it].IsActive() ) {
            const double energy = d[it].first;
            if ( energy > extraction_energy_threshold ) {
                data.push_back(CellDataRaw(it));
                data.back().SetEnergy(energy);
            }
        }
}

////////////////////////////////////////////////////////////////////////////////

#endif

}

#if 0
// --- Internal files ---
#include "Calorimeter.h"
#include "Cluster.h"
#include "Shower.h"

#include "ParticleSort.h"

using namespace std;

namespace Reco {

////////////////////////////////////////////////////////////////////////////////

void Calorimeter::ExtractData( const vector<CalorimeterParticle> &particles, vector<CellDataRaw> &external_calib_data_store )
{
  SetCalorimeterParticles( particles );
  external_calib_data_store.clear();
  ClearRealAmplitudes();
  for( vector<OneParticleResponse>::iterator it=manyparticles.begin();
                                             it!=manyparticles.end(); it++ )
  {
    it->AddToRealAmplitudes();
  }
  ExtractData(external_calib_data_store);
}

////////////////////////////////////////////////////////////////////////////////

void Calorimeter::TestEnergyDelivery(void)
{
  bool ok = true;
  if(manyparticles.size() == 0 ) return;
  cout << " WARNING!!! Calorimeter::TestEnergyDelivery  " << GetName() <<
                            " Using  real_amplitudes array NG = " << manyparticles.size() <<  endl;

  for( size_t i=0; i < NCells(); i++ )
  {
//    if( array_tmp_[i].first != 0. ) cerr <<  " cell " << i << " is not empty " << endl;
    array_tmp_[i] = pair< double, double> ( 0., 0.);
  }

  for( vector<OneParticleResponse>::iterator it=manyparticles.begin();
  						       it!=manyparticles.end(); it++ )
  {
    const std::vector<size_t> &gcells = it->GetCells();
//     const std::vector<std::pair<double,double> >   &egcells = it->GetCellsEnergy();
    const std::vector<double>    &egcells = it->GetCellsEnergy();
    for( size_t il=0; il < gcells.size(); il++ )
    {
      size_t jcell = gcells[il];
      array_tmp_[jcell].first += egcells[il];
    }
  }

  int ng=0;
  for( vector<OneParticleResponse>::iterator it=manyparticles.begin();
  						       it!=manyparticles.end(); it++ )
  {
    ng++;
    const std::vector<size_t> &gcells = it->GetCells();
//     const std::vector<std::pair<double,double> >   &egcells = it->GetCellsEnergy();
//     const std::vector<std::pair<double,double> >   &etcells = it->GetExpectedCellsEnergy();
    const std::vector<double >   &egcells = it->GetCellsEnergy();
    const std::vector<double >   &etcells = it->GetExpectedCellsEnergy();
    const std::vector<double >   &detcells = it->GetExpectedCellsDispEnergy();
    for( size_t il=0; il < gcells.size(); il++ )
    {
      size_t jcell = gcells[il];
      double de = real_amplitudes[jcell].first - array_tmp_[jcell].first;
      if( fabs(de) > 0.1 )
      {
        cerr << " KG= " << ng << " Problem in cell " << jcell << " Er " << real_amplitudes[jcell].first <<
	                                            " Es " << array_tmp_[jcell].first << endl;
        cerr << " E measured in cell  " << jcell << " E="  << real_amplitudes[jcell].first
                                      << " E wanted " <<  etcells[il]
                                      << " DispE wanted " <<  detcells[il]
                                      << " E obtained " <<  egcells[il]
                                       << " expected from all " <<  expected_amplitudes[jcell].first
                                         << " disp expected from all " <<  expected_amplitudes[jcell].second << endl;
        ok=false;
      }
    }
  }

  if( !ok ) exit(1);
  return;
}

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

void Calorimeter::CheckPrintManyParticlesResponse(void)
{
  for( vector<OneParticleResponse>::iterator it=manyparticles.begin();
  						       it!=manyparticles.end(); it++ )
  {
    it->CheckPrint();
  }

}

////////////////////////////////////////////////////////////////////////////////
vector<size_t > Calorimeter::GetGammasAroundCell( size_t ibcell ) const
{
  vector<size_t > cnt;
  if( manyparticles.size() <= 0 ) return cnt;
  for( vector<size_t>::const_iterator it=cells[ibcell].GetNeighbors().begin(); it!=cells[ibcell].GetNeighbors().end(); it++ )
  {
    for( size_t itg=0; itg != manyparticles.size(); itg++ )  // cycle over reconstructed particles
    {
      if( (*it) != manyparticles[itg].GetCells()[0] ) continue;
      cnt.push_back(itg);
    }
  }
  return cnt;
}
////////////////////////////////////////////////////////////////////////////////

} // namespace Reco

#endif

