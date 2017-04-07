/*
   $Source: /afs/cern.ch/user/g/gobbo/w0/coral/src/calorim/Reco/src/ReconstructionKolosovFindGamma.cc,v $
   $Date: 2010/06/24 21:51:05 $
   $Revision: 1.3 $
   -------------------------------------------------------------------------

   This file is part of cellular calorimeter reconstruction program.

   Authors:
     Vladimir  Kolosov   ( Kolosov@mx.ihep.su, Vladimir.Kolosov@cern.ch )
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

#include "ReconstructionKolosov.h"

#include "CellType.h"

namespace Reco {

////////////////////////////////////////////////////////////////////////////////

/*!
    \brief Check if cell with index cell_number is candidate for Gamma
    \param cell_number {Index of the cell to be checked}
    \param Nstand {Used to scale error of cell energy}
    \param CellThreshold {Cell energy threshhold}
    \return true if candidate for Gamma
    \callgraph
    \callergraph
    \remarks {
    Used options:
    tolerance_for_nearby_cells

    Hardcoded:
    bool debug : Enable/disable debuging output
    }
*/
bool ReconstructionKolosov::IsNicePeak0(const size_t& cell_number, const double& Nstand, const double& CellThreshold) const {
    bool debug = false;

    // Is it nice peak to consider it as gamma candidate?

    // Make sure that cell is active
    if( GetCalorimeter()->GetCells()[cell_number].IsActive()==false )
        return false; // Cell is passive

    // Check that cell is over threshold
    if( real_amplitudes[cell_number].first<=CellThreshold ) // Does cell have an amplitude?
        return false; // Amplitude too small.

    if ( debug )
        std::cout << GetCalorimeter()->GetName() << " A = " << real_amplitudes[cell_number].first
                  << " Sigma =" << sqrt(real_amplitudes[cell_number].second)*Nstand << std::endl;


    if( real_amplitudes[cell_number].first<=sqrt(real_amplitudes[cell_number].second)*Nstand ) // Can we believe in an amplitude?
        return false; // Amplitude too small.

    // Cycle on all cell's Neighbors.
    for (std::vector<size_t>::const_iterator it=GetCalorimeter()->GetCells()[cell_number].GetNeighbors().begin()+1;
         it!=GetCalorimeter()->GetCells()[cell_number].GetNeighbors().end(); it++) {
        // Check wether neighboring cell has higher amplitude
        if( real_amplitudes[cell_number].first < real_amplitudes[*it].first)
            return false;

        // Check wether neighboring cell has same amplitude
        if( real_amplitudes[cell_number].first == real_amplitudes[*it].first) {
            // See if it can still be considered center of a gamma cluster
            std::pair<int,int> mask = GetCalorimeter()->GetCells()[cell_number].NeighborMask(GetCalorimeter()->GetCells()[*it], GetCalorimeter()->GetOptions().tolerance_for_nearby_cells);
            if ( mask.first  > 0 )
                return false;
            else if( mask.first == 0 )
                if(mask.second > 0 )
                    return false;
        }
    }

    // If we reach here cell is a good candidate to be the center of a gamma cluster
    return true;
}

////////////////////////////////////////////////////////////////////////////////

 /*!
    \brief Try to find Gamma's
    \param CellThreshold {Cell energy threshhold}
    \callgraph
    \callergraph
    \remarks {
    Used options:
    particle_default

    Hardcoded:
    bool debug : Enable/disable debuging output,
    double Nstand : Used to check IsNicePeak0(...)
    }
  */

void ReconstructionKolosov::FindGamma0(const double& CellThreshold) {
    bool debug = false;

    if ( debug )
        std::cout << GetCalorimeter()->GetName() << " Just begin to find Gammas Threshold " << CellThreshold << std::endl;

    double Nstand = 2.;

    // Get distance between gamma detection and vertex position
    double dZ = GetCalorimeter()->GetPositionZ() - GetCalorimeter()->GetVertexPositionZ();
    // Set a default distance if calculation failed ????
    if ( dZ == 0. )
        dZ=6000.;

    // Loop over cells
    for(size_t it=0; it<GetCalorimeter()->NCells(); it++) {
        // Check if cell is candidate for a gamma cluster
        if( IsNicePeak0(it, Nstand, CellThreshold) == true ) {

            std::pair<double,double> eeg = real_amplitudes[it];

            // Get cell position
            double xxg = GetCalorimeter()->GetCells()[it].GetX();
            double yyg = GetCalorimeter()->GetCells()[it].GetY();
            double zzg = GetCalorimeter()->GetCells()[it].GetZ();

            // Don't we need vertex position as well to calculate angels ??????
            double ax = xxg/dZ;
            double ay = yyg/dZ;

            if ( debug )
                std::cout << " New gamma " << " E=" << eeg.first <<" X=" << xxg << " Y=" << yyg << " cell " << it << std::endl;

            // Get cell size
            double sx = GetCalorimeter()->GetCells()[it].GetCellType().GetSizeX();
            double sy = GetCalorimeter()->GetCells()[it].GetCellType().GetSizeY();
            double sz = GetCalorimeter()->GetCells()[it].GetCellType().GetSizeZ();

            // Genrate particle and calculate response (store in vector)
            manyparticles.push_back(OneParticleResponse(
                                        CalorimeterParticle(GetCalorimeter()->GetOptions().particle_default,
                                                            1, real_amplitudes[it].first,
                                                            xxg, yyg, zzg,
                                                            GetCalorimeter(),
                                                            real_amplitudes[it].first/2., sx/2., sy/2., sz/2.,   // errors: E,X,Y,Z
                                                            ax, ay),           // angles: X,Y
                                        GetCalorimeter(),
                                        GetCalorimeter()->GetCells()[it], real_amplitudes));

            if ( debug )
                std::cout << " X-check New gamma " << " E=" << manyparticles.back().GetParticle().GetE() << std::endl;
        }
    } // End loop over cells

    if ( debug )
        std::cout << GetCalorimeter()->GetName() << " FindGamma0 in total "
                  << manyparticles.size() <<" gamma candidates were found " << std::endl;
}

////////////////////////////////////////////////////////////////////////////////

} // namespace Reco

