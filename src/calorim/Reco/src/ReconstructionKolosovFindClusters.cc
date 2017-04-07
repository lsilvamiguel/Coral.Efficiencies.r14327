/*
   $Source: /afs/cern.ch/user/g/gobbo/w0/coral/src/calorim/Reco/src/ReconstructionKolosovFindClusters.cc,v $
   $Date: 2010/04/08 16:01:02 $
   $Revision: 1.2 $
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

#include "CellDataRaw.h"
#include "Cluster.h"

namespace Reco {

////////////////////////////////////////////////////////////////////////////////

void ReconstructionKolosov::AddCellToCluster(const size_t& cell_number, Cluster& cluster, std::vector<size_t>& in,
                                             const double& cell_porog, const double& Nstand) const {
    // Do we need to add the cell?
    if ( GetCalorimeter()->GetCells()[cell_number].IsActive()==false )
        return; // Cell is passive.

    if ( real_amplitudes[cell_number].first<= cell_porog ) // Does cell have an amplitude?
        return; // Amplitude too small.

    if ( real_amplitudes[cell_number].first<=sqrt(real_amplitudes[cell_number].second)*Nstand ) // Can we belive in amplitude?
        return; // Amplitude too small.

    if ( in.end()!=std::find(in.begin(),in.end(),cell_number) ) // Test that no clusters conatin the cell.
        return; // Cell already in cluster!

    // Add cell to cluster
    CellDataRaw signal(cell_number);
    signal.SetEnergy(real_amplitudes[cell_number].first);
    cluster.AddCell(signal);
    in.push_back(cell_number);

    // Cycle on all cell's Neighbors.
    for (std::vector<size_t>::const_iterator it=GetCalorimeter()->GetCells()[cell_number].GetNeighbors().begin()+1;
         it!=GetCalorimeter()->GetCells()[cell_number].GetNeighbors().end(); it++ )
        AddCellToCluster(*it,cluster,in,cell_porog,Nstand);
}

////////////////////////////////////////////////////////////////////////////////

void ReconstructionKolosov::FindClusters(std::vector<Cluster>& clusters,
                                         const double& cell_porog,
                                         const double& Nstand) const {
    clusters.clear();
    std::vector<size_t> in;  // List of cells already in a cluster.

    for (size_t it=0; it!=GetCalorimeter()->NCells(); it++) {
        Cluster cluster(GetCalorimeter());
        AddCellToCluster(it, cluster, in, cell_porog, Nstand);
        if ( cluster.Size()>0 )
            clusters.push_back(cluster);
    }
}

////////////////////////////////////////////////////////////////////////////////

} // namespace Reco

