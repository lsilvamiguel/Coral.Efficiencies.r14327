/*
   $Source: /afs/cern.ch/user/g/gobbo/w0/coral/src/calorim/Reco/src/Cluster.h,v $
   $Date: 2010/04/08 16:01:02 $
   $Revision: 1.8 $
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

#ifndef RecoCluster___include
#define RecoCluster___include

// --- STL ---
#include <ostream>
#include <vector>

namespace Reco {

////////////////////////////////////////////////////////////////////////////////

class Calorimeter;
class CellDataRaw;

/* \brief Cluster is a group of cells.
*/

class Cluster
{
  public:

    // //////////////////////////////////////////////////////////////
    //
    //                 Constructors and destructor
    //
    // //////////////////////////////////////////////////////////////

    /// Destructor
    virtual             ~Cluster                 (void) {}

    /// Construct Cluster with given properties.
                         Cluster                 (const Calorimeter* c);

    // //////////////////////////////////////////////////////////////
    //
    //                 Operators
    //
    // //////////////////////////////////////////////////////////////

    /// Print Cluster info to output stream
    friend std::ostream& operator <<           (std::ostream &o,const Cluster &c);

    // //////////////////////////////////////////////////////////////
    //
    // Methods
    //
    // //////////////////////////////////////////////////////////////

    void                 AddCell                 (const CellDataRaw &c) {cluster_data.push_back(c);}

    const std::vector<CellDataRaw>& GetCluster   (void) const {return cluster_data;}

          std::vector<CellDataRaw> &GetCluster  (void)       {return cluster_data;}

    /// \return Sum(Ai*Xi)/Sum(Ai)
    double               MeanX                   (void) const;

    /// \return Sum(Ai*Yi)/Sum(Ai)
    double               MeanY                   (void) const;

    /// \return Sum(Ai*Zi)/Sum(Ai)
    double               MeanZ                   (void) const;

    /// \return Sum(Ai*Xi*Xi)/Sum(Ai) - sqr(MeanX())
    double               VarX                    (void) const;

    /// \return Sum(Ai*Yi*Yi)/Sum(Ai) - sqr(MeanY())
    double               VarY                    (void) const;

    /// \return Sum(Ai*Zi*Zi)/Sum(Ai) - sqr(MeanZ())
    double               VarZ                    (void) const;

    /// \return Sum(Ai)
    double               AmplitudeTotal          (void) const;

    /// \return cell number with maximum amplitude
    size_t               GetMaxCell              (void) const;

    /// Clear cluster - remove all cells.
    void                 Clear                   (void) {cluster_data.clear();}

    size_t               Size                    (void) const {return cluster_data.size();}

    /// \return test whether all the blocks used for this cluster have the same size and are positioned regularly
    // epsilon allows for some imprecision
    bool                 IsXYRegular             (const double epsilon) const;

  private:

    // //////////////////////////////////////////////////////////////
    //
    //                 Data Members
    //
    // //////////////////////////////////////////////////////////////

    /// Every cluster refer to definite calorimeter
    const Calorimeter*       calorimeter;

    /// Cells of the cluster
    std::vector<CellDataRaw> cluster_data;

};

////////////////////////////////////////////////////////////////////////////////

} /// using namespace Reco
#endif //RecoCluster___include
