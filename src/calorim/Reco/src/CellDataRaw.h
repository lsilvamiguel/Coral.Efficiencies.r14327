/*
   $Source: /afs/cern.ch/user/g/gobbo/w0/coral/src/calorim/Reco/src/CellDataRaw.h,v $
   $Date: 2010/03/31 14:50:08 $
   $Revision: 1.2 $
   -------------------------------------------------------------------------

   This file is part of cellular calorimeter reconstruction program.

   Authors:
     Vladimir  Kolosov   ( Kolosov@mx.ihep.su )
     Alexander Zvyagin   ( Alexander.Zviagine@cern.ch, Zvyagin@mx.ihep.su )
     Denis     Murashev  ( Denis.Mourachev@cern.ch, Murashev@sirius.ihep.su )

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

#ifndef RecoCellDataRaw___include
#define RecoCellDataRaw___include

#include <cstddef>
#include <vector>

namespace Reco {

////////////////////////////////////////////////////////////////////////////////

class CellDataRaw {
    // =========================================================================
    // Constructors and destructor
    // =========================================================================

    public:

        /// Default constructor
                     CellDataRaw    (const size_t& cellIdx);

        /// Default constructor (also setting energy and time)
                     CellDataRaw    (const size_t& cellIdx,
                                     const double& energy,
                                     const double& time);

        /// Copy constructor for changing the cell index
                     CellDataRaw    (const size_t& cellIdx,
                                     const CellDataRaw& signal);

        /// Destructor
        virtual     ~CellDataRaw    () {}

    // ==========================================
    // Methods
    // ==========================================

    public:

        const size_t& GetCellIdx    ()                          const;

        const double& GetAmplitude  ()                          const;

        const double& GetEnergy     ()                          const;
        const double& GetEnergyErr  ()                          const;

        const double& GetTime       ()                          const;
        const double& GetTimeErr    ()                          const;


        const bool&   HasAmplitude  ()                          const;

        const bool&   HasEnergy     ()                          const;
        const bool&   HasEnergyErr  ()                          const;

        const bool&   HasTime       ()                          const;
        const bool&   HasTimeErr    ()                          const;


        void          SetAmplitude  (const double& amplitude);

        void          SetEnergy     (const double& energy);
        void          SetEnergyErr  (const double& energyErr);

        void          SetTime       (const double& time);
        void          SetTimeErr    (const double& timeErr);

    // ==========================================
    //  Attributes, data
    // ==========================================

    private:

        size_t fCellIdx;

        double fAmplitude;
        bool   fHasAmplitude;

        double fEnergy;
        double fEnergyErr;
        bool   fHasEnergy;
        bool   fHasEnergyErr;

        double fTime;
        double fTimeErr;
        bool   fHasTime;
        bool   fHasTimeErr;
};

} // namespace Reco

////////////////////////////////////////////////////////////////////////////////

#endif //RecoCellDataRaw___include
