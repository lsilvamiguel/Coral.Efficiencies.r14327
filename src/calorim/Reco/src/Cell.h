/*
   $Source: /afs/cern.ch/user/g/gobbo/w0/coral/src/calorim/Reco/src/Cell.h,v $
   $Date: 2010/08/27 16:56:58 $
   $Revision: 1.26 $
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

#ifndef RecoCell___include
#define RecoCell___include

#include <iostream>
#include <string>
#include <vector>

#include "CellType.h"
#include "Correction1D.h"

namespace Reco {

////////////////////////////////////////////////////////////////////////////////

/*! \brief Cell has box shape (the same shape as GEANT's "BOX ").
*/
class Cell
{
  //============================================================================
  // Constructors and destructor
  //============================================================================

  public:

    /// Destructor
    virtual            ~Cell                    (void) {}

    /*! \brief Base constructor

        \param ct cell type
        \param activity is false if this is not real cell (but just a material)
        \param x,y,z cell's position [mm]
	\param xstep x distance to next cell (center-to-center) [mm]
	\param ystep y distance to next cell (center-to-center) [mm]
    */
                        Cell                    (const CellType &ct, bool activity,
                                                 double x, double y, double z);

  //============================================================================
  // Operators
  //============================================================================

  public:

    /// Print Cell info to output stream
    friend std::ostream     &operator <<        (std::ostream &o, const Cell &cell);

  //============================================================================
  // Methods
  //============================================================================

  public:

    /// \returns List of neighbour cells.
    const std::vector<size_t> &GetNeighbors        (void) const {return fNeighbors;}

    /// Returns list of neighbour cells.
          std::vector<size_t> &GetNeighbors        (void)       {return fNeighbors;}

    /// \return Cell's type.
    const CellType     &GetCellType             (void) const {return *fCellType;}

    /// Set cell's XYZ-position.
    void                SetPosition             (double x, double y, double z) {fPosX=x;fPosY=y;fPosZ=z;}

    /// Set cell's X-position.
    void                SetX                    (double x) {fPosX=x;}

    /// Set cell's Y-position.
    void                SetY                    (double y) {fPosY=y;}

    /// Set cell's Z-position.
    void                SetZ                    (double z) {fPosZ=z;}

    /// X-position of cell center in mm.
    double              GetX                    (void) const {return fPosX;}

    /// Y-position of cell center in mm.
    double              GetY                    (void) const {return fPosY;}

    /// Z-position of cell center in mm.
    double              GetZ                    (void) const {return fPosZ;}

    /// X-position of left side of cell in mm (including padding)
    double              GetLeft                 (void) const {return fPosX - GetCellType().GetSizeX()/2;}

    /// X-position of right side of cell in mm (including padding)
    double              GetRight                (void) const {return fPosX + GetCellType().GetSizeX()/2;}

    /// Y-position of bottom of cell in mm (including padding)
    double              GetBottom               (void) const {return fPosY - GetCellType().GetSizeY()/2;}

    /// Y-position of top of cell in mm (including padding)
    double              GetTop                  (void) const {return fPosY + GetCellType().GetSizeY()/2;}

    /// Z-position of front (upstream side) of cell in mm
    double              GetFront                (void) const {return fPosZ - GetCellType().GetSizeZ()/2;}

    /// X-position of left side of active part of cell in mm (excluding padding)
    double              GetLeftActive           (void) const {return fPosX - GetCellType().GetTrueSizeX()/2;}

    /// X-position of right side of active part of cell in mm (excluding padding)
    double              GetRightActive          (void) const {return fPosX + GetCellType().GetTrueSizeX()/2;}

    /// Y-position of bottom of active part of cell in mm (excluding padding)
    double              GetBottomActive         (void) const {return fPosY - GetCellType().GetTrueSizeY()/2;}

    /// Y-position of top of active part of cell in mm (excluding padding)
    double              GetTopActive            (void) const {return fPosY + GetCellType().GetTrueSizeY()/2;}


    /// \return Is this cell left or right, up or down neighbor and how far.
    std::pair<int,int>  NeighborMask            (const Cell &cell,double eps=0.1) const;

    /// \return true if the cells are neighborhood.
    bool                IsNeighbor              (const Cell &cell,double eps=0.1) const;

    /// \return true if the point (x,y,z) is inside the cell (including padding)
    bool                InActiveArea            (const std::vector <double> &x) const;

    /// \return true if the point (x,y,z) is inside the cell (including padding)
    bool                InActiveArea            (double x, double y, double z) const;

    /// \return true if the point (x,y,z) is inside the active part of the cell (excluding padding)
    bool                InTrueActiveArea        (double x, double y, double z) const;

    /// \return distance from the point to the cell, negative if the point is inside the cell
    double              Distance                (const std::vector <double> &x) const;

    /// \return true if the cell is on Calorimeter boundary
    bool                IsOnBoundary            (void) const {return (fboundary_regions.size() > 0);}

    /// \return true if the point (x,y,z) is on cell's boundary region
    bool                InBoundaryRegion        (const std::vector <double> &x) const;

    /// \return true if the point (x,y,z) is on cell's boundary region
    bool                InBoundaryRegion        (double x, double y, double z) const;

    /// Set cell's boundary region (asuming point(x,y,z) outside the Calorimeter
    void                SetBoundaryRegion       (const std::vector <double> &x);

    /// Return true if cell is active.
    bool                IsActive                (void) const {return fIsActive;}

    /// \return energy-dependent correction (const)
    const Correction1D &GetEdepCorr             (void) const { return fEdepCorr; }

    /// \return energy-dependent correction (non-const)
    Correction1D       &GetEdepCorrNonConst     (void) { return fEdepCorr; }

    /// \return time-in-spill-dependent correction (const)
    const Correction1D &GetTiSdepCorr             (void) const { return fTiSdepCorr; }

    /// \return time-in-spill-dependent correction (non-const)
    Correction1D       &GetTiSdepCorrNonConst     (void) { return fTiSdepCorr; }

  private:
    /// Set cell's boundary region
    void                SetBoundaryRegion       (int quarter);

  //============================================================================
  // Attributes, data
  //============================================================================

  private:

    /// If is equal to 'false' than values of fAmplitude does not meen anything.
    bool                fIsActive;

    /// Cell type - cell size, cell material, ...
    const CellType     *fCellType;

    /// X position of the cell in mm (box center)
    double              fPosX;

    /// Y position of the cell in mm (box center)
    double              fPosY;

    /// Z position of the cell in mm (box center)
    double              fPosZ;

    /// Neighboures cells of the cell.
    std::vector<size_t> fNeighbors;

    /*! \brief The list of cell's parts(quarters) which are situated near Calorimeter boundary
                2 3
                0 1    implementation of quarters numbering schema

    */
    std::vector<size_t> fboundary_regions;

    /// energy-dependent calibration corrections
    Correction1D        fEdepCorr;

    /// time-in-spill-dependent correction
    Correction1D        fTiSdepCorr;
};

////////////////////////////////////////////////////////////////////////////////

} /// using namespace Reco

#endif //RecoCell___include
