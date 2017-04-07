/*
   $Source: /afs/cern.ch/user/g/gobbo/w0/coral/src/calorim/Reco/test/CalorimeterGAMS.h,v $ 
   $Date: 2000/07/25 15:33:10 $ 
   $Revision: 1.1 $ 
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

#ifndef CalorimeterGAMS___include
#define CalorimeterGAMS___include

#include "Calorimeter.h"

namespace Reco {

class CalorimeterGAMS : public Calorimeter
{
  private:
  
    class Hole
    {
      public:
        Hole(void) : Xmin(0), Xmax(0), Ymin(0), Ymax(0) {}
        Hole(size_t x_min,size_t x_max,size_t y_min,size_t y_max) :
            Xmin(x_min), Xmax(x_max), Ymin(y_min), Ymax(y_max) {}
        size_t Xmin, Xmax, Ymin, Ymax;   // [Xmin,Xmax][Ymin,Ymax]
        size_t SizeX(void) const {return (Xmax-Xmin+1);}
        size_t SizeY(void) const {return (Ymax-Ymin+1);}
        size_t Size (void) const {return SizeX()*SizeY();}
    };

  public:
  
    /// Destructor.
    virtual            ~CalorimeterGAMS         (void) {}

    /// Copy constructor.
                        CalorimeterGAMS         (const CalorimeterGAMS &gams) : Calorimeter(gams.GetName()) {*this=gams;}

    /// Base constructor.
                        CalorimeterGAMS         (const string &the_name,int type=0);

    /// \return Cell number (in list GetCells()) with point (x,y,z) inside or return -1 if there is no such cell.
    virtual int         FindCell                (double x,double y,double z) const;
    
    /// \return cell number on cell's coordinates.
    int                 GetGAMSCell             (size_t x,size_t y) const;

    /// \return cell number on cell's coordinates.
    int                 GetSADCell              (size_t x,size_t y) const;

    /// Assignment operator.
    CalorimeterGAMS    &operator =              (const CalorimeterGAMS &gams) {throw "CalorimeterGAMS::operator =    not implemented.";}

    size_t              GetSADFirstCell         (void) const {return SAD_cells;}

    /// GAMS size.
    size_t              GAMS_Nx, GAMS_Ny;    

    /// This is GAMS hole.
    Hole                GAMS_hole;

    /// SAD size. If Nx=0 SAD was not installed.
    size_t              SAD_Nx, SAD_Ny;
    
    /// This is SAD hole for beam.
    Hole                SAD_hole;
    
    /// First SAD cell number.
    size_t              SAD_cells;
};

} // namespace Reco

#endif // CalorimeterGAMS___include
