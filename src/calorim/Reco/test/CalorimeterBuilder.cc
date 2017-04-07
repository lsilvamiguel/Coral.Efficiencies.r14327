/*
   $Source: /afs/cern.ch/user/g/gobbo/w0/coral/src/calorim/Reco/test/CalorimeterBuilder.cc,v $ 
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

#include <stream.h>
#include "Reconstruction.h"

namespace Reco {

////////////////////////////////////////////////////////////////////////////////
  
void Calorimeter::AddCellsMatrix(vector<Cell> &cells,const CellType &type,size_t Nx,size_t Ny)
{
  if( Nx<=0 || Ny<=0 )
  {
    cerr << __PRETTY_FUNCTION__ << endl;
    cerr.form("Must Nx>0 and Ny>0.  Nx=%d Ny=%d\n",Nx,Ny);
    exit(1);
  }

  // OK. Lets start real code here.
  
  double pos_z = 0;
  for( size_t y=0; y<Ny; y++ )
  {
    double pos_y = (y-(Ny-1)/2.)*type.GetSizeY();

    for( size_t x=0; x<Nx; x++ )
    {
      double pos_x = (x-(Nx-1)/2.)*type.GetSizeX();
      cells.push_back( Cell(type,true,pos_x,pos_y,pos_z) );
    }
  }
}

////////////////////////////////////////////////////////////////////////////////
  
void Calorimeter::AddCellsMatrix(vector<Cell> &cells,const CellType &type,size_t Nx,size_t Ny,
                                 const vector <size_t> &holesXmin,const vector <size_t> &holesXmax,
                                 const vector <size_t> &holesYmin,const vector <size_t> &holesYmax)
{
  if( Nx<=0 || Ny<=0 )
  {
    cerr << __PRETTY_FUNCTION__ << endl;
    cerr.form("Must Nx>0 and Ny>0.  Nx=%d Ny=%d\n",Nx,Ny);
    exit(1);
  }

  if( holesXmin.size()!=holesXmax.size() ||
      holesXmin.size()!=holesYmin.size() ||
      holesXmin.size()!=holesYmax.size() )
  {
    cerr << __PRETTY_FUNCTION__ << endl;
    cerr.form("Bad hole sizes: %d %d %d %d\n",
      holesXmin.size(),holesXmax.size(),holesYmin.size(),holesYmax.size());
    exit(1);
  }

  // OK. Lets start real code here.
  
  double pos_z = 0;
  for( size_t y=0; y<Ny; y++ )
  {
    double pos_y = (y-(Ny-1)/2.)*type.GetSizeY();

    for( size_t x=0; x<Nx; x++ )
    {
      double pos_x = (x-(Nx-1)/2.)*type.GetSizeX();
      bool reality  = true;
      bool activity = true;
      for( size_t nhole=0; nhole<holesXmin.size(); nhole++ )
      {
        if( y >= holesYmin[nhole] && y <= holesYmax[nhole] &
            x >= holesXmin[nhole] && x <= holesXmax[nhole] ) {
          reality  = false;
          activity = false;
        }
      }
      if( reality )
        cells.push_back( Cell(type,activity,pos_x,pos_y,pos_z) );
    }
  }
}

////////////////////////////////////////////////////////////////////////////////

} // namespace Reco
