/*
   $Source: /afs/cern.ch/user/g/gobbo/w0/coral/src/calorim/Reco/src/SubSet.cc,v $
   $Date: 2010/05/27 16:12:39 $
   $Revision: 1.13 $
   -------------------------------------------------------------------------

   This file is part of cellular calorimeter reconstruction program.

   Authors:
     Vladimir  Kolosov   ( Vladimir.Kolosov@cern.ch, Kolosov@mx.ihep.su )

   Copyright(C): 1999-2004  V.Kolosov

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
#include <cmath>
#include <fstream>
#include <iostream>
#include <sstream>

// --- Internal files ----
#include "Calorimeter.h"

#include "CellType.h"
#include "Exception.h"

using namespace std;

namespace Reco {

////////////////////////////////////////////////////////////////////////////////

bool Calorimeter::SubSet::ReadConfig(const std::string &geom_file, Calorimeter &c)
{
  if( !geom_file.empty() )
  {
    ifstream f(geom_file.c_str());
    if( !f.is_open() )
      throw Exception("Calorimeter::SubSet::ReadConfig():  can not open file \"%s\"",geom_file.c_str());
    for( size_t line_n=1; true; line_n++ )
    {
      char line[555];

      f.getline( line, sizeof(line) );
      if( f.eof() )
        break;

      if( line[0]==0 || line[0]=='#' )
        continue;

      istringstream s(line);
      string opt;
      if( !(s >> opt) )
        throw Exception("Calorimeter::SubSet::ReadConfig(): bad format: \"%s\"",line);

      if( opt == "SubSet" )
      {
        string name;
        if( !(s >> name) )
          throw Exception("Calorimeter::SubSet::ReadConfig(): bad format: \"%s\"",line);
        if( name!=GetName() )
          continue;       // This is another SubSet
      }
      else
      if( opt == "add2subset_cells" )
      {
        string name;
        if( !(s >> name) )
          throw Exception("Calorimeter::SubSet::ReadConfig(): bad format: \"%s\"",line);

        if( name!=GetName() )
          continue;       // This is another SubSet

        int ncells, cell_index, step_index;

        if( !( s >> ncells
                 >> cell_index
                 >> step_index ) )
          throw Exception("Calorimeter::SubSet::ReadConfig(): bad format: \"%s\"",line);
        for (int i = 0; i < ncells; i++)
        {
          int icell = cell_index + i*step_index;
          if( icell >=0 && icell < (int)c.NCells() )
          {
            sub_set_cells.push_back( icell );
          }
          else
            throw Exception("Calorimeter::SubSet::ReadConfig(): bad cell index: %d ",icell);
        }
      }
      else
      if( opt == "add2subset_cellsxy" )
      {
        string name;
        if( !(s >> name) )
          throw Exception("Calorimeter::SubSet::ReadConfig(): bad format: \"%s\"",line);

        if( name!=GetName() )
          continue;       // This is another SubSet

        int nxcells, cell_x, step_x, nycells, cell_y, step_y;

        if( !( s >> cell_x
                 >> cell_y
                 >> nxcells
                 >> step_x
                 >> nycells
                 >> step_y  ) )
          throw Exception("Calorimeter::SubSet::ReadConfig(): bad format: \"%s\"",line);

        for (int ix = 0; ix < nxcells; ix++)
        {
          int x = cell_x + ix*step_x;
          for (int iy = 0; iy < nxcells; iy++)
          {
            int y = cell_y + iy*step_y;
            int icell = c.GetCellOfColumnRow(x,y);
            if( icell >=0 && icell < (int)c.NCells() )
            {
              sub_set_cells.push_back( icell );
            }
            else
              throw Exception("Calorimeter::SubSet::ReadConfig(): bad cell index: x= %d y= %d ",x,y);
          }
        }

      }
      else
        cerr << "Calorimeter::SubSet::ReadConfig():   Unknown keyword " << opt << endl;
        return false;
    }
  }
  return true;
}

////////////////////////////////////////////////////////////////////////////////

bool Calorimeter::SubSet::WriteConfig(const std::string &geom_file, Calorimeter &c) const
{
  ofstream out( geom_file.c_str(), ios::out );
  if( !out.is_open() )
  {
    cerr << " Calorimeter::SubSet::WriteConfig(): Can not open file  " << geom_file << endl;
  }
  else
  {
    out << "# SubSet " << GetName() << " of Calorimeter  " << c.GetName()  << endl;
    if( c.XYRegularGrid() )
    {
       for ( int i=0; i < (int)sub_set_cells.size(); i++)
       {
         int icell = sub_set_cells[i];
         int x = c.GetColumnOfCell( icell );
         int y = c.GetRowOfCell( icell);
         out << " add2subset_cells " << GetName() << " " << x << " " << y  << " " << 1 <<" "<< 1 <<" "<< 1 <<" "<< 1 << endl;
       }
    }
    else
    {
       for ( int i=0; i < (int)sub_set_cells.size(); i++)
       {
         int icell = sub_set_cells[i];
         out << " add2subset_cellsxy " << GetName() << " " << icell << " " << 1 << endl;
       }
    }

  }
  out.close();
  return true;
}

////////////////////////////////////////////////////////////////////////////////

double Calorimeter::SubSet::GetXmin( void ) const
{
  if( Size() == 0 ) return 0.;
  size_t it0 =GetCells()[0];
  double xs = calorimeter->GetCells()[it0].GetCellType().GetSizeX();
  double xc = calorimeter->GetCells()[it0].GetX();
  double xmin = xc-xs/2.;
  for( size_t iit=1; iit < Size(); iit++ )
  {
    size_t it = GetCells()[iit];
    xs = calorimeter->GetCells()[it].GetCellType().GetSizeX();
    xc = calorimeter->GetCells()[it].GetX();
    if( xmin > xc-xs/2. ) xmin = xc-xs/2.;
  }
  return xmin;
}

////////////////////////////////////////////////////////////////////////////////

double Calorimeter::SubSet::GetXmax( void ) const
{
  if( Size() == 0 ) return 0.;
  size_t it0 =GetCells()[0];
  double xs = calorimeter->GetCells()[it0].GetCellType().GetSizeX();
  double xc = calorimeter->GetCells()[it0].GetX();
  double xmax = xc+xs/2.;
  for( size_t iit=1; iit < Size(); iit++ )
  {
    size_t it = GetCells()[iit];
    xs = calorimeter->GetCells()[it].GetCellType().GetSizeX();
    xc = calorimeter->GetCells()[it].GetX();
    if( xmax < xc+xs/2. ) xmax = xc+xs/2.;
  }
  return xmax;
}

////////////////////////////////////////////////////////////////////////////////

double Calorimeter::SubSet::GetYmin( void ) const
{
  if( Size() == 0 ) return 0.;
  size_t it0 =GetCells()[0];
  double ys = calorimeter->GetCells()[it0].GetCellType().GetSizeY();
  double yc = calorimeter->GetCells()[it0].GetY();
  double ymin = yc-ys/2.;
  for( size_t iit=1; iit < Size(); iit++ )
  {
    size_t it = GetCells()[iit];
    ys = calorimeter->GetCells()[it].GetCellType().GetSizeY();
    yc = calorimeter->GetCells()[it].GetY();
    if( ymin > yc-ys/2. ) ymin = yc-ys/2.;
  }
  return ymin;
}

////////////////////////////////////////////////////////////////////////////////

double Calorimeter::SubSet::GetYmax( void ) const
{
  if( Size() == 0 ) return 0.;
  size_t it0 =GetCells()[0];
  double ys = calorimeter->GetCells()[it0].GetCellType().GetSizeY();
  double yc = calorimeter->GetCells()[it0].GetY();
  double ymax = yc+ys/2.;
  for( size_t iit=1; iit < Size(); iit++ )
  {
    size_t it = GetCells()[iit];
    ys = calorimeter->GetCells()[it].GetCellType().GetSizeY();
    yc = calorimeter->GetCells()[it].GetY();
    if( ymax < yc+ys/2. ) ymax = yc+ys/2.;
  }
  return ymax;
}

////////////////////////////////////////////////////////////////////////////////

int Calorimeter::SubSet::GetCellOfColumnRow(int x, int y) const
{
//   bool debug = true;
//   if( debug ) cout << " xy_structure_init_ " << xy_structure_init_ << endl;
  if( !xy_structure_init_ ) return -1;
//   if( debug ) cout << " fNcols_ " << fNcols_ <<  " fNrows_ " << fNrows_ << endl;
  if( x>=0 && x<fNcols_ && y>=0 && y<fNrows_ )
  {
//     if( debug ) cout << "  map_cell_xy_.size " << sizeof(map_cell_xy_) <<  " ask for element " << x+y*fNcols_ << endl;
//     if( debug ) cout << " Nu i ?? " << endl;
    int index = map_cell_xy_[x+y*fNcols_];
//     if( debug ) cout << " index " << index << endl;
//     if( debug ) cout << " sub_set_cells.size() " << sub_set_cells.size() << endl;
    if( index >= 0 ) return sub_set_cells[index];
  }
  return -1;
}

////////////////////////////////////////////////////////////////////////////////

bool  Calorimeter::SubSet::IsMyCell(int icell) const
{
  return ( GetCellIndex( icell) >= 0);
}

////////////////////////////////////////////////////////////////////////////////

int  Calorimeter::SubSet::GetCellIndex(int icell) const
{
  for ( int i=0; i< (int)Size(); i++ )
  {
    if( (int)GetCells()[i] == icell ) return i;
  }
  return -1;
}

////////////////////////////////////////////////////////////////////////////////

int  Calorimeter::SubSet::GetColumnOfCell(int icell) const
{
  int indxcell = GetCellIndex(icell);
  if( indxcell < 0 ) return -1;
  if( !xy_structure_init_ ) return -1;
  if( indxcell>=0 && indxcell < (int) Size() ) return map_xy_cell_[0][indxcell];
//   cerr << " Calorimeter::SubSet::GetColumnOfCell not implemented " << endl;
  return -1;
}

////////////////////////////////////////////////////////////////////////////////

int Calorimeter::SubSet::GetRowOfCell(int icell) const
{
  int indxcell = GetCellIndex(icell);
  if( indxcell < 0 ) return -1;
  if( !xy_structure_init_ ) return -1;
  if( indxcell>=0 && indxcell < (int)Size() ) return map_xy_cell_[1][indxcell];
//   cerr << " Calorimeter::SubSet::GetRowOfCell not implemented " << endl;
  return -1;
}

////////////////////////////////////////////////////////////////////////////////

bool Calorimeter::SubSet::InitXY(double tolerance)
{
//  bool debug = true;
  bool debug = false;
  if( debug )  cout << " Calorimeter " << calorimeter->GetName() << " subset" << GetName() << " InitXY started" << endl;

  if( Size() <= 0 )
  {
    cerr<< " Error InitXY: Empty subset " << " Calorimeter " << calorimeter->GetName() << " subset" << GetName() << endl;
    exit(1);
  }


  // Make sure that the z position of the front surface is identical.
  // (This must be checked for ECAL1, too, therefore we have to put it before
  // the regular grid test which aborts in the ECAL1 case.)
  {
    double zmin =  1e9;
    double zmax = -1e9;
    for( size_t iit=0; iit < Size(); iit++ ) {
      size_t it = GetCells()[iit];
      const Cell& cell = calorimeter->GetCells()[it];

      // z position of front edge
      const double zfront = cell.GetZ() - cell.GetCellType().GetSizeZ()/2.;
      if ( zfront < zmin )
        zmin = zfront;
      if ( zfront > zmax )
        zmax = zfront;
    }

    if ( zmax-zmin > 0.00001 )
      throw Exception("Calorimeter %s front surface z delta %gmm!  Fix detectors.dat!",
              GetName().c_str(), zmax-zmin);

    if ( debug )
      cout << "Calorimeter " << GetName() << " front surface z delta: " << zmax-zmin << endl;
  }

  // make sure that cells lie on a regular grid
  int i0             = GetCells()[0];
  xy_structure_init_ = true;
  fXStep_            = calorimeter->GetCells()[i0].GetCellType().GetStepX();
  fYStep_            = calorimeter->GetCells()[i0].GetCellType().GetStepY();
  const double XSize = calorimeter->GetCells()[i0].GetCellType().GetTrueSizeX();
  const double YSize = calorimeter->GetCells()[i0].GetCellType().GetTrueSizeY();
  fNcols_            = (int) rint( (GetXmax()-GetXmin()) / fXStep_ );
  fNrows_            = (int) rint( (GetYmax()-GetYmin()) / fYStep_ );

  const double xmin = GetXmin();
  const double ymin = GetYmin();
  for( size_t iit=0; iit < Size(); iit++ ) {
    size_t it = GetCells()[iit];
    const Cell& cell = calorimeter->GetCells()[it];
    // center of cell in calorimeter coordinates [mm]
    const double xc     = cell.GetX() - xmin;
    const double yc     = cell.GetY() - ymin;
    // cell position in integer coordinates
    const int     x     = (int) rint( (xc-XSize/2.) / fXStep_ );
    const int     y     = (int) rint( (yc-YSize/2.) / fYStep_ );
    // nominal (expected) center of cell in calo coords [mm]
    const double xc_nom = x*fXStep_ + XSize/2.;
    const double yc_nom = y*fYStep_ + YSize/2.;
    if( debug )
      cout << "InitXY:: Cell " << it
           << " dx =" << xc - xc_nom << " dy =" << yc - yc_nom << endl;

    if( fabs( xc - xc_nom ) > tolerance || fabs( yc - yc_nom ) > tolerance ) {
      cerr << "Calorimeter " << GetName() << " SubSet::InitXY Check failed at "
           << "cell index " << it << ", X=" << x << ", Y=" << y
           << ", xc=" << xc << ", yc=" << yc
           << ", nominal xc=" << xc_nom << ", nominal yc=" << yc_nom
           << ", nominal xstep=" << fXStep_ << ", nominal ystep=" << fYStep_
	   << endl;
      goto error;
    }
  }

  if ( debug )
    cout << "Calorimeter " << calorimeter->GetName()
         << " SubSet " << GetName() << " InitXY Check OK " << endl;

  map_cell_xy_.clear();
  map_xy_cell_[0].clear();
  map_xy_cell_[1].clear();
  for( int it=0; it!=fNcols_*fNrows_; it++ ) {
    map_cell_xy_.push_back(-1);
  }
  for( size_t it=0; it<Size(); it++ ) {
    map_xy_cell_[0].push_back(-1);
    map_xy_cell_[1].push_back(-1);
  }

  for( size_t iit=0; iit < Size(); iit++ ) {
    size_t it = GetCells()[iit];
    const Cell& cell = calorimeter->GetCells()[it];
    const int x = (int) rint( (cell.GetX() - xmin - fXStep_/2.) / fXStep_ );
    const int y = (int) rint( (cell.GetY() - ymin - fYStep_/2.) / fYStep_ );
    if( x < 0 || x >= fNcols_ || y < 0 || y >= fNrows_ ) {
      cerr << " ERROR xy setting: Cell outside SubSet " << GetName()
           << " in calorimeter " << calorimeter->GetName()
           << " description: x= " << x << ",  y= " << y << endl;
      goto error;
    }
    if ( map_cell_xy_[x + fNcols_*y] != -1 ) {
      cerr << " ERROR xy setting: Duplicate cell in SubSet " << GetName()
           << " in calorimeter " << calorimeter->GetName()
           << " description: x= " << x << ",  y= " << y << endl;
      goto error;
    }
    map_cell_xy_[x + fNcols_*y] = iit;
    map_xy_cell_[0][iit] = x;
    map_xy_cell_[1][iit] = y;
  }

  return true;

 error:

  xy_structure_init_ = false;
  fNcols_ = -1;
  fNrows_ = -1;
  fXStep_ = 0.;
  fYStep_ = 0.;
  return false;
}

////////////////////////////////////////////////////////////////////////////////

std::ostream& operator<<(std::ostream& o, Calorimeter::SubSet& s) {
    o << "SubSet " << s.GetName() << " of " << s.calorimeter->GetName() << ":" << endl;
    if ( s.xy_structure_init_ ) {
        o << "  Ncols=" << s.fNcols_ << ", Nrows=" << s.fNrows_
          << ", Xstep=" << s.fXStep_ << ", Ystep=" << s.fYStep_ << endl;
        o << "  size=" << s.Size() << endl;

        for (size_t i=0; i<s.Size(); i++) {
          o << "  " << i << ": x=" << s.map_xy_cell_[0][i] << ", y=" << s.map_xy_cell_[1][i] << endl;
        }
    } else {
        o << "  not initialized." << endl;
    }

    return o;
}

////////////////////////////////////////////////////////////////////////////////

} /// using namespace Reco
