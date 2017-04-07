/*
   $Source: /afs/cern.ch/user/g/gobbo/w0/coral/src/calorim/Reco/src/Cluster.cc,v $
   $Date: 2010/04/08 16:01:02 $
   $Revision: 1.10 $
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

#include <cmath>
#include <iostream>

#include "Cluster.h"

#include "Calorimeter.h"
#include "CellDataRaw.h"

using namespace std;

namespace Reco {

template <class T> T sqr(const T &x) {return x*x;}

////////////////////////////////////////////////////////////////////////////////

Cluster::Cluster(const Calorimeter* c) :
  calorimeter(c)
{}

////////////////////////////////////////////////////////////////////////////////

double Cluster::AmplitudeTotal(void) const
{
  double result = 0;
  for(vector<CellDataRaw>::const_iterator it=cluster_data.begin(); it!=cluster_data.end(); it++ )
    result += it->GetEnergy();
  return result;
}

////////////////////////////////////////////////////////////////////////////////

double Cluster::MeanX(void) const
{
  double result = 0;
  for(vector<CellDataRaw>::const_iterator it=cluster_data.begin(); it!=cluster_data.end(); it++ )
    result += it->GetEnergy() * calorimeter->GetCells()[it->GetCellIdx()].GetX();
  return result/AmplitudeTotal();
}

////////////////////////////////////////////////////////////////////////////////

double Cluster::MeanY(void) const
{
  double result = 0;
  for(vector<CellDataRaw>::const_iterator it=cluster_data.begin(); it!=cluster_data.end(); it++ )
    result += it->GetEnergy() * calorimeter->GetCells()[it->GetCellIdx()].GetY();
  return result/AmplitudeTotal();
}

////////////////////////////////////////////////////////////////////////////////

double Cluster::MeanZ(void) const
{
  double result = 0;
  for(vector<CellDataRaw>::const_iterator it=cluster_data.begin(); it!=cluster_data.end(); it++ )
    result += it->GetEnergy() * calorimeter->GetCells()[it->GetCellIdx()].GetZ();
  return result/AmplitudeTotal();
}

////////////////////////////////////////////////////////////////////////////////

double Cluster::VarX(void) const
{
  double result = 0;
  for(vector<CellDataRaw>::const_iterator it=cluster_data.begin(); it!=cluster_data.end(); it++ )
    result += it->GetEnergy() * sqr(calorimeter->GetCells()[it->GetCellIdx()].GetX());
  return result/AmplitudeTotal()-sqr(MeanX());
}

////////////////////////////////////////////////////////////////////////////////

double Cluster::VarY(void) const
{
  double result = 0;
  for(vector<CellDataRaw>::const_iterator it=cluster_data.begin(); it!=cluster_data.end(); it++ )
    result += it->GetEnergy() * sqr(calorimeter->GetCells()[it->GetCellIdx()].GetY());
  return result/AmplitudeTotal()-sqr(MeanX());
}

////////////////////////////////////////////////////////////////////////////////

double Cluster::VarZ(void) const
{
  double result = 0;
  for(vector<CellDataRaw>::const_iterator it=cluster_data.begin(); it!=cluster_data.end(); it++ )
    result += it->GetEnergy() * sqr(calorimeter->GetCells()[it->GetCellIdx()].GetZ());
  return result/AmplitudeTotal()-sqr(MeanZ());
}

////////////////////////////////////////////////////////////////////////////////

size_t Cluster::GetMaxCell(void) const
{
  double amp_max = 0;
  if( cluster_data.size() == 0 )
    throw Exception("GetMaxCell from Cluster of zero size!");

  size_t cell_max=size_t(-1);
  for(vector<CellDataRaw>::const_iterator it=cluster_data.begin(); it!=cluster_data.end(); it++ )
  {
    if (it->GetEnergy() >  amp_max )
    {
      amp_max = it->GetEnergy();
      cell_max=it->GetCellIdx();
    }
  }
  assert(cell_max!=size_t(-1));
  return cell_max;
}

////////////////////////////////////////////////////////////////////////////////

bool Cluster::IsXYRegular(const double epsilon) const {
    if( cluster_data.size() == 0 )
        throw Exception("IsXYRegular from cluster of zero size!");

    const Cell& firstCell = calorimeter->GetCells()[cluster_data.begin()->GetCellIdx()];
    const CellType& firstCellType = firstCell.GetCellType();

    for(vector<CellDataRaw>::const_iterator it=cluster_data.begin()+1; it!=cluster_data.end(); it++ ) {
        const Cell& cell = calorimeter->GetCells()[it->GetCellIdx()];
        const CellType& cellType = cell.GetCellType();

        if (fabs(cellType.GetStepX()-firstCellType.GetStepX())>epsilon)
            return false;
        if (fabs(cellType.GetStepY()-firstCellType.GetStepY())>epsilon)
            return false;

        const double stepX  = cellType.GetStepX();
        const int    countX = ((cell.GetX()-firstCell.GetX())/stepX)>=0. ? (int)((cell.GetX()-firstCell.GetX())/stepX+0.5) : (int)((cell.GetX()-firstCell.GetX())/stepX-0.5);
        if (fabs(cell.GetX()-countX*stepX-firstCell.GetX())>epsilon)
            return false;

        const double stepY  = cellType.GetStepY();
        const int    countY = ((cell.GetY()-firstCell.GetY())/stepY)>=0. ? (int)((cell.GetY()-firstCell.GetY())/stepY+0.5) : (int)((cell.GetY()-firstCell.GetY())/stepY-0.5);
        if (fabs(cell.GetY()-countY*stepY-firstCell.GetY())>epsilon)
            return false;
    }

    return true;
}

////////////////////////////////////////////////////////////////////////////////

ostream& operator << (ostream& o, const Cluster& c)
{
  o << "Cluster: Ncells=" << c.Size()
    << ", TotAmpl=" << c.AmplitudeTotal()
    << ", MeanXYZ=(" << c.MeanX() << "," << c.MeanY() << "," << c.MeanZ() << ")"
    << ", VarXYZ=(" << c.VarX() << "," << c.VarY() << "," << c.VarZ() << ")"
    << endl;
  return o;
}

////////////////////////////////////////////////////////////////////////////////
} // namespace Reco
