/*
   $Source: /afs/cern.ch/user/g/gobbo/w0/coral/src/calorim/Reco/src/Cell.cc,v $
   $Date: 2010/08/27 16:56:58 $
   $Revision: 1.23 $
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

#include <iostream>
#include <cmath>

#include "Cell.h"

#include "CellType.h"
#include "Exception.h"
#include "StatInfo.h"

using namespace std;

namespace Reco {

////////////////////////////////////////////////////////////////////////////////

Cell::Cell(const CellType &ct, bool activity, double x, double y, double z) :
  fIsActive     (activity),
  fCellType     (&ct),
  fPosX         (x),
  fPosY         (y),
  fPosZ         (z)
{
}

////////////////////////////////////////////////////////////////////////////////

bool Cell::InActiveArea(double x, double y, double z) const {
  if( fabs(x - fPosX ) > GetCellType().GetSizeX()/2 )
    return false;

  if( fabs(y - fPosY ) > GetCellType().GetSizeY()/2 )
    return false;

  if( fabs(z - fPosZ ) > GetCellType().GetSizeZ()/2 )
    return false;

  return true;
}

////////////////////////////////////////////////////////////////////////////////

bool Cell::InTrueActiveArea(double x, double y, double z) const {
  if( fabs(x - fPosX ) > GetCellType().GetTrueSizeX()/2 )
    return false;

  if( fabs(y - fPosY ) > GetCellType().GetTrueSizeY()/2 )
    return false;

  if( fabs(z - fPosZ ) > GetCellType().GetTrueSizeZ()/2 )
    return false;

  return true;
}

////////////////////////////////////////////////////////////////////////////////

bool Cell::InActiveArea(const vector <double> &x) const
{
  assert( x.size() == 3 );
  return InActiveArea(x[0], x[1], x[2]);
}

////////////////////////////////////////////////////////////////////////////////
/// Not precise function used in rare cases to find the nearest cell

double Cell::Distance(const vector <double> &x) const
{
  assert( x.size() == 3 );
  bool debug = false;
  if( debug ) cout << " Cell::Distance CellType " << fCellType << endl;
  if( debug ) cout << " Cell::Distance x " << x[0] << " sizex " << GetCellType().GetSizeX()/2. << endl;
  if( debug ) cout << " Cell::Distance y " << x[1] << " sizey " << GetCellType().GetSizeY()/2. << endl;
  if( debug ) cout << " Cell::Distance z " << x[2] << " sizey " << GetCellType().GetSizeZ()/2. << endl;
  double inside = -1.;
  double dx = x[0] - fPosX;
  if( debug ) cout << " dx " << dx << endl;
  double dy = x[1] - fPosY;
  if( debug ) cout << " dy " << dy << endl;
  double dz = x[2] - fPosZ;
  if( debug ) cout << " dz " << dz << endl;
  if( fabs( dx ) > GetCellType().GetSizeX()/2 )
  {
    inside = 1.;
    dx = fabs(dx) - GetCellType().GetSizeX()/2.;
  }
  if( debug ) cout << " dx " << dx << " inside " << inside << " sizex " << GetCellType().GetSizeX()/2. << endl;

  if( fabs(dy ) > GetCellType().GetSizeY()/2 )
  {
    inside = 1.;
    dy = fabs(dy) - GetCellType().GetSizeY()/2.;
  }
  if( debug ) cout << " dy " << dy << " inside " << inside <<  " sizey " << GetCellType().GetSizeY()/2. << endl;

  if( fabs(dz ) > GetCellType().GetSizeZ()/2 )
  {
    inside = 1.;
    dz = fabs(dz) - GetCellType().GetSizeZ()/2.;
  }
  if( debug ) cout << " dz " << dz << " inside " << inside << " sizez " << GetCellType().GetSizeZ()/2. <<  endl;

  double dd = dx*dx + dy*dy + dz*dz;
  if( debug ) cout << " dd " << dd << " inside " << inside << endl;
  dd = sqrt( fabs(dd) );
  if( debug ) cout << " dd " << dd*inside <<" OK!" << endl;

  return dd*inside;
}

////////////////////////////////////////////////////////////////////////////////

bool Cell::InBoundaryRegion(double x, double y, double z) const
{
  if( !IsOnBoundary()) return false;
  if( !InActiveArea(x,y,z) ) return false;
  if( fboundary_regions.size() ==4) return true;
  int dirx=0;
  int diry=0;
  if(x > fPosX ) dirx=1;
  if(y > fPosY ) diry=1;
  int in_zone = 2*diry+dirx;
  for(size_t it=0; it<fboundary_regions.size(); it++)
    if( (int)fboundary_regions[it] == in_zone) return true;

  return false;
}

////////////////////////////////////////////////////////////////////////////////

bool Cell::InBoundaryRegion(const vector <double> &x) const
{
  assert( x.size() == 3 );
  return InBoundaryRegion(x[0], x[1], x[2]);
}

////////////////////////////////////////////////////////////////////////////////

void Cell::SetBoundaryRegion( const vector <double> &x )
{
  assert( x.size() == 3 );
  double eps = 0.1;
  if( InActiveArea(x) ) return;  // Nothing to set
  int dirx=0;
  int diry=0;
  if(fabs(x[0]-fPosX) > eps && x[0] > fPosX ) dirx= 1;
  if(fabs(x[0]-fPosX) > eps && x[0] < fPosX ) dirx=-1;
  if(fabs(x[1]-fPosY) > eps && x[1] > fPosY ) diry= 1;
  if(fabs(x[1]-fPosY) > eps && x[1] < fPosY ) diry=-1;
  if(diry == 0 && dirx == 1)
  {
    SetBoundaryRegion(1);
    SetBoundaryRegion(3);
  }
  else if(diry == 0 && dirx == -1)
  {
    SetBoundaryRegion(0);
    SetBoundaryRegion(2);
  }
  else if(dirx == 0 && diry == 1)
  {
    SetBoundaryRegion(2);
    SetBoundaryRegion(3);
  }
  else if(dirx == 0 && diry == -1)
  {
    SetBoundaryRegion(0);
    SetBoundaryRegion(1);
  }
  else if(dirx != 0 && diry != 0)
  {
    if(dirx == -1) dirx=0;
    if(diry == -1) diry=0;
    SetBoundaryRegion(dirx+2*diry);
  }
  else
  {
    cout << " Cell::SetBoundaryRegion:: check code " << endl;
  }
}

////////////////////////////////////////////////////////////////////////////////

void Cell::SetBoundaryRegion(int quarter)
{
  if( quarter < 0 ||  quarter > 3 )
  {
     cout << " Wrong quarter number setting in Cell::SetBoundaryRegion(int quarter) " << quarter << endl;
  }
  for(size_t it=0; it < fboundary_regions.size(); it++)
  {
    if( (int)fboundary_regions[it] == quarter) return; // it is set already
  }
  fboundary_regions.push_back(quarter);

}

////////////////////////////////////////////////////////////////////////////////

bool Cell::IsNeighbor(Reco::Cell const &c,double eps) const
{
  if( fabs(fPosX-c.fPosX) > GetCellType().GetSizeX()/2+c.GetCellType().GetSizeX()/2+eps )
    return false;

  if( fabs(fPosY-c.fPosY) > GetCellType().GetSizeY()/2+c.GetCellType().GetSizeY()/2+eps )
    return false;

  if( fabs(fPosZ-c.fPosZ) > GetCellType().GetSizeZ()/2+c.GetCellType().GetSizeZ()/2+eps )
    return false;

  return true;
}

////////////////////////////////////////////////////////////////////////////////

pair <int,int> Cell::NeighborMask(Reco::Cell const &c,double accuracy) const
{
  double SX = GetCellType().GetSizeX();
  double SY = GetCellType().GetSizeY();

  double left   = fPosX - SX/2.;
  double right  = fPosX + SX/2.;
  double down   = fPosY - SY/2.;
  double up     = fPosY + SY/2.;

  double cfX = c.fPosX;
  double cSX = c.GetCellType().GetSizeX();
  double cfY = c.fPosY;
  double cSY = c.GetCellType().GetSizeY();

  double cleft  = cfX - cSX/2.;
  double cright = cfX + cSX/2.;
  double cdown  = cfY - cSY/2.;
  double cup    = cfY + cSY/2.;

  int bx = 0;
  int by = 0;

  if(left+accuracy  > cright ) bx = (int) (-1-(left+accuracy-cright)/SX);

  if(right-accuracy < cleft ) bx = (int) ( 1+(cleft+accuracy-right)/SX);

  if(down+accuracy > cup )    by = (int) (-1-(down+accuracy-cup)/SY);

  if(up-accuracy < cdown )    by = (int) ( 1+(cdown+accuracy-up)/SY);

  return pair<int,int>(bx,by);
}

////////////////////////////////////////////////////////////////////////////////

} // namespace Reco

