/*
   $Source: /afs/cern.ch/user/g/gobbo/w0/coral/src/calorim/Reco/src/CalorimeterParticle.cc,v $
   $Date: 2010/06/29 15:47:03 $
   $Revision: 1.14 $
   -------------------------------------------------------------------------

   This file is part of cellular calorimeter reconstruction program.

   Authors:
     Vladimir  Kolosov   ( Kolosov@mx.ihep.su )

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

#include "Calorimeter.h"
#include "CalorimeterParticle.h"
#include "Shower.h"

using namespace std;

namespace Reco {

////////////////////////////////////////////////////////////////////////////////

void CalorimeterParticle::Print(ostream &o) const
{
  o << " P=" << GetProb() <<
       " E=" << GetE() << " +/- " << GetEerr() <<
       " X=" << GetX() << " +/- " <<  GetXerr() <<
       " Y=" << GetY() << " +/- " <<  GetYerr() <<
       " Z=" << GetZ() << " +/- " <<  GetZerr() <<
       " AngleX=" << GetAngleX() << " AngleY=" << GetAngleY() <<
       " NCells " << main_hited_cells.size() <<
       " Time=" << time << " +/- " << time_err << endl;
}

////////////////////////////////////////////////////////////////////////////////

void CalorimeterParticle::SetClusterSize        (size_t ncells)
{
  if( cluster_size_in_cells==0 && cluster_data.size() == 0)
  {
    cluster_size_in_cells=ncells;
  }
  else
  {
    cerr << " WARNING !!! CalorimeterParticle::SetClusterSize may not be consistent ! It is set already!! " << endl;
    cluster_size_in_cells=ncells;
  }
}

////////////////////////////////////////////////////////////////////////////////

void   CalorimeterParticle::SetClusterData     (const vector < pair< size_t,double> > &cluster)
{
  if( cluster_data.size() != 0 )
  {
    cerr << " ERROR in CalorimeterParticle::SetClusterData to not empty cluster !!! " <<
            " Cluster cleaned by force!!! MIGHT BE FATAL !!! " << endl;
  }
  ReplaceClusterData ( cluster );
}

////////////////////////////////////////////////////////////////////////////////

void   CalorimeterParticle::ReplaceClusterData     (const vector < pair< size_t,double> > &cluster)
{
  bool debug = false;
//   bool debug = true;

  if( cluster_data.size() != 0 )
  {
//     cerr << " ERROR in CalorimeterParticle::SetClusterData to not empty cluster !!! " <<
//             " Cluster cleaned by force!!! MIGHT BE FATAL !!! " << endl;
    cluster_data.clear();
  }

  if( debug )
  {
    cout << " From " << GetCalorimeterName() <<
            " CalorimeterParticle::SetClusterData " << cluster.size() << endl;
    for( size_t it=0; it!=cluster.size(); it++)
      cout << " cell " << cluster[it].first << " energy " << cluster[it].second << endl;
  }
  vector <bool> sorted;
  for( size_t it=0; it!=cluster.size(); it++)  sorted.push_back(false);
  double esum=0;
  for( size_t it=0; it!=cluster.size(); it++)
  {
    esum += cluster[it].second;
    double emax=0;
    int itmax=-1;
    for( size_t it1=0; it1!=cluster.size(); it1++)
    {
      if( !sorted[it1] &&  cluster[it1].second > emax )
      {
         emax=cluster[it1].second;
         itmax=it1;
      }
    }

    if(itmax >= 0)
    {
      cluster_data.push_back(cluster[itmax]);
      sorted[itmax]=true;
    }
    else
      cerr << " CalorimeterParticle::SetClusterData  seems bad energy in clusters! "  << endl;

  }

  if( debug )
  {
    cout << " After Sorting Cluster size " << cluster_data.size() << endl;
    for( size_t it=0; it!=cluster_data.size(); it++)
      cout << " cell " << cluster_data[it].first << " energy " << cluster_data[it].second << endl;
  }

  if( cluster_data.size() != cluster.size() )
    cerr << " CalorimeterParticle::SetClusterData internal error or bad data to sort !! " << endl;

  double ecut=0;
  size_t ic=0;
  for(size_t nc=0; nc!=cluster_data.size(); nc++)
  {
    ecut += cluster_data[nc].second;
    ic++;
    if( ecut > 0.9*esum) break;
  }
  cluster_size_in_cells=ic;
}

////////////////////////////////////////////////////////////////////////////////

void   CalorimeterParticle::SetMiscInfo( MiscDataID data_id, double data )
{
  misc_data_[data_id] = data;
}

////////////////////////////////////////////////////////////////////////////////

pair < bool,double> CalorimeterParticle::GetMiscInfo( MiscDataID data_id ) const
{
  map<MiscDataID, double>::const_iterator e = misc_data_.find(data_id);
  if( e == misc_data_.end() ) return pair < bool,double>(false,0.);
  return pair < bool,double>(true,e->second);
}

////////////////////////////////////////////////////////////////////////////////

TLorentzVector   CalorimeterParticle::GetMomentum ( void ) const
{
  double p = e;
  if( mass_ > 0. )  p = sqrt((e-mass_)*(e+mass_));
  double d = sqrt(1.+angle[0]*angle[0]+angle[1]*angle[1]);
  double pz = p/d;
  return TLorentzVector( pz*angle[0] , pz*angle[1] , pz , e );
}

////////////////////////////////////////////////////////////////////////////////

const string& CalorimeterParticle::GetCalorimeterName (void) const {
  if ( calorimeter_name == "" && calorimeter )
    return calorimeter->GetName();
  else
    return calorimeter_name;
}

////////////////////////////////////////////////////////////////////////////////

double CalorimeterParticle::CalcZforTrack( double dxdz, double dydz, bool mip ) {

  // calculate shower midpoint
  const unsigned main = GetMainCells().at(0);  // index of main cell
  const double   zcal = calorimeter->GetPositionZ();
  if ( mip ) {
    return zcal + calorimeter->GetCells().at(main).GetZ();
  } else {
    double mid;
    if ( calorimeter->GetType() == Calorimeter::ElectroMagnetic ) {
      const double crit_energy = calorimeter->GetCells().at(main).GetCellType().GetCriticalEnergy();
      mid = ZmidShowerElectron( e, crit_energy );
      mid *= calorimeter->GetCells().at(main).GetCellType().GetRadiationLength();
    } else {
      mid = ZmidShowerHadronic( e );
      mid *= calorimeter->GetCells().at(main).GetCellType().GetNuclearLength();
    }

    // apply track angle to z calculation
    const double zcell = calorimeter->GetCells().at(main).GetFront();
    return zcal + zcell + mid / sqrt( 1. + dxdz*dxdz + dydz*dydz );
  }
}

} // namespace Reco
