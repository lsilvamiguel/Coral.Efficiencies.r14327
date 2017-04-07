/*
   $Source: /afs/cern.ch/user/g/gobbo/w0/coral/src/calorim/Reco/src/OneParticleResponse.cc,v $
   $Date: 2011/02/01 22:05:52 $
   $Revision: 1.45 $
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

#include <cassert>
#include <cmath>
#include <iostream>
#include <algorithm>   // for find()

#include "OneParticleResponse.h"

#include "mycernlib.h"
#include "Calorimeter.h"
#include "CellDataRaw.h"
#include "CellType.h"
#include "Shower.h"

using namespace std;

namespace Reco {

////////////////////////////////////////////////////////////////////////////////

OneParticleResponse::OneParticleResponse(const CalorimeterParticle& p, const Calorimeter* c,
                                         const Cell& cell, const vector<pair<double, double> >& realAmplitudes) :
  calorimeter   (c),
  fParticle     (p),
  e_overlap_    (0.)
{
  fEnergy_deposit = 0.;

  if( fParticle.GetE()<=0 )
    throw Exception("OneParticleResponse::OneParticleResponse1: negative energy: %g GeV",fParticle.GetE());

  // Get the list of cells, which are part of this particle (cluster)
  // and initilize vectors to hold information
  FindList(cell);

  // # warning OneParticleResponse:: Only for gamma in ECAL ExtendList() valid and depends on the cell size!!!

  // Extend list of cells by neigbors of current list
  if( fParticle.GetE() > calorimeter->GetOptions().ecut_extend_list )
    ExtendList();

  // Throw exception if there are no cells in cluster
  if( hited_cells.size()==0 )
    throw Exception("OneParticleResponse::OneParticleResponse1:  size of hited_cells vector =0 ");

  fParticle.SetHitedCell(hited_cells[0]);

  // Propagate particles z position dependent from particle id
  // and dimensions of the main cell to a z position within
  // the calorimeter cell
  PropagateParticleAtCorrectZ();

  // Calculate the expected cell response dependent from particle id position and energy
  CalculateExpectedResponse();

  // Loop over cluster cells to get real response
  for( int i =0; i !=(int)hited_cells.size(); i++ )
    {
      // Get global index
      int jcell = hited_cells[i];

      // Get real response
      e_real_[i]  = realAmplitudes[jcell].first;
      de_real_[i] = realAmplitudes[jcell].second;
    }
}

////////////////////////////////////////////////////////////////////////////////

OneParticleResponse::OneParticleResponse(const CalorimeterParticle &p,const Calorimeter* c,bool for_development ) :
  calorimeter (c),
  fParticle   (p),
  e_overlap_    (0.)
{

  bool debug=false;
  if(debug) cout << " OneParticleResponse::OneParticleResponse try to.. E= " << fParticle.GetE() <<
              " X=" << fParticle.GetX() << " Y=" << fParticle.GetY() << endl;

  fEnergy_deposit = 0.;

  cells_trace.clear();
  v_input_trace_x.clear();
  v_input_trace_y.clear();
  v_input_trace_z.clear();
  v_output_trace_x.clear();
  v_output_trace_y.clear();
  v_output_trace_z.clear();

  if( fParticle.GetE()<=0 )
    throw Exception("OneParticleResponse::OneParticleResponse2: negative energy: %g GeV",fParticle.GetE());

  // Get the list of cells, which are part of this particle (cluster)
  // and initilize vectors to hold information
  FindList();

  // Return if there are no cells in cluster
  if( hited_cells.empty() )
    return;

  // # warning OneParticleResponse:: Only for gamma in ECAL ExtendList() valid and depends on the cell size!!!

  // Extend list of cells by neigbors of current list
  if( fParticle.GetE() > calorimeter->GetOptions().ecut_extend_list )
    ExtendList();

  // Throw exception if there are no cells in cluster (cannot be if !hited_cells.empty() above)
  if( hited_cells.size()==0 )
    throw Exception("OneParticleResponse::OneParticleResponse2:  size of hited_cells vector =0 ");

  fParticle.SetHitedCell(hited_cells[0]);

  // Propagate particles z position dependent from particle id
  // and dimensions of the main cell to a z position within
  // the calorimeter cell
  PropagateParticleAtCorrectZ();

  if( for_development )                  // Development =  FMC
    CalculateExpectedResponseForFMC();
  else
    // Calculate the expected cell response dependent from particle id position and energy
    CalculateExpectedResponse();
}

////////////////////////////////////////////////////////////////////////////////

// # warning OneParticleResponse:: Could be done directly from cluster
OneParticleResponse::OneParticleResponse(const CalorimeterParticle &p, const Calorimeter* c,
                                                                    vector<CellDataRaw>& cluster_data) :
  calorimeter   (c),
  fParticle     (p),
  e_overlap_    (0.)
{
  fEnergy_deposit = 0.;
  if( fParticle.GetE()<=0 )
    throw Exception("OneParticleResponse::OneParticleResponse3: negative energy: %g GeV",fParticle.GetE());
  hited_cells.clear();
//   expected_one_particle_response_in_cells.clear();
//   real_one_particle_response_in_cells.clear();
  e_expected_.clear();
  de_expected_.clear();
  e_real_.clear();
  de_real_.clear();
  shower_fluctuations_one_particle_response.clear();
  fMaskCells.clear();
  de_one_particle_response_in_cells.clear();
  dx_one_particle_response_in_cells.clear();
  dy_one_particle_response_in_cells.clear();

// # warning OneParticleResponse:: Seems too heavy initialisation. Need to be improved.
  int p_main_cell;
  if(fParticle.GetMainCells().size() <=0 )
  {
//      p_main_cell = calorimeter->FindCell(fParticle.GetX(),fParticle.GetY(),fParticle.GetZ());
      p_main_cell = calorimeter->FindCell(fParticle);
  }
  else
      p_main_cell = fParticle.GetMainCells()[0];

// # warning OneParticleResponse:: !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! cluster_data.size() == 0
  static bool first_time=true;
  if( cluster_data.size() == 0 )
  {
    if( first_time )
    {
      cerr <<  calorimeter->GetName() <<
        " OneParticleResponse::OneParticleResponse3:  size of cluster_data vector =0 " << endl;
      cerr <<  calorimeter->GetName() <<
        " Sovershenneishii BARDAK !!!! " << endl;
      first_time=false;
    }
    return;
  }

  for(vector<CellDataRaw>::const_iterator it=cluster_data.begin(); it!=cluster_data.end(); it++ )
  {
    size_t icell = it->GetCellIdx();
    double e =it->GetEnergy();
    hited_cells.push_back(icell);

    fMaskCells.push_back(calorimeter->GetCells()[p_main_cell].NeighborMask(calorimeter->GetCells()[icell],
                             calorimeter->GetOptions().tolerance_for_nearby_cells));
//     expected_one_particle_response_in_cells   .push_back( pair<double,double>(0,0) );
//     real_one_particle_response_in_cells       .push_back( pair<double,double>(e,0) );
    e_expected_    .push_back( 0. );
    de_expected_   .push_back( 0. );
    e_real_        .push_back( e );
    de_real_       .push_back( 0. );
    shower_fluctuations_one_particle_response .push_back( pair<double,double>(0,0) );
    de_one_particle_response_in_cells.push_back(0);
    dx_one_particle_response_in_cells.push_back(0);
    dy_one_particle_response_in_cells.push_back(0);

  }

  if( hited_cells.size()==0 )
  {
    cerr <<  calorimeter->GetName() <<
     " OneParticleResponse::OneParticleResponse3:  size of hited_cells vector =0 " << endl;
    throw Exception("OneParticleResponse::OneParticleResponse3:  size of hited_cells vector =0 ");
  }
  PropagateParticleAtCorrectZ();
  CalculateExpectedResponse();
}

////////////////////////////////////////////////////////////////////////////////

OneParticleResponse::OneParticleResponse (const OneParticleResponse &h) :
    calorimeter(h.calorimeter),
    fParticle(h.fParticle),
    fEnergy_deposit(h.fEnergy_deposit),
    fXiSqr(h.fXiSqr),
    fNDF(h.fNDF),
    e_overlap_(h.e_overlap_),
    hited_cells(h.hited_cells),
//     expected_one_particle_response_in_cells(h.expected_one_particle_response_in_cells),
//     real_one_particle_response_in_cells(h.real_one_particle_response_in_cells),
    e_real_(h.e_real_),
    de_real_(h.de_real_),
    e_expected_(h.e_expected_),
    de_expected_(h.de_expected_),
    shower_fluctuations_one_particle_response(h.shower_fluctuations_one_particle_response),
    fMaskCells(h.fMaskCells),
    de_one_particle_response_in_cells(h.de_one_particle_response_in_cells),
    dx_one_particle_response_in_cells(h.dx_one_particle_response_in_cells),
    dy_one_particle_response_in_cells(h.dy_one_particle_response_in_cells),
    cells_trace(h.cells_trace),
    v_input_trace_x(h.v_input_trace_x),
    v_input_trace_y(h.v_input_trace_y),
    v_input_trace_z(h.v_input_trace_z),
    v_output_trace_x(h.v_output_trace_x),
    v_output_trace_y(h.v_output_trace_y),
    v_output_trace_z(h.v_output_trace_z)
{}

////////////////////////////////////////////////////////////////////////////////

OneParticleResponse &OneParticleResponse::operator =
                                     (const OneParticleResponse &h)
{
  if( &h!=this )
  {
//    cout << " Copy OneParticleResponse " << h.calorimeter->GetName() << endl;
    calorimeter = h.calorimeter;
    fParticle = h.fParticle;
    fEnergy_deposit = h.fEnergy_deposit;
    fXiSqr = h.fXiSqr;
    fNDF = h.fNDF;
    e_overlap_=h.e_overlap_;
    hited_cells = h.hited_cells;
//     expected_one_particle_response_in_cells = h.expected_one_particle_response_in_cells;
//     real_one_particle_response_in_cells = h.real_one_particle_response_in_cells;
    e_expected_ = h.e_expected_;
    de_expected_ = h.de_expected_;
    e_real_ = h.e_real_;
    de_real_ = h.de_real_;
    shower_fluctuations_one_particle_response = h.shower_fluctuations_one_particle_response;
    fMaskCells = h.fMaskCells;
    de_one_particle_response_in_cells = h.de_one_particle_response_in_cells;
    dx_one_particle_response_in_cells = h.dx_one_particle_response_in_cells;
    dy_one_particle_response_in_cells = h.dy_one_particle_response_in_cells;
    cells_trace = h.cells_trace;
    v_input_trace_x = h.v_input_trace_x;
    v_input_trace_y = h.v_input_trace_y;
    v_input_trace_z = h.v_input_trace_z;
    v_output_trace_x = h.v_output_trace_x;
    v_output_trace_y = h.v_output_trace_y;
    v_output_trace_z = h.v_output_trace_z;
  }
  return *this;
}

////////////////////////////////////////////////////////////////////////////////

void OneParticleResponse::PropagateParticleAtCorrectZ (void)
{
  // return if there are no cells in cluster
  if( hited_cells.size() == 0 )
    return;

  // Take first cell in cluster cells as main cell
  const Cell* p_main_cell = &calorimeter->GetCells()[hited_cells.front()];

  // Get the half Z Dimension of the main cell
  double SizeZ = p_main_cell->GetCellType().GetSizeZ()/2.;

  // Get the particle id
  CalorimeterParticle::ParticleID pid = fParticle.GetID();

  // Get the energy of the particle
  double E = fParticle.GetE();

  // Get z position dependent from particle id
  double zw;
  if ( 1 <= pid && pid <= 3 )   // gamma,electron,positron
    {
      double RadLeng = p_main_cell->GetCellType().GetRadiationLength();
      zw = ZmidShowerEM(E)*RadLeng;
    }
  else if ( pid == 8 || pid == 9 ) // pi+, pi-
    {
      double NuclLeng = p_main_cell->GetCellType().GetNuclearLength();
      zw = ZmidShowerHadronic(E)*NuclLeng;
    }
  else if ( pid == 0)
    {
      cerr << "Unknown particle in calorimeter!" << endl;
      zw = SizeZ;
    }
  else
    {
      cerr << "How do you know particle ID is " << pid << ", eh?  "
	   << "Using MC truth here?  Fix in source!" << endl;
      exit(1);
    }

  double znew = p_main_cell->GetZ() - SizeZ + zw;

  // Set new z position.  (x,y) position is determined by reconstruction and
  // must stay as it is!
  fParticle.SetZ(znew);
}

////////////////////////////////////////////////////////////////////////////////

#if 0 // suhl
pair< double, size_t> OneParticleResponse::CheckSinglePeakHypothesis (void) const
{
  if( hited_cells.empty() )
    return pair< double, size_t>(0.,0);

  const vector< size_t > &near = calorimeter->GetCells()[hited_cells[0]].GetNeighbors();
  double sum = 0.;
  size_t i =0;
  for( vector<size_t>::const_iterator it=hited_cells.begin(); it!=hited_cells.end(); it++ )
  {
    if( i == 0 ) continue;
    if( i >= near.size() ) break;

    double e = calorimeter->real_amplitudes[*it].first;
//     double de = calorimeter->real_amplitudes[*it].second -
//                                    real_one_particle_response_in_cells[i].second;
    double de = calorimeter->real_amplitudes[*it].second - de_real_[i];
    sum += (e*e)/de;
    i++;
  }
  return pair< double, size_t>(sum, near.size()-1);
}
#endif // suhl

////////////////////////////////////////////////////////////////////////////////

void OneParticleResponse::MonteRandomise (void)
{
  if( !calorimeter->GetOptions().add_shower_fluctuations_for_FMC )
    return;
  if( hited_cells.empty() )
    return;

  static double esum_before=0,esum_after=0;
  esum_before=0;
  esum_after=0;
//   cout << "Init esum_before = " << esum_before << "Init esum_after " << esum_after  << "PI " << M_PI <<endl;
  size_t i =0;
  for( vector<size_t>::iterator it=hited_cells.begin(); it!=hited_cells.end(); it++ ) {

      const double r1 = calorimeter->GetRandomGaus();

//       cout << " e " << i << "=" <<expected_one_particle_response_in_cells[i].first
//            << " SIGMA**2 " << expected_one_particle_response_in_cells[i].second << endl;

//       esum_before += expected_one_particle_response_in_cells[i].first;

      esum_before += e_expected_[i];
//       expected_one_particle_response_in_cells[i].first +=
//                       r1*sqrt(expected_one_particle_response_in_cells[i].second);
      e_expected_[i] +=
                      r1*sqrt(de_expected_[i]);
//       esum_after += expected_one_particle_response_in_cells[i].first;
      esum_after += e_expected_[i];
       i++;
  }

  i =0;

//   cout << "esum_before = " << esum_before << " esum_after " << esum_after  << endl;
//   if( !(esum_before > 0 && esum_after > 0) )
//      cout << "esum_before = " << esum_before << " esum_after " << esum_after  << endl;

  for( vector<size_t>::iterator it=hited_cells.begin(); it!=hited_cells.end(); it++ ) {
//       expected_one_particle_response_in_cells[i].first =
//               expected_one_particle_response_in_cells[i].first*esum_before/esum_after;
      e_expected_[i] = e_expected_[i]*esum_before/esum_after;
       i++;
  }
}

////////////////////////////////////////////////////////////////////////////////

/*!
    \brief Calculate expected amplitudes for cells in a cluster
    \callgraph
    \callergraph
    \remarks {
    Used options:
    rd_fiadc_digitization
    }
*/


void OneParticleResponse::CalculateExpectedResponse()
{
    std::vector<std::pair<double, double> > tmp;
    tmp.resize(calorimeter->NCells(), std::pair<double, double>(0., 0.));

    // Loop over hit cells to get real response
    for (size_t i=0; i<hited_cells.size(); i++) {
        // Get global index
        size_t jcell = hited_cells[i];

        // Get real response
        tmp[jcell].first  =  e_real_[i];
        tmp[jcell].second = de_real_[i];
    }

    CalculateExpectedResponse(tmp);
}

void OneParticleResponse::CalculateExpectedResponse(const std::vector<std::pair<double, double> >& realAmplitudes)
{

  // return if no cells in cluster
  if( hited_cells.empty() )
    return;

  // Get first cell as main cell
  const Cell*  p_main_cell = &(calorimeter->GetCells()[hited_cells.front()]);


  double etotal_expected(0.);

  // loop over cells
  for( int i =0; i !=(int)hited_cells.size(); i++ )
  {

    // Get cell index
    int jcell = hited_cells[i];

    // Get the cell dependent on index (assert that you are in valid range of indeces)
    assert( jcell>=0 && jcell < (int)calorimeter->GetCells().size() );
    const Cell*  pcell = &(calorimeter->GetCells()[jcell]);

    // Output of next function store results in expected_one_particle_response_in_cells[i] array
    CalculateParticleResponseInCell( i,pcell,p_main_cell,realAmplitudes, false);

    //  etotal_expected += expected_one_particle_response_in_cells[i].first;
    etotal_expected += e_expected_[i];
  }

  // We apply sparce energy cut to the expected ampitude which is not 100% correct We need to reconsider this solution later !??
  //  if( calorimeter->GetOptions().rd_fiadc_digitization )  it seems we have options doubling

  if( calorimeter->GetOptions().readout_sparsified_ )
  {
    etotal_expected = 0.;
//    double delta = calorimeter->GetOptions().rd_fiadc_sparce_delta;

    for( int i =0; i !=(int)hited_cells.size(); i++ )
    {
      int jcell = hited_cells[i];
//  it seems we have options doubling
//      double e_sparce_cut = delta*calorimeter->GetCellInfo( Calorimeter::CALIB, Calorimeter::OLD, jcell).GetMean();
      double e_sparce_cut = calorimeter->GetEnergyCutSparseMode()[jcell];
      double et = e_expected_[i];
      if( et < e_sparce_cut )
      {
        et=0.01;
      }
      etotal_expected += et;
    }
  }
  fEnergy_deposit=etotal_expected;
}

////////////////////////////////////////////////////////////////////////////////

void OneParticleResponse::CalculateExpectedResponseForFMC( void )
{
  bool debug=false;
//  bool debug=true;
  for( int i =0; i !=(int)hited_cells.size(); i++ )
  {
//     expected_one_particle_response_in_cells[i] = pair<double,double>(0.,0.);
    e_expected_[i] = 0.;
    de_expected_[i] = 0.;
  }
  if(debug) cout << " CalculateExpectedResponse in calorimeter " << calorimeter->GetName() << endl;
//  cout << " OneParticleResponse::CalculateExpectedResponse in Ncells= " << hited_cells.size() << endl;
  if( hited_cells.empty() )
    return;

// At this point we need to decide what kind of interactions occurs inside the detector to simulate response.
// Consider possible solutions:
//                  --   no interactions ( nothing for neutrals, MIP for charge particles)
//                  --   Electro-Magnetic shower starting at some point along Z direction
//                  --   Hadronic shower starting at some point along Z direction
//   For the moment this method is used only as FastMonteCarlo detector response simulation.
//   And must be somehow protected (for a while) from usage in reconstruction.

//  bool with_derivatives=false;  // As it is FMC method we dont need any derivatives calculation

  vector<double> l_trace;
  int npoints = cells_trace.size();
  if(debug) cout << " Points on trace = " << npoints << endl;
  for( int it=0; it != npoints; it++ )
  {
    double dx = v_input_trace_x[it]-v_input_trace_x[0];
    double dy = v_input_trace_y[it]-v_input_trace_y[0];
    double dz = v_input_trace_z[it]-v_input_trace_z[0];
    l_trace.push_back( sqrt(dx*dx+dy*dy+dz*dz) );
    if(debug) cout << " Input point " << it << " Ltrace " << l_trace.back() << endl;
  }
//  Shower shower;
  vector<size_t> cells_shower;
  vector<double> l_shower;
  vector<double> point_shower[3];
  vector<double> nl_shower_Lprofile;
  vector<double> rl_shower_Lprofile;
  vector<double> e_shower_Lprofile;

  double deltal=10;
//  double ang = fabs(fParticle.GetAngleX());
//  double angy = fabs(fParticle.GetAngleX());
//  if( angy > ang ) ang = angy;
//  if( ang > 0.001 )  deltal = 5/ang;


  l_shower.push_back(deltal);
  if( debug ) cout << " Add start shower point at l= " << deltal << endl;
  point_shower[0].push_back(v_input_trace_x[0]);
  point_shower[1].push_back(v_input_trace_y[0]);
  point_shower[2].push_back(v_input_trace_z[0]);
//  int last=0;

  double dx = v_output_trace_x[npoints-1]-v_input_trace_x[0];
  double dy = v_output_trace_y[npoints-1]-v_input_trace_y[0];
  double dz = v_output_trace_z[npoints-1]-v_input_trace_z[0];
  double llast =sqrt(dx*dx+dy*dy+dz*dz);  // Total shower length in the detector calculate it from Trace
  if( debug ) cout  << " Total shower length in the detector = " << llast << endl;

  l_trace.push_back(llast); // add final point to l_trace
  if(debug) cout << " Final Output point " << " Ltrace " << l_trace.back() << endl;
  double dax = dx/llast;
  double day = dy/llast;
  double daz = dz/llast;
  double lnext = 0.;

  for( int it=0; it != 100000; it++ )
  {
    lnext = l_shower.back()+deltal;                   //
    double xnext = point_shower[0].back()+dax*deltal; // using straigt line propagate at next point
    double ynext = point_shower[1].back()+day*deltal;
    double znext = point_shower[2].back()+daz*deltal;

    if(lnext > llast)
    {
//      if( debug ) cout << " Shower cutting finisfed OK lnext = " << lnext << " llast = " << llast << endl;
      break;
    }
    l_shower.push_back(lnext);
//    if( debug ) cout << " Add new shower point at l= " << lnext << endl;
    point_shower[0].push_back(xnext);
    point_shower[1].push_back(ynext);
    point_shower[2].push_back(znext);
  }

//???  cells_shower.push_back(cells_trace[0]);
  for( int it=0; it != (int)l_shower.size(); it++ )
  {
    int last = -1;
    for( int it1=0; it1 != (int)l_trace.size(); it1++ )
    {
      if( l_shower[it] > l_trace[it1]) continue;
      last = it1-1;
      break;
    }
    if( last >= 0 && last < (int)cells_trace.size() )
    {
      cells_shower.push_back(cells_trace[last]);
//      if( debug ) cout << "  At shower point at l= " << l_shower[it] << " cell =" << cells_shower.back() << endl;
    }
    else
    {
      cerr << " Internal ERROR in OneParticleResponse::CalculateExpectedResponseForDevelopment:" <<
                                       " wrong index = " << last << " at l_shower = " << l_shower[it] << endl;
      exit(1);
    }
  }

  if( debug ) cout << " Shower cutting finisfed  npoints " << cells_shower.size() << " check " << l_shower.size() << endl;

  for( size_t it=0; it != cells_shower.size(); it++ )
  {
    const Cell* p_cell = &(calorimeter->GetCells()[cells_shower[it]]);
    double RadLeng = p_cell->GetCellType().GetRadiationLength();
    double NucLeng = p_cell->GetCellType().GetNuclearLength();
    double rli = l_shower[it]/RadLeng;
    double nli = l_shower[it]/NucLeng;
    rl_shower_Lprofile.push_back(rli);
    nl_shower_Lprofile.push_back(nli);
  }


  for( int i =0; i !=(int)hited_cells.size(); i++ ) //check code
  {
//     if( expected_one_particle_response_in_cells[i].first != 0) cout << " Error in initialisation of CalculateExpectedResponse " << endl;
    if( e_expected_[i] != 0) cout << " Error in initialisation of CalculateExpectedResponse " << endl;
  }


  if( fParticle.GetID() == CalorimeterParticle::NEUTRINO ) // NO interactions
  {
//    return;
  }
  else if( fParticle.GetID() == CalorimeterParticle::GAMMA || fParticle.GetID() == CalorimeterParticle::POSITRON ||
                                                       fParticle.GetID() == CalorimeterParticle::ELECTRON)     // Electromagnetic shower
  {
    double e_shower=0.;
//     double l0_start_shower = 15.;
    double l0_start_shower = 0.;
    if( fabs( l0_start_shower ) > 0.00001 ) cerr << " Danger GAMES with shower functions!!! l0_start_shower = " <<l0_start_shower << endl;
    double l_start_shower = l0_start_shower + 1.5*(calorimeter->GetRandomFlat()-0.5);
    for( size_t it=0; it != cells_shower.size(); it++ )
    {
      double rl = 0.;
      double nl = 0.;
      if(it != 0)
      {
        rl = rl_shower_Lprofile[it-1];
        nl = nl_shower_Lprofile[it-1];
      }
      double rli = rl_shower_Lprofile[it];
//      double nli = nl_shower_Lprofile[it];
      double e = 0;

      e = cascade_EM(fParticle.GetE(),rl+l_start_shower,rli+l_start_shower);
//      e = cascade_EM(fParticle.GetE(),rl,rli);

      e_shower += e;

      e_shower_Lprofile.push_back(e);

      const Cell*  p_main_cell = &(calorimeter->GetCells()[cells_shower[it]]);
      double p_main_RadLeng = p_main_cell->GetCellType().GetRadiationLength();
//      double p_main_NucLeng = p_main_cell->GetCellType().GetNuclearLength();
      double xs = point_shower[0][it];
      double ys = point_shower[1][it];
      double zs = point_shower[2][it];

      double e_check=0;

      for( int i =0; i !=(int)hited_cells.size(); i++ )
      {
        int jcell = hited_cells[i];
        const Cell*  pcell = &(calorimeter->GetCells()[jcell]);
        double xc = pcell->GetX();
        double yc = pcell->GetY();
        double zc = pcell->GetZ();
        double sxc = pcell->GetCellType().GetSizeX();
        double syc = pcell->GetCellType().GetSizeY();
        double szc = pcell->GetCellType().GetSizeZ();
        double energy_in_cell = 0;
        double fluctuations_energy_in_cell = 0;

        if( zs < zc+szc/2 && zs > zc-szc/2)
	{
          energy_in_cell = ShowerLednev( p_main_RadLeng,e,xs-xc,ys-yc,sxc,syc,0.,0.); // Fluctuations calculation ???

  /* Add cell responce function  jcell - cell's index in calorimeter (energy_in_cell,xs,ys,zs)-shower hit
     Calorimeter::CellResponse function violate OO and here just for development  */

          energy_in_cell = calorimeter->CellResponse(jcell,energy_in_cell,xs,ys,zs);
        }

        e_check += energy_in_cell;
// ??        double  fraction = 1. - energy_in_cell/(0.9*fParticle.GetE());
// ??         fraction = fraction*fraction;
// ??         double fluctuations_energy_in_cell = 0.1*energy_in_cell*fraction;

//         expected_one_particle_response_in_cells[i].first += energy_in_cell;
//         expected_one_particle_response_in_cells[i].second += fluctuations_energy_in_cell;
        e_expected_[i] += energy_in_cell;
        de_expected_[i] += fluctuations_energy_in_cell;
      }

      if( debug ) cout << " Gamma Lprofile rl=" << rl << " rli=" << rli << " e =" <<
                                      e  << " e_check =" << e_check/2 << " Z= " << point_shower[2][it] << endl;

    }
    if( debug ) cout << " EM fRarticle E = " << fParticle.GetE() << " E Shower = " << e_shower << endl;
  }
  else if( fParticle.GetID() == CalorimeterParticle::MUON_PLUS || fParticle.GetID() == CalorimeterParticle::MUON_MINUS )  // MIP
  {
     double e = cascade_MIP(fParticle.GetE(),0,14.5);
     if( debug ) cout << " MIP fRarticle deltaE " << e/fParticle.GetE() << endl;
  }
  else                                                                                              // Hadronic shower
  {
     double e = cascade_Hadron(fParticle.GetE(),0,14.5);
     if( debug ) cout << " Hadronic fRarticle deltaE " << e/fParticle.GetE() << endl;
  }

    double etotal_expected =0;
    for( int i =0; i !=(int)hited_cells.size(); i++ )
    {
//        int jcell = hited_cells[i];
//        const Cell*  pcell = &(calorimeter->GetCells()[jcell]);
//         etotal_expected += expected_one_particle_response_in_cells[i].first;
        etotal_expected += e_expected_[i];
//        cout << " cell " << i << " # " << jcell << " Ampl " << expected_one_particle_response_in_cells[i].first << endl;
    }
    fEnergy_deposit=etotal_expected;
    if( debug ) cout << " OneParticleResponse::CalculateExpectedResponse E particle " << fParticle.GetE() << " Etotal expected " << etotal_expected << endl;
}

////////////////////////////////////////////////////////////////////////////////

void OneParticleResponse::CheckPrint( void ) const
{
   cout << " OneParticleResponse::CheckPrint is dummy !!! " << endl;
}

////////////////////////////////////////////////////////////////////////////////

/*!
    \brief Calculate expected amplitudes for pcell in a cluster
    \param cell_indx {Index of pcell}
    \param pcell {Calculate amplitude of this cell}
    \param p_main_cell {Main cell of particle/cluster}
    \param with_derivatives {Also calculate derivatives}
    \callgraph
    \callergraph
    \remarks {
    Used options:
    flag_fast_online

    Hardcoded:
    bool print_errors : Enable/disable error loging
    double fr : Fractions in fast online reconstruction
    }
*/


void
OneParticleResponse::CalculateParticleResponseInCell(int cell_indx,
								  const Cell* &pcell,
								  const Cell* &p_main_cell,
                                  const vector<pair<double, double> >& realAmplitudes,
								  bool with_derivatives)
{

  bool print_errors = false;

  if( with_derivatives )
    cout << " OneParticleResponse::CalculateParticleResponseInCell:: Sorry cant calculate derivatives for a while ! " << endl;

  // Get radiation/nuckear interaction length
  double RadLeng = p_main_cell->GetCellType().GetRadiationLength();
  double NucLeng = p_main_cell->GetCellType().GetNuclearLength();


  double  energy_in_cell = 0.;
  double  fluctuations_energy_in_cell=0.;
  double  shower_fluctuations_energy_in_cell=0.;
  double  detector_response_fluctuations_energy_in_cell=0.;
  double  parameters_uncertainty_fluctuations_energy_in_cell=0.;

  // Get cell global index
  int jcell = hited_cells[cell_indx];

  // Check if energy is set correct
  if( fParticle.GetE()<=0 )
    throw Exception("OneParticleResponse::ParticleResponseInCell: negative energy: %g GeV",fParticle.GetE());

  // Do nothing here if it is an neutino
  if( fParticle.GetID() == CalorimeterParticle::NEUTRINO ) {} // NO interactions

  // Particles with em showers
  else if( fParticle.GetID() == CalorimeterParticle::GAMMA || fParticle.GetID() == CalorimeterParticle::POSITRON || fParticle.GetID() == CalorimeterParticle::ELECTRON) // Electromagnetic shower
  {

    // Fast online reconstruction
    if (calorimeter->GetOptions().flag_fast_online)
    {
      // # warning " TODO improve  fast_online option in ParticleResponseInCell "

      // Get particle/cell position diffrence
      double dx = fParticle.GetX()-pcell->GetX();
      double dy = fParticle.GetY()-pcell->GetY();

      // Get cell size
      double sx = pcell->GetCellType().GetSizeX();
      double sy = pcell->GetCellType().GetSizeY();

      // Get minimum defiance
      double dmin = fabs(dx);
      if(fabs(dy) < fabs(dx))dmin = fabs(dy);

      // Determin energy fraction in cell
      double fr;
      // Particle position in cell
      if( fabs(dx) < sx && fabs(dy) < sy)
        fr = 0.8;
      // Particle position less than 1.5 cells apart
      else if(dmin < 1.5*sx)
        fr = 0.1;
      // Particle position less than 2 cells apart
      else if(dmin < 2*sx)
        fr = 0.05;
      // Particle position more than 2 cells apart
      else
        fr = 0.;

      // calculate expected energy
      energy_in_cell = fr*fParticle.GetE();

      // Calculate fluctuation
      // ??????????????????????
      double  fraction = 1. - energy_in_cell/(0.9*fParticle.GetE());
      fraction = fraction*fraction;
      fluctuations_energy_in_cell = 0.1*energy_in_cell*fraction;

    }
    // Offline reconstruction
    else
    {

      // Using the algorithm described in 'Shower seperation program for ECAL2' A.A. Lednev
      // Modification to this algorthim
      energy_in_cell = ShowerLednev( RadLeng,
                                     fParticle.GetE(),
                                     fParticle.GetX()-pcell->GetX(),
                                     fParticle.GetY()-pcell->GetY(),
                                     pcell->GetCellType().GetSizeX(),
                                     pcell->GetCellType().GetSizeY(),
                                     0.,0.);


      // Checking for errors
      if( energy_in_cell < 0. )
        {
          if( print_errors )
            cout << " ERROR!!! OneParticleResponse::ParticleResponseInCell:: E < 0  " << energy_in_cell << endl;
          // Can you simply invert the sign ??????
          energy_in_cell = fabs(energy_in_cell);
        }


      detector_response_fluctuations_energy_in_cell = calorimeter->EnergyDispersionInCellBasic(energy_in_cell,jcell);

      // Calculate fluctuation
      // ???????????????????? see above
      double  fract = energy_in_cell/(0.9*fParticle.GetE());
      double  fraction = (1.-fract)*(1.-fract);
      shower_fluctuations_energy_in_cell = 0.1*energy_in_cell*fraction;

      // In case of showers overlap we use direct energy fluctuation calculation using dX, dY estimation
      if( fabs( realAmplitudes[jcell].first - e_real_[cell_indx] ) > 0.01 )
        {
          // ??????????? not checking sign here but above !!!!!!!!!!!!!!!

          // Calculate maximal/minimal energy within position errors
          double emin = energy_in_cell;
          double emax = energy_in_cell;
          double e1010 = ShowerLednev( RadLeng, fParticle.GetE(),
                                       fParticle.GetX()-pcell->GetX()-fParticle.GetXerr(),
                                       fParticle.GetY()-pcell->GetY()-fParticle.GetYerr(),
                                       pcell->GetCellType().GetSizeX(),
                                       pcell->GetCellType().GetSizeY(),
                                       fParticle.GetAngleX(), fParticle.GetAngleY());
          if(e1010 > emax ) emax = e1010;
          if(e1010 < emin ) emin = e1010;
          double e1001 = ShowerLednev( RadLeng, fParticle.GetE(),
                                       fParticle.GetX()-pcell->GetX()-fParticle.GetXerr(),
                                       fParticle.GetY()-pcell->GetY()+fParticle.GetYerr(),
                                       pcell->GetCellType().GetSizeX(),
                                       pcell->GetCellType().GetSizeY(),
                                       fParticle.GetAngleX(), fParticle.GetAngleY());
          if(e1001 > emax ) emax = e1001;
          if(e1001 < emin ) emin = e1001;
          double e0110 = ShowerLednev( RadLeng, fParticle.GetE(),
                                       fParticle.GetX()-pcell->GetX()+fParticle.GetXerr(),
                                       fParticle.GetY()-pcell->GetY()-fParticle.GetYerr(),
                                       pcell->GetCellType().GetSizeX(),
                                       pcell->GetCellType().GetSizeY(),
                                       fParticle.GetAngleX(), fParticle.GetAngleY());
          if(e0110 > emax ) emax = e0110;
          if(e0110 < emin ) emin = e0110;
          double e0101 = ShowerLednev( RadLeng, fParticle.GetE(),
                                       fParticle.GetX()-pcell->GetX()+fParticle.GetXerr(),
                                       fParticle.GetY()-pcell->GetY()+fParticle.GetYerr(),
                                       pcell->GetCellType().GetSizeX(),
                                       pcell->GetCellType().GetSizeY(),
                                       fParticle.GetAngleX(), fParticle.GetAngleY());
          if(e0101 > emax ) emax = e0101;
          if(e0101 < emin ) emin = e0101;

          // Calculate possible fluctuation error
          double smima = (emax-emin)/2.;

          // Set fluctuation uncertainty
          parameters_uncertainty_fluctuations_energy_in_cell = smima*smima;

        }
      else
        {
          // Set default zero
          parameters_uncertainty_fluctuations_energy_in_cell = 0.;
        }

      //TODO   This formula for fluctuations_energy_in_cell is just to write something here.

      //TODO       fluctuations_energy_in_cell += pcell->GetCellType().Disp_phe_per_1GeV()*energy_in_cell;
      //   Add dispersion of photo-electron statistic.

      fluctuations_energy_in_cell = detector_response_fluctuations_energy_in_cell +
        shower_fluctuations_energy_in_cell +
        parameters_uncertainty_fluctuations_energy_in_cell;

    }

  }
  else if( fParticle.GetID() == CalorimeterParticle::MUON_PLUS || fParticle.GetID() == CalorimeterParticle::MUON_MINUS )                     // MIP
  {

    // # warning " TODO improve  muon response in ParticleResponseInCell "
    if(p_main_cell == pcell )
      {
        // !!!!!!!!!!!!!!!! hardcoded values for MIPs, where do they come from ?????????????
        energy_in_cell = 2.5;
        fluctuations_energy_in_cell = 0.5;
      }
    else
      {
        // MIPs only react in main cell
        energy_in_cell = 0.;
        fluctuations_energy_in_cell = 0.1;
      }
  }
  else                                                                        // Hadronic shower
  {

    // # warning " TODO improve  hadron response in ParticleResponseInCell "

    // Check interaction length
    if( NucLeng <= 0.0001 )
      {
        cerr << " Bug in  NucLeng settings !!! " << NucLeng << " Set it to 50." << endl;
        NucLeng = 50.;
      }

    // Using the algorithm described in 'Shower seperation program for ECAL2' A.A. Lednev
    // Modification to this algorthim
    energy_in_cell = ShowerLednev( NucLeng,
                                   fParticle.GetE(),
                                   fParticle.GetX()-pcell->GetX(),
                                   fParticle.GetY()-pcell->GetY(),
                                   pcell->GetCellType().GetSizeX(),
                                   pcell->GetCellType().GetSizeY(),
                                   0.,0.);
    // Checking energy
    if( energy_in_cell < 0. )
      {
        // ??????????? Just inverting value
        // ??????????? No output !!!!!!!!!!!!!
        energy_in_cell = fabs(energy_in_cell);
      }

    // Calculate fluctuation
    // ???????????????????? see above
    double  fraction = 1. - energy_in_cell/(0.9*fParticle.GetE());
    fraction = fraction*fraction;
    fluctuations_energy_in_cell = 0.1*energy_in_cell*fraction;

    //TODO   This formula for fluctuations_energy_in_cell is just to write something here.

    //TODO       fluctuations_energy_in_cell += pcell->GetCellType().Disp_phe_per_1GeV()*energy_in_cell;
    //   Add dispersion of photo-electron statistic.

  }

  // ?????????? why don't we use this common structures from the beginning

  // setting values in common structure
  e_expected_[cell_indx] = energy_in_cell;
  de_expected_[cell_indx] = fluctuations_energy_in_cell;

  // storing used fluctation in common structure
  shower_fluctuations_one_particle_response[cell_indx].first = shower_fluctuations_energy_in_cell;
  shower_fluctuations_one_particle_response[cell_indx].second = parameters_uncertainty_fluctuations_energy_in_cell;


}


////////////////////////////////////////////////////////////////////////////////

void OneParticleResponse::CalculateDerivativesExpectedResponse(void)
{

  // Check that there are cells in cluster
  if( hited_cells.empty() )
    return;

  // Calculate expected responses
  CalculateExpectedResponse();

  // Get main cell
  const Cell*  p_main_cell = &(calorimeter->GetCells()[hited_cells.front()]);

  // Get radiation length of main cell
  double RadLeng = p_main_cell->GetCellType().GetRadiationLength();

  // Iterate over cells in cluster
  for( int i =0; i !=(int)hited_cells.size(); i++ ) {

    // Get global cell index
    int jcell = hited_cells[i];

    // Get cell
    const Cell*  pcell = &(calorimeter->GetCells()[jcell]);

    // Calculate fraction of total particle energy
    // deposit in cell
    de_one_particle_response_in_cells[i] = e_expected_[i] / fParticle.GetE();

    // Fixed parameter variation !!!!!
    double dx = 0.005;

    // Calculate expected energy fluctuation in cell for
    // an cluster deviation of 0.005 in x and y direction
    // !!!!! Just using parameter variation to calculate derivative
    // ????? Maybe not only check positive variation
    dx_one_particle_response_in_cells[i] =
      (ShowerLednev(RadLeng,fParticle.GetE(),
                    fParticle.GetX()-pcell->GetX()+dx,
                    fParticle.GetY()-pcell->GetY(),
                    pcell->GetCellType().GetSizeX(),
                    pcell->GetCellType().GetSizeY(),
                    0.,0.)
       -ShowerLednev(RadLeng,fParticle.GetE(),
                     fParticle.GetX()-pcell->GetX(),
                     fParticle.GetY()-pcell->GetY(),
                     pcell->GetCellType().GetSizeX(),
                     pcell->GetCellType().GetSizeY(),
                     0.,0.)) / dx;

    dy_one_particle_response_in_cells[i] =
      (ShowerLednev(RadLeng,fParticle.GetE(),
                    fParticle.GetX()-pcell->GetX(),
                    fParticle.GetY()-pcell->GetY()+dx,
                    pcell->GetCellType().GetSizeX(),
                    pcell->GetCellType().GetSizeY(),
                    0.,0.)
       -ShowerLednev(RadLeng,fParticle.GetE(),
                     fParticle.GetX()-pcell->GetX(),
                     fParticle.GetY()-pcell->GetY(),
                     pcell->GetCellType().GetSizeX(),
                     pcell->GetCellType().GetSizeY(),
                     0.,0.)) / dx;
    }
}

////////////////////////////////////////////////////////////////////////////////

double OneParticleResponse::GetDDE( int cell)
{
  if( cell < 0 || cell >= (int)calorimeter->NCells() ) return 0.;
  if( de_one_particle_response_in_cells.size() <= 0 ) return 0.;
  for( size_t i=0; i< hited_cells.size(); i++ )
  {
    if( cell == (int)hited_cells[i] )
    {
      return de_one_particle_response_in_cells[i];
    }
  }
  return 0.;
}

////////////////////////////////////////////////////////////////////////////////

double OneParticleResponse::GetDDX( int cell)
{
  if( cell < 0 || cell >= (int)calorimeter->NCells() ) return 0.;
  if( dx_one_particle_response_in_cells.size() <= 0 ) return 0.;
  for( size_t i=0; i< hited_cells.size(); i++ )
  {
    if( cell == (int)hited_cells[i] )
    {
      return dx_one_particle_response_in_cells[i];
    }
  }
  return 0.;
}

////////////////////////////////////////////////////////////////////////////////

double OneParticleResponse::GetDDY( int cell)
{
  if( cell < 0 || cell >= (int)calorimeter->NCells() ) return 0.;
  if( dy_one_particle_response_in_cells.size() <= 0 ) return 0.;
  for( size_t i=0; i< hited_cells.size(); i++ )
  {
    if( cell == (int)hited_cells[i] )
    {
      return dy_one_particle_response_in_cells[i];
    }
  }
  return 0.;
}

////////////////////////////////////////////////////////////////////////////////

void OneParticleResponse::UpdateParticleParameters( double e, double se, double x, double sx, double y, double sy)
{
  fParticle.SetE( e , se );
  fParticle.SetX( x , sx );
  fParticle.SetY( y , sy );
}

////////////////////////////////////////////////////////////////////////////////

/*!
  \brief Calculate real response of particle in cell
  \param problems {Vector of cells with problems, in case of problems current cell added to vector}
  \return {true on success}
  \remarks{
  Hardcoded:
  bool debug : Debug output
  bool split_proportional : Enable/disable amplitude sharing
  }
 */

bool OneParticleResponse::CalculateRealResponse( vector < size_t > &problems,
                                                const vector<pair<double, double> >& realAmplitudes,
                                                const vector<pair<double, double> >& expectedAmplitudes)

{
  bool debug = false;

  bool ok = true;

  bool split_proportional = true;

  e_overlap_ = 0.;

  // Do only work if there are cells in cluster
  if( hited_cells.empty() )
    return ok;

  // loop over cells in cluster
  for( int i =0; i !=(int)hited_cells.size(); i++ )
    {

      // Get global cell index
      int jcell = hited_cells[i];


      if( split_proportional || realAmplitudes[jcell].first < 0.1 )
        {
          // Calculate real energy deposit derived from expected energy deposit
          e_real_[i] = e_expected_[i]*realAmplitudes[jcell].first/
            expectedAmplitudes[jcell].first;

          if( debug ) cout << " Got Ecell " << e_real_[i]  <<
                        " It was in total ErWall " <<  realAmplitudes[jcell].first <<
                        " et " <<  expectedAmplitudes[jcell].first <<
                        " wanted " <<  e_expected_[i] << endl;

        }
      else
        {
          // Real amplitude
          double ertot = realAmplitudes[jcell].first;
          // Expected amplitude
          double ettot = expectedAmplitudes[jcell].first;
          // Error of expected amplitude
          double dispettot = expectedAmplitudes[jcell].second;
          // Error of real amplitude derived for expected error
          double der = (ertot - ettot)/dispettot;

          // Calculate real energydeposit
          e_real_[i] = ( 1. + der ) * de_expected_[i];

          // There mustnot be negative energies
          if( e_real_[i] < 0. )
            {
              // Set enrgy to zero
              e_real_[i] = 0.;
              // Mark problem
              problems.push_back( jcell );
              ok = false;

            }

          // Real energy should not be smaller than real amplitude, Why???????????????????
          double dereal = realAmplitudes[jcell].first - e_real_[i];
          if( dereal < 0. )
            {
              // Set energy to amplitude
              e_real_[i] = realAmplitudes[jcell].first;
              // Mark problem (!!!! add assertion that cell is not marked twice ?????)
              problems.push_back( jcell );
              ok = false;
            }

       }

      double ec = e_real_[i];

      // Weights for error calculation (!!!! migrade to option)
      double weight_base = 1.; // Fluctuation in cell basis
      double weight_ro   = 1.; // Fluctuation in readout
      double weight_dig  = 1.; // Error by digitalisation

      // Calculate dispersion from bases and readout
      double e_disp = weight_base * calorimeter->EnergyDispersionInCellBasic(ec, jcell)
                      + weight_ro * calorimeter->EnergyDispersionInCellReadout(ec, jcell);

      // Rescale dispersion of bases and readout
      e_disp = calorimeter->ScaleDispersionInSparceMode(ec, e_disp, jcell);

      // Add error from digitalization
      e_disp += weight_dig*calorimeter->EnergyDispersionInCellDigitization(ec, jcell);

      // Calculate energy error
      de_real_[i]  = e_disp + shower_fluctuations_one_particle_response[i].second;
      if( fabs( expectedAmplitudes[jcell].second - de_expected_[i] ) > 0.0001 )
        de_real_[i] += expectedAmplitudes[jcell].second - de_expected_[i];


      // Energy error to e_overlap_ if larger than 0.001 (!!!! migrade to option)
      double dereal = realAmplitudes[jcell].first - e_real_[i];
      if( dereal > 0.001 ) e_overlap_ += dereal;

      // Just error checking/output !!!!!
      // TODO: Add reaktion on error
      if( dereal < -0.001 || e_real_[i] < -0.001)
      {
        const Cell &pcell = calorimeter->GetCells()[jcell];
        cerr << " Error in CalculateRealResponse " << calorimeter->GetName() <<" cell " << jcell << " dereal " << dereal << endl;
	   cout << " Egamma       " <<  fParticle.GetE() << endl;
	   cout << " Xgamma-Xcell " <<  fParticle.GetX()-pcell.GetX() << endl;
	   cout << " Ygamma-Ycell " <<  fParticle.GetY()-pcell.GetY() << endl;
	   cout << " Cell X Size  " <<  pcell.GetCellType().GetSizeX() << endl;
	   cout << " Cell Y Size  " <<  pcell.GetCellType().GetSizeY() << endl;
	   cout << " Xgamma Sigma " <<  fParticle.GetXerr() << endl;
	   cout << " Ygamma Sigma " <<  fParticle.GetYerr() << endl;
	   cout << " Ecell " << realAmplitudes[jcell].first <<
	                " e " << e_real_[i] <<
	                " Se " << sqrt(de_real_[i]) <<
	                " et " << e_expected_[i] <<
	                " Set " << sqrt(de_expected_[i]) << endl;

	  cout << "  Diff Disp " << expectedAmplitudes[jcell].second - de_expected_[i] <<
	           " Cell DispTot " <<  expectedAmplitudes[jcell].second  <<
	             " Cell Disp " << de_expected_[i]  << endl;
          cout << " E measured in cell  " << jcell << " E="  <<realAmplitudes[jcell].first
                                      << " E wanted " <<  e_expected_[i]
                                      << " expected from all " <<  expectedAmplitudes[jcell].first << endl;

         double ertot = realAmplitudes[jcell].first;
         double ettot = expectedAmplitudes[jcell].first;
         double dispettot = expectedAmplitudes[jcell].second;
         double der = (ertot - ettot)/dispettot;

         cout << " ertot " << ertot << " ettot " << ettot <<  " dispettot " << dispettot  << " der " << der << " result " <<
               e_expected_[i] + der*de_expected_[i] << endl;


      }

    }
  return ok;
}

////////////////////////////////////////////////////////////////////////////////

void OneParticleResponse::ChangeMain( void )
{

  double emax = 0;

  int imax = -1;

  // Search cell with maximal amplitude
  for( int i =0; i !=(int)hited_cells.size(); i++ )
  {
    if( e_real_[i] > emax )
    {
      emax = e_real_[i];
      imax = i;
    }
  }

  int jcell = hited_cells[imax];
  // ... tak razuvat'sya
  vector<size_t> mmm;
  mmm.push_back(jcell);

  // Set new main cell
  fParticle.SetMainCells(mmm);

  // Reset particles xy position
  double x0 = calorimeter->GetCells()[jcell].GetX();
  double y0 = calorimeter->GetCells()[jcell].GetY();
  double sizeX = calorimeter->GetCells()[jcell].GetCellType().GetSizeX();
  double sizeY = calorimeter->GetCells()[jcell].GetCellType().GetSizeY();
  double sigmax = sizeX/2.;
  double sigmay = sizeY/2.;
  fParticle.SetX( x0 , sigmax );
  fParticle.SetY( y0 , sigmay );

}

////////////////////////////////////////////////////////////////////////////////

void Calorimeter::InitPreshowerCorrections( void )
{
  if( options.preshower_length > 0. )
  {
    if( preshower_corrections_.size() == 0 )
    {
      for( size_t i=0; i<3000; i++ )
      {
        double e = 0.05 + double(i)*0.1;
        double cr = Reco::cascade_EM( e, 0., options.preshower_length)/e;
        preshower_corrections_.push_back( cr );
      }
    }
  }
}

////////////////////////////////////////////////////////////////////////////////

void OneParticleResponse::FindList(void)
{
  bool debug=false;

  if(hited_cells.size() > 0 )
  {
    cerr << " Call FindList(double dist) in second time ?! This is not expected " << endl;
    exit(1);
  }

  // Reset particle information
  ClearList();

  // Get main cell index
  int p_main_cell = calorimeter->FindCell(fParticle);
  // Return if failed
  if( p_main_cell==-1 )
    return;
  // Get main cell
  const Cell& cell  = calorimeter->GetCells()[p_main_cell];

  double dist_max = calorimeter->GetOptions().distmax_for_cells_search;
  if( dist_max <= 0 )
  {
    // Get neighboring cells if distmax not set
    hited_cells = cell.GetNeighbors();
  }
  else
  {
    // If distmax set add all cells with distance < distmax
    hited_cells.push_back(p_main_cell);
    for( int i=0; i< (int)calorimeter->NCells(); i++)
    {
      double xi[3],xo[3];
      double dist = TraceInCell(calorimeter->GetCells()[i],xi,xo);
      if( i != p_main_cell && dist < dist_max )
        hited_cells.push_back(i);
    }
  }

  // Generate vectors of the right size
  // Set fMaskCells (distances of hited_cells and cell)
  InitList(cell);

  // Loop over cluster cells to set traces
  for( vector<size_t>::iterator it=hited_cells.begin();
       it != hited_cells.end(); it++ )
  {
    // Get distance and trace
    double xi[3],xo[3];
    double dist = TraceInCell(calorimeter->GetCells()[*it],xi,xo);

    if(debug)
      cout << " TraceInCell cell " << *it << " dist = " << dist << endl;

    // Store trace if track inside this cell
    if( dist < 0)  // track is inside the cell
    {
      cells_trace.push_back(*it);
      v_input_trace_x.push_back(xi[0]);
      v_input_trace_y.push_back(xi[1]);
      v_input_trace_z.push_back(xi[2]);
      v_output_trace_x.push_back(xo[0]);
      v_output_trace_y.push_back(xo[1]);
      v_output_trace_z.push_back(xo[2]);
    }
  }



  if(debug)
  {
    cout << " OneParticleResponse:ExtendList(double) printout NCELLs=" << cells_trace.size() << endl;
    for( size_t it=0; it != cells_trace.size(); it++ )
    {
      cout << " Cell " << cells_trace[it] << " x_in " << v_input_trace_x[it] <<" y_in " << v_input_trace_y[it] <<
                                                                               " z_in " << v_input_trace_z[it] <<
                                        " x_out " <<  v_output_trace_x[it] <<" y_out " << v_output_trace_y[it] <<
				                                        " z_out " << v_output_trace_z[it] << endl;
    }
  }


  SortTrace();

  if( debug )  PrintTrace();

}

////////////////////////////////////////////////////////////////////////////////

void OneParticleResponse::FindList(const Cell& cell)
{
  // Reset particle information
  ClearList();

  // Add neighbors to cluster, which forms particle
  hited_cells = cell.GetNeighbors();

  // Generate vectors of the right size
  // Set fMaskCells (distances of hited_cells and cell)
  InitList(cell);

}


////////////////////////////////////////////////////////////////////////////////
// Reset the cluster information for particle                                 //
////////////////////////////////////////////////////////////////////////////////

void OneParticleResponse::ClearList(void)
{

  // empty the vector of cells in cluster (index)
  hited_cells.clear();

  // clear the energy deposition expectation of cells in cluster
  e_expected_.clear();
  de_expected_.clear();

  // clear the energy deposition of cells in cluster
  e_real_.clear();
  de_real_.clear();

  //
  shower_fluctuations_one_particle_response.clear();

  // clear mask for neighbour search
  fMaskCells.clear();

  // clear derivative vectors
  de_one_particle_response_in_cells.clear();
  dx_one_particle_response_in_cells.clear();
  dy_one_particle_response_in_cells.clear();

}

////////////////////////////////////////////////////////////////////////////////

void OneParticleResponse::InitList(const Cell &cell)
{

  // loop over cell indices contributing to this particle
  for( vector<size_t>::iterator it=hited_cells.begin(); it != hited_cells.end(); it++ )
  {

    // Get distance (dx,dy) between cell and cell
    // with index *it in units of the dimensions of cell
    fMaskCells.push_back(cell.NeighborMask(calorimeter->GetCells()[*it],
                                            calorimeter->GetOptions().tolerance_for_nearby_cells));

    // Set the expectation of energy contribution to (dummy) zero (for this cell)
    e_expected_   .push_back(0.);
    de_expected_  .push_back(0.);

    // Set the energy contribution to (dummy) zero (for this cell)
    e_real_       .push_back(0.);
    de_real_       .push_back(0.);

    // Set dummy for shower fluctuation
    shower_fluctuations_one_particle_response.push_back( pair<double,double>(0,0) );

    // Set dummy for derivatives
    de_one_particle_response_in_cells.push_back(0);
    dx_one_particle_response_in_cells.push_back(0);
    dy_one_particle_response_in_cells.push_back(0);
  }
}


///////////////////////////////////////////////////////////////////////////////
// Extend vector of hited_cells by adding neighbors of the cells which are   //
// already part of the cluster                                               //
///////////////////////////////////////////////////////////////////////////////

void OneParticleResponse::ExtendList(void)
{
  if( hited_cells.empty() )
    return;

  // need to make a copy because loop vector cannot be modified
  vector<size_t> dummy = hited_cells;

  // loop over cluster cells
  /// \todo: why is the first cell skipped?
  for( vector<size_t>::iterator it=dummy.begin()+1; it != dummy.end(); it++ )
  {
    for( vector<size_t>::const_iterator it1= calorimeter->GetCells()[*it].GetNeighbors().begin();
	 it1 != calorimeter->GetCells()[*it].GetNeighbors().end(); it1++ )
    {
      // Continue if current neighbor allready member of cluster
      if( std::find(hited_cells.begin(),hited_cells.end(),*it1) != hited_cells.end() )
        continue;

      // Add new cell to hited_cells and extend vectors
      hited_cells.push_back(*it1);
      fMaskCells.push_back( calorimeter->GetCells()[hited_cells[0]].NeighborMask(calorimeter->GetCells()[*it1],
                                                                            calorimeter->GetOptions().tolerance_for_nearby_cells));
      e_expected_   .push_back(0.);
      de_expected_  .push_back(0.);
      e_real_       .push_back(0.);
      de_real_      .push_back(0.);
      shower_fluctuations_one_particle_response .push_back( pair<double,double>(0,0) );
      de_one_particle_response_in_cells.push_back(0);
      dx_one_particle_response_in_cells.push_back(0);
      dy_one_particle_response_in_cells.push_back(0);
    }
  }
}


///////////////////////////////////////////////////////////////////////////////
// Recursively extend vector of hited_cells up to distmax from the track by  //
// adding neighbors of the cells which are already part of the cluster       //
///////////////////////////////////////////////////////////////////////////////

void OneParticleResponse::ExtendList(double distmax)
{
//  bool debug=true;
  bool debug=false;
  if(debug) cout << " OneParticleResponse:ExtendList(double) in " << calorimeter->GetName() << endl;
  if( hited_cells.empty() )
    return;

  // need to make a copy because loop vector cannot be modified
  vector<size_t> near = hited_cells;
  double xi[3],xo[3];

  for( vector<size_t>::iterator it=near.begin(); it != near.end(); it++ )
  {
    for( vector<size_t>::const_iterator it1= calorimeter->GetCells()[*it].GetNeighbors().begin();
	 it1 != calorimeter->GetCells()[*it].GetNeighbors().end(); it1++ )
    {
      // Continue if current neighbor allready member of cluster
      if( std::find(hited_cells.begin(),hited_cells.end(),*it1) != hited_cells.end() )
	continue;

      double dist = TraceInCell(calorimeter->GetCells()[*it1],xi,xo);
      if(debug) cout << " TraceInCell cell " << *it1 << " dist = " << dist << endl;

      if( dist > distmax ) continue;

      if( dist < 0)  // track is inside the cell
      {
        cells_trace.push_back(*it1);
        v_input_trace_x.push_back(xi[0]);
        v_input_trace_y.push_back(xi[1]);
        v_input_trace_z.push_back(xi[2]);
        v_output_trace_x.push_back(xo[0]);
        v_output_trace_y.push_back(xo[1]);
        v_output_trace_z.push_back(xo[2]);
      }

    // Add new cell to hited_cells
      hited_cells.push_back(*it1);
      fMaskCells.push_back( calorimeter->GetCells()[hited_cells[0]].NeighborMask(calorimeter->GetCells()[*it1],
                                                                              calorimeter->GetOptions().tolerance_for_nearby_cells));
      e_expected_   .push_back(0.);
      de_expected_  .push_back(0.);
      e_real_       .push_back(0.);
      de_real_      .push_back(0.);
      shower_fluctuations_one_particle_response .push_back( pair<double,double>(0,0) );
      de_one_particle_response_in_cells.push_back(0);
      dx_one_particle_response_in_cells.push_back(0);
      dy_one_particle_response_in_cells.push_back(0);
    }
  } // end for it

  if(hited_cells.size() != near.size()) ExtendList(distmax);
  if(debug) cout << " OneParticleResponse:ExtendList(double) hited_cells " << hited_cells.size() << endl;
  if(debug)
  {
    cout << " OneParticleResponse:ExtendList(double) printout NCELLs=" << cells_trace.size() << endl;
    for( size_t it=0; it != cells_trace.size(); it++ )
    {
      cout << " Cell " << cells_trace[it] << " x_in " << v_input_trace_x[it] <<
                                             " y_in " << v_input_trace_y[it] <<
					     " z_in " << v_input_trace_z[it] <<
                                             " x_out " << v_output_trace_x[it] <<
					     " y_out " << v_output_trace_y[it] <<
					     " z_out " << v_output_trace_z[it] << endl;
    }
  }
}

///////////////////////////////////////////////////////////////////////////////

void OneParticleResponse::CalculateChi2(void)
{
  double esum_real = 0.;
  double esum_expected = 0.;
  for( size_t it=0; it < hited_cells.size(); it++)
  {
//     esum_real += real_one_particle_response_in_cells[it].first;
//     esum_expected += expected_one_particle_response_in_cells[it].first;
    esum_real += e_real_[it];
    esum_expected += e_expected_[it];
  }
  double correction = esum_real/esum_expected;


  double chi2 = 0;
  int ndf = 0;
//   for( size_t it=0; it < expected_one_particle_response_in_cells.size(); it++)
  for( size_t it=0; it < e_expected_.size(); it++)
  {
//     pair <double,double> e_real = real_one_particle_response_in_cells[it];
//     pair <double,double> e_expected = expected_one_particle_response_in_cells[it];
//     e_expected.first = e_expected.first*correction;
//     double de = e_real.first - e_expected.first;
//     double we = 1./(e_real.second+e_expected.second);

    double de = e_real_[it] - e_expected_[it]*correction;
    double we = 1./(de_real_[it]+de_expected_[it]);
    double hi_add = de*de*we;
    chi2 += hi_add;
    ndf++;
  }
  fNDF = ndf;
  fXiSqr = chi2;
}

////////////////////////////////////////////////////////////////////////////////

double OneParticleResponse::TraceInCell (const Cell & pcell ,double *x_in,double *x_out)
{

  bool debug=false;

  // Measurment precision !!!!! assumed hardcoded
  double precision=0.000001;

  // Get Cell precision
  double xc = pcell.GetX();
  double yc = pcell.GetY();
  double zc = pcell.GetZ();

  // Get Cell size
  double sx = pcell.GetCellType().GetSizeX();
  double sy = pcell.GetCellType().GetSizeY();
  double sz = pcell.GetCellType().GetSizeZ();

  if(debug) {
    cout << " OneParticleResponse::TraceInCell Detector " << calorimeter->GetName()
         << " xc " << xc << " yc " << yc << " zc " << zc << " Ec " << endl
         << " ************************************************* " << endl
         << " Particle   E= " << fParticle.GetE() << " X=" << fParticle.GetX() << " Y=" << fParticle.GetY()
         << "  axp= " << fParticle.GetAngleX() << " ayp=" << fParticle.GetAngleY() << endl;

  }

  // GoTo Cell's RS.
  double xp = fParticle.GetX()-xc;
  double yp = fParticle.GetY()-yc;
  double zp = fParticle.GetZ()-zc;
  double axp = fParticle.GetAngleX();
  double ayp = fParticle.GetAngleY();

  double z_front = -sz/2;
  double z_back = +sz/2;

  // Calculate track position (x,y) at begiing/end of cell
  double x_front = axp*(z_front-zp) + xp;
  double y_front = ayp*(z_front-zp) + yp;
  double x_back  = axp*(z_back-zp)  + xp;
  double y_back  = ayp*(z_back-zp)  + yp;

  // Add deviation in x if outside precision
  double dax =0;
  if(fabs(y_back-y_front) > precision )
    dax= (x_back-x_front)/(y_back-y_front);

  // Add deviation in y if outside precision
  double day =0;
  if(fabs(x_back-x_front) > precision )
    day= (y_back-y_front)/(x_back-x_front);

  if(debug)
    cout << " x_front=" << x_front << " y_front=" << y_front <<
      " x_back=" << x_back << " y_back=" << y_back <<
      " dx/dy=" << dax << " dy/dx=" << day << endl;



  double xhit[3][10];
  int indx = -1;

  // Add hits with trace hit front of cell
  if( (x_front < sx/2+precision) && (x_front > -sx/2-precision) &&
      (y_front < sy/2+precision) && (y_front > -sy/2-precision) )
    {
      if(debug)
        cout << " OneParticleResponse::TraceInCell Hit front plane of cell " << endl
             << " At X=" << x_front+xc << " Y=" << y_front+yc << " Back At X=" << x_back+xc << " Back Y=" << y_back+yc << endl;

      // TODO: Check index to avoid out of range access
      indx++;
      xhit[0][indx]=x_front;
      xhit[1][indx]=y_front;
      xhit[2][indx]=z_front;

      if(debug) cout << " hit at front side z= " << xhit[2][indx]+zc << endl;
    }

  // Add hits with trace hit back of cell
  if( (x_back < sx/2+precision) && (x_back > -sx/2-precision) &&
      (y_back < sy/2+precision) && (y_back > -sy/2-precision) )
    {
      if(debug)
        cout << " OneParticleResponse::TraceInCell Hit back plane of cell " << endl
             << " At X=" << x_back+xc << " Y=" << y_back+yc << " Enter At X=" << x_front+xc << " Enter Y=" << y_front+yc <<endl;

      // TODO: Check index to avoid out of range access
      indx++;
      xhit[0][indx]=x_back;
      xhit[1][indx]=y_back;
      xhit[2][indx]=z_back;

      if(debug) cout << " hit at back side z= " << xhit[2][indx]+zc << endl;
    }

  // Try to recover from not found trace by extending the search conditions
  if(indx != 1) {

     double y_up = sy/2;
     double x_up = x_front + dax*(y_up-y_front);
     double z_up = z_front;
     if( fabs(ayp) > 0.00001 ) z_up += (y_up-y_front)/ayp;

     double y_down = -sy/2;
     double x_down = x_front + dax*(y_down-y_front);
     double z_down = z_front;
     if( fabs(ayp) > 0.00001 ) z_down += (y_down-y_front)/ayp;

     double x_right = sx/2;
     double y_right = y_front + day*(x_right-x_front);
     double z_right = z_front;
     if( fabs(axp) > 0.00001 ) z_right += (x_right-x_front)/axp;

     double x_left = -sx/2;
     double y_left = y_front + day*(x_left-x_front);
     double z_left = z_front;
     if( fabs(axp) > 0.00001 ) z_left += (x_left-x_front)/axp;

     if( ((x_up < sx/2+precision) && (x_up > -sx/2-precision)) &&
         ((x_front-x_up)*(x_back-x_up) < 0 && (y_front-y_up)*(y_back-y_up) < 0 ) )
       {
         if(debug) cout << " hit at top side x_up= " << x_up << " y_up=" << y_up << endl;
         indx++;
         if( indx > 1)
           {
             if(debug) cout << " Internal error ??? in OneParticleResponse::TraceInCell indx= " << indx << endl;
             //          exit(1);
           }
         xhit[0][indx]=x_up;
         xhit[1][indx]=y_up;
         xhit[2][indx]=z_up;
         if(debug) cout << " hit at top side z= " << xhit[2][indx]+zc << endl;
       }

     if( ((x_down < sx/2+precision) && (x_down > -sx/2-precision)) &&
         ((x_front-x_down)*(x_back-x_down) < 0 && (y_front-y_down)*(y_back-y_down) < 0 ) )
       {
         if(debug) cout << " hit at down side x_down= " << x_down << " y_down=" << y_down << endl;
         indx++;
         if( indx > 1)
           {
             if(debug) cout << " Internal error ??? in OneParticleResponse::TraceInCell indx= " << indx << endl;
             //          exit(1);
           }
         xhit[0][indx]=x_down;
         xhit[1][indx]=y_down;
         xhit[2][indx]=z_down;
         if(debug) cout << " hit at bottom side z= " << xhit[2][indx]+zc << endl;
       }

     if( ((y_left < sy/2+precision) && (y_left > -sy/2-precision)) &&
         ((x_front-x_left)*(x_back-x_left) < 0 && (y_front-y_left)*(y_back-y_left) < 0 ) )
       {
         if(debug) cout << " hit at left side x_left= " << x_left << " y_left=" << y_left << endl;
         indx++;
         if( indx > 1)
            {
              if(debug) cout << " Internal error ??? in OneParticleResponse::TraceInCell indx= " << indx << endl;
              //          exit(1);
            }
         xhit[0][indx]=x_left;
         xhit[1][indx]=y_left;
         xhit[2][indx]=z_left;
         if(debug) cout << " hit at left side z= " << xhit[2][indx]+zc << endl;
       }

     if( ((y_right < sy/2+precision) && (y_right > -sy/2-precision)) &&
         ((x_front-x_right)*(x_back-x_right) < 0 && (y_front-y_right)*(y_back-y_right) < 0 ) )
       {
         if(debug) cout << " hit at right side x_right= " << x_right << " y_right=" << y_right << endl;
         indx++;
         if( indx > 1)
           {
             if(debug) cout << " Internal error ??? in OneParticleResponse::TraceInCell indx= " << indx << endl;
             //          exit(1);
           }
         xhit[0][indx]=x_right;
         xhit[1][indx]=y_right;
         xhit[2][indx]=z_right;
         if(debug) cout << " hit at right side z= " << xhit[2][indx]+zc << endl;
       }
     if(debug) cout << " After u,bottom,left,right check indx= " << indx << endl;
  }


  if(indx == 1)
    {

      // Sort in/outgoing trace into structure (using z pos)
      int indx_first=0;
      int indx_last=1;
      if(xhit[2][indx_first] > xhit[2][indx_last])
      {
        indx_first=1;
        indx_last=0;
      }
      x_in[0]=xhit[0][indx_first];
      x_in[1]=xhit[1][indx_first];
      x_in[2]=xhit[2][indx_first];
      x_out[0]=xhit[0][indx_last];
      x_out[1]=xhit[1][indx_last];
      x_out[2]=xhit[2][indx_last];

      // Get minimal distance
      double dist = sx/2-fabs(x_in[0]+x_out[0])/2;
      double disty = sy/2-fabs(x_in[1]+x_out[1])/2;
      if(dist > disty ) dist=disty;

      // One should take reaktion on an internal error !!!!!!!!!!!!!
      if(dist < 0)
        cout << " Internal error in OneParticleResponse::TraceInCell dist_inside= " << dist << endl;

      // Get global position
      x_in[0] +=xc;
      x_in[1] +=yc;
      x_in[2] +=zc;
      x_out[0] +=xc;
      x_out[1] +=yc;
      x_out[2] +=zc;

      // Return distance (negative ???????)
      return -dist;

    }


  if(indx == 0)
    {
      cerr << " Internal error in OneParticleResponse::TraceInCell Finaly indx= " << indx << endl;
      cerr << " OneParticleResponse::TraceInCell Detector " << calorimeter->GetName() <<
        " E "  << fParticle.GetE() << " xc " << xc << " yc " << yc << " zc " << zc << endl;
      return precision;
      //      exit(1);
    }

  if(indx == -1)
  {
    double p_in[2],p_out[2];
    double x1[2],x2[2],x3[2],x4[2];
    x1[0]=-sx/2;
    x1[1]=-sy/2;

    x2[0]=-sx/2;
    x2[1]= sy/2;

    x3[0]= sx/2;
    x3[1]= sy/2;

    x4[0]= sx/2;
    x4[1]=-sy/2;

    p_in[0]= x_front;
    p_in[1]= y_front;

    p_out[0]= x_back;
    p_out[1]= y_back;

    x_in[0] =x_front;
    x_in[1] =y_front;
    x_in[2] =z_front;
    x_out[0] =x_front;
    x_out[1] =y_front;
    x_out[2] =z_front;

    double dist=distmin_LineRectang(p_in,p_out,x1,x2,x3,x4);
//     cout << " Near cell at dist= " << dist << endl;
    x_in[0] +=xc;
    x_in[1] +=yc;
    x_in[2] +=zc;
    x_out[0] +=xc;
    x_out[1] +=yc;
    x_out[2] +=zc;
    return dist;
  }
  else if( indx >= 2 && indx <= 3)
  {
    if(debug) cout << " Hit the corner?? " << endl;
    int indx_first=0;
    int indx_last=-1;
    for (int it=1; it< indx+1; it++)
    {
      double dx = xhit[0][it]-xhit[0][indx_first];
      double dy = xhit[1][it]-xhit[1][indx_first];
      double dz = xhit[2][it]-xhit[2][indx_first];
      double dist = sqrt( dx*dx +dy*dy + dz*dz);

      if(dist > precision)
      {
        indx_last=it;
        break;
      }
    }
    if(indx_last == -1)
    {
      cout << " Internal error in OneParticleResponse::TraceInCell Hit the corner but indx_last= " << indx_last << endl;
      exit(1);
    }

    if(xhit[2][indx_first] > xhit[2][indx_last])
    {
       indx_first=indx_last;
       indx_last=0;
    }
    x_in[0]=xhit[0][indx_first];
    x_in[1]=xhit[1][indx_first];
    x_in[2]=xhit[2][indx_first];
    x_out[0]=xhit[0][indx_last];
    x_out[1]=xhit[1][indx_last];
    x_out[2]=xhit[2][indx_last];

    double dist = sx/2-fabs(x_in[0]+x_out[0])/2;
    double disty = sy/2-fabs(x_in[1]+x_out[1])/2;
    if(dist > disty ) dist=disty;
    x_in[0] +=xc;
    x_in[1] +=yc;
    x_in[2] +=zc;
    x_out[0] +=xc;
    x_out[1] +=yc;
    x_out[2] +=zc;
    return -dist;
  }
  else if( indx > 3) cout << " Internal error in OneParticleResponse::TraceInCell Finaly indx= " << indx << endl;


  return double(indx);
}

///////////////////////////////////////////////////////////////////////////////

void OneParticleResponse::SortTrace(void)
{
  if(cells_trace.size() <= 1) return;
  int npoints = cells_trace.size();
  vector < int > sorted( npoints, -1 );
  vector < int > index ( npoints );

  for( int it1=0; it1 != npoints; it1++ )
  {
    double zmin = 100000000000.;
    int indx_min=-1;
    for( int it2=0; it2 != npoints; it2++ )
    {
      if(sorted[it2] >= 0) continue;
      if(v_input_trace_z[it2] >= zmin) continue;
      zmin = v_input_trace_z[it2];
      indx_min=it2;
    }
    sorted[indx_min]=1;
    index[it1]=indx_min;
  }

// Create Sorted Trace object
  vector<size_t> sorted_cells_trace;
  vector<double> sorted_v_input_trace[3];
  vector<double> sorted_v_output_trace[3];

  for( int it=0; it != npoints; it++ )
  {
    int indx = index[it];
    if( indx < 0 || indx >= npoints ) cerr << " Internal error in OneParticleResponse::SortTrace Wrong indx=" << indx << endl;
    sorted_cells_trace.push_back(cells_trace[indx]);
    sorted_v_input_trace[0].push_back(v_input_trace_x[indx]);
    sorted_v_input_trace[1].push_back(v_input_trace_y[indx]);
    sorted_v_input_trace[2].push_back(v_input_trace_z[indx]);
    sorted_v_output_trace[0].push_back(v_output_trace_x[indx]);
    sorted_v_output_trace[1].push_back(v_output_trace_y[indx]);
    sorted_v_output_trace[2].push_back(v_output_trace_z[indx]);
  }
  cells_trace=sorted_cells_trace;
  v_input_trace_x=sorted_v_input_trace[0];
  v_input_trace_y=sorted_v_input_trace[1];
  v_input_trace_z=sorted_v_input_trace[2];
  v_output_trace_x=sorted_v_output_trace[0];
  v_output_trace_y=sorted_v_output_trace[1];
  v_output_trace_z=sorted_v_output_trace[2];

//   PrintTrace();
}

///////////////////////////////////////////////////////////////////////////////

void OneParticleResponse::PrintTrace(void)
{
// So first let's take a look on it.
  cout << " OneParticleResponse::PrintTrace in Calorimeter " << calorimeter->GetName() <<  " Trace size " << cells_trace.size() <<endl;
  for( size_t it=0; it != cells_trace.size(); it++ )
  {
    double sx = calorimeter->GetCells()[cells_trace[it]].GetCellType().GetSizeX();
    double sy = calorimeter->GetCells()[cells_trace[it]].GetCellType().GetSizeY();
    if ( sx < sy )
      cout << " Cell X " << cells_trace[it] << " Z in " << v_input_trace_z[it] << " Z out " << v_output_trace_z[it] <<
                                                             " dZ=" << v_output_trace_z[it]-v_input_trace_z[it] << endl;
    else
      cout << " Cell Y " << cells_trace[it] << " Z in " << v_input_trace_z[it] << " Z out " << v_output_trace_z[it] <<
                                                             " dZ=" << v_output_trace_z[it]-v_input_trace_z[it] << endl;
  }

}

////////////////////////////////////////////////////////////////////////////////

std::ostream& operator<<( std::ostream& o, OneParticleResponse& h) {
    o << "OneParticleResponse:" << endl;

    o << "Hit cells (" << h.hited_cells.size() << "):";
    for (std::vector<size_t>::const_iterator it=h.hited_cells.begin(); it!=h.hited_cells.end(); it++)
        o << " " << *it;
    o << endl;

    if ( h.e_expected_.size() == h.de_expected_.size() ) {
        o << "Expected energy (" << h.e_expected_.size() << "):";
        for (size_t i=0; i<h.e_expected_.size(); i++)
            o << " " << h.e_expected_[i] << "(" << h.de_expected_[i] << ")";
        o << endl;
    } else
        o << "Size of expected energies does not correspond to the expected energy errors!" << endl;

    if ( h.e_real_.size() == h.de_real_.size() ) {
        o << "Real energy (" << h.e_real_.size() << "):";
        for (size_t i=0; i<h.e_real_.size(); i++)
            o << " " << h.e_real_[i] << "(" << h.de_real_[i] << ")";
        o << endl;
    } else
        o << "Size of real energies does not correspond to the real energy errors!" << endl;

    o << endl;

    return o;
}

} // namespace Reco
