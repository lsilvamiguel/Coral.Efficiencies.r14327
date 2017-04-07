/*
   $Source: /afs/cern.ch/user/g/gobbo/w0/coral/src/calorim/Reco/test/t1.cc,v $ 
   $Date: 2000/07/25 15:34:09 $ 
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

#include <iostream>
#include <cstdlib>
#include <cmath>
#include <cassert>

#include "Reconstruction.h"
#include "Calorimeter.h"
#include "Event.h"

using namespace Reco;

/*
  Noise information from a single cell.
*/
class CellNoise
{
  public:
    CellNoise(double m=0,double s=1) : mean(m), sigma(s) {}
    double mean,sigma;
};

/*
  Accumulated information about cell noise.
*/
class CellNoiseInfo
{
  public:
                        CellNoiseInfo           (const string &name="",size_t n_bins=3,float min=0,float max=1);
    void                Add                     (const CellDataR &cell);
    string              name;
    float               min,max;
    vector<float>       hist;
    vector<float>       momenta;                // 0,1,2,3...
};

inline
CellNoiseInfo::CellNoiseInfo(const string &name,size_t n_bins,float min_,float max_) :
  min(min_), max(max_), momenta(3)
//  bins( vector<float>() )
{
  
}

inline
void CellNoiseInfo::Add(const CellDataR &cell)
{
  float mom=1;
  for(size_t i=0; i<momenta.size(); i++ )
  {
    momenta[i] += mom;
    mom *= cell.GetAmplitude();
  }
}

void calibration_analysis(const vector<EventDataR>    &events,
                          vector<CellNoise>           &cells_noise,
                          vector<CellNoiseInfo>       &cells_noise_info)
{
  cells_noise.clear();
  cells_noise_info.clear();
  if( events.empty() )
    return;

  const size_t N = events.front().N();
  
  cells_noise_info = vector<CellNoiseInfo>(N);

  for( vector<EventDataR>::const_iterator it=events.begin(); it!=events.end(); it++ )
  {
    assert( it->N()==N );  // all events must have the same cells amount.
    for( size_t i=0; i<N; i++ )
      cells_noise_info[i].Add( (*it)[i] );
  }
  
  for( size_t i=0; i<N; i++ )
  {
    assert( cells_noise_info[i].momenta.size()>=2 );
    float
      n     = cells_noise_info[i].momenta[0],
      mean  = cells_noise_info[i].momenta[1]/n,
      sigma = sqrt( cells_noise_info[i].momenta[2]/n - mean*mean );
    
    cells_noise.push_back( CellNoise(mean,sigma) );
  }
}

int main( int argc, char *argv[] )
{
  //
  // Parameters initialisation
  //

  if( argc>3 )
  {
    cerr << "Usage: rec [Nevents_phys] [Nevents_cal]\n";
    return 1;
  }

  const size_t
    N_events_phys = argc>=2 ? atoi(argv[1]) : 1,
    N_events_cal  = argc>=3 ? atoi(argv[2]) : 0;
  cout.form("It will be generated %d physical and %d calibrated events.\n",
            N_events_phys,N_events_cal);

  Calorimeter GAMS("GAMS");

  EventDataR GAMS_ExtrData;
  vector<EventDataR> calibration_events;
  vector<Particle> GAMS_Monte_particles;
  vector<Particle> GAMSparticles;
  vector<CellNoise> cells_noise_gen, cells_noise_guessed;
  const size_t N_cells = GAMS.GetCells().size();

  if( N_events_cal>0 )
  {
    /////////////////////////////////////////////////////////
    // Create cells noise.
    /////////////////////////////////////////////////////////
    for( vector<Cell>::iterator it=GAMS.GetCells().begin(); it!=GAMS.GetCells().end(); it++ )
    {
      double
        u1=2*M_PI*drand48(), u2=sqrt(-2*log(drand48())), r1=sin(u1)*u2, r2=cos(u1)*u2,
        // r1 and r2 are gaussian distributed independent random numbers.
        mean  = fabs(r1)+1,
        sigma = fabs(r2)+0.1;
      cells_noise_gen.push_back( CellNoise(mean,sigma) );
    }
    /////////////////////////////////////////////////////////

    /////////////////////////////////////////////////////////
    // Create events of calibration.
    /////////////////////////////////////////////////////////

    calibration_events.clear();
    for( size_t i=0; i<N_events_cal; i++ )
    {
      calibration_events.push_back(EventDataR());

      // Create noises in all cells.
      for( size_t j=0; j<N_cells; j++ )
      {
        double
          u1=2*M_PI*drand48(), u2=sqrt(-2*log(drand48())), r1=sin(u1)*u2, r2=cos(u1)*u2;
        // r1,r2,r3,r4 are gaussian distributed independent random numbers.
        calibration_events.back().Add( CellDataR(j,cells_noise_gen[j].mean+r1*cells_noise_gen[j].sigma) );
      }
    }
    /////////////////////////////////////////////////////////

    vector<CellNoiseInfo> cells_noise_info;
    calibration_analysis(calibration_events,cells_noise_guessed,cells_noise_info);

    if( cells_noise_guessed.size()!=N_cells )
    {
      cerr.form("Fatal error!  cells_noise_guessed.size()=%d  but N_cells=%d\n",
                 cells_noise_guessed.size(),N_cells);
      exit(1);
    }
    assert( cells_noise_guessed.size()==cells_noise_gen.size()    );
    assert( cells_noise_guessed.size()==cells_noise_info.size()   );

    cout << "Result of calibartion reconstruction:\n";
    for( size_t i=0; i<N_cells; i++ )
      cout.form("Cell %4d noise:  generated=(%+9.4f,%+9.4f)   reconstructed=(%+9.4f,%+9.4f)\n",
                  i,
                  cells_noise_gen    [i].mean , cells_noise_gen    [i].sigma,
                  cells_noise_guessed[i].mean , cells_noise_guessed[i].sigma );
  }

  //
  // loop on events
  //
  
  for( size_t nevt=0; nevt!=N_events_phys; nevt++ )
  {
    // Physical event.
    GAMS_ExtrData.Clear();
    GAMS_Monte_particles.clear();
    cout << " Event " << nevt << endl;

    GAMS_Monte_particles.push_back(Particle(1,1000.,3.5,3.5,0.,
                                            0,2.,2.,5.,0.02,-0.02,1));
    GAMS_Monte_particles.push_back(Particle(1,500.,45.,40.,0.,
                                            0,2.,2.,5.,0.02,-0.02,1));
    GAMS_Monte_particles.push_back(Particle(1,700.,12.,-30.,0.,
                                            0,2.,2.,5.,0.02,-0.02,1));
    GAMS_Monte_particles.push_back(Particle(1,1000.,-40.,10.,0.,
                                            0,2.,2.,5.,0.02,-0.02,1));

    GAMS.MonteConstruction(GAMS_Monte_particles);
    GAMS.ExtractMonteData(GAMS_ExtrData,10.);
    cout << " DATA Size " << GAMS_ExtrData.N()  <<  endl;
    GAMS.Clear();
    GAMS.InsertData(GAMS_ExtrData.GetData());
    GAMS.Reconstruction(GAMSparticles);
//      GAMS.Reconstruction(GAMSparticles,cells_noise_guessed);
  }
  return 0;
}
