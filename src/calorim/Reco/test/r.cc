/*
   $Source: /afs/cern.ch/user/g/gobbo/w0/coral/src/calorim/Reco/test/r.cc,v $ 
   $Date: 2000/07/25 15:33:11 $ 
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
#include <ctime>
#include <cassert>

#include "Calorimeter.h"
#include "CalorimeterGAMS.h"

#include "TROOT.h"
#include "TFile.h"
TROOT ROOT("","");

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

int main( int argc, char *argv[] )
{
  ;{;
    // Random number initialisation.
    time_t t;
    time(&t);    // Get current time.
    srand48(t);  // Set seed for drand48() random number generator.
    // cout << "Random init test: " << drand48() << "\n";
    //HepRandom::setTheSeed((long)(100000000*drand48()));
  ;};

  if( argc>5 )
  {
    cerr << "Usage: gen [Nevents_phys] [Nevents_noise] [Nevents_ped] [Nevents_cal]\n";
    return 1;
  }

  const size_t
    Nevents_phys  = argc>=2 ? atoi(argv[1]) : 0,
    Nevents_noise = argc>=3 ? atoi(argv[2]) : 0,
    Nevents_ped   = argc>=4 ? atoi(argv[3]) : 0,
    Nevents_cal   = argc>=5 ? atoi(argv[4]) : 0,
    Nevents_total = Nevents_phys+Nevents_noise+Nevents_ped+Nevents_cal;
  cout.form("Nevents_phys=%d  Nevents_noise=%d  Nevents_ped=%d  Nevents_cal=%d\n",
             Nevents_phys, Nevents_noise, Nevents_ped, Nevents_cal);

  if( Nevents_total==0 )
    return 0;

  CalorimeterGAMS GAMS("GAMS",2000);

  EventDataR GAMS_ExtrData;
  vector<Particle> GAMS_Monte_particles;
  vector<Particle> GAMSparticles;
  vector<CellNoise> cells_noise_gen;
  const size_t N_cells = GAMS.GetCells().size();
  cout << "N cells = " << N_cells << "\n";

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

  TFile file("f.root","RECREATE","",9);

  Run                   run(1,Run::WRITE);

  Event                *event           = NULL;
  EventPhysics         *event_phys      = new EventPhysics;
  EventCalibration     *event_cal       = new EventCalibration;
  EventNoises          *event_noise     = new EventNoises;
  EventPedestals       *event_ped       = new EventPedestals;
  
  EventID event_id(0,run.ID());         // Event identification (event number + run number)

  for( size_t i=0; i<Nevents_total; i++ )       // Cycle for all events.
  {
    event_id++;
  
    double r = drand48();
    if( r<Nevents_phys/(double)Nevents_total )
    {
      // Physical event.

      event = event_phys;

      // Move to MonteConstruction GAMS_ExtrData.Clear();
      GAMS_Monte_particles.clear();

      GAMS_Monte_particles.push_back(Particle(Particle::GAMMA,1,1000.,3.5,3.5,0.,
                                              0,2.,2.,5.,0.02,-0.02));
      GAMS_Monte_particles.push_back(Particle(Particle::GAMMA,1,500.,45.,40.,0.,
                                              0,2.,2.,5.,0.02,-0.02));
      GAMS_Monte_particles.push_back(Particle(Particle::GAMMA,1,700.,12.,-30.,0.,
                                              0,2.,2.,5.,0.02,-0.02));
      GAMS_Monte_particles.push_back(Particle(Particle::GAMMA,1,1000.,-40.,10.,0.,
                                              0,2.,2.,5.,0.02,-0.02));

      GAMS.MonteConstruction(GAMS_Monte_particles,GAMS_ExtrData);
      // Move to MonteConstruction GAMS.ExtractMonteData(GAMS_ExtrData,10.);
      event->cells = GAMS_ExtrData.GetData();   // Simple copy. Slow, but it is all right.
    }
    else if( r<(Nevents_phys+Nevents_cal)/(double)Nevents_total )
    {
      // Event of calibration.
      event = event_cal;
      // put your code here.
    }
    else if( r<(Nevents_phys+Nevents_cal+Nevents_noise)/(double)Nevents_total )
    {
      // Events for cell noises determination.

      event = event_noise;
      
      event->cells.clear();

      // Create noises in all cells.
      for( size_t j=0; j<N_cells; j++ )
      {
        double
          u1=2*M_PI*drand48(), u2=sqrt(-2*log(drand48())), r1=sin(u1)*u2, r2=cos(u1)*u2;
        // r1,r2,r3,r4 are gaussian distributed independent random numbers.
        event->cells.push_back( CellDataR(j,cells_noise_gen[j].mean+r1*cells_noise_gen[j].sigma) );
      }
    }
    else
    {
      // Pedestal events
      event = event_ped;
      // put your code here.
    }

    *event = event_id;
    run.SaveEvent(event);
  }
  
  //////////////////////////////////////////////////////////////////
  // This is analysis of events with cells noises.
  //////////////////////////////////////////////////////////////////

  if( Nevents_noise>0 )
  {
    for( size_t i=0; run.GetEvent((Event*)event_noise,i); i++ )   // Cycle on all events in one run
      GAMS.InsertCellsNoiseData(event_noise->cells);

    cout << "Result of calibartion reconstruction:\n";
    for( size_t i=0; i<N_cells; i++ )
      cout.form("Cell %4d noise:  generated=(%+9.4f,%+9.4f)   reconstructed=(%+9.4f,%+9.4f)\n",
                  i,
                  cells_noise_gen    [i].mean , cells_noise_gen    [i].sigma,
                  GAMS.GetCellNoiseInfo(i).GetMean(), GAMS.GetCellNoiseInfo(i).GetSigma() );
  }


  GAMS.PrintCellsNoiseInfo();  

  file.Write();
  file.Close();
  return 0;
}
