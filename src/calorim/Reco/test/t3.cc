/*
   $Source: /afs/cern.ch/user/g/gobbo/w0/coral/src/calorim/Reco/test/t3.cc,v $ 
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
#include <ctime>
#include <cassert>

#include "Reconstruction.h"
#include "Calorimeter.h"
#include "Event.h"

#include "TROOT.h"
#include "TFile.h"
TROOT ROOT("","");

using namespace Reco;

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

  Calorimeter GAMS("GAMS");

  TFile file("f.root");

  Run                   run(1,Run::READ);

  Event                *event           = NULL;
  EventPhysics         *event_phys      = new EventPhysics;
  EventCalibration     *event_cal       = new EventCalibration;
  EventNoises          *event_noise     = new EventNoises;
  EventPedestals       *event_ped       = new EventPedestals;
  
//   EventID event_id(0,run.ID());
//   for( size_t i=0; i<Nevents_total; i++ )
//   {
//     double r = drand48();
//     if( r<Nevents_phys/(double)Nevents_total )
//     {
//       event = event_phys;
//     }
//     else if( r<(Nevents_phys+Nevents_cal)/(double)Nevents_total )
//     {
//       event = event_cal;
//     }
//     else if( r<(Nevents_phys+Nevents_cal+Nevents_noise)/(double)Nevents_total )
//     {
//       event = event_noise;
//     }
//     else
//     {
//       event = event_ped;
//     }
//     
//     *event = event_id++;
//     cout << event->ID() << endl;
//     run.SaveEvent(event);
//   }
// 
  for( size_t i=0; run.GetEvent((Event*)event_phys,i); i++ )
  {
    cout << event_phys << endl;
    event_phys->Print();
  }
  
// 
// 
// 
// 
// 
//   EventDataR GAMS_ExtrData;
//   vector<EventDataR> calibration_events;
//   vector<Particle> GAMS_Monte_particles;
//   vector<Particle> GAMSparticles;
//   vector<CellNoise> cells_noise_gen, cells_noise_guessed;
//   const size_t N_cells = GAMS.GetCells().size();
// 
//   if( N_events_cal>0 )
//   {
//     /////////////////////////////////////////////////////////
//     // Create cells noise.
//     /////////////////////////////////////////////////////////
//     for( vector<Cell>::iterator it=GAMS.GetCells().begin(); it!=GAMS.GetCells().end(); it++ )
//     {
//       double
//         u1=2*M_PI*drand48(), u2=sqrt(-2*log(drand48())), r1=sin(u1)*u2, r2=cos(u1)*u2,
//         // r1 and r2 are gaussian distributed independent random numbers.
//         mean  = fabs(r1)+1,
//         sigma = fabs(r2)+0.1;
//       cells_noise_gen.push_back( CellNoise(mean,sigma) );
//     }
//     /////////////////////////////////////////////////////////
// 
//     /////////////////////////////////////////////////////////
//     // Create events of calibration.
//     /////////////////////////////////////////////////////////
// 
//     calibration_events.clear();
//     for( size_t i=0; i<N_events_cal; i++ )
//     {
//       calibration_events.push_back(EventDataR());
// 
//       // Create noises in all cells.
//       for( size_t j=0; j<N_cells; j++ )
//       {
//         double
//           u1=2*M_PI*drand48(), u2=sqrt(-2*log(drand48())), r1=sin(u1)*u2, r2=cos(u1)*u2;
//         // r1,r2,r3,r4 are gaussian distributed independent random numbers.
//         calibration_events.back().Add( CellDataR(j,cells_noise_gen[j].mean+r1*cells_noise_gen[j].sigma) );
//       }
//     }
//     /////////////////////////////////////////////////////////
// 
//     vector<CellNoiseInfo> cells_noise_info;
//     calibration_analysis(calibration_events,cells_noise_guessed,cells_noise_info);
// 
//     if( cells_noise_guessed.size()!=N_cells )
//     {
//       cerr.form("Fatal error!  cells_noise_guessed.size()=%d  but N_cells=%d\n",
//                  cells_noise_guessed.size(),N_cells);
//       exit(1);
//     }
//     assert( cells_noise_guessed.size()==cells_noise_gen.size()    );
//     assert( cells_noise_guessed.size()==cells_noise_info.size()   );
// 
//     cout << "Result of calibartion reconstruction:\n";
//     for( size_t i=0; i<N_cells; i++ )
//       cout.form("Cell %4d noise:  generated=(%+9.4f,%+9.4f)   reconstructed=(%+9.4f,%+9.4f)\n",
//                   i,
//                   cells_noise_gen    [i].mean , cells_noise_gen    [i].sigma,
//                   cells_noise_guessed[i].mean , cells_noise_guessed[i].sigma );
//   }
// 
//   //
//   // loop on events
//   //
//   
//   for( size_t nevt=0; nevt!=N_events_phys; nevt++ )
//   {
//     // Physical event.
//     GAMS_ExtrData.Clear();
//     GAMS_Monte_particles.clear();
//     cout << " Event " << nevt << endl;
// 
//     GAMS_Monte_particles.push_back(Particle(1,1000.,3.5,3.5,0.,
//                                             0,2.,2.,5.,0.02,-0.02,1));
//     GAMS_Monte_particles.push_back(Particle(1,500.,45.,40.,0.,
//                                             0,2.,2.,5.,0.02,-0.02,1));
//     GAMS_Monte_particles.push_back(Particle(1,700.,12.,-30.,0.,
//                                             0,2.,2.,5.,0.02,-0.02,1));
//     GAMS_Monte_particles.push_back(Particle(1,1000.,-40.,10.,0.,
//                                             0,2.,2.,5.,0.02,-0.02,1));
// 
//     GAMS.MonteConstruction(GAMS_Monte_particles);
//     GAMS.ExtractMonteData(GAMS_ExtrData,10.);
//     cout << " DATA Size " << GAMS_ExtrData.N()  <<  endl;
//     GAMS.Clear();
//     GAMS.InsertData(GAMS_ExtrData.GetData());
//     GAMS.Reconstruction(GAMSparticles);
// //      GAMS.Reconstruction(GAMSparticles,cells_noise_guessed);
//   }
  file.Close();
  return 0;
}
