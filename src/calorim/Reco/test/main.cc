/*
   $Source: /afs/cern.ch/user/g/gobbo/w0/coral/src/calorim/Reco/test/main.cc,v $ 
   $Date: 2011/03/01 14:14:00 $ 
   $Revision: 1.8 $ 
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

#include <unistd.h>
#include <fcntl.h>

#include <iostream>
#include <cstdlib>
#include <algo.h>
#include <cmath>
#include <cstdio>

#include "TROOT.h"
#include "TSystem.h"
#include "TApplication.h"
#include "TCanvas.h"
#include "TFile.h"
#include "TH1.h"
#include "TH2.h"
TROOT ROOT("",""); // Global ROOT initialization.

#include "Calorimeter.h"
#include "DataBase.h"
#include "Calorimeter.h"
#include "CalorimeterGAMS.h"

using namespace Reco;

////////////////////////////////////////////////////////////////////////////////

int main( int argc, char *argv[] )
{
  try
  {
    { Calorimeter m("LG1A","calorimeters.dat"); }

    TApplication a("rec",&argc,argv);

    if( argc>2 )
    {
      cerr << "Usage: rec [Nevents_phys]\n";
      return 1;
    }

    const size_t N_events = argc==1 ? 1 : atoi(argv[1]);
    cout << "It will be generated " << N_events << " events.\n";

    TFile file("file.root","RECREATE","",9);

  // Calibration test 
    double coef = 0.5;
    CellInfo calib;  
    for( int i = 0; i < 10000 ; i++)
    {
      double
        u1=2*M_PI*drand48(), u2=sqrt(-2*log(drand48())), r1=sin(u1)*u2, r2=cos(u1)*u2,
        ss = 0.2*fabs(r2)+0.1,
        mm  = coef*ss*r1+2,
        ww  = 1./(coef*ss);
        ww  *= ww;
        calib.Add(mm/2.,ww);
    }

    double n,m,s;
    calib.Result(n,m,s);
    cout << " Calibration Stat " << n << " Mean " << m << " Sigma " << s << endl;

    CalorimeterGAMS GAMS("GAMS");
    
    DataBase db("DB",DataBase::WRITE|DataBase::READ|DataBase::CREATE);
    GAMS.PutInfoToDataBase(db);

    vector<CellDataRaw> data;
    vector<Calorimeter::Particle> GAMS_Monte_particles;
    vector<Calorimeter::Particle> GAMSparticles;
    vector<Calorimeter::Particle> GAMS_Answer_particles;

    TCanvas *canvas = new TCanvas(GAMS.GetName().c_str(),"Calorimeter display",500,500);

    //
    // loop calibration iterations
    //

    size_t max_iter = 5;

    for( size_t niter=0; niter <= max_iter-1; niter++ ) 
    {

    //
    // loop on events
    //

      for( size_t nevt=0; nevt!=N_events; nevt++ ) 
      {

        if( nevt%100==0 )
          cout << "Iter " << niter << "   Event " << nevt << endl;

        GAMS_Monte_particles.clear();
        //double eg  = 2.*(2.*drand48()-1.)+3., xg = 40.*(2.*drand48()-1.),
        //yg = 40.*(2.*drand48()-1.);
        // Electron calibtartion 10 GeV
        double eg  = 0.05*(2.*drand48()-1.)+10., xg = 20.*(2.*drand48()-1.),
        yg = 20.*(2.*drand48()-1.);
        // eg xg yg are uniformly distributed independent random numbers.

        GAMS_Monte_particles.push_back(Calorimeter::Particle(Calorimeter::Particle::GAMMA,1,
                                                  eg,xg,yg,0.,    // E,X,Y,Z
                                                  0,2.,2.,5.,     // Errors: E,X,Y,Z
                                                  0.02,-0.02));   // Angles: X,Y

        GAMS.MonteConstruction(GAMS_Monte_particles,data);


        // TODO: Print out function for Calorimeter Data
        //cout << " DATA Size " << data.size()  <<  endl;

        GAMS.InsertData(data);
        GAMS.Reconstruction(GAMSparticles);


  //  Some reconstruction test(Reconstruction case 3)

   /*
      GAMS_Answer_particles.clear();

      //GAMS_Answer_particles.push_back(Particle(1,1000.,3.5,3.5,0.,
      //                                        0,2.,2.,5.,0.02,-0.02,1));
      //GAMS_Answer_particles.push_back(Particle(1,500.,45.,40.,0.,
      //                                        0,2.,2.,5.,0.02,-0.02,1));
      //GAMS_Answer_particles.push_back(Particle(1,700.,12.,-30.,0.,
      //                                        0,2.,2.,5.,0.02,-0.02,1));
      GAMS_Answer_particles.push_back(Particle(1,990.,-40.,10.,0.,
                                              0,2.,2.,5.,0.02,-0.02,1));


      GAMS.Reconstruction(GAMS_Answer_particles);

   */
  //  Test of reconstruction in Calorimeter : Compare generated vector of
  //   particles with reconstructed one.
  //  
        if(niter == max_iter-1)
          GAMS.ReconstructionTest(GAMS_Monte_particles,GAMSparticles);

        GAMS.Draw(canvas,3);

        int ret = fcntl( fileno(stdin), F_SETFL, O_NONBLOCK);
        assert( ret != -1 );
        char c;
        do
          gSystem->ProcessEvents();
        while( -1==scanf("%c",&c) );
        ret = fcntl( fileno(stdin), F_SETFL, ~O_NONBLOCK);
        assert( ret != -1 );
      }

    // End of iteration.
//      GAMS.PutInfoToDataBase();

    }

    // End of session.
    file.Write();
    file.Close();

    cout << " That's all C++ lovers ..." << endl;
  }
  catch( std::exception &e )
  {
    cerr << "std::exception:\n" << e.what() << endl;
  }
  catch( ... )
  {
    cerr  << __FILE__ << " in " << __FUNCTION__ << " at " << __LINE__ << " Unknown exception\n";
  }

  return 0;
}

////////////////////////////////////////////////////////////////////////////////
