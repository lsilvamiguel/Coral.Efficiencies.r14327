#include <hash_map>

#include "Reco_config.h"
#include "CalorimeterMatrixSet.h"
#include "DataBase.h"

#if USE_Qt
#include <qapplication.h>
#include <qwidget.h>
#endif

#include "TROOT.h"
#include "TApplication.h"
#include "TSystem.h"
#include "TCanvas.h"
TROOT ROOT("","");

#include "GUICalorimeter.h"

using namespace Reco;

int main( int argc, char ** argv )
{
  try
  {
    #if USE_Qt
      QApplication app_Qt(argc,argv );
    #endif

      TApplication app_ROOT("gui",&argc,argv);

    Reco::CalorimeterMatrixSet GAMS("GAMS","GAMS-4pi.geom");
    GAMS.Init();

    hash_map<size_t,size_t> chan_to_cell;
    GAMS.ReadMap("GAMS-4pi.mapping",chan_to_cell);

    DataBase db("DB",DataBase::WRITE|DataBase::READ);
    if( db.WasCreated() )
      for( size_t i=0; i<GAMS.NCells(); i++ )
        GAMS.SetCellInfo(Calorimeter::CALIB,Calorimeter::OLD,i,CellInfo(2+i,1.+i/10.,1./(1+i)));
    else
      GAMS.GetInfoFromDataBase(db);

    GAMS.Print(cout,"GAMS: ");
    GAMS.CellsInfoPrint(Calorimeter::CALIB,Calorimeter::OLD);

    #if USE_Qt
      GAMS.CreateGUI(1);
      GAMS.GetGUI().show();
    #endif

      // new TCanvas("c","cccccccccccccccc",444,333); // Open ROOT canvas

    const double nominal_energy = 10;
    GAMS.SetCalibProc(true,nominal_energy);

    for( int i=1; i<=10000; i++ )
    {
      #if USE_Qt
        app_Qt.processEvents();
      #endif

        gSystem->ProcessEvents();

      vector<Calorimeter::Particle> particles_gen, particles_rec;
      double
        eg = 0.05*(2.*drand48()-1.)+nominal_energy,
        xg = GAMS.GetXmin()+(GAMS.GetXmax()-GAMS.GetXmin())*drand48(),
        yg = GAMS.GetYmin()+(GAMS.GetYmax()-GAMS.GetYmin())*drand48();
      if( eg<=0 )
        continue;
      particles_gen.push_back(Calorimeter::Particle(Calorimeter::Particle::GAMMA,1,
                                                eg,xg,yg,0.,    // E,X,Y,Z
                                                0,2.,2.,5.,     // Errors: E,X,Y,Z
                                                0.02,-0.02));   // Angles: X,Y

      vector<CellDataRaw> data;
      GAMS.MonteConstruction(particles_gen,data);
      GAMS.InsertData(data);
      GAMS.Reconstruction(particles_rec);

      if( i%1000==0 )
        GAMS.CellsInfoPrint(Calorimeter::CALIB,Calorimeter::NEW);
    }

    GAMS.PutInfoToDataBase(db);

    return 0;
  }
  catch( std::exception &e )
  {
    cerr << "Exception:\n" << e.what() << "\n";
  }
  catch( ... )
  {
    cerr << "Unknown exception:\n";
  }
  
  return 1;
}
