#include "Reco_config.h"
#include "Calorimeter.h"

#include "GUICalorimeter.h"

#if USE_Qt
#include <QtGui/QApplication>
#endif

#include "TApplication.h"
#include "TCanvas.h"
#include "TEnv.h"
#include "TSystem.h"

using namespace Reco;

class GuiCalorimeter : public Reco::Calorimeter {
    public:
        GuiCalorimeter(const std::string& name , const std::string& geom_file) : 
            Reco::Calorimeter(name, geom_file) {
            InitMemory();
        }

        virtual void Initialize() {
            options.mixed_blocks = true;
            options.SetOK();

            Reco::Calorimeter::Initialize();
        }
};

int main(int argc, char* argv[])
{
  //XInitThreads();

  try
  {
#if USE_Qt
    QApplication app_Qt(argc, argv);
#endif

    gEnv->SetValue("X11.XInitThread", 0);
    TApplication app_ROOT("gui", &argc, argv);

    GuiCalorimeter GAMS("GAMS","GAMS-4pi.geom");
    GAMS.Initialize();

    GAMS.CreateGUI();
#if USE_Qt
    GAMS.GetGUI().show();
#endif

    new TCanvas("c","cccccccccccccccc",444,333); // Open ROOT canvas

    while(1)
    {
#if USE_Qt
      app_Qt.processEvents();
#endif

      gSystem->ProcessEvents();
    }
  }
  catch( std::exception &e )
  {
    std::cerr << "Exception: " << e.what() << "\n";
  }
  catch( const char *s )
  {
    std::cerr << "Exception: " << s << "\n";
  }
  catch( ... )
  {
    std::cerr << "Unknown exception\n";
  }

  return 0;
}
