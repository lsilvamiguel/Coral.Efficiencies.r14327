// $Id: main.cc,v 1.53 2010/09/07 18:27:02 tnagel Exp $

#include "coral_config.h"

#include <TApplication.h>
#include <TEnv.h>
#include <TSystem.h>

#include "Coral.h"
#include "CsCalorimeter.h"
#include "CsEvent.h"
#include "CsGeom.h"
#include "CsInit.h"
#include "CsRegistrySing.h"
#include "Reco/GUICalorimeter.h"

int main(int argc, char* argv[]) {
    try {
#if USE_Qt
        QApplication app_Qt(argc, argv);
#endif

        gEnv->SetValue("X11.XInitThread", 0);
        TApplication app_ROOT("gui", NULL, NULL);

        // Package Initialization
        Coral::init(argc, argv);
        //CsInit::Instance(argc, argv);
        CsEvent* event = CsEvent::Instance();

        int nevt=0;

        // Loop on events
        while( event->getNextEvent() ) {
            nevt++;
            if((nevt)%5 == 0 )
                std::cout << "Event: " << nevt << std::endl;
            if( !(CsRegistrySing::Instance()->callEoeMethods()) ) {
                break;
            }

            for (size_t i=0; i<CsGeom::Instance()->getCalorimeters().size(); i++) {
                CsCalorimeter* calo = CsGeom::Instance()->getCalorimeters()[i];
                calo->ReconstructionTest();
                calo->CreateGUI();
#if USE_Qt
                calo->GetGUI().show();
#endif
            }

            while (1) {
#if USE_Qt
                app_Qt.processEvents();
#endif
                gSystem->ProcessEvents();
            }
        }

        // End session
        CsRegistrySing::Instance()->callEndMethods();
        CsErrLog::Instance()->dump( elDebugging );
    } catch(std::exception &e ) {
        std::cerr << "Exception:\n" << e.what() << std::endl;
    } catch(std::string &e) {
        std::cerr << "Exception:\n" << e << "\n";
    } catch(const char *e) {
        std::cerr << "Exception:\n" << e << "\n";
    } catch( ... ) {
        std::cerr << "Unknown exception!!\n";
    }

    return 0;
}
