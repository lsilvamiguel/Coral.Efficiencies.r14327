/*!
  \file    CsRichOneDisplay.cc
  \-----------------------
*/


#include <iostream>
#include <ostream>

#include <cstdio>
#include <cstdlib>

#include <unistd.h>
#include <fcntl.h>
#include <sys/types.h>
#include <sys/stat.h>

//---------------------------
#include <iostream>
#include <fstream>
#include <iomanip>
#include <sstream>
#include <string>
//----------------------------
#include "Coral.h"
#include "CsHistograms.h"
#include "CsRichOneDisplay.h"

#include "CsRCEDII.h"
#include "TROOT.h"
#include "TRint.h"
extern void InitGui();
TRint *theAppR;

using namespace std;

CsRichOneDisplay* CsRichOneDisplay::instance_ = 0;

CsRichOneDisplay* CsRichOneDisplay::Instance() {
  if( instance_ == 0 ) instance_ = new CsRichOneDisplay();
  return instance_;
}

CsRichOneDisplay::CsRichOneDisplay() {
  //  CsRegistry reg;
  //  reg.EOJRegistration( this );
  //  flag_ = 0;
}

CsRichOneDisplay::CsRichOneDisplay( const CsRichOneDisplay &richdisp ) {
  
  cout << "RICHONE : CsRichOne CopyConstructor" << endl;
  instance_ = instance_;
  
}

void CsRichOneDisplay::doRichOneDisplay() {

  static bool firstcall=true;

  if(firstcall){
    firstcall=false;
    VoidFuncPtr_t initfuncs[] = { InitGui, 0 };
    static TROOT ROOT("CsRoot","COMPASS",initfuncs);
    theAppR = new TRint("COMPASS RICH Event Display", NULL,NULL, 0);
    
    cout << " Init EventROOTDisplay RICH " << endl;
    
    graph_ = CsRCGraph::Init(initfuncs);
    IFdisplay_ = new CsRCIFtoDisplay();
    display_   = new CsRCDisplay(REmbedded);
    
  //@@--------------------------------------------
  }

  //----- EVENT DISPLAY ---------
  //-----------------------------
  //----- do the EventDisplay :
  //      ---------------------
  cout << " Next EventROOTDisplay RICH " << endl;
  
  //        ------------------------
  // Event Display stuff 
  
  CsRichOneDisplay* init = CsRichOneDisplay::Instance();

   if( (init->getDisplay())->GetContinueDisplay() ) {
    (init->getIFtoDisplay())->Trans( this );
    (init->getDisplay())->OpenEvent();
    if(!(init->getDisplay())->GetContinueDisplay() ) {
      (init->getDisplay())->CloseWindow();
    }
  }
  
}
