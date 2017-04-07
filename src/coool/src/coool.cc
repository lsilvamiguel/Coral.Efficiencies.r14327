#include <TROOT.h>
#include <TStyle.h>
#include <TApplication.h>
#include <TVirtualX.h>

#include "MainFrame.h"


TROOT root("GUI", "GUI test environement");

bool thr_flag = false;
//bool thr_flag = true;


int main(int argc, char **argv)
{

  gROOT->SetMacroPath(gSystem->Getenv("HOME"));
  gROOT->Macro("cooollogon.C");

  MainFrame *mainframe=0;
  try{
    char *root_argv[]={argv[0],NULL};
    int root_argc=1;
    TApplication theApp("App",&root_argc, root_argv);

    if (gROOT->IsBatch()) {
      fprintf(stderr, "%s: cannot run in batch mode\n", argv[0]);
      return 1;
    }

    
    mainframe=new MainFrame(gClient->GetRoot(),argc,argv, 400, 220);   

    theApp.Run();
  }   
  // Something is wrong...
  // Print error message and finish this data file.
      
  catch( const std::exception &e ) {
    std::cerr << "exception:\n" << e.what() << "\n";
  }
  catch( const char * s ) {
    std::cerr << "exception:\n" << s << "\n";
  }
  catch( ... ) {
    std::cerr << "Oops, unknown exception!\n";
  }
  
#ifndef __CINT__
  CS::Chip::PrintStatistics();
  CS::Exception::PrintStatistics();
#endif
  
  delete mainframe;
  exit(0);
}



