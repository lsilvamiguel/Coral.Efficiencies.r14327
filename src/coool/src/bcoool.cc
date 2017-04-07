
#include <math.h>

#include "Monitor.h"
#include "VariousSettings.h"

#include <map>
#include <string>

#include <TROOT.h>
#include <TApplication.h>
#include <TSystem.h>

TROOT root("GUI_BATCH", "GUI_BATCH  environement");

bool thr_flag = false;

namespace { // local routine
void usage(const char* pgm_name, const char* err_string) {

  if (err_string) std::cerr <<std::endl << err_string <<std::endl <<std::endl;
  std::cerr << "usage: " << gApplication->Argv(0)
       << " [-map map_file]"
       << " [-root root_file] [-rootdir root_dir]"
       << " [-text text_file] [-textdir text_dir]"
       << " [-group groups_file] [-nevent nb_of_events] [-evtspacing spacing of events]"
       << " [-cfg user_configuration_file]"
       << " [-trigmask (trigger mask binary value)]"
       << " [-evttypes (event types to read, given with that format:"
       << " -evttypes START_OF_RUN,START_OF_BURST,PHYSICS_EVENT,CALIBRATION_EVENT,...)]"
       << " [-notree] [-nocluster] [-nocalib] [-withref] [-withnofailedevt] [-experthistos]"
       << " [-cooolrc cooolrc file] [-CompMon cooolrc file (old name)]"
       << " [-workenv working environment string] [-geom detectors_geometry_file]"
       << " [-file_list file giving list of files (incompatible with data file argument)]"
       << " [data file or DATE source] [second data file] ...\n";
  exit(2);
}
}

int main(int argc, char **argv) {
  char *root_argv[]={argv[0],"-b",NULL};
  int root_argc=2;
  TApplication theApp("App", &root_argc, root_argv);

  gROOT->Macro("rootlogon.C");

  //for(int i=0; i<nfiles;i++) {

    std::string mapfile="/afs/cern.ch/compass/detector/maps";
    std::string rootfile="DEFAULT";
    std::string rootdir=".";
    std::string textfile="";
    std::string cooolrc ="";
    std::string textdir=".";
    std::string psfile;
    std::string datafile="/afs/cern.ch/compass/scratch/d72/neyret/cdr18009-22021.raw";
    std::list<std::string> datafilelist;
    std::string filelist="";
    std::string groupfile="/afs/cern.ch/compass/detector/monitor/groups.xml";
    std::string paramfile = "/afs/cern.ch/compass/detector/monitor/default_params";
//    std::string geomfile = "/afs/cern.ch/compass/detector/geometry/2002/detectors.22018.minus.dat";
    std::string geomfile = "/afs/cern.ch/compass/detector/geometry/2004/detectors.for.coool.dat";
    std::string configfile = "";
    std::string evttypes = "";
    const char* wkenv = "AFS";
    int nbofevt = 0;
    int evtspacing = 0;
    bool dotree = true;
    bool docluster = true;
    bool docalib = true;
    bool doshowref = false;
    bool dofailedevt = true;
    bool stopatendrun = false; 
    bool experthistos = false; 
    unsigned int trigmask = 0xffffffff;

    register int pargc = 1;


    // arguments parsing
    while (pargc < argc) {
      char *pargv, *pargvn;
      pargv =  argv[pargc];
      pargvn = 0;
      if ((pargc + 1) < argc) pargvn =  argv[pargc+1];

      if (pargv[0] == '-') {   // it is an option
        if (strcmp(pargv, "-map") == 0) {  // option to give Map file
          if (pargvn && (pargvn[0] != '-')) mapfile = pargvn;
            else usage(argv[0], "badly formed -map option");
          pargc+=2;
          continue;
        }
        if (strcmp(pargv, "-ps") == 0) {  // option to give ps output file
          if (pargvn && (pargvn[0] != '-')) psfile = pargvn;
            else usage(argv[0], "badly formed -ps option");
          pargc+=2;
          continue;
        }
        if (strcmp(pargv, "-root") == 0) {  // option to give root file
	  if (pargvn && (pargvn[0] != '-')) {
	    rootfile = pargvn;
	    pargc+=2;
	  }	    	
	  else usage(argv[0], "badly formed -root option");
	  continue;
        }
        if (strcmp(pargv, "-rootdir") == 0) {  // option to give output directory for root files
	  if (pargvn && (pargvn[0] != '-')) {
	    rootdir = pargvn;
	    pargc+=2;
	  }	    	
	  else usage(argv[0], "badly formed -rootdir option");
	  continue;
        }
        if (strcmp(pargv, "-text") == 0) {  // option to give output file name for text infos
	  if (pargvn && (pargvn[0] != '-')) {
	    textfile = pargvn;
	    pargc+=2;
	  }	    	
	  else usage(argv[0], "badly formed -text option");
	  continue;
        }
        if (strcmp(pargv, "-textdir") == 0) {  // option to give output directory for text infos
	  if (pargvn && (pargvn[0] != '-')) {
	    textdir = pargvn;
	    pargc+=2;
	  }	    	
	  else usage(argv[0], "badly formed -textdir option");
	  continue;
        }
        if (strcmp(pargv, "-group") == 0) {  // option to give groups file
          if (pargvn && (pargvn[0] != '-')) groupfile = pargvn;
	  else usage(argv[0], "badly formed -group option");
          pargc+=2;
          continue;
        }
        if (strcmp(pargv, "-cfg") == 0) {  // option to give user config file
          if (pargvn && (pargvn[0] != '-')) configfile = pargvn;
            else usage(argv[0], "badly formed -cfg option");
          pargc+=2;
          continue;
        }
        if (strcmp(pargv, "-geom") == 0) {  // option to give user config file
          if (pargvn && (pargvn[0] != '-')) geomfile = pargvn;
            else usage(argv[0], "badly formed -geom option");
          pargc+=2;
          continue;
        }
        if (strcmp(pargv, "-evttypes") == 0) {  // option to give event types to read
          if (pargvn && (pargvn[0] != '-')) evttypes = pargvn;
            else usage(argv[0], "badly formed -evttypes option");
          pargc+=2;
          continue;
        }
        if (strcmp(pargv, "-notree") == 0) {  // option to inhibate tree creation
          dotree = false;
          std::cout << "No tree created in the root file !\n";
          pargc++;
          continue;
        }
        if (strcmp(pargv, "-nocluster") == 0) {  // option to inhibate clustering
          docluster = false;
          dotree = false;  // because bcoool crashes if tree is filled without clustering (bug to find !)
          std::cout << "No clustering done ! Tree filling is removed too due to a coool bug.\n";
          pargc++;
          continue;
        }
        if (strcmp(pargv, "-nocalib") == 0) {  // option to inhibate clustering
          docalib = false;
          std::cout << "No calibration read !\n";
          pargc++;
          continue;
        }
        if (strcmp(pargv, "-withref") == 0) {  // option to enable reference plots
          doshowref = true;
          std::cout << "Reference plots will be added to the root file !\n";
          pargc++;
          continue;
        }
        if (strcmp(pargv, "-withfailedevt") == 0) {  // option to enable reading of evt badly decoded (for example, without trigger time)
          dofailedevt = true;
          std::cout << "Events partially decoded (without trigger time) will be read !\n";
          pargc++;
          continue;
        }
        if (strcmp(pargv, "-withnofailedevt") == 0) {  // option to disable reading of evt badly decoded (for example, without trigger time)
          dofailedevt = false;
          std::cout << "Events partially decoded (without trigger time) will be not read !\n";
          pargc++;
          continue;
        }
        if (strcmp(pargv, "-stopatendrun") == 0) {  // option to stop decoding if EOR event seen
          stopatendrun = true;
          std::cout << "bcoool will stop if an EOR event is seen !\n";
          pargc++;
          continue;
        }
        if (strcmp(pargv, "-experthistos") == 0) {  // option to activate expert histograms
          experthistos = true;
          std::cout << "bcoool will create expert histograms as well !\n";
          pargc++;
          continue;
        }
        if (strcmp(pargv, "-CompMon") == 0) {  // option to activate reading of .cooolrc file
          if (pargvn && (pargvn[0] != '-')) cooolrc = pargvn;
	  else usage(argv[0], "badly formed -CompMon option");
          pargc+=2;
          continue;
        }
        if (strcmp(pargv, "-cooolrc") == 0) {  // option to activate reading of .cooolrc file
          if (pargvn && (pargvn[0] != '-')) cooolrc = pargvn;
	  else usage(argv[0], "badly formed -cooolrc option");
          pargc+=2;
          continue;
        }
        if (strcmp(pargv, "-nevent") == 0) {  // option to give nb of events to read
          if (pargvn && (pargvn[0] != '-')) nbofevt = atoi(pargvn);
            else usage(argv[0], "badly formed -nb option");
          std::cout << nbofevt << " events to read\n";
          pargc+=2;
          continue;
        }
        if (strcmp(pargv, "-evtspacing") == 0) {  // option to give min event spacing
          if (pargvn && (pargvn[0] != '-')) evtspacing = atoi(pargvn);
            else usage(argv[0], "badly formed -evtspacing option");
          std::cout << "Event spacing value "<<nbofevt<<std::endl;
          pargc+=2;
          continue;
        }
        if (strcmp(pargv, "-trigmask") == 0) {  // option to give the trigger mask to use
          if (pargvn && (pargvn[0] != '-')) trigmask = strtoul(pargvn, 0, 0);
            else usage(argv[0], "badly formed -trigmask option");
          std::cout << nbofevt << " events to read\n";
          pargc+=2;
          continue;
        }
        if (strcmp(pargv, "-workenv") == 0) {  // option to give working environment string
          if (pargvn && (pargvn[0] != '-')) wkenv = pargvn;
            else usage(argv[0], "badly formed -workenv option");
          pargc+=2;
          continue;
        }
        if (strcmp(pargv, "-file_list") == 0) {  // option to give a file giving a list of data files
          if (pargvn && (pargvn[0] != '-')) filelist = pargvn;
            else usage(argv[0], "badly formed -file_list option");
          pargc+=2;
          continue;
        }
        if (strcmp(pargv, "--help") == 0) {  // option to request help
          usage(argv[0], 0);
        }
        // unknown option found
        char stmp[200];
        sprintf(stmp, "unknown option %s", pargv);
        usage(argv[0], stmp);
      }
      else {   // data file name
        // if ((pargc + 1) == argc) datafile = pargv;
        // else usage(argv[0],
        //           "data file name not at the end of command line or more than one data file name");
        datafilelist.push_back(pargv);
        if (filelist != "") {
          usage(argv[0], "data file argument and -file_list option incompatible, please remove one of those");
        }
        pargc++;
        continue;
      }
    }

    if (datafilelist.size() < 1) {    // no data file arg
      if (filelist == "") {
        datafilelist.push_back(datafile);  // if no data file arg and no file_list option

      } else {    // file_list given

        ifstream infl(filelist.c_str());
        if(! infl) {
          std::cerr<<"Error: can't open file_list file "<<filelist<<", exiting..."<<std::endl;
          exit(3);
        }
        while(1) {
          std::string fname;
          infl >> fname;
          if (infl.fail()) break;
          datafilelist.push_back(fname);
        }
        infl.close();
      }
    }

    Monitor *monitor=new Monitor();
    VariousSettings *settings = new VariousSettings(monitor);

    std::map<std::string,bool>& detintree=monitor->GetDetInTree();
    for (std::map<std::string,bool>::iterator dit=detintree.begin(); dit != detintree.end(); dit++) {
      dit->second=true;
    }
    if (cooolrc != "") settings->ReadPersoParam(cooolrc.c_str());
    monitor->DoTree(dotree);
    monitor->DoClustering(docluster);
    monitor->DoShowRef(doshowref);
    monitor->DoReadFailedEvt(dofailedevt);
    monitor->DoCalib(docalib);
    monitor->StopAtEndRun(stopatendrun);
    monitor->DoExpertHisto(experthistos);
    monitor->SetRootDir(rootdir.c_str());
    monitor->SetTextFile(textfile.c_str());
    monitor->SetTextDir(textdir.c_str());
    if (evttypes != "") monitor->CheckEvtTypes(evttypes.c_str());
    monitor->SetTriggerMask(trigmask);
    monitor->SetWkEnvStr(wkenv);
    monitor->Init(mapfile,rootfile,datafilelist,groupfile,geomfile);
    monitor->SetEventNumber(nbofevt);
    monitor->SetEventSpacing(evtspacing);

    if (paramfile != "") settings->LoadDefaultSettings(paramfile);

    if (configfile != "") {
     settings->LoadConfigSettings(configfile, 2, false);
    }

    monitor->Run();
    monitor->CloseRootFile();

    std::cout<<monitor->GetEventNumber()<<" events read" << std::endl;

    if(!psfile.empty()) {
      settings->GuiSettingsCreatePS(psfile);
      std::cout<<"ps file from "<<psfile<<" created"<<std::endl;
    }

    delete settings;
    delete monitor;
    gApplication->Terminate(0);
}






