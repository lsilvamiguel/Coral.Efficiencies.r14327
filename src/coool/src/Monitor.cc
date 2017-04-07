#include "Monitor.h"
#include "../expat/xmlparse/xmlparse.h"
#include "MonitorPanel.h"
#include <stdio.h>
#include <sys/resource.h>
#include <unistd.h>

#include "Exception.h"
#include "DaqEvent.h"

#include "Plane1V.h"
#include "PlaneVeto.h"
#include "PlaneTrigHodo.h"
#include "PlaneMISC.h"
#include "PlaneDriftChamber.h"
#include "PlaneDriftChamber2V.h"
#include "PlaneScifiJ.h"
#include "PlaneScifiG.h"
#include "PlaneMwpc.h"
#include "PlaneStrawTubes.h"
#include "PlaneECAL0.h"
#include "PlaneHCAL1.h"
#include "PlaneHCAL2.h"
#include "PlaneHCALT.h"
#include "PlaneHCAL2Sum.h"
#include "PlaneECAL1.h"
#include "PlaneFEM.h"
#include "PlaneECAL2.h"
#include "PlaneAPV.h"
#include "PlaneGEM.h"
#include "PlanePGEM.h"
#include "PlaneSili.h"
#include "PlaneRiAPV.h"
#include "PlaneMumega.h"
#include "PlaneMuonWallA.h"
#include "PlaneMuonWallB.h"
#include "PlanePMumega.h"
#include "PlaneBMS.h"
#include "PlaneScaler.h"
#include "PlaneRICH.h"
#include "PlaneRICH_MAPMT.h"
#include "PlaneCamera.h"
#include "PlaneGandalf.h"

#include "PlanePrimakoffHodo.h"
#include "PlaneBeamKiller.h"
#include "PlaneVBHodo.h"
#include "PlaneVBOX.h"
#include "PlaneCEDAR.h"
#include "PlaneRW.h"
#include "PlaneRPD_SADC.h"
#include "PlaneRPD_F1.h"
#include "PlaneSandwich.h"
#include "PlaneTrigger_SADC_F1.h"
#include "PlaneTcsPhase.h"

#include "GroupScifiJ.h"
#include "GroupScifiG.h"
#include "GroupTrigHodo.h"
#include "GroupDAQ.h"
#include "GroupMwpc.h"
#include "GroupMumega.h"
#include "GroupGEM.h"
#include "GroupPGEM.h"
#include "GroupSili.h"
#include "GroupMuonWallA.h"
#include "GroupStraws.h"
#include "GroupBeamKiller.h"
#include "GroupCEDAR.h"
#include "GroupRiAPV.h"
#include "GroupRICH_MAPMT.h"
#include "GroupRW.h"
#include "GroupRPD.h"
#include "GroupECAL1.h"
#include "GroupCamera.h"

#define THREADWAITTIME 50

#if USE_DATABASE == 1
#include "MySQLDB.h"
#define DBSERVER "wwwcompass.cern.ch"
#define DBUSER "anonymous"
#define DBNAME "runlb"
#endif

#include "Geometry.h"

#if __GNUC__ >= 3
__const unsigned short int *__ctype_b = *__ctype_b_loc ();
#endif


std::map<std::string,Plane*> Monitor::fDetMap;
std::vector<Plane*> Monitor::fPlanes;
std::map<std::string,Group*> Monitor::fGroupMap;
std::map<std::string,Plane*> Monitor::fDetMapRegistered;
std::map<std::string,Plane*> Monitor::fDetMapNotRegistered;

const char* Monitor::fAllPlanesName = "All Planes";

// workaround to the fact that this srcIDtoscan is deleted at each new run
std::set<CS::uint16>      srcID_to_scan_to_reload;


ClassImp(Monitor)



Monitor::Monitor() :
TObject(),fTh(0),fThreadRun(false), fClosingWindow(false), fRunClear(true), fTree(0),
  fDoTree(false), fDoTracking(false), fDoClustering(true),
  fDoTextOutput(false), fDoCalib(true), fStopAtEndRun(false),
  fReadFailedEvt(true), fExpertHistos(false), fStopGo(false),
  fRootFile(0), fRootFileFlushed(false), fNevent(0),
  fMinEvtspacing(0),
  flag_end(false), fThreadFinished(false), fIsInitialized(false),
  fControlPanel(0),
  fTMPattern(0xffffffff), fTMLowerOnly(true), fTMStrict(false), fWkEnvStr("AFS") {

  fEvtTypesOK.insert(CS::DaqEvent::PHYSICS_EVENT);
  fEvtTypesOK.insert(CS::DaqEvent::CALIBRATION_EVENT);

  fEvtTypes[CS::DaqEvent::START_OF_RUN]=("START_OF_RUN");
  fEvtTypes[CS::DaqEvent::END_OF_RUN]=("END_OF_RUN");
  fEvtTypes[CS::DaqEvent::START_OF_RUN_FILES]=("START_OF_RUN_FILES");
  fEvtTypes[CS::DaqEvent::END_OF_RUN_FILES]=("END_OF_RUN_FILES");
  fEvtTypes[CS::DaqEvent::START_OF_BURST]=("START_OF_BURST");
  fEvtTypes[CS::DaqEvent::END_OF_BURST]=("END_OF_BURST");
  fEvtTypes[CS::DaqEvent::PHYSICS_EVENT]=("PHYSICS_EVENT");
  fEvtTypes[CS::DaqEvent::CALIBRATION_EVENT]=("CALIBRATION_EVENT");
  fEvtTypes[CS::DaqEvent::END_OF_LINK]=("END_OF_LINK");
  fEvtTypes[CS::DaqEvent::EVENT_FORMAT_ERROR]=("EVENT_FORMAT_ERROR");

  fDetInTree["BMS"]=true;
  fDetInTree["CEDARs"]=true;
  fDetInTree["DriftChambers"]=true;
  fDetInTree["Camera"]=false;
  fDetInTree["ECAL1"]=true;
  fDetInTree["FEM"]=false;
  fDetInTree["ECAL2"]=true;
  fDetInTree["GEMs"]=true;
  fDetInTree["PixelGEMs"]=true;
  fDetInTree["HCAL1"]=true;
  fDetInTree["ECAL0"]=true;
  fDetInTree["HCAL2"]=true;
  fDetInTree["Hodoscopes"]=true;
  fDetInTree["Micromegas"]=false;
  fDetInTree["PixelMicromegas"]=true;
  fDetInTree["MuonWallA"]=true;
  fDetInTree["MuonWallB"]=true;
  fDetInTree["MWPC's"]=true;
  fDetInTree["RICH_APV"]=false;
  fDetInTree["RICH_BORA"]=false;
  fDetInTree["RICH_MAPMT"]=false;
  fDetInTree["RichWall"]=true;
  fDetInTree["RPD"]=false;
  fDetInTree["Sandwich"]=false;
  fDetInTree["Scalers"]=true;
  fDetInTree["SciFis(Germany)"]=true;
  fDetInTree["SciFis(Japan)"]=true;
  fDetInTree["SciFis(Warsaw)"]=false;
  fDetInTree["Silicons"]=false;
  fDetInTree["StrawTubes"]=true;
  fDetInTree["VetoBox"]=true;
  fDetInTree["W45"]=true;

#if USE_DATABASE == 1
  fDataBase = 0;
#endif

  fBenchDecoding=new TStopwatch();
  fBenchDecoding->Stop();
  fBenchClustering=new TStopwatch();
  fBenchClustering->Stop();
  fBenchTracking=new TStopwatch();
  fBenchTracking->Stop();
}


void Monitor::Init(std::string mapfile, std::string rootfile, std::string datafile, std::string groupfile, std::string geomfile) {

  std::list<std::string> datafilelist;
  datafilelist.push_back(datafile);
  Init(mapfile, rootfile, datafilelist, groupfile, geomfile);
}


void Monitor::Init(std::string mapfile, std::string rootfile, std::list<std::string>& datafilelist, std::string groupfile, std::string geomfile) {

  //CS::Exception::SetActivity(false);

  typedef std::list<std::string>::iterator LSiter;
  fNfiles=datafilelist.size();
  if (fNfiles < 1) {
    std::cerr<<"Monitor::Init: no datafile given, can't do anything !"<<std::endl;
    exit(1);
  }

  std::string& datafile1 = *datafilelist.begin();

  std::cout<<std::endl;
  std::cout<<std::endl;
  std::cout<<std::endl;
  std::cout<<"             CC     OO     OO     OO    L   "<<std::endl;
  std::cout<<"            C  C   O  O   O  O   O  O   L   "<<std::endl;
  std::cout<<"            C      O  O   O  O   O  O   L   "<<std::endl;
  std::cout<<"            C      O  O   O  O   O  O   L   "<<std::endl;
  std::cout<<"            C  C   O  O   O  O   O  O   L   "<<std::endl;
  std::cout<<"             CC     OO     OO     OO    LLLL"<<std::endl;
  std::cout<<std::endl;
  std::cout<<std::endl;
  std::cout<<"************* Starting monitoring session **************************"<<std::endl;
  std::cout<<"mapping file : "<<mapfile<<std::endl;
  std::cout<<"root file    : "<<rootfile<<std::endl;
  if (fNfiles == 1) {
    std::cout<<"data file    : "<<datafile1<<std::endl;
  } else {
    register int fcount = 1;
    for (LSiter ff = datafilelist.begin(); ff != datafilelist.end(); ff++) {
      std::cout<<"data file "<<fcount<<"  : "<<*ff<<std::endl;
      fcount++;
    }
  }
  std::cout<<"geometry     : "<<geomfile<<std::endl;
  std::cout<<"group file   : "<<groupfile<<std::endl;
  std::cout<<"working env  : "<<fWkEnvStr<<std::endl;
  std::cout<<"********************************************************************"<<std::endl;
  std::cout<<std::endl;

  // Just getting out the run number, and the list of detectors
#ifndef __CINT__
  fManager=new CS::DaqEventsManager;
  fManager->SetMapsDir(mapfile);
  fManager->AddDataSource(datafile1.c_str());

  std::cout<<std::endl;
  std::cout<<"Opening data file :"<<std::endl;
  std::cout<<"********************************************************************"<<std::endl;

  if( !fManager->ReadEvent() )
    {
      CS::Exception::PrintStatistics();
      exit(1);
    }

  CS::DaqEvent &event=fManager->GetEvent();
  fRunNumber=event.GetRunNumber();
  Monitor::fPlanes.resize(fManager->GetDetectors().size(),NULL);

  fRunTime = *localtime((time_t*)&((event.GetTime()).first));

  std::cout<<std::endl;
  std::cout<<"Run number "<<fRunNumber<< " start time "<<asctime(&fRunTime)<<std::endl;
  std::cout<<std::endl;

  // set working env in planes and groups
  Plane::setWkEnvStr(fWkEnvStr.c_str());
  Group::setWkEnvStr(fWkEnvStr.c_str());

  /// set expert histo flag
  Plane::setExpertHistoFlag(fExpertHistos);
  Group::setExpertHistoFlag(fExpertHistos);

  /// set do tracking flag
  Plane::setDoTrackingFlag(fDoTracking);
  Group::setDoTrackingFlag(fDoTracking);

  // creating all planes:
  FillDetList();

  // initialize database
#if USE_DATABASE == 1
  std::cout<<std::endl;
  std::cout<<std::endl;
  std::cout<<"Connect to MySQL database :"<<std::endl
	   <<"********************************************************************"<<std::endl;
  fDataBase = new MySQLDB(DBSERVER, DBUSER, "", DBNAME);
  fDataBase->ConnectDB();
  fDataBase->setSpecialPlace(fWkEnvStr);
  std::cout<<std::endl;
#endif

  // initialize geometry
  std::cout<<"Intialize geometry :"<<std::endl
	   <<"********************************************************************"<<std::endl;
  Geometry *geom=0;
  try {
    geom = new Geometry(geomfile.c_str());
  }
  catch(...) {
    std::cerr<<"WARNING Monitor::Init. Geometry problem, tracking deactivated"<<std::endl;
    geom=0;
  }
  std::cout<<std::endl;


  //global tree. will be branched by detectors in fDetectors
  if (rootfile == "DEFAULT") {

    // WARNING !!! will not work in case of a stream !!!
    rootfile = datafile1;
    rootfile.erase(0,rootfile.rfind("/")+1);
    if (rootfile.rfind(".") != (size_t)-1) rootfile.erase(rootfile.rfind("."),rootfile.size());
    std::string dir;
    if(!fRootDir.empty()) {
      dir=fRootDir + "/";
    }
    else {
      dir=getenv("HOME");
      dir += "/w0/";
    }
    rootfile = dir + rootfile;
    rootfile += ".root";
  }
  if (rootfile != "")
    fRootFile=new TFile(rootfile.c_str(),"recreate");
  printf("Root tree: %s\n",rootfile.c_str());
  if(fRootFile && fRootFile->IsOpen()) {
    if(fDoTree) {
      fTree=new TTree("T","COMPASS monitoring");
      fTree->Branch("runno",&fRunNumber,"runno/I",32000);
      fTree->Branch("spill",&fSpill,"spill/I",32000);
      fTree->Branch("evt",&fNevent,"evt/I",32000);
      fTree->Branch("evtinspl",&fEvtinspl,"evtinspl/I",32000);
      fTree->Branch("timesec",&fTimesec,"timesec/I",32000);
      fTree->Branch("timemic",&fTimemic,"timemic/I",32000);
      fTree->Branch("evtype",&fEvtType,"evtype/I",32000);
      fTree->Branch("pattern",&fPattern,"pattern/i",32000);
      fTree->Branch("attributes",&fAttributes,"attributes/I",32000);
      fTree->Branch("filterbits",&fFilterBits,"filterbits/I",32000);
    }
  }
  else {
    if (rootfile != "") {
      std::cout<<"WARNING : Monitor::Init : Cannot open root output file "
	       <<rootfile.c_str()<<std::endl;
    }
    fRootFile = 0;
  }

  // loop over all planes
  typedef std::map<std::string,Plane*>::iterator DI;
  for(DI i=Monitor::fDetMap.begin(); i!=Monitor::fDetMap.end(); i++ ) {
    Plane* plane = i->second;

    // reading calibrations
#if USE_DATABASE == 1
    plane->setDBpt(fDataBase);
#endif
    plane->OpenReference();
    if(fDoCalib)
      plane->ReadCalib(fRunTime);

    // association to a GeomPlane.
    if(geom && plane->NeedGeom())
      plane->SetGeometry( geom->GetPlane(plane->GetName()) );
    // plane initialization. must be done after geometry association
    // histograms are booked here. we create the good directory in
    // the root file and cd to it before calling the function

    if(fRootFile) {
      std::string dettype = plane->GetType();
      fRootFile->cd();

      // directory for this detector type
      if(!(fRootFile->GetKey(dettype.c_str())))
	fRootFile->mkdir(dettype.c_str());

      // directory for this detector
      TDirectory *subdir = dynamic_cast<TDirectory*>(fRootFile->Get(dettype.c_str()));
      subdir->mkdir(plane->GetName());

      std::string path = dettype;
      path += "/";
      path += plane->GetName();

      fRootFile->cd(path.c_str());
    }
    plane->Init(fTree);
    if(fRootFile) fRootFile->cd();
  }

  // load groups
  if (groupfile != "") LoadGroups(groupfile.c_str());

  // initialize groups
  for(std::map<std::string, Group*>::iterator ig=fGroupMap.begin(); ig!=fGroupMap.end(); ig++) {
    Group* group = ig->second;

#if USE_DATABASE == 1
    group->setDBpt(fDataBase);
#endif
    if(fRootFile) {
      std::string groupname = group->GetName();
      std::string grouptype = group->GetType();
      fRootFile->cd();
      if(groupname != fAllPlanesName) {
	std::string groupdirname = "GROUP_";
	groupdirname += groupname;

	// directory for this detector type
	TDirectory *subdir = (TDirectory*)fRootFile->Get(grouptype.c_str());
	if(!subdir) {
	  subdir = fRootFile->mkdir(grouptype.c_str());
	}
	// directory for this group
	subdir->mkdir(groupdirname.c_str());

	std::string path = grouptype;
	path += "/";
	path += groupdirname;
	fRootFile->cd(path.c_str());
      }
    }
    group->Init();
    if(fRootFile) fRootFile->cd();
  }

  // restarting decoding
  delete fManager;
  fManager=new CS::DaqEventsManager;
  for (LSiter ff = datafilelist.begin(); ff != datafilelist.end(); ff++) {
    fManager->AddDataSource((*ff).c_str());
  }
  fManager->SetMapsDir(mapfile);

  std::cout<<std::endl;
  std::cout<<"Starting decoding :"<<std::endl;
  std::cout<<"********************************************************************"<<std::endl;

  fManager->Print();
#endif // ifndef __CINT__

#if USE_DATABASE == 1
  fDataBase->DisconnectDB();
#endif

  // notify all groups and planes that a new run has started
  NewRun(fRunNumber);

  // workaround to the fact that this srcIDtoscan is deleted at each new run
  srcID_to_scan_to_reload = fManager->GetDaqOptions().GetSrcIDScan();

  std::cout <<"\nEvent types to read:";
  typedef std::set<int>::iterator evtIT;
  for (evtIT ii = fEvtTypesOK.begin(); ii != fEvtTypesOK.end(); ii++) {
    std::cout <<" "<<fEvtTypes[*ii];
  }
  std::cout<<std::endl;
  printf("Trigger mask used: 0x%x\n\n", fTMPattern);

  fIsInitialized=true;
}



Monitor::~Monitor() {
  volatile bool* thfinish = &fThreadFinished;

  if(fIsInitialized) {
    // stop thread
    if(fTh) {
      ThreadStopRequest();
      int ii;
      for (ii=0; ii<THREADWAITTIME; ii++) {  // we don't wait indefinitly...
        if (*thfinish) {
	  break;
	}
	gSystem->Sleep(4000);
      }
      ThreadStop();
    }

    // close rootfile
    if (fRootFile) {
      CloseRootFile();  // if it has not been already called by the thread
      fRootFile->Close();
      SafeDelete(fRootFile);
      fRootFile = 0;
      fTree = 0;
    }


#if USE_DATABASE == 1
    if (fDataBase) SafeDelete(fDataBase);
#endif

    std::cout<<"Time spent in decoding : \tRT "<<fBenchDecoding->RealTime()/static_cast<float>(fNevent)
	     <<" \tCPUT "<<fBenchDecoding->CpuTime()/static_cast<float>(fNevent)<<std::endl;
    std::cout<<"Time spent in clustering : \tRT "<<fBenchClustering->RealTime()/static_cast<float>(fNevent)
	     <<" \tCPUT "<<fBenchClustering->CpuTime()/static_cast<float>(fNevent)<<std::endl;
    std::cout<<"Time spent in tracking : \tRT "<<fBenchTracking->RealTime()/static_cast<float>(fNevent)
	     <<" \tCPUT "<<fBenchTracking->CpuTime()/static_cast<float>(fNevent)<<std::endl;

    //   exit(0); give crashes
  }
}


void Monitor::CloseRootFile() {

  if (fRootFileFlushed) return;  // this has already been done

  // protect this method if it is called by the 2 threads in the same time
  if (fThreadRun && thr_flag) TThread::Lock();

  if (!fRootFileFlushed) {

    typedef std::map<std::string,Plane*>::iterator PI;
    for(PI ip=fDetMap.begin(); ip!=fDetMap.end(); ip++) {
      (*ip).second->End();
    }

    TextOutput(fRunNumber);

    // Close properly the root file at the end of the run or at a break

    if (fRootFile) {
      fRootFile->cd();
      fRootFile->Write();
      std::cout<<"histograms saved in "<<fRootFile->GetName()<<std::endl;
      //     fRootFile->Close();
      //     SafeDelete(fTree);
      //     SafeDelete(fRootFile);
      //     fRootFile = 0;
      //     fTree = 0;
    }

    fRootFileFlushed = true;
  }

  if (fThreadRun && thr_flag) TThread::UnLock();
}


void Monitor::FillDetList() {

  std::cout<<std::endl;
  std::cout<<"Creating Plane objects :"<<std::endl
	   <<"********************************************************************"<<std::endl;

  NewGroup(fAllPlanesName,"ALL");
#ifndef __CINT__
  typedef CS::Chip::Maps::iterator MI;
  for(MI i=fManager->GetMaps().begin();i!=fManager->GetMaps().end();i++) {
    Plane* det;
    std::string detname=i->second->GetDetID().GetName();
    unsigned int detnum=i->second->GetDetID().GetNumber();
    det = NewDetector(detname.c_str(), detnum);
  }
  std::cerr<<std::endl;


  fManager->GetDaqOptions().Print();

#endif
}


Group* Monitor::NewGroup(const char* grpname, const char* grptype) {
  register Group* grp;

  // is this group already exist
  std::string sgrpname=grpname;
  std::map<std::string,Group*>::iterator i=Monitor::fGroupMap.find(sgrpname);
  if(i!=fGroupMap.end()) return i->second;

  if (grptype == 0) grptype = grpname;

  if (strncmp(grptype, "CEDAR", 5) == 0) { //CEDARs
    grp = new GroupCEDAR(grpname);
  }
  else if (strncmp(grptype, "HK", 2) == 0) { //Beam Killers
    grp = new GroupBeamKiller(grpname);
  }
  else if (strncmp(grptype, "H", 1) == 0) { // Trigger Hodoscope
    grp = new GroupTrigHodo(grpname);
  }
  else if (strncmp(grptype, "RA", 2) == 0) { // Rich APV
    grp = new GroupRiAPV(grpname);
  }
  else if (strncmp(grptype, "RM", 2) == 0) { // Rich MAPMT
    grp = new GroupRICH_MAPMT(grpname);
  }
  else if (strncmp(grptype, "FI", 2) == 0) { // SciFis japanese
    grp = new GroupScifiJ(grpname);
  }
  else if (strncmp(grptype, "SF", 2) == 0) { // SciFis german
    grp = new GroupScifiG(grpname);
  }
  else if (strncmp(grptype, "PS", 2) == 0) { // Mwpc's Star
    grp = new GroupMwpc(grpname);
  }
  else if (strncmp(grptype, "PA", 2) == 0) { // Mwpc's A
    grp = new GroupMwpc(grpname);
  }
  else if (strncmp(grptype, "PB", 2) == 0) { // Mwpc's B
    grp = new GroupMwpc(grpname);
  }
  else if (strncmp(grptype, "ST", 2) == 0) { // straws
    grp = new GroupStraws(grpname);
  }
  else if (strncmp(grptype, "DC", 2) == 0) { // drif chambers
    grp = new Group(grpname);  // no specific group for the moment
  }
  else if (strncmp(grptype, "MM", 2) == 0) { // micromegas
    grp = new GroupMumega(grpname);  // no specific group for the moment
  }
  else if (strncmp(grptype, "ALL", 3) == 0) { // groups all planes for future gui
    grp = new Group(grpname);  // no specific group for the moment
  }
  else if (strncmp(grptype, "DAQ", 3) == 0) { // display DAQ variables
    grp = new GroupDAQ(grpname);
  }
  else if (strncmp(grptype, "GEMtrack", 8) == 0) { // Group for GEM tracking
    grp = new GroupGEMTrack(grpname);
  }
  else if (strncmp(grptype, "GEM", 3) == 0) { // Group for 1 GEM station
    grp = new GroupGEM(grpname);
  }
  else if (strncmp(grptype, "PGEM", 4) == 0) { //Group for 1 PGEM station
    grp = new GroupPGEM(grpname);
  }
  else if (strncmp(grptype, "SILtrack", 8) == 0) { // Group for Silicon track
    grp = new GroupSiliTrack(grpname);
  }
  else if (strncmp(grptype, "SIL", 3) == 0) { // Group for 1 Silicon station
    grp = new GroupSili(grpname);
  }
  else if (strncmp(grptype, "MWA", 3) == 0) { // Group for Muon Wall A
    grp = new GroupMuonWallA(grpname);
  }
  else if (strncmp(grptype, "DR", 2) == 0) { // Group for RichWall
    grp = new GroupRichWall(grpname);
  }
  else if (strncmp(grptype, "RPD", 3) == 0) { // Group for RPD
    grp = new GroupRPD(grpname);
  }
  else if (strncmp(grptype, "EC", 2) == 0) { // Group for ECAL + FEM
    grp = new GroupECAL1(grpname);
  }
  else if (strncmp(grptype, "CA", 2) == 0) { // Group for Camera
    grp = new GroupCamera(grpname);
  }  
  else { // non recognized group
    std::cerr << "group type not recognized: " << grptype << " name: " << grpname<<std::endl;
    grp = new Group(grpname);
  }
  grp->SetType(grptype);

  Monitor::fGroupMap[sgrpname] = grp;
  std::cerr<<".";
  return grp;
}


Plane* Monitor::NewDetector(const char* detname,unsigned int detnum) {

  if(detnum >= fPlanes.size())
    throw CS::Exception("Monitor::NewDetector : DetID not valid");

  // this plane already exists
  std::string sdetname(detname);
  typedef std::map<std::string,Plane*>::iterator MI;
  MI i=Monitor::fDetMapRegistered.find(sdetname);
  if(i!=fDetMapRegistered.end()) return i->second;  // already registered
  MI inot=Monitor::fDetMapNotRegistered.find(sdetname);
  if(inot!=fDetMapNotRegistered.end()) return 0;  // already seen as not registered

  Plane *det=0;

#define  NCHAN fManager->GetMaps().GetWires(CS::DetID(detname))
  if(strncmp(detname, "TT",2)==0) { // Trigger time
    //AZ//fCatch2 = new PlaneCATCH2(detname,NCHAN,0,100);
    det=new Plane1V(detname,NCHAN,0,100);
    fManager->GetMaps().GetSrcIDs("TT.*",fManager->GetDaqOptions().GetSrcIDScan());
  }
  else if(strncmp(detname,"HMSC1",5) == 0) { // Special Trigger Hodoscope HMSC1 used for TCS phase, needed by APV based detectors,
    det=new PlaneMISC(detname,NCHAN,0,10000);  //  therefore always active
    fManager->GetMaps().GetSrcIDs("HM.*",fManager->GetDaqOptions().GetSrcIDScan());
  }
  //---------------------  ECAL0 ----------------------------------------
  
  else if(strncmp(detname, "EC00", 4 ) == 0)  { // ECAL0
    if((fDetInTree.find("ECAL0"))->second) {
      if ( strlen(detname) > 7 ) {
	if(strncmp(&detname[4], "F", 1 ) == 0) { // ECAL0 FEM
	  det=new PlaneECAL0(detname,2,2,2048,2048);
        } 
 	else  if(strncmp(&detname[4], "S", 1 ) == 0) { // ECAL0 SUM
	  det=new PlaneECAL0(detname,2,2,2048,2048);
        } 
        else { det=new PlaneECAL0(detname,33,27,2048,2048);  // main PLane 
        } 
        
        
      }  // if ( strlen(detname) > 7 )
      else {
	std::cerr << "Error in Monitor::NewDetector: detname "<<detname<<" is too short\n";
      }
      fManager->GetMaps().GetSrcIDs("EC00.*",fManager->GetDaqOptions().GetSrcIDScan());

    }  //if((fDetInTree.find("ECAL0"))->second)
  } //  if(strncmp(detname, "EC00", 4 ) == 0) { // ECAL0
  //-------------------------------------------------------------
  
  else if(strncmp(detname, "H", 1 ) == 0) {
    //-------------------------------------------------------------
    if(strncmp(detname, "HC01", 4 ) == 0) { // HCAL1
      if((fDetInTree.find("HCAL1"))->second) {
	if ( strlen(detname) > 7 ) {
	  if(strncmp(&detname[6], "T", 1 ) == 0) { // HCAL1 sum1 tdc
	    det=new PlaneHCALT(detname,8,8,-1024,2048);
	  }
	  else if(strncmp(&detname[6], "S", 1 ) == 0) { // HCAL1 sum2 adc
	    det=new PlaneHCAL1(detname,8,8,2048,2048);
	  }
	  else {
	    det=new PlaneHCAL1(detname,28,20,2048,2048);
	  }
	} else {
	  std::cerr << "Error in Monitor::NewDetector: detname "<<detname<<" is too short\n";
	}
        fManager->GetMaps().GetSrcIDs("HC01.*",fManager->GetDaqOptions().GetSrcIDScan());
      }
    }
    //-------------------------------------------------------------
   
    else if(strncmp(detname, "HC02", 4 ) == 0) { // HCAL2
      if((fDetInTree.find("HCAL2"))->second) {
	if ( strlen(detname) > 7 ) {
	  if(strncmp(&detname[6], "T", 1 ) == 0) { // HCAL2 sum1 tdc
	    det=new PlaneHCALT(detname,8,8,-1024,2048);
	  }
	  else if(strncmp(&detname[6], "S", 1 ) == 0) { // HCAL2 sum2 adc
	    det=new PlaneHCAL1(detname,8,8,2048,2048);
	  }
	  else {
	    det=new PlaneHCAL1(detname,22,10,2048,2048);
	  }
	} else {
	  std::cerr << "Error in Monitor::NewDetector: detname "<<detname<<" is too short\n";
	}
        fManager->GetMaps().GetSrcIDs("HC02.*",fManager->GetDaqOptions().GetSrcIDScan());
      }
    }
    else  { // Trigger Hodoscope
      if((fDetInTree.find("Hodoscopes"))->second) {
	if(strncmp("HP",detname,2) == 0)
	  det=new PlanePrimakoffHodo(detname,NCHAN,0,10000);
	else
	  if(strncmp("HK",detname,2) == 0)
	    det=new PlaneBeamKiller(detname,NCHAN,0,10000);
	  else
	    if(strncmp("HM01P",detname,5) == 0) {
	      // count TDC and ADC channels
	      std::set<CS::uint16> srcIds;
	      int nAdcChan(0);
	      int nTdcChan(0);
	      fManager->GetMaps().GetSrcIDs(detname, srcIds);
	      for (std::set<CS::uint16>::iterator srcIdIt=srcIds.begin(); srcIdIt!=srcIds.end(); srcIdIt++) {
		CS::Chip::DataID minDataId = CS::Chip::CreateDataID5(*srcIdIt,     0, 0, 0, 0);
		CS::Chip::Maps::const_iterator minMapIt = fManager->GetMaps().lower_bound(minDataId);
		CS::Chip::DataID maxDataId = CS::Chip::CreateDataID5((*srcIdIt)+1, 0, 0, 0, 0);
		CS::Chip::Maps::const_iterator maxMapIt = fManager->GetMaps().lower_bound(maxDataId);
		while (minMapIt!=maxMapIt) {
		  if (minMapIt->second->GetDetID().GetName()==detname) {
		    if (dynamic_cast<CS::ChipF1::Digit*>(minMapIt->second)   != 0) nTdcChan++;
		    if (dynamic_cast<CS::ChipSADC::Digit*>(minMapIt->second) != 0) nAdcChan++;
		  }
		  minMapIt++;
		}
	      }
	      // at the moment there are two entries per F1 chip in the daq maps
	      det=new PlaneTrigger_SADC_F1(detname, nAdcChan, nTdcChan/2, 0, 10000);
	    } else
	      if(strncmp("HMSC",detname,4))
		det=new PlaneTrigHodo(detname,NCHAN,0,10000);
	      else
		det=new PlaneMISC(detname,NCHAN,0,10000);

	fManager->GetMaps().GetSrcIDs("HI.*",fManager->GetDaqOptions().GetSrcIDScan());
	fManager->GetMaps().GetSrcIDs("HM.*",fManager->GetDaqOptions().GetSrcIDScan());
	fManager->GetMaps().GetSrcIDs("HL.*",fManager->GetDaqOptions().GetSrcIDScan());
	fManager->GetMaps().GetSrcIDs("HO.*",fManager->GetDaqOptions().GetSrcIDScan());
	fManager->GetMaps().GetSrcIDs("HP.*",fManager->GetDaqOptions().GetSrcIDScan());
	fManager->GetMaps().GetSrcIDs("HK.*",fManager->GetDaqOptions().GetSrcIDScan());
      }
    }
  }
  else if(strncmp(detname, "CE", 2 ) == 0) { // CEDARs
    if((fDetInTree.find("CEDARs"))->second) {
      det=new PlaneCEDAR(detname/*,NCHAN,0,10000*/); // uncomment to make the old version running
      fManager->GetMaps().GetSrcIDs("CE.*",fManager->GetDaqOptions().GetSrcIDScan());
    }
  }
  else if(strncmp(detname, "VT", 2 ) == 0) { // Veto Hodoscopes
    if((fDetInTree.find("Hodoscopes"))->second) {
      det=new PlaneVeto(detname,NCHAN,0,10000);
      fManager->GetMaps().GetSrcIDs("VT.*",fManager->GetDaqOptions().GetSrcIDScan());
    }
  }
  else if(strncmp(detname, "RI", 2 ) == 0) { // RICH
    if((fDetInTree.find("RICH_BORA"))->second) {
      det=new PlaneRICH(detname,16,72,72,2048,2048);
      fManager->GetMaps().GetSrcIDs("RI.*",fManager->GetDaqOptions().GetSrcIDScan());
    }
  }
  else if(strncmp(detname, "RM", 2 ) == 0) { // RICH_MAPMT
    if((fDetInTree.find("RICH_MAPMT"))->second) {
      //     det=new PlaneRICH_MAPMT(detname,48,48,-2000,1000);
      //      det=new PlaneRICH_MAPMT(detname,48,48,-30000,30000);
      det=new PlaneRICH_MAPMT(detname,48,48,-3000,3000);
      fManager->GetMaps().GetSrcIDs("RM.*",fManager->GetDaqOptions().GetSrcIDScan());
    }
  }
  else if(strncmp(detname, "RA", 2 ) == 0) { // Rich APV
    if((fDetInTree.find("RICH_APV"))->second) {
      //      PlaneRiAPV* riapv = new PlaneRiAPV(detname,NCHAN,500,500);
      PlaneRiAPV* riapv = new PlaneRiAPV(detname,6144,500,500);
      riapv -> SetRunMaps(&Monitor::fManager->GetMaps());
      det = riapv;
      fManager->GetMaps().GetSrcIDs("RA.*",fManager->GetDaqOptions().GetSrcIDScan());
    }
  }
  else if(strncmp(detname, "BM", 2 ) == 0) { // BMS
    if((fDetInTree.find("BMS"))->second) {
      det=new PlaneBMS(detname,NCHAN,0,20000);
      fManager->GetMaps().GetSrcIDs("BM.*",fManager->GetDaqOptions().GetSrcIDScan());
    }
  }

  else if(strncmp(detname, "VBX", 3 ) == 0) { // VetoBox - LeadGlass Detector
    if((fDetInTree.find("VetoBox"))->second) {   
      det=new PlaneVBOX(detname,12,8,2048,2048);
      fManager->GetMaps().GetSrcIDs("VBX1.*",fManager->GetDaqOptions().GetSrcIDScan());
    }
  }
  else if(strncmp(detname, "VB", 2 ) == 0) { // VetoBox - inner and front hodoscopes
    if((fDetInTree.find("VetoBox"))->second) {
      det=new PlaneVBHodo(detname,NCHAN,0,10000);
      fManager->GetMaps().GetSrcIDs("VB0.*",fManager->GetDaqOptions().GetSrcIDScan());
    }
  }

  else if(strncmp(detname, "EC01", 4 ) == 0) { // ECAL1
    if((fDetInTree.find("ECAL1"))->second) {
      //   	if((fDetInTree.find("FEM"))->second) {
      //       		if (strncmp(detname, "EC01FEM", 7 ) == 0) {  //ECAL1 FEMs
      //        	det=new PlaneFEM(detname,1,8,2048,2048);
      //                }
      //       }
      if(strncmp(detname, "EC01P00", 7 ) == 0) {  //ECAL1 GAMS
	det=new PlaneECAL1(detname,44,24,2048,2048);
      }
      else if(strncmp(detname, "EC01P01", 7 ) == 0) {  //ECAL1 MAINZ
	det=new PlaneECAL1(detname,22,26,2048,2048);
      }
      else if(strncmp(detname, "EC01P02", 7 ) == 0) {  //ECAL1 OLGA
	det=new PlaneECAL1(detname,16,20,2048,2048);
      }
      else if(strncmp(detname, "EC01P1S", 7 ) == 0) {  //ECAL1 trigger
	det=new PlaneECAL1(detname,44,24,2048,2048);
      }
      else if(strncmp(detname, "EC01P1T", 7 ) == 0) {  //ECAL1 trigger
	det=new Plane1V(detname,NCHAN,0,10000);
      }
      fManager->GetMaps().GetSrcIDs("EC01.*",fManager->GetDaqOptions().GetSrcIDScan());
    } else {    }
    if((fDetInTree.find("FEM"))->second) {
      if (strncmp(detname, "EC01FEM", 7 ) == 0) {  //ECAL1 FEMs
	det=new PlaneFEM(detname,1,8,2048,2048);
      }
    }     
  }

  else if(strncmp(detname, "EC02", 4 ) == 0) { // ECAL
    if((fDetInTree.find("ECAL2"))->second) {
      //      det=new PlaneHCAL1(detname,64,32,2048,2048);
      det=new PlaneECAL2(detname,64,48,2048,2048);
      fManager->GetMaps().GetSrcIDs("EC02.*",fManager->GetDaqOptions().GetSrcIDScan());
    } else {
      //      std::cerr << "Error in Monitor::NewDetector: detname "<<detname<<" is too short\n";
    }
  }

  else if(strncmp(detname, "SADC", 4 ) == 0) { // ECAL SADC
    if(0) {   // not treated up to now
      //       det=new PlaneSADC (does not exist...) (detname);
      fManager->GetMaps().GetSrcIDs("SADC.*",fManager->GetDaqOptions().GetSrcIDScan());
    } else {
      //      std::cerr << "Error in Monitor::NewDetector: detname "<<detname<<" is too short\n";
    }
  }

  else if(strncmp(detname, "RP01Q", 5 ) == 0) { // RPD SADC
    if((fDetInTree.find("RPD"))->second) {
      det=new PlaneRPD_SADC(detname);
      fManager->GetMaps().GetSrcIDs("RP01Q.*",fManager->GetDaqOptions().GetSrcIDScan());
    }
  }
  else if(strncmp(detname, "RP01T", 5 ) == 0) { // RPD F1
    if((fDetInTree.find("RPD"))->second) {
      det=new PlaneRPD_F1(detname,NCHAN,0,100);
      fManager->GetMaps().GetSrcIDs("RP01T.*",fManager->GetDaqOptions().GetSrcIDScan());
    }
  }

  //   else if(strncmp(detname, "RP01T", 5 ) == 0) { // RPD TDC
  //     if((fDetInTree.find("RPD"))->second) {
  //       det=new Plane1V(detname,NCHAN,0,60000);
  //       fManager->GetMaps().GetSrcIDs("RP01T.*",fManager->GetDaqOptions().GetSrcIDScan());
  //     }
  //   }

  else if(strncmp(detname, "SW", 2 ) == 0) { // Sandwich
    if((fDetInTree.find("Sandwich"))->second) {
      det=new PlaneSandwich(detname);
      fManager->GetMaps().GetSrcIDs("SW.*",fManager->GetDaqOptions().GetSrcIDScan());
    }
  }

  else if(strncmp(detname, "FI", 2 ) == 0) { // SciFis

    std::string s(detname);
    //int pos=s.find("0");
    char nums[3];

    s.copy(nums,2,2);
    nums[2] = '\000';

    int num=atoi(nums);

    switch (num) {
    case 1:
    case 2:
    case 3:
    case 4:
      if((fDetInTree.find("SciFis(Japan)"))->second) {
	det=new PlaneScifiJ(detname,NCHAN,-8400,500);
        fManager->GetMaps().GetSrcIDs("FI0[1-4].*",fManager->GetDaqOptions().GetSrcIDScan());
      }
      break;
    case 5:
    case 6:
    case 7:
    case 8:
      if((fDetInTree.find("SciFis(Germany)"))->second) {
	if(s.find("6V")!=string::npos) {
	  PlaneScifiJ *p=new PlaneScifiJ(detname,NCHAN,-7000,10000);
	  //?//p->SetHRTime(true);
          fManager->GetMaps().GetSrcIDs("FI06V.*",fManager->GetDaqOptions().GetSrcIDScan());
	  det = p;
	}
	else {
	  det=new PlaneScifiG(detname,NCHAN,-7000,10000);
	  fManager->GetMaps().GetSrcIDs("FI0[5-8].*",fManager->GetDaqOptions().GetSrcIDScan());
	}
      }
      break;
    case 12:
    case 13:
    case 14:
    case 15:
      {
	if((fDetInTree.find("SciFis(Germany)"))->second) {
	  PlaneScifiJ *p=new PlaneScifiJ(detname,NCHAN,-7000,10000);
	  //?//p->SetHRTime(true);
	  fManager->GetMaps().GetSrcIDs("FI15.*",fManager->GetDaqOptions().GetSrcIDScan());
	  det = p;
	  break;
	}
      }
      break;
    case 35:
      if((fDetInTree.find("SciFis(Japan)"))->second) {
	det=new PlaneScifiJ(detname,NCHAN,-8500,1500);
	fManager->GetMaps().GetSrcIDs("FI35.*",fManager->GetDaqOptions().GetSrcIDScan());
      }
      break;
    case 55:
      if((fDetInTree.find("SciFis(Warsaw)"))->second) {
	det=new PlaneScifiJ(detname,NCHAN,-8400,500);
        fManager->GetMaps().GetSrcIDs("FI55.*",fManager->GetDaqOptions().GetSrcIDScan());
      }
      break;
    default:
      break;
    }
  }
  else if(strncmp(detname, "PS", 2 ) == 0) { // Mwpc's Star
    if((fDetInTree.find("MWPC's"))->second) {
      det=new PlaneMwpc(detname,800,-1050,2000);
      fManager->GetMaps().GetSrcIDs("PS.*",fManager->GetDaqOptions().GetSrcIDScan());
    }
  }
  else if(strncmp(detname, "PA", 2 ) == 0) { // Mwpc's A
    if((fDetInTree.find("MWPC's"))->second) {
      det=new PlaneMwpc(detname,800,-1050,2000);
      fManager->GetMaps().GetSrcIDs("PA.*",fManager->GetDaqOptions().GetSrcIDScan());
    }
  }
  else if(strncmp(detname, "PB", 2 ) == 0) { // Mwpc's B
    if((fDetInTree.find("MWPC's"))->second) {
      det=new PlaneMwpc(detname,800,-1050,2000);
      fManager->GetMaps().GetSrcIDs("PB.*",fManager->GetDaqOptions().GetSrcIDScan());
    }
  }
  else if( strncmp(detname, "DW", 2 ) == 0 ) { // W45
    if((fDetInTree.find("W45"))->second) {
      det=new Plane1V(detname,800,-8000,10000);
      fManager->GetMaps().GetSrcIDs("DW.*",fManager->GetDaqOptions().GetSrcIDScan());
    }
  }
  else if(strncmp(detname, "MB", 2 ) == 0) { // Muon Walls
    if((fDetInTree.find("MuonWallB"))->second) {
      det=new PlaneMuonWallB(detname,800,-9000,2000);
      fManager->GetMaps().GetSrcIDs("MB.*",fManager->GetDaqOptions().GetSrcIDScan());
    }
  }
  else if(strncmp(detname, "MA", 2 ) == 0) { // Muon Walls
    if((fDetInTree.find("MuonWallA"))->second) {
      det=new PlaneMuonWallA(detname,800,-9000,2000);
      fManager->GetMaps().GetSrcIDs("MA.*",fManager->GetDaqOptions().GetSrcIDScan());
    }
  }
  else if(strncmp(detname, "DC", 2 ) == 0) { // drift chambers
    if((fDetInTree.find("DriftChambers"))->second) {
      if(strncmp(detname, "DC05", 4 ) == 0) { // dc05 uses planedriftchamber2v
	det=new PlaneDriftChamber2V(detname,NCHAN,-11500,2500);
	fManager->GetMaps().GetSrcIDs("DC.*",fManager->GetDaqOptions().GetSrcIDScan());
      }
      else
	det=new PlaneDriftChamber(detname,NCHAN,-11500,2500);
      fManager->GetMaps().GetSrcIDs("DC.*",fManager->GetDaqOptions().GetSrcIDScan());
    }
  }
  else if(strncmp(detname, "MM", 2 ) == 0) { // micromegas
    if((fDetInTree.find("Micromegas"))->second) {
      if(fRunNumber < 14000) {
	// xlsat was used in 2001
	if(strncmp(detname, "MM01X1", 6) == 0)
	  det=new PlaneMumega(detname,1152,-10000,2500,1);
	else
	  det=new PlaneMumega(detname,1024,-10000,2500);
      }
      else
	det=new PlaneMumega(detname,1024,-10000,2500);
      fManager->GetMaps().GetSrcIDs("MM.*",fManager->GetDaqOptions().GetSrcIDScan());
    }
  }
  else if(strncmp(detname, "MP", 2) == 0) { //PixelMumegas
    if((fDetInTree.find("PixelMicromegas"))->second) {
      
      int val=0;
      if (strncmp(&detname[4],"P",1)==0) 
	val=1;
      else if (strncmp(&detname[4],"M",1)==0)
	val=2;
      PlanePMumega *pmumega = new PlanePMumega(detname,NCHAN,500,500,val); /* distinguish between pixel and strip readout */
      pmumega -> SetRunMaps(&Monitor::fManager->GetMaps());
      det = pmumega;
      fManager->GetMaps().GetSrcIDs("MP.*",fManager->GetDaqOptions().GetSrcIDScan());
    }
  }
  else if(strncmp(detname, "ST", 2 ) == 0) { // Straw tubes
    if((fDetInTree.find("StrawTubes"))->second) {
      det = new PlaneStrawTubes(detname, -12250, 1250);
      fManager->GetMaps().GetSrcIDs("ST.*",fManager->GetDaqOptions().GetSrcIDScan());
    }
  }
  else if(strncmp(detname, "GM", 2 ) == 0) { // GEM's
    if((fDetInTree.find("GEMs"))->second) {
      PlaneGEM *gem = new PlaneGEM(detname,NCHAN,500,500);
      gem -> SetRunMaps(&Monitor::fManager->GetMaps());
      det = gem;
      fManager->GetMaps().GetSrcIDs("GM.*",fManager->GetDaqOptions().GetSrcIDScan());
    }
  }
  else if(strncmp(detname, "GP", 2) == 0) { //PixelGEMs
    if((fDetInTree.find("PixelGEMs"))->second) {
      PlanePGEM *pgem = new PlanePGEM(detname,NCHAN,500,500,(strncmp(&detname[4],"P",1)==0)/*distinguish between pixel and strip readout*/);
      pgem -> SetRunMaps(&Monitor::fManager->GetMaps());
      det = pgem;
      fManager->GetMaps().GetSrcIDs("GP.*",fManager->GetDaqOptions().GetSrcIDScan());
    }
  }
  else if(strncmp(detname, "SI", 2 ) == 0) { // Silicon's
    if((fDetInTree.find("Silicons"))->second) {
      PlaneSili *sili = new PlaneSili(detname,NCHAN,500,500);
      sili -> SetRunMaps(&Monitor::fManager->GetMaps());
      det = sili;
      fManager->GetMaps().GetSrcIDs("SI.*",fManager->GetDaqOptions().GetSrcIDScan());
    }
  }
  else if(strncmp(detname, "SC", 2 ) == 0) { // Scalers
    if((fDetInTree.find("Scalers"))->second) {
      det=new PlaneScaler(detname,1);
      fManager->GetMaps().GetSrcIDs("SC.*",fManager->GetDaqOptions().GetSrcIDScan());
    }
  }
  else if(strncmp(detname, "DR", 2 ) == 0) { // RICH Wall
    if((fDetInTree.find("RichWall"))->second) {

      if( (strncmp(detname, "DR01X1", 6 ) == 0) || (strncmp(detname, "DR01X2", 6 ) == 0) ||
	  (strncmp(detname, "DR02X1", 6 ) == 0) || (strncmp(detname, "DR02X2", 6 ) == 0)   ){

	det=new PlaneRichWall(detname,592,-3000,3000);
	fManager->GetMaps().GetSrcIDs("RW.*",fManager->GetDaqOptions().GetSrcIDScan());

      }else if( (strncmp(detname, "DR01Y1", 6 ) == 0) || (strncmp(detname, "DR01Y2", 6 ) == 0) ||
		(strncmp(detname, "DR02Y1", 6 ) == 0) || (strncmp(detname, "DR02Y2", 6 ) == 0)   ){

	det=new PlaneRichWall(detname,416,-3000,3000);
	fManager->GetMaps().GetSrcIDs("RW.*",fManager->GetDaqOptions().GetSrcIDScan());

      }
    }
  }
  else if(strncmp(detname, "TCSphase", 8 ) == 0) { // TCS phase
    // should always be activate as it is needed for time reconstruction of
    // all APV detectors plus ECALS
    if (NCHAN != 1) {
      std::cerr << "Monitor::NewDetector: Exactly one channel must be mapped onto TCSphase (" << NCHAN << " actually are). This is a major flaw in your mapping!" << std::endl;
      exit(1);
    }
    det = new PlaneTcsPhase(detname);
    fManager->GetMaps().GetSrcIDs("TCSphase",fManager->GetDaqOptions().GetSrcIDScan());
  }
  else if(strncmp(detname, "CA", 2 ) == 0) { // Camera detector
    if((fDetInTree.find("Camera"))->second) {
      det=new PlaneCamera(detname);
      fManager->GetMaps().GetSrcIDs("CA.*",fManager->GetDaqOptions().GetSrcIDScan());
    }
  }
  else {
    std::cout<<"Detector type not recognized: "<<detname<<". See Monitor::NewDetector"<<std::endl;
    det=new Plane1V(detname,NCHAN,0,60000);
    fManager->GetMaps().GetSrcIDs(detname,fManager->GetDaqOptions().GetSrcIDScan());
  }
  Monitor::fPlanes[detnum]=det;

  // #warning "Straw-specific code?"
  // yes but no more useful as straw scheme has changed
  // anyway it does not give problems and I don't remember what I did so lets keep it. Damien
  if(det) {
    det->SetTree(fTree);
    // change due to PlaneStrawTubes which change the detname given ("ST03Y1ua" -> "ST03Y1__")
    if ( ! fDetMap[det->GetName()]) GetGroup(fAllPlanesName)->Add(det);
    Monitor::fDetMap[det->GetName()]=det;
    Monitor::fDetMapRegistered[detname]=det;
    //    std::cerr<<"detname entree: "<<detname<<" sortie: "<<det->GetName()<<std::endl;
    std::cerr<<".";
  } else {
    Monitor::fDetMapNotRegistered[detname]=0;
  }

  return det;
}


void Monitor::Run() {

#ifndef __CINT__
  try {
    fLastDecodedEventnum=0;
    while(!flag_end) {

      // checking if the stopgo button is true (suspend events reading)
      // wait in this case
      if (fStopGo) {
        gSystem->Sleep(1000);
        continue;
      }


      // resetting all planes at the beginning of the event
      typedef std::map<std::string,Plane*>::iterator DI;
      for(DI i=Monitor::fDetMap.begin(); i!=Monitor::fDetMap.end(); i++ ) {
	(*i).second->Reset();
      }

      // do we need to reset groups ?
      if(Plane::ModifiedSemaphore()) {
	Plane::ModifiedAcknowledge();
	for(std::map<std::string,Group*>::iterator ig=fGroupMap.begin(); ig!=fGroupMap.end(); ig++) {
	  ig->second->ResetIfNeeded();
	}
	for(std::map<std::string,Plane*>::iterator ip=fDetMap.begin(); ip!=fDetMap.end(); ip++)
	  ip->second->Modified(0);
      }

      if( !fManager->ReadEvent() ) {
        ThreadStopRequest();
        break;
      }

      CS::DaqEvent &event = fManager->GetEvent();

      if ( fStopAtEndRun && (event.GetType() == 2)) {
	std::cerr << "Begin of run event and StopAtEndRun set -> end of data read-out for run " << fRunNumber <<std::endl;
        ThreadStopRequest();
	break;
      }

      // get some information about this event from the decoding library
      fSpill=event.GetBurstNumber();
      fEvtinrun=event.GetEventNumberInRun();
      fEvtinspl=event.GetEventNumberInBurst();
      fTimesec=event.GetTime().first;
      fTimemic=event.GetTime().second;
      fEvtType=event.GetType();
      //      std::cerr << " event number = "<<fEvtinrun<<"\n";
      //       CS::DaqEvent::Header head=event.GetHeader();
      //       fPattern=head.typeAttribute[1];
      fPattern=event.GetTrigger();    // for new version of decoding lib (new DATE v5)
      //       fAttributes=head.typeAttribute[0];
      // for new version of decoding lib (new DATE v5)
      //       fAttributes= head.GetVersion()<0xffff ? head.GetHeaderOld().typeAttribute[1] : head.GetHeader36().event_trigger_pattern[0];
      fAttributes= 0;
      // std::cerr<<"spill "<<fSpill<<" evtnb "<<fEvtinrun<<" fTimesec "<<fTimesec<<" fEvtType "<<fEvtType<<std::endl;

      if ( event.GetRunNumber() != fRunNumber ) {
	// run number has changed...

	TextOutput(fRunNumber);

        if ( fStopAtEndRun ) {
	  std::cerr << "New run seen and StopAtEndRun set -> end of data read-out for run " << fRunNumber <<std::endl;
          ThreadStopRequest();
	  break;
        }

	register unsigned int oldrunnumber = fRunNumber;
	fRunNumber = event.GetRunNumber();
        Monitor::fPlanes.resize(fManager->GetDetectors().size(),NULL);
	std::cerr << "New run number " << fRunNumber << " (last one was "<<oldrunnumber<<")"<<std::endl;

	NewRun(fRunNumber);
	// workaround to the fact that this srcIDtoscan is deleted at each new run
        fManager->GetDaqOptions().GetSrcIDScan() = srcID_to_scan_to_reload;
      }

      std::set<int>::iterator tok = fEvtTypesOK.find(fEvtType);
      if(tok == fEvtTypesOK.end()) continue;

      // evt type is accepted

      unsigned int tPattern(fPattern);
      if (fTMLowerOnly)
        tPattern &= 0xffff;
      bool accepted(true);
      if (fTMStrict)
        accepted = (tPattern==fTMPattern);
      else
        accepted = (tPattern&fTMPattern);
      if(!accepted && fPattern)
        continue;

      // trigger mask is accepted...

      // throw away events that are too close
#warning  Is that normal ?: fAttributes=fEvtinspl - fLastDecodedEventnum;
      // fAttributes is filled above with head.typeAttribute[0] ????
      fAttributes=fEvtinspl - fLastDecodedEventnum;
      if (fLastDecodedEventnum > fEvtinspl) // we are in a new spill
	fLastDecodedEventnum=0;
      if (fEvtinspl - fLastDecodedEventnum < fMinEvtspacing)
	continue;
      else
	fLastDecodedEventnum=fEvtinspl;

      //...decode
      fBenchDecoding->Start(false);
      const bool decoded = fManager->DecodeEvent();
      fBenchDecoding->Stop();
      if( !decoded )
	{
	  if (! fReadFailedEvt) {
	    std::cerr << "Event "<<fEvtinrun<<": decoding has failed! skipped\n";
	    continue;
	  } else {
	    //           std::cerr << "Event "<<fEvtinrun<<": decoding has failed! using it anyway\n";
	  }
	}

      // at last, count this event

      //    std::cout<<"\n\n"<<fNevent<<" events read" << std::endl;

      fNevent++;
      if((fNevent%1000)==0) {
	std::cout<<"\n\n"<<fNevent<<" events read" << std::endl;
      }
      if ( fNeventMax != 0 && fNevent >= fNeventMax) ThreadStopRequest();

      //... set online filter data...
      CS::OnlineFilter& of=event.GetOnlineFilter();
      fFilterBits=0;
      if (of.IsMonitoring()) fFilterBits |= 0x00000001;
      if (of.IsAccepted()) fFilterBits |= 0x00000002;

      // Pass digits from the decoding library to coool
      for(CS::Chip::Digits::const_iterator it=fManager->GetEventDigits().begin();
          it!=fManager->GetEventDigits().end(); it++ ) {
	register unsigned int detnum = it->first.GetNumber();
	if(detnum < fPlanes.size())
	  if(fPlanes[detnum])
	    fPlanes[detnum]->StoreDigit(it->second);
      }


      //...Planes analysis
      fBenchClustering->Start(false);
      for(DI di=Monitor::fDetMap.begin(); di!=Monitor::fDetMap.end(); di++ ) {
	// digitize
	//	std::cout<<"Detector analysed: "<<di->second->GetName()<<std::endl;
	di->second->UpdateRateCounter();
	di->second->EndEvent(event);
	// clusterize
	if(fDoClustering) di->second->Clusterize();
      }
      fBenchClustering->Stop();

      //...Groups analysis
      fBenchTracking->Start(false);
      // std::cout<<"TSV: Event: "<<fNevent<<" \n"; // TSV, don't remove!!!
      typedef std::map<std::string,Group*>::iterator GI;
      for(GI gi=Monitor::fGroupMap.begin(); gi!=Monitor::fGroupMap.end(); gi++ ) {
	(*gi).second->UpdateRateCounter();
	(*gi).second->EndEvent(event);
	if(fDoTracking) gi->second->Tracking();
      }
      fBenchTracking->Stop();

      if(fThreadRun)
	if (thr_flag) TThread::Lock();
      if (fTree) fTree->Fill();
      if(fThreadRun)
	if (thr_flag) TThread::UnLock();
    }
  }
  catch( CS::DaqEvent::ExceptionEndOfStream ) {
    // This is the normal exit from the loop
    std::cout<<"End of file "<<std::endl;
    //     if(fThreadRun)
    //       ThreadStop();
    //     else return;
    ThreadStopRequest();
    return;
  }

  // Something is wrong...
  // Print error message and finish this data file.
  catch(CS::DaqError &e) {
    //std::cerr<<e.what()<<std::endl;
    ThreadStopRequest();
    return;
  }

  catch( const std::exception &e ) {
    std::cerr << "exception:\n" << e.what() << "\n";
    ThreadStopRequest();
    return;
  }
  catch( const char * s ) {
    std::cerr << "exception:\n" << s << "\n";
    exit(1);
  }
  catch( ... ) {
    std::cerr << "Oops, unknown exception!\n";
    exit(1);
  }
  ThreadStopRequest();
#endif

  CS::Exception::PrintStatistics();
  CS::Chip::PrintStatistics();
  gSystem->Sleep(1000);
}


void Monitor::ResetTree() {

  if (thr_flag) TThread::Lock();
  delete fTree;
  fTree=new TTree("T","COMPASS monitoring");
  if (thr_flag) TThread::UnLock();
}

// #ifndef __CINT__
// void Monitor::Decode(CS::Chip& chip) {
//
//   try {
//     // Loop over all data lines in given chip.
//     // Note that argument 'maps' for constructor of CS::Chip::DataIter class
//     // is _optional_, you may skip it. Read more in the documentation page.
//
//     CS::Chip::Digits digits;
//     chip.Decode(Monitor::fManager->GetMaps(),digits,fManager->GetDaqOptions());
//
//     typedef CS::Chip::Digits::const_iterator DI;
//     for(DI it=digits.begin(); it!=digits.end(); it++ ) {
//       unsigned int detnum = it->first.GetNumber();
//
//       // Pass the digit for trigger time decoding.
//       CS::ChipF1::GetTT().DataAdd(*it->second);
//
//       if(detnum < fPlanes.size())
// 	if(fPlanes[detnum])
// 	  fPlanes[detnum]->StoreDigit(it->second);
//
//     }
//   }
//   catch( CS::Exception &e ) {
//     //e.Print();
//   }
//   catch( ... ) {
//     std::cerr << "Unknown error in decoding!\n";
//   }
// }
// #endif


int Monitor::ThreadStart() {

  fThreadRun=true;
  if(!fTh){
    std::cout<<"Creating data treatment thread"<<std::endl;
    gSystem->Sleep(1000);
    fTh= new TThread("memberfunction",
		     (void(*) (void *))&Thread0,
		     (void*) this);
    std::cout<<"Starting data treatment thread"<<std::endl;
    gSystem->Sleep(1000);
    fTh->Run();
    return 0;
  }
  return 1;
}


void Monitor::ThreadStopRequest() {
  flag_end = true;
}


int Monitor::ThreadStop() {
  // stop all active threads
  ThreadStopRequest();
  if(fThreadRun){
    gSystem->Sleep(1000); // just to be sure...
    TThread::Delete(fTh);
    //delete fTh;
    //    fTh=0;
    fThreadRun=false;    // aborting flag
    return 0;
  }
  return 1;
}


void Monitor::Thread0(void* arg) {
  // thread function which calls user specified action Func0

  Monitor *th=(Monitor*)arg;
  TThread::SetCancelOn();
  TThread::SetCancelDeferred();
  int meid=TThread::SelfId(); // get pthread id
  std::cout << "\nThread 0, id:" <<meid<< " is running..\n"<<std::endl;
  TThread::CancelPoint();
  th->Run(); // call the user defined threaded function
  //  th->CloseRootFile();
  if (th->fClosingWindow == false) th->CloseRootFile();
  th->fThreadFinished = true;

  //   th->ThreadStop();
}


void Monitor::ResetHistos() {
  fNevent = 0;

  typedef std::map<std::string,Plane*>::iterator PI;
  for(PI ip=fDetMap.begin(); ip!=fDetMap.end(); ip++)
    (*ip).second->ResetHistograms();

  typedef std::map<std::string,Group*>::iterator GI;
  for(GI ig=fGroupMap.begin(); ig!=fGroupMap.end(); ig++)
    (*ig).second->ResetHistograms();

  std::cout << "\nAll histograms have been cleared"<<std::endl;
}


void Monitor::NewRun(int runnumber) {
  fBenchDecoding->Reset();
  fBenchClustering->Reset();
  fBenchTracking->Reset();
  if(fRunClear) { fNevent = 0; }

  if (thr_flag) TThread::Lock();

  typedef std::map<std::string,Plane*>::iterator PI;
  for(PI ip=fDetMap.begin(); ip!=fDetMap.end(); ip++) {
    (*ip).second->NewRun(runnumber);
  }
  typedef std::map<std::string,Group*>::iterator GI;
  for(GI ig=fGroupMap.begin(); ig!=fGroupMap.end(); ig++) {
    (*ig).second->NewRun(runnumber);
  }

  if (fRunClear) ResetHistos();

  if (thr_flag) TThread::UnLock();
}


TH1* Monitor::GetHisto(const std::string hName) {
  register TH1* th1 = 0;

  for (register unsigned int ii=0; ii<fPlanes.size(); ii++) {
    register Plane* pl = fPlanes[ii];
    if (pl) {
      th1 = pl->GetHisto(hName);
      if (th1) return th1;
    }
  }

  for(std::map<std::string,Group*>::iterator ig=fGroupMap.begin(); ig!=fGroupMap.end(); ig++) {
    register Group* gr = ig->second;
    if (gr) {
      th1 = gr->GetHisto(hName);
      if (th1) return th1;
    }
  }

  return 0;
}


void Monitor::TextOutput(int runnumber) {

  if ((!fDoTextOutput) && fTextDir.empty() && fTextFile.empty()) return;
  if(fTextDir.empty()) fTextDir += ".";

  // open output file
  if ( fTextFile.empty() || (fTextFile == "" )) {
    char outfilename[100];
    sprintf(outfilename,"%s/coool_output_%d.dat",fTextDir.c_str(),runnumber);
    fTextFile = outfilename;
  }

  ofstream out(fTextFile.c_str());
  if(!out) {
    std::cerr<<"WARNING cannot open text output file "<<fTextFile<<std::endl;
    return;
  }
  else
    std::cout<<"text output saved in "<<fTextFile<<std::endl;

  typedef std::map<std::string,Plane*>::iterator PI;
  for(PI ip=fDetMap.begin(); ip!=fDetMap.end(); ip++) {
    out<<"BEGIN Plane "<<((*ip).second->GetName())<<std::endl;
    (*ip).second->TextOutput(out);
    out<<"END Plane "<<((*ip).second->GetName())<<std::endl<<std::endl;
  }
}


void Monitor::CheckEvtTypes(const char* evttypesstr) {
  typedef std::map<int,std::string>::iterator evtIT;

  const char *evtp = evttypesstr;
  const char *evtmp;
  std::string ledet;
  std::set<int> evttypesok;

  while ((evtmp = strsep((char **)&evtp, ","))) {
    ledet = evtmp;
    for (evtIT ii = fEvtTypes.begin(); ii != fEvtTypes.end(); ii++) {
      if (ledet == ii->second) { evttypesok.insert(ii->first); }
    }
  }
  CheckEvtTypes(evttypesok);
}



namespace MonitorFile {  // private namespace for the groups description decoding

  class UserData
  {
  public:
    UserData(void)
    { Clear(); }
    void    Clear(void)
    { monitor = 0; read_mode = 0; ClearGroup();
    }
    void    ClearGroup(void)
    { groupname = ""; grouptype = ""; group = 0; ClearPlane();
    }
    void    ClearPlane(void)
    { planename = ""; planetype = "";
    }
    Monitor* monitor;
    int      read_mode;
    std::string   groupname;
    std::string   grouptype;
    Group*   group;
    std::string   planename;
    std::string   planetype;
  };



  // XML_StartElementHandler
  void startElement(void *userData, const char *name, const char **atts)
  {
    UserData &d = *reinterpret_cast<UserData*>(userData);
    std::string stname(name);

    for( register int i=0; atts[i]!=0; i+=2 ) {
      if (atts[i+1]==0 )
	{
#ifndef __CINT__
	  CS::Exception("XML:startElement(): an attribute is zero!!");
#endif
	  return;
	}
    }

    if (stname == "groupcontents") {
      //    register int i;
      d.Clear();
      d.read_mode = 1;
      //     for (i=0; atts[i]!=0; i+=2 ) {
      //       std::string opt(atts[i]), val(atts[i+1]);
      //       if (opt == "version") vers = atof(atts[i+1]);
      //     }
      return;
    }

    if (d.read_mode == 1 && stname == "group") {
      register int i;
      d.ClearGroup();
      for (i=0; atts[i]!=0; i+=2 ) {
	std::string opt(atts[i]), val(atts[i+1]);
	if (opt == "name") d.groupname = atts[i+1];
	if (opt == "type") d.grouptype = atts[i+1];
      }
      register Group* grp = 0;
      if (d.groupname != "") {
	grp = d.monitor->GetGroup(d.groupname.c_str());
	if (grp) d.group = grp;
	else {
	  grp = d.monitor->NewGroup(d.groupname.c_str(), d.grouptype.c_str());
	  d.group = grp;
	}
      }
      return;
    }

    if (d.read_mode == 1 && stname == "plane") {
      register int i;
      d.ClearPlane();
      if (d.group == 0) return;
      for (i=0; atts[i]!=0; i+=2 ) {
	std::string opt(atts[i]), val(atts[i+1]);
	if (opt == "name") d.planename = atts[i+1];
	if (opt == "type") d.planetype = atts[i+1];
      }
      register Plane* pln = 0;
      if (d.planename != "") {
	pln = d.monitor->GetPlane(d.planename.c_str());
	if (pln) d.group->Add(pln);
	//      else {
	//  std::cerr << "Plane "<<d.planename<<" not known in Group "<<d.groupname<<std::endl;
	//}
      }
      return;
    }
  }


  /* XML_EndElementHandler */
  void endElement(void *userData, const char *name)
  {
    UserData &d = *reinterpret_cast<UserData*>(userData);
    std::string stname(name);
    if (stname == "groupcontents") d.read_mode = 0;
    // std::cerr << "endElement  name=" << name <<std::endl;
  }


  /* XML_DefaultHandler */
  void default_handler(void *userData,const XML_Char *s,int len)
  {
    //   UserData &d = *reinterpret_cast<UserData*>(userData);
    //   char buf[500];
    //   for (register int i=0; i<len; i++) buf[i]=s[i];
    //   buf[len]=0;
    // std::cerr << "*** default_handler  len=" << len << " xml_char= " << buf <<std::endl;
    // std::cerr << "******\n";
    // Search and remove empty strings
    //   for( int j=0; j<len; j++ )
    //     if( s[j]!=' ' && s[j]!='	' && s[j]!='\n' )
    //       std::cerr << "    mot vide\n";
    // Empty string. Ignore it.
    return;
  }

}  // end of private namespace




Bool_t Monitor::LoadGroups(const char* filename)
{
  using namespace MonitorFile;
  XML_Parser parser = XML_ParserCreate(NULL);
  UserData d;

  const unsigned buf_size=1000000;
  char *buf = new char[buf_size];
  XML_SetUserData                       (parser, &d);
  XML_SetElementHandler                 (parser, startElement, endElement);
  //   XML_SetCharacterDataHandler           (parser, characterdata_handler);
  //   XML_SetProcessingInstructionHandler   (parser, instruction_handler);
  XML_SetDefaultHandler                 (parser, default_handler);
  //   XML_SetCommentHandler                 (parser, comment_handler);


  ifstream in(filename);
  if (!in) {
#ifndef __CINT__
    CS::Exception("Can not read groups settings file\n");
#endif
    return false;
  }

  std::cout<<std::endl;
  std::cout<<"Creating Group objects :"<<std::endl;
  std::cout<<"********************************************************************";
  std::cout<<std::endl;

  d.monitor = this;
  try
    {
      int done=0;
      do
	{
	  in.read(buf,buf_size);
	  if( ((unsigned int)in.gcount())==buf_size )
#ifndef __CINT__
	    throw CS::Exception("Monitor::LoadGroups(): too short buffer.");
#endif

	  done = ((unsigned int)in.gcount()) < buf_size;
	  // std::cerr << "appel XML_Parse: in.gcount()=" << in.gcount() << " done=" <<done <<std::endl;
	  if( !XML_Parse(parser, buf, in.gcount(), done) )
#ifndef __CINT__
	    CS::Exception("Monitor::LoadGroups(): Error: %s  at line %d\n",
			  XML_ErrorString(XML_GetErrorCode(parser)),
			  XML_GetCurrentLineNumber(parser)).Print();
#endif
	} while (!done);
    }
  catch(...)
    {
      delete [] buf;
      throw;
    }

  delete [] buf;
  in.close();
  XML_ParserFree(parser);
  std::cout<<std::endl;
  return true;
}


void Monitor::ControlPanel(const TGWindow *p, const TGWindow *main) {
  if (!fControlPanel) {
    fControlPanel = new MonitorPanel(p, main, 100, 100,
				     this, fTMPattern, fTMLowerOnly, fTMStrict);

    typedef std::set<int>::iterator IT;
    for(IT it = fEvtTypesOK.begin(); it!= fEvtTypesOK.end(); it++)
      fControlPanel->CheckEvtType(*it);
  }
}








