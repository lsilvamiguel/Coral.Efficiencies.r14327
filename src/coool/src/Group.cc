#include "TThread.h"

#include "Group.h"
#include "GroupPanel.h"
#include "GroupMumega.h"
#include "GroupMwpc.h"
#include "GroupScifiG.h"
#include "GroupScifiJ.h"
#include "GroupGEM.h"
#include "GroupPGEM.h"
#include "GroupSili.h"
#include "GroupDAQ.h"
#include "GroupTrigHodo.h"
#include "GroupMuonWallA.h"
#include "GroupStraws.h"
#include "GroupBeamKiller.h"
#include "GroupCEDAR.h"
#include "GroupRiAPV.h"
#include "GroupRICH_MAPMT.h"
#include "GroupRW.h"
#include "GroupRPD.h"
#include "GroupCamera.h"
ClassImp(Group);

using namespace std;


std::string Group::defaultRefFileName = "/afs/cern.ch/compass/detector/monitor/References/reference.root";
std::string Group::forcedRefFileName = "";
std::string Group::wkEnvStr = "AFS";
bool Group::fExpertHistos = false;
bool Group::fDoTracking = false;

Group::~Group() {
  fHistList.clear();
}

void Group::ResetHistograms() {

  fRateCounter = 0;
  if (thr_flag) TThread::Lock();
  for(unsigned int i=0; i<fHistList.size(); i++)
    fHistList[i]->Reset();
  if (thr_flag) TThread::UnLock();
}

void Group::ResetIfNeeded() {

  unsigned i=0;
  for(; i<fPlanes.size(); i++) {
    if(fPlanes[i]->IsModified())
      break;
  }
  if(i<fPlanes.size()) ResetHistograms();
}

void Group::Add(const Plane* plane) {
  if(plane) {
    bool isin=false;
    for(size_t i=0; i<fPlanes.size(); i++) {
      if((*fPlanes[i])==(*plane)) {
	isin=true;
	break;
      }
    }
    if(isin)
      std::cout<<"Plane "<<plane->GetName()<<" already in group "<<fName<<std::endl;
    else {
//       std::cout<<"Plane "<<plane->GetName()<<" added to group "<<fName<<std::endl;
      fPlanes.push_back(plane);
    }
  }
  else std::cerr<<"Group "<<fName<<": Plane does not exist !"<<std::endl;
}

const char* Group::GetType() const {

  register const char *grpname = GetName();
  register const char *grptype = fType.c_str();
  if (fType == "") grptype = grpname;

// cerr<<"Group::GetType: grpname "<<grpname<<" grptype "<<grptype<<endl;

  if(strncmp(grptype, "MM",2)==0) {
    return "MM";
  } else if(strncmp(grptype, "FI",2)==0) {
    return "FJ";
  } else if(strncmp(grptype, "DAQ",3)==0) {
    return "DAQ";
  } else if(strncmp(grptype, "GM",2)==0) {
    return "GM";
  } else if(strncmp(grptype, "GEM",3)==0) {
    return "GM";
  } else if(strncmp(grptype, "GP",2)==0) {
    return "GP";
  } else if(strncmp(grptype, "PGEM",4)==0) {
    return "GP";
  } else if(strncmp(grptype, "SIL",3)==0) {
    return "SI";
  } else if(strncmp(grptype, "PS",2)==0) {
    return "MWPC";
  } else if(strncmp(grptype, "PA",2)==0) {
    return "MWPC";
  } else if(strncmp(grptype, "PB",2)==0) {
    return "MWPC";
  } else if(strncmp(grptype, "SF",2)==0) {
    return "FG";
  } else if(strncmp(grptype, "HK",2)==0) {
    return "HK";
  } else if(strncmp(grptype, "H",1)==0) {
    return "HS";
  } else if(strncmp(grptype, "MWA",3)==0) {
    return "MA";
  } else if(strncmp(grptype, "ST",2)==0) {
    return "Straw";
  } else if(strncmp(grptype, "CEDAR",5)==0) {
    return "CEDAR";
  } else if(strncmp(grptype, "RA",2)==0) {
    return "RichAPV";
  } else if(strncmp(grptype, "RM",2)==0) {
    return "RichMAPMT";
  } else if(strncmp(grptype, "DR",2)==0) {
    return "DR";
  } else if(strncmp(grptype, "ALL",3)==0) {
    return "UNKNOWN";
  } else if(strncmp(grptype, "DC",2)==0) {
    return "DC";
  } else if(strncmp(grptype, "RPD",3)==0) {
    return "RPD";
  } else if(strncmp(grptype, "SW",2)==0) {
    return "SW";
  } else if(strncmp(grptype, "EC",2)==0) {
    return "EC";
  } else if(strncmp(grptype, "DW",2)==0) {
    return "W45";
  } else if(strncmp(grptype, "CA",2)==0) {
		return "Camera";
  }  

cerr << "Group::GetType: entering in old part of the method, strange... (grpname "
     << grpname << "  grptype "<<grptype<<")" << endl;


  if( NULL!=dynamic_cast<const GroupMumega*>(this) )
    return "MM";

  if( NULL!=dynamic_cast<const GroupScifiJ*>(this) )
    return "FJ";

  if( NULL!=dynamic_cast<const GroupDAQ*>(this) )
    return "DAQ";

  if( NULL!=dynamic_cast<const GroupGEM*>(this) )
    return "GEM";

  if( NULL!=dynamic_cast<const GroupGEMTrack*>(this) )
    return "GEMTrack";

  if( NULL!=dynamic_cast<const GroupPGEM*>(this) )
    return "GP";

  if( NULL!=dynamic_cast<const GroupSili*>(this) )
    return "Silicon";

  if( NULL!=dynamic_cast<const GroupSiliTrack*>(this) )
    return "SiliconTrack";

  if( NULL!=dynamic_cast<const GroupMwpc*>(this) )
    return "MWPC";

  if( NULL!=dynamic_cast<const GroupScifiG*>(this) )
    return "FG";

  if( NULL!=dynamic_cast<const GroupTrigHodo*>(this) )
    return "HS";

  if( NULL!=dynamic_cast<const GroupMuonWallA*>(this) )
    return "MW1";

  if( NULL!=dynamic_cast<const GroupStraws*>(this) )
    return "Straw";

  if( NULL!=dynamic_cast<const GroupBeamKiller*>(this) )
    return "HK";
  if( NULL!=dynamic_cast<const GroupCEDAR*>(this) )
    return "CEDAR";
  if( NULL!=dynamic_cast<const GroupRiAPV*>(this) )
    return "RichAPV";
  if( NULL!=dynamic_cast<const GroupRICH_MAPMT*>(this) )
    return "RichMAPMT";
  if( NULL!=dynamic_cast<const GroupRichWall*>(this) )
    return "RichWall";
  if( NULL!=dynamic_cast<const GroupRPD*>(this) )
    return "RPD";
  if( NULL!=dynamic_cast<const GroupCamera*>(this) )
    return "CA";
  return "UNKNOWN";
}


void Group::AddHistogram(TH1* hist) {
  fHistList.push_back(hist);
}


TH1* Group::GetHisto( const std::string hName ) {
  for( unsigned int i = 0; i< fHistList.size(); i++ )
  if( std::string(fHistList[i]->GetName()) == hName ||
      std::string(fHistList[i]->GetName()) == fName+hName ) return fHistList[i];

//  std::cout << "Plane::GetHisto - ERROR: cannot find histogram \"" << hName << "\".\n";
  return 0;
}


#if USE_DATABASE == 1
const char* Group::GetRefDirNameDB() {
  std::string letype(GetType());
  if(fDataBase) {
    return fDataBase->getRootRefFile(letype.c_str(), wkEnvStr.c_str());
  }
  else return 0;
}
#endif


const char* Group::GetRefDirName() {
  if (forcedRefFileName != "") return forcedRefFileName.c_str();
#if USE_DATABASE == 1
  const char* ref = GetRefDirNameDB();
  if(ref) {
    return ref;
  }
#endif
  return defaultRefFileName.c_str();
}


void Group::OpenReference() {
//  fReferenceDirectory=defaultRefdir->LoadReferenceFile(fRefDirName.c_str());
  if (TH1F_Ref::GetShowRef())
    fReferenceDirectory=ReferenceDirectory::GiveReferenceDir(GetRefDirName())->GiveReferenceFile();
}


void Group::ControlPanel(const TGWindow *p, const TGWindow *main) {

// std::cerr<<"Group::ControlPanel: original group panel"<<std::endl;
//  if (!fControlPanel) fControlPanel = new GroupPanel(p, main, 100, 100, this);
}

