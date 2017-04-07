
#include <iostream>

#include "Plane.h"
#include "GeomPlaneMumega.h"
#include "GeomPlanePMumega.h"
#include "GeomPlaneMwpc.h"
#include "GeomPlaneGEM.h"
#include "GeomPlanePGEM.h"
#include "GeomPlaneSili.h"
#include "GeomPlaneScifiJ.h"
#include "GeomPlaneMuonWallA.h"
#include "GeomPlaneScifiG.h"
#include "PlaneMwpc.h"
#include "PlaneHCALT.h"
#include "PlaneRICH.h"
#include "PlaneRICH_MAPMT.h"
#include "PlaneRiAPV.h"
#include "PlaneDriftChamber.h"
#include "PlaneMumega.h"
#include "PlanePMumega.h"
#include "PlaneScaler.h"
#include "PlaneDriftChamberPanel.h"
#include "PlaneMuonWallA.h"
#include "PlaneMuonWallB.h"
#include "PlaneScifiG.h"
#include "PlaneGEM.h"
#include "PlanePGEM.h"
#include "PlaneMwpc.h"
#include "PlaneScifiJ.h"
#include "PlaneHCAL1.h"
#include "PlaneECAL1.h"
#include "PlaneFEM.h"
#include "PlaneECAL2.h"
#include "PlaneSili.h"
#include "PlaneHCAL2.h"
#include "PlaneStrawTubes.h"
#include "PlaneBMS.h"
#include "PlaneTrigHodo.h"
#include "PlaneVeto.h"
#include "PlaneRW.h"
#include "PlaneRPD_SADC.h"
#include "PlaneRPD_F1.h"
#include "PlaneTrigger_SADC_F1.h"


ClassImp(Plane);

using namespace std;

bool Plane::sModified = false;
bool Plane::fExpertHistos = false;
bool Plane::fDoTracking = false;
std::string Plane::defaultRefFileName = "/afs/cern.ch/compass/detector/monitor/References/reference.root";
std::string Plane::forcedRefFileName = "";
std::string Plane::wkEnvStr = "AFS";

std::map<std::string,const Plane*> Plane::sInstances;


Plane::Plane(const char* name, TTree *tree)
  : TObject(), fReferenceDirectory(0), fRateCounter(0),
  fModified(false), fNeedGeom(false),fGeom(0),
  fName(name),fIsInTree(false),
  fNhits(0), fControlPanel(0), fTMPattern(0), fUseCalib(false),
  fTree(tree) {

  // store a pointer to this instance in static map for access
  // from other planes.

  Plane::sInstances[fName] = this;
}


Plane::~Plane() {

  std::cerr<<"deleting plane "<<fName<<std::endl;

  if (thr_flag) TThread::Lock();
  for(register unsigned i=0; i<fHistList.size(); i++)
    delete fHistList[i];
  if (thr_flag) TThread::UnLock();

  fHistList.clear();

  for(register unsigned i=0; i<fVariables.size(); i++)
    delete fVariables[i];
  fVariables.clear();
}


void Plane::Reset() {

  for(register unsigned int i=0;i<fVariables.size();i++) {
    fVariables[i]->Reset();
  }
  fNhitsKept=0;
  fNhits=0;

  lDigits.clear();
}


void Plane::ResetHistograms() {

  fRateCounter = 0;

  if (thr_flag) TThread::Lock();
  for(register unsigned int i=0; i<fHistList.size(); i++) {
    fHistList[i]->Reset();
  }
  if (thr_flag) TThread::UnLock();
}


bool Plane::operator== (const Plane& plane) const  {
  return fName==std::string(plane.GetName());
}


Variable* Plane::AddVariable(const char* name, int nbins, float min, float max,
			     int maxsize)  {
  fVariables.push_back(new Variable(name,nbins, min, max, maxsize));
  return fVariables.back();
}

void Plane::AddVariable(Variable* var){
  fVariables.push_back(var);
}


void Plane::AddHistogram(TH1* hist) {
  fHistList.push_back(hist);
}

//__________________________________________________________________
Variable* Plane::GetVariable( const std::string vName )
{
  for( unsigned int i = 0; i< fVariables.size(); i++ )
  if( fVariables[i]->GetName() == vName ) return fVariables[i];

  std::cout << "Plane::GeVariable - ERROR: cannot find variable \"" << vName << "\".\n";
  return 0;
}
//__________________________________________________________________
TH1* Plane::GetHisto( const std::string hName )
{
  for( unsigned int i = 0; i< fHistList.size(); i++ )
  if( std::string(fHistList[i]->GetName()) == hName ||
      std::string(fHistList[i]->GetName()) == fName+hName ) return fHistList[i];

//  std::cout << "Plane::GetHisto - ERROR: cannot find histogram \"" << hName << "\".\n";
  return 0;
}

//__________________________________________________________________
void Plane::SetGeometry(GeomPlane* geom) {
  GeomPlaneMumega *mmgeom = dynamic_cast<GeomPlaneMumega*> (geom);
  if(mmgeom) {
    PlaneMumega *mm = dynamic_cast<PlaneMumega*> (this);
    if(mm) {
      fGeom = geom;
      return;
    }
    else
      std::cerr<<"WARNING Plane::SetGeometry. PlaneMumega must be associated to GeomPlaneMumega !"<<std::endl;
  }
  GeomPlanePMumega *pmmgeom = dynamic_cast<GeomPlanePMumega*> (geom);
  if(pmmgeom) {
    PlanePMumega *pmm = dynamic_cast<PlanePMumega*> (this);
    if(pmm) {
      fGeom = geom;
      return;
    }
    else
      std::cerr<<"WARNING Plane::SetGeometry. PlanePMumega must be associated to GeomPlaneMumega !"<<std::endl;
  }
  GeomPlaneMwpc *mwpcgeom = dynamic_cast<GeomPlaneMwpc*> (geom);
  if(mwpcgeom) {
    PlaneMwpc *mwpc = dynamic_cast<PlaneMwpc*> (this);
    if(mwpc) {
      fGeom = geom;
      return;
    }
    else
      std::cerr<<"WARNING Plane::SetGeometry. PlaneMwpc must be associated to GeomPlaneMwpc !"<<std::endl;
  }
  GeomPlaneGEM *gemgeom = dynamic_cast<GeomPlaneGEM*> (geom);
  if(gemgeom) {
    PlaneGEM *gem = dynamic_cast<PlaneGEM*> (this);
    if(gem) {
      fGeom = geom;
      return;
    }
    else
      std::cerr<<"WARNING Plane::SetGeometry. PlaneGEM must be associated to GeomPlaneGEM !"<<std::endl;
  }
  GeomPlanePGEM *pgemgeom = dynamic_cast<GeomPlanePGEM*> (geom);
  if(pgemgeom) {
    PlanePGEM *pgem = dynamic_cast<PlanePGEM*> (this);
    if(pgem) {
      fGeom = geom;
      return;
    }
    else
      std::cerr<<"WARNING Plane::SetGeometry. PlanePGEM must be associated to GeomPlanePGEM !"<<std::endl;
  }
  GeomPlaneSili *silgeom = dynamic_cast<GeomPlaneSili*> (geom);
  if(silgeom) {
    PlaneSili *sil = dynamic_cast<PlaneSili*> (this);
    if(sil) {
      fGeom = geom;
      return;
    }
    else
      std::cerr<<"WARNING Plane::SetGeometry. PlaneSili must be associated to GeomPlaneSili !"<<std::endl;
  }
  GeomPlaneScifiJ *scijgeom = dynamic_cast<GeomPlaneScifiJ*> (geom);
  if(scijgeom) {
    PlaneScifiJ *scij = dynamic_cast<PlaneScifiJ*> (this);
    if(scij) {
      fGeom = geom;
      return;
    }
    else
      std::cerr<<"WARNING Plane::SetGeometry. PlaneScifiJ must be associated to GeomPlaneScifiJ !"<<std::endl;
  }

  GeomPlaneMuonWallA *mwageom = dynamic_cast<GeomPlaneMuonWallA*> (geom);
  if(mwageom) {
    PlaneMuonWallA *mwa = dynamic_cast<PlaneMuonWallA*> (this);
    if(mwa) {
      fGeom = geom;
      return;
    }
    else
      std::cerr<<"WARNING Plane::SetGeometry. PlaneMuonWallA must be associated to GeomPlaneMuonWallA !"<<std::endl;
  }
  GeomPlaneScifiG *sciggeom = dynamic_cast<GeomPlaneScifiG*> (geom);
  if(sciggeom) {
    PlaneScifiG *scig = dynamic_cast<PlaneScifiG*> (this);
    if(scig) {
      fGeom = geom;
      return;
    }
    else
      std::cerr<<"WARNING Plane::SetGeometry. PlaneScifiG must be associated to GeomPlaneScifiG !"<<std::endl;
  }
}


ostream& operator<<(ostream& out, const std::set<int>& dataset) {

  // this function compresses a set of integers :
  // 1,2,3,10,11,12,15 is written 1-3 10-12 15
  // used for missing/noisy channels output

  if(!out) return out;

  typedef std::set<int>::const_iterator IM;
  int last = -2;
  bool sep=false;
  for(IM im=dataset.begin(); im!=dataset.end(); im++) {
    if(*im == last+1) {
      if(!sep) {
	out<<"-";
	sep=true;
      }
    }
    else {
      if(sep) {
	out<<last;
	sep=false;
      }
      out<<" "<<(*im);
    }
    last = *im;
  }
  if(sep) out<<last;
  return out;
}


const char* Plane::GetType() const {

  register const char *detname = GetName();

  if(strncmp(detname, "MM",2)==0) {
    return "MM";
  } else if(strncmp(detname, "MP",2)==0) {
    return "MP";
  } else if(strncmp(detname, "PS",2)==0) {
    return "MWPC";
  } else if(strncmp(detname, "PA",2)==0) {
    return "MWPC";
  } else if(strncmp(detname, "PB",2)==0) {
    return "MWPC";
  } else if(strncmp(detname, "FI",2)==0) {
    std::string s(detname);
    char nums[3];
    s.copy(nums,2,2);
    nums[2] = '\000';
    register int num=atoi(nums);
    switch (num) {
     case 1:
     case 2:
     case 3:
     case 4:
      return "FJ";
      break;
     case 5:
     case 6:
     case 7:
     case 8:
     case 15:
      return "FG";
      break;
     case 35:
      return "FJ";
      break;
     case 55:
      return "FG";
      break;
    }
  } else if(strncmp(detname, "GM",2)==0) {
    return "GM";
  } else if(strncmp(detname, "GP",2)==0) {
    return "GP";
  } else if(strncmp(detname, "DC",2)==0) {
    return "DC";
  } else if(strncmp(detname, "ST",2)==0) {
    return "ST";
  } else if(strncmp(detname, "BM",2)==0) {
    return "BM";
  } else if(strncmp(detname, "DW",2)==0) {
    return "DW";
  } else if(strncmp(detname, "MA",2)==0) {
    return "MA";
  } else if(strncmp(detname, "MB",2)==0) {
    return "MB";
  } else if(strncmp(detname, "HB",2)==0) {
    return "HS";
  } else if(strncmp(detname, "HF",2)==0) {
    return "HS";
  } else if(strncmp(detname, "HI",2)==0) {
    return "HS";
  } else if(strncmp(detname, "HM",2)==0) {
    return "HS";
  } else if(strncmp(detname, "HL",2)==0) {
    return "HS";
  } else if(strncmp(detname, "HO",2)==0) {
    return "HS";
  } else if(strncmp(detname, "HP",2)==0) {
    return "HS";
  } else if(strncmp(detname, "HQ",2)==0) {
    return "HS";
  } else if(strncmp(detname, "HK",2)==0) {
    return "HK";
  } else if(strncmp(detname, "SC",2)==0) {
    return "SC";
  } else if(strncmp(detname, "HC",2)==0) {
    return "HCAL";
  } else if(strncmp(detname, "SI",2)==0) {
    return "SI";
  } else if(strncmp(detname, "RI",2)==0) {
    return "RICH";
  } else if(strncmp(detname, "VT",2)==0) {
    return "VETO";
  } else if(strncmp(detname, "RA",2)==0) {
    return "RichAPV";
  } else if(strncmp(detname, "DR",2)==0) {
    return "DR";
  } else if(strncmp(detname, "RM",2)==0) {
    return "RichMAPMT";
  } else if(strncmp(detname, "EC",2)==0) {
    return "ECAL";
  } else if(strncmp(detname, "TT",2)==0) {
    return "TT";
  } else if(strncmp(detname, "CE",2)==0) {
    return "CEDAR";
  } else if(strncmp(detname, "VB",2)==0) {
    return "VB";
  } else if(strncmp(detname, "DC",2)==0) {
    return "DC";
  } else if(strncmp(detname, "RP",2)==0) {
    return "RPD";
  } else if(strncmp(detname, "SW",2)==0) {
    return "SW";
  } else if(strncmp(detname, "TCSphase", 8) == 0) {
    return "TCSphase";
  }

cerr << "Plane::GetType: entering in old part of the method, strange... (detname "
     << detname << ")"<<endl;

  const PlaneMumega *mm = dynamic_cast<const PlaneMumega*>(this);
  const PlaneScifiJ *fj = dynamic_cast<const PlaneScifiJ*>(this);
  const PlaneGEM    *gm = dynamic_cast<const PlaneGEM*>(this);
  const PlaneSili   *si = dynamic_cast<const PlaneSili*>(this);
  const PlaneBMS    *bm = dynamic_cast<const PlaneBMS*>(this);
  const PlaneMwpc   *mwpc = dynamic_cast<const PlaneMwpc*>(this);
  const PlaneDriftChamber *dc = dynamic_cast<const PlaneDriftChamber*>(this);
  const PlaneMuonWallA *mwa = dynamic_cast<const PlaneMuonWallA*>(this);
  const PlaneMuonWallB *mwb = dynamic_cast<const PlaneMuonWallB*>(this);
  const PlaneStrawTubes *st = dynamic_cast<const PlaneStrawTubes*>(this);
  const PlaneScifiG *fg = dynamic_cast<const PlaneScifiG*>(this);
  const PlaneTrigHodo *hs = dynamic_cast<const PlaneTrigHodo*>(this);
  const PlaneECAL1 *ecal1 = dynamic_cast<const PlaneECAL1*>(this);
  const PlaneFEM *fem = dynamic_cast<const PlaneFEM*>(this);
  const PlaneECAL2 *ecal2 = dynamic_cast<const PlaneECAL2*>(this);
  const PlaneHCAL1 *hcal1 = dynamic_cast<const PlaneHCAL1*>(this);
  const PlaneHCAL2 *hcal2 = dynamic_cast<const PlaneHCAL2*>(this);
  const PlaneHCALT *hcalt = dynamic_cast<const PlaneHCALT*>(this);
  const PlaneRICH *rich = dynamic_cast<const PlaneRICH*>(this);
  const PlaneRICH_MAPMT *rich_mapmt = dynamic_cast<const PlaneRICH_MAPMT*>(this);
  const PlaneRiAPV *riapv = dynamic_cast<const PlaneRiAPV*>(this);
  const PlaneRCath *rcath = dynamic_cast<const PlaneRCath*>(this);
  const PlaneScaler *sc = dynamic_cast<const PlaneScaler*>(this);
  const PlaneVeto *veto = dynamic_cast<const PlaneVeto*>(this);
  const PlaneRichWall *richwall = dynamic_cast<const PlaneRichWall*>(this);
  const Plane1V *plane1V = dynamic_cast<const Plane1V*>(this);

  if(mm)
    return "MM";
  else if(fj)
    return "FJ";
  else if(gm)
    return "GM";
  else if(si)
    return "SI";
  else if(bm)
    return "BM";
  else if(mwpc)
    return "MWPC";
  else if(dc)
    return "DC";
  else if(mwa)
    return "MA";
  else if(mwb)
    return "MB";
  else if(st)
    return "ST";
  else if(fg)
    return "FG";
  else if(hs)
    return "HS";
  else if(ecal1)
    return "ECAL";
  else if(fem)
    return "ECAL";
  else if(ecal2)
    return "ECAL";
  else if(hcal1)
    return "HCAL";
  else if(hcal2)
    return "HCAL";
  else if(hcalt)
    return "HCAL";
  else if(rich)
    return "RICH";
  else if(rich_mapmt)
    return "RM";
  else if(rcath)
    return "RICH";
  else if(riapv)
    return "RichAPV";
  else if(sc)
    return "SC";
  else if(veto)
    return "VETO";
  else if(richwall)
    return "DR";
  else if(plane1V)
    return "PLANE1V";
  else
    return "UNKNOWN";
}


#if USE_DATABASE == 1
const char* Plane::GetRefDirNameDB() {
  std::string letype(GetType());
  if(letype == "UNKNOWN") return 0;
  if(fDataBase) {
    return fDataBase->getRootRefFile(letype.c_str(), wkEnvStr.c_str());
  }
  else return 0;
}
#endif


const char* Plane::GetRefDirName() {
  if (forcedRefFileName != "") return forcedRefFileName.c_str();
#if USE_DATABASE == 1
  const char* ref = GetRefDirNameDB();
  if(ref) {
    return ref;
  }
#endif
  return defaultRefFileName.c_str();
}


void Plane::OpenReference() {
//  fReferenceDirectory=defaultRefdir->LoadReferenceFile(fRefDirName.c_str());
  if (TH1F_Ref::GetShowRef())
    fReferenceDirectory=ReferenceDirectory::GiveReferenceDir(GetRefDirName())->GiveReferenceFile();
}


const Plane* Plane::GetPlane(const char* name) const {
  std::map<std::string, const Plane*>::iterator iplane = Plane::sInstances.find(name);
  if(iplane != Plane::sInstances.end() ) return iplane->second;
  else return 0;
}
