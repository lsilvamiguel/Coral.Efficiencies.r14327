#include <memory.h>

#include <CsBeam.h>
#include <Reco/CalorimeterParticle.h>
#include "../src/evmc/pdgParticles.h"

#include "CsDst_04.h"
#include "runflag.h"

using namespace std;
using namespace CLHEP;

CsDstCluster_04::CsDstCluster_04()
{
  type = CLUSTER_TYPE_04;
  dstVersion = 4;
  u = v = w = time = analog = 0.;
  memset(cov,0,sizeof(cov));
  hasTime = hasAnalog = 0;
  detectorPosition = 0;
}

void CsDstCluster_04::Init(const CsCluster& cl)
{
  u                = cl.getU();
  v                = cl.getV();
  w                = cl.getW();
  double t = 0.;
  hasTime          = (cl.getTime(t)?1:0);
  time             = (float32)t;
  hasAnalog        = (cl.getAnalog(t)?1:0);
  analog           = (float32)t;

  HepMatrix covm = cl.getCov();
  cov[0]       = covm(1,1);
  cov[1]       = covm(1,2);
  cov[2]       = covm(2,2);
  cov[3]       = covm(1,3);
  cov[4]       = covm(2,3);
  cov[5]       = covm(3,3);

  map<CsDetector*,int,less<CsDetector*> > detmap = 
    CsEvent::Instance()->getDetectorMap();
  if( ! detmap.empty() )  // protection
    detectorPosition = detmap[ cl.getDetsList().front() ];
  else 
    detectorPosition = -1;
}

CsCluster* CsDstCluster_04::GetCluster() const
{
  HepMatrix covm(3,3,0); 
  covm(1,1) = cov[0];
  covm(1,2) = cov[1];
  covm(2,2) = cov[2];
  covm(1,3) = cov[3];
  covm(2,3) = cov[4];
  covm(3,3) = cov[5];
  CsCluster* cluster = new CsCluster(u,v,w,covm);
  if(hasTime) cluster->setTime(double(time));
  if(hasAnalog) cluster->setAnalog(double(analog));
  vector<CsDetector*> detvec = CsEvent::Instance()->getDetectorVector();
  if(detectorPosition != -1) cluster->addDet(* detvec[detectorPosition]);
  return cluster;
}

bool CsDstCluster_04::Load(std::istream& file)
{
  file.read((char*)&u,sizeof(u));
  file.read((char*)&v,sizeof(v));
  file.read((char*)&w,sizeof(w));
  file.read((char*)&detectorPosition,sizeof(detectorPosition));
  file.read((char*)&hasTime,sizeof(hasTime));
  if(hasTime) file.read((char*)&time,sizeof(time));
  file.read((char*)&hasAnalog,sizeof(hasAnalog));
  if(hasAnalog) file.read((char*)&analog,sizeof(analog));
  if(file.fail()) 
    {
      cerr << "CsDstCluster_04::Load(istream): Input stream reading error." << endl;
      return false;
    }
  return true;
}

bool CsDstCluster_04::Load(FlatFile& file)
{
  if(!file.Read((char*)&u,sizeof(u))) return false;
  if(!file.Read((char*)&v,sizeof(v))) return false;
  if(!file.Read((char*)&w,sizeof(w))) return false;
  if(!file.Read((char*)&detectorPosition,sizeof(detectorPosition))) return false;
  if(!file.Read((char*)&hasTime,sizeof(hasTime))) return false;
  if(hasTime) if(!file.Read((char*)&time,sizeof(time))) return false;
  if(!file.Read((char*)&hasAnalog,sizeof(hasAnalog))) return false;
  if(hasAnalog) if(!file.Read((char*)&analog,sizeof(analog))) return false;
  return true;
}

bool CsDstCluster_04::IsEquil(const CsDstObject* obj) const
{
  if(type != obj->type) return false;
  const CsDstCluster_04* cl = (const CsDstCluster_04*)obj;
  if(u != cl->u) return false;
  if(v != cl->v) return false;
  if(w != cl->w) return false;
  if(hasTime != cl->hasTime) return false;
  if(hasTime && time != cl->time) return false;
  if(hasAnalog != cl->hasAnalog) return false;
  if(hasAnalog && analog != cl->analog) return false;
  if(detectorPosition != cl->detectorPosition) return false;
  for(int i = 0; i < 6; i++)
    if(cov[i] != cl->cov[i]) return false;
  return true;  
}

void CsDstCluster_04::Print(std::ostream& msgout) const
{
  char MSG_STR[1024];
  sprintf(MSG_STR,"\t(U,V,W)=(%.1f, %.1f, %.1f); Det.position: %d;", 
	      u, v, w, detectorPosition);
  msgout << MSG_STR;
  if(hasAnalog) msgout << " Analog: " << analog;
  if(hasTime) msgout << " Time: " << time;
  msgout << "\n";
}

void CsDstCluster_04::GetBuffer(std::vector<uint8>& buffer) const
{
  AddValueToVector(buffer,&u,sizeof(u));
  AddValueToVector(buffer,&v,sizeof(v));
  AddValueToVector(buffer,&w,sizeof(w));
  AddValueToVector(buffer,&detectorPosition,sizeof(detectorPosition));
  AddValueToVector(buffer,&hasTime,sizeof(hasTime));
  if(hasTime) AddValueToVector(buffer,&time,sizeof(time));
  AddValueToVector(buffer,&hasAnalog,sizeof(hasAnalog));
  if(hasAnalog) AddValueToVector(buffer,&analog,sizeof(analog));
}

CsDstHelix_04::CsDstHelix_04()
{
  type = HELIX_TYPE_04;
  dstVersion = 4;
  x = y = z = dxdz = dydz = cop = 0.;
  memset(cov,0,sizeof(cov));
}

void CsDstHelix_04::Init(const CsHelix& helix)
{
  CsHelix& hl = const_cast<CsHelix&>(helix);
  x    = hl.getX();
  y    = hl.getY();
  z    = hl.getZ();
  dxdz = hl.getDXDZ();
  dydz = hl.getDYDZ();
  cop  = hl.getCop();
  for(int i = 0; i < 15; i++)
    cov[i] = hl.getCov()[i];
}

CsHelix* CsDstHelix_04::GetHelix() const
{
  double covm[15];
  for(int i = 0; i < 15; i++)
    covm[i] = cov[i];
  CsHelix* helix = new CsHelix(x, y, z, dxdz, dydz, cop, covm);
  return helix;
}

bool CsDstHelix_04::Load(std::istream& file)
{
  file.read((char*)&x,sizeof(x));
  file.read((char*)&y,sizeof(y));
  file.read((char*)&z,sizeof(z));
  file.read((char*)&dxdz,sizeof(dxdz));
  file.read((char*)&dydz,sizeof(dydz));
  file.read((char*)&cop,sizeof(cop));
  file.read((char*)&cov[0],sizeof(cov));
  if(file.fail())
    {
      cerr << "CssHelix04::Load(): Input stream reading error." << endl;
      return false;
    }
  return true;
}

bool CsDstHelix_04::Load(FlatFile& file)
{
  if(!file.Read((char*)&x,sizeof(x))) return false;
  if(!file.Read((char*)&y,sizeof(y))) return false;
  if(!file.Read((char*)&z,sizeof(z))) return false;
  if(!file.Read((char*)&dxdz,sizeof(dxdz))) return false;
  if(!file.Read((char*)&dydz,sizeof(dydz))) return false;
  if(!file.Read((char*)&cop,sizeof(cop))) return false;
  if(!file.Read((char*)&cov[0],sizeof(cov))) return false;
  return true;
}

bool CsDstHelix_04::IsEquil(const CsDstObject* obj) const
{
  if(type != obj->type) return false;
  const CsDstHelix_04* hl = (const CsDstHelix_04*)obj;
  if(x != hl->x) return false;
  if(y != hl->y) return false;
  if(z != hl->z) return false;
  if(dxdz != hl->dxdz)
    {
      if(fabs(1. - dxdz/hl->dxdz) > 0.1)
	{
	  cout << "CsDstHelix_04: very big dxdz difference: "
	       << dxdz << " and " << hl->dxdz << endl;
	}
    }
  if(dydz != hl->dydz)
    {
      if(fabs(1. - dydz/hl->dydz) > 0.1)
	{
	  cout << "CsDstHelix_04: very big dxdz difference: "
	       << dydz << " and " << hl->dydz << endl;
	}
    }
  if(cop != hl->cop) return false;
  for(int i = 0; i < 15; i++)
    if(cov[i] != hl->cov[i]) return false;
  return true;
}

void CsDstHelix_04::Print(std::ostream& msgout) const
{
  char MSG_STR[1024];
  sprintf(MSG_STR,"\t(X,Y,Z)=(%.1f, %.1f, %.1f); DX/DZ = %f; DY/DZ = %f; cop = %f;\n",
	      x, y, z, dxdz, dydz, cop);
  msgout << MSG_STR;
  msgout << "\t\tCov.Matrix:\n";
  for(int j = 0; j < 15; j++)
    {
      if(j%5 == 0) msgout << "\t\t";
      sprintf(MSG_STR,"%8.2f ",  cov[j]);
      msgout << MSG_STR;
      if(j%5 == 4) msgout << "\n";
    }
}

void CsDstHelix_04::GetBuffer(std::vector<uint8>& buffer) const
{
  AddValueToVector(buffer,&x,sizeof(x));
  AddValueToVector(buffer,&y,sizeof(y));
  AddValueToVector(buffer,&z,sizeof(z));
  AddValueToVector(buffer,&dxdz,sizeof(dxdz));
  AddValueToVector(buffer,&dydz,sizeof(dydz));
  AddValueToVector(buffer,&cop,sizeof(cop));
  AddValueToVector(buffer,cov,sizeof(cov));
}

CsDstTrack_04::CsDstTrack_04()
{
  type = TRACK_TYPE_04;
  dstVersion = 4;
  chi2 = XX0 = 0.;
  flags = 0;
  time = 9999.0;    // default values if the track
  timeError = -1.;  // has not time information
  memset(rich1Probs,0,sizeof(rich1Probs));
  memset(expDets,0,sizeof(expDets));
  memset(firDets,0,sizeof(firDets));
}

void CsDstTrack_04::Init(const CsTrack& track)
{
  CsTrack& tr = const_cast<CsTrack&>(track);

  chi2 = tr.getChi2();
  if(tr.hasMeanTime())
    {
      time = tr.getMeanTime();
      timeError = tr.getMeanTimeError();
    }
  XX0 = tr.getXX0();
  if( tr.hasRich1Probs() ) 
    {
      const double* _rich1Probs = tr.getRich1Probs();
      for(unsigned i=0; i<sizeof(rich1Probs)/sizeof(rich1Probs[0]); i++ ) 
	rich1Probs[i] = _rich1Probs[i];
    }
  memcpy(expDets, tr.getExpectedDetsBitmap(),sizeof(expDets));
  memcpy(firDets, tr.getFiredDetsBitmap(),sizeof(firDets));
  const vector<CsHelix> helices = tr.getHelices();
  double zmin = 1E6;
  if(!helices.empty())
    {
      for(uint32 i = 0; i < helices.size(); i++)
	{
	  double z = helices[i].getZ();
	  if(z < zmin)
	    zmin = z;
	}
    }
  if(zmin < 0.)
    {
      uint32 test = 0;
      for(int i = 0; i < CSTRACK_MAPSIZE; i++)
	test += firDets[i];
      if(test == 0 && chi2 == 99)
	{
	  flags = 1; // it's beam track
	}
    }
}

CsTrack* CsDstTrack_04::GetTrack() const
{
  CsTrack* track = NULL;
  if((flags&1) != 0)
    track = new CsBeam;
  else
    track = new CsTrack;

  track->setChi2( chi2 );

  track->setMeanTime( time, timeError );

  track->setXX0( XX0 );
  double r1p[sizeof(rich1Probs)/sizeof(rich1Probs[0])]; 
  memset(r1p,0,sizeof(r1p));
  bool hasRich1Prob = false;
  for(unsigned i = 0; i < sizeof(rich1Probs)/sizeof(rich1Probs[0]); i++)
    {
      r1p[i] = rich1Probs[i];
      if(r1p[i]) hasRich1Prob = true;
    }
  if(hasRich1Prob)
    track->setRich1Probs(r1p);

  track->setFiredDetsBitmap(firDets);
  track->setExpectedDetsBitmap(expDets);

  return track;
}

bool CsDstTrack_04::Load(std::istream& file)
{
  file.read((char*)&chi2, sizeof(chi2));
  file.read((char*)&time, sizeof(time));
  file.read((char*)&timeError, sizeof(timeError));
  file.read((char*)&XX0, sizeof(XX0));

  //VD061218
  //Before 2006 production dsts were not saved with full CSTRACK_RICHDATASIZE rich probs
  //
  memset(rich1Probs,0,sizeof(rich1Probs));  // reset probs because fewer might be read!!!
  if ( runGT46692 ) {
	  file.read((char*)&rich1Probs[0],sizeof(rich1Probs));
  } else {
	  file.read((char*)&rich1Probs[0],sizeof(rich1Probs[0])*15); //read only first 15 probs!!!
  }
  //VD061218
  file.read((char*)&flags,sizeof(flags));
  memset(expDets,0,sizeof(expDets));
  memset(firDets,0,sizeof(firDets));
  unsigned expDetsSize(sizeof(expDets));
  unsigned firDetsSize(sizeof(firDets));
  if( !runGT32832 ) { // xxxDets size was [10] before 2004 
	  expDetsSize = sizeof(expDets[0])*10;
	  firDetsSize = sizeof(firDets[0])*10;
  }
  file.read((char*)&expDets[0], expDetsSize);
  file.read((char*)&firDets[0], firDetsSize);
  if(file.fail())
    {
      cerr << "CsDstTrack_04::Load(): Input stream reading error." << endl;
      return false;
    }
  return true;
}

bool CsDstTrack_04::Load(FlatFile& file)
{
  if(!file.Read((char*)&chi2, sizeof(chi2))) return false;
  if(!file.Read((char*)&time, sizeof(time))) return false;
  if(!file.Read((char*)&timeError, sizeof(timeError))) return false;
  if(!file.Read((char*)&XX0, sizeof(XX0))) return false;
  //VD061218
  //Before 2006 production dsts were not saved with full CSTRACK_RICHDATASIZE rich probs
  //
  memset(rich1Probs,0,sizeof(rich1Probs));  // reset probs because fewer might be read!!!
  if ( runGT46692 ) {
	  if(!file.Read((char*)&rich1Probs[0],sizeof(rich1Probs))) return false;
  } else {
	  if(!file.Read((char*)&rich1Probs[0],sizeof(rich1Probs[0])*15)) return false; //read only first 15 probs!!!
  }
  //VD061218
  if(!file.Read((char*)&flags,sizeof(flags))) return false;
  memset(firDets,0,sizeof(firDets));
  unsigned expDetsSize(sizeof(expDets));
  unsigned firDetsSize(sizeof(firDets));
  if( !runGT32832 ) { // xxxDets size was [10] before 2004
	  expDetsSize = sizeof(expDets[0])*10;
	  firDetsSize = sizeof(firDets[0])*10;
  }
  if(!file.Read((char*)&expDets[0],expDetsSize)) return false;
  if(!file.Read((char*)&firDets[0],firDetsSize)) return false;
  return true;
}

bool CsDstTrack_04::IsEquil(const CsDstObject* obj) const
{
  if(type != obj->type) return false;
  const CsDstTrack_04* tr = (const CsDstTrack_04*)obj;
  if(chi2 != tr->chi2)
    {
      if(fabs(chi2 - tr->chi2) > 0.01)
	{
	  cout << "CsDstTrack_04: Very big chi2 difference: "
	       << chi2 << " and " << tr->chi2 << endl;
	}
    }
  if(time != tr->time) return false;
  if(timeError != tr->timeError) return false;
  if(XX0 != tr->XX0) return false;
  for(unsigned i = 0; i < sizeof(rich1Probs)/sizeof(rich1Probs[0]); i++)
    if(rich1Probs[i] != tr->rich1Probs[i]) return false;
  if(flags != tr->flags) return false;
  for(unsigned i = 0; i < sizeof(expDets)/sizeof(expDets[0]); i++)
    if(expDets[i] != tr->expDets[i]) return false;
  for(unsigned i = 0; i < sizeof(firDets)/sizeof(firDets[0]); i++)
    if(firDets[i] != tr->firDets[i]) return false;
  return true;
}

void CsDstTrack_04::Print(std::ostream& msgout) const
{
  char MSG_STR[1024];
  msgout << (((flags&1) == 1)?"Beam track:\n":"Track:\n");
  sprintf(MSG_STR,"\tchi2 = %.1g; time = %f; time.error = %f; n.rad.length = %f;\n",
	      chi2, time, timeError, XX0);
  msgout << MSG_STR;
  msgout << "\t\tRICH1 infos on this track:\n";
  for(unsigned j = 0; j < sizeof(rich1Probs)/sizeof(rich1Probs[0]); j++)
    {
      if(j%5 == 0) msgout << "\t\t";
      sprintf(MSG_STR,"%f ",  rich1Probs[j]);
      msgout << MSG_STR;
      if(j%5 == 4) msgout << "\n";
    }
  msgout << "\t\tBit map of expected fired dets:\n";
  for(unsigned j = 0; j < sizeof(expDets)/sizeof(expDets[0]); j++)
    {
      if(j%5 == 0) msgout << "\t\t";
      sprintf(MSG_STR,"%8X ",  expDets[j]);
      msgout << MSG_STR;
      if(j%5 == 4) msgout << "\n";
    }
  msgout << "\t\tBit map of fired fired dets:\n";
  for(unsigned j = 0; j < sizeof(firDets)/sizeof(firDets[0]); j++)
    {
      if(j%5 == 0) msgout << "\t\t";
      sprintf(MSG_STR,"%8X ",  firDets[j]);
      msgout << MSG_STR;
      if(j%5 == 4) msgout << "\n";
    }
}

void CsDstTrack_04::GetBuffer(std::vector<uint8>& buffer) const
{
  AddValueToVector(buffer,&chi2, sizeof(chi2));
  AddValueToVector(buffer,&time, sizeof(time));
  AddValueToVector(buffer,&timeError, sizeof(timeError));
  AddValueToVector(buffer,&XX0, sizeof(XX0));
  AddValueToVector(buffer,rich1Probs,sizeof(rich1Probs));
  AddValueToVector(buffer,&flags,sizeof(flags));
  unsigned expDetsSize(sizeof(expDets));
  unsigned firDetsSize(sizeof(firDets));
  if( runGT32832 ) { // xxxDets size was [10] before 2004 
	  expDetsSize = sizeof(expDets[0])*10;
	  firDetsSize = sizeof(firDets[0])*10;
  }
  AddValueToVector(buffer,expDets,expDetsSize);
  AddValueToVector(buffer,firDets,firDetsSize);
}


CsDstVertex_04::CsDstVertex_04()
{
  type = VERTEX_TYPE_04;
  dstVersion = 4;
  x = y = z = chi2 = 0.;
  isPrimary = 0;
}

void CsDstVertex_04::Init(const CsVertex& vt)
{
  x = vt.getX();
  y = vt.getY();
  z = vt.getZ();
  chi2 = vt.getChi2();
  isPrimary = (vt.isPrimary()?1:0);
}

CsVertex* CsDstVertex_04::GetVertex() const
{
  CsVertex* vertex = new CsVertex(x,y,z);
  vertex->setChi2(chi2);
  vertex->setType(isPrimary == 1);
  return vertex;  
}

bool CsDstVertex_04::Load(std::istream& file)
{
  file.read((char*)&x, sizeof(x));
  file.read((char*)&y, sizeof(y));
  file.read((char*)&z, sizeof(z));
  file.read((char*)&chi2, sizeof(chi2));
  file.read((char*)&isPrimary, sizeof(isPrimary));
  if(file.fail()) return false;
  return true;
}

bool CsDstVertex_04::Load(FlatFile& file)
{
  if(!file.Read((char*)&x, sizeof(x))) return false;
  if(!file.Read((char*)&y, sizeof(y))) return false;
  if(!file.Read((char*)&z, sizeof(z))) return false;
  if(!file.Read((char*)&chi2, sizeof(chi2))) return false;
  if(!file.Read((char*)&isPrimary, sizeof(isPrimary))) return false;
  return true;
}

bool CsDstVertex_04::IsEquil(const CsDstObject* obj) const
{
  if(type != obj->type) return false;
  const CsDstVertex_04* vt = (const CsDstVertex_04*)obj;
  if(chi2 != vt->chi2) return false;
  if(x != vt->x) return false;
  if(y != vt->y) return false;
  if(z != vt->z) return false;
  if(isPrimary != vt->isPrimary) return false;
  return true;
}

void CsDstVertex_04::Print(std::ostream& msgout) const
{
  char MSG_STR[256];
  sprintf(MSG_STR,"\tIsPrimary = \"%3s\"; (X,Y,Z)=(%.1f,%.1f,%.1f); chi2 = %f\n",
	      (isPrimary?"yes":"no"), x, y, z, chi2);
  msgout << MSG_STR;
}

void CsDstVertex_04::GetBuffer(std::vector<uint8>& buffer) const
{
  AddValueToVector(buffer, &x, sizeof(x));
  AddValueToVector(buffer, &y, sizeof(y));
  AddValueToVector(buffer, &z, sizeof(z));
  AddValueToVector(buffer, &chi2, sizeof(chi2));
  AddValueToVector(buffer, &isPrimary, sizeof(isPrimary));
}

CsDstCalobject_04::CsDstCalobject_04()
{
  type = CALOBJ_TYPE_04;
  dstVersion = 4;
  x = y = z = e = 0.;
  xErr = yErr = zErr = eErr = 0.;
}

void CsDstCalobject_04::Init(const Reco::CalorimeterParticle& co)
{
  e = co.GetE();
  x = co.GetX();
  y = co.GetY();
  z = co.GetZ();
  eErr = co.GetEerr();
  xErr = co.GetXerr();
  yErr = co.GetYerr();
  zErr = co.GetZerr();
}

Reco::CalorimeterParticle* CsDstCalobject_04::GetCalobject() const
{
  return new Reco::CalorimeterParticle( Reco::CalorimeterParticle::GAMMA, 
					0, e, x, y, z,
                                        NULL,
					eErr, xErr, yErr, zErr, 0, 0 );
}
 
bool CsDstCalobject_04::Load(std::istream& file)
{
  file.read((char*)&x, sizeof(x));
  file.read((char*)&y, sizeof(y));
  file.read((char*)&z, sizeof(z));
  file.read((char*)&e, sizeof(e));
  file.read((char*)&xErr, sizeof(xErr));
  file.read((char*)&yErr, sizeof(yErr));
  file.read((char*)&zErr, sizeof(zErr));
  file.read((char*)&eErr, sizeof(eErr));
  if(file.fail()) return false;
  return true;
}

bool CsDstCalobject_04::Load(FlatFile& file)
{
  if(!file.Read((char*)&x, sizeof(x))) return false;
  if(!file.Read((char*)&y, sizeof(y))) return false;
  if(!file.Read((char*)&z, sizeof(z))) return false;
  if(!file.Read((char*)&e, sizeof(e))) return false;
  if(!file.Read((char*)&xErr, sizeof(xErr))) return false;
  if(!file.Read((char*)&yErr, sizeof(yErr))) return false;
  if(!file.Read((char*)&zErr, sizeof(zErr))) return false;
  if(!file.Read((char*)&eErr, sizeof(eErr))) return false;
  return true;
}

bool CsDstCalobject_04::IsEquil(const CsDstObject* obj) const
{
  if(type != obj->type) return false;
  const CsDstCalobject_04* co = (const CsDstCalobject_04*)obj;
  if((finite(x) || finite(co->x)) && x != co->x) return false;
  if((finite(y) || finite(co->y)) && y != co->y) return false;
  if((finite(z) || finite(co->z)) && z != co->z) return false;
  if(e != co->e) return false;
  if(xErr != co->xErr) return false;
  if(yErr != co->yErr) return false;
  if(zErr != co->zErr) return false;
  if(eErr != co->eErr) return false;
  return true;
}

void CsDstCalobject_04::Print(std::ostream& msgout) const
{
  char MSG_STR[256];
  sprintf(MSG_STR,"\t(X,Y,Z)=(%.1f, %.1f, %.1f); (Xerr,Yerr,Zerr)=(%.1f, %.1f, %.1f);\n"
	      "\t                  E = %f, Eerr = %f\n",
	      x, y, z, xErr, yErr,zErr, e, eErr);
  msgout << MSG_STR;
}

void CsDstCalobject_04::GetBuffer(std::vector<uint8>& buffer) const
{
  AddValueToVector(buffer, &x, sizeof(x));
  AddValueToVector(buffer, &y, sizeof(y));
  AddValueToVector(buffer, &z, sizeof(z));
  AddValueToVector(buffer, &e, sizeof(e));
  AddValueToVector(buffer, &xErr, sizeof(xErr));
  AddValueToVector(buffer, &yErr, sizeof(yErr));
  AddValueToVector(buffer, &zErr, sizeof(zErr));
  AddValueToVector(buffer, &eErr, sizeof(eErr));
}

CsDstParticle_04::CsDstParticle_04()
{
  type = PARTICLE_TYPE_04;
  dstVersion = 4;
  partType = 0;
  PDGid = 0;
}

void CsDstParticle_04::Init(const CsParticle& pat)
{
  CsParticle& pt = const_cast<CsParticle&>(pat);
  partType = pt.getType();
  string name = pt.getName();
  PDGid = 0;
  for( int i=0; i<nParticles; i++ ) 
    {
      if( name == PDGpart[i].name ) 
	{
	  PDGid  = PDGpart[i].number;
	}
    }
}

CsParticle* CsDstParticle_04::GetParticle() const
{
  CsParticle* particle = new CsParticle();
  particle->setType(partType);
  particle->setName("");
  for( int i = 0; i < nParticles; i++ ) 
    {
      if( PDGid == PDGpart[i].number ) 
	{
	  particle->setName( PDGpart[i].name );
	}
    }
  return particle;
}

bool CsDstParticle_04::Load(std::istream& file)
{
  file.read((char*)&partType, sizeof(partType));
  file.read((char*)&PDGid, sizeof(PDGid));
  if(file.fail()) return false;
  return true;
}

bool CsDstParticle_04::Load(FlatFile& file)
{
  if(!file.Read((char*)&partType, sizeof(partType))) return false;
  if(!file.Read((char*)&PDGid, sizeof(PDGid))) return false;
  return true;
}

bool CsDstParticle_04::IsEquil(const CsDstObject* obj) const
{
  if(type != obj->type) return false;
  const CsDstParticle_04* pt = (const CsDstParticle_04*)obj;
  if(partType != pt->partType) return false;
  if(PDGid != pt->PDGid) return false;
  return true;  
}

void CsDstParticle_04::Print(std::ostream& msgout) const
{
  char MSG_STR[256];
  sprintf(MSG_STR,"\tParticle type = %d; PDG_id: %d (%s)\n",
	      partType, PDGid, GetPartName().c_str());
  msgout << MSG_STR;
}

void CsDstParticle_04::GetBuffer(std::vector<uint8>& buffer) const
{
  AddValueToVector(buffer, &partType, sizeof(partType));
  AddValueToVector(buffer, &PDGid, sizeof(PDGid));
}

std::string CsDstParticle_04::GetPartName() const
{
  string name = "";
  for( int i = 0; i < nParticles; i++ ) 
    {
      if( PDGid == PDGpart[i].number ) 
	{
	  name = string( PDGpart[i].name );
	}
    }
  return name;
}

void CsDstEventHeader_04::Init(const CsRecoEvent& event)
{
  runNumber    = event.getEventHeader()[1];
  eventNumber  = event.getEventHeader()[2]; // see CsEvent.cc:1359
  eventInBurst = event.getEventHeader()[4];
  burstNumber  = event.getEventHeader()[3];
  timeInSec    = event.getEventHeader()[6];
  timeInUSec   = event.getEventHeader()[7];
  errorCode    = event.getEventHeader()[8];
  triggerMask  = event.getEventHeader()[9];
}

void CsDstEventHeader_04::GetHeader(CsRecoEvent& event) const
{
  event.setEventHeader(runNumber,1);
  event.setEventHeader(eventNumber,2);
  event.setEventHeader(eventInBurst,4);
  event.setEventHeader(burstNumber,3);
  event.setEventHeader(timeInSec,6);
  event.setEventHeader(timeInUSec,7);
  event.setEventHeader(errorCode,8);
  event.setEventHeader(triggerMask,9);
}

bool CsDstEventHeader_04::Load(std::istream& file)
{
  file.read((char*)&eventNumber,   sizeof(eventNumber));
  file.read((char*)&eventInBurst,  sizeof(eventInBurst));
  file.read((char*)&burstNumber,   sizeof(burstNumber));
  file.read((char*)&timeInSec,     sizeof(timeInSec));
  file.read((char*)&timeInUSec,    sizeof(timeInUSec));
  file.read((char*)&errorCode,     sizeof(errorCode));
  file.read((char*)&triggerMask,   sizeof(triggerMask));
  file.read((char*)&extVars[0],    sizeof(extVars));
  if(file.fail()) return false;
  return true;  
}

bool CsDstEventHeader_04::Load(FlatFile& file)
{
  if(!file.Read((char*)&eventNumber,   sizeof(eventNumber)))  return false;
  if(!file.Read((char*)&eventInBurst,  sizeof(eventInBurst))) return false;
  if(!file.Read((char*)&burstNumber,   sizeof(burstNumber)))  return false;
  if(!file.Read((char*)&timeInSec,     sizeof(timeInSec)))    return false;
  if(!file.Read((char*)&timeInUSec,    sizeof(timeInUSec)))   return false;
  if(!file.Read((char*)&errorCode,     sizeof(errorCode)))    return false;
  if(!file.Read((char*)&triggerMask,   sizeof(triggerMask)))  return false;
  if(!file.Read((char*)&extVars[0],    sizeof(extVars)))      return false;
  return true;  
}

bool CsDstEventHeader_04::IsEquil(const CsDstObject* obj) const
{
  if(type != obj->type) return false;
  const CsDstEventHeader_04* header = (const CsDstEventHeader_04*)obj;
  if(eventNumber != header->eventNumber) return false;
  if(eventInBurst != header->eventInBurst) return false;
  if(burstNumber != header->burstNumber) return false;
  if(timeInSec != header->timeInSec) return false;
  if(timeInUSec != header->timeInUSec) return false;
  if(errorCode != header->errorCode) return false;
  if(triggerMask != header->triggerMask) return false;
  for(int i = 0; i < 3; i++)
    if(extVars[i] != header->extVars[i]) return false;
  return true;  
}

void CsDstEventHeader_04::Print(std::ostream& msgout) const
{
  msgout << "\tEvent number:             " << eventNumber                   << "\n";
  msgout << "\tNumber of event in burst: " << eventInBurst                  << "\n";
  msgout << "\tBurst number:             " << burstNumber                   << "\n";
  msgout << "\tEvent time:               " << CsTime(timeInSec, timeInUSec) << "\n";
  msgout << "\tError code:               " << errorCode                     << "\n";
  msgout << "\tTrigger mask:             " << triggerMask                   << "\n";
  msgout << "\tExternal variable [0]:    " << extVars[0]                    << "\n";
  msgout << "\tExternal variable [1]:    " << extVars[1]                    << "\n";
  msgout << "\tExternal variable [2]:    " << extVars[2]                    << "\n";
}

void CsDstEventHeader_04::GetBuffer(std::vector<uint8>& buffer) const
{
  AddValueToVector(buffer, &eventNumber,   sizeof(eventNumber));
  AddValueToVector(buffer, &eventInBurst,  sizeof(eventInBurst));
  AddValueToVector(buffer, &burstNumber,   sizeof(burstNumber));
  AddValueToVector(buffer, &timeInSec,     sizeof(timeInSec));
  AddValueToVector(buffer, &timeInUSec,    sizeof(timeInUSec));
  AddValueToVector(buffer, &errorCode,     sizeof(errorCode));
  AddValueToVector(buffer, &triggerMask,   sizeof(triggerMask));
  AddValueToVector(buffer, extVars,        sizeof(extVars));
}

CsDstEvent_04::CsDstEvent_04(int fat)
{
  type = DST_EVENT_TYPE_04;
  dstVersion = 4;
  CsDstEvent::extVars.resize(3);
  //  nReserved = 1; // Only TimeInSpill information reserved at present.
  fatness = fat;
  vta.clear();
  tca.clear();
  ptca.clear();
  triggerTime = 0;
  nVtxToTrk = nTrkToClu = nPartTkCo = nClusters = nHelices =
    nTracks = nVertices = nCalobjects = nParticles = nTrkParAtVtx =
    nVtxErrMatrices = nHodoData = 0;
  memset(scalers,0,sizeof(scalers));
  vtxToTrk.clear();
  trkToClu.clear();
  partTkCo.clear();
  hodoData.clear();
  clusters.clear();
  helices.clear();
  tracks.clear();
  vertices.clear();
  calobjects.clear();
  particles.clear();
  trkParAtVtx.clear();
  vtxErrMatrices.clear();
  objects.clear();
  reserved.clear();
}

CsDstEvent_04::~CsDstEvent_04()
{
  for(unsigned i = 0; i < objects.size(); i++)
    delete objects[i];
}

void CsDstEvent_04::Init(const CsRecoEvent& ev)
{

  CsRecoEvent& event = const_cast<CsRecoEvent&>(ev);
  runNumber    = event.getEventHeader()[1];
  eventNumber  = event.getEventHeader()[2]; // see CsEvent.cc:1359
  eventInBurst = event.getEventHeader()[4];
  burstNumber  = event.getEventHeader()[3];
  timeInSec    = event.getEventHeader()[6];
  timeInUSec   = event.getEventHeader()[7];
  errorCode    = event.getEventHeader()[8];
  triggerMask  = event.getEventHeader()[9];

  vector<CsParticle*> _particles = event.getParticles();
  list<CsVertex*>     _vertices  = event.getVertices();
  list<CsCluster*>    _clusters  = event.getClusters();
    
  int n_Particles  = _particles.size();
  int n_Vertices   = _vertices.size();
  int n_Clusters   = _clusters.size();
  int n_Tracks     = 0;
  int n_Helices    = 0; 
  int n_Calobjects = 0;

  // event data

  vector<float32> calo_cluster_info;


  const int* _scalers = event.getScalers();
  
  for( unsigned int i=0; i<8; i++ ) 
    {
      scalers[i] = _scalers[i];
    }

  triggerTime = event.getTriggerTime();

  vector<unsigned int> _hodoData = event.getHodoData();
  nHodoData = _hodoData.size();
  hodoData.resize( nHodoData );
  for( int i=0; i<nHodoData; i++ ) 
    {
      hodoData[i] = _hodoData[i];
    }

  // particles

  particles.resize( n_Particles );

  for( unsigned int i=0; i < _particles.size(); i++ ) 
    {
      particles[i] = new CsDstParticle_04;
      particles[i]->Init( *_particles[i] );

      const CsTrack* trk = _particles[i]->getTrack();
      if( trk != NULL ) 
	{
	  n_Tracks++;
	  //BG20021015. Now save only 1st helix
	  if( (trk->getHelices()).size() > 0 ) 
	    {
	      n_Helices++;
	    }
	  //BG20021015. End
	}
      n_Calobjects += (_particles[i]->getCalObjects()).size(); 
    }

  int jh = 0;  // helices counter
  int jt = 0;  // tracks counter
  int jv = 0;  // vertices counter
  int jc = 0;  // clusters counter
  int jo = 0;  // calobjects counter
  
  helices.resize( n_Helices );
  tracks.resize( n_Tracks ); 
  vertices.resize( n_Vertices );
  calobjects.resize( n_Calobjects );

  short* ttoc = NULL;
  map< CsCluster*, int, less<CsCluster*> > mc;
  // clusters only if fatness=10
  if( fatness == 10 ) 
    {
      clusters.resize( n_Clusters );
      ttoc = new short[n_Tracks*(n_Clusters+2)]; // it's more than enough...
      for( list<CsCluster*>::iterator i=_clusters.begin();
	   i!=_clusters.end(); i++ ) 
	{
	  clusters[jc] = new CsDstCluster_04;
	  clusters[jc]->Init( *(*i) );
	  mc[ *i ] = jc;
	  jc ++;
	}
    }

  short* vtot = new short[n_Vertices*(n_Tracks+2)]; // it's more than enough...
  int countttoc = 0;
  map< const CsTrack*, int, less<const CsTrack*> > mt;

  short* pto = new short[n_Particles*7]; 
  int countpto = 0;

  for( unsigned int i=0; i<_particles.size(); i++ ) 
    {
      const CsTrack* trk = _particles[i]->getTrack();
      vector<Reco::CalorimeterParticle*> calob = 
	_particles[i]->getCalObjects(); 

      pto[countpto++] = i;

      if( trk != NULL ) 
	{
	  vector<CsHelix> _helices = trk->getHelices();
	  //BG20021015. Now save only 1st helix
	  helices[jh] = new CsDstHelix_04;
	  helices[jh++]->Init( _helices[0] );
	  //BG20021015. End
	  tracks[jt] = new CsDstTrack_04;
	  tracks[jt]->Init( *trk );
	  if( trk->IsBeamTrack() ) 
	    tracks[jt]->SetAsBeam();
	  mt[ trk ] = jt;
	  if( fatness == 10 ) 
	    {
	      ttoc[countttoc++] = jt;
	      list<CsCluster*> myclusters = trk->getClusters();
	      ttoc[countttoc++] = myclusters.size();
	      for( list<CsCluster*>::iterator j=myclusters.begin();
		   j!=myclusters.end(); j++ ) 
		{
		  ttoc[countttoc++] = mc[ *j ];
		}
	    }
	  pto[countpto++] = jt;
	  jt++;
	}
      else 
	{
	  pto[countpto++] = -1;
	}

      pto[countpto++] = calob.size();
      if( ! calob.empty() ) 
	{
	  for( unsigned int j=0; j<calob.size(); j++ ) 
	    {
	      calo_cluster_info.push_back(calob[j]->GetClusterSize());
	      if(calob[j]->HasTime())
		calo_cluster_info.push_back(calob[j]->GetTime());
	      else
		calo_cluster_info.push_back(1.e+7);
	      calobjects[jo] = new CsDstCalobject_04;
	      calobjects[jo]->Init(*calob[j]);
	      pto[countpto++] = jo;
	      jo++;
	    }
	}
    }

  int countvtot = 0;
  
  for( list<CsVertex*>::iterator i=_vertices.begin(); i!=_vertices.end(); i++ ) 
    {
      vertices[jv] = new CsDstVertex_04;
      vertices[jv]->Init( *(*i) );
      vtot[countvtot++] = jv;
      list<CsTrack*> myTks = (*i)->getTracks();
      vtot[countvtot++] = myTks.size();
      for( list<CsTrack*>::iterator j=myTks.begin(); j!=myTks.end(); j++ ) 
	{
	  CsTrack* trk  = (*j);
	  vtot[countvtot++] = mt[ trk ];
	  Cs3Vector trkparatv(0, 0, 0);
	  CsDst3Vector* v = new CsDst3Vector;
	  (*i)->getPar( trk, trkparatv );
	  v->Init(trkparatv);
	  trkParAtVtx.push_back(v);
	}
      vector<HepMatrix*> mtx = (*i)->getCov();
      for( unsigned int j=0; j<mtx.size(); j++ ) 
	{
	  CsDst3Matrix* m = new CsDst3Matrix;
	  m->Init(*(mtx[j]));
	  vtxErrMatrices.push_back( m );
	}
      jv++;  
    }

  if( fatness == 10 ) 
    {
      trkToClu.resize( countttoc );
      for( int i=0; i<countttoc; i++ ) 
	{
	  trkToClu[i] = ttoc[i];
	}
    }

  vtxToTrk.resize( countvtot );
  for( int i=0; i<countvtot; i++ ) 
    {
      vtxToTrk[i] = vtot[i];
    }

  partTkCo.resize( countpto );
  for( int i=0; i<countpto; i++ ) 
    {
      partTkCo[i] = pto[i];
    }

  delete[] vtot;
  if( fatness == 10 ) delete[] ttoc;
  delete[] pto;

  nVtxToTrk       = vtxToTrk.size();
  nTrkToClu       = trkToClu.size();
  nPartTkCo       = partTkCo.size();
  nClusters       = clusters.size();
  nHelices        = helices.size();
  nTracks         = tracks.size();
  nVertices       = vertices.size();
  nCalobjects     = calobjects.size();
  nParticles      = particles.size();
  nTrkParAtVtx    = trkParAtVtx.size();
  nVtxErrMatrices = vtxErrMatrices.size();
  objects.resize(nClusters+nHelices+nTracks+nVertices+
		 nCalobjects+nParticles+nTrkParAtVtx+nVtxErrMatrices);
  int n = 0;
  for(int i = 0; i < nClusters; i++, n++)
    objects[n] = clusters[i];
  for(int i = 0; i < nHelices; i++, n++)
    objects[n] = helices[i];
  for(int i = 0; i < nTracks; i++, n++)
    objects[n] = tracks[i];
  for(int i = 0; i < nVertices; i++, n++)
    objects[n] = vertices[i];
  for(int i = 0; i < nCalobjects; i++, n++)
    objects[n] = calobjects[i];
  for(int i = 0; i < nParticles; i++, n++)
    objects[n] = particles[i];
  for(int i = 0; i < nTrkParAtVtx; i++, n++)
    objects[n] = trkParAtVtx[i];
  for(int i = 0; i < nVtxErrMatrices; i++, n++)
    objects[n] = vtxErrMatrices[i];
  if(!SetAssociations())
    {
      cerr << "CsDstEvent_04::Init: Association error." << endl;
      exit(1); // should not be ever
    }
  nReserved = 1+calo_cluster_info.size();
  reserved.resize(nReserved);
  reserved[0] = (float32)ev.getTimeInSpill();
  for(unsigned i = 0; i < calo_cluster_info.size(); i++)
    {
      reserved[1+i] = calo_cluster_info[i];
    }
}

void CsDstEvent_04::GetEvent(CsRecoEvent& event) const
{

  event.setEventHeader(runNumber,1);
  event.setEventHeader(eventNumber,2);
  event.setEventHeader(eventInBurst,4);
  event.setEventHeader(burstNumber,3);
  event.setEventHeader(timeInSec,6);
  event.setEventHeader(timeInUSec,7);
  event.setEventHeader(errorCode,8);
  event.setEventHeader(triggerMask,9);

  vector<CsCluster*>   _clusters;
  vector<CsTrack*>     _tracks;
  vector<CsHelix*>     _helices;
  vector<CsVertex*>    _vertices;
  vector<CsParticle*>  _particles;
  vector<Reco::CalorimeterParticle*> _calobjs;

  event.clear();

  event.setTimeInSpill(double(reserved[0]));

  vector<float32> calo_cluster_info;
  if(nReserved-1 != 0)
    {
      calo_cluster_info.resize(nReserved-1);
      for(unsigned i = 0; i < calo_cluster_info.size(); i++)
	calo_cluster_info[i] = reserved[i-1];
    }
  else
    {
      calo_cluster_info.resize(nCalobjects*2);
      for(int i = 0; i < nCalobjects; i++)
	{
	  calo_cluster_info[2*i] = 0;
	  calo_cluster_info[2*i+1] = 1e+7;
	}
    }
  event.setTriggerTime(triggerTime);
  event.setScalers(scalers);

  for(int i = 0; i < nHodoData; i++)
    event.addHodoDatum(hodoData[i]);

  for(int i = 0; i < nClusters; i++)
    {
      _clusters.push_back(clusters[i]->GetCluster());
    }

  for(int i = 0; i < nHelices; i++)
    {
      _helices.push_back(helices[i]->GetHelix());
    }

  if(nTracks)
    {
      list<CsTrack*>      trs;
      for(int i = 0; i < nTracks; i++)
	{
	  trs.push_back(tracks[i]->GetTrack());
	}
      event.setTracks(trs);
      // important REMOVE LOCAL TRACKS: CsRecoEvent made a copy...
      for(list<CsTrack*>::iterator i = trs.begin(); i != trs.end(); i++)
	delete *i;
      list<CsTrack*>& trs1 = event.getTracks();
      int n = 0;
      for(list<CsTrack*>::iterator i = trs1.begin(); i != trs1.end(); i++, n++)
	{
	  _tracks.push_back(*i);
	  (*i)->addHelix(*(const CsHelix*)_helices[n]); // for DST4 only 1 helix per track
	}
    }

  // Delete all helices:
  for(unsigned i = 0; i < helices.size(); i++)
    delete _helices[i];

  for(int i = 0; i < nVertices; i++)
    {
      _vertices.push_back(vertices[i]->GetVertex());
      for(unsigned j = 0; j < vta[i].cov.size(); j++)
	{
	  _vertices[i]->addCov(vtxErrMatrices[vta[i].cov[j]]->GetHepMatrix());
	}
      for(unsigned j = 0; j < vta[i].trackParams.size(); j++)
	{
	  _tracks[vta[i].tracks[j]]->addVertex(_vertices[vta[i].vertex]);
	  _vertices[vta[i].vertex]->addTrack(_tracks[vta[i].tracks[j]]);
	  const float64* par = trkParAtVtx[vta[i].trackParams[j]]->GetVector();
	  Cs3Vector vec(par[0],par[1],par[2]);
	  _vertices[i]->addTrackAtVertex(_tracks[vta[i].tracks[j]],vec);
	}
    }

  for(int i = 0; i < nCalobjects; i++)
    {
      Reco::CalorimeterParticle* calPat = calobjects[i]->GetCalobject();
      calPat->SetClusterSize((size_t)calo_cluster_info[2*i]);
      if(calo_cluster_info[2*i+1] != 1e+7)
	calPat->SetTime(calo_cluster_info[2*i+1]);
      _calobjs.push_back(calPat);
    }

  if(fatness == 10)
    {
      for(unsigned i = 0; i < tca.size(); i++)
	{
	  for(unsigned j = 0; j < tca[i].clusters.size(); j++)
	    _tracks[tca[i].track]->addCluster(*_clusters[tca[i].clusters[j]]);
	}
    }
  else // Add dummy clusters list (it's neccessary for testing DST only)
    {
      for(unsigned i = 0; i < _tracks.size(); i++)
	{
	  int nCl = tracks[i]->GetNClusters();
	  list<CsCluster*> cl_list;
	  cl_list.resize(nCl);
	  _tracks[i]->setClusters(cl_list);
	}
    }

  if(_calobjs.size())
    {
      vector<Reco::CalorimeterParticle> cobjs;
      for(unsigned i = 0; i < _calobjs.size(); i++)
	{
	  cobjs.push_back(*_calobjs[i]);
  	  delete _calobjs[i];
	}
      event.setCalObjsVector(cobjs);
      _calobjs.clear();
      _calobjs = event.getCalObjs();
    }

  for(unsigned i = 0; i < ptca.size(); i++)
    {
      if(ptca[i].track != -1)
	_particles.push_back(new CsParticle(_tracks[ptca[i].track]));
      else
	_particles.push_back(new CsParticle);
      _particles[i]->setType(particles[i]->GetPartType());
      _particles[i]->setName(particles[i]->GetPartName());
      for(unsigned j = 0; j < ptca[i].calos.size(); j++)
	_particles[i]->addCalobj(_calobjs[ptca[i].calos[j]]);
    }

  if(_particles.size())
    {
      event.setParticles(_particles);
    }

  if(_clusters.size())
    {
      for(unsigned i = 0; i < _clusters.size(); i++)
	event.addCluster(*(_clusters[i]));
    }

  if(_vertices.size())
    {
      list<CsVertex*> vts;
      for(unsigned i = 0; i < _vertices.size(); i++)
	vts.push_back(_vertices[i]);
      event.setVertices(vts);
    }
}

bool CsDstEvent_04::Load(std::istream& file)
{
  uint32 hType = 0;
  file.read((char*)&hType,         sizeof(hType));
  if(hType != EVENT_HEADER_TYPE_04)
    {
      cerr << "CsDstEvent_04::Load(istream): DST reading error." << endl;
      return false;
    }
  file.read((char*)&eventNumber,   sizeof(eventNumber));
  file.read((char*)&eventInBurst,  sizeof(eventInBurst));
  file.read((char*)&burstNumber,   sizeof(burstNumber));
  file.read((char*)&timeInSec,     sizeof(timeInSec));
  file.read((char*)&timeInUSec,    sizeof(timeInUSec));
  file.read((char*)&errorCode,     sizeof(errorCode));
  file.read((char*)&triggerMask,   sizeof(triggerMask));
  file.read((char*)&extVars[0],    sizeof(extVars[0])*extVars.size());

  file.read((char*)&nVtxToTrk,sizeof(nVtxToTrk));
  file.read((char*)&nTrkToClu,sizeof(nTrkToClu));
  file.read((char*)&nPartTkCo,sizeof(nPartTkCo));
  file.read((char*)&nClusters,sizeof(nClusters));
  file.read((char*)&nHelices,sizeof(nHelices));
  file.read((char*)&nTracks,sizeof(nTracks));
  file.read((char*)&nVertices,sizeof(nVertices));
  file.read((char*)&nCalobjects,sizeof(nCalobjects));
  file.read((char*)&nParticles,sizeof(nParticles));
  file.read((char*)&nTrkParAtVtx, sizeof(nTrkParAtVtx));
  file.read((char*)&nVtxErrMatrices, sizeof(nVtxErrMatrices));
  file.read((char*)&scalers[0], sizeof(scalers));
  file.read((char*)&nHodoData, sizeof(nHodoData));
  file.read((char*)&triggerTime, sizeof(triggerTime));
  vtxToTrk.resize(nVtxToTrk);
  trkToClu.resize(nTrkToClu);
  partTkCo.resize(nPartTkCo);
  file.read((char*)&fatness,sizeof(fatness));
  file.read((char*)&vtxToTrk[0], nVtxToTrk*sizeof(int16));
  file.read((char*)&trkToClu[0], nTrkToClu*sizeof(int16));
  file.read((char*)&partTkCo[0], nPartTkCo*sizeof(int16));
  hodoData.resize(nHodoData);
  file.read((char*)&hodoData[0], nHodoData*sizeof(uint32));
  clusters.resize(nClusters);
  helices.resize(nHelices);
  tracks.resize(nTracks);
  vertices.resize(nVertices);
  calobjects.resize(nCalobjects);
  particles.resize(nParticles);
  trkParAtVtx.resize(nTrkParAtVtx);
  vtxErrMatrices.resize(nVtxErrMatrices);
  objects.clear();
  for(int i = 0; i < nClusters; i++)
    {
      clusters[i] = new CsDstCluster_04;
      objects.push_back(clusters[i]);
      if(!clusters[i]->Load(file)) return false;
    }
  for(int i = 0; i < nHelices; i++)
    {
      helices[i] = new CsDstHelix_04;
      objects.push_back(helices[i]);
      if(!helices[i]->Load(file)) return false;
    }
  for(int i = 0; i < nTracks; i++)
    {
      tracks[i] = new CsDstTrack_04;
      objects.push_back(tracks[i]);
      if(!tracks[i]->Load(file)) return false;
    }
  for(int i = 0; i < nVertices; i++)
    {
      vertices[i] = new CsDstVertex_04;
      objects.push_back(vertices[i]);
      if(!vertices[i]->Load(file)) return false;
    }
  for(int i = 0; i < nCalobjects; i++)
    {
      calobjects[i] = new CsDstCalobject_04;
      objects.push_back(calobjects[i]);
      if(!calobjects[i]->Load(file)) return false;
    }
  for(int i = 0; i < nParticles; i++)
    {
      particles[i] = new CsDstParticle_04;
      objects.push_back(particles[i]);
      if(!particles[i]->Load(file)) return false;
    }
  for(int i = 0; i < nTrkParAtVtx; i++)
    {
      trkParAtVtx[i] = new CsDst3Vector;
      objects.push_back(trkParAtVtx[i]);
      if(!trkParAtVtx[i]->Load(file)) return false;
    }
  for(int i = 0; i < nVtxErrMatrices; i++)
    {
      vtxErrMatrices[i] = new CsDst3Matrix;
      objects.push_back(vtxErrMatrices[i]);
      if(!vtxErrMatrices[i]->Load(file)) return false;
    }
  file.read((char*)&nReserved,sizeof(nReserved));
  reserved.resize(nReserved);
  file.read((char*)&reserved[0], nReserved*sizeof(float32));
  if(file.fail()) return false;
  return SetAssociations();
}

bool CsDstEvent_04::Load(FlatFile& file)
{
  uint32 hType = 0;
  if(!file.Read((char*)&hType,          sizeof(hType)))                    return false;
  if(hType != EVENT_HEADER_TYPE_04)
    {
      cerr << "CsDstEvent_04::Load(): DST reading error." << endl;
      return false;
    }
  if(!file.Read((char*)&eventNumber,   sizeof(eventNumber)))               return false;
  if(!file.Read((char*)&eventInBurst,  sizeof(eventInBurst)))              return false;
  if(!file.Read((char*)&burstNumber,   sizeof(burstNumber)))               return false;
  if(!file.Read((char*)&timeInSec,     sizeof(timeInSec)))                 return false;
  if(!file.Read((char*)&timeInUSec,    sizeof(timeInUSec)))                return false;
  if(!file.Read((char*)&errorCode,     sizeof(errorCode)))                 return false;
  if(!file.Read((char*)&triggerMask,   sizeof(triggerMask)))               return false;
  if(!file.Read((char*)&extVars[0],    sizeof(extVars[0])*extVars.size())) return false;

  if(!file.Read((char*)&nVtxToTrk,sizeof(nVtxToTrk))) return false;
  if(!file.Read((char*)&nTrkToClu,sizeof(nTrkToClu))) return false;
  if(!file.Read((char*)&nPartTkCo,sizeof(nPartTkCo))) return false;
  if(!file.Read((char*)&nClusters,sizeof(nClusters))) return false;
  if(!file.Read((char*)&nHelices,sizeof(nHelices))) return false;
  if(!file.Read((char*)&nTracks,sizeof(nTracks))) return false;
  if(!file.Read((char*)&nVertices,sizeof(nVertices))) return false;
  if(!file.Read((char*)&nCalobjects,sizeof(nCalobjects))) return false;
  if(!file.Read((char*)&nParticles,sizeof(nParticles))) return false;
  if(!file.Read((char*)&nTrkParAtVtx, sizeof(nTrkParAtVtx))) return false;
  if(!file.Read((char*)&nVtxErrMatrices, sizeof(nVtxErrMatrices))) return false;
  if(!file.Read((char*)&scalers[0], sizeof(scalers))) return false;
  if(!file.Read((char*)&nHodoData, sizeof(nHodoData))) return false;
  if(!file.Read((char*)&triggerTime, sizeof(triggerTime))) return false;
  vtxToTrk.resize(nVtxToTrk);
  trkToClu.resize(nTrkToClu);
  partTkCo.resize(nPartTkCo);
  if(!file.Read((char*)&fatness,sizeof(fatness))) return false;
  if(!file.Read((char*)&vtxToTrk[0], nVtxToTrk*sizeof(int16))) return false;
  if(!file.Read((char*)&trkToClu[0], nTrkToClu*sizeof(int16))) return false;
  if(!file.Read((char*)&partTkCo[0], nPartTkCo*sizeof(int16))) return false;
  hodoData.resize(nHodoData);
  if(!file.Read((char*)&hodoData[0], nHodoData*sizeof(uint32))) return false;
  clusters.resize(nClusters);
  helices.resize(nHelices);
  tracks.resize(nTracks);
  vertices.resize(nVertices);
  calobjects.resize(nCalobjects);
  particles.resize(nParticles);
  trkParAtVtx.resize(nTrkParAtVtx);
  vtxErrMatrices.resize(nVtxErrMatrices);
  objects.clear();
  for(int i = 0; i < nClusters; i++)
    {
      clusters[i] = new CsDstCluster_04;
      objects.push_back(clusters[i]);
      if(!clusters[i]->Load(file)) return false;
    }
  for(int i = 0; i < nHelices; i++)
    {
      helices[i] = new CsDstHelix_04;
      objects.push_back(helices[i]);
      if(!helices[i]->Load(file)) return false;
    }
  for(int i = 0; i < nTracks; i++)
    {
      tracks[i] = new CsDstTrack_04;
      objects.push_back(tracks[i]);
      if(!tracks[i]->Load(file)) return false;
    }
  for(int i = 0; i < nVertices; i++)
    {
      vertices[i] = new CsDstVertex_04;
      objects.push_back(vertices[i]);
      if(!vertices[i]->Load(file)) return false;
    }
  for(int i = 0; i < nCalobjects; i++)
    {
      calobjects[i] = new CsDstCalobject_04;
      objects.push_back(calobjects[i]);
      if(!calobjects[i]->Load(file)) return false;
    }
  for(int i = 0; i < nParticles; i++)
    {
      particles[i] = new CsDstParticle_04;
      objects.push_back(particles[i]);
      if(!particles[i]->Load(file)) return false;
    }
  for(int i = 0; i < nTrkParAtVtx; i++)
    {
      trkParAtVtx[i] = new CsDst3Vector;
      objects.push_back(trkParAtVtx[i]);
      if(!trkParAtVtx[i]->Load(file)) return false;
    }
  for(int i = 0; i < nVtxErrMatrices; i++)
    {
      vtxErrMatrices[i] = new CsDst3Matrix;
      objects.push_back(vtxErrMatrices[i]);
      if(!vtxErrMatrices[i]->Load(file)) return false;
    }
  if(!file.Read((char*)&nReserved,sizeof(nReserved))) return false;
  reserved.resize(nReserved);
  if(!file.Read((char*)&reserved[0], nReserved*sizeof(float32))) return false;
  return SetAssociations();  
}

bool CsDstEvent_04::IsEquil(const CsDstObject* obj) const
{
  if(type != obj->type) return false;
  const CsDstEvent_04* event = (const CsDstEvent_04*)obj;

  if(eventNumber != event->eventNumber) return false;
  if(eventInBurst != event->eventInBurst) return false;
  if(burstNumber != event->burstNumber) return false;
  if(timeInSec != event->timeInSec) return false;
  if(timeInUSec != event->timeInUSec) return false;
  if(errorCode != event->errorCode) return false;
  if(triggerMask != event->triggerMask) return false;
  for(int i = 0; i < 3; i++)
    if(extVars[i] != event->extVars[i]) return false;

  if(fatness != event->fatness) return false;
  if(nVtxToTrk != event->nVtxToTrk) return false;
  if(nTrkToClu != event->nTrkToClu) return false;
  if(nPartTkCo != event->nPartTkCo) return false;
  if(nClusters != event->nClusters) return false;
  if(nHelices != event->nHelices) return false;
  if(nTracks != event->nTracks) return false;
  if(nVertices != event->nVertices) return false;
  if(nCalobjects != event->nCalobjects) return false;
  if(nParticles != event->nParticles) return false;
  if(nTrkParAtVtx != event->nTrkParAtVtx) return false;
  if(nVtxErrMatrices != event->nVtxErrMatrices) return false;
  if(nHodoData != event->nHodoData) return false;
  for(int i = 0; i < 8; i++)
    if(scalers[i] != event->scalers[i]) return false;
  for(int i = 0; i < nHodoData; i++)
    if(hodoData[i] != event->hodoData[i]) return false;
  for(int i = 0; i < nVtxToTrk; i++)
    {
      if(vtxToTrk[i] != event->vtxToTrk[i]) 
	return false;
    }
  for(int i = 0; i < nTrkToClu; i++)
    {
      if(trkToClu[i] != event->trkToClu[i]) 
	return false;
    }
  for(int i = 0; i < nPartTkCo; i++)
    {
      if(partTkCo[i] != event->partTkCo[i]) 
	return false;
    }
  for(unsigned i = 0; i < objects.size(); i++)
    {
      if(!objects[i]->IsEquil(event->objects[i]))
	{
	  cout << "Object: " << i << endl;
	  objects[i]->Print(cout);
	  cout << "event->Object: " << endl;
	  event->objects[i]->Print(cout);
	  cout << "*******************" << endl;
	  return false;
	}
    }
  if(nReserved != event->nReserved) return false;
  for(int i = 0; i < nReserved; i++)
    if(reserved[i] != event->reserved[i]) return false;
  return true;
}

void CsDstEvent_04::Print(std::ostream& msgout) const
{
  char MSG_STR[1024];
  msgout << "\tEvent number:             " << eventNumber                   << "\n";
  msgout << "\tNumber of event in burst: " << eventInBurst                  << "\n";
  msgout << "\tBurst number:             " << burstNumber                   << "\n";
  msgout << "\tEvent time:               " << CsTime(timeInSec, timeInUSec) << "\n";
  msgout << "\tError code:               " << errorCode                     << "\n";
  msgout << "\tTrigger mask:             " << triggerMask                   << "\n";
  msgout << "\tExternal variable [0]:    " << extVars[0]                    << "\n";
  msgout << "\tExternal variable [1]:    " << extVars[1]                    << "\n";
  msgout << "\tExternal variable [2]:    " << extVars[2]                    << "\n";

  msgout << "Fatness:               " << (int)fatness << "\n";
  msgout << "Number of clusters:    " << nClusters << "\n";
  msgout << "Number of helices:     " << nHelices << "\n";
  msgout << "Number of tracks:      " << nTracks << "\n";
  msgout << "Number of vertices:    " << nVertices << "\n";
  msgout << "Number of cal.objects: " << nCalobjects << "\n";
  msgout << "Number of particles:   " << nParticles << "\n";
  msgout << "Number of TrkParAtVtx objects: " << nTrkParAtVtx << "\n";
  msgout << "Number of error matrices:      " << nVtxErrMatrices << "\n";

  for(int i = 0; i < nClusters; i++)
    {
      msgout << "Cluster # " << i << "\n";
      clusters[i]->Print(msgout);
    }
  for(int i = 0; i < nHelices; i++)
    {
      msgout << "Helix # " << i << "\n";
      helices[i]->Print(msgout);
    }
  for(int i = 0; i < nTracks; i++)
    {
      msgout << "Track # " << i << "\n";
      tracks[i]->Print(msgout);
    }
  for(int i = 0; i < nVertices; i++)
    {
      msgout << "Vertex # " << i << "\n";
      vertices[i]->Print(msgout);
    }
  for(int i = 0; i < nCalobjects; i++)
    {
      msgout << "Calorimeter object # " << i << "\n";
      calobjects[i]->Print(msgout);
    }
  for(int i = 0; i < nParticles; i++)
    {
      msgout << "Particle # " << i << "\n";
      particles[i]->Print(msgout);
    }
  for(int i = 0; i < nTrkParAtVtx; i++)
    {
      msgout << " Particle\'s track # " << i << " at vertex parameters: " << "\n";
      trkParAtVtx[i]->Print(msgout);
    }
  for(int i = 0; i < nVtxErrMatrices; i++)
    {
      msgout << "Covariant error matrix for vertex # " << i << "\n";
      vtxErrMatrices[i]->Print(msgout);
    }
  msgout << "Scaler\'s data:\n{ ";
  for(int i = 0; i < 8; i++)
    {
      sprintf(MSG_STR,"%8d ", scalers[i]); 
      msgout << MSG_STR;
    }
  msgout << "}\n";
  msgout << "\n Hodos data:";
  for(int i = 0; i < nHodoData; i++)
    {
      if(i%8 == 0) msgout << "\n";
      sprintf(MSG_STR,"%08X ", hodoData[i]);
      msgout << MSG_STR;
    }
  msgout << "\n";
  if(nVtxToTrk != 0)
    {
      msgout << "Vertex to tracks:\n  ";
      for(unsigned i = 0; i < vta.size(); i++)
	{
	  msgout << "\t* Vertex # " << vta[i].vertex << ": ";
	  msgout << "number of tracks: " << vta[i].tracks.size() << ", ";
	  msgout << "track\'s numbers: ";
	  for(unsigned j = 0; j < vta[i].tracks.size(); j++)
	    {
	      if(j != 0) msgout << ", ";
	      msgout << vta[i].tracks[j];
	    }
	  msgout << "\n";
	}
      msgout << "\n";
    }
  if(nTrkToClu != 0)
    {
      msgout << "Track to clusters:\n  ";
      for(unsigned i = 0; i < tca.size(); i++)
	{
	  msgout << "\t* Track # " << tca[i].track << ": ";
	  msgout << "number of clusters: " << tca[i].clusters.size() << ", ";
	  msgout << "cluster\'s numbers: ";
	  for(unsigned j = 0; j < tca[i].clusters.size(); j++)
	    {
	      if(j != 0) msgout << ", ";
	      msgout << tca[i].clusters[j];
	    }
	  msgout << "\n";
	}
      msgout << "\n";
    }
  if(nPartTkCo != 0)
    {
      msgout << "\tParticle\'s tracks to calorimeter:\n  ";
      for(unsigned i = 0; i < ptca.size(); i++)
	{
	  msgout << "\t* Particle # " << ptca[i].particle << ": ";
	  if(ptca[i].track != -1)
	    msgout << "track # " << ptca[i].track << ", ";
	  else
	    msgout << "no track, ";
	  msgout << "number of cal.objects: ";
	  msgout << ptca[i].calos.size();
	  if(ptca[i].calos.size())
	    {
	      msgout << ", cal.object\'s numbers: ";
	      for(unsigned j = 0; j < ptca[i].calos.size(); j++)
		{
		  if(j != 0) msgout << ", ";
		  msgout << ptca[i].calos[j];
		}
	    }
	  msgout << "\n";
	}
      msgout << "\n";
    }
  if(nReserved)
    {
      msgout << "Reserved words: \nTime in spill: " << reserved[0]
	     << endl;
    }
}

void CsDstEvent_04::GetBuffer(std::vector<uint8>& buffer) const
{
  uint32 hType = EVENT_HEADER_TYPE_04;
  AddValueToVector(buffer, &hType,         sizeof(hType));
  AddValueToVector(buffer, &eventNumber,   sizeof(eventNumber));
  AddValueToVector(buffer, &eventInBurst,  sizeof(eventInBurst));
  AddValueToVector(buffer, &burstNumber,   sizeof(burstNumber));
  AddValueToVector(buffer, &timeInSec,     sizeof(timeInSec));
  AddValueToVector(buffer, &timeInUSec,    sizeof(timeInUSec));
  AddValueToVector(buffer, &errorCode,     sizeof(errorCode));
  AddValueToVector(buffer, &triggerMask,   sizeof(triggerMask));
  AddValueToVector(buffer, &extVars[0],    sizeof(extVars[0])*extVars.size());

  AddValueToVector(buffer, &nVtxToTrk,  sizeof(nVtxToTrk));
  AddValueToVector(buffer, &nTrkToClu,  sizeof(nTrkToClu));
  AddValueToVector(buffer, &nPartTkCo,  sizeof(nPartTkCo));
  AddValueToVector(buffer, &nClusters,  sizeof(nClusters));
  AddValueToVector(buffer, &nHelices,   sizeof(nHelices));
  AddValueToVector(buffer, &nTracks,    sizeof(nTracks));
  AddValueToVector(buffer, &nVertices,  sizeof(nVertices));
  AddValueToVector(buffer, &nCalobjects,sizeof(nCalobjects));
  AddValueToVector(buffer, &nParticles, sizeof(nParticles));
  AddValueToVector(buffer, &nTrkParAtVtx, sizeof(nTrkParAtVtx));
  AddValueToVector(buffer, &nVtxErrMatrices, sizeof(nVtxErrMatrices));
  AddValueToVector(buffer, scalers, sizeof(scalers));
  AddValueToVector(buffer, &nHodoData, sizeof(nHodoData));
  AddValueToVector(buffer, &triggerTime, sizeof(triggerTime));
  AddValueToVector(buffer, &fatness,sizeof(fatness));
  AddValueToVector(buffer, &vtxToTrk[0], nVtxToTrk*sizeof(int16));
  AddValueToVector(buffer, &trkToClu[0], nTrkToClu*sizeof(int16));
  AddValueToVector(buffer, &partTkCo[0], nPartTkCo*sizeof(int16));
  AddValueToVector(buffer, &hodoData[0], nHodoData*sizeof(uint32));
  for(unsigned i = 0; i < objects.size(); i++)
    {
      objects[i]->GetBuffer(buffer);
    }
  AddValueToVector(buffer, &nReserved, sizeof(nReserved));
  AddValueToVector(buffer, &reserved[0], nReserved*sizeof(float32));
}

bool CsDstEvent_04::SetAssociations()
{
  vta.clear();
  vta.resize(nVertices);
  int nVert = 0;
  int kp = 0;
  int km = 0;
  for(int i = 0; i < nVtxToTrk; i++, nVert++)
    {
      vta[nVert].vertex = vtxToTrk[i++];
      int nTr = vtxToTrk[i];
      int nm  = nTr*2+1;
      for(int j = 0; j < nTr; j++)
	{
	  i++;
	  vta[nVert].tracks.push_back(vtxToTrk[i]);
	  vta[nVert].trackParams.push_back(kp++);
	}
      vta[nVert].cov.resize(nm);
      for(int j = 0; j < nm; j++)
	vta[nVert].cov[j] = km++;
    }
  if(nVert != nVertices)
    {
      cerr << "CsDstEvent_04: vertex association decoding error." << endl;
      return false;
    }
  if(fatness == 10)
    {
      tca.clear();
      tca.resize(nTracks);
      int nTr = 0;
      for(int i = 0; i < nTrkToClu; i++, nTr++)
	{
	  tca[nTr].track = trkToClu[i++];
	  int nClu = trkToClu[i];
	  for(int j = 0; j < nClu; j++)
	    {
	      i++;
	      tca[nTr].clusters.push_back(trkToClu[i]);
	    }
	}
      if(nTr != nTracks)
	{
	  cerr << "CsDstEvent_04: track association decoding error." << endl;
	  return false;
	}
    }
  int nPat = 0;
  ptca.clear();
  ptca.resize(nParticles);
  for(int i = 0; i < nPartTkCo; i++, nPat++)
    {
      ptca[nPat].particle = partTkCo[i++];
      ptca[nPat].track = partTkCo[i++];
      int nCo = partTkCo[i];
      for(int j = 0; j < nCo; j++)
	{
	  i++;
	  ptca[nPat].calos.push_back(partTkCo[i]);
	}
    }
  if(nPat != nParticles)
    {
      cerr << "CsDstEvent_04: particle association decoding error." << endl;
      return false;
    }
  return true;
}

CsDstChunk_04::CsDstChunk_04()
{
  type = DST_CHUNK_TYPE_04;
  dstVersion = 4;
  chunkName = "";
  TBNames = "";
  extVarNames.clear(); 
  extVarNames.resize(3);  
}

void CsDstChunk_04::Init(const char* _chunkName )
{
  chunkName = _chunkName;
}

CsDstEvent* CsDstChunk_04::GetNextEvent(std::istream& file) const
{
  CsDstEvent* event = new CsDstEvent_04(0);
  if(!event->Load(file))
    {
      delete event;
      return NULL;
    }
  return event;  
}

CsDstEvent* CsDstChunk_04::GetNextEvent(FlatFile& file) const
{
  CsDstEvent* event = new CsDstEvent_04(0);
  if(!event->Load(file))
    {
      delete event;
      return NULL;
    }
  return event;  
}

CsDstEvent* CsDstChunk_04::GetNextEvent(std::vector<uint8>& buffer) const
{
  VectorBuffer buf(buffer);
  istream istr(&buf);
  CsDstEvent* event = new CsDstEvent_04(0);
  if(!event->Load(istr))
    {
      delete event;
      return NULL;
    }
  return event;
}

bool CsDstChunk_04::Load(std::istream& file)
{
  file.read((char*)&type,sizeof(type));
  if(type != DST_CHUNK_TYPE_04)
    {
      cerr << "CsDstChunk_04::Load: incorrect chunk type: " << type << endl;
      return false;
    }
  file.read((char*)&dstVersion,sizeof(dstVersion));
  if(dstVersion != 4)
    {
      cerr << "CsDstChunk_04::Load: incorrect DST version: " << dstVersion << endl;
      return false;
    }
  uint32 length = 0;
  file.read((char*)&length,sizeof(length));
  char* str = new char[length];
  file.read(str,length);
  chunkName = string(str);
  delete str;
  length = 0;
  file.read((char*)&length,sizeof(length));
  str = new char[length];
  file.read(str,length);
  TBNames = string(str);
  delete str;
  char var[12];
  for(int i = 0; i < 3; i++)
    {
      file.read(var,12);
      extVarNames[i] = string(var);
    }
  if(file.fail()) return false;
  return true;
}

bool CsDstChunk_04::Load(FlatFile& file)
{
  if(!file.Read((char*)&type,sizeof(type))) return false;
  if(type != DST_CHUNK_TYPE_04)
    {
      cerr << "CsDstChunk_04::Load: incorrect chunk type: " << type << endl;
      return false;
    }
  if(!file.Read((char*)&dstVersion,sizeof(dstVersion))) return false;
  if(dstVersion != 4)
    {
      cerr << "CsDstChunk_04::Load: incorrect DST version: " << dstVersion << endl;
      return false;
    }
  uint32 length = 0;
  if(!file.Read((char*)&length,sizeof(length))) return false;
  char* str = new char[length];
  if(!file.Read(str,length))
    {
      delete str;
      return false;
    }
  chunkName = string(str);
  delete str;
  length = 0;
  if(!file.Read((char*)&length,sizeof(length))) return false;
  str = new char[length];
  if(!file.Read(str,length))
    {
      delete str;
      return false;
    }
  TBNames = string(str);
  delete str;
  char var[12];
  for(int i = 0; i < 3; i++)
    {
      if(!file.Read(var,12)) return false;
      extVarNames[i] = string(var);
    }
  return true;
}

bool CsDstChunk_04::IsEquil(const CsDstObject* obj) const
{
  const CsDstChunk_04* cont = (const CsDstChunk_04*)obj;
  if(chunkName != cont->chunkName) return false;
  if(TBNames != cont->TBNames) return false;
  for(int i = 0; i < 3; i++)
    if(extVarNames[i] != cont->extVarNames[i]) return false;
  return true;
}

void CsDstChunk_04::Print(std::ostream& msgout) const
{
  msgout << "DST chunk:            " << chunkName << "\n";
  msgout << "TB names: \n{\n" << TBNames << "\n}\n";
  if(extVarNames[0][0])
    msgout << "External variable #0: " << extVarNames[0] << "\n";
  if(extVarNames[1][0])
    msgout << "External variable #1: " << extVarNames[1] << "\n";
  if(extVarNames[2][0])
    msgout << "External variable #2: " << extVarNames[2] << "\n";
}

void CsDstChunk_04::GetBuffer(std::vector<uint8>& buffer) const
{
  AddValueToVector(buffer, &type, sizeof(type));
  AddValueToVector(buffer, &dstVersion, sizeof(dstVersion));
  uint32 length = chunkName.size()+1;
  AddValueToVector(buffer, &length, sizeof(length));
  AddValueToVector(buffer, chunkName.c_str(), length);
  length = TBNames.size()+1;
  AddValueToVector(buffer, &length, sizeof(length));
  AddValueToVector(buffer, TBNames.c_str(), length);
  char vars[3][12];
  memset(vars,0,sizeof(vars));
  for(int i = 0; i < 3; i++)
    {
      strcpy(vars[i],extVarNames[i].c_str());
    }
  AddValueToVector(buffer, vars, sizeof(vars));
}

void CsDstChunk_04::SetExtVarName(int n, const char* name)
{
  if(n < 0 || n > 2) return;
  char var[12];
  strncpy(var,name,11);
  var[11] = 0;
  extVarNames[n] = string(var);
}
