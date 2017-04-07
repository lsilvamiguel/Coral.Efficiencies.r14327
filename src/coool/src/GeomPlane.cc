
#include <cassert>
#include <cstring>
#include <iostream>
#include "TMath.h"
#include "GeomPlane.h"

#define max(a,b) (((a) > (b)) ? (a) : (b))

//ClassImp(GeomPlane);

const double GeomPlane::sqrt12 = sqrt(12.);

GeomPlane:: GeomPlane (int id,const char* name,int nwir,
		       double x,double y,double z,
		       double dx,double dy,double dz,
		       double ang,double pitch)
  : fId(id),fName(name),fNwir(nwir),fX(x),fY(y),fZ(z),
    fDx(dx),fDy(dy),fDz(dz),fAngle(ang), fHalfsize(nwir*pitch/2.) {

  char *gname = new char[strlen(name)+1];
  sprintf(gname, "%s", name);
 //fNameC used for SI reading the spacial resolution in GeomPlaneSili
  fNameC = gname;

  fPitchs[0] = pitch;
  
  fCos = cos(fAngle*TMath::Pi()/180.);
  fSin = sin(fAngle*TMath::Pi()/180.);
}

GeomPlane::~GeomPlane() {
}

void GeomPlane::AddSubPlane(int id, const char *name, Int_t nwires, 
			    double x, double y, double z, 
			    double dx, double dy, double dz,
			    double angle, double inpitch) {

  typedef std::map<int,double>::iterator IP;
  
  if(fName!=name)
    throw "GeomPlane::AddSubPlane : different names !";

  if( fabs(angle - fAngle) > 0.001 )
    throw "GeomPlane::AddSubPlane : different angles !";
  
  // CONTINUITY TEST (in WRS) 
  //  double tolerance=fPitchs[0]/10.;
  //  double tolerance=inpitch/10.;  // like in Coral (DN 19/5/2004)
  double tolerance=max(inpitch, fPitchs[0])/10.*6;  // to remove continuity test pb...
  
  
  // old half size of the active zone (WRS)
  //    double halfxs = 0;
  //    int lastw=0;
  //    double lastp=0.;
  //    for(IP i=fPitchs.begin(); i!=fPitchs.end(); i++) {
  //      halfxs += (i->first-lastw)*lastp;
  //      lastw = i->first;
  //      lastp = i->second;
  //    }
  //    halfxs += (fNwir - lastw)*lastp; 
  //    halfxs = halfxs/2.;
  
  // vector from old detector center to subdetector center in WRS
  double cscy =  GetCos()*(y-fY) + GetSin()*(z-fZ);
  double cscz = -GetSin()*(y-fY) + GetCos()*(z-fZ);

  double dist = sqrt((fabs(cscy)-(fHalfsize+nwires*inpitch/2.))
		     *(fabs(cscy)-(fHalfsize+nwires*inpitch/2.))
		     + cscz*cscz);  

  // for GP, the subdetector is perpendicular, see above
  if( fabs(dist)>tolerance ) {
    std::cerr<<"GeomPlane::AddSubPlane : error dist "<<dist<<" > tolerance "<<tolerance<<std::endl;
    std::cerr<<"                       : name "<<name<<" x "<<x<<" y "<<y<<" z "<<z<<std::endl;
    std::cerr<<"                       : angle "<<angle<<std::endl;
    std::cerr<<"                       : pitch "<<fPitchs[0]<<" inpitch "<<inpitch<<std::endl;
    throw "GeomPlane::AddSubPlane : Continuity test failed";
  }

  // storing subdet pitch, calculating center pos in WRS

  double cmwrsy =  GetCos()*fY + GetSin()*fZ;

  if(cscy>0) {
    fPitchs[fNwir] = inpitch;
    cmwrsy += nwires*inpitch/2.;
  }
  else {
    std::map<int, double> newpitch;
    newpitch[0] = inpitch;
    
    for(IP i=fPitchs.begin(); i!=fPitchs.end(); i++) {
      newpitch[i->first + nwires] = i->second;
    }
    fPitchs = newpitch;
    cmwrsy -= nwires*inpitch/2.;
  }
  fNwir += nwires;

  //new size
  fDy += dy;
  fHalfsize += nwires*inpitch/2.;

  // new center pos in MRS
  double cmwrsz = -GetSin()*fY + GetCos()*fZ;
  fY = GetCos()*cmwrsy - GetSin()*cmwrsz;
  fZ = GetSin()*cmwrsy + GetCos()*cmwrsz;
} 

double GeomPlane::Wire2Pos(double wire) const {
  int lastch=0;
  double lastp=0.;
  double dist=0.;
  //double firstpitch=0.;
  for(std::map<int,double>::const_iterator ipitch=fPitchs.begin(); ipitch!=fPitchs.end(); ipitch++ ) {
    //if(!firstpitch) firstpitch=ipitch->second;
    if(wire < ipitch->first) break;
    dist += (ipitch->first-lastch)*lastp;
    lastp = ipitch->second;
    lastch = ipitch->first;
  }
  //dist += (wire-lastch+0.5)*lastp - firstpitch/2.;
  dist += (wire-lastch+0.5)*lastp - fHalfsize;

  return dist;
}

double GeomPlane::Pitch(double wire) const {

  double pitch = 0;
  for(std::map<int,double>::const_iterator ipitch=fPitchs.begin(); ipitch!=fPitchs.end(); ipitch++ ) {
    if(ipitch->first > wire) 
      break;
    pitch = ipitch->second;
  }  
  return pitch;
}

double GeomPlane::Resolution(double wire) const {
  return Pitch(wire)/sqrt12;
}

std::set<CCluster1>& GeomPlane::Clusterize(const std::map<int,CDigit1>& digits) {

  fClusters.clear();

  for (std::map<int,CDigit1>::const_iterator id = digits.begin(); id!=digits.end(); ) {

    const CDigit1& digit = id->second;
    
    int firstchan = digit.ch;
    int lastchan = firstchan;
    int ndigits=1;
    int delta=0;
    do {
      id++;
      if (id!=digits.end()) {
        const CDigit1& htmp = id->second;
        delta = htmp.ch - lastchan;
        if (delta <= 1) {
	  ndigits++;
          lastchan = htmp.ch;
        }
      }
    } while ((id!=digits.end()) && (delta<=2));
    
    assert(lastchan>=firstchan);
    double wire = (firstchan+lastchan)/2.;
    double pos = Wire2Pos(wire);
    double res = Resolution(wire);
    
    // this is not optimal : cluster data is data of the first digit...
    fClusters.insert(CCluster1(GetID(),pos,ndigits,res,digit.dt));    
  }
  return fClusters;
}

