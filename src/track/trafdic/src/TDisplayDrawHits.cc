// $Id: TDisplayDrawHits.cc 14069 2015-09-17 20:44:46Z lsilva $

/*!
  \brief Draw hits (THit objects) on a projection plane in the event display.
  \arg \a ifl: if !=0, no actual drawing, wait instead for a mouse click locating a hit and return its index (in "TEv::vecHit").

  - The hits' position along the coordinate axis measured by their parent detector is plotted, as opposed to the event display plane.
  - Drift hits are corrected for various effects (signal propagation, ...), as resulting for the reco processing.
  - The silent partner of coalesced hits is not drawn.
  - GP P-pixel hits are drawn in all of 4 proj.: X, Y and 2 extras, defined by "TraF" options.
  - Upon option, MP M-pixel hits are drawn in one extra proj., viz. along their v axis.
  - MC data: colour code (blue=default,green=ghost, e.g. mirror in drift).
  - Get the angle of the current event display proj. from "TSetup::vecProj".
*/

#include <cfloat>
#include <iostream>

#include "TDisplay.h"
#include "TSetup.h"
#include "TEv.h"
#include "TOpt.h"
#include "TConstants.h"
#include "higz.h"

using namespace std;

int TDisplay::DrawHits(int ifl)
{
  const TSetup& setup = TSetup::Ref();
  TEv& ev = TEv::Ref(); const vector<THit>& vHit = ev.vHit();

  float zc, x, y;
  ISMK(20);

  float px, py; if (ifl!=0) {
    int istat, nt; IRQLC(1,20,istat,nt,px,py);
    if (istat==0) return -1;
  }

  //        ********** ANGLE of CURRENT DISPLAYED PROJECTION **********
  int iAngle = (int)floor(setup.vProj()[Proj]/10.+.5);
  double Angle = iAngle*M_PI/180;

  int zone, j_picked; float dist, dist_min;
  for (zone = 0, j_picked = -1, dist_min = FLT_MAX;
       zone<(int)setup.vIplFirst().size(); zone++) {

    //      *************** LOOP ON ZONES ***************

    //                 ***** PIXELGEM/MMs *****
    // GP P-pixels could a priori be projected on any proj. In TraFDic, they are
    // available on 4 different proj.'s: 0 and 90 deg., and 2 others, which
    // +/-45 deg. defaults can be overriden by option on a per zone basis.
    double aPixProjs[] = {0,90,45,-45};
    if (TOpt::iCut[29]&1<<zone) { // Zones where special pixelGEM projection
      aPixProjs[2] = TOpt::dCut[82]; aPixProjs[3] = TOpt::dCut[83];
    }
    int drawPPixDet, i03, vProj, iproj; for (i03=drawPPixDet = 0; i03<=3; i03++)
      if (abs(floor(aPixProjs[i03]+.5)-iAngle)<=1) {
	drawPPixDet = 1; break;
      }
    // MP M-pixels, upon option: project them on one extra proj.: their v-axis.
    for (iproj = 0, vProj = -1; iproj<(int)setup.vProj().size(); iproj++) {
      int jAngle = (int)floor(setup.vProj()[iproj]/10.+.5);
      if (abs(abs(iAngle-jAngle)-90)<=1) {
	vProj = iproj;
	break; // Breaking here prevents the -90 proj., which contains only BMS, from being retained.
      }
    }

    for (int j = 0; j<(int)vHit.size(); j++) {

      // *********************************************************************
      // *************** LOOP ON ALL HITS ZONE AFTER ZONE ***************
      // *********************************************************************

      const THit &h = vHit[j];
      if (h.Status==-4) continue;  // Skip silent component of coalesced hit
      const TPlane &p = setup.vPlane(h.IPlane);
      if (h.IPlane<setup.vIplFirst()[zone]) continue;
      if (TOpt::ReMode[1]==2) {// THit's are ordered by increasing det. index
	
      }
      else if (setup.vIplLast()[zone]<h.IPlane) continue;
      if ((TOpt::Graph[5]&0x1) &&  // (Optionally) Skip switched OFF planes
	  p.IFlag==0) continue;
      const TDetect &d = setup.vDetect(p.IDetRef);
      int isPixDet = 0; // Note: better not rely on "TPlane::IFlag&0x30" to...
      // ...single out pixels, for we want them also when turned off: "IFlag=0".
      if (d.IType==29) isPixDet = 1; // GP P-pixel
      else if ((TOpt::Graph[5]&0x2) && d.IType==32)
	isPixDet = 2;                    // (upon option) MP M-pixel
      if (p.IProj!=Proj &&              // Skip all but current proj. except...
	  (isPixDet!=1 || !drawPPixDet) &&     // ...P-pixels
	  (isPixDet!=2 || p.IProj!=vProj))     // ...(optionally) M-pixels
	continue;

      x = d.X(0); zc = d.XR(2);
      if (TOpt::Graph[8]==0) {     // Simplistic cluster representation
	ISMKSC(0.40+0.0008*zc);
      }
      else {                       // ...else acount for cluster size
	// Nota bene: this is not properly taken care of for pixelGEM's 
	int clsiz = min(5, int(h.sDigits().size()));
	ISMKSC(0.40 + TOpt::Graph[8]*0.10*(clsiz-1));
      }

      if (!isPixDet) y = h.U;                      // ***** HIT COORDINATE *****
      else {
	double a = atan2(d.Sa,d.Ca)-Angle; y = h.U*cos(a)-h.V*sin(a);
      }

      if (ifl==0) {  // ********** DRAWING MODE **********
	ISPMCI(4); ISPLCI(4);
	if (ev.IsMC() &&// MC: Special color for hits not associated to MC track
	    h.IKine<0 && h.IKin2<0) ISPMCI(109);
	if (h.PtrClusMirr()!=NULL){ // Special case: drift detector hit
	  float xx[2]; float yy[2];
	  xx[0]=xx[1]=x; 
	  yy[0]=y-h.DeltaR;
	  yy[1]=y+h.DeltaR;
	  if (!ev.IsMC() || (ev.IsMC() && h.DRsignMC == -1)) IPMm(xx[0],yy[0]); 
	  if (!ev.IsMC() || (ev.IsMC() && h.DRsignMC == +1)) IPMm(xx[1],yy[1]);
	  IPLm(2,xx,yy); // connecting line
	}
	else IPMm(x,y);
      }
      else {       // ********** PICKING UP MODE **********
	if ((dist = sqrt((x-px)*(x-px)+(y-py)*(y-py)))<dist_min) { 
	  dist_min = dist; j_picked = j;
	}
      }
    } // End of loop over hits
  } // End loop on zones

  if (ifl==0) {
    IUWK(0,1); ISPLCI(1);
  }
  return j_picked;
}
