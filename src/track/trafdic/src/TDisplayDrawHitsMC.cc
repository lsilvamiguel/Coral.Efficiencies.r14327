// $Id: TDisplayDrawHitsMC.cc 13508 2012-07-10 17:21:50Z ybedfer $

/*!
 \brief Draw MC hits (THitMC objects) on a projection plane in the event display.
 \arg \a ifl: if !=0, no actual drawing, wait instead for a mouse click locating a MC hit and return its index (in "TEv::vecHitMC").

 - The MC coord's (from the MC truth bank, in MRS) are projected onto the axis of the TDetector the hit pertains to, so that they can be directly compared to their smeared counterpart.
  - GP P-pixel hits are drawn in all of 4 proj.: X, Y and 2 extras, defined by "TraF" options.
  - Upon option, MP M-pixel hits are drawn in one extra proj., viz. along their v axis.
 - Colour code = THitMC::Origin: magenta if original hit, cyan otherwise.
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
#include <math.h>

using namespace std;

int TDisplay::DrawHitsMC(int ifl)
{
  float px, py; if (ifl!=0) {
    int istat, nt; IRQLC(1,20,istat,nt,px,py);
    if (istat==0) return -1;
  }

  const TSetup &setup = TSetup::Ref();
   TEv &ev = TEv::Ref(); const vector<THitMC> &vHitMC = ev.vHitMC();

  //        ********** ANGLE of CURRENT DISPLAYED PROJECTION **********
  int iAngle = (int)floor(setup.vProj()[Proj]/10.+.5);
  double Angle = iAngle*M_PI/180, cosA = cos(Angle), sinA = sin(Angle);

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

    for (int j = 0; j<(int)vHitMC.size(); j++) {

      // *********************************************************************
      // *************** LOOP ON ALL HITS ZONE AFTER ZONE ***************
      // *********************************************************************

      const THitMC &h = vHitMC[j];
      const TPlane &p = setup.Id2Plane(h.IDet);
      if (p.IPlane<setup.vIplFirst()[zone]) continue;
      if (setup.vIplLast()[zone]<p.IPlane) continue;
      const TDetect &d = h.DetRef();

      float x ,uv;
      if (ev.vKine(h.IKine).Q()== 0) {  // Neutral
	ISPMCI(7); ISMK(1);
      }
      else {
	if (h.IOrig!=0) ISPMCI(110);
	else            ISPMCI(6);
	ISMK(30);
      }

      x = h.Xyz(0);
      //         ********** PROJECT TO CURRENT VIEW PLANE **********
      int isPixDet = 0; // Note: better not rely on "TPlane::IFlag&0x30" to...
      // ...single out pixels, for we want them also when turned off: "IFlag=0".
      if (d.IType==29) isPixDet = 1; // GP P-pixel
      else if ((TOpt::Graph[5]&0x2) && d.IType==32)
	isPixDet = 2;                    // (upon option) MP M-pixel
      if (p.IProj==Proj) {
	//      ***** DETECTOR'S PROJECTION == CURRENT VIEW PLANE *****
	// The MC coord's (MRS) are projected onto the detector's measurement
	// axis (WRS), rather than onto the proj. plane of the event display, so
	// that they can be directly compared to their smeared counterparts.
	ISMKSC(0.8); uv = d.Ca*h.Xyz(1) + d.Sa*h.Xyz(2);
      }
      else if ((isPixDet==1 && drawPPixDet) ||     // ...P-pixels
               (isPixDet==2 && p.IProj==vProj)) {  // ...(optionally) M-pixels
	// The MC coord's (MRS) are projected onto the proj. plane of the event
	// display.
	ISMKSC(0.8); uv = cos(Angle)*h.Xyz(1) + sin(Angle)*h.Xyz(2);
      }
      else if (DrOpt[4]!=-1) {  // ***** (OPTIONALLY) DRAW ALL PROJ. *****
	ISMKSC(0.5); uv = cosA*h.Xyz(1) + sinA*h.Xyz(2);
      }
      else continue; // Skip this MC hit

      if (ifl==0) IPMm(x,uv);
      else {
	if ((dist = sqrt((x-px)*(x-px)+(uv-py)*(uv-py)))<dist_min) { 
	  dist_min = dist; j_picked = j;
	}
      }
    } // End of loop over hits
  } // End loop over zones

  if (ifl==0) IUWK(0,1);

  return j_picked;
}
