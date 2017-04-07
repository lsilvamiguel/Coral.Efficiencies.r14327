// $Id: TDisplayDrawMCTracks.cc 13508 2012-07-10 17:21:50Z ybedfer $

/*!
  \brief Draw MC tracks (TKine objects) on a projection plane in the event display.
  \arg \a ifl: if !=0, no actual drawing, wait instead for a mouse click locating a MC track and return its index (in "TEv::vecKine").

  - 2 modes of play:
  i)  Polyline: through only (by default) those hits pertaining to proj. plane.
  ii) Extrap.: extrapolating from helix at vertex of origin.
  - In Polyline mode, P-pixel hits are taken into account in all of 4 proj.: X, Y and 2 extras, defined by "TraF" options.
  - In Polyline mode but only upon option, M-pixels are taken into account in one extra proj., viz. their v-axis.
  - MC data: colour code (grey=pileup,magenta otherwise), label w/ associated reco'd track (if any) index in addition to MC track index.
  - Get the angle of the current event display proj. from "TSetup::vecProj".
  - Draw the beam track even if it has no, MC, hits and hence cannot be delimited.
  - Options (cpp) for skipping off-time or pile-up tracks.
*/

#include <cfloat>
#include <iostream>

#include "TDisplay.h"
#include "TSetup.h"
#include "TEv.h"
#include "TOpt.h"
#include "TConstants.h"
#include "TAlgo.h"
#include "higz.h"
#include <math.h>

using namespace std;

int TDisplay::DrawMCTracks(int ifl)
{
  if (DrOpt[9]<0) { // ***** MC TRACK DRAWING IS OFF *****
    if (ifl!=0) 
      cout<<"TDisplay::DrawMCTracks ==> Track can not be picked up, as drawing if OFF."<<endl; 
    return -1;
  }
  
  const TSetup &setup = TSetup::Ref();
  TEv &ev = TEv::Ref();
  const vector<TKine> &vKine  = ev.vKine();
  const vector<THitMC> &vHitMC = ev.vHitMC();

  float px, py; if (ifl!=0) {
    int nt, istat; IRQLC(1,20,istat,nt,px,py);
    if (istat==0) return -1;
  } else {
  }

  //        ********** ANGLE of CURRENTLY DISPLAYED PROJECTION **********
  int iAngle = (int)floor(setup.vProj()[Proj]/10.+.5);
  double Angle = iAngle*M_PI/180, cosA = cos(Angle), sinA = sin(Angle);

  //         ********** RE-ORDER RECONSTRUCTION ZONES **********
  int iZone, nZones = (int)setup.vIplLast().size();
  int iZone2Zone[nZones], iplLasts[nZones];
  for (iZone = 0; iZone<nZones; iZone++) {
    iZone2Zone[iZone] = iZone; iplLasts[iZone] = setup.vIplLast()[iZone];
  }
  for (iZone = 0; iZone<nZones; iZone++) {
    int last = iplLasts[iZone], zone = iZone2Zone[iZone];
    for (int iZp = iZone+1; iZp<nZones; iZp++) if (iplLasts[iZp]<last) {
      iplLasts[iZone] = iplLasts[iZp];     iplLasts[iZp] = last;
      last = iplLasts[iZone];
      iZone2Zone[iZone] = iZone2Zone[iZp]; iZone2Zone[iZp] = zone;
      zone = iZone2Zone[iZone];
    }
  }

  //                 ***** PIXELGEM/MMs *****
  // GP P-pixels could a priori be taken into account in any proj. In TraFDic,
  // they are available on 4 different proj.'s: 0 and 90 deg., and 2 others,
  // which +/-45 deg. defaults can be overriden by option on a per zone basis.
  int zone; unsigned int pixZones;
  for (zone = 0, pixZones = 0; zone<(int)setup.vIplFirst().size(); zone++) {
    double aPixProjs[] = {0,90,45,-45};
    if (TOpt::iCut[29]&1<<zone)  { // Zones where special pixelGEM projection
      aPixProjs[2] = TOpt::dCut[82]; aPixProjs[3] = TOpt::dCut[83];
    }
    int i03; for (i03 = 0; i03<=3; i03++)
      if (abs(floor(aPixProjs[i03]+.5)-iAngle)<=1) {
	pixZones |= 1<<zone; break;
      }
  }
  // MP M-pixels, upon option: take them into account in one extra proj.: their
  // v-axis.
  int iproj, vProj; for (iproj = 0, vProj = -1; iproj<(int)setup.vProj().size();
			 iproj++) {
    int jAngle = (int)floor(setup.vProj()[iproj]/10.+.5);
    if (abs(abs(iAngle-jAngle)-90)<=1) {
      vProj = iproj;
      break; // Breaking here prevents the -90 proj., which contains only BMS, from being retained.
    }
  }

  float uv, x[2], y[2];
  int ik, id_picked; float dist, dist_min;
  for (ik = 0, id_picked = -1, dist_min = FLT_MAX; ik<(int)vKine.size(); ik++) {

    // ***********************************************************************
    // ************************* LOOP OVER MC TRACKS *************************
    // ***********************************************************************

    const TKine &k = vKine[ik];
    if (k.Q()==0) continue;                               // ***** SKIP NEUTRALS
    bool isPileup = k.isPileup();        // ***** PILE-UPS *****
    //#define TDisplay_SKIP_PILEUP
#ifdef TDisplay_SKIP_PILEUP
    if (isPileup) continue;                 // ***** (UPON DEFINE) SKIP PILE-UPS
#endif
    const vector<int> &vHitMCRef = k.vHitMCRef();
    vector<int>::const_iterator j;
#define TDisplay_SKIP_OFFTIME    // ***** (UPON DEFINE) SKIP FAR OFF-TIME TRACKS
#ifdef TDisplay_SKIP_OFFTIME
    // The idea is to draw only those tracks that stand reasonable chance to
    // be reconstructed.
    if (vHitMCRef.empty()) {
      // We still draw all !pileups, even if not reconstructible, for
      // completeness' sake.
      if (isPileup) continue;
    }
    else {
      //         ********** OFF-TIME? **********
      // - Determine track's time from its component hits.
      // - Derive time gate from detector type, whether VSAT, SAT, Drift or LAT.
      // - Only (some, cf. infra) tracking detectors (and hence, hits thereof)
      //  are considered. (Pixel)GEM are excluded, given that their timing
      //  can be smeared, when the "AmpCorrelationMC" feature is activated.
      // VSAT = SI(21), FI(22)
      // SAT  = MM(27)
      // LAT  = P(1,2), DC(15), ST(11:14,17), DW(16), MB(18)
      double sTime, tMnMn, tMxMx; int timed, s1;
      for (j = vHitMCRef.begin(), sTime=tMnMn=tMxMx = 0, s1=timed = 0;
		       j!=vHitMCRef.end(); j++) {
	const THitMC &h = vHitMC[*j]; const TDetect &d = h.DetRef();
	int t = d.IType, reliable = 1; static double tMn, tMx;
	if      (t==21 || t==22) { tMn = -20; tMx = 20; }  // VSAT
	else if (t==27)          { tMn = -40; tMx = 40; }  // SAT
	else if (t==1 || t==2)   { tMn = -60; tMx = 60; }  // MWPC
	else if (t==16)          { tMn = -20; tMx = 80; }  // DW
	else if (11<=t && t<=18) { tMn = -10; tMx = 60; }  // DC, DR, ST, MB
	else if (26<=t && t<=29) { 
	  reliable = 0;            tMn = -40; tMx = 40;    // GM, GP
	}
	if (reliable) {
	  if (!timed) { // Reset average time, now that reliable info at hand
	    sTime = 0; s1= 0; timed = 1; 
	  }
	  sTime += h.DeltaT; s1++;
	}
	else if (!timed) {// Accept unreliable info only if nothing else at hand
	  sTime += h.DeltaT; s1++;
	}
	if (tMn<tMnMn) tMnMn = tMn; if (tMx>tMxMx) tMxMx = tMx;
	//#  define TDisplay_DEBUG_MCTracks
#  ifdef TDisplay_DEBUG_MCTracks
	printf("%3d: %s %7.2f => %s\n",ik,d.Name.c_str(),h.DeltaT,
	       fabs(h.DeltaT)>tGate?"OUT":"IN");
#  endif
      }
      int inTime = -1; if (s1) {
	sTime /= s1; inTime = tMnMn<sTime && sTime<tMxMx ? 1 : 0;
      }
      if      (inTime==0) continue;    // Off-time => Skip
      else if (inTime==-1 &&
	       isPileup) continue;     // Undefined => Skip only if it's pileup
    }
#endif

    int color = isPileup ? 106 : 6; // MC track color
    char c[10]; string s, s1;       // Track label string

    if (k.sTrackID().size()==0) {     // ***** TRACK LABEL *****
      sprintf(c,"%u",ik); s = c;
    }
    else {
      sprintf(c,"%u(",ik); s = c;
      for (set<unsigned int>::iterator is = k.sTrackID().begin();
	   is!=k.sTrackID().end(); is++) {
	sprintf(c,"%u,",(*is)); s1 = c; s += s1;
      }
      s.replace(s.rfind(","),1,")");
    }
    float tga(0);

    //           *************** DRAW TRACK ***************
   
    if (DrOpt[5]<0) {      //   ********** POLYLINE **********
      float xx[2], yy[2]; int npt = 0;
      int nPats = vHitMCRef.size()/32+1; unsigned int hitsPats[nPats];
      memset((void*)hitsPats,0,nPats*sizeof(unsigned int));
      int iNext; do {
	int iref; static double xLowest;
	for (iref = 0, iNext = -1; iref<(int)vHitMCRef.size(); iref++) {
	  //                                            *****  LOOP OVER MC HITS
	  int bref = 1<<iref%32, jref = iref/32; if (hitsPats[jref]&bref)
	    continue;
	  int ihit = vHitMCRef[iref]; const THitMC &h = vHitMC[ihit];
	  if (h.IOrig!=0) {                   // ***** SKIP NON-ORIGINAL MC HITS
	    hitsPats[jref] |= bref; continue;
	  }
	  const TDetect &d =// Note: "Det2Ignore" skipped, cf. "TEv::GetMCInfo".
	    h.DetRef(); 
	  const TPlane &p = setup.Id2Plane(h.IDet);
	  int iZone; for (iZone = 0, zone = -1; iZone<nZones; iZone++) {
	    if (iplLasts[iZone]<p.IPlane) continue;
	    if (p.IPlane<=iplLasts[iZone]) { zone = iZone2Zone[iZone]; break; }
	  }
	  if (zone==-1) continue;

	  int isPixDet = 0; if (d.IType==29) isPixDet = 1; // GP P-pixel
	  else if ((TOpt::Graph[5]&0x2) && d.IType==32)
	    isPixDet = 2;                    // (upon option) MP M-pixel
	  if (DrOpt[4]!=-1 && p.IProj!=Proj &&
	      (isPixDet!=1  || !(pixZones&1<<zone)) &&
	      (isPixDet!=2  || p.IProj!=vProj)) {
	    hitsPats[jref] |= bref; continue;
	  }
	  if (iNext==-1 || d.X(0)<xLowest) {
	    xLowest = d.X(0); iNext = iref;
	  }
	}
	if (iNext>=0) {
	  hitsPats[iNext/32] |= 1<<iNext%32;
	  int ihit = vHitMCRef[iNext]; const THitMC &h = vHitMC[ihit];
	  const TDetect &d = h.DetRef();
	  const TPlane &p = setup.Id2Plane(h.IDet);
	  int iZone; for (iZone = 0, zone = -1; iZone<nZones; iZone++) {
	    if (iplLasts[iZone]<p.IPlane) continue;
	    if (p.IPlane<=iplLasts[iZone]) { zone = iZone2Zone[iZone]; break; }
	  }
	  if (zone==-1) {
  printf("TDisplay::DrawMCTracks: TKine #%d, TDetect #%d %s not in any zone\n",
	    ik,h.IDet,d.Name.c_str()); continue;
	  } 

	  xx[npt] = h.Xyz(0);
	  bool isPixDet = d.IType==29 || d.IType==32;
	  if (p.IProj==Proj ||                    // ***** CURRENT PROJECTION...
	      isPixDet && (pixZones&1<<zone)) {
	    if (p.IProj==Proj)  // Plane is in currently displayed proj. =>
	      // project MC MRS coord's onto detector plane axis.
	      uv = d.Ca*h.Xyz(1)+d.Sa*h.Xyz(2);
	    else                // Else, hence necessarily pixel, =>
	      // project MC MRS coord's onto proj. plane.
	      uv = cosA*h.Xyz(1)+sinA*h.Xyz(2);
	  }
	  else if (DrOpt[4]==-1) {         // ***** ...OR ALL PROJECTIONS OPTION
	    uv = cosA*h.Xyz(1)+sinA*h.Xyz(2);   // Project to current view plane
	  }
	  else continue;
	  yy[npt] = uv; 
	
	  if (++npt==2) {
	    if (ifl==0) {                            // ***** DRAW TRACK SEGMENT
	      ISLWSC(3); ISPLCI(color); ISPMCI(color); IPLm(2,xx,yy);
	    }
	    else {                                         // ***** LOCATOR MODE
	      if ((dist = TAlgo::Point2Line(px,py,xx,yy))<dist_min) {
		dist_min = dist; id_picked = ik;
	      }
	    }
	    // redraw hits
	    //IPMm(xx[0],yy[0]); IPMm(xx[1],yy[1]); 
	    if (xx[1]-xx[0]!=0) tga = (yy[1]-yy[0])/(xx[1]-xx[0]);
	    xx[0] = xx[1]; yy[0] = yy[1]; npt = 1;
	  }
	}
      } while (iNext>=0); // End of loop over hits
      if (ifl==0 && // In drawing mode and...
	  npt>1)    // ...current track is w/in currently displayed canvas
	DrawString(s,xx[0],yy[0],color,tga);
    }
    else {           // ********** PROPAGATION TRAJECTORY **********

      if (k.IVtxOrig < 0) continue;
      bool beam = false; if (k.P(0)<0) beam = true;
      const TVtxMC &v = ev.vVtxMC(k.IVtxOrig);
      THlx H0, H1, Hrot;
      H0(0) = v.V(0); H0(1) = v.V(1); H0(2) = v.V(2);
      H0(3) = k.P(1)/k.P(0); H0(4) = k.P(2)/k.P(0); H0(5) = k.Pinv();
      if (beam) H0(5)*=-1; // As beam goes backward in time in GEANT

      double xmax = -100000., xmin =  100000.;
      for (j = vHitMCRef.begin(); j!=vHitMCRef.end(); j++) {
	//                     ********** LOOP OVER MC HITS **********
	const THitMC &h =vHitMC[*j];
	if (h.IOrig!=0) continue;                // ***** SKIP !ORIGINAL MC HITS
	if (h.Xyz(0)>xmax) xmax = h.Xyz(0);
	if (h.Xyz(0)<xmin) xmin = h.Xyz(0);
      }
      if (beam) {
	if (H0(0)<Xvf[0]) continue;
	if (xmin >Xvf[1]) { // This is because of MC hits...
	  xmin = Xvf[0]+1;  // ...let's still draw the track. 
	}
      }
      else {
	if(H0(0)>Xvf[1]) continue;
	if(xmax <Xvf[0]) continue;
      }
      H0*=0; // Reset cov matrix for faster extrapolation
      double xl,xr,x1;
      if (beam) {
	xl = xmin  > double(Xvf[0]) ? xmin : double(Xvf[0]);
	xr = H0(0) < double(Xvf[1]) ? H0(0): double(Xvf[1]);
	if(xmin ==  100000.) xl = Xvf[0];
	if(H0(0) > Xvf[1]){ 
	  H1(0)=xr;  
	  if(!H0.Extrapolate(H1,false)) continue; 
	  H0=H1;
	}
      }
      else {
	xl = H0(0) > double(Xvf[0]) ? H0(0) : double(Xvf[0]);
	xr = xmax  < double(Xvf[1]) ? xmax  : double(Xvf[1]);
	if(xmax == -100000.) xr = Xvf[1];
	if(H0(0) < Xvf[0]){ 
	  H1(0)=xl;  
	  if(!H0.Extrapolate(H1,false)) continue; 
	  H0=H1;
	}
      }

      double step = (Xvf[1]-Xvf[0])/300.;

      const int nstep = 1000;
      float xx[nstep],yy[nstep];
      
      H0.Rotate(cosA,sinA,Hrot);
      int istep=0;
      xx[istep]=Hrot(0); yy[istep]=Hrot(1);
      while(1){
	if(beam){
	  x1=H0(0)-step; x1=x1 > xl ? x1 : xl;
	} else {
	  x1=H0(0)+step; x1=x1 < xr ? x1 : xr;
	}
	H1(0)=x1;
	if( !H0.Extrapolate(H1,false)) break;
	H1.Rotate(cosA,sinA,Hrot);
	if(istep+1 >= nstep) {
	  cout<<"TDisplay::DrawMCTracks ==> Wrong number of steps: Xl = "
	      << xl << "  Xr = " << xr << " X1 = " << x1 << "  step = " << step << endl;
	  break;
	}
	if( (Xvf[0] < Hrot(0) && Hrot(0) < Xvf[1]) && (Yvf[0] < Hrot(1) && Hrot(1) < Yvf[1])){
	  xx[++istep]=Hrot(0); yy[istep]=Hrot(1); // store point
	}
	if(beam) {
	  if(x1==xl) break;
	} else {
	  if(x1==xr) break;
	}
	H0=H1;
      }
      if(ifl == 0) {
	ISLN(1);
       ISMK(1);
       ISPLCI(color);ISPMCI(color);
       IPLm(istep+1,xx,yy); // draw track ( long polyline)
      }
      else {
	for (int l = 0; l<istep; l+=2) {
	 x[0]=xx[l]; x[1]=xx[l+1]; y[0]=yy[l]; y[1]=yy[l+1];
	 if ((dist = TAlgo::Point2Line(px,py,x,y))<dist_min){
	   dist_min=dist; id_picked = ik;
	 }
       }
     }

     if(istep > 1 && (xx[istep]-xx[istep-1]) != 0) 
       tga = (yy[istep]-yy[istep-1])/(xx[istep]-xx[istep-1]);
     DrawString(s, xx[istep],yy[istep], color, tga);
  
     
    }
  }     // End of loop over Kine tracks
  if (ifl==0) {
    ISPLCI(1); IUWK(0,1);
  }
  return id_picked;
}
