// $Id: TDisplayDrawTracks.cc 14080 2015-10-26 14:03:47Z lsilva $

/*!
  \brief Draw reco'ed tracks (TTrack objects) on a projection plane in the event display.
  \arg \a ifl: if !=0, no actual drawing, wait instead for a mouse click locating a track and return its index ("TTrack::Id").

  - 3 modes of play:
  i)   Polyline: through only (by default) those hits pertaining to proj. plane.
  ii)  Smooth: through all hits, via projections of the smoothed (at hit's abscissa) track.
  iii) Extrap.: extrapolating from helix at first (most upstream) measured point.
  - In Polyline mode, P-pixel hits are taken into account in all of 4 proj.: X, Y and 2 extras, defined by "TraF" options.
  - In Polyline mode but only upon option, M-pixels are taken into account in one extra proj., viz. their v-axis.
  - Upon option, extrapolate tracks to Z of pVertex.
  - MC data: colour code (gold/red=perfect/valid MC-association,brown=ghost), label w/ MC track index in addition to reco'd track index.
  - Get the angle of the current event display proj. from "TSetup::vecProj".
  - Upon option, delete tracks and display instead the hinging points for the bridging over the muFilter.
*/

//  The working of this method is conditioned by DrOpt options
// - DrOpt[5]: -1 = Polyline, 0 = Smooth (=polyline through ''smoothed'' points,
//             +1 = Extrapolation
// - DrOpt[13]: >0 = Propagate to vertex (conditioned in turn by DrOpt[5]<=0).

#include <cfloat>
#include <iostream>

#include "CsEvent.h"
#include "CsTrack.h"
#include "CsVertex.h"
#include "TDisplay.h"
#include "TEv.h"
#include "TOpt.h"
#include "THit.h"
#include "TPlane.h"
#include "TConstants.h"
#include "TAlgo.h"
#include "higz.h"

using namespace std;

/*
 iii) 2nd KINE track index.
*/

int TDisplay::DrawTracks(int ifl)
{  
  if (DrOpt[10]<0) {        // ***** TRACK DRAWING IS OFF *****
    if (ifl!= 0)
      cout<<"\nTDisplay::DrawTracks ==> Track can not be picked up, as drawing is OFF.\n"; 
    return -1;
  }

  TEv &ev = TEv::Ref(); const TSetup &setup = TSetup::Ref();

  if (TOpt::Graph[2]==-1) { // Special setting of Graph[2]

    //  **********************************************************************
    //  ********* DRAW the HINGING POINTS for BRIDGING OVER muFILTER *********
    //  **********************************************************************

 printf("\n****************************\n");
 printf("Special Display Mode: DRAWING HINGING POINTS for BRIDGING OVER WALL\n");
 printf("  => All reco'd tracks are erased. Tracking is upset beyond repair\n");
 printf("****************************\n");


    //  The bridging over the muFilter (in order to connect 0x4 and 0x8 track
    // segments) is based on hinging points, at which 0x4 and 0x8 extrapolations
    // are required to meet. (Note: The definition of the muFilter meant here is
    // an extended one, including H/ECAL2). The simple hinging point connection
    // is followed by a more rigorous fitting, using the material maps.)
    //  The hinging points are derived from the material maps. There are 3,
    // corresponding to 3 different regions defined in the yz plane. For each,
    // an helix is extrapoled through the material maps.
    //  The present special display is meant to give a graphical account of this
    // derivation: the 3 helices used draw as TTrack's (w/ a ``smooth'' point
    // defined @ the hinging point).
    static bool first = true; if (first) { // Else TTrack's already exist. 
      first = false;
      list<TTrack> &lTrack = const_cast<list<TTrack>&>(ev.lTrack());
      lTrack.clear(); // So that there is no confusion, neither on the display, nor in the mind of the user who would have triggered this special mode by chance.

      THlx Hs[12]; // 12 THlx, for 4 TTrack's w/ each Hfirst, Hlast, Hsmooth
      setup.ComputeHinges(0,Hs);

      int nRegions = TOpt::Calo[30]>1 /* !=(double)0 => ECAL2 exists */ ? 4 : 3;
      for (int region = 0; region<nRegions; region++) {
	TTrack t; t.Type = 0xc;
	t.Hfirst = Hs[3*region]; t.Hlast = Hs[3*region+1];
	lTrack.push_back(t);
	printf("Creating dummy TTrack #%d in region where mu-Absober consists in ",
	       lTrack.back().Id);
	switch (region) {
	case 0: printf("ECAL2\n"); break;
	case 1: printf("ECAL2+HCAL2\n"); break;
	case 2: printf("ECAL2+HCAL2+WALL\n"); break;
	case 3: printf("HCAL2+WALL\n");
	}
	// The 3rd helix output by TSetup::ComputeHinges is @ the hinging point.
	// Make it a smooth helix of the dummy track, so that a dot is drawn
	// at that abscissa on the event display.
	lTrack.back().lHsm.push_back(Hs[3*region+2]);
	// Enable printing smoothed helices in TTrack::Print, so that the
	// numerical value of the abscissa of the hinging point is also printed.
	TOpt::Print[0] |= 0x8;
      }
    }
  }

  float px, py; if (ifl!=0) {
    int nt, istat;
    IRQLC(1,20,istat,nt,px,py);
    if (istat==0) return -1;
  }
  else {
    px=py = 0;
  }

  if (ev.lTrack().empty()) return -1;

  CsVertex *bestVertex = 0;   //     ***** EXTRAPOLATE TO VERTEX ? *****
  if (DrOpt[5]>0 && DrOpt[13]>0) {
    const list<CsVertex*> &vertices = ev.ptrEvt()->getVertices();
    list<CsVertex*>::const_iterator iv;
    for (iv = vertices.begin(); iv!=vertices.end(); iv++) {
      CsVertex *v = *iv; if (!v->isBestVertex()) continue;
      if (!v->isPrimary()) {// Case should have been rejected supra
 printf("\n\a*** Inconsistency: Best vertex (Z = %.2f cm) not a primary\a\n",
	v->getZ()); continue;
      }
      if (bestVertex) { // Best vertex should be unique
 printf("\n\a*** Inconsistency: Event w/ 2 best vertices: Z = %.2f,%.2f cm\a\n",
	bestVertex->getZ(),v->getZ()); continue;
      }
      bestVertex = v;
    }
  }

  list<TTrack>::reverse_iterator it;
  list<TTrack> &lTrack = const_cast<list<TTrack>&>(ev.lTrack());
  float dist, dist_min; int id_picked;
  for (it = lTrack.rbegin(), dist_min = FLT_MAX, id_picked = -1;
       it !=lTrack.rend(); it++) {

    // ********************************************************************
    // ************************* LOOP OVER TRACKS *************************
    // ********************************************************************

    const TTrack &t = *it;
    if ((int)t.lHitPat.size()<TOpt::Graph[2])
      continue;         // (Optionally) Limit on number of hits


    int color(0);     // ***** TRACK COLOR AND POLYMARKER STYLE *****
    if (ev.IsMC()) {
      if (t.IKine==-1) {                        // Not associated to MC
	color=111; ISLWSC(1);
      }
      else {                                    // Associated to MC
	ISLWSC(3);
	if (t.NHits==t.NHsame) color = 112;       // Golden track
	else                   color = 2;
      }
    }
    else {
      color = 2; ISLWSC(3);                     // Real data tracks are red
    }
    ISMK(20);

    if (DrOpt[5]<=0 || t.Hfirst.empty()) {
      // ***** POLYLINE: IF REQUESTED OR IF THE TRACK WAS NOT FITTED *****

      bool smflg = // Can track be drawn as polyline through "smoothed" points?
	t.mHef.size() && t.mHeb.size() && t.mHuf.size() && t.mHub.size();
      bool smooth = // Draw "smoothed" polyline if it's required in addition
	DrOpt[5]==0 && smflg;

      if (!DrawPolyLine(smooth,color,ifl,px,py,dist,t)) continue;
    }
    else { // ***** TRACK'S HELIX PROPAGATION TRAJECTORY ***** 

      if (!DrawExtrapolate(bestVertex,color,ifl,px,py,dist,t)) continue;
    } // End drawing mode if/else 

    if (ifl &&    //             ***** LOCATING MODE *****
	dist>=0 && // Granted when Draw(PolyLine|Extraplate) returns true, yet, to be on the safe side...
	dist<dist_min) {
      dist_min = dist; id_picked = t.Id;
    }

    // Put dots where track parameters had been measured
    {
      list<THlx>lTmp(t.lHsm);
      lTmp.push_back(t.Hfirst);
      lTmp.push_back(t.Hlast);
      list<THlx>::const_iterator ihsm;

      THlx  H, Hrot;
      for (ihsm = lTmp.begin(); ihsm!=lTmp.end(); ihsm++) {
	H = (*ihsm);
	double ca,sa;
	double a = 0.1*setup.vProj()[Proj]*M_PI/180.; // approximate current proj. angle 
	ca=cos(a), sa=sin(a);
	H.Rotate(ca,sa,Hrot);
	ISPMCI(color);
	ISMK(25);
	IPMm(Hrot(0),Hrot(1));
	ISMK(20);
      } // End of loop over smoothed track pars.
    }
  } // End of loop over tracks

  if (ifl==0) {
    ISPLCI(1);
    IUWK(0,1); ISLWSC(1);
  }

  return id_picked;

}
// ===========================================================================
// =========================       DrawPolyLine      =========================
// ===========================================================================
bool TDisplay::DrawPolyLine(bool smooth, int color,
			    int ifl, float px, float py, float &dist,
			    const TTrack &t)
{
  TEv &ev = TEv::Ref(); const TSetup &setup = TSetup::Ref();

  //        ********** ANGLE of CURRENTLY DISPLAYED PROJECTION **********
  int iAngle = (int)floor(setup.vProj()[Proj]/10.+.5);
  double Angle = iAngle*M_PI/180, cosA = cos(Angle), sinA = sin(Angle);

  ISLN(1);

  float x[2], y[2];

  int npts, nsgmnts, iZone, nZones = setup.vIplLast().size(); float tga(0);
  for (iZone=npts=nsgmnts = 0, dist = -1; iZone<nZones; iZone++) {
    int zone = iZone2Zone[iZone];
    //      *************** LOOP ON ZONES ***************

    //                 ***** PIXELGEM/MMs *****
    // Pixels could be projected onto any projection. In TraFDic's proj.
    // search, they are projected on 4 different proj's: 0 and 90 degrees,
    // and 2 others, which +/-45 deg. defaults can be overriden by option on
    // a per zone basis.
    double aPixProjs[] = {0,90,45,-45};
    if (TOpt::iCut[29]&1<<zone) { // Zones where special pixelGEM projection
      aPixProjs[2] = TOpt::dCut[82]; aPixProjs[3] = TOpt::dCut[83];
    }
    int drawPPixDet, i03, vProj, iproj;
    for (i03=drawPPixDet = 0; i03<=3; i03++)
      if (abs(floor(aPixProjs[i03]+.5)-iAngle)<=1) {
	drawPPixDet = 1; break;
      }
    for (iproj = 0, vProj = -1; iproj<(int)setup.vProj().size(); iproj++) {
      int jAngle = (int)floor(setup.vProj()[iproj]/10.+.5);
      if (abs(abs(iAngle-jAngle)-90)<=1) {
	vProj = iproj;
	break; // Breaking here prevents the -90 proj., which contains only BMS, from being retained.
      }
    }

    list<int>::const_iterator ipl, ihp;
    for (ipl = t.lPRef().begin(), ihp = t.lHPat().begin();
	 ipl!=t.lPRef().end(); ipl++, ihp++) {

      // ****************************************************************
      //  ****** LOOP OVER TRACK'S PLANE REFERENCES ZONE AFTER ZONE *****
      // ****************************************************************

      if (*ihp<=-2) continue; // Out of active area or  plane is off

      const TPlane &p = setup.vPlane(*ipl);
      if (iZone && *ipl<=iplLasts[iZone-1]) continue;
      if (iplLasts[iZone]<*ipl) break;
      const TDetect &d = setup.iPlane2Detect(*ipl);
	
      const THit *h = 0;
      if (smooth) {          // ***** POLYLINE THROUGH SMOOTHED TRACK PARAMETERS
	THlx  H, Hrot; t.GetSmoothed(H,(*ipl));
	double ca, sa;
	if (p.IProj==Proj && d.IType!=29) { // Plane is in current proj. =>
	  // drawing more appropriate to evaluate departure from hit
	  ca = d.Ca; sa = d.Sa;
	}
	else
	  ca = cosA, sa = sinA;
	H.Rotate(ca,sa,Hrot); x[npts] = Hrot(0); y[npts] = Hrot(1);
      }
      else {                                      // ***** POLYLINE THROUGH HITS
	int isPixDet = 0; if (d.IType==29) isPixDet = 1; // GP P-pixel
	else if ((TOpt::Graph[5]&0x2) && d.IType==32)
	  isPixDet = 2;                    // (upon option) MP M-pixel
	if (p.IProj!=Proj &&         // Skip all but current proj. except...
	    (isPixDet!=1 || !drawPPixDet) && // ...P-pixels
	    (isPixDet!=2 || p.IProj!=vProj)) // ...(upon option) M-pixels
	  continue;
	if (*ihp<0) continue;              // Missing hit
	h = &ev.vHit(*ihp);
	if (!isPixDet) y[npts] = h->U;                   // ***** HIT COORDINATE
	else {
	  double a = atan2(d.Sa,d.Ca)-Angle;
	  y[npts] = h->U*cos(a)-h->V*sin(a);
	}
	x[npts] = d.X(0);
      }

      if (++npts==2) {  //          ***** DRAW TRACK SEGMENT *****
	nsgmnts++;
	if (ifl==0) {                                // ***** TRACK DRAWING MODE
	  ISPLCI(color); ISPMCI(color); IPLm(2,x,y);
	  if (smooth){
	    ISMK(31); IPMm(x[0],y[0]); IPMm(x[1],y[1]); ISMK(20);
	  }
	  else { // Redraw hits (only in hit connection mode)
	    if (h && (h->IKine==t.IKsame && h->IKine>=0 ||
		      h->IKin2==t.IKsame && h->IKin2>=0)) {
	      IPMm(x[0],y[0]); IPMm(x[1],y[1]);
	    }
	  }
	}
	else {                                             // ***** LOCATOR MODE
	  float d = TAlgo::Point2Line(px,py,x,y);
	  if (dist<0 || d<dist) dist = d;
	}
	if (x[1]-x[0]!=0) tga = (y[1]-y[0])/(x[1]-x[0]);
	x[0] = x[1]; y[0] = y[1]; npts = 1;
      }
    } // End of loop over plane references
  } // End loop on zones
  if (nsgmnts) {
    if (ifl==0) {
      char c[10]; if (t.IKine<0) sprintf(c,"%u",t.Id);      // ***** TRACK LABEL
      else                       sprintf(c,"%u(%u)",t.Id,t.IKine);
      string s(c); DrawString(s,x[0],y[0],color,tga);// Draw labels at track end
    }
    return true;
  }
  return false;
}
// ===========================================================================
// =========================      DrawExtrapolate    =========================
// ===========================================================================
bool TDisplay::DrawExtrapolate(CsVertex *bestVertex, int color,
			       int ifl, float px, float py, float &dist,
			       const TTrack &t)
{
  TEv &ev = TEv::Ref(); const TSetup &setup = TSetup::Ref();

  float x[2], y[2];

  float tga(0);

  THlx Hfirst = t.Hfirst, Hlast  = t.Hlast;

  if (DrOpt[13]>0) {   //        ***** EXTRAPOLATE TO VERTEX *****
    //#define TDisplay_DEBUG_VERTEX
#if defined DEBUG && defined TDisplay_DEBUG_VERTEX   // ***** INTERACTIVE VERTEX
    // Vertex coordinates have to be entered manually via, e.g., gdb
    static double xVertex = 0; static int doVertex = 0;
    static int ids[16] = {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
    if (doVertex) {
      int iId, ok; for (iId=ok = 0; iId<16; iId++)
		     if (t.Id==ids[iId]) ok = 1;
      if (ok) {
	THlx Hextr;
	if (t.InGroup(4)) {
	  Hextr(0) = xVertex+10;
	  if (Hlast.Extrapolate(Hextr))  Hlast  = Hextr;
	}
	else {
	  Hextr(0) = xVertex-10;
	  if (Hfirst.Extrapolate(Hextr)) Hfirst = Hextr;
	}
      }
    }
#else               // ***** MC: (upon option) TAKE VERTEX from MC TRUTH BANK...
    if (ev.IsMC() && TOpt::Graph[9]) {
      if (Hfirst.with_mom() &&               // Track w/ momentum...
	  // ...let's restrict ourself to 0x1 or 0x10. Other tracks could
	  // indeed be associated to the vertex. Yet, in order not to clog
	  // the display w/ too much drawing...
	  (t.InGroup(0) || t.InGroup(4))) {
	THlx Hextr;
	Hextr(0) = ev.vVtxMC(0).V(0); // pVertex assumed to be 1st one
	if (t.InGroup(4)) {
	  if (Hlast.Extrapolate(Hextr))  Hlast  = Hextr;
	}
	else {
	  if (Hfirst.Extrapolate(Hextr)) Hfirst = Hextr;
	}
      }
    }
    else {                        // ***** ...else: TAKE VERTEX from CsEvent
      if (bestVertex) {// If best pV exists: does current track belong to it
	const list<CsTrack*> vTracks = bestVertex->getTracks();
	list<CsTrack*>::const_iterator icst;
	bool inVertex; for (icst = vTracks.begin(), inVertex = false;
			    icst!=vTracks.end(); icst++) {
	  const CsTrack *cst = *icst;
	  map<int, int, less<int> >::const_iterator im =
	    ev.mapCsTrId2TrId.find(cst->getId());
	  if (im==ev.mapCsTrId2TrId.end())
	    CsErrLog::msg(elFatal,__FILE__,__LINE__,
  "CsTrack %d has no counterpart in TraFFic's event object",cst->getId());
	    if (im->second==(int)t.Id) { inVertex = true; break; }
	}
	if (inVertex) {                   // ***** CURRENT TRACK in BEST pVERTEX
	  //  Draw its extrapolation to Z of pVertex, using helix derived from
	  // tracking (as opposed to vertexing). Which is not always the most
	  // interesting drawing...
	  THlx Hextr; Hextr(0) = bestVertex->getZ()/10;
	  if (Hlast(0)<Hextr(0) || Hextr(0)<Hfirst(0)) { // Standard case...
	    if (Hlast(0)<Hextr(0)) {
	      if (Hlast.Extrapolate(Hextr))  Hlast  = Hextr;
	    }
	    else {
	      if (Hfirst.Extrapolate(Hextr)) Hfirst = Hextr;
	    }
	  }
	  // ...else: A track can overlap a vertex and yet belong to it,
	  // provided its last (resp. first) point is not too much downstream (
	  // resp. upstream) of it (tolerance is typically 20 cm)
	  // => In such a case: no extrapolation.
	}
      }
    }
#endif
  }

  if (Hfirst(0)>Xvf[1] || Hlast(0)<Xvf[0]) return false;

  //        ********** ANGLE of CURRENTLY DISPLAYED PROJECTION **********
  int iAngle = (int)floor(setup.vProj()[Proj]/10.+.5);
  double Angle = iAngle*M_PI/180, cosA = cos(Angle), sinA = sin(Angle);

  THlx H0 = Hfirst, H1, Hrot;
  H0*=0; // reset cov matrix for faster extrapolation
  double xl = Hfirst(0) > double(Xvf[0]) ? Hfirst(0) : double(Xvf[0]);
  double xr = Hlast (0) < double(Xvf[1]) ? Hlast (0) : double(Xvf[1]);
  double x1;

  if(Xvf[0] > H0(0)){ H1(0)=xl;  H0.Extrapolate(H1,true); H0=H1;}

  double step = (Xvf[1]-Xvf[0])/300.;
  const int NSTEP_MAX = 1000;
  float xx[NSTEP_MAX],yy[NSTEP_MAX];
      
  H0.Rotate(cosA,sinA,Hrot);
  int nstep=0;
  xx[nstep]=Hrot(0); yy[nstep]=Hrot(1);
  bool first = true;
  while(1){
    x1=H0(0)+step; x1=x1 < xr ? x1 : xr;
    H1(0)=x1;
    H0.Extrapolate(H1,true);
    H1.Rotate(cosA,sinA,Hrot);
    if (nstep+1>=NSTEP_MAX) {
      cout<<"TDisplay::DrawTracks ==> Wrong number of steps: Xl = "
	  << xl << "  Xr = " << xr << " X1 = " << x1 << "  step = " << step << endl;
      break;
    }
    if( (Xvf[0] < Hrot(0) && Hrot(0) < Xvf[1]) && (Yvf[0] < Hrot(1) && Hrot(1) < Yvf[1])){
      xx[++nstep]=Hrot(0); yy[nstep]=Hrot(1); // store point
    }
    if(x1==xr) break;
    H0=H1;
  } // end of while() loop

  if (!nstep) return false;

  if (ifl==0) {
    ISLN(1);
    ISMK(1);
    ISPLCI(color); ISPMCI(color);
    IPLm(nstep+1,xx,yy); // draw track ( long polyline)
    // Draw labels at track ends
    if(nstep > 1 && (xx[nstep]-xx[nstep-1]) != 0) 
      tga = (yy[nstep]-yy[nstep-1])/(xx[nstep]-xx[nstep-1]);
    char c[10]; if (t.IKine<0) sprintf(c,"%u",t.Id);    // ***** TRACK LABEL
    else                       sprintf(c,"%u(%u)",t.Id,t.IKine);
    string s(c); DrawString(s,xx[nstep],yy[nstep],color,tga);

  }
  else {
    dist = -1;
    for (int k = 0; k <nstep; k += 2) {
      x[0]=xx[k]; x[1]=xx[k+1];
      y[0]=yy[k]; y[1]=yy[k+1];
      float d = TAlgo::Point2Line(px,py,x,y);
      if (dist<0 || d<dist) dist = d;
    }
  }
  return true;
}
