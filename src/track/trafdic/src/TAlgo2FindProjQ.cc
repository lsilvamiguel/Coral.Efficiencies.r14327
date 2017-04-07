// $Id: TAlgo2FindProjQ.cc 14094 2015-11-06 15:28:48Z lsilva $

/*!
 Preliminary Pattern Recognition 
 in one projection using "pivot planes" method

 Input:
  igr       - Detectors' group, i.e. zone or set of zones.
             if   igr<5 igr == zone#
	     else special case of paraxial tracks
  X0        - Reference X coordinate

  npl       - Number of planes
  xpl [ipl] - X abscissa
  tols[ipl] - Residual tolerance
  fh  [ipl] - First hit index (in yhit[])
  lh  [ipl] - After-the-last hit index (in yhit[])
              If there are no hits on the plane, fh[ipl]==lh[ipl]
  idpl[ipl] - Type

  nhits     - Total number of hits
  yhit[ih]  - Hit coordinate in detector reference system (DRS)

  mode      - Mode of processing: expected to be iteration# in
             the iterative pattern finding algorithm.

 In/output
  Ntk       - In:  Max possible number of tracks (just for checking))
              Out: Number of found tracks
 Output:
  hits[i]   - Array of hit references (used in PrePattern only for debugging):
             [ih0,ih1... ihn,-1,ih0,ih1...ihn,-1,.....]
              \            /    \           /
               Hits of tr.1      Hits of tr.2  ...

  Y0_trk[i] - Track intercept at X0
  Yp_trk[i] - Track slope
  Y2_trk[1] - Quadratic term
*/

/*
  (Version specific to "lattice" alternative.
  Changes w/ respect to "traffic/TAlgoFindProj.cc":
  i) MACRO's "CURVATURE_TOLERANCE", "DRIFTS_SUBDUED".
  ii) argument "mode", used to tune option parameters...
  ...for the maximum slope to be considered,
  ...for the number of hits to qualify a track,
  depending upon iteration number in the iterative pattern finding algorithm.)
  iii) MACRO "QUADRATIC_PROJ"
  iv) Allows first, "igr", argument to stand for special groups (of
  detectors) used for paraxial track finding.
*/

#include <iostream>
#include <cstdio>
#include <cmath>
#include "cfortran.h"
#include "CsErrLog.h"
#include "TOpt.h"
#include "TSetup.h"
#include "TAlgo2.h"

PROTOCCALLSFSUB6(SORTZV,sortzv,PFLOAT,PINT,INT,INT,INT,INT)
#define SORTZV(A1,A2,A3,A4,A5,A6)  CCALLSFSUB6(SORTZV,sortzv,PFLOAT,PINT,INT,INT,INT,INT,A1,A2,A3,A4,A5,A6)

int TAlgo2::FindProjQ(int igr,float X0,
		      int npl,float xpl[],float tols[],int fh[],int lh[], 
		      int nhits,float yhit[],
		      int &Ntk,int hits[],float Y0_trk[],float Yp_trk[],
		      int mode,float *Y2_trk,int *idpl,int *href)
{
  int zone; // ***** DETECTORS GROUP -> SPECTROMETER ZONE *****
  if      (igr<5)  zone = igr;
  else {
    switch (igr) {
    case 5:  zone = -1; break; // Special case of paraxial PR in Z proj.
    case 8:  zone = 5; break;  // Zones 0x4 and 0x8 at once
#if defined PR_LOW_P_FORWARD || PR_2ZONE_ZPROJ
    case 9:
    case 10: zone = -1; break; // Special case of PR in Z spanning 2 zones
#endif
    default: zone = igr-6;
    }
  }

  const int maxcomb = 10000; // ***** ATTRIBUTES of COMBINATIONS *****
  static float Ysl [maxcomb], Y0[maxcomb], Nh[maxcomb];
  static int  ind  [maxcomb];
  //#define CURVATURE_TOLERANCE
#ifdef CURVATURE_TOLERANCE
  static float x0s [maxcomb], tol0s[maxcomb], dtoldxs[maxcomb], dxs[maxcomb];
  float tol0, dtoldx;
#endif

#if defined DRIFTS_SUBDUED || defined QUADRATIC_PROJ
  // Min distance between two drift planes that are not associated: the idea is
  // to give hits less weight when they are associated in a pair from
  // staggered drifts planes. (For if one falls w/in route, so does the other.)
  // A lower limit is set by the MB case, where interplane distance is ~2.8 cm.
  // The upper limit is not much constrained.
  static float deltax = 3.0;
#endif
 
  // ***** HIT FLAGS
  short int ifl[nhits]; for (int i = 0; i<nhits; i++) ifl[i]=0;

  //         ********** REQUIREMENTS on PR PARAMETERS **********
  float nHitsMn =                               // ***** MINIMUM # of HITS *****
    mode>0 ? TOpt::iPRpar[igr*10+3] : TOpt::iPRpar[igr*10];
  if (igr==3 &&                       // Special case of group #3 when...
      nhits<TOpt::iCut[24] &&           // ...few hits and...
      mode<=1) nHitsMn = 2;             // ...routes are at their narrowest
  if (zone==0 && fabs(nHitsMn-3)<.01 && // Special case of zone #0...in order to
      // ...catch MM01+DC01 tracks or decays downstream of MM01.
      // In order to limit #combinations, expected to be high in that case,...
      mode>=0) // ...exclude the case of <0 iterations
    nHitsMn -= .50;
#ifdef PR_LOW_P_FORWARD
  else if (igr==9 && fabs(nHitsMn-8)<.01) // Special case of low P forward PR
    // In case "nHitsMn" is required elsewhere (not the case as of 05/06) to
    // be as high as 8: relax somewhat the constraint:
    nHitsMn -= .50;    // E.g. 2 in zone 0x1 + 3 DCs or STs *1.75
#endif
  if (mode<0 && igr==1 && TOpt::iCut[31]) nHitsMn += 1;

  float angleMx;                         // ***** MAXIMUM ANGLE *****
  // Cure "dstar...2002.04...1" Evt#6 (has to be merged w/ case igr==2?)
  static float angleMxp;
  if      (igr==0) {
    angleMxp = TOpt::dPRpar[igr*10+1] * (1+mode/2.); angleMx = angleMxp*2.;
  }
  else if (igr==1)
    angleMx = TOpt::dPRpar[igr*10+1] * (1+mode/3.);
  else
    angleMx = TOpt::dPRpar[igr*10+1];  // If igr==2, it's updated infra.

  double dxTarget = X0 - TSetup::Ref().TargetCenter[0];
  double beamTAngle = TSetup::Ref().BeamTAngle;

  int ipl, nFound, nUsed, k, k0, kl, kr;
  double yp, y0, y1, dy, y, qual;
  double x, x0; float x1, dx; const float xSM1 = 350, xSM2 = 1800;
  int ipl0, ipl1;   // pivot planes
  int ncomb = 0;
  const int gem = 26;

  //        ***** VSAT SPECIFICS *****
  int scifi, stripGP, pixGP, pixMP;
  bool vSAT0x3Enhanced = 1<<(mode>=0 ? mode : 0)&TOpt::iCut[28];
  if (vSAT0x3Enhanced ||// Condition VSAT enhanced w/ iteration#("mode")..
      // ...This in case of groups 0x3, and hence excluding, in particular:
      // - Special group "igr"==5.
      // - All groups "igr>5", even when "zone<=1". For they do not go
      //  through an iterative PR.
      1<igr && igr<10) {
    scifi = 22; stripGP = 28; pixGP = 29;
    // For the pixelated pieces of MPs, we have to foresee 2 cases: those setups
    // where all micromegas are MPs or where redundancy in the beam region is
    // provided by SIs and the 2012 DVCS one, where only 01X and 01Y are MPs and
    // the pixelated pieces of these are the only detectors in the beam region
    // apart from FI04. In the latter case, we can have at most 2 hits per
    // projection, and therefore these need be all (MP as well as SI) enhanced.
    if (TOpt::ReMode[48]) pixMP = 32;
    else                  pixMP = -1;
  }
  else {
    scifi=stripGP=pixGP=pixMP = -1;
  }
#ifdef PR_LOW_P_FORWARD
  if (igr==9)   // Special case of low P forward PR: there scifis play...
    scifi = -1; // ...the role of SAT (as opposed to VSAT) => no enhancement
#endif
  static bool scifisOnly[maxcomb];

  int iplPG = -1;    //        ***** ZONE 0x1 SPECIFICS: GP01X/Y *****
  // Good resolution => narrow routes. A route based on a SI pivot + something
  // else as the 2nd pivot, will seldom catch the 2nd, intervening, SI. In order
  // to make for this effect, tolerance about SIs has already been enlarged, cf.
  // "TEv::PrePattern2". This presumably solves most cases. For the case GP01X/Y
  // (w/ even worse resolution): we decide to base the route on the sole SIs.
  // Need for this to extend the track search beyond the 2nd pivot down to GP01X/Y => single the latter out.
  float xuPA04 = 0, xdPA05 = 0;// ***** ZONE 0x4 SPECIFICS... *****
  //                 ***** ... I) LARGE ANGLE TRACKS *****
  //                 ***** ...II) 2-HIT TRACKS       *****
  //  Theses special searches are restricted...
  //   ... In case (I), to candidates starting far upstream, hence upstream of
  //      PA04 (i.e. 2nd PA in zone 0x4).
  //   ... In case (II) to candidates ending in the 1st 3 PA stations, hence in
  //      or upstream of PA05  (i.e. 3rd PA in zone 0x4).
  // To implement these restrictions, we determine here "xuPA04" = upstream of
  // PA04 and "xdPA05" = downstrean of PA05, by scanning the list of detectors
  // in "idpl" and counting the # of encountered PAs. This is fine, except when
  // any of the PAs has been turned off software-wise, as e.g. in the case of a
  // coral job aimed at doing detector studies on PAs. If it is so, we decide to
  // cancel the special searches: "xu/dPA04/5" = 0.
  if (igr==0) {
    for (ipl = 0; ipl<npl; ipl++) if (idpl[ipl]==29) iplPG = ipl;
  }
  if (igr==2) {
    int nPAs; for (ipl=nPAs = 0; ipl<npl; ipl++) {
      int detType = idpl[ipl];
      if (detType!=1 && detType!=26/*GM*/ && detType!=28 && detType!=29/*GP*/)
	break;
      if (detType!=1) continue;
      if    (++nPAs==2) xuPA04 = xpl[ipl]-1;
      else if (nPAs==3) { xdPA05 = xpl[ipl]+1; break; }
    }
    if (nPAs<3) { xuPA04=xdPA05 = 0; } // If no 3 PAs in [SM2,HO03]: reset
  }

#ifdef PR_ENHANCE_CENTRAL_GEMs
  // Enhance GEM w/ central zone activated (cf. "TEv::PrePattern2").
  int iplCentralGEM = -1; float u0CentralGEM = Y0_trk[0]; if (igr==1) {
    for (ipl = 0; ipl<npl; ipl++) {
      if (!(idpl[ipl]&0x200)) continue;
      iplCentralGEM = ipl;
      idpl[ipl] %= 0x200;  // Get rid of the tag => plane later be ID'd as a GEM
      if (idpl[ipl]!=26)   // Check that it's indeed a GEM
	CsErrLog::msg(elFatal,__FILE__,__LINE__,
	  "Plane %d tagged as DZfree GEM has type %d != 26!",ipl,idpl[ipl]);
    }
  }
#endif

  for (ipl0 = 0; ipl0<npl-1; ipl0++) { // ***** LOOP ON 1ST PIVOTS *****
    if (fh[ipl0]==lh[ipl0]) continue;         // No hits
    int idpl0 = idpl[ipl0];
    bool scifi_0 = idpl0==scifi || idpl0==stripGP || idpl0==pixGP ||
      idpl0==pixMP;
    bool si_0 = idpl0==21;
    bool mb_0 = idpl0==18 && nhits<2*TOpt::iCut[24] && mode==0;

    x0 = xpl[ipl0]; float angMx;
    if (igr==2) {  // Special case of tracks bent to large angles by SM2
      angMx = angleMx;
      if (x0<xuPA04) angMx *= 2; if (mode>0)   angMx /= 2;
    }
    else angMx = angleMx;

    float dx_mn;                // ***** MINIMUM DISTANCE TO 2ND PIVOT *****
    if (12<=idpl0 && idpl0<=14 && // ...special case of outer|vertical straws...
	mode>=2)                  // ...when in a later iteration
      dx_mn = 5;     // ... relax condition to catch ST03X1:2 + ST03Y1:2
    else
      dx_mn = 15;

    for (ipl1 = npl-1; ipl1>ipl0; ipl1--) { // ***** LOOP ON 2ND PIVOTS *****
      if (fh[ipl1]==lh[ipl1]) continue;     // No hits
      int idpl1 = idpl[ipl1]; bool scifi_01 = scifi_0 &&
	(idpl1==scifi || idpl1==stripGP || idpl1==pixGP) &&
	igr!=5; // Special case of paraxial PR in Z proj.
      bool noneButScifi = scifi_01 && // If both pivots are scifis...
	// ...track is expected to travel in scifis throughout, except...
	igr!=4; // ...in scifi/Si telescope: this flag to speed up processing
      int iplLast = si_0 && idpl1==21 && iplPG!=-1 ? iplPG : ipl1;
      // DVCS 2012:
      bool tooCrowded = // DVCS 2012: Low redundancy in VSAT of zone 0x1
	// => Have then to retain all combinations of 1*MP +1*FI => Can yield
	// high number of tracks in the proj. search and then overflow the space
	// search, if MPs are noisy or not time-cut (because of lack of time
	// calibration). This situation is expected to be flagged by
	// "ReMode[48]", and is checked on the # of hits per plane.
	igr==0 && (lh[ipl1]-fh[ipl1])*(lh[ipl0]-fh[ipl0])>100 &&
	TOpt::ReMode[48]>1;

      x1 = xpl[ipl1]; dx = x1-x0;
      if (dx<dx_mn) continue;               // Less than minimum allowed spacing

      // Min # of hits in zone 0x4: 2 is found to improve efficiency but is also
      // comsuming very much CPU => restrict 2 to:
      // - tracks ending in the 1st 3 PA stations,
      // - tracks in, or close to, the beam region. 
      float nHsMn; if (igr==2 && x1>xdPA05 &&
		       !(x1>3000 && (idpl1==22 || idpl1==26 || idpl1==43))) {
	if (nHitsMn<2.1)   nHsMn = 3;
        else               nHsMn = nHitsMn;
      }
      else if (mb_0 && idpl1==18) 
      	nHsMn = 2;
      else
	nHsMn = nHitsMn;
      nHsMn -= .01;  // Since we are dealing with float max value

#ifdef CURVATURE_TOLERANCE
      tol0 = tols[ipl0]; dtoldx = (tols[ipl1]-tol0)/dx;
#endif
      for (int ih0 = fh[ipl0]; ih0<lh[ipl0]; ih0++){
	//                        ***** LOOP HITS in 1ST PIVOT  *****
	y0 = yhit[ih0];
	for(int ih1 = fh[ipl1]; ih1<lh[ipl1]; ih1++){
	  //                        ***** LOOP ON HITS in 2ND PIVOT *****
	  y1 = yhit[ih1]; yp = (y1-y0)/dx;
#ifdef PR_ENHANCE_CENTRAL_GEMs
	  if (scifi_0 &&                // Treat central part of GEM as a scifi.
	      fabs(y1-u0CentralGEM)<5.0 &&
	      ipl1==iplCentralGEM)      // Never happens if "igr==5" (cf. supra)
	    scifi_01 = true;
	  // Do not set "noneButScifi" for this "ih1" might still be outside DZ.
#endif

	  //                        ***** ANGLE CUT *****

	  float angle;
	  if (igr==4) angle = fabs(yp-beamTAngle);
	  else {
#define TARGET_POINTING_PROJ 2
#ifdef TARGET_POINTING_PROJ
#  if TARGET_POINTING_PROJ == 2
	    angle = fabs((y0 + yp*(X0-x0))/dxTarget-yp);
#  else
	    angle = fabs(yp);
#  endif
#else
	    angle = fabs(yp);
#endif
	  }
	  if (angle>angMx) continue;

	  // ***** FIND HITS in INTERMEDIATE PLANES for DEFINED ROUTE *****
	  float xPrv, nabscissa; int hrefPrv; float nDCSTs; int nSATs;
	  if (mode<0 && // Early iter. w/ tight PR requirements...
	      // ...in zones 0x3: let's disregard tracks in LAT, which are
	      // typically low momenta, hence much impacted by MS or field and
	      // therefore expected to deviate much from straight line.
	      zone<=1) nDCSTs = 0;
	  else nDCSTs = -1;
#if defined PR_LOW_P_FORWARD || defined PR_2ZONE_ZPROJ
	  // Special case of PR in Z spanning 2 zones: # of hits in zone 0x1
	  static int nZone0x1; int nZ0x1Mn = TOpt::iPRpar[10*igr+2];
#endif
#ifdef PR_2ZONE_ZPROJ
	  static double s1, sx, sx2, sy, sxy, sy2;
	  if (igr==10) {
	    s1=sx=sx2=sy=sxy=sy2 = 0;
	  }
#endif
	  for (ipl = ipl0, nFound=nSATs = 0, nabscissa=qual = 0,
		 xPrv = X0-1000, hrefPrv = -1; ipl<=iplLast; ipl++) {
	    if (fh[ipl]==lh[ipl]) continue;
	    int type = idpl[ipl];
	    if (noneButScifi &&
		type!=scifi && type!=stripGP && type!=pixGP)
	      //&& type!=pixMP: No! No pixelated MP in the 2nd pivot, given that
	      // if they are ever enhanced, it's in the 2012 DVCS context w/
	      // only MP01X/Y, which in turn only allows ipl0=MP+ipl1=FI proj. 
	      continue;

	    x = xpl[ipl]; y = y0+yp*(x-x0);

	    float tol = tols[ipl];
#ifdef CURVATURE_TOLERANCE
	    float xx0 = x-x0;
	    float dx_lim = xx0>dx ? dx : xx0;  // Limit the range of xx0/dx
	    float tolp = tol0+dtoldx*dx_lim;
	    tol = tolp>tol ? tolp : tol;
#endif
	    if (y<yhit[fh[ipl]]  -tol ||
		y>yhit[lh[ipl]-1]+tol) continue; // "y" is too far from hits on this plane  

	    // find hit, closest to track position "y"
	    kl = fh[ipl]; kr = lh[ipl]-1; while (kr-kl>1) {
	      k0 = (kl+kr)>>1;
	      if   (y<=yhit[k0]) kr=k0;
	      else               kl=k0;
	    }
	    if   (yhit[kr]-y<y-yhit[kl]) k0=kr;
	    else                         k0=kl;
	    if (fabs(y-yhit[k0])>tol)  continue;

#if defined PR_LOW_P_FORWARD || defined PR_2ZONE_ZPROJ
	    if (igr>=9) {         // Special case of PR in Z spanning 2 zones
	      if (x>xSM1) {
		if (nFound<nZ0x1Mn) break;     // Require min. #hits in zone 0x1
	      }
	      else nZone0x1 = nFound+1;        // # hits in zone 0x1
#  ifdef PR_2ZONE_ZPROJ
	      if (igr==10) {        // Special case of 2-zone proj. in Z
		float yi = yhit[k0], w2 = 1/tol/tol;
		s1 += w2; sx += x*w2; sx2 += x*x*w2;
		sy += yi*w2; sxy += x*yi*w2; sy2 += yi*yi*w2;
	      }
#  endif
	    }
#endif
	    nFound++;

	    float iabscissa;
	    if (igr==5 && // Special case of paraxial PR in Z proj.
		type==22) {
	      // Scifis only enhancement. So as to fulfill typical requirement
	      // of 6 w/, e.g. 4 out 5 FI03::07
	      if (scifi_0) iabscissa = 1.50;
	      else         iabscissa = 1.;
	    }
#ifdef DRIFTS_SUBDUED
	    else if (x<xPrv+deltax) {   // No twice dble-layered drifts
	      if (12<=type && type<=15)      // DCs, Outer/Vertical !ST04 straws
		iabscissa = .75;
	      else if (11==type)             // Inner XUV !ST04 straws
		iabscissa = .20;
	      else if (type==16 || type==18) // W45s, MBs
		iabscissa = .20;
	      else if (type==17)             // ST04, DRs
		iabscissa = .75;		
	      else    // Virtual XY GEM made out of UV by amplitude correlation
		iabscissa = 1.;
	    }
#endif
	    //#  define SOLICIT_GEMs
#ifdef SOLICIT_GEMs
	    else if (type==gem) iabscissa = .75;
#endif
	    else iabscissa = 1.;
	    if (// type==26 /* this is ensured by construction */ &&
		href[k0]==hrefPrv) iabscissa -= .50;

	    qual += iabscissa-fabs(y-yhit[k0])/tol; nabscissa += iabscissa;
	    if (nDCSTs>=0) {  // Special handling of LAT in early iteration
	      if (11<=type && type<=15) nDCSTs += iabscissa;
	      if (type==26 /* GM */ || type==27 /* MM */) nSATs++;
	    } 

	    hrefPrv = href[k0]; xPrv = x;
	  } // End of loop over planes 

	  //                                    ***** CUT on # of HITS *****
	  if (nabscissa<nHsMn && !scifi_01) continue;
	  if (nDCSTs>3 && // Special handling of LAT in early iteration
	      // We disregard combinations of quasi-exclusively DC/ST02/3... Yet
	      // allowing still for tracks going through the overlap w/ SAT.
	      nSATs<2) continue;

	  if (igr==0) {  // group 0: Stricter angle cut for not so good comb
	    if (angle>angleMxp &&
		(nabscissa<nHsMn+1 || // E.g. case 2MM+2DC(==3.75)>3.5 passes
		 qual/nFound<.75)) continue;
	    if (scifi_01 && nabscissa<nHsMn && tooCrowded) {
	      // => Let's then severly cut on angle.
	      if (angle>.02) continue;
	    }
	  }
	  else if (zone==1) {// Zone0x2, not scifis: extra hit to bridge RICH...
	    if (x1-x0>750 && !scifi_0 && nabscissa<nHsMn+1 &&
		TOpt::iCut[5]) continue; // ...upon option
	  }
	  else if (igr==5) { // Special case of paraxial PR in Z proj....
	    if      (x1<xSM1) continue;  // ...require several zones
	    else if (xSM1<x0 && x1<xSM2) continue;
	    else if (xSM2<x0) continue;
	  }

	  // ***** STORE COMBINATION *****
	  Ysl[ncomb] = yp; Y0 [ncomb] = y0 + yp*(X0-x0);
#ifdef PR_2ZONE_ZPROJ
	  if (igr>=10) {          // Special case of 2-zone proj. in Z
	    double alpha = (s1*sxy-sx*sy)/(s1*sx2-sx*sx), b = sy/s1-alpha*sx/s1;
	    float chi2 = (alpha*alpha*sx2+b*b*s1+sy2+
			  2*(alpha*b*sx-alpha*sxy-b*sy))/(nFound-2);
	    if (chi2>.25) continue;
	    Nh [ncomb] = nFound-nZone0x1-10*chi2;
	  }
	  else
#endif
	    Nh [ncomb] = nFound+qual/nFound;

	  if (scifi_01) {
#ifndef M_PI
#define M_PI   3.14159265358979323846  /* pi ! Compilation error w/out this! */
#endif
	    Nh[ncomb] += atan2(1,fabs(yp))/M_PI;
	    scifisOnly[ncomb] = true;
	  }
	  else {
#ifdef TARGET_POINTING_PROJ
	    Nh[ncomb] += .5-5*fabs(Y0[ncomb]/dxTarget-yp);	    
#endif
	    scifisOnly[ncomb] = false;
	  }
	  ind[ncomb] = ncomb+1;
#ifdef CURVATURE_TOLERANCE
	  x0s[ncomb] = x0; tol0s[ncomb] = tol0; dtoldxs[ncomb] = dtoldx;
	  dxs[ncomb] = dx;
#endif
	  if (++ncomb>=maxcomb) {
	    CsErrLog::msg(elError,__FILE__,__LINE__,	
			  "Too many combinations (%d) in group %d",ncomb,igr);
	    goto sorting;
	  }
	  
	} //end of loop over 2-d pivot plane hits 
      }//end of loop over 1-st pivot plane hits
 
    } // end of loop over 2-d pivot planes
  } // end of loop over 1-st pivot planes
  
  // sort Nh-wise

 sorting:
  SORTZV(Nh[0],ind[0],ncomb,1,1,0);

  //     ***** MAXIMUM # COMMON HITS (but for scifis tracks) *****
  int nUsedMx = mode>0 ? TOpt::iPRpar[igr*10+7] : TOpt::iPRpar[igr*10+6];
  int nUsedMX;  // Case of tracks w/ many hits
  if (mode>0)
    nUsedMX = TOpt::iPRpar[igr*10+9] ? TOpt::iPRpar[igr*10+9] : nUsedMx;
  else
    nUsedMX = TOpt::iPRpar[igr*10+8] ? TOpt::iPRpar[igr*10+8] : nUsedMx;


  int m; int nt=0, nh=0;
  int hh[npl];
  for(int l=0; l < ncomb; l++){ // loop over combinations starting from Nh max 
    m=ind[l]-1;

#ifdef CURVATURE_TOLERANCE
    x0 = x0s[m]; tol0 = tol0s[m]; dtoldx = dtoldxs[m]; dx = dxs[m];
#endif

    float s1 = 0, sx = 0, sx2 = 0, sy = 0, sxy = 0, sy2 = 0;
#ifdef QUADRATIC_PROJ
    static int   nabsc;                // "static" in order to ... 
    static float sx3, sx4, sx2y, xprv; // ...prevent compiler from warning
    static float xFirst, xLast;
    if (Y2_trk) { sx3=sx4=sx2y = 0; nabsc = 0; xprv = -1000; }
#endif
    bool scfOnly = scifisOnly[m];
    // Starting w/ 2008, there can be 3 VSAT planes per projection, e.g.
    // FI05+GP02P1+GP02P2 => "nUMx" = 3 in order to collect all combinations.
    int nUMx = scfOnly ? 3 : (Nh[m]>nHitsMn+1 ? nUsedMX : nUsedMx);
    for (ipl=nFound=nUsed = 0; ipl<npl; ipl++) { 
      // count hits on _all_ planes in defined direction
      if (fh[ipl]==lh[ipl]) continue;
#ifdef PR_ENHANCE_CENTRAL_GEMs
      if (scfOnly &&
	  idpl[ipl]!=scifi && idpl[ipl]!=stripGP && idpl[ipl]!=pixGP && 
	  idpl[ipl]!=pixMP && ipl!=iplCentralGEM) continue;
#else
      if (scfOnly &&
	  idpl[ipl]!=scifi && idpl[ipl]!=stripGP && idpl[ipl]!=pixGP
	  && idpl[ipl]!=pixMP) continue;
#endif

      float xX0 = xpl[ipl]-X0; y=Y0[m]+Ysl[m]*xX0;

      float tol = tols[ipl];
#ifdef CURVATURE_TOLERANCE
      float xx0 = xpl[ipl]-x0;
      float dx_lim = xx0>dx ? dx : xx0;   // Limit the range of xx0/dx
      float tolp = tol0+dtoldx*dx_lim;
      tol = tolp>tol ? tolp : tol;
#endif

      if (y<yhit[fh[ipl]]  -tol ||
	  y>yhit[lh[ipl]-1]+tol) continue; // "y" is too far from hits

      // find hit, closest to track position "y"
      kl=fh[ipl]; kr=lh[ipl]-1;
      while((kr-kl) > 1){
	k0=(kl+kr)>>1;
	if(y <= yhit[k0]) kr=k0;
	else              kl=k0;
      }
      if((yhit[kr]-y) < (y-yhit[kl])) k0=kr;
      else                            k0=kl;
      
      if (fabs(y-yhit[k0])<tol) {   // ***** W/IN TOLERANCE: STORE HIT
	hh[nFound++] = k0;
	if (ifl[k0]!=0) nUsed++; if (nUsed>nUMx) goto next_comb;
	float xi = xX0, yi = yhit[k0];
	s1 += 1; sx += xi; xi *= xX0; sx2 += xi;
	sy += yi; sxy += xX0*yi; sy2 += yi*yi;
#ifdef QUADRATIC_PROJ
	if (Y2_trk) {
	  sx2y += xi*yi; xi *= xX0; sx3 += xi; xi *= xX0; sx4 += xi;
	  if (xpl[ipl]>xprv+deltax) nabsc++;
	  if (nFound==1) xFirst = xX0; xLast = xX0;
	}
#endif 
      } 
    } // end of loop over planes 

    for (int i = 0; i<nFound; i++) {            // ***** STORE FOUND HITS *****
      hits[nh++] = hh[i]; ifl[hh[i]] = 1; // Flag hit as used
    }
    hits[nh++] = -1;   // End of hit list

    //store found track parameters
    //Y0_trk[nt]=Y0[m];
    //Yp_trk[nt]=Ysl[m];
#ifdef QUADRATIC_PROJ
    if (Y2_trk && nabsc>3 && xLast-xFirst>1000) {
      float B = sx3*s1-sx2*sx, C = sx2*s1-sx*sx, A = sx4*s1-sx2*sx2;
      float det = A*C-B*B;
      if (fabs(det)>1.e-3*A*A) {
	float E = sx2y*s1-sx2*sy, F = sxy*s1-sx*sy;
	float a = (E*C-B*F)/det, b = (A*F-B*E)/det;
	Y2_trk[nt] = a; Yp_trk[nt] = b;
	Y0_trk[nt] = (sy-a*sx2-b*sx)/s1;
      }
      else {
	Y2_trk[nt] = 0.; 
	float alpha = (s1*sxy-sx*sy)/(s1*sx2-sx*sx);
	Y0_trk[nt] = sy/s1-alpha*sx/s1;
	Yp_trk[nt] = alpha;
      }
    }
    else
#endif
    {
      float alpha = (s1*sxy-sx*sy)/(s1*sx2-sx*sx), b = sy/s1-alpha*sx/s1;
      Y0_trk[nt] = b; Yp_trk[nt] = alpha;
    }

    if(++nt >= Ntk)  {
      CsErrLog::msg(elError,__FILE__,__LINE__,
		    "Too many tracks (=%d) in group %d",nt,igr);
      return 1;
    }

  next_comb:;
  } // end of loop over combination

  if (TOpt::Print[5])
    printf("  N comb. = %5u   N found tracks = %4u \n", ncomb, nt);

  Ntk = nt;  return 0;
}
