//   $Id: TTrackQNewtonFit.cc 13534 2012-07-18 00:05:11Z ybedfer $

/*!

  Interface to Quasi-Newton == fit through already found hits, using either
track propagator or Dico and in that case possibly fast fit, depending upon
option flags.
  Same i/o as QuickKF, except that it does not determine covariance.

 \param dir  fit direction
 \param mode defines how to fit 

 \verbatim
 Input:

   dir =  1 - fit in downstream direction (resulting helix is Hlast )
   dir = -1 - fit in   upstream direction (NOT implemented!!)
 
   mode =  0 - Straight line fit
   mode =  1 - Starts from known (to 1st approx.) helix 
   mode =  3 - Starts from known (to 1st approx.) helix w/ 1/p == 1/pxz
   mode =  2 - Special fit w/, as much as possible, fixed momentum while
              restricting oneself to zones 0x3. It's been specially designed to
              meet the needs of "TEv::BackTrackSAS". If too few space points in
	      zones 0x3, the fit domain is extended to include zone 0x4 and the
	      momentum is freed.
   mode =  9   Special fit, also designed for "TEv::BackTrackSAS", where still
              restricting oneself to zones 0x3, if enough space points, but this
	      time, free momentum in any case. It is to be used after some 0x1
	      hits have already been picked up, when these hits are SIs, whose
	      good resolution is expected to be too demanding for a fixed
	      momentum QN fit.
   mode = 10 - Special fit, also designed for "TEv::BackTrackSAS", where the fit
              domain is restricted to hits upstream of RICH and momentum fixed.
   mode =  4 - Init values in Haux (as opposed to Hfirst)
   mode =  5 - Init values in Haux w/ 1/p == 1/pxz
   mode =  6 - Special fit w/ fixed Y vertex. It's specially designed for
              "TEv::BackTrackZ2". Init values in Haux.
   mode =  8 - Special fit designed for "TEv::BackTrackZ1", where the fit
              domain is restricted to hits upstream of RICH.
   mode =  7 - Init values in Haux and update all guests
\endverbatim

*/

#include <iostream>
#include <iomanip>
#include "CsErrLog.h"
#include "TOpt.h"
#include "TEv.h"
#include "TLattice.h"

using namespace std;

TLattice *Lat = NULL;

#ifdef EXTENSIVE_DICO
#  define HINI Hfirst
#  define HOUT Haux
#else
#  define HINI Hfirst
#  define HOUT Hfirst
#endif

bool TTrack::QNewtonFit(int dir, int mode)
{
#ifdef QNewton_ENABLED
  const TSetup &setup = TSetup::Ref();
  TEv &ev = TEv::Ref(); const vector<THit>& vHit = ev.vHit();
  
  if (dir!=1) CsErrLog::mes(elFatal,"Backward Fit not allowed!");
  if (lPlnRef.size()!=lHitPat.size())
    CsErrLog::mes(elFatal,"lPlnRef.size() != lHitPat.size()");

  if (!Lat) Lat = new TLattice(1);
  static float resols[NROWS]; Lat->sigmas = resols;

  Track  ft; Parameters *par = &ft.par;
  // # of DoF's + 1 when it comes to actual fit,
  //            + 0 when simply extrapolating
  static unsigned int nh0; // static to avoid warning "might be uninitialized.."
  static unsigned int nh1;

  // *************** FILL STARTING VALUES ***************

  bool all_guests = false;
  switch (mode) {
  case 0:                   // Straight line fit (no momentum)
    HOUT(1,1)=HOUT(2,2) = 50.*50.; HOUT(3,3)=HOUT(4,4) = 0.5*0.5;
    HOUT(5,5)=1.;                              // Big error on q/P
    par->pxzi = 0;                             // Set infinite momentum
    //par->x=par->y=par->tx=par->ty = 0.;
    ft.par_pattern = 0x1b; nh0 = 4; nh1 = 1;
    break;
  case 1:                   // Fit with momentum
    par->x = HINI(1)-HINI(0)*HINI(3); par->tx = HINI(3);
    par->y = HINI(2)-HINI(0)*HINI(4); par->ty = HINI(4);
    par->pxzi = HINI(5)*sqrt(1+par->ty*par->ty/(1+par->tx*par->tx));
    if (Lat->OutDico(par) &&    // Probably outside Dico on view of init pars..
	(TOpt::ReMode[14]&0x4)) // ...while using Dico...
      return false;             // ...=> return (and try KF...)
    HOUT *= 100.;             // Scale up diagonal cov. matrix
    ft.par_pattern = 0x1f; nh0 = 6; nh1 = 6;
    break;
  case 2:                   // Fixed momentum
  case 9:
    if (!HINI.with_mom())
      CsErrLog::mes(elFatal,"Fit with fixed momentum, but momentum not known");
    par->x=par->y=par->tx = 0;
    par->ty = HINI(4);
    par->pxzi = HINI(5)*sqrt(1+par->ty*par->ty/(1+par->tx*par->tx));
    HOUT(1,1)=HOUT(2,2) = 50.*50.;
    HOUT(3,3)=HOUT(4,4) = 0.5*0.5;
    par->pxzi = HINI(5);
    HOUT(5,5) = 1.;
    ft.par_pattern = 0x1b; nh0 = 4; nh1 = 1;
    all_guests = true;
    if (mode==9) ft.par_pattern = 0x1f;
    break;
  case 3:                   // Expecting HINI to contain 1/pxz, NOT 1/p
    par->x = HINI(1)-HINI(0)*HINI(3); par->tx = HINI(3);
    par->y = HINI(2)-HINI(0)*HINI(4); par->ty = HINI(4);
    par->pxzi = HINI(5);
    if (Lat->OutDico(par) &&    // Probably outside Dico on view of init pars..
	(TOpt::ReMode[14]&0x4)) // ...while using Dico...
      return false;             // ...=> return (and try KF...)
    HOUT *= 100.;             // Scale up diagonal cov. matrix
    ft.par_pattern = 0x1f; nh0 = 6; nh1 = 6;
    break;
#  ifdef EXTENSIVE_DICO
  case 7:
    all_guests = true;
  case 4:                   // Init values (w/ 1/p) in Haux
    par->x = Haux(1)-Haux(0)*Haux(3); par->tx = Haux(3);
    par->y = Haux(2)-Haux(0)*Haux(4); par->ty = Haux(4);
    par->pxzi = Haux(5)*sqrt(1+par->ty*par->ty/(1+par->tx*par->tx));
    if (Lat->FarOutDico(par) && // Far outside Dico on view of init pars..
	(TOpt::ReMode[14]&0x4)) // ...while using Dico...
      return false;             // ...=> return (and try KF...)
    HOUT *= 100.;             // Scale up diagonal cov. matrix
    ft.par_pattern = 0x1f; nh0 = 6; nh1 = 6;
    break;
  case 5:                   // Init values (w/ 1/pxz) in Haux
    par->x = Haux(1)-Haux(0)*Haux(3); par->tx = Haux(3);
    par->y = Haux(2)-Haux(0)*Haux(4); par->ty = Haux(4);
    par->pxzi = Haux(5);
    if (Lat->OutDico(par) &&    // Probably outside Dico on view of init pars..
	(TOpt::ReMode[14]&0x4)) // ...while using Dico...
      return false;             // ...=> return (and try KF...)
    HOUT *= 100.;             // Scale up diagonal cov. matrix
    ft.par_pattern = 0x1f; nh0 = 6; nh1 = 6;
    break;
  case 8:
  case 6:     // Fixed Y vertex: Init values (w/ 1/p) in Haux??
    par->x = Haux(1)-Haux(0)*Haux(3); par->tx = Haux(3);
    par->y = Haux(2)-Haux(0)*Haux(4); par->ty = Haux(4);
    par->pxzi = Haux(5);
    if (Lat->FarOutDico(par) && // Far outside Dico on view of init pars..
	(TOpt::ReMode[14]&0x4)) // ...while using Dico...
      return false;             // ...=> return (and try KF...)
    HOUT *= 100.;             // Scale up diagonal cov. matrix
    ft.par_pattern = 0x1e; nh0 = 3 /* case FI07XY+FI08XY */; nh1 = 1;
    all_guests = true;
    break;
  case 10:    // Fixed momentum: w/ 1/px in Haux;
    par->x = HINI(1)-HINI(0)*HINI(3); par->tx = HINI(3);
    par->y = HINI(2)-HINI(0)*HINI(4); par->ty = HINI(4);
    par->pxzi = HINI(5)*sqrt(1+par->ty*par->ty/(1+par->tx*par->tx));
    if (Lat->FarOutDico(par) && // Far outside Dico on view of init pars..
	(TOpt::ReMode[14]&0x4)) // ...while using Dico...
      return false;             // ...=> return (and try KF...)
    HOUT *= 100.;             // Scale up diagonal cov. matrix
    ft.par_pattern = 0x1b; nh0 = 3 /* case FI05XY+GM03XY|UV */; nh1 = 1;
    all_guests = true;
    break;
#  endif
  default:
    CsErrLog::mes(elFatal,"Unexpected mode switch value."); 
  }
  
  // *************** BUILD STRUCTURE "Track" ***************

  double chi2max = 0;
  Track_Hits *hits = &ft.hits; memset(&hits->hp,0,sizeof(hits->hp));
  list<int>::iterator ih = lHitPat.begin();

  if (mode==2 || mode==9) { // ***** SPECIAL CASE OF FIT w/ FIXED MOMENTUM *****
    // This is for backtracking SM2 track in view of extending it upstream of
    // SM1.
    // - Since it seems that the full fit, including the part downstream of
    //  SM2 is difficult
    //       => TRY AND STRIP THE ZONE #2 TAIL.
    //  Since the momentum is already defined we do not need it. Except if track
    //  is so poor in hits in zone #1 that it does not even have 2
    //  ``space points'', e.g. track w/ FI05XorY + FI06XYV
    // - Never mind redefinig "NDics". Probably not useful anymore.
    int nSpacePts = 0, projs = 0, nProjs = 0; const TStation *sPrv = 0;
    bool firstZ2 = true; NDics = 0; while (ih!=lHitPat.end()) {
      if (*ih>=0) {
	const THit &h = vHit[*ih]; int ipl = h.IPlane; 
	const TPlane &p = setup.vPlane(ipl);
#  ifdef EXTENSIVE_DICO
	if (ipl<setup.vIplLast()[1]) {
	  const TStation *&s = p.Station; if (s!=sPrv) {
	    sPrv = s; if (nProjs>=2) { nSpacePts++; nProjs=projs = 0; }
	  }
	  int jproj = 1<<p.IProj;
	  if ((jproj&projs)==0)    { projs |= jproj;    nProjs++; }
	}
	else if (firstZ2) {
	  firstZ2 = false; if (nProjs>=2) nSpacePts++;
	  if (nSpacePts<2) {  // Case not enough ``space points'': full TTrack,
	    // including 0x2 segment is to be fitted => allow some versatility
	    ft.par_pattern = 0x1f; nh0 = 6; nh1 = 5; // => relax momentum
	    // (N.B.: the case where there's not enough info in the bending
	    // dimensions for the free momentum fit to be precise (e.g. when
	    // having in view to extrapolate the track upstream of its reco
	    // zones) is not dealt with here and should be pinned down before
	    // entering "QNewtonFit" and taken care of elsewhere.)
	  }
	}
#  endif
	if (ipl<=setup.vIplLast()[1] || nSpacePts<2) {
	  int idx; if ((idx = Lat->dico_idx[ipl])>=0) {
	    hits->h[idx] = h.U; hits->hp[idx/32] |= 1<<(idx%32); NDics++;
	    resols[idx] = h.SigU;
	    if      (p.IFlag&0x10) {        // P-pixel planes of XY kind
	      // Y follows immediately X in TLattice's "dicoDets"
	      int jdx = idx+1;
	      hits->h[jdx] = h.V; hits->hp[jdx/32] |= 1<<(jdx%32); NDics++;
	      resols[jdx] = h.SigV;
	    }
	    else if (p.IFlag&0x20) {        // P-pixel planes of UV kind
	      int jdx = idx+1;
	      hits->h[jdx] = -h.V; hits->hp[jdx/32] |= 1<<(jdx%32); NDics++;
	      resols[jdx] = h.SigV;
	    }
	  }
	}
      }
      ih++;
    }
  }
  else if (mode==8 || mode==10) {
    //   ***** SPECIAL CASE W/ FIXED Y VERTEX|MOMENTUM, UPSTREAM of RICH *****
    // We want the fit not to be disturbed by muliscattering in RICH
    // => strip tail
    NDics = 0; while (ih!=lHitPat.end()) {
      if (*ih>=0) {
	const THit &h = vHit[*ih]; int ipl = h.IPlane; 
	if (ipl<setup.vIplLast()[1]) {
	  const TDetect &d = setup.iPlane2Detect(ipl);
	  if (d.X(0)<900) {  // Require upstream of RICH
	    int idx; if ((idx = Lat->dico_idx[ipl])>=0) {
	      hits->h[idx] = h.U; hits->hp[idx/32] |= 1<<(idx%32); NDics++;
	      resols[idx] = h.SigU;
	      if      (setup.vPlane(ipl).IFlag&0x10) {// P-pixel planes, XY kind
		int jdx = idx+1;
		hits->h[jdx] = h.V; hits->hp[jdx/32] |= 1<<(jdx%32); NDics++;
		resols[jdx] = h.SigV;
	      }
	      else if (setup.vPlane(ipl).IFlag&0x20) {// P-pixel planes, UV kind
		int jdx = idx+1;
		hits->h[jdx] = -h.V; hits->hp[jdx/32] |= 1<<(jdx%32); NDics++;
		resols[jdx] = h.SigV;	    
	      }
	    }
	  }
	}
      }
      ih++;
    }
  }
  else {                                            // ***** STANDARD CASE *****
    NDics = 0; while (ih!=lHitPat.end()) {        // Loop over the hits
      if (*ih>=0) {                                   // Hit on current plane
	const THit &h = vHit[*ih]; int ipl = h.IPlane, idx;
	if ((idx = Lat->dico_idx[ipl])>=0) {          // Dico is not universal
	  hits->h[idx] = h.U; hits->hp[idx/32] |= 1<<(idx%32); NDics++;
	  resols[idx] = h.SigU;
	  if      (setup.vPlane(ipl).IFlag&0x10) {    // P-pixel planes, XY kind
	    int jdx = idx+1;
	    hits->h[jdx] = h.V; hits->hp[jdx/32] |= 1<<(jdx%32); NDics++;
	    resols[jdx] = h.SigV;
	  }
	  else if (setup.vPlane(ipl).IFlag&0x20) {    // P-pixel planes, UV kind
	    int jdx = idx+1;
	    hits->h[jdx] = -h.V; hits->hp[jdx/32] |= 1<<(jdx%32); NDics++;
	    resols[jdx] = h.SigV;	    
	  }
	}
	//else  if (idx==-2) break;   // Why not? Have to check...
      }
      ih++;
    } // End of loop over planes and hits
  }
  if (NDics<=nh0) return false;

  //                      *************** FIT ***************

  //#define DEBUG_DICO
#  ifdef DEBUG_DICO 
  Track ot = ft;
#  endif     
  int status = Lat->Chi2Fit(&ft);
  if (!status) {
    if (Lat->FarOutDico(par) && // Far outside Dico on view of final pars..
	(TOpt::ReMode[14]&0x4)) // ...while using Dico...
      return false;
    if (ft.chi2/(NDics-nh1)<TOpt::dCut[9]) {

      // ********** CONVERT BACK Track-> Helix **********

      // HOUT(0) set @ Dico ref. plane =   target
      float dZStart = Lat->dZStart; //   + dZ over which a straight line extrapolation is enforced  
      HOUT(0) = setup.TargetCenter[0]+Lat->dZStart;
      // Extrapolate from Dico definition plane(= target) to ref. plane, hence
      // in straight track line
      HOUT(1) = par->x+dZStart*par->tx;
      HOUT(2) = par->y+dZStart*par->ty;
      HOUT(3) = par->tx; HOUT(4) = par->ty;
      HOUT(5) = par->pxzi/sqrt(1+par->ty*par->ty/(1+par->tx*par->tx));
      // Do not change Hlast, even its momentum component
      Chi2tot = ft.chi2;

#  ifdef EXTENSIVE_DICO
      if (nh1>=5) {      // ***** FULL FIT *****
	IFit = 0x8; // Set fit type bit pattern to == QN (it overwrites chi2tot)
	Chi2aux = Chi2tot;  // Saveguard Chi2tot for it might be overwritten
      }
#  endif

      //#define DEBUG_DICO
#  ifdef DEBUG_DICO 
      if (TOpt::Print[5]>0) {
	int idx;
	for (idx = 0; idx<NROWS; idx++) {
	  if (ot.hits.hp[idx/32]&(1<<idx%32))
	    printf("%8.2f ",ft.hits.h[idx]-ot.hits.h[idx]);
	  else
	    printf("         ");
	  if (idx%8==7) printf("\n");
	}
	printf("\n");
      }
#  endif     
#  ifdef EXTENSIVE_DICO
      // ********** GUESTIMATES **********

      int igroup, ipl0, guests_groups = GuestsGroups; GuestsGroups = 0;
      for (igroup= 0, ipl0 = -1; igroup<3; igroup++) {
	if (InGroup(igroup) ||    // Guests in the sole groups where hits
	    all_guests) {         // Guests everywhere, in view of extrapolating
	  if (ipl0<0) ipl0 = setup.vIplFirst()[igroup];
	  if ((1<<igroup)&guests_groups) {
	    // ***** UPDATE *****
	    for (int ipl = setup.vIplFirst()[igroup];
		 ipl<=setup.vIplLast()[igroup]; ipl++) {
	      int idx = Lat->dico_idx[ipl];
	      if      (idx>=0) vGuests[ipl-ipl0] = hits->h[idx];
	      else if (idx==-2) break;
	    }
	  }
	  else {
	    if (igroup==0 && guests_groups) {
	      vGuests.erase(vGuests.begin(),vGuests.end()); // ***** RESET *****
	      guests_groups = 0;
	    }
	    // ***** INITIALIZE *****
	    vGuests.reserve(TLattice::IplLast);
	    for (int ipl = setup.vIplFirst()[igroup];
		 ipl<=setup.vIplLast()[igroup]; ipl++) {
	      int idx = Lat->dico_idx[ipl];
	      if      (idx>=0) vGuests.push_back(hits->h[idx]);
	      else if (idx==-1) vGuests.push_back(0);
	      else if (idx==-2) break;
	    }
	    guests_groups |= 1<<igroup;
	  }
	}
      }
      GuestsGroups = guests_groups;
#  endif
    }
    else return false;
  }
  else
    return false;
#endif
  return true;
}
