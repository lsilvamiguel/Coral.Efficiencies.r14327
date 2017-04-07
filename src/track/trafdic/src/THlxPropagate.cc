// $Id: THlxPropagate.cc 13148 2011-12-28 16:55:25Z kbicker $

/*!
  Extrapolate helix w/ Multi-Scattering from material maps AND from detectors outside maps (Note that the latter is not taken into account in standard "THlx::Extrapolate" and may be signficant, e.g. in DRs equipped w/ lead).
*/

#include <float.h>
#include "CsErrLog.h"
#include "TOpt.h"
#include "TSetup.h"
#include "TDisplay.h"
#include "THlx.h"
#include "CsGeom.h"

using namespace std;

// The method extrapolates
//  - Through the full (up to targeted end point) thickness of the material maps encountered
//  - Outside these maps, from detector to detector, checking the track is traversing them,
//    I) either in their sole active area,
//   II) or anywhere w/in the extent of their outer frame, including the central dead zone.
//   This makes a difference for detector
//      - massive,
//      - w/ an empty dead zone,
// hence particularly so for the 2006/7 version of the DRs, w/ their thick layer of lead. Here we use version (I) in addition for ST04, MBs, MAs and Hodoscopes.
//   This->THlx is most probably the result of a fit that ends up at a particular detector. That detector (in practice, the one in the neighbourhood of the starting point) is skipped, on the ground that it has been taken into account in the fit. The targeted end point, if a detector, is, on the contrary, taken into account (so that we have a consistent scheme over successive executions of fit and propagation).

bool THlx::Propagate(THlx &Hout) const
{

  if (&Hout==this) {
    cout<<"THlx::Extrapolate ==> Output helix must differ from 'this' helix"<<endl;
    assert(false);
  }

  double x = Hpar[0], xStop = Hout.Hpar[0];

  if (fabs(xStop-x)<FLT_EPSILON) { // Trivial case: extrapolate to same X
    Hout = *this; Hout.path = fabs(xStop-x); Hout.radLenFr = 0;
    return true;
  }

  if (TOpt::ReMode[20]>0 && with_mom()) { // Material maps are ON, track w/ P
    CsMaterialMap *map = CsGeom::Instance()->getCsMaterialMap();
    if (!map) CsErrLog::mes(elFatal,
 "Extrapolation with material map requested (cf. \"TraF ReMode[20]\"), while no material map avalable!");

    if (map->usingROOTGeometry()) {
      const TSetup &setup = TSetup::Ref();
      THlx Hin;
      float  Len, RadLen, Step;
      bool forward = xStop>x;

      int iter = 0;
      double tot_len(0), tot_rad_len(0), tot_eloss(0); Hout = *this;
      do {     // ***** ...=> ITERATIVE EXTRAPOLATION THROUGH MAP
	// Get media properties and recommended step
	map->getRadLength(Hout,forward,RadLen,Step);
#define THlx_DEBUG_PROPAGATE
#ifdef THlx_DEBUG_PROPAGATE
	// Draw point (if requested) where material map was called
	if (TDisplay::Ptr()) { // TDisplay object exists
	  if (TDisplay::Ref().Mmaptraj()) {
	    double r[] = {Hout(0),Hout(1),Hout(2)};
	    TDisplay::Ref().point(r,110);
	  }
	}
#endif
	// Prepare next step
	  
	if (forward) x+Step<xStop ? x += Step : x = xStop;
	else         x+Step>xStop ? x += Step : x = xStop;

	Hin = Hout;
	if (!Hin.Extrap(x,Hout))     // Do extrapolation
	  return false;                          // Give up if in error

	Len = Hout.Path(); tot_len += Len; // Update trajectory path and...
	tot_rad_len += Len/RadLen;         // ...fraction of rad. length

	Hout.AddNoise(Len,RadLen);         // Add multiple scattering
	if (TOpt::ReMode[20]>1) {          // Add energy losses
	  float eloss = map->getdE(Hin,Len); tot_eloss += eloss;
	  if (TOpt::ELossStraggling)       // Propagate straggling in cov. matrix
	    Hout.Hcov[14] = Hout.Hpar[5]*Hout.Hpar[5]*Hout.Hpar[5]*Hout.Hpar[5]
	      * ( Hin.Hcov[14]/Hin.Hpar[5]/Hin.Hpar[5]/Hin.Hpar[5]/Hin.Hpar[5] +
		  map->getdEStraggling(Hin,Len) );
	  if (forward) Hout.Hpar[5] /= 1-eloss/Hout.Mom();
	  else         Hout.Hpar[5] /= 1+eloss/Hout.Mom();
	}
	x = Hout.Hpar[0];
	if (++iter>5000) { 
	  cout<<"THlx::Propagate ==> # of iteration exceeds 5000."<<endl;
	  cout<<"Last step was @ Z = "<<x<<"  Z final is "<<xStop<<endl;
	  if (!Extrap(x,Hout)) return false;
	  x = Hout.Hpar[0];
	  break;
	}
      } while ((forward && x < xStop) || (!forward && x > xStop)); // End of iteration loop
      
      // Save total length passed and rad. len. fraction 
      Hout.path =     tot_len;
      Hout.radLenFr = tot_rad_len;
      Hout.eLoss =    tot_eloss;

      return true;
    } else { // Material maps are ON and track w/ P, not using ROOT geometry

      CsMaterialMap *map = CsGeom::Instance()->getCsMaterialMap();
      if (!map) CsErrLog::mes(elFatal,
			      "Extrapolation with material map requested (cf. \"TraF ReMode[20]\"), while no material map avalable!");

      const TSetup &setup = TSetup::Ref();
      THlx Hin; float  Len, RadLen, Step;
      bool forward = xStop>x;

      double margin = // 1 micron of margin. Expected to be less than the distance between any 2 detectors and yet larger than the rounding error
	forward ? .0001 : -.0001;

      //   ***** INIT SEARCH for NEXT DET *****
      const vector<TDetect> &dets = setup.vDetect(); static const TDetect *det;
      int nDets = dets.size(), idet; double xDet;
      if (forward) {
	for (idet = 0,       xDet = xStop+2*margin; idet<nDets; idet++) {
	  det = &dets[idet]; xDet = det->X(0);
	  if (xDet-margin<x) continue; // => Skip detector at starting point
	  break;
	}
      }
      else {
	for (idet = nDets-1, xDet = xStop+2*margin; idet>=0;    idet--) {
	  det = &dets[idet]; xDet = det->X(0);
	  if (xDet-margin>x) continue; // => Skip detector at starting point
	  break;
	}
      }


      //   ***** INIT SEARCH for NEXT MAP *****
      const vector<float> &mapBorders = setup.vMaterialBorders();
      int nMaps = mapBorders.size(), imap;
      double xMap0, xMap1; if (forward) {
	imap = 0;       xMap0 = mapBorders[imap++]; xMap1 = mapBorders[imap++];
      }
      else {
	imap = nMaps-1; xMap0 = mapBorders[imap--]; xMap1 = mapBorders[imap--];
      }

      double tot_len(0), tot_rad_len(0); Hout = *this;
      do {
	bool winMap; if (forward) {                      // ***** UPDATE NEXT MAP?
	  while (xMap1<=x && imap<nMaps-1) {
	    xMap0 = mapBorders[imap++]; xMap1 = mapBorders[imap++];
	  }
	  winMap = xMap0<x && x<xMap1;
	}
	else {
	  while (xMap1>=x && imap>0) {
	    xMap0 = mapBorders[imap--]; xMap1 = mapBorders[imap--];
	  }
	  winMap = xMap1<x && x<xMap0;
	}
	if (winMap) {                                        // ***** W/IN MAP?...
	  int iter(0); do {     // ***** ...=> ITERATIVE EXTRAPOLATION THROUGH MAP
	    // Get media properties and recommended step
	    map->getRadLength(Hout,forward,RadLen,Step);
#ifdef THlx_DEBUG_PROPAGATE
	    // Draw point (if requested) where material map was called
	    if (TDisplay::Ptr()) { // TDisplay object exists
	      if (TDisplay::Ref().Mmaptraj()) {
		double r[] = {Hout(0),Hout(1),Hout(2)};
		TDisplay::Ref().point(r,110);
	      }
	    }
#endif
	    // Prepare next step
	  
	    if (forward) x+Step<xMap1 ? x += Step : x = xMap1;
	    else         x+Step>xMap1 ? x += Step : x = xMap1;

	    Hin = Hout;
	    if (!Hin.Extrap(x,Hout))     // Do extrapolation
	      return false;                          // Give up if in error

	    Len = Hout.Path(); tot_len += Len; // Update trajectory path and...
	    tot_rad_len += Len/RadLen;         // ...fraction of rad. length

	    Hout.AddNoise(Len,RadLen);         // Add multiple scattering
	    if (TOpt::ReMode[20]>1) {          // Add energy losses 
	      if ( TOpt::ELossStraggling )     // Propagate sigma from energy loss straggling to covariance matrix
		Hout.Hcov[14] = Hout.Hpar[5]*Hout.Hpar[5]*Hout.Hpar[5]*Hout.Hpar[5]
		  * ( Hin.Hcov[14]/Hin.Hpar[5]/Hin.Hpar[5]/Hin.Hpar[5]/Hin.Hpar[5] + map->getdEStraggling(Hin,Len) );
	      if ( forward )
		Hout.Hpar[5] /= 1-map->getdE(Hin, Len)/Hout.Mom();
	      else
		Hout.Hpar[5] /= 1+map->getdE(Hin, Len)/Hout.Mom();
	    }
	    x = Hout.Hpar[0];
	    if (++iter>5000) { 
	      cout<<"THlx::Propagate ==> # of iteration exceeds 5000."<<endl;
	      cout<<"Last step was @ Z = "<<x<<"  Z final is "<<xMap1<<endl;
	      if (!Extrap(xMap1,Hout)) return false;
	      x = Hout.Hpar[0];
	      break;
	    }
	  } while (fabs(x-xMap1)>FLT_EPSILON); // End of iteration loop
	  x = xMap1+margin; // Beyond map's end: next loop won't re-enter it
	}
	else {
	  //  ***** SEARCH NEXT GOAL: DET, MAP or END POINT ***** 
	  double xNext; bool isDet; if (forward) {
	    while (xDet<x && idet<nDets) {          // ***** UPDATE NEXT DETECTOR?
	      det = &dets[idet++]; xDet = det->X(0);
	    }
	    if (idet==nDets)          // End of det list: next must be "xStop" ...
	      xDet = xStop+2*margin;  // ...make so that condition infra fulfilled
	    xNext = xDet; isDet = true;
	    if (xNext>xMap0)        { xNext = xMap0; isDet = false; }
	    if (xNext>xStop+margin) { xNext = xStop; isDet = false; }
	  }
	  else {
	    while (xDet>x && idet>=0) {
	      det = &dets[idet--]; xDet = det->X(0);
	    }
	    if (idet<0)               // End of det list: next must be "xStop" ...
	      xDet = xStop+2*margin;  // ...make so that condition infra fulfilled
	    xNext = xDet; isDet = true;
	    if (xNext<xMap0)        { xNext = xMap0; isDet = false; }
	    if (xNext<xStop+margin) { xNext = xStop; isDet = false; }
	  }

	  // ***** EXTRAPOLATE to EITHER NEXT DET or NEXT MAP or END POINT *****
	  Hin = Hout;
	  if (!Hin.Extrap(xNext,Hout))   // Do extrapolation
	    return false;                            // Give up if in error

	  tot_len += Hout.Path();              // Update trajectory path

	  if (isDet) {  // ***** IF IT'S NEXT DET: ADD MULTI-SCATTERING *****
	    int type = det->IType;
	    if (type!=14 /* DR or ST04 */ && type!=18 /* MB */ &&
		type<39  /* MA and H */ &&  // If detector w/ massive dead zone...
		det->InFrame (Hout(1),Hout(2)) ||  // ...and w/in frame
		det->InActive(Hout(1),Hout(2))) {  // or w/in active area
	      double Len    = det->Siz(0)/Hout.DirCos(1); // Det. thickness
	      double RadLen = det->RadLen;                // Det. rad. len.
	      tot_rad_len += Len/RadLen;       // Update fraction of rad. length
	      Hout.AddNoise(Len, RadLen);      // Add multiple scattering
	      // No energy loss in detectors: no info about their stopping power.
	    }
	  }

	  x = xNext+margin;
	
	}
      } while ((forward&&x<xStop) || ((!forward)&&x>xStop));

      // Save total length passed and rad. len. fraction 
      Hout.path     = tot_len;
      Hout.radLenFr = tot_rad_len;

      return true;

    }
  }
  else { // Use of material map is OFF or helix without momentum.

    return(this->Extrap(xStop, Hout));
    
  }
  return(false);
}








