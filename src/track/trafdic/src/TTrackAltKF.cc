// $Id: TTrackAltKF.cc 13148 2011-12-28 16:55:25Z kbicker $

#include <iostream>
#include <iomanip>
#include "CsErrLog.h"
#include "TEv.h"
#include "TSetup.h"
#include "TTrack.h"
#include "THit.h"

using namespace std;

/*!
  Alternative full Kalman fits (cf. "TTrack::FullKF"). So far, only one alternative: "TTRack::KFwExtraMS", w/ extra MultiScattering.

 \param extraRadLenFr Extra MS added at the hinging point of the muFilter

 \verbatim
 - Backward direction => Resulting helix is Hfirst.
 - Various simplifying assumptions.
   - Argument track has momentum.
   - Use of material map has been requested.
 - Bypass all consistency checks (everthing's been already checked by the time
  the alternative method is called upon)
 - No use of drift info
 \endverbatim
*/

bool TTrack::KFwExtraMS(double extraRadLenFr)
{
  if (TOpt::ReMode[20]<=0) { // I.e. material maps not used...
    // ...the whole method makes no sense => Inconsistency => Exit
    printf("TTrack::KFwExtraMS: Inconsistency: Material maps not used, cf. \"TraF ReMode[20]\"\n"); exit(1);
  }

  const TSetup &setup = TSetup::Ref();
  TEv &ev = TEv::Ref(); const vector<THit>& vHit = ev.vHit();

  //      ***** INITIALIZATION of HELICES *****
  mHeb.clear(); mHub.clear(); // Backward direction => backward helices
  THlx Hext, Hupd, H(Hlast);
  H *= 100;                   // Diagonalize and scale up cov. matrix

  double tot_rad_len(0); Chi2tot = 0; bool extraAdded = false;
  list<int>::reverse_iterator rp = lPlnRef.rbegin() ;
  list<int>::reverse_iterator rh = lHitPat.rbegin() ;
  while (rp!=lPlnRef.rend()) { // ***** LOOP OVER PLANES, HITS *****
    const TDetect &d = setup.iPlane2Detect(*rp);

    //                ***** EXTRAPOLATE *****
    bool ok = true; if (!extraAdded && d.X(0)<setup.MuFilterHinges[1]) {
      Hext(0) = setup.MuFilterHinges[1]; if ((ok = H.Extrapolate(Hext,true))) {
	Hext.AddNoise(1,1/extraRadLenFr); H = Hext;
	tot_rad_len += Hext.RadLenFr() + extraRadLenFr;
	extraAdded = true;
      }
    }
    if (ok) {
      Hext(0) = d.X(0); ok = H.Extrapolate(Hext,true);
    }
    if (!ok) {
      CsErrLog::msg(elError,__FILE__,__LINE__,
 "Error extrapolating track #%d backward to \"%s\"",Id,d.Name.c_str());
      return false;
    }
    tot_rad_len += Hext.RadLenFr();       // Update radiation length

    // ***** ADD SCATTERING ON DETECTOR if not done in extrapolation... *****
    // ...i.e. if map is OFF or we are outside any material map
    if (!setup.InMaterialMap(Hext(0))) {
      if (d.InMassive(Hext(1),Hext(2))) {
	// If inside detector massive area, i.e. active zone + possibly central
	// dead zone if massive (case of e.g. MM, as opposed to ST, DR, MB).
	double Len    = d.Siz(0)/Hext.DirCos(1); // Thickness along path
	double RadLen = d.RadLen;                // Det. rad. len.
	tot_rad_len += Len/RadLen;
	Hext.AddNoise(Len,RadLen);
      }
    }  

    mHeb[*rp] = Hext;     // Store extrapolated helix

    //            ***** UPDATE *****
    if (*rh<0) Hupd = Hext;  // No hit associated w/ this track on current plane
    else {
      const THit &h = vHit[*rh]; THit hit = h;
      double chi2 = Hext.HitChi2(hit);
      mChi2[*rp] = chi2;     // Store Chi2 increment into the map
      Chi2tot += chi2;
      Hext.Update(hit,Hupd); // <--------- Update
    }

    mHub[*rp] = Hupd;   // Store updated helix

    H = Hupd; 
    rp++; rh++;
  } // End of loop over planes and hits

  Hfirst = H; radLenFr = tot_rad_len;
  bFitDone = true; // Set fit flag
  return true;
}









