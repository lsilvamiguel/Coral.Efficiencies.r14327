// $Id: TTrackFullKF.cc 14069 2015-09-17 20:44:46Z lsilva $

#include <iostream>
#include <iomanip>
#include "TEv.h"
#include "TSetup.h"
#include "TTrack.h"
#include "THit.h"

using namespace std;

/*!
 Kalman fit through already found hits (no hits re-assignment).
 - Material traversed is taken into account, via ROOTGeometry or a combination of material maps and TDetect::RadLen, as specified by options.
 - Accumulated fraction of rad. length and energy loss are stored into "radLenFr" and "eLoss" data member. Nota bene that e-loss is only computed while ROOTGeometry is enabled or in special mat.maps, and never while traversing discrete detectors.
 - Upon option, refitting while recycling the initial part of the previous fit.
  NOTA BENE: The chi2 corresponding to the initial part of the refit (that which is a mere recycling) is taken from "TTrack::mChi2" which holds the results of the lastly performed fit.
   => Refit has to be performed in the same direction as the latter.
   => This should be changed in the future, by storing 2 chi2 maps, one for each direction.
 - The cov. matrix is scaled up prior to the fitting computation proper, so as to forget about the a priori setting of the input parameter vector. One can work around the covariance upscaling feature in order to fix a particular parameter, e.g. q/P, by assigning to the corresponding diagonal term of the matrix a very small (yet finite) value, before calling "FullKF".
 - One can also bypass the upscaling, when fitting the extension of a track previously fitted: in that case, old information is not reprocessed and one des not have to forget about the previous fit to it.
 - In some cases, the KF can be sent far away from the true trajectory, particularly when the track is refitted backward starting from the MA/HG02 detectors. As explained in "wwwcompass.cern.ch/twiki/bin/view/DataReconstruction/MuWallAReconstruction". It then takes into account irrelevant MultiScattering and Energy-Loss. A solution to this problem is to derive MS and EL from the trajectory found in the forward fit and stored in the map of mHef helices.

 \param dir  fit direction

 \verbatim
 Input:

   dir =  1 - fit in  forward direction (resulting helix is Hlast )
   dir = -1 - fit in backward direction (resulting helix is Hfirst)

 \endverbatim
 \param ipl  If >=0, only a refit is performed, starting from plane <ipl>. If <ipl> turns out to lies outside track's range, the first (or last) hit w/in range is chosen instead.
*/

/*
  Version specific to "lattice" alternative.
  Changes w/ respect to "traffic/sources/TTrackFullKF.cc".
    i) Noise is not added if track passes through an empty dead zone, as is the case for, e.g., DRs and their thick layer of Pb in 2006/7.
*/

bool TTrack::FullKF(int dir, int ipl, bool scaleUp)
{
  if (TOpt::ReMode[24] &&  // ***** UPON OPTION...
      (dir>0 && !Hfirst.with_mom() || dir<0 && !Hlast.with_mom()))
      //                      ***** MOMENTUM-LESS TRACKS => SIMPLE FIT
    return this->QuickKF(dir,0);

  //         ********** INITIALIZATIONS, CHECKS **********
  const TSetup &setup = TSetup::Ref(); TEv &ev = TEv::Ref();
  const vector<THit>& vHit = ev.vHit(); bool ok = true, print = false;
  if (dir!=1 && dir!=-1) {
    cout<<"TTrack::FullKF ==> wrong direction argument : "<<dir<<endl;
    assert(false);
  }
  if (lPlnRef.size()!=lHitPat.size()) {
    cout<<"TTrack::FullKF ==> lPlnRef.size() != lHitPat.size() "<<endl;
    assert(false);
  }
  //#define KF_DEBUG_KF_PATH
#ifdef KF_DEBUG_KF_PATH
  // The KF can be sent off trajectory in its earlier steps (when the covariance
  // is still unconstrained): this patch allows to follow its path graphically.
  if (TOpt::Graph[6]&0x8) {
    // Switch ON debug drawing of calls to field and material.
    TDisplay::Ref().SetMmaptraj(true); TDisplay::Ref().SetRkutraj(true);
  }
#endif
  if (dir<0) {              // ***** BACKWARD: REVERSE HIT LIST, SWAP END POINTS
    lPlnRef.reverse(); lHitPat.reverse(); swap(Hfirst,Hlast);
  }
  list<int>::iterator ip = lPlnRef.begin(), ih = lHitPat.begin(); THlx H;
  Chi2tot = 0;
  bool all = ipl<0 ||
		 // Can be that partial fit be requested at first(last) plane 
		 dir>0 && ipl<=lPlnRef.front() || dir<0 && ipl>=lPlnRef.front();
  if (all) {                                 // ***** ALL-OUT FIT: CLEAR MAPS...
    mChi2.clear();
    if (dir>0) { mHef.clear(); mHuf.clear(); }
    else {       mHeb.clear(); mHub.clear(); }
    H = Hfirst;
    if (scaleUp) H *= 100;           // ...DIAGONALIZE AND SCALE UP COV. MATRIX
  }
  else {                         // ***** PARTIAL REFIT: RE-INITIALIZE ITERATORS
    const map<int,THlx> *mH;
    if (dir>0) {
      mH = &mHuf; map<int,double>::iterator idx = mChi2.begin();
      while (ip!=lPlnRef.end() && *ip!=ipl) {
	if (*ih>=0) { Chi2tot += (*idx).second; idx++; }
	ip++; ih++;
      }
    }
    else {
      mH = &mHub; map<int,double>::reverse_iterator idx = mChi2.rbegin();
      while (ip!=lPlnRef.end() && *ip!=ipl) {
	if (*ih>=0) { Chi2tot += (*idx).second; idx++; }
	ip++; ih++;
      }
    }
    list<int>::iterator jp = ip; jp--; int jpl = *jp;
    if (ip==lPlnRef.end() && // Case: starting point is outside track
      // (This happens when FullKF is called for refitting after the hit
      // orginally last in the "dir" direction has been cancelled and the list
      // of plane references has been consequently cancelled (note that the case
      // of an originally first hit cancelled is dealt with supra).)
	(dir>0 && ipl<jpl || dir<0 && ipl>jpl))
      CsErrLog::msg(elFatal,__FILE__,__LINE__,"Track #%d: refitting from "
		    "plane #%d, which is not inserted in track",Id,ipl);
    map<int,THlx>::const_iterator im = mH->find(jpl);
    if (im==mH->end())
      CsErrLog::msg(elFatal,__FILE__,__LINE__,"Track #%d: refitting from "
		    "plane #%d, which is not in track helix map",Id,jpl);
    H = (*im).second;
  }
  unsigned int muWABackFit = 0;// Backward refit starting from muWallA sub-zone?
  // => Backward MWA fit: recylcing radlength/eloss info derived in forward pass
  if (all && dir<0 && TOpt::ReMode[20] &&
      // For radlength and eloss, we need rely on the info in the mHef helices
      mHef.size()!=0) {
    if (setup.InMuWallA(*ip)) {
      list<int>::iterator jp = ip, jh = ih; int ok;
      for (jp++, jh++; jp!=lPlnRef.end(); jp++, jh++) {
	if (*jh>=0) {
	  if (setup.InMuWallA(*jp)) muWABackFit = 0x3;
	  break;
	}
      }
    }
  }

  if (print) cout<<endl<<endl;
  THlx Hext, Hupd;
  int nstep(0); double tot_len(0), tot_rad_len(0), tot_eloss(0);
  while (ip!=lPlnRef.end()) {
    //    ********** LOOP OVER PLANES and HITS (or HIT PLACEHOLDERS) **********
    nstep++;
    if (print)
      cout<<"-----> pl "<<(*ip)<<"  proj = "<<setup.vPlane(*ip).IProj<<endl;
    
    const TDetect &det = setup.iPlane2Detect(*ip);

    //                                                         ***** EXTRAPOLATE
    Hext(0) = det.X(0);
    if (!muWABackFit)
      ok = H.Extrapolate(Hext,true);  // Extrapolate (with use of material map)
    else {
      ok = H.Extrapolate(Hext,false); // Special Backward MWA fit
      if (ip!=lPlnRef.begin()) {
	list<int>::iterator jp = ip;
	map<int, THlx>::const_iterator im = mHef.find(*jp);
	if (im!=mHef.end()) Hext.SetqP((*im).second(5));
	if (jp!=lPlnRef.end()) {
	  jp--; im = mHef.find(*jp); if (im!=mHef.end()) {
	    const THlx &hef = (*im).second; Hext.SetPath(hef.Path());
	    double RadLen = hef.RadLenFr();
	    Hext.SetRadLenFr(RadLen); Hext.SetELoss(hef.ELoss());
	    if (RadLen) // If "RadLen" finte!
	      Hext.AddNoise(Hext.Path(),RadLen);
	  }
	}
      }
      if (*ih>=0 && !setup.InMuWallA(*ip)) {// Exiting special Backward MWA fit?
	float sa = setup.vDetect(*ip).Sa;
	if (sa<.05) muWABackFit &= 0x2;
	else        muWABackFit &= 0x1;
      }
    }
    if (!ok) {                                            // ***** EXIT IF FAILS
      cout<<"TTrack::FullKF() ==> Extrapolate error for track ID = "
	  <<Id<<"  dir = "<<dir<<" at step # "<<nstep<<endl;
      if (dir<0) { // Restore before exiting...
	// ...The idea being that, even if the present fit fails, this->TTrack
	// can still be usefull, and then better if the calling program can
	// rely on the TTrack's attributes to reconstitute a valid track out of
	// it and fit it.
	lPlnRef.reverse(); lHitPat.reverse(); swap(Hfirst,Hlast);
      }
      return false;
    }
    if (print) Hext.Print("Hext");

    // Cumulate pass, fraction of rad. length and ELoss
    tot_len     += Hext.Path();
    tot_rad_len += Hext.RadLenFr();
    tot_eloss   += Hext.ELoss();

    // Add scattering on detector if it's not yet done inside extrapolation
    // function, i.e. if map is OFF or map is ON but we are out of material map.
    if (TOpt::ReMode[20]<=0 ||
	TOpt::ReMode[20]>0 && !setup.InMaterialMap(Hext(0))) {
      if (det.InMassive(Hext(1),Hext(2))) {
	// If inside detector massive area, i.e. active zone + possibly central
	// dead zone if massive (case of e.g. MM, as opposed to ST, DR, MB).
	double Len    = det.Siz(0)/Hext.DirCos(1); // det. thickness
	double RadLen = det.RadLen;                // det. rad. len.
	tot_len     += Len;
	tot_rad_len += Len/RadLen;
	Hext.AddNoise(Len, RadLen);
	Hext.SetRadLenFr(Hext.RadLenFr()+Len/RadLen);
	if (print) Hext.Print("Hext + noise");
      }
      else if (print) cout<<"   track is out of plane. No noise added\n";
    }

    if (Hext.RadLenFr()<0) {
      cout<<"TTrack::FullKF() ==> Extrapolating from "<<H(0)<<" to "<<Hext(0)
	  <<" gives negative X/X0 = "<<Hext.RadLenFr()<<" (tot_rad_len = "<< tot_rad_len <<")"<<endl;
      assert(false);
    }

    //   *************** STORE EXTRAPOLATED HELIX ***************
    if (dir==1) mHef[*ip] = Hext;
    else        mHeb[*ip] = Hext;

    //           *************** UPDATE ***************
    if (*ih<0)                      // No associated hit
      Hupd = Hext;
    else {
      const THit &h = vHit[*ih];    // Associated hit

      double chi2 = Hext.HitChi2(h);
      mChi2[*ip] = chi2; // Store chi2 increment
      Chi2tot += chi2;
      if (print) {
	cout<<"     hit # "<<*ih
	    <<" det = "<<sqrt(1./(h.G0*h.G2 - h.G1*h.G1))
	    <<" chi2 "<<chi2<<endl;
      }
      Hext.Update(h,Hupd);        // <--------- Update
    }
    if (print) Hupd.Print("Hupdated");
    if (print) cout<<endl;

    //      *************** STORE UPDATED HELIX ***************
    if (dir==1) mHuf[*ip] = Hupd;
    else        mHub[*ip] = Hupd;

    H = Hupd; 
    ip++; ih++;  // Next plane, next hit

  } // End of loop over planes and hits

  Hlast = H;                       // ********** UPDATE
  if (ipl<0) {
    this->radLenFr = tot_rad_len;
    this->eLoss = tot_eloss;
    // Set fit flags
    if (dir==1) fFitDone = true;
    else        bFitDone = true;
  }
  if (dir<0) {// ***** BACKWARD: RESTORE HIT LIST ORDERING, SWAP BACK END POINTS
    lPlnRef.reverse(); lHitPat.reverse(); swap(Hfirst,Hlast);
  }

#ifdef KF_DEBUG_KF_PATH
  if (TOpt::Graph[6]&0x8) {
    TDisplay::Ref().SetMmaptraj(false); TDisplay::Ref().SetRkutraj(false); // Switch OFF debug drawing
  }
#endif
  return true;
}
