// $Id: TEvTracksFit.cc 13308 2012-03-21 19:53:38Z ybedfer $

/*! 
  \brief Basic method for fitting the list of tracks

  - The method is brought into play by option "TraF ReMode[13] 0", cf.
   "CsTrafficFitting".
  - The very most basic working is obtain when "TraF ReMode[11|12] 0".
  - Otherwise some features of TEv::TracksFit2" are included, viz.:
    - Beam momentum assigned a fix value, whether it's final (case of hadron
     beam, w/o BMS) or temporary (later to be updated by BMS). This is based on
     options dCut[4,5], iCut[15], conditioned by the later being !=0:
        P=i[15]/d[4], d(1/P) = sqrt(d[5]). 
    - In MC case, provided "ReMode[25]&0x1": a posteriori BMS simulation.
    - Fringe-field fit based on P inherited from "BridgeSegments2".
    - Erasing fragments left over by "BridgeSegments2".
    - etc...
*/

#include <iostream>
#include "TEv.h"
#include "TOpt.h"
#include "TDisplay.h"

using namespace std;

void TEv::TracksFit()
{

  bool reMode2 = // Tracking (PR, Bridging) track alla TrafDic (ReMode==2)
    TOpt::ReMode[11]==2 && TOpt::ReMode[12]==2;

  //        ******************** BEAM TRACKS ********************
  // I.e., in fact, beam telescope
  //  - Either a standard beam track, i.e. restricted to the 0x10 zone of the
  //   beam telescope alone,
  //  - Or a beam track bridged over the target (i.e. 0x11 type). Whether it's
  //   over a target dipole ("MagFlag1[0]>2") or not does not make much
  //   difference: the dipole strength is too low for a precise determination
  //   of a beam particle's momentum (expected to be of the order of 100 GeV).
  if (reMode2) BeamsFit2(); // Special handling of beams alla TraFDic

  bool ret;
  list<TTrack>::iterator it = listTrack.begin(); 
  while (it != listTrack.end()) {
    //  ************************* LOOP OVER TRACKS *************************

    if ((*it).NGroups()==0) {
      cout<<"TEv::TracksFit ==> track with NGroups == 0"<<endl;
      assert(false);
    }

    int dofit = 1; if (reMode2) { // Special processings alla TraFDic

      TTrack &t = *it;
      if (t.Type==0x10 || t.Type==0x11) {  // Beams: already processed supra.
	it++; continue;
      }

      if (t.IMark>0 && !(t.Type&0x2)) {
	//                 ***** DOWMSTREAM SEGMENT OF YOKE TRACK *****
	t.Hfirst(5,5) = 1.E-20;
      }
      else if (t.IMark<0) //   ***** TRACK PIECE LEFT OVER FROM BRIDGING *****
	// - One expects to catch here the single-zone segments already merged
	//  w/ (appended to) another track and awaiting the validation of the
	//  merging for their fate to be decided.
	// - "BridgeSegments2" leaves over such track segments in the case of
	//  the bridging over the (downstream) muon absorber (and only in that
	//  case). The reason being that there, the simplistic "TTrack::QuickKF"
	//  method is used which cannot do a decent job in solving the difficult
	//  task of bridging over the thickness of the absorber.
	dofit = 0;
      else if (t.NGroups()==1) {
	// ********** SINGLE SEGMENT, NON-BEAM, NON SM2-YOKE, TRACK **********
	// - Either magnet-OFF data (I), fringe-field track (II) or else
	//  fragments of reconstruction (III).
	// - In all cases, assigning momentum allows to account for multi-
	//  scattering in the fit, and, technically, enables smoothing points,
	//  i.e. best track estimators, that are only computed, at any given
	//  point, by the FullKF.
	//  I) FullKF and later on smoothing can be obtained via option
	//   "ReMode[24]". The momentum can only be the beam line tune, which is
	//   expected to be specified by "iCut[15], dCut[4] and [5]".
	// II) The bending through fringe-field allows to determine momentum,
	//   albeit w/ limited precision. The starting value for the momentum is
	//   retrieved, if availbale, from the dico fit. The whole thing is
	//   assumed to be of little interest in case (I), and therefore ignored
	//   when "ReMode[24]".
	// III) No action explicitly intended. When "ReMode[24]", the fragments
	//   get assigned the momentum of the beam line tune.
	if (TOpt::ReMode[24]==0 ||// But do fit still, if smoothing required or...
	    TOpt::ReMode[23]) {   // ...if explictly requested by option...
	  t.Hfirst(5) = TOpt::iCut[15]/TOpt::dCut[4];
	  t.Hfirst(5,5) = 1.E-20; // ...variance ~= 0, to fix momentum
	}
#ifdef FRINGE_FIELD_KF
	// - Require a finite starting value that cannot but have been assigned
	//  by "TTrack::IFit&0x8", cf. "BridgeSegments2".
	else if ((t.IFit&0x8) &&                       // If QN fit and...
		 (TOpt::dCut[66]<0 ||
		  fabs(1/t.Haux(5))<TOpt::dCut[66]) && // ...low enough momentum
		 t.InGroup(0)) {                       // ...and in zone 0x1
	  t.Hfirst(5) = t.Haux(5);        // => KF full fit
	  t.Hfirst(5,5) = 1.E-4;          // ...w/ large initial momemtum error
	}
#endif
	else {
	  t.Hfirst(5) = 0.; dofit = 0;   // ...Else nothing...and reset P (to be safe)
	}
      }
      else {// ********** SOME SPECIAL CASES of MULTI SEGMENT TRACKS **********

	if (t.Type==0xc) {                  // Mu-wall momentum-less
	  t.Hfirst(5) = TOpt::iCut[15]/TOpt::dCut[4];
	  t.Hfirst(5,5) = 1.E-20;             // ...variance ~=0 to fix momentum
	}
	else if (!t.Hfirst.with_mom())
	  t.Hfirst = t.Haux;                // QN only so far
      }
    } // End special processing alla TraFDic

    if (!dofit) { (*it).UseHitTime(); it++; continue; }

    // Insert misses planes
    //
    // It is needed as missed planes (it's where hit is not found) may add
    // mult. scattering.
    // To be done before final fit.
    //
    if(TOpt::ReMode[21] == 0) 
      (*it).InsertMissedPlanes();

    if ((*it).Hfirst.empty() ||      // If track has not yet been fitted
	!(*it).Hfirst.with_mom()) {  // ...or is momentum-less
      it++; continue;                // => Skip
    }

    //
    // Full (with calculation of smoothed helices) Kalman fit forward and backward
    //
    if (!(*it).Hfirst(5,5)) {   // q/P^2 component of the cov. matrix =0?...
      (*it).Hfirst(5,5) = 1.E-4;// ...initialize it.
      CsErrLog::msg(elError,__FILE__,__LINE__,// Yet, shouldn't happen=> warning
		    "Evt #%d Trk #%d(0x%x) has null d2(q/P)=0 while q/P finite",
		    event,(*it).Id,(*it).Type);
    }
    ret = (*it).FullKF(1); // -----> forward
    if(!ret) { // forward  fit has failed. Erase track. Next one.
      cout<<"TEv::TracksFit() ==> Forward Kalman fit 1 failed.  Track ID "
	  <<(*it).Id<<" removed"<<endl;
      listTrack.erase(it++); 
      goto next_track;
    }

    ret = (*it).FullKF(-1); // <------ backward
    if(!ret) { // backward  fit has failed. Erase track. Next one
      cout<<"TEv::TracksFit() ==> Backward Kalman fit 1 failed. Track ID "
	  <<(*it).Id<<" removed"<<endl;
      listTrack.erase(it++);
      goto next_track;
    }

    //
    // Track "refiner"
    //
    if(TOpt::ReMode[22] == 0) goto next_track; // skip all track improvement procedures
    ret=(*it).Refine(1);
    if(!ret) {
      cout<<"TEv::TracksFit() ==> Refining  failed. Track ID "
	  <<(*it).Id<<" removed"<<endl;
      listTrack.erase(it++); 
      goto next_track;
    }
    
    //
    // Fit backward/upward again
    //
    ret = (*it).FullKF(1); // -----> forward
    if(!ret) { // forward  fit has failed. Erase track. Next one.
      cout<<"TEv::TracksFit() ==> Forward Kalman fit 2 failed.  Track ID "
	  <<(*it).Id<<" removed"<<endl;
      listTrack.erase(it++); 
      goto next_track;
    }

    if(TOpt::Graph[6] > 0){ // switch ON debug drawing of calls to field and maperial maps
      TDisplay::Ref().SetMmaptraj(true); TDisplay::Ref().SetRkutraj(true);
    }

    ret = (*it).FullKF(-1); // <------ backward
    if(!ret) { // backward  fit has failed. Erase track. Next one
      cout<<"TEv::TracksFit() ==> Backward Kalman fit 2 failed. Track ID "
	  <<(*it).Id<<" removed"<<endl;
      listTrack.erase(it++);
      goto next_track;
    }
    TDisplay::Ref().SetMmaptraj(false); TDisplay::Ref().SetRkutraj(false); // switch OFF debug drawing


    // Cut
    if(TOpt::dCut[17] != 0. && ((*it).Chi2tot/(*it).NHits) > TOpt::dCut[17]){ // fit OK. Check Chi2       
      listTrack.erase(it++);     // Bad Chi2. Erase.
      goto next_track;
    }
    
    (*it).UseHitTime();     // Time measurement

    //Calculate "smoothed" helices (if requested) 

    for(int i = 0; i < int(sizeof(TOpt::SmoothPos)/sizeof(double)); i++){ // loop over "smoothing" points
      double x = TOpt::SmoothPos[i];
      if(x == 0.) break;
      THlx Hsm;
      if((*it).GetSmoothed(Hsm,x) < 0) continue;
      (*it).lHsm.push_back(Hsm); // store smoothed
    }// end of loop over "smoothing" points

    it++; // next track
  next_track:;
  } // end of loop over tracks



  // In case of MC event:
  // - assosiate found track with MC track.
  // - extrap. tracks upwards to prim. vertex if requested,
  //   (for tracks starting befor SM1)
  if(! this->isMC ) return;

  if(TOpt::Graph > 0){ // switch on debug drawing of calls to field and maperial maps
    TDisplay::Ref().SetMmaptraj(true); TDisplay::Ref().SetRkutraj(false);
  }

  for(it = listTrack.begin(); it != listTrack.end(); it++) {  // loop over tracks
    (*it).FindKine();
    if(TOpt::ReMode[9] == 1 && (*it).InGroup(0) == 1){
      if( (*it).Hfirst.with_mom() ) { // track with momentum
	THlx Hextr;
	Hextr(0)=vVtxMC(0).V(0);   // coord of prim. vertx. (assumed to be the first one)
	if((*it).Hfirst.Extrapolate(Hextr)) (*it).Hfirst = Hextr;
      }
    }
  }
  
  TDisplay::Ref().SetMmaptraj(false); TDisplay::Ref().SetRkutraj(false); // switch off debug drawing


}
