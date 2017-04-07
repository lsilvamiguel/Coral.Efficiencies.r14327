// $Id: TTrackPrint.cc 14094 2015-11-06 15:28:48Z lsilva $

#include "TEv.h"
#include "TSetup.h"
#include "TTrack.h"

using namespace std;

// Track information printouts
// level 1 - brief 

/*
  Version specific to "lattice" alternative.
  Changes w/ respect to "traffic/TrackPrint":
   i) Handle, lattice-specific, THit's w/ 2 TTrack references (IKine+IKin2)
  ii) Don't display DC flags (they are not used in TraFDic)
 iii) Coalesced mirrors, cf. lattice's "ImportClusters", are tagged as such.
  iv) Include new method "UpdateProjs", which re-builds "sProj" from hits list (
     "sProj" is not updated by TraFDic throughout the construction of the list.)
   v) Print which-side info for split channels of straws (except ST04), MAs and
     stripped pixelGEM/MMs.
  vi) Print v-coord of pixels.
*/

void TTrack::Print(int level, string comment) const
{

  //       ***** SOME X-CHECKS *****
  if (lPlnRef.size()-lHitPat.size()!=0)
    cout<<"TTrack::Print ==> inconsistency in list sizes : "
	<<lPlnRef.size()<<" ; "<<lHitPat.size()<<endl;

  const TEv &ev = TEv::Ref();
  const TSetup &setup = TSetup::Ref();

  if (level==0) {   // ***** SHORT PRINT MODE *****
    printf("\nTrk ID: %d",Id);
    if (IMark>0)      printf(", [%d pass]",IMark);
    if (Associate>=0) printf(", Associate = %d",Associate);
    cout<<endl;
  }

  if (level==1) {   // ***** PRINT "HIT PATTERN" *****
    if(insPlaneDone && fFitDone && bFitDone) { // i.e. if missed planes are inserted and smoother helices are there
      THlx H; double chi2;
      list<int>::const_iterator ipl = lPlnRef.begin();
      list<int>::const_iterator ihp = lHitPat.begin(); 
      while(ipl != lPlnRef.end()){ // loop over track's plane references
	const TDetect& d = setup.iPlane2Detect(*ipl);
	const TPlane&  p = setup.vPlane(*ipl);

	chi2 = GetSmoothed(H,(*ipl));

	if (p.IFlag && d.InActive(H(1),H(2)))
	  cout<<"*"; // in active area 
	else
	  cout<<"-"; // in dead zone or outside or plane is OFF
	ipl++; ihp++;
      } // end of loop over track's plane references
      cout<<" : expected"<<endl;

      ipl = lPlnRef.begin();
      ihp = lHitPat.begin(); 
      while(ipl != lPlnRef.end()){ // loop over track's plane references
	if((*ihp) >=  0 ) cout<<"*"; // hit found
	if((*ihp) == -1 ) cout<<"o"; // no  hit found, while has to be there
	if((*ihp) <= -2 ) cout<<"-"; // not expected to be found (out of active area or plane is off)
	ipl++; ihp++;
      } // End of loop over track's plane references
      cout<<" :    found"<<endl;
      cout<<endl;
    } else {
      if(!insPlaneDone) cout<<"missed planes insertion had not been done"<<endl;
      if(!fFitDone)     cout<<"forward  Kalman fit had not been done"<<endl;
      if(!bFitDone)     cout<<"backward Kalman fit had not been done"<<endl;
    }
  }

  if (level==2 || level>=3) {     // *************** HIT INFO ***************

    list<int>::const_iterator ipr = lPlnRef.begin();
    list<int>::const_iterator ihp = lHitPat.begin();
    printf("          Plane           U (cm)   hit #   dChi2 %%  MC track        Time \n"); 

    for (; ipr!=lPlnRef.end(); ipr++, ihp++) {

      int ipl = *ipr;      // ********** LOOP ON PLANE REFERENCES **********
      //if (*ipr<57) continue; if (*ipr>66) break;
      // ***** CHI2 INCREMENTS *****
      char str5[8] = "  ---  ";
      map<int,double>::const_iterator ichi;
      ichi = mChi2.find(ipl);
      if (ichi!=mChi2.end()) sprintf(str5,"%7.1f", 100.*(*ichi).second/Chi2tot);

      //                                   ***** HIT INDICES *****
      char str6[9] = "--------"; // enough for 2 track-IDs + a "+" + margin
      const TDetect &d = setup.iPlane2Detect(ipl);
      if (d.IType==30) continue; // Special case of dummy MP (cf. TSetup::Init)
      const TPlane &p = setup.vPlane(ipl);
      float ang = 0.1*setup.vProj()[p.IProj];
      if (ev.IsMC() && *ihp>=0) {
	const THit &h = ev.vHit(*ihp);
	int ikin = h.IKine, ikin2 = h.IKin2;
	if (ikin>=0 || ikin2>=0) {
	  if      (ikin>=0 && ikin2>=0)
	    sprintf(str6,"%3d+%3d%c",ikin2,ikin, h.IOrig?'*':' ');
	  else if (ikin==-1 && ikin2>=0)
	    sprintf(str6," ??+%3d%c",      ikin2,h.IOrig?'*':' ');
	  else {
	    CsCluster *c = h.PtrClus(), *m = c->getMirrorCluster();
	    bool coalesce = false; if (m) {
	      double uc = c->getU()/10, um = m->getU()/10;
	      const TPlane  &p = setup.vPlane(h.IPlane);
	      const TDetect &d = setup.vDetect(p.IDetRef);
	      if (fabs(um-uc)<TOpt::dCut[27]*d.Resol)// Close enough to coalesce
		sprintf(str6,"mir&%3d%c",ikin,h.IOrig?'*':' ');
	    }
	    if (!coalesce)
	      sprintf(str6,"%4d%c",ikin,h.IOrig?'*':' ');
	  }
	}
	else if (ikin==-1) sprintf(str6,"bgr.");
	else sprintf(str6,"mir.%3d%c",-2-ikin,h.IOrig?'*':' ');
      }

      char strU[10] = "     ---", strV[10] = "     ---";
      char str7[10] = "     ---";       // ***** TIME *****
      char whichSide = ' ';
      if (*ihp>=0) {
	const THit &h = ev.vHit(*ihp);
	if (h.SigT>0) sprintf(str7,"%8.2f", h.Time);
	else          sprintf(str7,"%s","     ...");
	sprintf(strU,"%8.2f",h.U);
	if (d.IType==29 || d.IType==32) // ***** V-COORD of PIXELGEM/MMs *****
	  sprintf(strV,"%8.2f",h.V);
	//                                   ***** WHICH-SIDE INFO *****
	if (d.IType==11 || d.IType==12 || // Inner straws, except ST04 (waiting for X-check which-side sign consistency...)
	    d.IType==39 || d.IType==40) { // MAs
	  double side = h.PtrClus()->getDigitsList().front()->getData()[1];
	  if      (side>.5)  whichSide = d.Name[4]=='Y'? '-' : '+';
	  else if (side<-.5) whichSide = d.Name[4]=='Y'? '+' : '-';
	}
	else if (d.IType==28 || d.IType==31) {    // Stripped pixelGEM/MMs
	  const vector<double> &info = h.PtrClus()->getAllAnalogDataErrors();
	  if (info.size()<4) whichSide = '?';
	  else {
	    // Note that only those strips that cross the pixelised central
	    // piece yield an unambiguous info per se. For all others, the
	    // cluster might sit close to the border between the 2 sides. The
	    // unambiguous strips correspond to channels [88,167], cf.
	    // "THit::Print()". Note also that, in any case, part of the <0 side
	    // may correspond to a >0 V coordinate because of the detector's
	    // offset w.r.t. spectrometer axis.
	    if      (info[3]> 0.1)     whichSide = '+';
	    else if (info[3]<-0.1)     whichSide = '-';
	    else /* if == (double)0 */ whichSide = ' ';
	  }
	}
      }

      printf("%3u[%s] %6.1f  %s%c  %4d   %s   %9s     %s\n",
	     ipl,d.Name.c_str(),ang,strU,whichSide,(*ihp),str5,str6,str7);
      if (d.IType==29 || d.IType==32) // PixelGEM/MM
	printf("              %6.1f  %s\n",ang+90,strV);

      if (level>=3) {       // ***** PRINT "SMOOTHER" HELICES *****
	unsigned int doPrint = ~(level-3);
	map<int, THlx>::const_iterator im;
	if (doPrint&0x1) {
	  im = mHef.find(ipl);
	  if(im != mHef.end()){
	    (*im).second.Print("Hef");
	  } else {
	    cout<<"    no Hef for plane "<<ipl<<endl;
	  }
	}
	if (doPrint&0x2) {
	  im = mHeb.find(ipl);
	  if(im != mHeb.end()){
	    (*im).second.Print("Heb");
	  } else {
	    cout<<"    no Heb for plane "<<ipl<<endl;
	  }
	}
	if (doPrint&0x1) {
	  im = mHuf.find(ipl);
	  if(im != mHuf.end()){
	    (*im).second.Print("Huf");
	  } else {
	    cout<<"    no Huf for plane "<<ipl<<endl;
	  }
	}
	if (doPrint&0x2) {
	  im = mHub.find(ipl);
	  if(im != mHub.end()){ 
	    (*im).second.Print("Hub");
	  } else {
	    cout<<"    no Hub for plane "<<ipl<<endl;
	  }
	}
      }

    } // End of loop over plane references for hits

    if (TOpt::Print[0]&0x8) {   // ********** "SMOOTHED" HELICES **********
      list<THlx>::const_iterator ihsm;
      for (ihsm = lHsm.begin(); ihsm!=lHsm.end(); ihsm++)
	(*ihsm).Print("Smoothed helix");
    }
  } // End long print w/ hit info and smoothed helices
}

void TTrack::UpdateProjs()
{
  // Re-built set of projections (which TraFDic does not look after in the
  // course of track 's construction)
  const TSetup &setup = TSetup::Ref();
  sProj.clear();
  list<int>::const_iterator ipr = lPlnRef.begin();
  list<int>::const_iterator ihp = lHitPat.begin();
  for (; ipr!=lPlnRef.end(); ipr++, ihp++) {
    if (*ihp<0) continue;
    sProj.insert(setup.vPlane(*ipr).IProj);
  }
}
