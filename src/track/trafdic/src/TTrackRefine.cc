// $Id: TTrackRefine.cc 13148 2011-12-28 16:55:25Z kbicker $


#include <iostream>
#include <iomanip>
#include "TH1.h"
#define TTr_smthPulls_BINS 1  // Binning of smoothed pulls...
#if TTr_smthPulls_BINS == 1   // ... ==1 means TProfile
#  include "TProfile.h"
#else
#  include "TH2.h"
#endif
#include "CsHistograms.h"
#include "TOpt.h"
#include "TEv.h"
#include "TSetup.h"
#include "TTrack.h"
#include "THit.h"

using namespace std;

/*!
  \brief Hits Maps: Found and Expected (i.e. when in active area)...
  ...AND smoothed pulls
*/

/*
  Version specific to "lattice" alternative.
  Changes w/ respect to "traffic/TracksFit":
   i) Cancel all what concerns DCs.
  ii) No longer mu Wall 1.
 iii) Histogramming.
  iv) Add smoothed pulls
   v) PixelGEM/MMs. 
*/

bool TTrack::Refine(int iter)
{

  // Following algorithms require smoothed helices
  if(!fFitDone || !bFitDone){
    cout<<"TTrack::Refine ==> FullKF function must be called before. Track ID "<<Id<<endl;
    return(false);
  } 

  const TSetup &setup = TSetup::Ref();
  const TEv    &ev    = TEv::Ref();
  static CsHist1D *expctd[5], *found[5];
#if TTr_smthPulls_BINS == 1   // ... ==1 means TProfile
  static TProfile *smthPulls[5];
#else
  static TH2D *smthPulls[5];
#endif


  // ************************* INITIALISATION *************************
  static int book = (TOpt::Hist[1] ? 1 : 0) | (TOpt::Hist[3] ? 2 : 0),
    fill = book;
  if (book){
    book = 0;
    CsHistograms::SetCurrentPath("/Traffic/BitMaps");

    TH1D *det_types[4];
    //char hName[] = "Expected0", hTitle[] = "Expected hits  -  Z0";
    char hName[] = "smthPulls0";
    char hTitle[] = "#DeltaU(hit-smoothed track)/#sigmaU  -  Z0";

    int igr_mx = (int)setup.vIplFirst().size();
    for (int igr = 0; igr<igr_mx; igr++) {

      //                     ********** BOOKING BIT MAP HISTOs **********
      double Xmin = setup.vIplFirst()[igr]-.5, Xmax = setup.vIplLast()[igr]+.5; 
      int NbinsX = (int)(Xmax-Xmin);
      sprintf(hName, "Expected%d",igr);
      sprintf(hTitle,"Expected hits  -  Z%d",igr);
      expctd[igr] = new CsHist1D(hName,hTitle,NbinsX,Xmin,Xmax);
      sprintf(hName, "Found%d",igr);
      sprintf(hTitle,"Found hits  -  Z%d",igr);
      found[igr]  = new CsHist1D(hName,hTitle,NbinsX,Xmin,Xmax);
      sprintf(hName, "DetTyp%d",igr);
      sprintf(hTitle,"Det Types  -  Z%d",igr);
      det_types[igr] = new TH1D(hName,hTitle,NbinsX,Xmin,Xmax);
      int ipl0 = setup.vIplFirst()[igr], ipl;
      for (ipl = ipl0; ipl<=setup.vIplLast()[igr]; ipl++) {
	const TPlane &p = setup.vPlane(ipl);
	const TDetect &d = setup.vDetect(p.IDetRef);
	double type = 0;
	if      (d.Name.find("FI")==0) type = .1;
	else if (d.Name.find("SI")==0) type = .2;
	else if (d.Name.find("MP")==0) type = .25;
	else if (d.Name.find("MM")==0) type = .3;
	else if (d.Name.find("GP")==0) type = .25;
	else if (d.Name.find("GM")==0) type = .3;
	else if (d.Name.find("DC")==0) type = .4;
	else if (d.Name.find("DW")==0) type = .4;
	else if (d.Name.find("DR")==0) type = .4;
	else if (d.Name.find("ST")==0) {
	  if (d.Name.c_str()[7]=='b')  type = .5;
	  else                         type = .55;
	}
	else if (d.Name.find("P") ==0) type = .6;
	else if (d.Name.find("M") ==0) type = .7;
	else if (d.Name.find("H") ==0) type = .8;
	det_types[igr]->SetBinContent(ipl-ipl0+1,type);
      }

      if (TOpt::Hist[3]) {   // ********** BOOKING SMOOTHED PULLS **********
	sprintf(hName,"smthPulls%d",igr);
	sprintf(hTitle,"#DeltaU(hit-smoothed track)/#sigmaU  -  Z%d",igr);
#if TTr_smthPulls_BINS == 1   // ... ==1 means TProfile
	smthPulls[igr] = new TProfile(hName,hTitle,NbinsX,Xmin,Xmax,"S");
#else
	smthPulls[igr] = new TProfile(hName,hTitle,NbinsX,Xmin,Xmax,
				      TTr_smthPulls_BINS,-5,5);
#endif
      }
    }
    CsHistograms::SetCurrentPath("/");
  }  

  //          ***** FILL HITS MAPS and/or SMOOTHED PULLS?... *****
  int doFill = NDFs>8 && Chi2tot/(NDFs-5)<3 &&      // ***** ...HIGH QUALITY
    SigmaTime>0 &&
    fabs(MeanTime-ev.GetEventTime())/SigmaTime<3 &&         // ***** ... IN-TIME
    Pinv() ? fill : 0;                                      // ***** w/ MOMENTUM

  THlx H; double chi2;
  list<int>::iterator ipl = lPlnRef.begin(), ihp = lHitPat.begin(); 
  while (ipl!=lPlnRef.end()){
    //      ********** LOOP OVER TRACK'S PLANE REFERENCES **********
    const TDetect &d = setup.iPlane2Detect(*ipl);
    const TPlane &p = setup.vPlane(*ipl);

    //     ***** SET "NOT-IN-ACTIVE-AREA" and "PLANE-OFF" FLAGS *****

    chi2 = GetSmoothed(H,*ipl,false); // Get "smoothed" track parameters

    if (p.IFlag==0) { // Plane is OFF
      if (*ihp>=0){
	cout<<"TTrack::Refine() ==> Hit from switched off plane had been associated to the track! "<<endl;
	assert(false);
      }
      if (*ihp==-1) *ihp = -3;
    }
    if (!d.InActive(H(1),H(2))) {                    // ***** NOT IN ACTIVE AREA
      if (*ihp==-1)      // At the moment only "no-hit" flag -1 
	*ihp = -2;       // Changed to -2 (out of active area of the plane) 
    }
    else {                                               // ***** IN ACTIVE AREA
      if (doFill) {
	//                 ***** HIGH QUALITY, W/IN ACTIVE AREA... *****
	int igr; if (*ipl<setup.vIplFirst()[0]) igr = 4;
	else {
	  igr = 0; while (*ipl>setup.vIplLast()[igr]) igr++;
	}
	if (doFill&0x1) expctd[igr]->Fill(*ipl);   // Expected hits map
	if (*ihp>=0) {                       // ***** HIT FOUND on CURRENT PLANE
	  if (doFill&0x1) found[igr]->Fill(*ipl);  // Found hits map
	  if (doFill&0x2) {                        // Smoothed pulls
	    const THit &h = ev.vHit(*ihp);
	    smthPulls[igr]->Fill(*ipl,(h.U-H(1)*d.Ca-H(2)*d.Sa)/d.Resol);
	  }
	}
      }
    }

    // check
    if (*ihp>=0) {
      const THit& h = ev.vHit(*ihp);
      if (h.IPlane!=(*ipl)) {
	cout<<"TTrack::Refine() ==> Track --> Plane <-- hit  mismatch.  Det. name "<<d.Name<<endl;
      }
    }

    ipl++; ihp++;
  } // end of loop over track's plane references
 
  return true;
}
