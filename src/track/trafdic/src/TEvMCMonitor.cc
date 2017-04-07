// $Id: TEvMCMonitor.cc 14069 2015-09-17 20:44:46Z lsilva $

#include <iostream>
#include "TH1.h"
#include "TH2.h"
#include "TProfile.h"
#include "CsHistograms.h"
#include "CsMCUtils.h"
#include "CsMCTrkHit.h"
#include "CsGeom.h"
#include "TOpt.h"
#include "Traffic.h"
#include "TSetup.h"
#include "TEv.h"

using namespace std;
using namespace CLHEP;

/*! 
  \brief Traffic method for monitoring MC reconstruction.

  Probably obsolete: much better is achievable by analyzing mDSTs.
*/

#define M_p   .93827231

void TEv::MCMonitor()
{
  if(TOpt::Hist[1] == 0) return;
  if(! IsMC()) return;
  
  const TSetup& setup = TSetup::Ref();
  CsMCUtils utils;

  //----------------------------------------

  // ********** HISTOGRAMS **********

  static CsHist1D *tPMCall, *tPMCsel, *tPMCse2;   // MC distributions
  static CsHist1D *tPall, *tPrec, *tPbrg, *tPrbrg, *tPsrbrg, *tPfake, *tPfbrg;
  static CsHist1D *tPsrbr2;
  static CsHist1D *tNall, *tNfake, *tNsel, *tNzon1;
  static CsHist2D *tEff;
  static CsHist1D *tCorCl, *tRNhit, *tRMCh, *tTrChi, *tNcl;
  static TProfile *tdppv, *tdphv, *tdPpv, *tdppt, *tdppl, *tdppL;
  //#define DEBUG_MC_MONITOR
#ifdef DEBUG_MC_MONITOR
  static TH1D *dball, *dbsel;
#endif
  static CsHist1D *pulls_v[5], *pulls_t[5];

  static CsHist1D* res[250];
  static CsHist2D *used[250], *mused[250];
  static CsHist2D *unused[250], *unmused[250];

  // ***** Scattered muon *****
  static TH1D *tQ2MC, *tQ2sel, *tQ2rec, *tQ2s2l; static TH2D *tQ2r2c;
  static TH1D *txBMC, *txBsel, *txBrec, *txBs2l; static TH2D *txBr2c;
  static CsHist1D *tyBMC, *tyBsel, *tyBrec, *tyBs2l; static CsHist2D *tyBr2c;
  static TH1D *tthMC, *tthsel, *tthrec, *tths2l; static TH2D *tthr2c;

  // ********** OTHER STATIC INITIALISED UPON FIRST **********
  // ID of SCATTERED mu, or OTHERWISE SPECIAL PARTICLE
  static int scttrdMuID = 1;// ID of scattered mu in AROMA and COMGeanted PYTHIA
  static int incdntMuID = 0;// ID of incident mu in AROMA and PYTHIA

  static bool first = true;
  if(first){
    first = false;
    CsHistograms::SetCurrentPath("/Traffic/MCmonitor");
    tPMCall= new CsHist1D("tPMCall","All MC tracks",                 200,0,200);
    tPMCsel= new CsHist1D("tPMCsel","Reconstructible tracks",        200,0,200);
    tPMCse2= new CsHist1D("tPMCse2","Reconstructible SAS",           200,0,200);
    tPall  = new CsHist1D("tPall",  "All tracks",                    200,0,200);
    tPbrg  = new CsHist1D("tPbrg",  "All Bridged tracks",            200,0,200);
    tPrec  = new CsHist1D("tPrec",  "Genuine tracks",                200,0,200);
    tPrbrg = new CsHist1D("tPrbrg", "Genuine Bridged tracks",        200,0,200);
    tPsrbrg= new CsHist1D("tPsrbrg","Reconstructible Genuine tracks",200,0,200);
    tPsrbr2= new CsHist1D("tPsrbr2","Reconstructible Genuine SAS",   200,0,200);
    tPfake = new CsHist1D("tPfake" ,"Fake tracks",                   200,0,200);
    tPfbrg = new CsHist1D("tPfbrg" ,"Fake Bridged tracks",           200,0,200);

    tEff = new CsHist2D("tEff","Reconstruction #epsilon vs. Tracks in Zone 1",10,0,1+1.e-6,5,0,30);

    tNall  = new CsHist1D("tNall",  "Number of tracks with P in event",30,0,30);
    tNfake = new CsHist1D("tNfake", "Number of fakes  with P in event",30,0,30);
    tNsel  = new CsHist1D("tNsel",  "Reconstructible tracks per event",30,0,30);
    tNzon1 = new CsHist1D("tNzon1", "Tracks in zone 1 per event",      30,0,30);

    tNcl   = new CsHist1D("tNcl",   "N clusters / track", 100, 0., 100);
    tCorCl = new CsHist1D("tCorCl", "N cl. same / N cl" , 100, 0., 1.25);
    tRNhit = new CsHist1D("tRNhit", "N cl. / max N cl." , 100, 0., 1.25);
    tRMCh  = new CsHist1D("tRMCh",  "N cl. / N MC hits" , 100, 0., 1.25);

    tTrChi = new CsHist1D("tTrChi", "Track Chi2/NDF-5", 100, 0., 20.);

    tdppv = new TProfile("tdppv","#Deltap/p vs. 1/p  -  Vtx"  ,50,-2,2,"S");
    tdPpv = new TProfile("tdPpv","#Deltap/p vs. 1/p  -  Vtx"  ,50,-2,2,"S");
    tdphv = new TProfile("tdphv","#Deltap/p vs. #hits  -  Vtx",50,.5,100.5,"S");
    tdppt = new TProfile("tdppt","#Deltap/p vs. 1/p  -  Trk"  ,50,-2,2,"S");
    tdppl = new TProfile("tdppl","pT - pV vs. 1/p"  ,    50,-2,2,"S");
    tdppL = new TProfile("tdppL","pTMC - pVMC vs. 1/p"  ,50,-2,2,"S");
#ifdef DEBUG_MC_MONITOR
    dball = new TH1D("dball","All Spectro tracks",150,-100,200);
    dbsel = new TH1D("dbsel","Sel Spectro tracks",150,-100,200);
#endif
    
    pulls_v[0] = new CsHist1D("PullY_v",   "Y pull  -  Vtx",    40, -10., 10.);
    pulls_v[1] = new CsHist1D("PullZ_v",   "Z pull  -  Vtx",    40, -10., 10.);
    pulls_v[2] = new CsHist1D("PulldYdX_v","dY/dX pull  -  Vtx",40, -10., 10.);
    pulls_v[3] = new CsHist1D("PulldZdX_v","dZ/dX pull  -  Vtx",40, -10., 10.);
    pulls_v[4] = new CsHist1D("PullCop_v", "q/P pull  -  Vtx",  40, -10., 10.);

    pulls_t[0] = new CsHist1D("PullY_f",   "Y pull  -  Trk",    40, -10., 10.);
    pulls_t[1] = new CsHist1D("PullZ_f",   "Z pull  -  Trk",    40, -10., 10.);
    pulls_t[2] = new CsHist1D("PulldYdX_f","dY/dX pull  -  Trk",40, -10., 10.);
    pulls_t[3] = new CsHist1D("PulldZdX_f","dZ/dX pull  -  Trk",40, -10., 10.);
    pulls_t[4] = new CsHist1D("PullCop_f", "q/P pull  -  Trk",  40, -10., 10.);

    CsHistograms::SetCurrentPath("/");

    if (TOpt::Hist[1]&0x20) {
      // MC hit - Cluster residuals
      CsHistograms::SetCurrentPath("/Traffic/MCmonitor/resid");
      char nam[7]; char tit[60];
      char name[] = "tMUSMM1_999", title[] = "Unmused SMM1_999";
      for(int ipl=0; ipl < int(setup.vPlane().size()); ipl++){ // loop over det. planes
	if(ipl >= 250) continue;
	const TPlane&  pl = setup.vPlane(ipl);
	const TDetect&  d = setup.iPlane2Detect(ipl);
	sprintf(nam,"res%03u",ipl);
	sprintf(tit,"%s  ID %3u  prj. %i",d.Name.c_str(), d.IDet, pl.IProj);
	res[ipl] = new CsHist1D( nam, tit, 100, -10., 10. );

	sprintf(name,"tu%s_%d",d.PtrDet()->getName(),d.PtrDet()->getUnit());
	sprintf(title,"used %s_%d",d.PtrDet()->getName(),d.PtrDet()->getUnit());
	used[ipl] = new CsHist2D(name,title,25,-d.Siz(1)/2,+d.Siz(1)/2,
			                    25,-d.Siz(2)/2,+d.Siz(2)/2);
	sprintf(name,"tMu%s_%d",d.PtrDet()->getName(),d.PtrDet()->getUnit());
	sprintf(title,"mused %s_%d",d.PtrDet()->getName(),d.PtrDet()->getUnit());
	mused[ipl] = new CsHist2D(name,title,25,-d.Siz(1)/2,+d.Siz(1)/2,
		  	                     25,-d.Siz(2)/2,+d.Siz(2)/2);

	sprintf(name,"tU%s_%d",d.PtrDet()->getName(),d.PtrDet()->getUnit());
	sprintf(title,"Unused %s_%d",d.PtrDet()->getName(),d.PtrDet()->getUnit());
	unused[ipl] = new CsHist2D(name,title,25,-d.Siz(1)/2,+d.Siz(1)/2,
			                      25,-d.Siz(2)/2,+d.Siz(2)/2);
	sprintf(name,"tMU%s_%d",d.PtrDet()->getName(),d.PtrDet()->getUnit());
	sprintf(title,"Unmused %s_%d",d.PtrDet()->getName(),d.PtrDet()->getUnit());
	unmused[ipl] = new CsHist2D(name,title,25,-d.Siz(1)/2,+d.Siz(1)/2,
		  	                       25,-d.Siz(2)/2,+d.Siz(2)/2);
      }
      CsHistograms::SetCurrentPath("/");
    }

    if (TOpt::Hist[1]&4) {

      // ***** Scattered muon *****

      CsHistograms::SetCurrentPath("/Traffic/MCmonitor");
      double sq4_10; sq4_10 = sqrt(sqrt(10.0));
      double Q2bin; int bin; double Q2bins[24];
      for (bin = 23, Q2bin = 100/sq4_10; bin>=0; bin--) {
	Q2bins[bin] = Q2bin; Q2bin /= sq4_10;
      }
      tQ2MC  = new TH1D("tQ2MC", "Q^{2} - MC",                  23,Q2bins);
      tQ2sel = new TH1D("tQ2sel","Q^{2}_{MC} - Reconstructible",23,Q2bins);
      tQ2rec = new TH1D("tQ2rec","Q^{2}_{MC} - Reconstructed",  23,Q2bins);
      tQ2s2l = new TH1D("tQ2s2l","Q^{2}_{MC} - R2constructible",23,Q2bins);
      tQ2r2c = new TH2D("tQ2r2c","Q^{2}_{MC} - R2constructed",  23,Q2bins,4,-.5,3.5);
      double xBbin; double xBbins[26];
      for (bin = 25, xBbin =   1; bin>=0; bin--) {
	xBbins[bin] = xBbin; xBbin /= sq4_10;
      }
      txBMC  = new TH1D("txBMC", "x_{B} - MC",                  25,xBbins);
      txBsel = new TH1D("txBsel","x_{B}^{MC} - Reconstructible",25,xBbins);
      txBrec = new TH1D("txBrec","x_{B}^{MC} - Reconstructed",  25,xBbins);
      txBs2l = new TH1D("txBs2l","x_{B}^{MC} - R2constructible",25,xBbins);
      txBr2c = new TH2D("txBr2c","x_{B}^{MC} - R2constructed",  25,xBbins,4,-.5,3.5);
      tyBMC  = new CsHist1D("tyBMC", "y_{B} - MC",                  15,.0,.8);
      tyBsel = new CsHist1D("tyBsel","y_{B}^{MC} - Reconstructible",15,.0,.8);
      tyBrec = new CsHist1D("tyBrec","y_{B}^{MC} - Reconstructed",  15,.0,.8);
      tyBs2l = new CsHist1D("tyBs2l","y_{B}^{MC} - R2constructible",15,.0,.8);
      tyBr2c = new CsHist2D("tyBr2c","y_{B}^{MC} - R2constructed",  15,.0,.8,4,-.5,3.5);
      double thbin; double thbins[11];
      for (bin = 0, thbin = 1e-6; bin<=10; bin++) {
	thbins[bin] = thbin; thbin = (bin+1)*(bin+1)*(bin+1)*1e-5;
      }
      tthMC  = new TH1D("tthMC", "#Theta - MC",                  10,thbins);
      tthsel = new TH1D("tthsel","#Theta_{MC} - Reconstructible",10,thbins);
      tthrec = new TH1D("tthrec","#Theta_{MC} - Reconstructed",  10,thbins);
      tths2l = new TH1D("tths2l","#Theta_{MC} - R2constructible",10,thbins);
      tthr2c = new TH2D("tthr2c","#Theta_{MC} - R2constructed",  10,thbins,4,-.5,3.5);
      CsHistograms::SetCurrentPath("/");
    }

    // ***** RESET ID of INCIDENT and SCATTERED MUON... *****
    // ...(or any SPECIAL PARTICLE one want to MONITOR INDEPENDENTLY)
    if (TOpt::Hist[10]) incdntMuID = TOpt::Hist[10]-1;
    if (TOpt::Hist[11]) scttrdMuID = TOpt::Hist[11]-1;
  } // End of booking (and initialisation) block

  //----------------------------------------

  {
    // ******************** MOMENTUM RESOLUTION  and PULLS ********************

    // coord of prim. vertx. (assumed to be the first one)
    double vX = vVtxMC(0).V(0), vY = vVtxMC(0).V(1), vZ = vVtxMC(0).V(2);

    list<TTrack>::iterator it;
    for(it = listTrack.begin(); it != listTrack.end(); it++){
      // loop over TTracks
      TTrack& tr = *it;
      if (tr.Pinv() == 0) continue;                   // No momentum
      if (tr.Type==0x10) continue;                    // Beam
      if (tr.Type==0x1)  continue;                    // Fringe field
#ifdef DEBUG_MC_MONITOR
      dball->Fill(1/tr.Pinv());
#endif
      if (tr.IKine < 0)   continue;                   // No MC association
      if (tr.Chi2tot/(tr.NDFs-5)>10) continue;        // Bad chi2
      TKine &t = vecKine[tr.IKine];
      const TVtxMC &vtx = vVtxMC(t.IVtxOrig);
      if (fabs(vtx.V(0)-vX)>1) continue;              // Primary track
#ifdef DEBUG_MC_MONITOR
      dbsel->Fill(1/tr.Pinv());
#endif
      double Pmc  = 1/ t.Pinv();
      double vYp = t.P(1)/t.P(0), vZp = t.P(2)/t.P(0); 

      // ********** RECONSTRUCTED vs. GENERATED @ MC VERTEX **********
      THlx Hextr; Hextr(0) = vX;
      if (tr.Hfirst.Extrapolate(Hextr,true)) { // extrapolate to MC vertex
	double dpp   = (1/Hextr(5)-Pmc)/Pmc;
	double Ppull = (Hextr(5)-1/Pmc)/sqrt(Hextr(5,5));
	tdppv->Fill(1/Pmc,           dpp);
	tdphv->Fill(double(tr.NHits),dpp);
	tdPpv->Fill(1/Pmc,Ppull);

	pulls_v[0] ->Fill( (Hextr(1)-vY   ) / sqrt(Hextr(1,1)) );
	pulls_v[1] ->Fill( (Hextr(2)-vZ   ) / sqrt(Hextr(2,2)) );
	pulls_v[2] ->Fill( (Hextr(3)-vYp  ) / sqrt(Hextr(3,3)) );
	pulls_v[3] ->Fill( (Hextr(4)-vZp  ) / sqrt(Hextr(4,4)) );
	pulls_v[4] ->Fill( (Hextr(5)-1/Pmc) / sqrt(Hextr(5,5)) );
      }

      // ********** RECONSTRUCTED vs. GENERATED @ FIRST POINT **********
      const vector<int>& HitsMC = t.vHitMCRef();
      vector<int>::const_iterator iH;
      for (iH = HitsMC.begin(); iH!=HitsMC.end(); iH++) {
	if (fabs(vecHitMC[*iH].Xyz(0)-tr.Hfirst(0))<0.5 &&
	    vecHitMC[*iH].IOrig==0) {
	  double YMC = vecHitMC[*iH].Xyz(1), ZMC = vecHitMC[*iH].Xyz(2);
	  pulls_t[0]->Fill((tr.Hfirst(1)-YMC)/ sqrt(tr.Hfirst(1,1)) );
	  pulls_t[1]->Fill((tr.Hfirst(2)-ZMC)/ sqrt(tr.Hfirst(2,2)) );
	  const CsMCTrkHit* HitMC = dynamic_cast<const CsMCTrkHit*> (vecHitMC[ *iH ].PtrHit());
	  if (HitMC==NULL) continue;
	  pulls_t[2]->Fill((tr.Hfirst(3)-HitMC->getP().x()/HitMC->getP().z())/
			   sqrt(tr.Hfirst(3,3)));
	  pulls_t[3]->Fill((tr.Hfirst(4)-HitMC->getP().y()/HitMC->getP().z())/
			   sqrt(tr.Hfirst(4,4)));
	  pulls_t[4]->Fill((fabs(tr.Pinv())   -1/HitMC->getP().mag()        )/
			   sqrt(tr.Hfirst(5,5)));
	  tdppt->Fill(1/Pmc,(fabs(1/tr.Pinv())-HitMC->getP().mag())/fabs(Pmc));
	  tdppl->Fill(1/Pmc,(fabs(1/tr.Pinv())-fabs(1/Hextr(5)))/fabs(Pmc));
	  tdppL->Fill(1/Pmc,(HitMC->getP().mag()-fabs(Pmc))/fabs(Pmc));
	  break;
	}
      }
    }
  }

  if (setup.vIplFirst().size()<5) // The rest of the routine implies #zones>=5
    return;

  // *****************************************************************
  // ******************** EFFICIENCY HISTOGRAMS...********************
  // *****************************************************************

  static double Q2, xB, nu, yB, theta; // ***** scattered mu *****
  int n_selected = 0, n_reconstructed = 0, n_reconstructible = 0, n_zon1 = 0;
  int n_select2d = 0, n_reconstruct2d = 0;
  for (int ik = 0; ik<(int)vecKine.size(); ik++) {

    //         ********** LOOP over MC TRACKS **********

    TKine& t = vecKine[ik];
    //                             ***** SKIP INCIDENT PARTICLE (w/ Px<0) 
    // (Rely on it to have a negative momentum, rather than on the
    // incident particle ID retrieved from option, which could be wrong.)
    if (t.P(0)<0) continue;
    if (t.Q()==0) continue;     // ***** SKIP NEUTRALS
    if (t.isPileup()) continue; // ***** SKIP PILEUP/HALO
    
    tPMCall->Fill(fabs(1./t.Pinv())); // All PRIMARY CHARGED MC tracks
    //#define MCMon_DEBUG
#if defined DEBUG || defined MCMon_DEBUG
    unsigned int  sctype = 0; static int sctypes[8] = {0,0,0,0,0,0,0,0};
#endif

    int nh1(0), nh2(0), nh3(0), nh4(0);
    if(!t.vHitMCRef().empty()){
      double r = t.vHitRef().size()/t.vHitMCRef().size();
      tRMCh->Fill(r);
    }
    int ihr, iDetPrv;
#define MCMon_EXCLUDE_DC1ONLY_TRACKS
#ifdef MCMon_EXCLUDE_DC1ONLY_TRACKS
#warning TODO: Revisit MCMon_EXCLUDE_DC1ONLY_TRACKS
    int nhDC1 = 0;
#endif
    for (ihr = 0, iDetPrv = -1; ihr<int(t.vHitMCRef().size()); ihr++) {

      // ***** LOOP OVER CLUSTERS *****

      THitMC &hMC = vecHitMC[t.vHitMCRef()[ihr]];
      if (hMC.IKine!=ik)
	CsErrLog::msg(elError,__FILE__,__LINE__,
		      "MC track %d -> MC Hit (#%d->%d) association mismatch.",
		      ik,ihr,hMC.IKine);
      if (hMC.IOrig!=0) continue; // Skip hits from secondaries (delta rays,...)
      const TDetect &d = hMC.DetRef(); if (d.IDet==iDetPrv) continue;
      iDetPrv = d.IDet;
      const TPlane &p = setup.Id2Plane(hMC.IDet); int ipl = p.IPlane;
      if      (ipl<setup.vIplFirst()[0]) continue;
      else if (ipl<setup.vIplFirst()[1]) {
#ifdef MCMon_EXCLUDE_DC1ONLY_TRACKS
	if (d.IType==15) nhDC1++;
#endif
	/*  */                           nh1++;
      }
      else if (ipl<setup.vIplFirst()[2]) nh2++;
      else if (ipl<setup.vIplFirst()[3]) nh3++;
      else                               nh4++;
    } // End of loop over clusters
#ifdef MCMon_EXCLUDE_DC1ONLY_TRACKS
    if (nhDC1==nh1)                      nh1 = 0;
#endif

    int reconscmuctible = 0;    // Scattered muon reconstructibility
    bool scttrdMu = ik==scttrdMuID; if (scttrdMu) {

      // ***** SCATTERED mu (or otherwise SPECIAL PARTICLE) *****

      // Compute some kinematical variables...
      TKine &mui = vecKine[incdntMuID];     // Incident muon + M_PI!
      const HepLorentzVector &k = mui.ptrTrk->getP(), &kp = t.ptrTrk->getP();
      nu = k.e()-kp.e(); yB = nu/k.e();
      theta = -k.angle(kp)+M_PI;
      float stheta2 = sin(theta/2);
      Q2 = 4*k.e()*kp.e()*stheta2*stheta2; xB = Q2/2/M_p/nu;
      //printf("Q2,xB,y = %.4g %.4g %.3f\n",Q2,xB,yB);
      if (TOpt::Hist[1]&4) { // ... if in addition histo requested for sc mu
	tQ2MC->Fill(Q2); txBMC->Fill(xB); tyBMC->Fill(yB); tthMC->Fill(theta);
	if (nh2>=8) {
	  tQ2sel->Fill(Q2); txBsel->Fill(xB); tyBsel->Fill(yB); tthsel->Fill(theta);
	  reconscmuctible |= 0x1;
	}
	if (nh3>=4 && nh4>=TOpt::iCut[12]) {
	  tQ2s2l->Fill(Q2); txBs2l->Fill(xB); tyBs2l->Fill(yB); tths2l->Fill(theta);
	  reconscmuctible |= 0x2;
	}
      }
    }

    int nante = TOpt::iCut[1] ? TOpt::iCut[1] : 8;
    int npost = TOpt::iCut[2] ? TOpt::iCut[2] : 8;

    if (nh1>=6) n_zon1++;

    // ***** RECONSTRUCTIBILITY: ? HITS ANTE and ? POST SM1(2)... *****
    //                   ... or SCATTERED MUON (cf. supra)

    bool reconstructible, reconstructibl2;
    int nant2 = TOpt::iCut[3] ? TOpt::iCut[3] : 5;
    int npos2 = TOpt::iCut[4] ? TOpt::iCut[4] : 6;
    {
      static bool first = true;
      if (first) {
	CsErrLog::msg
	  (elInfo,__FILE__,__LINE__,
	   "TraFDic: LAS Reconstructibility == %d ANTE %d POST",nante,npost);
	CsErrLog::msg
	  (elInfo,__FILE__,__LINE__,
	   "TraFDic: SAS Reconstructibility == %d ANTE %d POST",nant2,npos2);
	first = false;
      }
      reconstructible = nh1>=nante && nh2>=npost;
      reconstructibl2 = nh2>=nant2 && nh3>=npos2;
    }

    if (reconstructible) {
      tPMCsel->Fill(fabs(1./t.Pinv()));
      if (fabs(1./t.Pinv())>1.0 && nh2>=8 && ik!=scttrdMuID) {
	n_reconstructible++; 
	Traffic::Ref().n_selected_MC++;         // 1-value efficiency estimate
	n_selected++;
      }
    }
    if (reconstructibl2) {
      tPMCse2->Fill(fabs(1./t.Pinv()));
      n_select2d++;
    }

    // look if this MC track had been reconstructed ?
    //if( t.sTrackID().empty() ) goto nextMCtrack;  // skip if not

    set<unsigned int, less<unsigned int> >::iterator iID;
    for(iID = t.sTrackID().begin(); iID != t.sTrackID().end(); iID++){

      // ***** LOOP OVER FOUND TRACK IDs *****

      list<TTrack>::iterator it;
      for(it = listTrack.begin(); it != listTrack.end(); it++){

	// loop over TTracks

	TTrack& tr = (*it);
	if(tr.Id != (*iID))    continue;  // not the ID we need
	tRNhit->Fill( double(tr.NHits)/double(nh1+nh2+nh3) ); // N hits / MC N hits ratio

	if(tr.Pinv() == 0)     continue;  // momentum is not measured

	// ********** ALL RECONSTRUCTED (w/ MOMENTUM) TRACKS **********

	tPrec->Fill(fabs(1./t.Pinv()));      // Fill with |P|, taken from TKine

#if defined DEBUG || defined MCMon_DEBUG
	if (scttrdMu)  sctype = tr.Type&0x1 | (tr.Type&0xc)>>1;
#endif
	if ((reconscmuctible&0x2) &&    // ***** SCATTERED mu IN SAS *****
	    tr.InGroup(1) && tr.InGroup(2)) {
	  int type = tr.Type&0x1 | (tr.Type&0x8)>>2;
	  tQ2r2c->Fill(Q2,(float)type); txBr2c->Fill(xB,   (float)type);
	  tyBr2c->Fill(yB,(float)type); tthr2c->Fill(theta,(float)type);
	}

	if (tr.NGroups()==1) continue;       // Exclude fringe reconstruction

	//         ***** BRIDGED TRACKS *****

	tPrbrg->Fill(fabs(1./t.Pinv()));

	//           ***** SA SPECTRO *****

	if (reconstructibl2 &&               // SA Sp Reconstructibility
	    tr.InGroup(2)) {
	    tPsrbr2->Fill(fabs(1./t.Pinv())); n_reconstruct2d++;
	}

	if (tr.InGroup(0)==0) continue;      // Exclude SA only tracks

	//             ***** LA SPECTRO *****

	if (reconstructible) {               // LA Sp Reconstructibility
	  tPsrbrg->Fill(fabs(1./t.Pinv()));
	  if (fabs(1./t.Pinv())>1.0 && nh2>=8 && ik!=scttrdMuID) { 
	    Traffic::Ref().n_reconstructed_MC++;  // Average Efficiency
	    n_reconstructed++;
	  }
	}
	if (reconscmuctible&0x1) {        // ***** SCATTERED mu *****
	  tQ2rec->Fill(Q2); txBrec->Fill(xB);
	  tyBrec->Fill(yB); tthrec->Fill(theta);
	}
	goto nextMCtrack; // take only the first corresponding TTrack (as best tracks are first)
      } // end of loop over TTracks
    } // end of loop over IDs
  nextMCtrack:;
#if defined DEBUG || defined MCMon_DEBUG
    if ((TOpt::Hist[1]&4) && scttrdMu) {
      if (sctype<0 || 7<sctype )
	CsErrLog::msg
	  (elFatal,__FILE__,__LINE__,"MCMon_DEBUG S µ type = %d!",sctype);
      sctypes[sctype]++;
      printf("Sµ  (0x1|0x4|0x8): 0x%x  (",sctype);
      for (int isct = 0; isct<8; isct++) printf("  %4d",sctypes[isct]);
      printf(")\n");
    }
#endif
  }// end of loop over MC tracks

  // ***** RECONSTRUCTION EFF. per EVENT *****

  tEff->Fill((float)n_reconstructed/n_selected,n_zon1+0.5);
  if (TOpt::Print[8]) {
    printf("Eff's (current, average, 2nd) %f %f   %f\n",
	   (float)n_reconstructed/n_selected,
	   (float)Traffic::Ref().n_reconstructed_MC/Traffic::Ref().n_selected_MC,
	   (float)n_reconstruct2d/n_select2d);
  }
  tNsel->Fill(n_reconstructible+0.5); tNzon1->Fill(n_zon1+0.5);

  //----------------------------------------

  //         ********** TRACK QUALITY HISTOGRAMS **********

  int Nall(0), Nfake(0);
  list<TTrack>::iterator it;
  for (it = listTrack.begin(); it!=listTrack.end(); it++) { // Loop on TTracks
    TTrack &t = *it;
    tNcl->Fill(double(t.NHits));
    double r = double(t.NHsame)/double(t.NHits);
    tCorCl->Fill(r);
    if (t.Pinv()!=0) {
      Nall++; tPall->Fill(fabs(1./t.Pinv()));
      if (t.NGroups()>1) tPbrg->Fill(fabs(1./t.Pinv()));
    }
    if (t.IKine>=0 &&  r<TOpt::dCut[10] &&  // Check RD<->MC association.
	(t.HasShower!=1 || !TOpt::dCut[81])) {// RICHWall shower is an exception
      int pid = vKine(t.IKine).PDGid();       // Check it's a decaying pion
      int nVtcs = vKine(t.IKine).vVtxMCRef().size();
      if (abs(pid)!=211 || nVtcs!=1) {
	CsErrLog::msg(elError,__FILE__,__LINE__,
	  "Event #%d Track #%d: %.1f%%(<cut=%.1f%%) hits from same MC (none of which from a RICHWall shower)...\n"
	  "  ...while associated to MC #%d w/ PDG ID = %d and #decay vertices = %d\n",
		      event,t.Id,r*100,TOpt::dCut[10]*100,t.IKine,pid,nVtcs);
      }
    }
    if (t.IKine<0) { // Considered to be fake track
      if (t.Pinv()!=0) {
	Nfake++; tPfake->Fill(fabs(1./t.Pinv()));
	if (t.NGroups()>1) tPfbrg->Fill(fabs(1./t.Pinv()));
      }
    }
    tTrChi->Fill(t.Chi2tot/t.NDFs);
  } // End of loop on tracks

  tNall ->Fill(Nall +0.5);
  tNfake->Fill(Nfake+0.5);

  //----------------------------------------

  // MC hit - Cluster residuals
  if (TOpt::Hist[1]&0x20) {
    set<int, less<int> >::iterator icl;
    for(int ipl=0; ipl < int(setup.vPlane().size()); ipl++){ // loop over det. planes
      if(ipl >= 250) continue;
      const TPlane&  pl = setup.vPlane(ipl);
      const TDetect&  d = setup.iPlane2Detect(ipl);
      for(int ih = 0; ih < int(pl.vHitMCRef().size()); ih++){ // loop over MC hits of the det. plane
	const THitMC& h = vHitMC(pl.vHitMCRef()[ih]);
	if(h.DetRef().IDet != d.IDet) {
	  cout<<"TEv::Monitor() ==> Mismatch: h.DetRef().IDet  = "<<h.DetRef().IDet<<" d.IDet = "<< d.IDet<<endl;
	}
	double uv=d.Ca*h.Xyz(1) + d.Sa*h.Xyz(2); // project MC hit to current det. WRS
	for(icl = h.sHitRef().begin(); icl != h.sHitRef().end(); icl++){ // loop over corresponding clusters
	  const THit& cl = vHit(*icl);
	  if( ! utils.isAGoodCluster(cl.PtrClus()) ) continue; // skip mirror hit
	  double dU = (cl.U - uv)/cl.SigU;
	  res[ipl]->Fill(dU);
	  if (cl.sTrackID().empty()) {
	    unused[ipl]->Fill(h.Xyz(1),h.Xyz(2));
	    if (cl.IKine==1 || cl.IKin2==1)
	      unmused[ipl]->Fill(h.Xyz(1),h.Xyz(2));
	  }
	  else {
	    used[ipl]->Fill(h.Xyz(1),h.Xyz(2));
	    if (cl.IKine==1 || cl.IKin2==1)
	      mused[ipl]->Fill(h.Xyz(1),h.Xyz(2));
	  }
	}
      }
    }
  } // MC hit - Cluster residuals
}
