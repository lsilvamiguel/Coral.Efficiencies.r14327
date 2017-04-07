// $Id: FindPrimary.cc 14069 2015-09-17 20:44:46Z lsilva $

/*!
   \file    FindPrimary.cc
   \brief   Compass Vertex Pattern Class.
   \author  Alexandre Korzenev
   \version $Revision: 14069 $
   \date    $Date: 2015-09-17 22:44:46 +0200 (Thu, 17 Sep 2015) $
*/

#include "CsAverPattern.h"
#include "CsGeant3.h"
#include "THlx.h"
#include "CsEvent.h"
#include "CsGeom.h"

using namespace std;
using namespace CLHEP;

extern void lindist(THlx *HTr1, THlx *HTr2, double xmin[], double *dist);

namespace {
  inline void Print1(int, TMtx&, THlx*, double, double);
}

bool CsAverPattern::FindPrimary( THlx** HTrMom, int* vIact, CsTrack** TrkRef, int tn )
{
  const int NTRK_MAX=32;
  double DCA[NTRK_MAX],ZCA[NTRK_MAX],COS[NTRK_MAX];
  
  int Schema = 0;
  int i, nTrksInV;
  double v12[3]={0,0,0},d12=0;
  TMtx Xn(3);
  double zT = CsGeom::Instance()->getTargetCenter()/10;// Center of target in cm

  nTrksInV = 1;     // Initialize "nTrksInV" w/ one count, for the beam track.
  for (i=0;i<NTRK_MAX;i++) {
    DCA[i]= 999999; ZCA[i]=-999999; COS[i]= 1;
    if (vIact[i]!=5) vIact[i] = 0;
  }
  CsEvent*  event  = CsEvent::Instance();
    
  //          ***** HISTOGRAM BOOKING *****
  static CsHist1F *hX, *hY, *hZ, *hDZ;
  static CsHist2F *hXY, *hZT, *hZTmm;
  static CsHist1S *hNTrk;
  static CsHist1F *hMeanTime, *hMeanTimeS, *hDiffTimeM;
  static CsHist1F *hHlxDist;
  static bool first = true; if (hist_ && first) {
    const double Zint = 5000;
    CsHistograms::SetCurrentPath("/CsAverPattern/Distribution");
    hX = new CsHist1F("PX","Prefilter2: X", 100,-40       ,40);
    hY = new CsHist1F("PY","Prefilter2: Y", 100,-40       ,40);
    hZ = new CsHist1F("PZ","Prefilter2: Z", 500,10*zT-Zint/2,10*zT+Zint );
    hZTmm = new CsHist2F( "hZTmm", "Z/Theta correlation of mu/mu' intersection", 100,10*zT-Zint,10*zT+Zint, 100,0,0.03 );
    hZT = new CsHist2F( "hZT", "Z/Theta correlation of intersection for all tracks", 100,10*zT-Zint,10*zT+Zint, 100,0,0.03 );
    hXY = new CsHist2F( "hXY", "XY distr. of prim. vertex", 100,-40,40, 100,-40,40 );
    hNTrk = new CsHist1S( "hNTrk", "N of tracks in prim. vertex before fit", 15,0,15 );
    hMeanTime  = new CsHist1F( "hMeanTime" , "Mean Time of tracks in prim. vertex", 50,-10,10 );
    hMeanTimeS = new CsHist1F( "hMeanTimeS", "Mean Time of tracks in prim. vertex (in \\sigma's)", 50,-5,5 );
    hDiffTimeM = new CsHist1F( "hDiffTimeM", "(t_{mu}-t_{mu'})/(\\sigma_{mu}^{2}+\\sigma_{mu'}^{2})^{1/2}", 50,-5,5 );
    hHlxDist = new CsHist1F("hHlxDist","Helix Distance (TimePrimCut*AcceptTill)",100,0,50);
    CsHistograms::SetCurrentPath("/");
    first = false;
  }

  double zBLast = // Abscissa of end of beam track
    const_cast<CsHelix&>(TrkRef[0]->getHelices()[1])(0)/10;

  bool downstreamTgt =  // ***** INTERACTION DOWNSTREAM OF PHYSICS TARGET *****
    //  If beam tracks ends up downstream of target...
    // ...These are tracks built by TraFDic (upon option), made up of a beam
    // telescope segment and a zone 0x1 segment, and interpreted as the track
    // of an incident particle which undergoes an interaction far downstream
    // of the physics target. They need a special processing.
    // Initial guess for vertex abscissa is no longer the target abscissa...
    zBLast>zT;

  if (Print_[0]) cout<<endl;
  double incMom = HTrMom[0]->Mom();
  for (i = 1; i<tn; i++) {
    //     ********** 1ST LOOP ON ALL OUTGOING TRACKS  **********

    if (Print_[0])
      cout<<"CsAP::FindPrim==> i = "<<i<<";   SPEC="<<vIact[i]<<endl;

    if (TrkRef[i]->hasMeanTime() && TrkRef[0]->hasMeanTime()) {
      //                                                          ***** TIME CUT
      double t   = TrkRef[i]->getMeanTime(), tb = TrkRef[0]->getMeanTime();
      double dT  = TrkRef[i]->getMeanTimeError();
      double dTb = TrkRef[0]->getMeanTimeError();
      double td = t-tb, dTd = sqrt(dT*dT+dTb*dTb);
      if (fabs(td)/dTd>TimePrimCut_) {
	if (vIact[i]==5) {
	  if (Print_[0])
	    printf("CsAP::FindPrim==> Mu' off time: %.3f Time cut: t/dT = %.1f (%.1f+/-%.1f - %.1f+/-%.1f) /%.2f > %.2f\n",1/(*HTrMom[i])(5),td,t,dT,tb,dTb,dTd,TimePrimCut_);
	  return false;
	}
	if (Print_[0])
	  printf("CsAP::FindPrim==> %.3f Time cut: t/dT = %.1f (%.1f+/-%.1f - %.1f+/-%.1f) /%.2f > %.2f\n",1/(*HTrMom[i])(5),td,t,dT,tb,dTb,dTd,TimePrimCut_);
	vIact[i] = 1; continue;
      }
    }

    //                                              ***** SMALL ENOUGH MOMENTUM?
    double mom = HTrMom[i]->Mom();
    double dp = sqrt((*HTrMom[i])(5,5))*mom*mom;
    if (vIact[i]!=5 && mom-3*dp>AcceptTill_*incMom) {
      // Note: No momentum cut on scattered mu!
      if (Print_[0])
	printf("CsAP::FindPrim==> Momentum cut: P-3dP=%.2f-3*%.2f > %.2f\n",
	       mom,dp,incMom);
      vIact[i] = 1; continue;
    }

    //                        ********** CDA **********
    //  I) The search for CDA is bounded (restricted to the gap between the beam
    // and 2ndary tracks (w/ some tolerance outside),
    // II) CDA is requisite.
    //  This has the drawback to suppress interactions in the very forward
    // direction (because then the uncertainty on CDA tends to infinite).
    //  E.g. low Q2 mu' or Primakoff. For the former, there would be a last
    // chance to recue the vertex: building it 1st on the 2ndary hadrons and
    // associating mu' to it a posteriori.
    //  It would be probably better not to put any requirement on CDA, but
    // rather to cut on the distance of the 2ndary track w.r.t. beam =
    //  - distance @ CDA if CDA turns out to lie w/in the gap,
    //  - minimum distance w/in gap otherwise.
    // Or, in other terms, define CDA as closest distance w/in gap. This, at
    // least for the scattered mu and in case there is no other vertex. For the
    // former case, alternatively, we could have an immediate 2nd step, mu'
    // dedicated, rescuing it if it's consistent w/ the vertex otherwise found.
    THlx Dum1 = *HTrMom[i], Dum2 = *HTrMom[0]; double s3;
    bool ok; if (downstreamTgt) {
      // ***** SPECIAL CASE of an INTERACTION DOWNSTREAM OF PHYSICS TARGET *****
      // Initial guess for CDA is downstream end of beam track. This sets also
      // the upstream limit of the CDA search, w/ some tolerance allowing the (
      // 1st estimate of the) vertex to actually lie w/in track. The tolerance
      // on the downstream limit is relaxed w.r.t. standard case. These loose
      // tolerances because tracking may not have been able (it's not designed
      // for) to tell the incident track from the outgoing one. 
      ok = Dum1.FindCDA(Dum2,zBLast,zBLast-30,Dum1(0)+20);
    }
    else
      // ***** STANDARD CASE *****
      // Initial guess for CDA is target center. The allowed range goes from end
      // of beam track to beginning of spectro track. W/ some margin allowing
      // CDA to actually fall w/in the span of either of these 2 tracks.
      ok = Dum1.FindCDA(Dum2,zT,zBLast-20,Dum1(0)+20);
    if (!ok) {
      // Look for closest distance approach up to abscissa of 1st helix of
      // scattered track, i.e. do not allow the latter to have clusters upstream
      // of CDA. Still allow for 10 cm margin to account for track imprecision.
      if (vIact[i]==5) {
	if (Print_[0])
	  printf("CsAP::FindPrim==> Mu': No intersection -> give up!\n");
	return false;
      }
      if (Print_[0]) printf("CsAP::FindPrim==> No intersection!\n");
      d12 = 1e10; s3 = 0; vIact[i] = 1;
    }
    else {
      d12 = Dum1.Dist( Dum2 );
      v12[0] = ( Dum1(0) + Dum2(0) ) / 2;  // Z
      v12[1] = ( Dum1(1) + Dum2(1) ) / 2;  // X
      v12[2] = ( Dum1(2) + Dum2(2) ) / 2;  // Y
      ZCA[i] = v12[0];

      // Error on distance between two helices
      if (DDistUseMMap_) { // (upon option) Take material into account in track uncertainties
	Dum2(0) = v12[0]; HTrMom[0]->Extrapolate(Dum2);
	Dum1(0) = v12[0]; HTrMom[i]->Extrapolate(Dum1);
      }		
      s3 = HelDist_ / d12 * sqrt( (Dum1(1)-Dum2(1))*(Dum1(1)-Dum2(1))*(Dum1(1,1)+Dum2(1,1)) +
				  (Dum1(2)-Dum2(2))*(Dum1(2)-Dum2(2))*(Dum1(2,2)+Dum2(2,2)) +
				2*(Dum1(1)-Dum2(1))*(Dum1(2)-Dum2(2))*(Dum1(1,2)+Dum2(1,2)));
      DCA[i] = d12/s3;
      COS[i] = Dum1.DirCos(1)*Dum2.DirCos(1)+Dum1.DirCos(2)*Dum2.DirCos(2)+Dum1.DirCos(3)*Dum2.DirCos(3);

      if (hist_) hHlxDist->Fill(d12);

      if (Print_[0])
	printf("CsAP::FindPrim==> Helixdist: z=%5.1g, mom=%#5.3g; dist=%#5.3g, cut:%#5.3g\n",v12[0],1/(*HTrMom[i])(5),d12,s3);
    }
  }
  if (Print_[0]) printf("\n");

  //  ***** E CONSERVATION *****
  // The aim is 2-fold:
  // - Prevent crude violation of the energy conservation by discarding tracks
  //  w/ a momentum larger than the (upper bound on the) missing momentum of the
  //  scattered muon. These are most probably accidentally coincident beam or
  //  near halo tracks unduly associated to the pVertex (for the current
  //  particular choice of scattered muon).
  // - Take into account energy (or rather longitudinal P) conservation in the
  //  evaluation of pVertices (which is based on #tracks and the crude estimate
  //  of chi2 available at the PR step) by giving a malus to those that violate
  //  conservation.
  double dp0 = sqrt((*HTrMom[0])(5,5))*incMom*incMom; // Uncertainty on the incident momentum
  double upperMom = -1; // Upper bound on momentum, given inc. and sc. tracks 
  for (i = 1; i<tn; i++) {

    //           ***** CUT on TRANSVERSE DISTANCE and Z *****

    if (vIact[i]!=0 && vIact[i]!=5) continue;

    if (DCA[i]<1) {                                      // ***** TRANSVERSE CUT
      // Note: No Z cut (w.r.t. target) is applied. Neither on the scattered mu,
      // because Z can be very imprecise (even infinitely imprecise) at low Q2,
      // nor on any other track, because then vertices far from target center
      // would be eliminated. This has a drawback: tracks unrelated to the
      // pVertex can be unduly associated to it. The work-around it is to
      // exclude tracks w/ Z of CA far from target center in the guessing the
      // initial value of Z of vertex (if other tracks, closer to target exist).
      if (Print_[0])	
	printf("CsAP::FindPrim==> Accepted. mom=%#7.4g, ZCA=%#6.3g,  SPEC=%d\n",
	       HTrMom[i]->Mom(),ZCA[i],vIact[i]);
    }
    else {  
      if (vIact[i]==5) {
	if (Print_[0])
	  printf("CsAP::FindPrim==> Mu' rejected (P = %.2f) -> Give up!\n",
		 HTrMom[i]->Mom());
	return false;
      }
      vIact[i] = 1;
      if (Print_[0])
	printf("CsAP::FindPrim==> Rejected, P = %.2f\n",HTrMom[i]->Mom());
    }
    if (vIact[i]==5) {   //  ***** PREPARE FOR E CONSERVATION *****
      // Let's have 4sigmas for the incident momentum. In order to make for the
      // case where the indident momentum is not determined by a BMS reco and
      // its uncertainty is equated to the beam spread divided by sqrt(12). 
      double mom = HTrMom[i]->Mom(), dp = sqrt((*HTrMom[i])(5,5))*mom*mom;
      upperMom = incMom+4*dp0-mom+3*dp;
    }
  }

  double sD2s, sPz, sDP2;// Evaluation of pVertices (for later selecting best...
  // ...pVertex in view of determining the event time before refitting). "sD2s"
  // collects the sum of squared distances, "sPz" the total z momentum and
  // "sDP2" the uncertainty on total momentum (easier to compute than, and to
  // serve as an estimate of, the uncertainty on the total z momentum).
  for (i = 1, sD2s = 0, sPz = -HTrMom[0]->DirCos(1)*HTrMom[0]->Mom(),
	 sDP2 = dp0*dp0; i<tn; i++) {
    if (vIact[i]!=0 && vIact[i]!=5) continue;
    double mom = HTrMom[i]->Mom();              // ***** ENERGY CONSERVATION CUT
    double dp = sqrt((*HTrMom[i])(5,5))*mom*mom;
    if (vIact[i]!=5 && upperMom>0 && mom-3*dp>upperMom) {
      printf("CsAP::FindPrim==> %.3f E cut: p-3dp(=%.2f) > p_max(=%.2f)\n",
	     1/(*HTrMom[i])(5),mom-3*dp,upperMom);
      if (Print_[0]) cout<<"CsAP::FindPrim==> Rejected. mom="<<HTrMom[i]->Mom()
			 <<",  SPEC="<<vIact[i]<<endl;
      vIact[i] = 1; continue;
    }
    nTrksInV++; sD2s += DCA[i]; sPz += HTrMom[i]->DirCos(1)*mom;
    sDP2 += (*HTrMom[i])(5,5)*mom*mom*mom*mom;
    if (hist_) {
      hZT->Fill(10*ZCA[i],acos(COS[i]));
      if (TrkRef[i]->hasMeanTime()) {
	double time = TrkRef[i]->getMeanTime();
	double timeERR = TrkRef[i]->getMeanTimeError();
	hMeanTime->Fill(time); hMeanTimeS->Fill(time/timeERR);
      }
    }
  }

  //   ***** DETERMINE THE INITIAL VERTEX POSITION: 2 SCHEMAS *****
  // - mu' takes precedence, depending upon scattering angle:
  //  0) Large: Some average Z of CDA, cf. infra.
  //  1) Small: Z of CDA of mu'
  // - No mu' <-> Schema #0.
  double zAv;    int nAvs;   // Average Z
  double zGuess; int nGuess; // Guess for Z of vertex (in schema #0): excluding tracks far from target center, at low angle, not bridged over SM1
  for (i = 1, zAv=zGuess = 0, nAvs=nGuess = 0; i<tn; i++) {
    // Note: Why aren't these blocks included in the loops on all outgoing tracks(already 2!) supra?
    if (vIact[i]!=0 && vIact[i]!=5) continue;
    double d12 = DCA[i], zca = ZCA[i], cca = COS[i];

    //    ***** LAST SELECTION CUT: TRY and REMOVE TERTIARIES *****
    // - This is a quick and dirty hack for rejecting some tertiary tracks...
    // ...and preventing them from upsetting the vertex.
    // - I(Y.B.) had understood that the procedure in "CsKF::doFitting" that
    //  filters out tracks incompatible w/ the vertex in the making was good
    //  enough, able not only to identify (worst chi2 increment) and reject fake
    //  secondaries, but also to put the vertex fit back on track afterward.
    // - This turns out not to be true: a fake track, w/ Z of CA far from Z of
    //  true vertex, can upset definitively the vertex fit: e.g.
    //  "/castor/cern.ch/user/g/guskov/Primakoff_t55/zebradat.0.fz", Evt #3.
    // - In order to fix these cases, the best would be to revisit 
    //  "CsKF::doFitting" in order precisely to put the fit back on track.
    // - But a more complete reshuffling of the method is otherwise needed, that
    //  would also deal w/ the case of re-interactions.
    // => Therefore, in the meanwhile, let's try and reject far-fetched
    //  tertiaries here in "CsAP::FindPrimary".
    // - This is achieved here based on the diff between suspected track's
    //  CA and average CA of all other tracks.
    // - The cut is very loose and should have no impact on cases where there
    //  is a re-interaction close to the primary one. It rather copes w/ the
    //  cases of a V0 decay or a gamma conversion.
    int j, tne; double ZCAe; // Mean Z of CA, bounded by Z of target
    for (j = 1, tne = 0, ZCAe = 0; j<tn; j++)
      if (j!=i && (vIact[j]==0 || vIact[j]==5)) { ZCAe += ZCA[j]; tne++; }
    if (tne) ZCAe /= tne;
    if (vIact[i]!=5 /* Exclude mu' */ && tne && fabs(ZCA[i]-ZCAe)>300) {
      double dp[4] = {.01,.01,.0001,.0001}, dZdp[4], dZ; int ip; bool ok;
      for (ip = 0, ok = true; ip<4; ip++) {
	THlx Dum1 = *HTrMom[i], Dum2 = *HTrMom[0]; Dum1(ip+1) += dp[ip];
	ok &= Dum1.FindCDA(Dum2,ZCA[i],ZCA[i]-100,ZCA[i]+100);
	if (ok) dZdp[ip] = ((Dum1(0)+Dum2(0))/2-ZCA[i])/dp[ip];
      }
      if (ok) {
	THlx &Dum1 = *HTrMom[i];
	for (ip = 0, dZ = 0; ip<4; ip++) {
	  int ip1 = ip+1; dZ += dZdp[ip]*dZdp[ip]*Dum1(ip1,ip1);
	  for (int jp = ip1; jp<4; jp++)
	    dZ += 2*dZdp[ip]*dZdp[jp]*Dum1(ip1,jp+1);
	}
	dZ = sqrt(dZ);
	if (ZCA[i]-ZCAe-300>6*dZ) {
	  vIact[i] = 1;
	  if (Print_[0])
	    printf("CsAP::FindPrim==> Rejected. mom=%#7.4g,"
		   " ZCAi:%#6.3g+/-%#6.3g << |ZCA|j!=i:%.2f\n",
		   HTrMom[i]->Mom(),ZCA[i],dZ,ZCAe);
	  continue;
	}
      }
    }

    if (vIact[i]==5) {                                    // ********** CASE MU'
      if (d12>1)  // Note: Never happens!! If "d12>1", "vIact" is set=1 supra!
	CsErrLog::msg(elFatal,__FILE__,__LINE__,
		      "CsAP::FindPrim==> Inconsistency: mu' accepted,"
		      " yet too far from bram (dist/cut=%.2f)\n",d12);
      else {                                // **** MU': EITHER OF THE 2 SCHEMAS
	if (cca>CosPrimMM_)
	  Schema = 0;
	else {
	  Schema = 1; Xn(1) = zca;
	}
	if(Print_[0]) cout<<endl<<"CsAP::FindPrim==> Angle(mu,mu') = "<<acos(cca)*180/3.141593<<endl;
	
	if (hist_) hZTmm->Fill(10*zca,fabs(acos(cca)));
      }
      if (hist_ && TrkRef[i]->hasMeanTime() ) {
	double time_B    = TrkRef[0]->getMeanTime();
	double timeERR_B = TrkRef[0]->getMeanTimeError();
	double time_M    = TrkRef[i]->getMeanTime();
	double timeERR_M = TrkRef[i]->getMeanTimeError();
	hDiffTimeM->Fill( (time_B-time_M) / sqrt( timeERR_B*timeERR_B - timeERR_M*timeERR_M ) );
      }
    } // End case mu'
    if (nTrksInV==2) // This correspond to case where only 1 2ndary...
      // - The following does not then override "Schema==1", but confirms it.
      // - It supplies a definite "Xn(1)" to vertices w/ single, non LAS, track.
      Xn(1) = zca; 
    if ((*HTrMom[i])(0)<350 &&   // If not a pure SAS far from pVertex track...
	cca<CosPrim_ &&          // ...at reasonably large angle
	fabs(zT-zca)<LinDist_) { // ...at reasonable distance from target center
      zGuess += zca; nGuess++;
    }
    zAv += zca; nAvs++;
  } // End loop on secondary tracks.
  if (!Schema) {
    if (nGuess) Xn(1) = zGuess/nGuess;
    else        Xn(1) = zAv/nAvs;
  }

  if (nTrksInV==1) {  // ***** UNIQUE TRACK (i.e. ONLY BEAM) => GIVE-UP *****
    statistics_[4]++;  // for statistics 
    return false;
  }

  int count; for (i=count = 0; i<tn; i++)
    if (vIact[i]==0 || vIact[i]==5) count ++;
  if (Print_[0]) cout<<"CsAP::FindPrim==> _______ End of loop over "<<tn<<
                       "("<<count<<")"<<" tracks _______\n";
  if (count<2) return false;

  //     ***** EXTRAPOLATE to FIRST GUESS for VERTEX' Z *****
  THlx hDum; HTrMom[0]->Extrap(Xn(1),hDum); Xn(2) = hDum(1); Xn(3) = hDum(2);

  //        ********** INSTANTIATE and STORE VERTEX **********
  CsVertex *vrt = new CsVertex( 10 * Xn(2), 10 * Xn(3), 10 * Xn(1) );
  HepMatrix &Cn = *new HepMatrix(3,3);
  Cn[0][0] = CnSigmas_[0];  Cn[0][1] = 0;  Cn[0][2] = 0; // sigma=0.2cm
  Cn[1][0] = 0;  Cn[1][1] = CnSigmas_[1];  Cn[1][2] = 0; // sigma=0.2cm
  Cn[2][0] = 0;  Cn[2][1] = 0;  Cn[2][2] = CnSigmas_[2]; // sigma=5cm
  vrt-> addCov(&Cn); vrt-> setType(true);
  //   ***** TEMPORARY VALUE FOR CHI2 *****
  // - Will be used to select best one among pVertices that have same #tracks.
  // - Single out those pVertices that violate momentum conservation and
  //  tag them by giving them a <0 chi2.
  if (sPz>4*sqrt(sDP2) /* 4-sigma cut */) sD2s *= -1; vrt-> setChi2(sD2s);
  //   ***** LIST of TRACKS: specials come 1st, w/ incident track very 1st.
  for (i = 0; i<tn; i++) if (vIact[i]==5) vrt->addTrack( TrkRef[i] );
  for (i = 0; i<tn; i++) if (vIact[i]==0) vrt->addTrack( TrkRef[i] );
  vrts_.push_back( vrt );

  if (Print_[0]) Print1(nTrksInV,Xn,HTrMom[0],TrkRef[0]->getMeanTime(),sD2s);
  statistics_[8] += vrt->getNTracks(); // for statistics
   
  if (hist_) {  // ********** HISTOGRAMING **********
    hX->Fill(10*Xn(2)); hY->Fill(10*Xn(3)); hZ->Fill(10*Xn(1));
    hXY->Fill(10*Xn(2),10*Xn(3));
    int NTRK = vrt->getNTracks(); hNTrk->Fill( NTRK );

    if (event->isAMonteCarloEvent()) {  // ***** MC HISTOS *****
      list<CsMCVertex*>mcvertices = CsGeant3::Instance()->getMCVertices();
      list<CsMCTrack*> mctracks   = CsGeant3::Instance()->getMCTracks();
      double XXMC = mcvertices.front()->getX();
      double XYMC = mcvertices.front()->getY();
      double XZMC = mcvertices.front()->getZ();
      hDeltaPref[0]->Fill(XXMC-10*Xn(2));
      hDeltaPref[1]->Fill(XYMC-10*Xn(3));
      hDeltaPref[2]->Fill(XZMC-10*Xn(1));      
    }
  }
  return true;
}


////////////// PRINTS. ARE LOCAL FOR THIS FILE. /////////////////////

namespace {

  inline void Print1(int nTracks, TMtx& Xn, THlx* HTrMom, double time, double chi2) {
    printf("\nCsAP::FindPrim==> pVertex %.2f ns, z,x,y = %.2f,%.2f,%.2f, %d Tracks, chi2 = %.2f\n",time,Xn(1),Xn(2),Xn(3),nTracks,chi2);
    if (CsEvent::Instance()->isAMonteCarloEvent()) {
      CsGeant3* geant3 = CsGeant3::Instance();
      list<CsMCVertex*>mcvertices = geant3->getMCVertices();
      list<CsMCTrack*> mctracks   = geant3->getMCTracks();
      double XXMC=mcvertices.front()->getX()/10;
      double XYMC=mcvertices.front()->getY()/10;
      double XZMC=mcvertices.front()->getZ()/10;
      cout<<"CsAP::FindPrim==> MC:         z="<<XZMC  <<", x="<<XXMC  <<", y="<<XYMC  <<endl;
      cout<<"CsAP::FindPrim==> MC-PATTERN: z="<<setw(8)<<XZMC-Xn(1)<<", x="
	  <<setw(8)<<XXMC-Xn(2)<<", y="<<setw(8)<<XYMC-Xn(3)<<endl;
      cout<<endl;
      
      list<CsMCTrack*>::iterator itrk;
      for(itrk=mctracks.begin(); itrk!=mctracks.end(); itrk++){
	if( (*itrk)->getGnum() == 1 ){ 
	  cout<<"CsAP::FindPrim==> E_MC-E_BEAM = "
	      <<( (*itrk)->getE() - fabs( 1/(*HTrMom)(5) ) )<<endl;
	  break;      
	}
      }
      cout<<endl;
    }
  }
}
