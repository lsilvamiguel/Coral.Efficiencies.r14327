/*!
  \file    Kalman.cc
  \brief   Kalman Filter w/ elimination of worst chi2 contributors.
  \author  A.Korzenev, C.Ulvegren 
  \version $Revision: 14069 $
  \date    $Date: 2015-09-17 22:44:46 +0200 (Thu, 17 Sep 2015) $

  - Elimnation is of all contributing tracks w/ KF-chi2 > "CsKalmanFitting CUTS" and of worst contributors until total chi2 < "CsKalmanFitting Chi2VCut".
  - The elination procedure is not fully satisfactory, because the global fit can be biased because of inital guess too far away from true value. Several rescue procedures are tried: still not fully satisfactory.
  
*/

//#define Kalman_DEBUG_ELoss 1
#ifdef Kalman_DEBUG_ELoss
#include "TH2.h"
#endif
#include "TMath.h"
#include "coral_config.h"
#include "CsErrLog.h"
#include "CsKalmanFitting.h"
#include "CsVTrack.h"
#include "CsGeom.h"
#include "CsEvent.h"
#include "CsGeant3.h"
#include "CsOpt.h"
#include "CsMCTrack.h"
#include "CsGeom.h"
#include "CsStopwatch.h"

using namespace std;

TMtx inv3x3(TMtx* M, int& ierr);
extern void thesort(double a[], int in[], int n, int imod);
extern void helixdist(THlx *H1, THlx *H2, double a,double b,double x0,
	       double ver[], double *dist, THlx *H1v, THlx *H2v);

namespace {
  inline void Print1( int );
  inline void Print2( list<CsVTrack>& );
  inline void Print3( list<CsVTrack>&,double);
  inline void Print4( CsVTrack &, double &, TMtx &, const char *);
  inline void Print6( double, double*, int );
  inline void Print7( TMtx &, TMtx &, int, int, double );
  inline void Print8( TMtx &, TMtx &, int, int, double );
}


void CsKalmanFitting::Kalman(TMtx &Xn,TMtx &Cn,list<CsVTrack> &vtrks,
			     double &chi2)
{
  int iact = vtrks.size(); if (iact<2) return;
  const int NTRK_MAX = 32; int iend = 1, Niter = 0 /* # filtering steps */;
  double Chi2trk[NTRK_MAX], Chi2tot=-1;
  TMtx Xnew(3), Cnew(3,3), C0(Cn), X0(Xn);
  list<CsVTrack>::iterator iTr;
  if (Print_[0]) Print1(iact);

  if (Schema_ == 0) {      //******** INVERSE KALMAN FILTER *********
    if (Print_[1]) Print2( vtrks );
    Globalfit(Cn,Xn,vtrks); Smoother(C0,Cn,X0,Xn,Chi2tot,vtrks);
    if (Print_[1]) { CovMatrix(Cn,vtrks); Print3(vtrks,Chi2tot/(2*iact-3)); }
    if (Print_[0] || Print_[1])
      cout<<endl<<"CsKF::Kalman==> GLOBAL FIT:  #tracks ="<<iact
	  <<" chi2="<<Chi2tot<<" chi2/NDF="<<Chi2tot/(2*iact-3)
	  <<", x,y,z="<<Xn(1)<<","<<Xn(2)<<","<<Xn(3)<<endl<<endl;

    for (; iact>1 && iend!=0 ; Niter++) {

      int it, ntNS, TrkPntr[NTRK_MAX];
      for (int i = 0; i<(int)vtrks.size(); i++) TrkPntr[i] = 0;
      for (iTr = vtrks.begin(), it=ntNS = 0; iTr!=vtrks.end(); iTr++, it++) {

	// ***** LOOP on CsVTrack's ***** 

	if (!(*iTr).primary()) {            // ***** TRACK ALREADY EXCLUDED: SKIP
	  if (Print_[2]) printf("CsKF::Kalman==> Trk %2d\n",it); continue;
	}
	else if (Print_[2]) printf("CsKF::Kalman==> Trk %2d",it);
	InvFilter(Cn,Xn,*iTr,Chi2trk[ntNS],&Xnew,&Cnew);     // ***** INVERSE KF
	if (Chi2trk[ntNS]<0) Chi2trk[ntNS] = 1000;
	if (Print_[2]) Print4(*iTr,Chi2trk[ntNS],Xnew," IKF");

	if ((*iTr).getStatus()!=5)          // ***** IF Non-Special CsVTrack...
	  TrkPntr[ntNS++] = it;                 //     ...STORE REFERENCE
      }
      if (ntNS==0) break;

      thesort(Chi2trk,TrkPntr,ntNS,1);      // ***** SORT ACCORDING to Chi2trk
      int vNDF = 2*iact-3; float vChi2 = Chi2tot/vNDF;
      if (Print_[2]) Print6(vChi2,Chi2trk,ntNS+1);

      bool vChi2Out = vChi2>Chi2VCut_;
      if ((Chi2trk[ntNS-1]>InvChi2Trk_) ||  // ***** IF WORST'S Chi2trk >>...
      	  vChi2Out) {                       //  ...OR TOTAL CHI2 >>...
	if ((Print_[2]) && Chi2trk[ntNS-1]<=InvChi2Trk_)
	  cout<<"CsKF::Kalman==> vChi2 "<<vChi2<<">"<<Chi2VCut_<<endl;
	bool ok = false;
#define Kalman_ReFIT_WNW  // Refit w/ worse and w/o next to worst
	// This is found to improve (very) marginally the reco efficiency for
	// D0 (80.21->80.26) D* (69.03->69.10) Vertexing (95.84->95.89)
	// while alleviating, somewhat, the vertex fakes rate (15.80->15.79)
	// (All measured on AROMA D* w/ 2002 setup and ideal reco)
	// => And is therefore enabled
#ifdef Kalman_ReFIT_WNW
	if (ntNS>1 && Chi2trk[ntNS-1]<Chi2trk[ntNS-2]*2 &&
	    (ntNS==2 || Chi2trk[ntNS-3]<InvChi2Trk_)) {
	  // If Next to Worst pretty bad AND NNWorst (if exists) OK
	  // Reevaluate thoroughly the 2 cases (viz. Worst and NWorst)
	  double chi2s[2] = {0,0}; int iOut, iter;
	  for (iter=iOut = 0; iter<=2; iter++) {
	    for (it = 0, iTr = vtrks.begin(); it<TrkPntr[ntNS-1-iOut]; it++)
	      iTr++;
	    (*iTr).setStatus(2);
	    InvFilter(Cn,Xn,*iTr,chi2s[iOut],&Xnew,&Cnew); // ...inverse KF
	    resetTrks(Xnew,vtrks); Globalfit(Cnew,Xnew,vtrks);
	    Smoother(C0,Cnew,X0,Xnew,chi2s[iOut],vtrks);
	    if (Print_[2]) {
	      cout<<"\nCsKF::Kalman==> Refitting w/o track #"
		  <<TrkPntr[ntNS-1-iOut]<<endl; CovMatrix(Cnew,vtrks );
	      Print3(vtrks,chi2s[iOut]/(2*(iact-1)-3));
	    }
	    if (ok) break; // Success, iter: don't reset status and CsVTrack's
	    (*iTr).setStatus(0); resetTrks(X0,vtrks);
	    if      (iter==0) iOut = 1;      // Prepare for 2nd iter
	    else if (iter==1) {              // Prepare for 3rd iter...
	      if (chi2s[1]<chi2s[0]) iOut = 1; else iOut = 0;
	      //                          ...if chi2 OK retain best
	      vChi2Out = chi2s[iOut]/(2*iact-5)>Chi2VCut_;
	      if (vChi2Out) break;      //...else give-up
	      else ok = true;
	    }
	  }
	  if (ok) {                     //...if chi2 OK: update
	    Chi2tot = chi2s[iOut]; iact--; Cn = Cnew; Xn = Xnew;
	  }
	}
#endif
	if (!ok) {
	  if (Print_[2])
	    cout<<"CsKF::Kalman==> Track with Chi2="<<Chi2trk[ntNS-1]<<" OUT\n";
	  for (it = 0, iTr = vtrks.begin(); it<TrkPntr[ntNS-1]; it++) iTr++;
	  (*iTr).setStatus(2); iact--;             // ... remove worst...
	  InvFilter(Cn,Xn,*iTr,Chi2trk[ntNS-1],&Xnew,&Cnew); // ...inverse KF
	  Cn = Cnew; Xn = Xnew; Smoother(C0,Cn,X0,Xn,Chi2tot,vtrks );
	}
      }
      else iend = 0;
      
    }
#define Kalman_ReFIT_WOB  // Refit w/o best
#ifdef Kalman_ReFIT_WOB
    if (Chi2tot/(2*iact-3)>4 && // Bad vertex chi2 &&...
	// ...suspicious configuration. We have in view cases where a track from
	// a V0, or not related to the primary interaction, hijacks the pVertex
	// and prevents truly secondaries from being associated.
	//  The hijacking can only take place if few true secondaries are
	// present, otherwise these would tip the scales in their favour.
	// Therefore considering: 4 or 6 tracks in all in pVertex, and 2 or 3
	// tracks eventually retained. E.g. beam + mu' + 2 V0-decays, one of
	// which is the hijacker + 1 or 2 true 2ndaries. Or: beam + diffracted
	// pion + 2*2 V0-decays, one of which is the hijacker.
	2<=iact && iact<=3 &&
	iact+1<=(int)vtrks.size() && (int)vtrks.size()<=6) {
      // Alternative working variables:
      TMtx Ca(C0), Xa(X0); double Chi2alt; int jact, Njter;
      const int NPATS = NTRK_MAX/32+1; unsigned int statusPats[NPATS]; int it;
      for (iTr = vtrks.begin(), it = 0, jact = vtrks.size(),
	     memset(statusPats,0,NPATS*sizeof(unsigned int)); iTr!=vtrks.end();
	   iTr++, it++) {
	if      ((*iTr).getStatus()==2) {// Backup result of original fit
	  statusPats[it/32] |= 1<<it%32;
	  (*iTr).setStatus(0);
	}
	else if ((*iTr).getStatus()==0) {// Exclude hijacker (not special and formerly associated)
	  (*iTr).setStatus(2); jact--;
	  if (Print_[0] || Print_[1]) {
	    printf("\n**************************************************\n");
	    printf("CsKF::Kalman==> Refitting w/o track @ %.3f GeV\n",
		   1/(*iTr).getCop());
	  }
	}
      }
      if(Print_[1]) Print2( vtrks );
      Globalfit(Ca,Xa,vtrks); Smoother(C0,Ca,X0,Xa,Chi2alt,vtrks);
      if (Print_[1]) { CovMatrix(Ca,vtrks); Print3(vtrks,Chi2alt/(2*jact-3)); }
      if (Print_[0] || Print_[1])
	cout<<" (again)       GLOBAL FIT:  #tracks ="<<vtrks.size()
	    <<" Chi2tot="<<Chi2alt<<", x,y,z="<<Xa(1)<<","<<Xa(2)<<","<<Xa(3)
	    <<endl<<endl;
      for (Njter = 0, iend = 1; jact>1 && iend!=0 ;
	   Njter++) {
	int it, ntNS, TrkPntr[NTRK_MAX];
	for (int i = 0; i<(int)vtrks.size(); i++) TrkPntr[i] = 0;
	for (iTr = vtrks.begin(), it=ntNS = 0; iTr!=vtrks.end(); iTr++, it++) {
	  // ***** LOOP on CsVTrack's ***** 
	  if (!(*iTr).primary()) {         // ***** TRACK ALREADY EXCLUDED: SKIP
	    if (Print_[2]) printf("CsKF::Kalman==> Trk %2d\n",it); continue;
	  }
	  else if (Print_[2]) printf("CsKF::Kalman==> Trk %2d",it);
	  InvFilter(Ca,Xa,*iTr,Chi2trk[ntNS],&Xnew,&Cnew);   // ***** INVERSE KF
	  if (Chi2trk[ntNS]<0) Chi2trk[ntNS] = 1000;
	  if (Print_[2]) Print4(*iTr,Chi2trk[ntNS],Xnew," IKF");
	  if ((*iTr).getStatus()!=5)         // ***** IF Non-Special CsVTrack...
	    TrkPntr[ntNS++] = it;                 //     ...STORE REFERENCE
	}
	if (ntNS==0) break;

	thesort(Chi2trk,TrkPntr,ntNS,1);      // ***** SORT ACCORDING to Chi2trk
	int vNDF = 2*jact-3; float vChi2 = Chi2alt/vNDF;
	if (Print_[2]) Print6(vChi2,Chi2trk,ntNS+1);

	bool vChi2Out = vChi2>Chi2VCut_;
	if ((Chi2trk[ntNS-1]>InvChi2Trk_) ||  // ***** IF WORST'S Chi2trk >>...
	    vChi2Out) {                       //  ...OR TOTAL CHI2 >>...
	  if ((Print_[2]) && Chi2trk[ntNS-1]<=InvChi2Trk_)
	    cout<<"CsKF::Kalman==> vChi2 "<<vChi2<<">"<<Chi2VCut_<<endl;
	  if (Print_[2])
	    cout<<"CsKF::Kalman==> Track with Chi2="<<Chi2trk[ntNS-1]<<" OUT\n";
	  for (it = 0, iTr = vtrks.begin(); it<TrkPntr[ntNS-1]; it++) iTr++;
	  (*iTr).setStatus(2); jact--;             // ... remove worst...
	  InvFilter(Ca,Xa,*iTr,Chi2trk[ntNS-1],&Xnew,&Cnew); // ...inverse KF
	  Ca = Cnew; Xa = Xnew; Smoother(C0,Ca,X0,Xa,Chi2alt,vtrks );
	}
	else iend = 0;
      }
      if (jact>iact && Chi2alt/(2*jact-3)<Chi2tot/(2*iact-3) ||
	  jact==iact && Chi2alt/(2*jact-3)<2) {
	if (Print_[0] || Print_[1]) {
	  printf("CsKF::Kalman==> Retaining the outcome of the refit\n");
	  printf("\n**************************************************\n");
	}
	Cn = Ca; Xn = Xa; Chi2tot = Chi2alt; iact = jact;
      }
      else {
	if (Print_[0] || Print_[1]) {
	  printf("CsKF::Kalman==> Restoring the outcome of the original fit\n");
	  printf("**************************************************\n\n");
	}
	for (iTr = vtrks.begin(), it = 0; iTr!=vtrks.end(); iTr++, it++) {
	  if ((*iTr).getStatus()==5) continue;
	  if (statusPats[it/32]&1<<it%32) (*iTr).setStatus(2);
	  else                            (*iTr).setStatus(0);
	}
      }
    }
#endif
    CovMatrix( Cn, vtrks );
    
  }
  else if (Schema_==1) {        //******** DIRECT KALMAN FILTER *********
    
    iact = 0;
    for(iTr=vtrks.begin(); iTr != vtrks.end(); iTr++ )
      if( (*iTr).getStatus() != 5 ) (*iTr).setStatus( 2 );
      else iact++;
    
    //****** MU MU' FIT ************ 
    Cn = C0; Xn = X0;
    
    Chi2tot = 0;
    for(iTr=vtrks.begin(); iTr != vtrks.end(); iTr++ ) {
      double Chi2 = 0;
      
      //------ DIRECT KALMAN FILTER ---------- 

      DirFilter( Cn, Xn, (*iTr), Chi2 );
      if (Print_[2]) Print4(*iTr,Chi2,Xn,"CsKF::Kalman==> DKF");

      if( (*iTr).getStatus() == 5 ) {

        X0 = Xn; C0 = Cn;
        Chi2tot += Chi2;
        if(Print_[2]) cout<<"CsKF::Kalman==> SPECIAL TRACK."<<endl;

      } else if( Chi2 < DirChi2Trk_ ) {

        (*iTr).setStatus( 0 );
        X0 = Xn; C0 = Cn;
        Chi2tot += Chi2;
        iact++;
        if(Print_[2]) cout<<"CsKF::Kalman==> OK."<<endl;
	
      } else {

        Xn = X0; Cn = C0;
        if(Print_[2]) cout<<"CsKF::Kalman==> BIG CHI2, EXCLUDED."<<endl;
      }
      
    }

    if(iact < 2) { chi2=1e9; return; }
    Smoother( C0,Cn, X0,Xn, Chi2tot, vtrks );
    Niter++;
  }
  
  if(Print_[0]) { 
    if( vtrks.front().getStatus()==5 ) Print7( Xn, Cn, iact, Niter, Chi2tot );
    else                               Print8( Xn, Cn, iact, Niter, Chi2tot );
  }
  
  // ***** ERASE OUT DISCARDED TRACKS *****
  for (iTr = vtrks.begin(); iTr != vtrks.end(); iTr++ )
    if( (*iTr).getStatus() == 2 ){
      vtrks.erase( iTr );
      iTr = vtrks.begin();
    }

  chi2 = Chi2tot; ///(2*iact-3);

  return;
}


bool CsKalmanFitting::Globalfit(TMtx &C0, TMtx &X0, list<CsVTrack> &vVertTrk )
{ 
  TMtx Xn(3), Cn(3,3);
  TMtx C0inv(3,3);
  TMtx SUM1(3), SUM2(3,3), Cninv(3,3), Xd(3), Rnk(5), Chi2(1,1);
  int ierr3 = 0;
  
  for (int i=1; i<4; i++){
    SUM1(i)=0;
    for (int j=1; j<4; j++)
      SUM2(i,j) = 0;
  }

  if(Print_[1]){
    cout<<endl;
    cout<<"CsKF::Globalfit==> ########## INSIDE OF GLOBAL FIT ##########"<<endl;
    X0.Print("X0");
    C0.Print("C0");
  }
  
  C0inv = inv3x3( &C0,ierr3);
  if(ierr3 != 0) cout<<"CsKF::Globalfit==> Error inverting Cninv ierr="<<ierr3<<endl;
  
  list<CsVTrack>::iterator iTrack = vVertTrk.begin();
  for(; iTrack != vVertTrk.end(); iTrack++){
    CsVTrack &it = (*iTrack);
    if(it.primary()){
      SUM1 += it.AkT *(it.GkB *(it.Pk - it.Cke));
      SUM2 += it.AkT *(it.GkB * it.Ak);
    }
  }
  
  Cninv =  C0inv + SUM2;
  Cn = inv3x3(&Cninv,ierr3);
  if(ierr3 != 0) cout<<"CsKF::Globalfit==> Error inverting Cninv ierr="<<ierr3<<endl;
 
  Xn = Cn *(C0inv * X0 + SUM1);
  
  X0 = Xn;
  C0 = Cn;
  
  if( ierr3 !=0 ) return false;
  return true;
}


bool CsKalmanFitting::Smoother(TMtx &C0, TMtx &Cn, TMtx &X0, TMtx &Xn, double &Chi2tot, list<CsVTrack> &vVertTrk )
{
  int ierr3 = 0;
  TMtx Chi2(1,1), Rnk(5), Cninv(3,3);
  TMtx C0inv(3,3), Xd(3);
  
  Cninv = inv3x3( &Cn, ierr3 );
  if(ierr3 != 0) cout<<"CsKF::Smoother==> Error inverting Cninv ierr="<<ierr3<<endl;

  Chi2(1,1) = 0;
  list<CsVTrack>::iterator iTrack = vVertTrk.begin();
  for(; iTrack != vVertTrk.end(); iTrack++ ){
    CsVTrack &it = (*iTrack);
    if(it.primary()){
      it.Qnk = it.Wk *(it.BkT *( it.Gk *(it.Pk - it.Cke - (it.Ak * Xn))));
      it.Pnk = it.Cke + (it.Ak * Xn) + (it.Bk * it.Qnk);
      Rnk = it.Pk - it.Pnk;
      Chi2 += Rnk.t() * (it.Gk * Rnk);
      //it.Enk = Cn * (it.AkT * (it.Gk * (it.Bk * it.Wk)));
      //it.Enk *= -1;
      //it.Dnk = it.Wk + it.Enk.t() * (Cninv * it.Enk);
    }
  }
  Chi2tot = Chi2(1,1);
  
  
  C0inv = inv3x3( &C0,ierr3);
  if(ierr3 != 0) cout<<"CsKF::Globalfit==> Error inverting Cninv ierr="<<ierr3<<endl;
  
  Xd = X0 - Xn;
  Chi2 = Xd.t() * (C0inv * Xd);
  Chi2tot += Chi2(1,1);
  
  
  if( ierr3!=0 ) return false;
  return true;
}


bool CsKalmanFitting::CovMatrix( TMtx &Cn, list<CsVTrack> &vVertTrk )
{
  int ierr3 = 0;
  TMtx Rnk(5), Cninv(3,3);
  
  Cninv = inv3x3( &Cn, ierr3 );
  if(ierr3 != 0) cout<<"CsKF::Smoother==> Error inverting Cninv ierr="<<ierr3<<endl;
  
  list<CsVTrack>::iterator iTrack = vVertTrk.begin();
  for(; iTrack != vVertTrk.end(); iTrack++ ){
    CsVTrack &it = (*iTrack);
    if(it.primary()){
      Rnk = it.Pk - it.Pnk;
      it.Enk = Cn * (it.AkT * (it.Gk * (it.Bk * it.Wk)));
      it.Enk *= -1;
      it.Dnk = it.Wk + it.Enk.t() * (Cninv * it.Enk);
    }
  }
  
  
  
  if( ierr3!=0 ) return false;
  return true;
}



void CsKalmanFitting::setTrks(TMtx &XVtx, list<CsTrack*> &trks,     // Inputs...
			      bool isPrimary,map<CsTrack*,bool> &specials,// ...
			      list<CsVTrack> &vtrks)                // Output
{
  // ********** BUILD THE CsVTrack's FROM ARGUMENTS "XVtx" AND "trks" **********

  // Output CsVTrack's are flagged as special following arg "specials" 

  THlx hlx; CsVTrack vtrk;
  TMtx Qke(3), Hke(5),   Winv(3,3),  COV(5,5), COV_REF(5,5), COV_MS(5,5);
  int ierr3 = 0, ierr5 = 0; double xDum = 0; TMtx MDum(5);

  double Z_Target = CsGeom::Instance()->getTargetCenter(); // in mm

  // ***** DEFINE "CurrentRefPlane_/Beam_", VALID for THIS VERY VERTEX *****
  //   This determination encompasses several different versions, depending
  // upon options "RefPlane", "RefBeam", "RefMargin". (Some combinations of
  // options are meant for backward, or otherwise, compatibility):
  //   I) Original, Alex' version:
  //      - RefBeam > 1e6 => No reference for beam, which helix is furthermore
  //       considered @ 1st measured point (probably because it was overlooked
  //       that for beam back helix is more appropriate than front one).
  //      - RefMargin = 0, i.e. vertex can come as close as possible to
  //       ref. plane.
  //      - RefPlane specified (e.g., 500mm)
  //  II) No ref. plane version:
  //      - RefBeam, RefPlane > 1e6: No ref. plane (beam's helix @ LAST point).
  // III) Fixed beam ref. planes: RefBeam != 0.
  //   In addition:
  //      - If RefMargin != 0, (spectro) ref. plane is shifted from its
  //       base value so that vertex<RefPlane-RefMargin, for outside,
  //       downstream of, target vertices.
  //   The, base, reference is expected to lie outside target's field, so that
  // derivative w.r.t. momentum take full account of this field. 
  //   There is an alternative for what concerns beam reference:
  //      - One may argue that it has to be @ a same distance from vertex as
  //       that of spectro tracks. In addition it has not to fulfill the
  //       outside target's field condition mentioned supra, for its momentum
  //       is not correlated with the rest of its helix parameters. This is
  //       provided for by option (IV).
  //      - Alternatively (Empirically it TENDS to be the best solution) a
  //       fixed beam ref. gives a lever arm to tune the beam position @
  //       the vertex when the latter is far from its last measured point.
  //  IV) Variable beam ref.: RefBeam_ == 0. Beam is referenced @ a
  //     plane symmetrical of the ref. plane defined for spectro tracks w.r.t.
  //     vertex.
  if (RefPlane_<1e6) {
    if (XVtx(1)<RefPlane_-RefMargin_) CurrentRefPlane_ = RefPlane_;
    else if (RefMargin_)              CurrentRefPlane_ = XVtx(1)+RefMargin_;
    else                              CurrentRefPlane_ = 1.1e6;
    if (RefBeam_<1e6) {
      if (RefBeam_==0) CurrentRefBeam_ = 2*XVtx(1)-CurrentRefPlane_;  // (IV)
      else             CurrentRefBeam_ = RefBeam_;                    // (III)
      const vector<CsHelix> &tpar = trks.front()->getHelices(); 
      if (CurrentRefBeam_<tpar[0].getZ()+10) CurrentRefBeam_ = 1.1e6;
      BeamHelixFirst_ = false;
    }
    else {
      CurrentRefBeam_ = 1.1e6; BeamHelixFirst_ = true;                // (I)
    }
  }
  else {                                                              // (II)
    CurrentRefPlane_ = 1.1e6; CurrentRefBeam_ = 1.1e6; BeamHelixFirst_ = false;
  }

  list<CsTrack*>::iterator it; bool hasBeam;
  for (it = trks.begin(), hasBeam = false; it!=trks.end(); it++) {

    // ********** LOOP ON CsTracks's IN ARG "trks" **********

    const vector<CsHelix> &tpar = (*it)->getHelices(); 
    bool isBeam = tpar.front().getZ()<Z_Target;

    // ***** EXTRAPOLATE to some REF. PLANE *****
    bool doExtrap(false); double zRef(CurrentRefBeam_);
    if (isBeam) {                         // BEAM... 
      if (CurrentRefBeam_<1e6) {
	doExtrap = true; zRef = CurrentRefBeam_;
      }
      else {
	if (!BeamHelixFirst_)
	  hlx.ImportHelix(tpar.back());           // ...LAST measured plane 
	else // (FOR BACKWARD COMPATIBILITY and although not the best choice...)
	  hlx.ImportHelix(tpar.front());          // (I) ...1ST measured plane
      }
    }
    else {                                // SPECTRO...
      if (CurrentRefPlane_<1e6 &&
	  CurrentRefPlane_<tpar[0].getZ()) {   // ...Do NOT extrapolate forwward
	zRef = CurrentRefPlane_; doExtrap = true; 
      }
      else hlx.ImportHelix(tpar[0]);           // ...1ST  measured plane
    }
    if (doExtrap) {
      THlx hlxD;
      if (isBeam) hlxD.ImportHelix(tpar.back()); // BEAM   : LAST measured plane
      else	  hlxD.ImportHelix(tpar[0]);     // SPECTRO: 1ST  measured plane
      hlx(0) = zRef; if (!hlxD.Extrapolate(hlx))
	CsErrLog::msg(elError,__FILE__,__LINE__,
		      "Error extrapolating from %.3f -> Ref. Plane @ %.3f cm",
		      hlxD(0),zRef);
    }

    // Covar matrix @ ref. plane, as determined by track reco
    hlx.Get(xDum,vtrk.Pk,COV_REF);
    GetDerivs( &hlx, XVtx, vtrk.Ak, vtrk.Bk, Qke, Hke, COV); // COV is dummy
    bool muWoBMS = false; if (isBeam) {
      // The off-diagonal components involving the 5th coordinate should be and
      // should have been kept equal to zero throughout the extrapolation to the
      // reference plane, because there is no (computed) correlation between the
      // momentum, determined in the BMS (or assigned from option), and the
      // other DoFs, determined in the scifi/Si telescope. But the imprecision
      // in the extrapolation breaks the rule. We restore it here, in particular
      // in order to keep the momentum fixed throughout the vertexing when it
      // has been assigned a (a priori round number) value from option. 
      for (int k = 1; k<=4; k++) { COV_REF(k,5)=COV_REF(5,k) = 0; }
      // Determine if beam has BMS. If not, and if the uncertainty on its
      // momentum is large, it means that it's a track to which a momentum
      // has been assigned via option. Changing the value of this momentum,
      // either via energy loss or by imprecise extrapolation would have the
      // adverse effect of modifying the a priori round number that has been
      // assigned and hence of blurring this remarkable signature.
      const list<CsZone*> &zones = (*it)->getZones();
      list<CsZone*>::const_iterator iz; int hasBMS;
      for (iz = zones.begin(), hasBMS = 0; iz!=zones.end(); iz++)
	if ((*iz)->getName()=="BMSzone") { hasBMS = 1; break; }
      if (!hasBMS) {
	const CsHelix &trackingHelix = (*it)->getHelices().front();
	double r = trackingHelix.getCop(), sr2 =  trackingHelix.getCov()[14];
	if (fabs(sqrt(sr2)/r)>.02) {// E.g. P=160+/-15/sqrt(12) GeV yields ~2.7%
	  muWoBMS = true;
	  vtrk.Pk(5)=Qke(3)=Hke(5) = hlx(5);
	  // Momentum (5th coord) should be invariant except for energy loss and
	  // energy loss is meaningless in the present case.
	  for (int k = 1; k<=3; k++) vtrk.Ak(5,k) = 0;
	  for (int k = 1; k<=2; k++) vtrk.Bk(5,k) = 0; vtrk.Bk(5,3) = 1;
	  COV_REF(5,5) = 1e-15;
	}
      }
    }

    if (Agm_) { // ***** AUGMENTING Track's error TO MAKE FOR MS *****
      THlx HelRef(hlx); HelRef *= 0.00001*0.00001;
      THlx HelVtx; HelVtx(0) = XVtx(1); if (!HelRef.Extrapolate(HelVtx))
	CsErrLog::msg(elError,__FILE__,__LINE__,
		      "Error Extrapolating: Ref. plane @ %.1f -> Vertex @ %.1f",
		      HelRef(0),HelVtx(0));
      HelVtx.Get(xDum,MDum,COV_MS);
      for(int k = 1; k<=5; k++) { COV_MS(k,5)=COV_MS(5,k) = 0; }
      COV_MS(5,5) = 0;
    } else {
      //COV_MS *= 0;
      for(int k=1;k<6;k++) for(int l=1;l<6;l++) COV_MS(k,l)=0;
    }
    COV = COV_REF + COV_MS;

    // ***** COMPUTE ALL OF CsVTrack MATRICES *****
    vtrk.Gk  = COV.i5(ierr5);
    if (ierr5!=0) {
      CsErrLog::msg(elError,__FILE__,__LINE__,
		    "V inversion error %d: Wrong Track's Cov Matrix",ierr5);
      COV.Print("CsKalmanFitting::setTrks:"); cout<<endl;
    }

    vtrk.Cke = Hke - vtrk.Ak*XVtx - vtrk.Bk*Qke;
    if (muWoBMS) vtrk.Cke(5) = 0;
    vtrk.AkT = vtrk.Ak.t(); vtrk.BkT = vtrk.Bk.t();
    Winv   = vtrk.BkT * (vtrk.Gk * vtrk.Bk);
    vtrk.Wk  = inv3x3(&Winv,ierr3);
    if (ierr3!=0) CsErrLog::msg(elError,__FILE__,__LINE__,
				"Error inverting Winv: %d",ierr3);
    vtrk.GkB = vtrk.Gk - vtrk.Gk * (vtrk.Bk * (vtrk.Wk * ( vtrk.BkT * vtrk.Gk )));

    vtrk.setChi2(-1); vtrk.setTrkRef(*it);

    if (isPrimary &&

	// ***** PRIMARY: FLAG SPECIAL CsVTrack's ACCORDING TO "specials" *****

	specials.find(*it)!=specials.end()) {
      if (isBeam) {                                     // => BEAM MUON
	vtrk.setStatus(5); hasBeam = true;
      }
      else if (Specials_ && it==(++trks.begin())) {     // => SCATTERED MUON
	if (hasBeam) {
	  vtrk.setStatus(5);
	  // Set Schema_=IKF (Remnant of a time where DKF was still in use,
	  //                  Retained here to be on the safe side)
	  Schema_ = 0;
	}
	else CsErrLog::mes(elFatal,"Beam track is not the first!");
      }
      else vtrk.setStatus(0);                           // => ORDINARY CsVtrack
    }
    else vtrk.setStatus(0);                   //*** ORDINARY CsTrack => CsVTrack

    vtrks.push_back( vtrk );
  }

 
  return;
}

void CsKalmanFitting::resetTrks(TMtx &XVtx,             // Input
				list<CsVTrack> &vtrks)  // In/Output
{
  // ********** RESET THE CsVTrack's GIVEN ARGUMENT VERTEX XVtx **********

  THlx hlx;
  TMtx Qke(3), Hke(5),   Winv(3,3),  COV(5,5), COV_REF(5,5), COV_MS(5,5);
  int ierr3 = 0, ierr5 = 0; double xDum = 0; TMtx MDum(5);

  double Z_Target = CsGeom::Instance()->getTargetCenter(); // in mm

  // ***** RE-DEFINE "CurrentRefPlane_/Beam_ (Cf. "setTrks" supra) *****
  if (RefPlane_<1e6) {
    if (XVtx(1)<RefPlane_-RefMargin_) CurrentRefPlane_ = RefPlane_;
    else if (RefMargin_)              CurrentRefPlane_ = XVtx(1)+RefMargin_;
    else                              CurrentRefPlane_ = 1.1e6;
    if (RefBeam_<1e6) {
      if (RefBeam_==0) CurrentRefBeam_ = 2*XVtx(1)-CurrentRefPlane_;  // (IV)
      else             CurrentRefBeam_ = RefBeam_;                    // (III)
      const vector<CsHelix> &tpar =
	vtrks.front().getAssociatedTrk()->getHelices();
      if (CurrentRefBeam_<tpar[0].getZ()+10) CurrentRefBeam_ = 1.1e6;
    }
  }

  list<CsVTrack>::iterator ivt;
  for (ivt = vtrks.begin(); ivt!=vtrks.end(); ivt++) {
    if ((*ivt).getStatus()==2) continue;

    // ********** LOOP ON, ACTIVE, CsVTracks's IN ARG "vtrks" **********

    CsVTrack &vtrk = *ivt;
    const vector<CsHelix> &tpar = vtrk.getAssociatedTrk()->getHelices();
    bool isBeam = tpar.front().getZ()<Z_Target;

    // ***** EXTRAPOLATE to some REF. PLANE *****
    bool doExtrap(false); double zRef(CurrentRefBeam_);
    if (isBeam) {                         // BEAM... 
      if (CurrentRefBeam_<1e6) {
	doExtrap = true; zRef = CurrentRefBeam_;
      }
      else {
	if (!BeamHelixFirst_)
	  hlx.ImportHelix(tpar.back());           // ...LAST measured plane 
	else // (FOR BACKWARD COMPATIBILITY and although not the best choice...)
	  hlx.ImportHelix(tpar.front());          // (I) ...1ST measured plane
      }
    }
    else {                                // SPECTRO...
      if (CurrentRefPlane_<1e6 &&
	  CurrentRefPlane_<tpar[0].getZ()) {   // ...Do NOT extrapolate forwward
	zRef = CurrentRefPlane_; doExtrap = true; 
      }
      else hlx.ImportHelix(tpar[0]);           // ...1ST  measured plane
    }
    if (doExtrap) {
      THlx hlxD;
      if (isBeam) hlxD.ImportHelix(tpar.back()); // BEAM   : LAST measured plane
      else	  hlxD.ImportHelix(tpar[0]);     // SPECTRO: 1ST  measured plane
      hlx(0) = zRef; if (!hlxD.Extrapolate(hlx))
	CsErrLog::msg(elError,__FILE__,__LINE__,
		      "Error extrapolating from %.3f -> Ref. Plane @ %.3f cm",
		      hlxD(0),zRef);
    }

    // Covar matrix @ ref. plane, as determined by track reco
    hlx.Get(xDum,vtrk.Pk,COV_REF);
    GetDerivs( &hlx, XVtx, vtrk.Ak, vtrk.Bk, Qke, Hke, COV); // COV is dummy
    bool muWoBMS = false; if (isBeam) {
      // The off-diagonal components involving the 5th coordinate should be and
      // should have been kept equal to zero throughout the extrapolation to the
      // reference plane, because there is no (computed) correlation between the
      // momentum, determined in the BMS (or assigned from option), and the
      // other DoFs, determined in the scifi/Si telescope. But the imprecision
      // in the extrapolation breaks the rule. We restore it here, in particular
      // in order to keep the momentum fixed throughout the vertexing when it
      // has been assigned a (a priori round number) value from option. 
      for (int k = 1; k<=4; k++) { COV_REF(k,5)=COV_REF(5,k) = 0; }
      // Determine if beam has BMS. If not, and if the uncertainty on its
      // momentum is large, it means that it's a muon track to which a momentum
      // has been assigned via option. Changing the value of this momentum,
      // either via energy loss or by imprecise extrapolation would have the
      // adverse effect of modifying the a priori round number that has been
      // assigned and hence of blurring this remarkable signature.
      const CsTrack *cst = (*ivt).getAssociatedTrk();
      const list<CsZone*> &zones = cst->getZones();
      list<CsZone*>::const_iterator iz; int hasBMS;
      for (iz = zones.begin(), hasBMS = 0; iz!=zones.end(); iz++)
	if ((*iz)->getName()=="BMSzone") { hasBMS = 1; break; }
      if (!hasBMS) {
	const CsHelix &trackingHelix = cst->getHelices().front();
	double r = trackingHelix.getCop(), sr2 =  trackingHelix.getCov()[14];
	if (fabs(sqrt(sr2)/r)>.02) {// E.g. P=160+/-15/sqrt(12) GeV yields ~2.7%
	  muWoBMS = true;
	  vtrk.Pk(5)=Qke(3)=Hke(5) = hlx(5);
	  // Momentum (5th coord) should be invariant except for energy loss and
	  // energy loss is meaningless in the present case.
	  for (int k = 1; k<=3; k++) vtrk.Ak(5,k) = 0;
	  for (int k = 1; k<=2; k++) vtrk.Bk(5,k) = 0; vtrk.Bk(5,3) = 1;
	  COV_REF(5,5) = 1e-15;
	}
      }
    }

    if (Agm_) { // ***** AUGMENTING Track's error TO MAKE FOR MS *****
      THlx HelRef(hlx); HelRef *= 0.00001*0.00001;
      THlx HelVtx; HelVtx(0) = XVtx(1); if (!HelRef.Extrapolate(HelVtx))
	CsErrLog::msg(elError,__FILE__,__LINE__,
		      "Error Extrapolating: Ref. plane @ %.1f -> Vertex @ %.1f",
		      HelRef(0),HelVtx(0));
      HelVtx.Get(xDum,MDum,COV_MS);
      for(int k = 1; k<=5; k++) { COV_MS(k,5)=COV_MS(5,k) = 0; }
      COV_MS(5,5) = 0;
    } else {
      //COV_MS *= 0;
      for(int k=1;k<6;k++) for(int l=1;l<6;l++) COV_MS(k,l)=0;
    }
    COV = COV_REF + COV_MS;

    // ***** COMPUTE ALL OF CsVTrack MATRICES *****
    vtrk.Gk  = COV.i5(ierr5);
    if (ierr5!=0) {
      CsErrLog::msg(elError,__FILE__,__LINE__,
		    "V inversion error %d: Wrong Track's Cov Matrix",ierr5);
      COV.Print("CsKalmanFitting::setTrks:"); cout<<endl;
    }

    vtrk.Cke = Hke - vtrk.Ak*XVtx - vtrk.Bk*Qke;
    if (muWoBMS) vtrk.Cke(5) = 0;
    vtrk.AkT = vtrk.Ak.t(); vtrk.BkT = vtrk.Bk.t();
    Winv   = vtrk.BkT * (vtrk.Gk * vtrk.Bk);
    vtrk.Wk  = inv3x3(&Winv,ierr3);
    if (ierr3!=0) CsErrLog::msg(elError,__FILE__,__LINE__,
				"Error inverting Winv: %d",ierr3);
    vtrk.GkB = vtrk.Gk - vtrk.Gk * (vtrk.Bk * (vtrk.Wk * ( vtrk.BkT * vtrk.Gk )));
    vtrk.setChi2( -1 );
  }
  return;
}



bool CsKalmanFitting::InvFilter(TMtx &Cn, TMtx &Xn, CsVTrack &it, double &Chi2trk, TMtx *Xnk, TMtx *Cnk){
 
  TMtx Cnkinv(3,3), Rnk(5), Xnd(3), Chi2kS(1,1);
  TMtx Cglobinv(3,3);
  int ierr3=0;
  
  Rnk = it.Pk - it.Pnk;
  
  Cglobinv = inv3x3(&Cn,ierr3);
  if(ierr3 != 0) cout<<"CsKF::InvFilter==> Error inverting Cn ierr="<<ierr3<<endl;
  
  Cnkinv = Cglobinv - it.AkT *(it.GkB * it.Ak);
  (*Cnk) = inv3x3(&Cnkinv,ierr3);
  if(ierr3 != 0) cout<<"CsKF::InvFilter==> Error inverting Cnkinv ierr="<<ierr3<<endl;
  
  (*Xnk) = (*Cnk) *(Cglobinv * Xn - it.AkT *(it.GkB *(it.Pk - it.Cke)));
  Xnd = Xn - (*Xnk);
  Chi2kS = Xnd.t() *(Cnkinv * Xnd) + Rnk.t()* (it.Gk * Rnk);
  Chi2trk = Chi2kS(1,1);
  it.setChi2(Chi2trk);
  
  if(ierr3 !=0 ) return false;
  return true;
}



bool CsKalmanFitting::DirFilter(TMtx &C0, TMtx &X0, CsVTrack &it, double &Chi2trk ){

  TMtx Ck(3,3), CkInv(3,3), Rk(5), Xkd(3), Xk(3), Chi2k(1);
  TMtx C0Inv(3,3);
  int ierr3=0;

  C0Inv = inv3x3(&C0,ierr3);
  if(ierr3 != 0) cout<<"CsKF::DirFilter==> Error inverting Cn ierr="<<ierr3<<endl;

  CkInv = C0Inv + (it.AkT * it.GkB * it.Ak);
  Ck = inv3x3(&CkInv,ierr3);
  if(ierr3 != 0) cout<<"CsKF::DirFilter==> Error inverting Cnkinv ierr="<<ierr3<<endl;

  Xk = Ck *(C0Inv * X0 + it.AkT *(it.GkB *(it.Pk - it.Cke)));
  it.Qnk = it.Wk * it.BkT * it.Gk * ( it.Pk - it.Cke - (it.Ak * Xk) );
  it.Pnk = it.Cke + (it.Ak * Xk) + (it.Bk * it.Qnk);
  Rk = it.Pk - it.Pnk;

  Xkd = Xk - X0;
  Chi2k = Xkd.t() *(C0Inv * Xkd) + Rk.t()* (it.Gk * Rk);

  Chi2trk = Chi2k(1);
  it.setChi2(Chi2trk);

  C0 = Ck;
  X0 = Xk;

  if(ierr3!=0) return false;
  return true;
}




void CsKalmanFitting::GetDerivs(THlx *HTrck,TMtx &Xe,TMtx &Ake,
			       TMtx &Bke,TMtx &Qke, TMtx &Hke, TMtx &COV)
{
  double xDum;
  THlx HelExp,HelExpSh,HelExpSh2,HelRefSh,HelRefSh2;
  TMtx  HkeSh(5), HkeSh2(5), Vdum(5,5);

  THlx &HelRef = *HTrck;

  //const double dX[3] = {1,0.1,0.1};             // steps for dH/dX
  const double dQ[3] = {0.001, 0.001, 0.01};     // steps for dH/dQ

  const double dX[3] = {2.,0.02,0.02};             // steps for dH/dX
  //const double dQ[3] = {0.0005, 0.0005, 0.005};     // steps for dH/dQ

  HelExp(0)=Xe(1);
  if( !HelRef.Extrapolate(HelExp) ) {
    cout<<endl<<"CsKF::GetDerivs=> Wrong track seting."
        <<" You try to extrapolate from "<<HelRef(0)<<" to "
        <<HelExp(0)<<endl<<endl;
  }
  
  // Get covarience matrix COV. HkeSh is used as dummy
  if( Agm_ ) HelExp.Get(xDum,HkeSh,COV);

  // Qk[3] (dX2/dX1, dX3/dX1, q/P) at Xe
 
  Qke(1)=HelExp(3);
  Qke(2)=HelExp(4);
  Qke(3)=HelExp(5);
 
  // Set X2,X3 to X2e, X3e and trace back to the reference plane
 
  HelRefSh(0)=HelRef(0);
  HelExp(1)=Xe(2);
  HelExp(2)=Xe(3);
 
  // HelExp is now defined in Xe with slopes and etc. from Qke
 
  if( !HelExp.Extrapolate(HelRefSh,0) ) {    // trace to the ref. plane
    cout<<endl<<"CsKF::GetDerivs=> Wrong track seting."
        <<" You try to extrapolate from "<<HelExp(0)<<" to "
        <<HelRefSh(0)<<endl<<endl;
  }
  HelRefSh.Get(xDum,Hke,Vdum);
 
  HelExpSh=HelExp;
  HelExpSh(0)+=dX[0];
  HelRefSh(0)=HelRef(0);
  if( !HelExpSh.Extrapolate(HelRefSh,0) ) {
    cout<<endl<<"CsKF::GetDerivs=> Wrong track seting."
        <<" You try to extrapolate from "<<HelExpSh(0)<<" to "
	<<HelRefSh(0)<<endl<<endl;
  }
  HelRefSh.Get(xDum,HkeSh,Vdum);
 
  Ake(1,1)=(HkeSh(1)-Hke(1))/dX[0];
  Ake(2,1)=(HkeSh(2)-Hke(2))/dX[0];
  Ake(3,1)=(HkeSh(3)-Hke(3))/dX[0];
  Ake(4,1)=(HkeSh(4)-Hke(4))/dX[0];
  Ake(5,1)=(HkeSh(5)-Hke(5))/dX[0];
 
  //dHk/dX2 at Xe
 
  HelExpSh=HelExp;
  HelExpSh(1)+=dX[1];
  HelRefSh(0)=HelRef(0);
  if( !HelExpSh.Extrapolate(HelRefSh,0) ) {
    cout<<endl<<"CsKF::GetDerivs=> Wrong track seting."
        <<" You try to extrapolate from "<<HelExpSh(0)<<" to "
        <<HelRefSh(0)<<endl<<endl;
  }
  HelRefSh.Get(xDum,HkeSh,Vdum);
 
  Ake(1,2)=(HkeSh(1)-Hke(1))/dX[1];
  Ake(2,2)=(HkeSh(2)-Hke(2))/dX[1];
  Ake(3,2)=(HkeSh(3)-Hke(3))/dX[1];
  Ake(4,2)=(HkeSh(4)-Hke(4))/dX[1];
  Ake(5,2)=(HkeSh(5)-Hke(5))/dX[1];
 
  //dHk/dX3 at Xe
 
  HelExpSh=HelExp;
  HelExpSh(2)+=dX[2];
  HelRefSh(0)=HelRef(0);
  if( !HelExpSh.Extrapolate(HelRefSh,0) ) {
    cout<<endl<<"CsKF::GetDerivs=> Wrong track seting."
        <<" You try to extrapolate from "<<HelExpSh(0)<<" to "
        <<HelRefSh(0)<<endl<<endl;
  }
  HelRefSh.Get(xDum,HkeSh,Vdum);
 
  Ake(1,3)=(HkeSh(1)-Hke(1))/dX[2];
  Ake(2,3)=(HkeSh(2)-Hke(2))/dX[2];
  Ake(3,3)=(HkeSh(3)-Hke(3))/dX[2];
  Ake(4,3)=(HkeSh(4)-Hke(4))/dX[2];
  Ake(5,3)=(HkeSh(5)-Hke(5))/dX[2];
 
  //dHk/dQ1 at Xe
 
  HelExpSh=HelExp;
  HelExpSh(3)+=dQ[0];
  HelRefSh(0)=HelRef(0);
  if( !HelExpSh.Extrapolate(HelRefSh,0) ) {
    cout<<endl<<"CsKF::GetDerivs=> Wrong track seting."
        <<" You try to extrapolate from "<<HelExpSh(0)<<" to "
        <<HelRefSh(0)<<endl<<endl;
  }
  HelRefSh.Get(xDum,HkeSh,Vdum);
 
  Bke(1,1)=(HkeSh(1)-Hke(1))/dQ[0];
  Bke(2,1)=(HkeSh(2)-Hke(2))/dQ[0];
  Bke(3,1)=(HkeSh(3)-Hke(3))/dQ[0];
  Bke(4,1)=(HkeSh(4)-Hke(4))/dQ[0];
  Bke(5,1)=(HkeSh(5)-Hke(5))/dQ[0];
 
  //dHk/dQ2 at Xe
 
  HelExpSh=HelExp;
  HelExpSh(4)+=dQ[1];
  HelRefSh(0)=HelRef(0);
  if( !HelExpSh.Extrapolate(HelRefSh,0) ) {
    cout<<endl<<"CsKF::GetDerivs=> Wrong track seting."
        <<" You try to extrapolate from "<<HelExpSh(0)<<" to "
        <<HelRefSh(0)<<endl<<endl;
  }
  HelRefSh.Get(xDum,HkeSh,Vdum);
 
  Bke(1,2)=(HkeSh(1)-Hke(1))/dQ[1];
  Bke(2,2)=(HkeSh(2)-Hke(2))/dQ[1];
  Bke(3,2)=(HkeSh(3)-Hke(3))/dQ[1];
  Bke(4,2)=(HkeSh(4)-Hke(4))/dQ[1];
  Bke(5,2)=(HkeSh(5)-Hke(5))/dQ[1];
 
  //dHk/dQ3 at Xe
 
  HelExpSh=HelExp;
  HelExpSh(5)+=HelExpSh(5)*dQ[2];
  HelRefSh(0)=HelRef(0);
  if( !HelExpSh.Extrapolate(HelRefSh,0) ) {
    cout<<endl<<"CsKF::GetDerivs=> Wrong track seting."
        <<" You try to extrapolate from "<<HelExpSh(0)<<" to "
        <<HelRefSh(0)<<endl<<endl;
  }
  HelRefSh.Get(xDum,HkeSh,Vdum);
 
  Bke(1,3)=(HkeSh(1)-Hke(1))/(HelExp(5)*dQ[2]);
  Bke(2,3)=(HkeSh(2)-Hke(2))/(HelExp(5)*dQ[2]);
  Bke(3,3)=(HkeSh(3)-Hke(3))/(HelExp(5)*dQ[2]);
  Bke(4,3)=(HkeSh(4)-Hke(4))/(HelExp(5)*dQ[2]);
  Bke(5,3)=(HkeSh(5)-Hke(5))/(HelExp(5)*dQ[2]);
};
 


inline TMtx inv3x3(TMtx *M, int &ierr)
{ 
  TMtx Res(3,3);
  double Deter;
 
  ierr = 0;
 
  Res(1,1)=  (*M)(2,2)*(*M)(3,3) - (*M)(2,3)*(*M)(3,2);
  Res(1,2)= -(*M)(2,1)*(*M)(3,3) + (*M)(2,3)*(*M)(3,1);
  Res(1,3)=  (*M)(2,1)*(*M)(3,2) - (*M)(2,2)*(*M)(3,1);
  Res(2,1)= -(*M)(1,2)*(*M)(3,3) + (*M)(1,3)*(*M)(3,2);
  Res(2,2)=  (*M)(1,1)*(*M)(3,3) - (*M)(1,3)*(*M)(3,1);
  Res(2,3)= -(*M)(1,1)*(*M)(3,2) + (*M)(1,2)*(*M)(3,1);
  Res(3,1)=  (*M)(1,2)*(*M)(2,3) - (*M)(1,3)*(*M)(2,2);
  Res(3,2)= -(*M)(1,1)*(*M)(2,3) + (*M)(1,3)*(*M)(2,1);
  Res(3,3)=  (*M)(1,1)*(*M)(2,2) - (*M)(1,2)*(*M)(2,1);
 
  Deter = (*M)(1,1)*Res(1,1) + (*M)(1,2)*Res(1,2) + (*M)(1,3)*Res(1,3);
 
  if( Deter !=0 ){
    for (int i=1; i<4; i++){
      for (int j=1; j<4; j++){
        Res(i,j) /= Deter;
      }
    }
  }
  else{
    cout<<"CsKF::inv3x3==> Error inversing 3x3 matrix "<<endl;
    ierr = 1;
  }
  return(Res);
}



void CsKalmanFitting::ELossCorrection(TMtx& Xn, list<CsVTrack> &vtrks)
{
  //  At this point, only part of energy loss in target may have been accounted
  // for: that part between the beginning of the target ELoss material
  // map and the reference plane.
  //  Now we are about to fill that possible gap between that ref plane and
  // vertex. But, upon define ("Kalman_ELoss_FROM_1ST"), we do it in fact from
  // scratch, resuming extrapolation from 1st measured point. This in order
  // to be able to compare, in MC, energy loss between 1st measured point
  // and vertex:
  //  - as generated
  //  - vs. as reconstructed,
  // by including the 1st->ref part in the variable to be histogrammed
  // for the reconstructed term. (This is only used when "Kalman_DEBUG_ELoss"
  // is defined.)

  // ********** HISTOGRAM BOOKING **********

  static bool first = true;
  static CsHist1F *hEL,*hEL20,*hEL_MC,*hEL20_MC;
  static CsHist2F *h0,*h1,*h2,*h3,*h4,*h5,*h6,*h7; 

#ifdef Kalman_DEBUG_ELoss

  // ***** EXTRA HISTOs FOR DEBUG SPECIAL HISTOGRAMMING *****

  double zs[9] = {-1000,-700,-350,0,300,900,1500,2100,2700};
  double ps[4] = {60,100,130,160};
  // TH2D's for bEAM, AS A F(Z)
  static TH2D *p1_b, *p2_b, *p3_b, *p4_b, *p5_b, *p6_b;
  // TH2D's for SCATTERED mUONS, AS A F(Z), for p SLICES
  static TH2D *p1_m[4], *p2_m[4], *p3_m[4], *p4_m[4], *p5_m[4], *p6_m[4];
  // TH2D's for ``PIONS'' and sECONDARIES (AS A F(p), for Z slices)
  static TH2D *p1, *p2[5], *p3[5], *p4[5], *p5[5], *p6[5];
  static TH2D *p1_s[9], *p2_s[9], *p3_s[9], *p4_s[9];
  // TH2D of EXTRAPOLATION ERROR for BEAM
  static TH1D *ee_b;
  // CORRESPONDING DISTRIBUTIONS
  static TH1D *d1_b, *d2_b, *d3_b, *d4_b, *d6_b;
  static TH1D *d1_m[4], *d2_m[4], *d3_m[4], *d4_m[4], *d6_m[4];
  static TH1D *d1, *d2[5], *d3[5], *d4[5], *d6[5], *d1_s[9], *d2_s[9], *d3_s[9];
#  if Kalman_DEBUG_ELoss > 1
  int iTrk = 0;
#  endif
#endif
  if( hist_ && first ) {
    CsHistograms::SetCurrentPath("/CsKalmanFitting/ELoss");
    hEL = new CsHist1F( "hEL", "Energy losses", 100,0,0.25);
    hEL20 = new CsHist1F( "hEL20", "Energy losses", 100,0,0.25);
    if( CsEvent::Instance()->isAMonteCarloEvent() ) {
      
      hEL_MC   = new CsHist1F("hEL_MC"  , "Energy losses (MC)", 100,0,0.25);
      hEL20_MC = new CsHist1F("hEL20_MC", "Energy losses (MC)", 100,0,0.25);
      
      h0= new CsHist2F("h0","(P_{MC}^{V}-P_{TR}^{ref})/P_{MC} vs P_{MC}^{V}"    , 20 ,0,20 , 50 ,-0.1 ,0.1  );
      h1= new CsHist2F("h1","(P_{MC}^{V}-P_{V}^{V})/P_{MC} vs P_{MC}^{V}"       , 20 ,0,20 , 50 ,-0.1 ,0.1  );
      h2= new CsHist2F("h2","(P_{MC}^{V}-P_{V}^{ref}) /P_{MC} vs P_{MC}^{V}"    , 20 ,0,20 , 50 ,-0.1 ,0.1  );
      h3= new CsHist2F("h3","(P_{MC}^{V}-P_{MC}^{ref})/P_{MC} vs P_{MC}^{V}"    , 20 ,0,20 , 50 ,-0.1 ,0.1  );
      h4= new CsHist2F("h4","(P_{MC}^{ref}-P_{TR}^{ref})/P_{MC} vs P_{MC}^{ref}", 20 ,0,20 , 50 ,-0.1 ,0.1  );
      h5= new CsHist2F("h5","(P_{MC}^{ref}-P_{V}^{ref})/P_{MC} vs P_{MC}^{ref}" , 20 ,0,20 , 50 ,-0.1 ,0.1  );
      
      h6= new CsHist2F("h6","(P_{MC}^{V}-P_{TR}^{ref})/P_{MC} vs 1/P_{MC}^{V}"    , 30 ,-1.5,1.5 , 50 ,-0.1 ,0.1  );
      h7= new CsHist2F("h7","(P_{MC}^{V}-P_{V}^{V})/P_{MC} vs 1/P_{MC}^{V}"       , 30 ,-1.5,1.5 , 50 ,-0.1 ,0.1  );

#ifdef Kalman_DEBUG_ELoss
      // TH2D's AS A F(Z) for bEAM...
      p1_b = new
	TH2D("p1_b","P_{MC}^{last}-P_{T}^{last} vs. z_{MC} - BEAM",
		 15,-1100,400,80,-2.,1.2);
      p2_b = new
	TH2D("p2_b","P_{MC}^{V}-P_{V}^{ref} vs. z_{MC} - BEAM",
		 15,-1100,400,80,-2.,1.2);
      p3_b = new
	TH2D("p3_b","P_{MC}^{V}-P_{V}^{V} vs. z_{MC} - BEAM",
		 15,-1100,400,80,-2.,1.2);
      p4_b = new
	TH2D("p4_b","P_{MC}^{V}-P_{MC}^{last} vs. z_{MC} - BEAM",
		 15,-1100,400,60,-.25,.05);
      p5_b = new
	TH2D("p5_b","P_{V}^{ref}-P_{T}^{last} vs. z_{MC} - BEAM",
		 15,-1100,400,40,-.15,.05);
      p6_b = new
	TH2D("p6_b","E Loss vs. z_{MC} - BEAM",
		 15,-1100,400,60,-.25,.05);
      // ...CORRESPONDING DISTRIBUTIONS
      d1_b = new
	TH1D("d1_b","(P_{MC}^{1st}-P_{T}^{1st})/P_{MC} - BEAM", 50,-.05,.05);
      d2_b = new
	TH1D("d2_b","(P_{MC}^{V}-P_{V}^{ref})/P_{MC}  - BEAM",  50,-.05,.05);
      d3_b = new
	TH1D("d3_b","(P_{MC}^{V}-P_{V}^{V})/P_{MC}  - BEAM",    50,-.05,.05);
      d4_b = new
	TH1D("d4_b","(P_{MC}^{V}-P_{MC}^{last})/P_{MC}  - BEAM",50,-.025,.025);
      d6_b = new
	TH1D("d6_b","(E Loss)/P_{MC}  - BEAM",                  50,-.025,.025);
      // Extrapolation error
      ee_b = new TH1D("ee_b","Extrapolation Error  -  BEAM",50,0,1);
      // TH2D's AS A F(Z) for mUONS...
      char name[] = "p1_s0";
      char title[] =
	"(P_{MC}^{1st}-P_{T}^{1st})/P_{MC} vs. P_{MC} - #mu',-1000<Z<-700  ";
      char pRange[] = ",130<p<160";
      for (int pBin = 0; pBin<4; pBin++) {
	if (pBin==0) strcpy(pRange,"");
	else sprintf(pRange,",%.0f<p<%.0f",ps[pBin-1],ps[pBin]);
	sprintf(name,"p1_m%d",pBin);
	sprintf(title,"P_{MC}^{1st}-P_{T}^{1st} vs. z_{MC}- #mu'%s", pRange);
	p1_m[pBin] = new TH2D(name,title,15,-1100,400,80,-.6,1.);
	sprintf(name,"p2_m%d",pBin);
	sprintf(title,"P_{MC}^{V}-P_{V}^{ref} vs. z_{MC} - #mu'%s", pRange);
	p2_m[pBin] = new TH2D(name,title,15,-1100,400,80,-.6,1.);
	sprintf(name,"p3_m%d",pBin);
	sprintf(title,"P_{MC}^{V}-P_{V}^{V} vs. z_{MC} - #mu'%s",   pRange);
	p3_m[pBin] = new TH2D(name,title,15,-1100,400,80,-.6,1.);
	sprintf(name,"p4_m%d",pBin);
	sprintf(title,"P_{MC}^{V}-P_{MC}^{1st} vs. z_{MC} - #mu'%s",pRange);
	p4_m[pBin] = new TH2D(name,title,15,-1100,400,60,-.05,.25);
	sprintf(name,"p5_m%d",pBin);
	sprintf(title,"P_{V}^{ref}-P_{T}^{1st} vs. z_{MC} - #mu'%s",pRange);
	p5_m[pBin] = new TH2D(name,title,15,-1100,400,40,-.05,.15);
	sprintf(name,"p6_m%d",pBin);
	sprintf(title,"E Loss vs. z_{MC} - #mu'%s",pRange);
	p6_m[pBin] = new TH2D(name,title,15,-1100,400,60,-.05,.25);
	// ...CORRESPONDING DISTRIBUTIONS
	sprintf(name,"d1_m%d",pBin);
	sprintf(title,"(P_{MC}^{1st}-P_{T}^{1st})/P_{MC}  - #mu'%s",pRange);
	d1_m[pBin] = new TH1D(name,title,50,-.05,.05);
	sprintf(name,"d2_m%d",pBin);
	sprintf(title,"(P_{MC}^{V}-P_{V}^{ref})/P_{MC}  - #mu'%s",  pRange);
	d2_m[pBin] = new TH1D(name,title,50,-.05,.05);
	sprintf(name,"d3_m%d",pBin);
	sprintf(title,"(P_{MC}^{V}-P_{V}^{V})/P_{MC}  - #mu'%s",    pRange);
	d3_m[pBin] = new TH1D(name,title,50,-.05,.05);
	sprintf(name,"d4_m%d",pBin);
	sprintf(title,"(P_{MC}^{V}-P_{MC}^{1st})/P_{MC}  - #mu'%s", pRange);
	d4_m[pBin] = new TH1D(name,title,50,-.025,.025);
	sprintf(name,"d6_m%d",pBin);
	sprintf(title,"(E Loss)/P_{MC}  - #mu'%s",                  pRange);
	d6_m[pBin] = new TH1D(name,title,50,-.025,.025);
      }
      // TH2D and DISTRIBUTION for ``PIONS'' @ 1ST MEASURED POINT
      double pbins[15] = {0,1, 2, 3, 4,5,6,7,8,
			   10,12,14,16,
			  20,24};
      p1   = new
	TH2D("p1"  ,"P_{MC}^{1st}-P_{T}^{1st} vs. P_{MC}",14,pbins,60,-.05,.25);
      d1   = new
	TH1D("d1"  ,"(P_{MC}^{1st}-P_{T}^{1st})/P_{MC}",        50,-.1,.1);
      // TH2D's AS a F(P), for Z SLICES, FOR ``PIONS'' and sECONDARIES,...
      // ...and CORRESPONDING DISTRIBUTIONS
      char zRange[] = " - -1000<Z<-700";
      for (int zBin = 0; zBin<9; zBin++) {
	// SECONDARIES...
	if (zBin==0) strcpy(zRange,"");
	else sprintf(zRange," - %.0f<z<%.0f",zs[zBin-1],zs[zBin]);
	sprintf(name,"p1_s%d",zBin);
	sprintf(title,
		"P_{MC}^{1st}-P_{T}^{1st} vs. P_{MC}%s",zRange);
	p1_s[zBin]   = new TH2D(name,title,14,pbins,60,-.05,.25);
	sprintf(name,"p2_s%d",zBin);
	sprintf(title,
		"P_{MC}^{V}-P_{V}^{ref} vs. P_{MC}%s",zRange);
	p2_s[zBin]   = new TH2D(name,title,14,pbins,60,-.05,.25);
	sprintf(name,"p3_s%d",zBin);
	sprintf(title,
		"P_{MC}^{V}-P_{V}^{V} vs. P_{MC}%s",zRange);
	p3_s[zBin]   = new TH2D(name,title,14,pbins,60,-.05,.25);
	sprintf(name,"p4_s%d",zBin);
	sprintf(title,
		"P_{MC}^{V}-P_{MC}^{1st} vs. P_{MC}%s",zRange);
	p4_s[zBin]   = new TH2D(name,title,14,pbins,60,.0,.15);
	sprintf(name,"d1_s%d",zBin);
	sprintf(title,
		"(P_{MC}^{1st}-P_{T}^{1st})/P_{MC}%s",zRange);
	d1_s[zBin]   = new TH1D(name,title,50,-.1,.1);
	sprintf(name,"d2_s%d",zBin);
	sprintf(title,
		"(P_{MC}^{V}-P_{V}^{ref})/P_{MC}%s",zRange);
	d2_s[zBin]   = new TH1D(name,title,50,-.1,.1);
	sprintf(name,"d3_s%d",zBin);
	sprintf(title,
		"(P_{MC}^{V}-P_{V}^{V})/P_{MC}%s",zRange);
	d3_s[zBin]   = new TH1D(name,title,50,-.1,.1);

	if (zBin>4) continue;                    // ...PRIMARY ``PIONS''
	sprintf(name,"p2%d",zBin);
	sprintf(title,
		"P_{MC}^{V}-P_{V}^{ref} vs. P_{MC}%s",zRange);
	p2[zBin]     = new TH2D(name,title,14,pbins,60,-.05,.25);
	sprintf(name,"p3%d",zBin);
	sprintf(title,
		"P_{MC}^{V}-P_{V}^{V} vs. P_{MC}%s",zRange);
	p3[zBin]     = new TH2D(name,title,14,pbins,60,-.05,.25);
	sprintf(name,"p4%d",zBin);
	sprintf(title,
		"P_{MC}^{V}-P_{MC}^{1st} vs. P_{MC}%s",zRange);
	p4[zBin]     = new TH2D(name,title,14,pbins,60,.0,.15);
	sprintf(name,"p5%d",zBin);
	sprintf(title,
		"P_{V}^{ref}-P_{T}^{1st} vs. P_{MC}%s",zRange);
	p5[zBin]     = new TH2D(name,title,14,pbins,60,-.05,.05);
	sprintf(name,"p6%d",zBin);
	sprintf(title,
		"E Loss vs. P_{MC}%s",zRange);
	p6[zBin]     = new TH2D(name,title,14,pbins,60,.0,.15);
	sprintf(name,"d2%d",zBin);
	sprintf(title,
		"(P_{MC}^{V}-P_{V}^{ref})/P_{MC}%s",zRange);
	d2[zBin]     = new TH1D(name,title,50,-.1,.1);
	sprintf(name,"d3%d",zBin);
	sprintf(title,
		"(P_{MC}^{V}-P_{V}^{V})/P_{MC}%s",zRange);
	d3[zBin]     = new TH1D(name,title,50,-.1,.1);
	sprintf(name,"d4%d",zBin);
	sprintf(title,
		"(P_{MC}^{V}-P_{MC}^{1st})/P_{MC}%s",zRange);
	d4[zBin]     = new TH1D(name,title,50,-.05,.05);
	sprintf(name,"d6%d",zBin);
	sprintf(title,
		"(E Loss)/P_{MC}%s",zRange);
	d6[zBin]     = new TH1D(name,title,50,-.05,.05);
      }
#endif

    }
    CsHistograms::SetCurrentPath("/");
    first = false;
  }

  double Z_Target = CsGeom::Instance()->getTargetCenter(); // in mm
  THlx hlxR, hlxV;
  list<CsVTrack>::iterator ivt; for(ivt = vtrks.begin(); ivt!=vtrks.end();
				    ivt++ ) {

    const CsHelix *trackingHelix =
      &((*ivt).getAssociatedTrk()->getHelices().front());
    double zTrk = trackingHelix->getZ()/10;
    bool beamTrk = zTrk<Z_Target/10; if (beamTrk) {
      // Determine if it has BMS. If not, and if the uncertainty on its momentum
      // is large: skip it. It means that it's a muon track, as opposed to a
      // hadron one, for which, given the large spread in momentum of the muon
      // beam, BMS is required to get a momentum precise enough for the energy
      // loss to be meaningful. Apart from being meaningless, the energy loss
      // would have the adverse effect of modifying the a priori round number
      // that has been assigned to the momentum (via option) and hence of
      // blurring this remarkable signature.
      const CsTrack *cst = (*ivt).getAssociatedTrk(); 
      const list<CsZone*> &zones = cst->getZones();
      list<CsZone*>::const_iterator iz; int hasBMS;
      for (iz = zones.begin(), hasBMS = 0; iz!=zones.end(); iz++)
	if ((*iz)->getName()=="BMSzone") { hasBMS = 1; break; }
      if (!hasBMS) {
	double r = trackingHelix->getCop(), sr2 =  trackingHelix->getCov()[14];
	if (fabs(sqrt(sr2)/r)>.02) {// E.g. P=160+/-15/sqrt(12) GeV yields ~2.7%
	  CsVTrack &vt = *ivt; if (fabs(vt.getCop()-r)>.000002)
	    CsErrLog::msg(elError,__FILE__,__LINE__,"Beam track w/o BMS: Vertex fit doesn't leave P invariant: %.2f -> %.2f GeV",
			  1/r,1/vt.getCop());
	  continue;
	}
      }
      if (!BeamHelixFirst_) {
	trackingHelix = &((*ivt).getAssociatedTrk()->getHelices().back());
	zTrk = trackingHelix->getZ()/10;
      }
    }
    double rTrk = trackingHelix->getCop(), pTrk = fabs(1/rTrk);
    double xVtx, yVtx;
    

    //#define Kalman_ELoss_FROM_1ST
#ifdef Kalman_ELoss_FROM_1ST

    // Resuming Eloss calculation from 1st measured point (as opposed
    // to ref plane) in order to be able to histogram the whole quantity
    // (Does not affect the actual energy loss correction to the momentum).

    double el; double zRef;      // ***** EXTRAPOLATION to some REF. PLANE *****
    if   (beamTrk) zRef = CurrentRefBeam_;
    else           zRef = CurrentRefPlane_;
    if (zRef<zTrk) {   // This holds for beam too: cf. "setTrks"
      hlxR(0) = zTrk;
      hlxR(1) = trackingHelix->getX()/10; hlxR(2) = trackingHelix->getY()/10;
      hlxR(3) = trackingHelix->getDXDZ(); hlxR(4) = trackingHelix->getDYDZ();
      hlxR(5) = rTrk;
      hlxV(0) = zRef;
      if (!hlxR.Extrapolate(hlxV)) {
	CsErrLog::msg(elError,__FILE__,__LINE__,
		      "Error Extrapolating from %.2f to %.2f",hlxR(0),hlxV(0));
      }
      el = hlxV.Mom() - hlxR.Mom();
    }
    else el = 0;
#endif

    // ********** ENERGY LOSS FROM 1ST or REF. PLANE TO VERTEX **********

    if      (beamTrk && CurrentRefBeam_<1e6) hlxR(0) = CurrentRefBeam_;
    else if (!beamTrk && CurrentRefPlane_<1e6 &&
	     CurrentRefPlane_<zTrk)          hlxR(0) = CurrentRefPlane_;
    else                                     hlxR(0) = zTrk;
    hlxR(1) = (*ivt).Pnk(1); hlxR(2) = (*ivt).Pnk(2);
    hlxR(3) = (*ivt).Pnk(3); hlxR(4) = (*ivt).Pnk(4);
    hlxR(5) = (*ivt).Pnk(5);
    hlxV(0) = Xn(1);
    if (!hlxR.Extrapolate(hlxV)) {
      CsErrLog::msg
	(elError,__FILE__,__LINE__,
	 "Wrong track setting: trying to extrapolate from %.2f to %.2f",
	 hlxR(0),hlxV(0));
    }

    double EL = hlxV.Mom() - hlxR.Mom();

    // ***** UPDATE MOMENTUM (q/P in fact) *****

    (*ivt).Qnk(3) = ((*ivt).Qnk(3) < 0 ? -1. : 1. )/( (*ivt).getMom() + EL );

#ifdef Kalman_ELoss_FROM_1ST
    EL += el;    // Total energy loss
#endif
    (*ivt).setELoss( EL );
    
    if( hist_ ){                 // *************** HISTOGRAMING ***************

      hEL->Fill( (*ivt).getELoss() );
      if( (*ivt).getMom() < 20 ) hEL20->Fill( (*ivt).getELoss() );
      
      // **********************************************************************
      //                  ********** MC **********
      // **********************************************************************

      if (CsEvent::Instance()->isAMonteCarloEvent()){

	CsTrack *trk = (*ivt).getAssociatedTrk(); 

#ifdef Kalman_DEBUG_ELoss
	// VERTEX'S (as opposed to track's) ATTRIBUTES...
	static int vertexType;   // =0: Unknowm, =+1: Primary, =-1: Secondary
	static int zBin; static double zMC;
	if (ivt==vtrks.begin()) { // ...Init them upon 1st track in vertex
	  vertexType = 0;
#  if Kalman_DEBUG_ELoss > 1
	  iTrk = 0;
#  endif
	}
#endif

	// ***** ASSOCIATED MC *****
	CsMCTrack* trkMC = const_cast<CsMCTrack*>(trk->getAssociatedMCTrack());
	if (trkMC) {
	  const list<CsMCHit*>& hits = trkMC->getMCHits();

#ifdef Kalman_DEBUG_ELoss
	  // 1ST TRACK: RETRIEVE some VERTEX' ATTRIBUTES and remember them
	  if (vertexType==0) {
	    if (zTrk<Z_Target/10) {
	      vertexType = 1;   //                   ***** PRIMARY VERTEX? *****
	      zMC = CsEvent::Instance()->getMCVertices().front()->getZ(); 
	    }
	    else {
	      const CsMCVertex *vtxMC = trkMC->getInVertex();
	      if (vtxMC!=NULL) {
		vertexType = -1;
		zMC = vtxMC->getZ();
	      }
	    }
	    if (vertexType) {
	      //                                           ***** z BINNING *****
	      if      (zMC<zs[0]) zBin = 0;
	      else {
		for (zBin = 1; zBin<9; zBin++) {
		  if (zMC<zs[zBin]) break;
		}
	      }
	    }
	  }
#endif

	  // MC and TRACKING MOMENTA @ VERTEX & @ 1ST (or last) MEASURED POINT
          double pV_MCX = trkMC->getPX(), pV_MCY = trkMC->getPY();
          double pV_MCZ = trkMC->getPZ();
          double pV_MC  = sqrt( pV_MCX*pV_MCX + pV_MCY*pV_MCY + pV_MCZ*pV_MCZ );
	  double pR_MC, pR_T;
	  if (zTrk<Z_Target/10) {  // First track == beam track
	    // COMGeant has beam track gain energy! => Correct for this
	    pR_MC = 2*pV_MC-hits.back()->getP().mag(); // (refine: use .back()!)
	    pR_T = fabs( 1/trk->getHelices().back().getCop() );
	  }
	  else {
	    pR_MC = hits.front()->getP().mag();
	    pR_T = fabs( 1/trk->getHelices().front().getCop() );
	  }
	  // ``VERTEXING'' MOMENTA @ VERTEX (i.e. corrected), @ 1ST (last) POINT
	  double pR_V = fabs( 1/(*ivt).Pnk(5) );
	  double pV_V = (*ivt).getMom();
	  
	  int Q = trkMC->getParticle()->getCharge();
	  
	  h0->Fill( pV_MC, ( pV_MC - pR_T )/pV_MC );
	  h1->Fill( pV_MC, ( pV_MC - pV_V )/pV_MC );
	  h2->Fill( pV_MC, ( pV_MC - pR_V )/pV_MC );
	  h3->Fill( pV_MC, ( pV_MC - pR_MC)/pV_MC );
	  h4->Fill( pV_MC, ( pR_MC - pR_T )/pR_MC );
	  h5->Fill( pV_MC, ( pR_MC - pR_V )/pR_MC );
	  
	  h6->Fill( Q/pV_MC, ( pV_MC - pR_T )/pV_MC );
	  h7->Fill( Q/pV_MC, ( pV_MC - pV_V )/pV_MC );
	  
	  hEL_MC->Fill( pV_MC - pR_MC );
	  if( pV_MC < 20 ) hEL20_MC->Fill( pV_MC - pR_MC );

#ifdef Kalman_DEBUG_ELoss

	  // ***** EXTRA HISTOs FOR DEBUG SPECIAL HISTOGRAMMING *****

	  if (vertexType) {  // If vertexType has yet been determined
	    // Some more histograms: The 3 quantities that really matters are
	    // the errors (i.e. MC - reco'd)
	    //   i) @ 1st measured point, where reco means tracking,
	    //     in order to single out that part of E loss, or any other
	    //     error generating phenomenon for that matter, that does
	    //     not arise from target => 1 TH2D
	    //  ii) @ vertex, where reco means vertexing, before correction.
	    //     => 4 TH2D's for 4 different bins in the target thickness
	    // iii) @ vertex, after correction => 4 TH2D's
	    //  In addition, the cases of the beam track and scatterd muon are
	    // dissociated from the rest.
	    //  Also dissociated: secondaries, for they are expected to be less
	    // realiable, and we don't want them to blur the picture.
	    if (vertexType==1) {
	      if      (zTrk<Z_Target/10) {                       // beam
		p1_b->Fill        (   zMC, pR_MC - pR_T );
		p2_b->Fill        (   zMC, pV_MC - pR_V );
		p3_b->Fill        (   zMC, pV_MC - pV_V );
		p4_b->Fill        (   zMC, pV_MC - pR_MC );
		p5_b->Fill        (   zMC, pR_V  - pR_T );
		p6_b->Fill        (   zMC, EL );
		d1_b->Fill        ( ( pR_MC - pR_T )/pR_MC );
		d2_b->Fill        ( ( pV_MC - pR_V )/pV_MC );
		d3_b->Fill        ( ( pV_MC - pV_V )/pV_MC );
		d4_b->Fill        ( ( pV_MC - pR_MC)/pV_MC );
		d6_b->Fill        ( ( EL )/pV_MC );
		ee_b->Fill        ( extrapError );
	      }
	      else if ((*ivt).getStatus()==5) {                  // ...mu'...
		//            ***** p BINNING *****
		int pBin = 0;
		if (pV_MC>ps[0]) {
		  for (pBin = 1; pBin<4; pBin++) {
		    if (pV_MC<ps[pBin]) break;
		  }
		}
		p1_m[0]->Fill     (   zMC, pR_MC - pR_T );
		p2_m[0]->Fill     (   zMC, pV_MC - pR_V );
		p3_m[0]->Fill     (   zMC, pV_MC - pV_V );
		p4_m[0]->Fill     (   zMC, pV_MC - pR_MC );
		p5_m[0]->Fill     (   zMC, pR_V  - pR_T );
		p6_m[0]->Fill     (   zMC, EL );
		d1_m[0]->Fill     ( ( pR_MC - pR_T )/pR_MC );
		d2_m[0]->Fill     ( ( pV_MC - pR_V )/pV_MC );
		d3_m[0]->Fill     ( ( pV_MC - pV_V )/pV_MC );
		d4_m[0]->Fill     ( ( pV_MC - pR_MC)/pV_MC );
		d6_m[0]->Fill     ( ( EL )/pV_MC );
		if (0<pBin && pBin<4) {
		  p1_m[pBin]->Fill(   zMC, pR_MC - pR_T );
		  p2_m[pBin]->Fill(   zMC, pV_MC - pR_V );
		  p3_m[pBin]->Fill(   zMC, pV_MC - pV_V );
		  p4_m[pBin]->Fill(   zMC, pV_MC - pR_MC );
		  p5_m[pBin]->Fill(   zMC, pR_V  - pR_T );
		  p6_m[pBin]->Fill(   zMC, EL );
		  d1_m[pBin]->Fill( ( pR_MC - pR_T )/pR_MC );
		  d2_m[pBin]->Fill( ( pV_MC - pR_V )/pV_MC );
		  d3_m[pBin]->Fill( ( pV_MC - pV_V )/pV_MC );
		  d4_m[pBin]->Fill  ( ( pV_MC - pR_MC)/pV_MC );
		  d6_m[pBin]->Fill  ( ( EL )/pV_MC );
		}
	      }
	      else {                                  // ...primary ``pions''...
		p1->Fill          ( pR_MC, pR_MC - pR_T );
		p2[0]->Fill       ( pV_MC, pV_MC - pR_V );
		p3[0]->Fill       ( pV_MC, pV_MC - pV_V );
		p4[0]->Fill       ( pV_MC, pV_MC - pR_MC );
		p5[0]->Fill       ( pV_MC, pR_V  - pR_T );
		p6[0]->Fill       ( pV_MC, EL );
		d1->Fill          ( ( pR_MC - pR_T )/pR_MC );
		d2[0]->Fill       ( ( pV_MC - pR_V )/pV_MC );
		d3[0]->Fill       ( ( pV_MC - pV_V )/pV_MC );
		d4[0]->Fill       ( ( pV_MC - pR_MC)/pV_MC );
		d6[0]->Fill       ( ( EL )/pV_MC );
		if (0<zBin && zBin<5) {
		  p2[zBin]->Fill  ( pV_MC, pV_MC - pR_V );
		  p3[zBin]->Fill  ( pV_MC, pV_MC - pV_V );
		  p4[zBin]->Fill  ( pV_MC, pV_MC - pR_MC );
		  p5[zBin]->Fill  ( pV_MC, pR_V  - pR_T );
		  p6[zBin]->Fill  ( pV_MC, EL );
		  d2[zBin]->Fill  ( ( pV_MC - pR_V )/pV_MC );
		  d3[zBin]->Fill  ( ( pV_MC - pV_V )/pV_MC );
		  d4[zBin]->Fill  ( ( pV_MC - pR_MC)/pV_MC );
		  d6[zBin]->Fill  ( ( EL )/pV_MC );
		}
	      }
	    }
	    else {
	      p1_s[0]->Fill       ( pR_MC, pR_MC - pR_T );
	      p2_s[0]->Fill       ( pV_MC, pV_MC - pR_V );
	      p3_s[0]->Fill       ( pV_MC, pV_MC - pV_V );	      
	      p4_s[0]->Fill       ( pV_MC, pV_MC - pR_MC );
	      d1_s[0]->Fill       ( ( pR_MC - pR_T )/pR_MC );
	      d2_s[0]->Fill       ( ( pV_MC - pR_V )/pV_MC );
	      d3_s[0]->Fill       ( ( pV_MC - pV_V )/pV_MC );	      
	      if (0<zBin && zBin<9) {                       // ...secondaries
		p1_s[zBin]->Fill  ( pR_MC, pR_MC - pR_T );
		p2_s[zBin]->Fill  ( pV_MC, pV_MC - pR_V );
		p3_s[zBin]->Fill  ( pV_MC, pV_MC - pV_V );	      
		p4_s[zBin]->Fill  ( pV_MC, pV_MC - pR_MC );
		d1_s[zBin]->Fill  ( ( pR_MC - pR_T )/pR_MC );
		d2_s[zBin]->Fill  ( ( pV_MC - pR_V )/pV_MC );
		d3_s[zBin]->Fill  ( ( pV_MC - pV_V )/pV_MC );	      
	      }
	    }
#  if Kalman_DEBUG_ELoss > 2
	    printf("Evt %d Trk %d  %s  %.1f(%d)  %.3f  %.3f  %.3f\n",
		   CsEvent::Instance()->getEventNumberInRun(),iTrk,
		   zTrk<0?"beam":(ivt->getStatus()==5?"mu' ":"    "),
		   zMC,zBin,
		   (pR_MC-pR_T)/pR_MC,(pV_MC-pR_V)/pV_MC,(pV_MC-pV_V)/pV_MC);
	    iTrk++;
#  endif
	  }  // End of detailed diagnosis (enabled when "vertexType" is defined)
#endif
        }  // End associated MC track found
      }  // End MC diagnosis block
    }  //  End histogramming block

  }  // End loop on tracks in vertex

  return;
}










////////////// PRINTS. ARE LOCAL FOR THIS FILE. ///////////////////// 

namespace {

  inline void Print1( int iact )
  {
    cout<<"CsKF::Kalman==> IACT="<<iact<<endl;
    cout<<"CsKF::Kalman==> ****************** KALMAN FILTER ******************"<<endl;
  }

  inline void Print2( list<CsVTrack>& vtrks )
  {
    cout<<endl<<"CsKF::Kalman==> ------- GLOBAL FIT -------"<<endl;
    cout<<"CsKF::Kalman==> TRACKS BEFORE GL FIT:"<<endl;
    list<CsVTrack>::iterator iTr = vtrks.begin();
    for(int k=0; iTr != vtrks.end(); iTr++,k++){
      if( (*iTr).primary() ){
        vector<CsHelix> tpar=(*iTr).getAssociatedTrk()->getHelices();
        cout<<"CsKF::Kalman==> "<<k<<", "<<setprecision(3)<<setw(4)
            <<tpar[0].getZ()/10<<" cm, dX/dZ="
            <<setw(10)<<tpar[0].getDXDZ()<<", dY/dZ="
            <<setw(10)<<tpar[0].getDYDZ()<<", Mom="
            <<setw(5)<<1/tpar[0].getCop()<<"+-"<<sqrt( tpar[0](5,5) )/tpar[0].getCop()/tpar[0].getCop()<<endl;
      }
    }
  }

  inline void Print3(list<CsVTrack> &vtrks, double vChi2)
  {
    cout<<endl<<"CsKF::Kalman==> AFTER GL FIT: chi2="<<vChi2<<endl;
    list<CsVTrack>::iterator iTr = vtrks.begin();
    for(int k=0; iTr != vtrks.end(); iTr++,k++){
      if( (*iTr).primary() ){
        vector<CsHelix> tpar=(*iTr).getAssociatedTrk()->getHelices();
        cout<<"CsKF::Kalman==> "<<k<<", "<<setprecision(3)<<setw(4)
            <<tpar[0].getZ()/10<<" cm, dX/dZ="
            <<setw(10)<<(*iTr).getDXDZ()<<", dY/dZ="
            <<setw(10)<<(*iTr).getDYDZ()<<", Mom="
            <<setw(5)<<(*iTr).getMom()<<"+-"<<sqrt( (*iTr).Dnk(3,3) )/(*iTr).getCop()/(*iTr).getCop()<<endl;
      }
    }
  }

  inline void Print4(CsVTrack &vtrk, double &Chi2trk, TMtx &Xnew,
		     const char *intro = 0)
  {
    if (intro) printf("%s ",intro);

    const CsTrack *trk = vtrk.getAssociatedTrk(); 
    const vector<CsHelix> &vect= trk->getHelices();
    double MOM = 1/vect[0].getCop(), mom = 1/vtrk.getCop();
    printf("chi2Incr=%7.3f z=%7.2f",Chi2trk,Xnew(1));
    if (vtrk.getStatus()==5) printf(" P:\033[34m%#7.4g->%#7.4g\033[m",MOM,mom);
    else                     printf(" P:%#7.4g->%#7.4g",MOM,mom);

    // ********** ---- MC - MC - MC - MC ---- **********
    // - Print particle's name
    // - And its mother's name
    // - And it's radial distance to primary vertex if it's the grand-daughter
    //   - of the pVertex, e.g. K|pi<-D0
    //   - or of a 0-lifetime grand-daughter of the pVertex, e.g. K|pi<-D0<-D*.
    if (CsEvent::Instance()->isAMonteCarloEvent()) {
      // (One need const_cast here because CsTrack::getParticle is not const
      // (as of 2004/07))
      CsMCTrack* trkMC = const_cast<CsMCTrack*>(trk->getAssociatedMCTrack());
      if (trkMC) {
	printf(" (%-8s",trkMC->getParticle()->getName());
	CsMCVertex *mVtx = const_cast<CsMCVertex*>(trkMC->getInVertex());
        if (mVtx->getGnum()==1) printf(" from primary)\n");
        else {
	  CsMCTrack *mTrk = const_cast<CsMCTrack*>(mVtx->getInTrack());
	  if (mTrk) {
	    double xm = mVtx->getX(), ym = mVtx->getY(), D = -1;
	    CsMCVertex *gmVtx = const_cast<CsMCVertex*>(mTrk->getInVertex());
	    double xgm = gmVtx->getX(), ygm = gmVtx->getY();
	    if (gmVtx->getGnum()==1)
	      D = sqrt((xm-xgm)*(xm-xgm)+(ym-ygm)*(ym-ygm));
	    else {
	      CsMCTrack *gmTrk = const_cast<CsMCTrack*>(gmVtx->getInTrack());
	      if (gmTrk) {
		CsMCVertex *ggmVtx = const_cast<CsMCVertex*>(gmTrk->getInVertex());
		double xggm = gmVtx->getX(), yggm = gmVtx->getY();
		if (ggmVtx->getGnum()==1 &&
		    sqrt((xgm-xggm)*(xgm-xggm)+(ygm-yggm)*(ygm-yggm))<1e-5)
		  // If grand-mother is from primary and readily decays...
		  D = sqrt((xm-xgm)*(xm-xgm)+(ym-ygm)*(ym-ygm));
	      }
	    }
	    if (D>=0) printf(" from %d(%.3f)=%s)\n",mVtx->getGnum(),D,
			     mTrk->getParticle()->getName());
	    else      printf(" from %d=%s)\n",mVtx->getGnum(),
			     mTrk->getParticle()->getName());
	  }
	  else
	    printf(" from %d=?!?)\n",mVtx->getGnum());
	}
      }
    }
    else printf("\n");
    //------ END OF MC ---------------- 

  }

  inline void Print6( double vChi2, double* Chi2trk, int iact )
  {
    printf("CsKF::Kalman==> chi2/NDF=%6.3f, sorted tracks=",vChi2);
    for (int j = 0; j<iact-1; j++) cout<<Chi2trk[j]<<", ";
    cout<<endl<<endl;
  }

  inline void Print7( TMtx &Xn, TMtx &Cn, int iact, int Niter, double Chi2tot_ )
  {
    cout<<"CsKF::Kalman==> *********** END of FIT OF PRIMARY ************"<<endl;
    cout<<setprecision(4);
    double XXMC(0),XYMC(0),XZMC(0);

    cout<<endl;
    cout<<"CsKF::Kalman==> Niter = "<<Niter<<", iact = "<<iact
        <<", chi2= "<<Chi2tot_<<", chi2/NDF = "<<Chi2tot_/(2*iact-3)<<endl;
    cout<<"CsKF::Kalman==> Sig X = "<<sqrt(Cn(2,2))
                     <<", Sig Y = "<<sqrt(Cn(3,3))<<", Sig Z = "<<sqrt(Cn(1,1))<<endl;

    if( CsEvent::Instance()->isAMonteCarloEvent() ) {
      CsGeant3* geant3 = CsGeant3::Instance();
      list<CsMCVertex*>mcvertices = geant3->getMCVertices();

      XXMC=mcvertices.front()->getX()/10;
      XYMC=mcvertices.front()->getY()/10;
      XZMC=mcvertices.front()->getZ()/10;

      cout<<endl;
      cout<<"CsKF::Kalman==> Vertex Pull = "<<( XZMC-Xn(1) ) / sqrt(Cn(1,1))<<endl;
      cout<<"CsKF::Kalman==> MC:        z="<<setw(8)<<XZMC       <<", x="
          <<setw(8)<<XXMC <<", y="<<setw(8)<<XYMC <<endl;
    }
    cout<<"CsKF::Kalman==> KALMAN:    z="<<setw(8)<<Xn(1)     <<", x="
        <<setw(8)<<Xn(2)<<", y="<<setw(8)<<Xn(3)<<endl;
    if( CsEvent::Instance()->isAMonteCarloEvent() ) {
      cout<<"CsKF::Kalman==> MC-KALMAN: z="<<setw(8)<<XZMC-Xn(1)<<", x="
          <<setw(8)<<XXMC-Xn(2)<<", y="<<setw(8)<<XYMC-Xn(3)<<endl;
    }

    cout<<endl;    
    Cn.Print(" FINAL Cn ");
  }

  inline void Print8( TMtx &Xn, TMtx &Cn, int iact, int Niter, double Chi2tot_ )
  {
    cout<<"CsKF::Kalman==> ************ END of FIT of SECONDARY ************"<<endl;
    cout<<setprecision(4);
    //double XXMC,XYMC,XZMC;
    if( CsEvent::Instance()->isAMonteCarloEvent() ) {
      CsGeant3* geant3 = CsGeant3::Instance();
      list<CsMCVertex*>mcvertices = geant3->getMCVertices();

      //XXMC=mcvertices.front()->getX()/10;
      //XYMC=mcvertices.front()->getY()/10;
      //XZMC=mcvertices.front()->getZ()/10;

      cout<<endl;
      cout<<"CsKF::Kalman==> Niter = "<<Niter<<", iact = "<<iact
          <<", Chi2tot_ = "<<Chi2tot_<<endl;
      cout<<"CsKF::Kalman==> Chi2tot_/NDF = "<<Chi2tot_/(2*iact-3)
          <<",  Sigma Z = "<<sqrt(Cn(1,1))<<endl;
      cout<<endl;
      /* 
      cout<<"CsKF::Kalman==> MC:        z="<<setw(8)<<XZMC       <<", x=" 
          <<setw(8)<<XXMC <<", y="<<setw(8)<<XYMC <<endl; 
      */
    }
    cout<<"CsKF::Kalman==> KALMAN:    z="<<setw(8)<<Xn(1)     <<", x="
        <<setw(8)<<Xn(2)<<", y="<<setw(8)<<Xn(3)<<endl;
    /* 
    if( CsEvent::Instance()->isAMonteCarloEvent() ) { 
      cout<<"CsKF::Kalman==> MC-KALMAN: z="<<setw(8)<<XZMC-Xn(1)<<", x=" 
          <<setw(8)<<XXMC-Xn(2)<<", y="<<setw(8)<<XYMC-Xn(3)<<endl; 
    } 
    
    cout<<endl;
    cout<<"CsKF::Kalman==> Vertex Pull = "<<( XZMC-Xn(1) ) / sqrt(Cn(1,1))<<endl;
    */
    cout<<endl;
    Cn.Print(" FINAL Cn ");
  }

}
