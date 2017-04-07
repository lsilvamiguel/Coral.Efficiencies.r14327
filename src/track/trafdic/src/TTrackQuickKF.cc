// $Id: TTrackQuickKF.cc 13148 2011-12-28 16:55:25Z kbicker $

#include <iostream>
#include <iomanip>
#include "TEv.h"
#include "TSetup.h"
#include "TTrack.h"
#include "THit.h"

using namespace std;

/*!
  Quick (no mult. scattering corrections,  no hit re-assignment here) Kalman fit through already found hits.

  \param dir  fit direction
  \param mode defines how to fit 

  \verbatim
  Input:

  dir =  1 - Fit in downstream direction (resulting helix is Hlast )
  dir = -1 - Fit in   upstream direction (resulting helix is Hfirst)
 
  mode =  0 - Straight line fit
  mode =  1 - Start from known helix  (first or last, depending on direction), scale up covariance matrix and do the fit
  mode =  2 - start from fully uncertain (except of momentum) helix (for debug only)

\endverbatim

*/
bool TTrack::QuickKF(int dir, int mode)
{

  const TSetup& setup = TSetup::Ref();
  TEv& ev = TEv::Ref();
  const vector<THit>& vHit = ev.vHit();
  
  bool ok=true; bool print=false;

  if(dir != 1 && dir != -1){
    cout<<"TTrack::QuickKF ==> wrong direction argument : "<<dir<<endl;
    assert(false);
  }
  if(lPlnRef.size() != lHitPat.size()){
    cout<<"TTrack::QuickKF ==> lPlnRef.size() != lHitPat.size() "<<endl;
    assert(false);
  }

  if(dir == -1){
    lPlnRef.reverse();
    lHitPat.reverse();
    swap(Hfirst,Hlast);
  }

  THlx Hext, Hupd, H;

  list<int>::iterator ip = lPlnRef.begin() ;
  list<int>::iterator ih = lHitPat.begin() ;
  Chi2tot = 0;

  // Initial helix
  if(mode == 0){ // straight line fit (no momentum)
    H(0) = setup.iPlane2Detect(*ip).X(0); // X of initial helix
    H(1)=H(2)=H(3)=H(4)=0;
    H(1,1)=H(2,2)=50.*50.;
    H(3,3)=H(4,4)=0.5*0.5;
    H(5)=0;    // set infinit momentum
    H(5,5)=1.; // set some big error on q/P
  }
  if(mode == 1){  // fit with momentum
    if( !Hfirst.with_mom() ){
      cout<<"TTrack::QuickKF ==> fit with momentum, but initial momentum is not known"<<endl;
      return(false);
    }
    H=Hfirst;
    H*=100.;     // diagonalize and scale up cov. matrix
  }

  if(mode == 2){ // refit with known momentum (for debug only)
    if( !Hfirst.with_mom() ){
      cout<<"TTrack::QuickKF ==> refit with momentum, but initial momentum not known"<<endl;
      assert(false);
    }
    H(0) = setup.iPlane2Detect(*ip).X(0); // X of initial helix
    H(1)=H(2)=H(3)=H(4)=0;
    H(1,1)=H(2,2)=50.*50.;
    H(3,3)=H(4,4)=0.5*0.5;
    H(5)=Hfirst(5);
    H(5,5)=1.;
  }
  
  if(print) cout<<endl<<endl;

  int nstep(0);
  while(ip != lPlnRef.end() && ih != lHitPat.end()){ // loop over the planes and hits
    nstep++;
    if(print) cout<<"-----> pl "<<(*ip)<<"  proj = "<<setup.vPlane(*ip).IProj<<endl;
    
    // Extrapolate

    Hext(0) = setup.iPlane2Detect(*ip).X(0);
    ok = H.Extrap(Hext(0), Hext);
    if(!ok) {
      //cout<<"TTrack::QuickKF() ==> Extrapolate error for track ID = "
      //  <<Id<<"  dir = "<<dir<<" mode = "<<mode<<" at step # "<<nstep<<endl;
      return(false);
    }
    if(print)Hext.Print("Hext");
    
    // Update

    if((*ih) < 0) { // no hits assosiated with this track on the plane
      Hupd=Hext;
    } else {
      const THit& h = vHit[(*ih)];
      double chi2 = Hext.HitChi2(h); // contribution to Chi2
      Chi2tot+=chi2;
      if(print) {
	cout<<"     hit # "<<(*ih)<<" det = "<<sqrt(1./(h.G0*h.G2 - h.G1*h.G1)) <<" chi2 "<<chi2<<endl;
      }
      Hext.Update(h,Hupd);          // update
    }
    if(print)Hupd.Print("Hupdated");
    if(print) cout<<endl;

    H=Hupd; ip++; ih++;                  // Next plane and next hit
  } // end of loop over planes and hits

  Hlast = H;

  if(dir == -1){
    lPlnRef.reverse();
    lHitPat.reverse();
    swap(Hfirst,Hlast);
  }

  return(true);
}













