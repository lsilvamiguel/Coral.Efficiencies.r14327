// $Id: TTrackGetSmoothed.cc 13355 2012-04-14 23:58:47Z ybedfer $

#include <iostream>
#include <iomanip>
#include "TEv.h"
#include "TSetup.h"
#include "TTrack.h"
#include "THit.h"
#ifdef COMPASS_USE_OSPACE_STD
# include <ospace/std/algorithm>
#else
# include <algorithm>  
#endif

using namespace std;

/*!

  Calculate "smoothed" track parameters
  on the plane # ipl and returns Chi2/ndf.

 \param Hout  - output smoothed helix
 \param ipl - input: plane number where smoothed
  track parameters has to be calculated
 \param flg - input: if == false, measurement on the current plane
 is excluded. If there is no hit on the plane ipl, 
 this flag has no effect.

*/

double TTrack::GetSmoothed(THlx& Hout, int ipl, bool flg) const
{

  bool print(false);

  if(mHef.size() == 0 || mHeb.size() == 0 || mHuf.size() == 0 || mHub.size() == 0){
    cout<<"TTrack::GetSmoothed(THlx&, int) ==> FullKF function must be called before for this track: "<<Id<<endl;
    return(-1);
  } 

  if((3*mHef.size() - mHeb.size() - mHuf.size() - mHub.size()) != 0){
    cout<<"TTrack::GetSmoothed(THlx&, int) ==> Track "<<Id<<" has different mH sizes : "
	<<mHef.size()<<" "<<mHeb.size()<<" "<<mHuf.size()<<" "<<mHub.size()<<endl;
    return(-1);
  } 

  if(ipl < lPlnRef.front() || ipl > lPlnRef.back()){
    cout<<"TTrack::GetSmoothed(THlx&, int) ==> can't smooth track ID "<<Id<<" on plane # "<<ipl
	<<" while track's first plane: "<< lPlnRef.front()
	<<" and last plane: "<<lPlnRef.back()<<endl;
    return(-1);
  }

  if(find(lPlnRef.begin(), lPlnRef.end(), ipl) ==  lPlnRef.end()){
    cout<<"TTrack::GetSmoothed(THlx&, int) ==> plane # "<<ipl<<" IS NOT in planes' list of this track : "<<Id<<endl;
    return(-1);
  }
  
  list<int>::const_iterator ipr = lPlnRef.begin(); 
  list<int>::const_iterator ihp = lHitPat.begin(); 
  while(ipr != lPlnRef.end() && (*ipr) != ipl) {ipr++; ihp++;}
  int ihit = (*ihp); // hit # on the plane "ipl" (if any)


  double chi2=0;
  map<int, THlx>::const_iterator im;

  im = mHef.find(ipl);
  if(im == mHef.end()){
    cout<<"TTrack::GetSmoothed(THlx&, int) ==> for track "<<Id<<
      " no extrapolated forward helix on plane # "<<ipl<<endl;
    return(-1);
  }
  im = mHeb.find(ipl);
  if(im == mHeb.end()){
    cout<<"TTrack::GetSmoothed(THlx&, int) ==> for track "<<Id<<
      " no extrapolated backward helix on plane # "<<ipl<<endl;
    return(-1);
  }
  im = mHuf.find(ipl);
  if(im == mHuf.end()){
    cout<<"TTrack::GetSmoothed(THlx&, int) ==> for track "<<Id<<
      " no updated forward helix on plane # "<<ipl<<endl;
    return(-1);
  }
  im = mHub.find(ipl);
  if(im == mHub.end()){
    cout<<"TTrack::GetSmoothed(THlx&, int) ==> for track "<<Id<<
      " no updated backward helix on plane # "<<ipl<<endl;
    return(-1);
  }

  //
  // Update "left" helix by "right" helix
  //

  THlx Hleft, Hright;
  im = mHef.find(ipl);
  Hleft = (*im).second;    // extrapolated forward
  if(ihit < 0 || flg == false) { // if hit is not there or "exclude hit" mode 
    im = mHeb.find(ipl);
    Hright = (*im).second; // extrapolated backward
  } else {                       // with hit
    im = mHub.find(ipl);
    Hright = (*im).second; // updated backward
  }

  if(print) {
    cout<<endl<<endl<<" ---- Track "<<Id<<" smoothing on the plane "<<ipl<<" with hit "<<ihit<<" ----- "<<endl;
    Hright.Print("Right helix");
    Hleft. Print("Left helix");
    (Hright-Hleft).Print("Right helix - Left helix");
  }
  

  Hleft.Update(Hright, Hout, chi2); // update helix by helix 

  if(print) {
    Hout.  Print("Smoothed helix");
    cout<<"Chi2/5 = "<<(chi2/5.)<<endl;
  }


  return(chi2/5.);
}

/*!

  Calculate "smoothed" track parameters
  at given X and returns Chi2/ndf.
  If track do not cross plane X, returns -1;

 \param Hout  - output smoothed helix
 \param X - input: coordinate along the beam, where smoothed 
  track parameters has to be calculated
*/
double TTrack::GetSmoothed(THlx& Hout, double X) const
{
  bool print(false);

  // track do not cross plane X
  if(X < this->Hfirst(0) || X > this->Hlast(0)) {
    if(print) cout<<"TTrack::GetSmoothed(THlx&, double) ==> X = "<<X
		  <<" not between "<<Hfirst(0)<<" and "<<Hlast(0)<<endl;
    return(-1);
  }

  if(mHef.size() == 0 || mHeb.size() == 0 || mHuf.size() == 0 || mHub.size() == 0){
    cerr<<"FATAL: TTrack::GetSmoothed(THlx&, double) ==> FullKF function must be called before for this track: "<<Id<<endl;
    exit(1);
  } 

  if((3*mHef.size() - mHeb.size() - mHuf.size() - mHub.size()) != 0){
    cerr<<"FATAL: TTrack::GetSmoothed(THlx&, double) ==> Track "<<Id<<" has different mH sizes : "
	<<mHef.size()<<" "<<mHeb.size()<<" "<<mHuf.size()<<" "<<mHub.size()<<endl;
    exit(1);
  } 

  
  map<int, THlx>::const_iterator imf, imb;
  imf = this->mHuf.begin();
  imb = this->mHub.begin();
  
  imb++; // advance "updated-backward" iterator
  if(imb == mHub.end()){
    cerr<<"FATAL: TTrack::GetSmoothed(THlx&, double) ==>  abnormal mHub size :"<<mHub.size()<<endl;
    exit(1);
  }

  bool found=false;
  // tolerance to possible difference in helix H(0) if to come to a plane from left and from right
  // (to fix the case when X is exactly detector position)
  const double tol(1.E-4); 
  THlx Hl,Hr;
  do{
    Hl = (*imf).second;
    Hr = (*imb).second;
    if(Hl(0)-tol < X && X < Hr(0)+tol) {
      found=true;
      break;
    }
    imf++; imb++;
  } while (imf!=mHuf.end() && imb!=mHub.end());
  
  if(!found){
    cerr.precision(16);
    cerr<<"FATAL: TTrack::GetSmoothed(THlx&, double) ==> Trk #"<<Id<<" Can't find closest to X = "<<X<<"[cm] helices"<<endl;
    cerr<<"                                       Xfirst = "<<Hfirst(0)<<" Xlast = "<<Hlast(0)<<endl;
    imf = this->mHuf.begin();
    imb = this->mHub.begin();
    imb++; 
    do{
      Hl = (*imf).second;
      Hr = (*imb).second;
      cerr<<"                                       Xhuf = "<<Hl(0)<<" Xhub = "<<Hr(0)<<endl;
      imf++; imb++;
    } while (imf!=mHuf.end() && imb!=mHub.end());
    exit(1);
  }

  if(print) {
    cout<<endl<<endl<<" ---- Track "<<Id<<" smoothing at "<<X<<" ----- "<<endl;
    Hl.Print("Left  helix");
    Hr.Print("Right helix");
  }

  // Extrapolate to requested X
  THlx Hel,Her;
  Hel(0)=X, Her(0)=X;
  
  Hl.Extrapolate(Hel);
  Hr.Extrapolate(Her);

  if(print) {
    Hel.Print("Extrapolated left  helix");
    Her.Print("Extrapolated right helix");
  }

  // get smoothed helix
  double chi2;
  Hel.Update(Her, Hout, chi2); // update helix by helix 

  if(print) {
    Hout.Print("Smoothed helix");
    cout<<"Chi2/5 = "<<(chi2/5.)<<endl;
  }


  return(chi2/5.);
    
  

}


