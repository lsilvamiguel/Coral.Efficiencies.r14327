/*!
//
// Alignment check histograms
//
// input:
//
//        igr   - detector's group #
//        npl   - number of planes
//        idpl[ipl] - index of plane #ipl in TSetup::vDetect
//        iprj[ipl] - projection index of plane # ipl
//        cosa[ipl], sina[ipl] - rotation of plane # ipl
//        tol [ipl] - residual tolerance of plane
//        fh  [ipl] - first hit number (in Uhit[]) on the plane
//        lh  [ipl] - after-the-last hit number (in Uhit[]) on the plane
//                    if there are no nits on the plane, fh[ipl]=lh[ipl]
//        
//        nhits     - total number of hits
//        Uhit[ih]  - hit coordinate in detector reference system (DRS)
*/

#include <iostream>
#include <stdio.h>
#include <math.h>
#include "TSetup.h"
#include "TOpt.h"
#include "TAlgo.h"
#include "CsHistograms.h"
#include "TConstants.h"
#include "cfortran.h"
#include "TROOT.h"
#include "TH1.h"
#ifdef COMPASS_USE_OSPACE_STD
# include <ospace/std/set>
#else
# include <set>
#endif
# include "TSetup.h"

using namespace std;

void TAlgo::Alignment(int igr, int npl, int idpl[], int iprj[], float cosa[], float sina[], 
		      float tol[], int fh[], int lh[], int nhit, float Uhit[])
{
  
  if(TOpt::Hist[1] == 0) return; // all histograms are off
  if(TOpt::Hist[3] == 0) return; // alignment histograms are OFF
  if(npl <= 4) {
    cout<<"TAlgoAlignment() ==> Not enouph planes in group # "<<igr<<endl;
    return;
  }
  CsHistograms::SetCurrentPath("/Traffic/Alignment");

  string np1L, np1R, np2L, np2R;
  int ip1L(-1),ip1R(-1),ip2L(-1),ip2R(-1); // 4 base planes (out of npl)
  bool print(false);
  // base planes selection


  if(igr == 4){
    np1L="FI01X1"; np1R="FI02X1";
    np2L="FI01Y1"; np2R="FI02Y1";
  }

  if(igr == 0){
    np1L="FI03X1"; np1R="FI04X1";
    np2L="FI03U1"; np2R="FI04U1";
  }

  if(igr == 1){
    np1L="FI05X1"; np1R="FI06X1";
    np2L="FI05Y1"; np2R="FI06Y1";
    //np1L="GM01X1"; np1R="GM06X1";
    //np2L="GM01Y1"; np2R="GM06Y1";
  }

  if(igr == 2){
    np1L="PA03U1"; np1R="PA06U1";
    np2L="PA03V1"; np2R="PA06V1";
    //np1L="GM08X1"; np1R="GM09X1";
    //np2L="GM08Y1"; np2R="GM09Y1";
  }

  if(igr == 3){
    np1L="PB01U1"; np1R="PB05U1";
    np2L="PB02V1"; np2R="PB06V1";
  }


   for(int ii=0; ii < npl; ii++){
     const string& name = TSetup::Ref().vDetect(idpl[ii]).Name;
     if(name.find(np1L) == 0) ip1L=ii;
     if(name.find(np1R) == 0) ip1R=ii;
     if(name.find(np2L) == 0) ip2L=ii;
     if(name.find(np2R) == 0) ip2R=ii;
   }

 // check if base planes are correctly selected
 if(ip1L == ip1R || ip2L == ip2R) {
   cout<<"TAlgoAlignment() ==> Wrong base planes selection for group # "<<igr
       <<" ("<<npl<<" planes)"<<endl;
   cout<<"ip1L = "<<ip1L<<"  ip1R = "<<ip1R<<"  ip2L = "<<ip2L<<"  ip2R = "<<ip2R<<endl;
   assert(false);
 }

 if(ip1L < 0 || ip1R < 0 || ip2L < 0 || ip2R < 0 ) {
   cout<<"TAlgoAlignment() ==> No base planes selection for group # "<<igr
       <<" ("<<npl<<" planes)"<<endl;
   cout<<"ip1L = "<<ip1L<<"  ip1R = "<<ip1R<<"  ip2L = "<<ip2L<<"  ip2R = "<<ip2R<<endl;
   assert(false);
 } 

 if(ip1L >= npl || ip1R >= npl || ip2L >= npl || ip2R >= npl ){
   cout<<"TAlgoAlignment() ==> Wrong base planes selection for group # "<<igr
       <<" ("<<npl<<" planes)"<<endl;
   cout<<"ip1L = "<<ip1L<<"  ip1R = "<<ip1R<<"  ip2L = "<<ip2L<<"  ip2R = "<<ip2R<<endl;
   assert(false);
 }


 if((iprj[ip1L] != iprj[ip1R]) ||
    (iprj[ip2L] != iprj[ip2R]) ||
    (iprj[ip1L] == iprj[ip2L])) {
   cout<<"TAlgoAlignment() ==> Wrong projections selection for base planes. The group # "<<igr<<endl;
   for(int ii=0; ii < npl; ii++){
     cout<<iprj[ii]<<" "<<TSetup::Ref().vDetect(idpl[ii]).Name<<" ("<<ii<<")";
     if(ii == ip1R || ii == ip2R || ii == ip1L || ii == ip2L) cout<<" <----";
     cout<<endl;
   }
   cout<<endl;
   assert(false);
 }

 // book histograms 
 // (explicit use of ROOT just because of need "FindObject" function
 char id[20], tit[90];
 static set<int> grps;
 if(grps.find(igr) == grps.end()){ // "first call" blocks
   grps.insert(igr);
   cout<<"TAlgoAlignment() ==> group "<<igr<<" base planes : "
       <<TSetup::Ref().vDetect(idpl[ip1L]).Name<<TSetup::Ref().vDetect(idpl[ip2L]).Name
       <<TSetup::Ref().vDetect(idpl[ip1R]).Name<<TSetup::Ref().vDetect(idpl[ip2R]).Name<<endl;
   for(int i = 0; i < npl; i++){
     string tbnam = TSetup::Ref().vDetect(idpl[i]).Name;
     double pitch = fabs(TSetup::Ref().vDetect(idpl[i]).Pitch);
     sprintf(id, "alig_%u_%02u",igr,i);
     sprintf(tit,"%s alignment",tbnam.c_str()); 
     new TH1D(id,tit,100, -100*pitch, 100*pitch);
     sprintf(id, "alig_%u_%02u_z",igr,i);
     new TH1D(id,tit,100, -10*pitch, 10*pitch);
     sprintf(id, "alig_%u_%02u_zz",igr,i);
     new TH1D(id,tit,100, -2*pitch, 2*pitch);
     cout<<"... book histogram '"<<tit<<"'";
     if(i == ip1R || i == ip2R || i == ip1L || i == ip2L) cout<<" <----";
     cout<<endl;
   }
   cout<<endl;
 }// end "if"

 if(fh[ip1L] == lh[ip1L] ||
    fh[ip1R] == lh[ip1R] ||
    fh[ip2L] == lh[ip2L] ||
    fh[ip2R] == lh[ip2R]) { // must be at least 1 hit on each of base planes
   if(print) {
     cout<<"TAlgoAlignment() ==> No hits on one of base planes of group # "<<igr<<endl;
     cout<<"Plane ID - N hits : "
	 <<TSetup::Ref().vDetect(idpl[ip1L]).IDet<<"-"<<lh[ip1L]-fh[ip1L]<<" "
	 <<TSetup::Ref().vDetect(idpl[ip1R]).IDet<<"-"<<lh[ip1R]-fh[ip1R]<<" "
	 <<TSetup::Ref().vDetect(idpl[ip2L]).IDet<<"-"<<lh[ip2L]-fh[ip2L]<<" "
	 <<TSetup::Ref().vDetect(idpl[ip2R]).IDet<<"-"<<lh[ip2R]-fh[ip2R]<<" "<<endl;
   }
   return;
 } 

 if((lh[ip1L] - fh[ip1L]) > 1 ||
    (lh[ip1R] - fh[ip1R]) > 1 ||
    (lh[ip2L] - fh[ip2L]) > 1 ||
    (lh[ip2R] - fh[ip2R]) > 1  ) { // if more then 1 hit on base planes, check distances
   int ipl(-1);
   for(int i=0; i < 4; i++){ // loop over base planes
     if(i == 0) ipl = ip1L;
     if(i == 1) ipl = ip1R;
     if(i == 2) ipl = ip2L;
     if(i == 3) ipl = ip2R;
     double pitch = fabs(TSetup::Ref().vDetect(idpl[ipl]).Pitch);
     for(int ihit = fh[ipl]; ihit < lh[ipl]-1; ihit++){
       if(fabs(Uhit[ihit]-Uhit[ihit+1]) < 50*pitch) return; // skip if base plane hits are too close
     }
   }
 } 

 const TSetup& setup = TSetup::Ref();
 double xpl[npl]; 
 double x0=setup.vDetect(idpl[ip1L]).X(0);
 for(int i=0; i < npl;  i++)  {
   xpl[i] = setup.vDetect(idpl[i]).X(0)-x0;
 }

 double dx1 = xpl[ip1R] - xpl[ip1L];
 double dx2 = xpl[ip2R] - xpl[ip2L];
 if( dx1 == 0. || dx2 == 0.){
   cout<<"TAlgoAlignment() ==> Something is wrong in base planes selection for the group # "<<igr<<endl;
   assert(false);
 }

 double yp1, y01, yp2, y02, Y0, YP, Z0, ZP, y, yy,zz, dU;
 int kl,kr,k0, ntrk(0);
 double c1=cosa[ip1L]; double c2=cosa[ip2L]; double s1=sina[ip1L]; double s2=sina[ip2L];
 double det=c1*s2-c2*s1;

 for(int i1L = fh[ip1L]; i1L < lh[ip1L]; i1L++){ // loop over hits on base plane 1-Left
   for(int i1R = fh[ip1R]; i1R < lh[ip1R]; i1R++){ // loop over hits on base plane 1-Right
     //proj. # 1 track candidate
     yp1 = (Uhit[i1R]-Uhit[i1L])/dx1; // slope
     y01 = Uhit[i1L]-yp1*xpl[ip1L]; // offset at x0 

     for(int i2L = fh[ip2L]; i2L < lh[ip2L]; i2L++){ // loop over hits on base plane 2-Left
       for(int i2R = fh[ip2R]; i2R < lh[ip2R]; i2R++){ // loop over hits on base plane 2-Right
	 //proj. # 2 track candidate
	 yp2 = (Uhit[i2R]-Uhit[i2L])/dx2; // slope
	 y02 = Uhit[i2L]-yp2*xpl[ip2L]; // offset at x0
	 
	 // space track candidate
	 Y0= (y01*s2-s1*y02)/det;
	 Z0=-(y01*c2-c1*y02)/det;
	 YP= (yp1*s2-s1*yp2)/det;
	 ZP=-(yp1*c2-c1*yp2)/det;
	 ntrk++;
	 
	 for(int ipl=0; ipl < npl; ipl++){ // loop over ALL planes
	   if(fh[ipl]==lh[ipl]) continue; // no hits
	   // extrapolate to xpl[ipl]
	   yy=Y0+YP*xpl[ipl];
	   zz=Z0+ZP*xpl[ipl];
	   // rotate
	   y = yy*cosa[ipl] + zz*sina[ipl];

	   // fill residuals histogram
	   sprintf(id, "alig_%u_%02u",igr,ipl);
	   TH1D* h1 = ((TH1D*)gROOT->FindObject(id));	   
	   sprintf(id, "alig_%u_%02u_z",igr,ipl);
	   TH1D* h2 = ((TH1D*)gROOT->FindObject(id));	   
	   sprintf(id, "alig_%u_%02u_zz",igr,ipl);
	   TH1D* h3 = ((TH1D*)gROOT->FindObject(id));	   
	   for(int ih = fh[ipl]; ih < lh[ipl]; ih++) { // loop over hits on ipl
	     dU = y - double(Uhit[ih]); // residual: "Expected coordinate - Hit coordinate"
	     h1->Fill(dU);
	     h2->Fill(dU);
	     h3->Fill(dU);
	   } // end of loop over hits on ipl

	 } // end of loop over planes 
	 
       } // end of loop over hits on base plane 2-Right
     } // end of loop over hits on base plane 2-Left
   } // end of loop over hits on base plane 1-Right
 } // end of loop over hits on base plane 1-Left

 if(print) cout<<"alignment for group "<<igr<<" use "<<ntrk<<" space track candidates"<<endl;
 CsHistograms::SetCurrentPath("/");

}













