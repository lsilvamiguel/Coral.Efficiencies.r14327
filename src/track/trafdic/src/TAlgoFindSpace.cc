/*!
//
// Preliminary pattern recognition 
// of tracks in space using found projection tracks
//
// input:
//
//        igr   - detector's group #
//        X0    - reference X coordinate
//        npl   - number of planes
//        idpl[ipl] - index of plane #ipl in TSetup::vDetect
//        iprj[ipl] - projection index of plane # ipl
//        cosa[ipl], sina[ipl] - rotation of plane # ipl
//        tol [ipl] - residual tolerance of plane
//        sigt[ipl] - time resolution of plane squared (< 0 means "do not measure time")
//        fh  [ipl] - first hit number (in Uhit[]) on the plane
//        lh  [ipl] - after-the-last hit number (in Uhit[]) on the plane
//                    if there are no nits on the plane, fh[ipl]=lh[ipl]
//        
//        nhits     - total number of hits
//        Uhit[ih]  - hit coordinate in detector reference system (DRS)
//        Thit[ih]  - Time of hit
//        nproj     - number of projections where tracks had been searched
//        Ntrk_prj[ip] - number of tracks in proj ip
//        cos_prj[ip], sin_prj[ip] - angle of proj. ip
//        Y0_prj_[ip][it] - offset of track it in proj. ip
//        Yp_prj[ip][it]  - slope of track it in proj. ip
//        tT_prj_[ip][it] - time of track it in proj. ip
//        sT_prj[ip][it]  - time error of track it in proj. ip (< 0 means "time is not known")
//
//       
// output:
//
//       Ntk       - number of found tracks
//       hits[i]   - array of hit references:
//         [ih0,ih1... ihn,-1,ih0,ih1...ihn,-1,.....]  where ih* = [0,N)
//          \            /    \           /
//           hit of tr.1       hits of tr.2  ...
//
*/
#include <stdio.h>
#include <math.h>
#include <iostream>
#include <set>
#include "TSetup.h"
#include "TOpt.h"
#include "TAlgo.h"
#include "CsHistograms.h"
#include "TConstants.h"
#include "cfortran.h"


PROTOCCALLSFSUB6(SORTZV,sortzv,PFLOAT,PINT,INT,INT,INT,INT)
#define SORTZV(A1,A2,A3,A4,A5,A6)  CCALLSFSUB6(SORTZV,sortzv,PFLOAT,PINT,INT,INT,INT,INT,A1,A2,A3,A4,A5,A6)

/*
#ifdef TIME_CHECK
#include "CsStopwatch.h"
extern double t_algo1, t_algo2;
#endif
*/

int TAlgo::FindSpace(int igr, int ipass, float X0, int npl, int idpl[], int iprj[], float cosa[], float sina[], 
		     float tol[], float sigt[], int fh[], int lh[], int nhit, float Uhit[], float Thit[],
		     int nproj, int Ntk_prj[], float cos_prj[], float sin_prj[], 
		     float **Y0_prj_, float **Yp_prj_, float **tT_prj_, float **sT_prj_, int& Ntk, int hits[])
{
  float (*Y0_prj)[Ntk] = (float(*)[Ntk]) Y0_prj_;
  float (*Yp_prj)[Ntk] = (float(*)[Ntk]) Yp_prj_;
  float (*tT_prj)[Ntk] = (float(*)[Ntk]) tT_prj_;
  float (*sT_prj)[Ntk] = (float(*)[Ntk]) sT_prj_;

  const TSetup& setup = TSetup::Ref();

  const int maxcomb = 70000; // maximum number of combinations
  static float Ysl [maxcomb];
  static float Y0  [maxcomb];
  static float Zsl [maxcomb];
  static float Z0  [maxcomb];
  static float Qt  [maxcomb];
  static float T   [maxcomb];
  static float S   [maxcomb];
  static int  ind  [maxcomb];
  
  float xpl [npl];
  
  short int ifl [nhit];   // flag of hit
  for(int i=0; i < nhit; i++) ifl[i]=0;
 
  for(int i=0; i < npl;  i++) {
    xpl[i] = setup.vDetect(idpl[i]).X(0) - X0;
  }
  
  int ipl,nfound,nused,k,k0,kl,kr;
  double yp,y0,zp,z0,y,z,qual,yy,zz;
  double c1,c2,s1,s2,y1,yp1,y2,yp2,det;
  double e1, e2, t1, t2;
  bool usetime(false);
  int nt1,nt2;
  char id[20];
  std::set<int> snprj;

  int ncomb=0;
  for(int ip1=0; ip1 < nproj-1; ip1++){ // loop over 1-st base proj.
    nt1=Ntk_prj[ip1];
    c1=cos_prj[ip1]; s1=sin_prj[ip1];
  
     
    for(int ip2=ip1+1; ip2 < nproj; ip2++){ // loop over 2-d base proj. 
      nt2=Ntk_prj[ip2];
      c2=cos_prj[ip2]; float s2=sin_prj[ip2];
      det=c1*s2-c2*s1;

      for(int it1=0; it1 < nt1; it1++){ // loop over tracks on 1-st base proj.
	y1 =Y0_prj[ip1][it1];
	yp1=Yp_prj[ip1][it1];
	t1 =tT_prj[ip1][it1];
	e1 =sT_prj[ip1][it1];
		
	for(int it2=0; it2 < nt2; it2++){ // loop over tracks on 2-d base proj.
	  y2 =Y0_prj[ip2][it2];
	  yp2=Yp_prj[ip2][it2];

	  t2 =tT_prj[ip2][it2];
	  e2 =sT_prj[ip2][it2];
	  

	  if(e1 > 0 && e2 > 0) usetime = true;
	  else                 usetime = false;

	  // cut on chi2 of differece of proj. track times
	  if(usetime) {
	    if((t1-t2)*(t1-t2)/(e1 + e2) > TOpt::dCut[20]) continue;
	  }

	  // space track candidate params
	  y0= (y1*s2-s1*y2)/det;
	  z0=-(y1*c2-c1*y2)/det;
	  yp= (yp1*s2-s1*yp2)/det;
	  zp=-(yp1*c2-c1*yp2)/det;

	  //	  std::cout<<"  pr1 = "<<ip1<<" y1= "<<y1<<" yp1 = "<<yp1
	  //    <<"  pr2 = "<<ip2<<" y2= "<<y2<<" yp2 = "<<yp2<<std::endl
	  //     <<"  y0  = "<<y0<<"\t z0 = "<<z0<<"\t yp = "<<yp<<"\t zp = "<<zp<<std::endl;
	  
	  nfound=0;
	  qual=0;
	  snprj.clear();
	  int igap=0;
	  int maxgap = 0;
          for(ipl=0; ipl < npl; ipl++){ // count hits on _all_ planes in defined direction
	    yy=y0+yp*xpl[ipl];
	    zz=z0+zp*xpl[ipl];
	    if(!setup.vDetect(idpl[ipl]).InActive(yy,zz)) continue; // not in active area
	    igap++;

	    if(fh[ipl]==lh[ipl]) continue; // no hits on plane

	    y = yy*cosa[ipl] + zz*sina[ipl]; // rotate
	    if(y < (Uhit[fh[ipl]]  -tol[ipl]) ||
	       y > (Uhit[lh[ipl]-1]+tol[ipl])) continue; // "y" is too far from hits on this plane 

	    // find hit, closest to track position "y"
	    kl=fh[ipl]; kr=lh[ipl]-1;
	    while((kr-kl) > 1){
	      k0=(kl+kr)>>1;
	      if(y <= Uhit[k0]) kr=k0;
	      else              kl=k0;
	    }
	    if((Uhit[kr]-y) < (y-Uhit[kl])) k0=kr;
	    else                            k0=kl;
	    
	    // Cuts chi2 of difference of the hit time and proj. track  times
	    if(usetime && sigt[ipl] > 0){
	      if((t1-Thit[k0])*(t1-Thit[k0])/(e1 + sigt[ipl]) > TOpt::dCut[20]) continue;
	      if((t2-Thit[k0])*(t2-Thit[k0])/(e2 + sigt[ipl]) > TOpt::dCut[20]) continue;
	    }

	    if(fabs(y-Uhit[k0]) < tol[ipl] && ifl[k0] == 0) {
	      nfound++;
	      snprj.insert(iprj[ipl]);
	      qual+=1.0-fabs(y-Uhit[k0])/tol[ipl];
	      if(igap > maxgap) maxgap = igap;
	      igap=0; // reset missed hits counter
	    }
	  } // end of loop over planes 
	  
	  // Space track must contains at least 1 hit of "stereo" (3-d) projection
	  if(snprj.size() < 3) continue;

	  // Too many missed-in-the-row hits
	  if(TOpt::iCut[igr] != 0 && maxgap > TOpt::iCut[igr]) continue;

	  // Cut on N hits n space
	  if(nfound < TOpt::iPRpar[igr*10 + 2]) continue; // check next pair of proj. tracks
	  
	  // Store combination
	  Ysl[ncomb] = yp;
	  Y0 [ncomb] = y0;
	  Zsl[ncomb] = zp;
	  Z0 [ncomb] = z0;
	  Qt [ncomb] = nfound+qual/nfound; // Number of hits + quality factor (less then 1)
	  if(usetime){
	    S[ncomb] = 1./(1./e1 + 1./e2);
	    T[ncomb] = (t1/e1 + t2/e2) * S[ncomb];
	  } else {
	    T[ncomb] =  0;
	    S[ncomb] = -1;
	  }
	  ind[ncomb] = ncomb+1;
	  if( ++ncomb >= maxcomb) {
	    std::cout<<"TEvPrePatRec == > too many combinations "<<ncomb<<std::endl;
	    goto sorting;
	  }
	  
	} // end of loop over tracks on 2-d base proj.
      } // end of loop over tracks on 1-st base proj.
      
    } // end of loop over 2-d base proj
  } // end of loop over 1-st base proj
  
  // sort Qt-wise

 sorting:
  SORTZV(Qt[0],ind[0],ncomb,1,1,0);

  int m; int nt=0, nh=0;
  
  int hh[npl];
  for(int l=0; l < ncomb; l++){ // loop over combinations starting from Qt max 
    m=ind[l]-1;

    nfound = 0; nused = 0;
    
    snprj.clear();
    int igap = 0; 
    int maxgap = 0;
    for(ipl=0; ipl < npl; ipl++){ // count hits on _all_ planes in defined direction
      yy=Y0[m]+Ysl[m]*xpl[ipl];
      zz=Z0[m]+Zsl[m]*xpl[ipl];
      if(!setup.vDetect(idpl[ipl]).InActive(yy,zz)) continue; // not in active area
      igap++;
      
      if(fh[ipl]==lh[ipl]) continue; // no hits on the plane

      y = yy*cosa[ipl] + zz*sina[ipl]; // rotate
      if(y < (Uhit[fh[ipl]]  -tol[ipl]) ||
	 y > (Uhit[lh[ipl]-1]+tol[ipl])) continue; // "y" is too far from hits on this plane  
      
      // find hit, closest to track position "y"
      kl=fh[ipl]; kr=lh[ipl]-1;
      while((kr-kl) > 1){
	k0=(kl+kr)>>1;
	if(y <= Uhit[k0]) kr=k0;
	else              kl=k0;
      }
      if((Uhit[kr]-y) < (y-Uhit[kl])) k0=kr;
      else                            k0=kl;
      
      // Cuts chi2 of difference of the hit time and base proj. track mean times
      if(S[m] > 0 && sigt[ipl] > 0){
	if((T[m]-Thit[k0])*(T[m]-Thit[k0])/(S[m] + sigt[ipl]) > TOpt::dCut[20]) continue;
      }
 
     if(fabs(y-Uhit[k0]) < tol[ipl]) {
	hh[nfound++]=k0;          // strore hit
	snprj.insert(iprj[ipl]);  // store projection index
	if(ifl[k0] !=0) nused++;
	// cut
	if(nused > TOpt::iPRpar[igr*10 + 1]) goto next_comb;
	if(igap > maxgap) maxgap = igap;
	igap=0; // reset missed hits counter
     } 
    } // end of loop over planes 
    
    // Space track must contains at least 1 hit of "stereo" (3-d) projection
    if(snprj.size() < 3) goto next_comb;

    // Too many missed-in-the-row hits
    if(TOpt::iCut[igr] != 0 && maxgap > TOpt::iCut[igr]) goto next_comb;

    //store found hits
    for(int i=0; i < nfound; i++) {
      hits[nh++]=hh[i];
      ifl[hh[i]]=1;  // label hit like used
    }
    hits[nh++]=-1; //end of hit list

    if(++nt >= TConstants_NTtrack_max)  {
      if(TOpt::Print[0] != 0) 
	std::cout<<"TEvPrePatRec == > (space) too many tracks "<<nt<<std::endl;
      return(2);
    }

  next_comb:;
    //std::cout<<"icomb = "<<l<<"out of "<<ncomb<<" igr = "<<igr<<" nhits = "<<nfound<<" maxgap = "<<maxgap<<std::endl;
  } // end of loop over combination

  if(TOpt::Print[5]){
    printf("TAlgoFindSpace ==> in group %3u  N comb. = %5u   N found tracks = %4u \n", igr, ncomb, nt);
  }
  Ntk=nt;  return(0);
}














