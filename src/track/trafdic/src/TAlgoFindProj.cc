#include <iostream>
#include <stdio.h>
#include <math.h>
#include "TOpt.h"
#include "TAlgo.h"
#include "cfortran.h"

PROTOCCALLSFSUB6(SORTZV,sortzv,PFLOAT,PINT,INT,INT,INT,INT)
#define SORTZV(A1,A2,A3,A4,A5,A6)  CCALLSFSUB6(SORTZV,sortzv,PFLOAT,PINT,INT,INT,INT,INT,A1,A2,A3,A4,A5,A6)

/*!
  Preliminary pattern recognition 
  in one projection using "pivot planes" method.

  If detectors measure the time, it is used for filtering out
  wrong combinations.

  input:
        igr   - detector's group #
        X0    - reference X coordinate
        npl   - number of planes
        xpl [ipl] - x coodinates of the plane #ipl
        tol [ipl] - residual tolerance of plane
        sigt[ipl] - time resolution on the plane squared (< 0 means "do not use time")
        fh  [ipl] - first hit number (in yhit[]) on the plane
        lh  [ipl] - after-the-last hit number (in yhit[]) on the plane
                    if there are no nits on the plane, fh[ipl]=lh[ipl]
        
        nhits    - total number of hits
        yhit[ih] - hit coordinate in detector reference system (WRS)
        thit[ih] - time of hit

        Ntk      - max possible number of tracks (just for checking))
       
  output:

       Ntk       - number of found tracks
       hits[i]   - array of hit references:
         [ih0,ih1... ihn,-1,ih0,ih1...ihn,-1,.....]  where ih* = [0,N)
          \            /    \           /
           hit of tr.1       hits of tr.2  ...

      Y0_trk[i] - track offsets at X0
      Yp_trk[i] - track slopes

      T_trk[i]  - track mean time (or 0)
      sT_trk[i] - track time sigma (or -1 if not known)
*/

  int TAlgo::FindProj(int igr, int ipass, float X0, int npl, float xpl[], float tol[], float sigt[], int fh[], int lh[], 
		      int nhit, float yhit[], float thit[], int& Ntk, int hits[], float Y0_trk[], float Yp_trk[],
		      float T_trk[], float sT_trk[]) 
{
  const int maxcomb = 70000; // maximum number of combinations
  static float Ysl [maxcomb];
  static float Y0  [maxcomb];
  static float Nh  [maxcomb];
  static int  ind  [maxcomb];
  
  short int ifl [nhit];   // flag of hit
  for(int i=0; i < nhit; i++) ifl[i]=0; 

  int ipl,nfound,nused,k,k0,kl,kr;
  double yp,y0,y1,dy,y,qual,s0,s1,t,ss,w0,w1,t0(0),t1(0);
  double x0,x1,dx;
  bool usetime;

  // pivot planes
  int ipl0;
  int ipl1;
  
  int ncomb=0;

  for(ipl0=0; ipl0 < npl-1; ipl0++){ // loop over 1-st pivot planes 
    if(fh[ipl0]==lh[ipl0]) continue; // no hits

    for(ipl1=npl-1; ipl1 > ipl0; ipl1--){ // loop over 2-d pivot planes
      if(fh[ipl1]==lh[ipl1]) continue; // no hits
      
      x0=xpl[ipl0];
      x1=xpl[ipl1];
      dx=x1-x0;
      if(dx == 0.) continue;
      s0=sigt[ipl0];
      s1=sigt[ipl1];
      if(s0 > 0 && s1 > 0) usetime = true;
      else                 usetime = false;

      for(int ih0=fh[ipl0]; ih0 < lh[ipl0]; ih0++){ // loop over hits on 1-st pivot plane
	y0=yhit[ih0];
	if(usetime) t0=thit[ih0];
	
	for(int ih1=fh[ipl1]; ih1 < lh[ipl1]; ih1++){ // loop over hits on 2-d pivot plane
	  y1=yhit[ih1];
	  if(usetime) t1=thit[ih1];

	  yp=(y1-y0)/dx;
	  if(fabs(yp) > TOpt::dPRpar[igr*10 + 1]) continue;
	  
	  // cut on chi2 of differece of pivot plane hit times
	  if(usetime) {
	    if((t0-t1)*(t0-t1)/(s0+s1) > TOpt::dCut[20]) continue;
	  }
	  nfound=0;
	  qual=0;
	  
	  for(ipl=0; ipl < npl; ipl++){ // count hits on _all_ planes in defined direction
	    if(fh[ipl]==lh[ipl]) continue;
	    y=y0+yp*(xpl[ipl]-x0);

	    if(y < (yhit[fh[ipl]]  -tol[ipl]) ||
	       y > (yhit[lh[ipl]-1]+tol[ipl])) continue; // "y" is too far from hits on this plane  


	    // find hit, closest to track position "y"
	    kl=fh[ipl]; kr=lh[ipl]-1;
	    while((kr-kl) > 1){
	      k0=(kl+kr)>>1;
	      if(y <= yhit[k0]) kr=k0;
	      else              kl=k0;
	    }
	    if((yhit[kr]-y) < (y-yhit[kl])) k0=kr;
	    else                            k0=kl;
	    
	    // Cuts chi2 of difference of the hit time and pivot plane hits times
	    if(usetime && sigt[ipl] > 0){
	      if((t0-thit[k0])*(t0-thit[k0])/(s0 + sigt[ipl]) > TOpt::dCut[20]) continue;
	      if((t1-thit[k0])*(t1-thit[k0])/(s1 + sigt[ipl]) > TOpt::dCut[20]) continue;
	    }
	    
	    if(fabs(y-yhit[k0]) < tol[ipl] && ifl[k0] == 0) {
	      nfound++;
	      qual+=1.0-fabs(y-yhit[k0])/tol[ipl];
	    }
	  } // end of loop over planes 
	  
	  // Cut
	  if(nfound < TOpt::iPRpar[igr*10 + 0]) continue; // check next pair of hits

	  // Store combination
	  Ysl[ncomb] = yp;
	  Y0 [ncomb] = y0 + yp*(X0-x0); // offset at X0
	  Nh [ncomb] = nfound+qual/nfound;;
	  ind[ncomb] = ncomb+1;
	  if( ++ncomb >= maxcomb) {
	    if(TOpt::Print[0] != 0) 
	      std::cout<<"TAlgoFindProj == > too many combinations ("<<ncomb<<") in group "<<igr<<std::endl;
	    goto sorting;
	  }
	  
	} //end of loop over 2-d pivot plane hits 
      }//end of loop over 1-st pivot plane hits
 
    } // end of loop over 2-d pivot planes
  } // end of loop over 1-st pivot planes
  
  // sort Nh-wise

 sorting:
  SORTZV(Nh[0],ind[0],ncomb,1,1,0);

  int m; int nt=0, nh=0;
  float tw(0),w(0); 

  int   hh[npl];
  int   hp[npl];
  for(int l=0; l < ncomb; l++){ // loop over combinations starting from Nh max 
    m=ind[l]-1;

    nfound = 0; nused = 0;
    for(ipl=0; ipl < npl; ipl++){ // count hits on _all_ planes in defined direction
      if(fh[ipl]==lh[ipl]) continue;
      y=Y0[m]+Ysl[m]*(xpl[ipl]-X0);

      if(y < (yhit[fh[ipl]]  -tol[ipl]) ||
	 y > (yhit[lh[ipl]-1]+tol[ipl])) continue; // "y" is too far from hits on this plane  

      // find hit, closest to track position "y"
      kl=fh[ipl]; kr=lh[ipl]-1;
      while((kr-kl) > 1){
	k0=(kl+kr)>>1;
	if(y <= yhit[k0]) kr=k0;
	else              kl=k0;
      }
      if((yhit[kr]-y) < (y-yhit[kl])) k0=kr;
      else                            k0=kl;
      
      if(fabs(y-yhit[k0]) < tol[ipl]) {
	// strore hit
	hp[nfound]  = ipl;
	hh[nfound++]=k0; 
	if(ifl[k0] !=0) nused++;
	// cut
	if(nused > TOpt::iPRpar[igr*10 + 1]) goto next_comb;
      } 
    } // end of loop over planes 

    
    //store found hits
    int j,k;
    tw=0; w=0;
    for(int i=0; i < nfound; i++) { // loop over hits
      j = hh[i]; // hit # 
      k = hp[i]; // pl. #
      hits[nh++]=j;
      ifl[j]=1;  // label hit like used
      // for mean track time
      if(sigt[k] > 0){
	w  += 1./ sigt[k];
	tw += thit[j] / sigt[k];
	//std::cout<<thit[j]<<" +- "<<sigt[k]<<std::endl;
      }
    } // end of loop over found hits

    //if( w != 0.) std::cout<<"  nhit = "<<nfound<<" time = "<<tw/w 
    //<<" sigma = "<< sqrt(1./w)<<" tw = "<<tw<<" w = "<<w<<std::endl;

    hits[nh++]=-1; //end of hit list

    //store found track parameters
    Y0_trk[nt]=Y0[m];
    Yp_trk[nt]=Ysl[m];
    // store time with error (if any)
    if(w !=0.) {
      T_trk [nt]=  tw/w;
      sT_trk[nt]=  1./w;
    } else {
      T_trk [nt]=  0.;
      sT_trk[nt]= -1.;
    }

    if(++nt >= Ntk)  {
      if(TOpt::Print[0] != 0) 
	std::cout<<"TAlgoFindProj == > too many tracks ("<<nt<<") in group "<<igr<<std::endl; 
      return(1);
    }

  next_comb:;
  } // end of loop over combination

  if(TOpt::Print[5]){
    printf("  N comb. = %5u   N found tracks = %4u \n", ncomb, nt);
  }
  Ntk=nt;  return(0);
}














