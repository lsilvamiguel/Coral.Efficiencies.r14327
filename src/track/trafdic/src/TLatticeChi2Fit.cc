// $Id: TLatticeChi2Fit.cc 13148 2011-12-28 16:55:25Z kbicker $

/*!
  \file   TLatticeChi2Fit.cc
  \brief  Tracks chi2 fit using QNewton and interpolation
  \author Yann Bedfer
  */


/*
 *   Chi2Fit: Determine the track-parameter from the hit coordinates
 *
 *         INPUT:  - pointer to 'Track' structure with the following
 *                   fields defined:
 *                   ->hit_pattern : bit-pattern telling which hits are valid
 *                   ->hit[i]      : hits
 *                   ->par_pattern : 5-bit pattern telling which par's to fit
 *                   ->par         : start-values for the parameters
 *
 *        OUTPUT:  ->hit_pattern : hit-pattern for fitted track
 *                 ->hit[i]      : hits of fitted track
 *                 ->par         : fitted parameters
 */

#include <stdio.h>
#include <math.h>
#include <cstring>   // for memset()

#include "TOpt.h"
#include "TLattice.h"

using namespace std;

extern TLattice *Lat;

// dummy functions, not yet implemented
void get_coords_via_tcp(Track *t) {}
void trafo_track_to_z0(Track *t) {}

float compute_chi_square(Track_Hits *,Track_Hits *);

void TLattice::get_coords_auto(Track *t)
{
  if (TOpt::ReMode[14]&2) {
    int keep = 0;
    Parameters p_keep;
    float      z_keep;
    if (t->z !=0.) {
      /* if the z-co-ord is not 0, use the track-back routine to make it 0. */
      /* to be safe: keep the original parameters */
      p_keep = t->par; z_keep = t->z;
      keep = 1; trafo_track_to_z0(t);
    }
    Lat->get_coords(t);
    if (keep) {
     t->par = p_keep; t->z   = z_keep;
    }
  }    
  else get_coords_via_tcp(t);
}


static float delta[5] = { .25,.5,0.001,0.0001,0.0002 };  /* to give ?mm */

//float sigma[NROWS] = {.01,.01,.01,.01,.01,.01,.01,.01,.01,.01,.01,.01,.01,.01,.01,.01,.01,.01,.01,.01,.01,.01,.01,.01,.01,.01,.01,.01,.01,.01,.01,.01,.01,.01,.01,.01,.01,.01,.01,.01,.01,.01,.01,.01,.01,.01,.01,.01,.01,.01,.01,.01,.01,.01,.01,.01,.01,.01,.01,.01,.01,.01,.01,.01}; /* resolutions */


int TLattice::Chi2Fit(Track *t)
{
  int ret;

  /*#define PXZ*/
#ifdef PXZ
  {
    float tx = t->par.tx, ty = t->par.ty;
    t->par.pinv *= sqrt(1+ty*ty/(1+tx*tx));
  }
#endif

#define FAST_FIT
#ifdef FAST_FIT
  /* switch for fast tracking */
      
  if (TOpt::ReMode[14]&4) ret = Lat->FastChi2Fit(t,1);
  else                    ret = Lat->Chi2Fit_(t);
#else
  ret = Lat->Chi2Fit_(t);
#endif

  if (TOpt::Print[4]>1) if (ret<0) 
    cout<<"** TLattice:Chi2Fit: Error "<<ret<<endl;
#ifdef PXZ
  {
    float tx = t->par.tx, ty = t->par.ty;
    t->par.pinv /= sqrt(1+ty*ty/(1+tx*tx));
  }
#endif
  return ret;
}

/********** (SLOW) CHI SQUARE FIT HERE **********/
int TLattice::Chi2Fit_(Track *t)
{

  Track p1;  /* local copy of track-structure */
  Track_Hits hits_p, hits_m, hits_0;
  float dfdx[NROWS][5];

  float *dd, *pp, chi2;
  float a[5][5],bb[5],b[5];
  int i, j, k, iter;
  
  int old_grid_point = -1, grid_checks = 0;

  memset(dfdx,0,sizeof(dfdx)); memset(a,0,sizeof(a));
  dd = delta; pp = &(p1.par.x); hits_0 = t->hits;
  for (iter=0; iter<10; iter++) {

    /***** LOOP OVER MANY ITERATIONS *****/

    /*
     * For each iteration, compute the lattice-grid-point ONCE
     * to avoid discontinuities when calculating derivatives.
     *  In case of linear interpolation it is also an easy convergence 
     * check by asking whether the grid-point changed.
     *  Anyhow, the calculation of derivatives should be made
     * DIRECTLY from the lattice, which is more accurate and faster.
     */
       
    global_grid_point = Lat->get_point(&t->par);
    if (TOpt::Print[4]>4)
      printf("iteration %d gridpoint %d\n",iter,global_grid_point);

    if (iter && (TOpt::ReMode[14]&2)) {
      if (global_grid_point==old_grid_point) {
	grid_checks++;
	if (grid_checks>1) {
	  if (TOpt::Print[4]>4) printf("Converged after %d iter's!\n",iter);
	  break;
	}
      }
      else grid_checks = 0;
    }
    
    old_grid_point = global_grid_point; 
    
    /***** FILL b WITH FIRST DERIVATIVE, a WITH SECOND DERIVATIVE *****/

    for(i=0;i<5;i++)
      if ((t->par_pattern)&(1<<i)) {
	int j, k, l, m;
	unsigned int hp;
	float chi_p, chi_m;

	/* ADD delta STEP TO PARAMETER */ 

	p1 = *t; pp[i] += dd[i]; Lat->get_coords_auto(&p1); hits_p = p1.hits;

	/* SUBTRACT delta STEP FROM PARAMETER */ 

	p1 = *t; pp[i] -= dd[i]; Lat->get_coords_auto(&p1); hits_m = p1.hits;

	/* COMPUTE CHI2 FOR plus AND minus*/

	chi_p = compute_chi_square( &hits_0 , &hits_p );
	chi_m = compute_chi_square( &hits_0 , &hits_m );

	/* COMPUTE (NEGATIVE) FIRST DERIVATIVE */

	b[i] = (chi_m - chi_p)/dd[i]/2.;

	/* NOW COMPUTE MATRIX OF DERIVATIVES */

	for (l=k = 0; l<NRW32; l++) {
	  hp = hits_0.hp[l];
	  for(m = 0; m<32; m++) {
	    if (hp&(1<<m))
	      dfdx[k][i] = (hits_p.h[k]-hits_m.h[k])/dd[i]/Lat->sigmas[k];
	    k++;
	  }
	}
	hp = hits_0.hp[NRW32];
	for(m = 0; m<NROWS%32; m++) {
	  if (hp&(1<<m))
	    dfdx[k][i] = (hits_p.h[k]-hits_m.h[k])/dd[i]/Lat->sigmas[k];
	  k++;
	}

	for (j = 0; j<=i; j++) {
	  float d2 = 0.;
	  for (l=k = 0; l<NRW32; l++) {
	    hp = hits_0.hp[l];
	    for (m = 0; m<32; m++) {
	      if (hp&(1<<m)) d2 += dfdx[k][i]*dfdx[k][j];
	      k++;
	    }
	  }
	  hp = hits_0.hp[NRW32];
	  for (m = 0; m<NROWS%32; m++) {
	    if (hp&(1<<m)) d2 += dfdx[k][i]*dfdx[k][j];
	    k++;
	  }
	  a[i][j] = a[j][i] = 0.5*d2;
	}
      }

    /* SOLVE THE LINEAR EQUATIONS */

    for (i = 0; i<4; i++)
      if (t->par_pattern&(1<<i))
	for(j = i+1; j<5; j++)
	  if (t->par_pattern&(1<<j)) {
	    float f;
	    if (a[i][i]==0.) return -2;
	    f = a[i][j]/a[i][i];
	    for (k = i+1; k<5; k++) 
	      if (t->par_pattern&(1<<k)) a[k][j] -= f*a[k][i];
	    b[j] -= f*b[i];
	  }

    for (i = 4; i>0; i--)
      if (t->par_pattern&(1<<i))
	for (j = i-1; j>=0; j--)
	  if (t->par_pattern&(1<<j)) {
	    if (a[i][i]==0.) return -3;
	    b[j] -= a[i][j]/a[i][i]*b[i];
	  }

    for (i = 0; i<5; i++)
      if (t->par_pattern&(1<<i)) {
	if (a[i][i]==0.) return -4;
	b[i] /= a[i][i];
      }

    /* APPLY PARAMETER CHANGE */

    p1 = *t;
    for (i = 0; i<5; i++) if (t->par_pattern&(1<<i)) pp[i] += b[i];

    t->par = p1.par;

  } /* for iter */


  /* COMPUTE CHI2 FOR RESULTING PARAMETERS */ 

  p1 = *t;
      
  Lat->get_coords_auto(&p1); hits_p = p1.hits;
  chi2 = compute_chi_square(&hits_0,&hits_p);

  *t = p1;
  t->chi2 = chi2;

  return 0;    
}

/***** COMPUTE CHI2 FOR HITS h1 VERSUS REFERENCE HITS h0 *****/

float compute_chi_square(Track_Hits *h0,Track_Hits *h1)
{
  float chi2, diff;
  int   l, m, k, hp;

  for (l=k = 0, chi2 = 0.; l<NRW32; l++) {
    hp = h0->hp[l];
    for (m = 0; m<32; m++) {
      if (hp&(1<<m)) {
	diff = (h0->h[k]-h1->h[k])/Lat->sigmas[k];
	chi2 += diff*diff;
      }
      k++;
    }
  }
  hp = h0->hp[NRW32];
  for (m = 0; m<NROWS%32; m++) {
    if (hp&(1<<m)) {
      diff = (h0->h[k]-h1->h[k])/Lat->sigmas[k];
      chi2 += diff*diff;
    }
    k++;
  }
  return chi2;
}
