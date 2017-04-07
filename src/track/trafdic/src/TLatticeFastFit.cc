// $Id: TLatticeFastFit.cc 13148 2011-12-28 16:55:25Z kbicker $

/*!
  \file   TLatticeFastFit.cc
  \brief  Tracks chi2 fit using QNewton and interpolation in a lattice
  \author Yann Bedfer
  */

/*
 * FastFit: determine the track-parameters from the hit coordinates
 *          for a list of tracks
 *
 *         INPUT:  - pointer to a list of tracks where the following
 *                  fields are used:
 *                   ->hit[i]   : hits
 *                   ->par      : these parameters are overwritten
 *                               by the callback-routine "get_par_from_gpar"
 *                               but they may contain fixed parameters
 *                               which are then used as is
 *                   ->any other fields evaluated by the callback-routine
 *                     (e.g. the vertex-z-parameter)
 *                 - start parameter for the generalized parameters
 *                 - deltas for that parameters used for finite-difference
 *                  computation
 *                 - a pointer to a callback routine that computes
 *                  the track-parameters from the generalized parameters
 *                 - the maximum number of iterations
 *
 *        OUTPUT:  ->gpar contains the fit-result for the generalized
 *                   parameters
 *                 
 *        RETURN VALUE: 0 if o.k, negative values are the following errors:
 *                 -2,-3,-4 : different types of singularities when
 *                            solving the linear equations
 *                 -5 : the parameters went outside the allowd range
 */

#include <math.h>

#include <cstring>   // for memset()

#include "TLattice.h"

using namespace std;

extern TLattice *Lat;
// extern int jter;     /* Working around CoRe NaN bug*/

/* this should be the most gereral fit we can do:
   - an arbitrary set of (generalized) parameters
   - an arbitrary number of tracks
   the above are connected by a callback-routine calculating the
   individual track-parameters from the generalized parameters
   
   this replaces the traditional single-track fit (Chi2Fit),
   the double-track fit (vertex_fit) and serves as a new type of
   kinematical 4-track fit
   
   now that we compute coordinate derivatives directly from the lattice,
   the drawback is that we can derive only with respect to the primary
   parameters (that make up the lattice). If the parameters for the fit
   are different (e.g. if applying an elastic constraint or if doing the
   fit for z!=0 ), we have to compute a matrix transforming the primary
   derivatives to the ones we are interested in. This matrix is calculated
   using finite differences and the callback-routine. Since usually
   this matrix contains many 0-entries, there is special optimization
   to skip the 0's so to limit the amount of cpu-time needed for
   the additional matrix multiplication
*/
     
typedef struct {
  float h[NROWS];            /* Coordinates */
  float sigma[NROWS];        /* Resolution */
  int index[NROWS];
  int n_activ;               /* Number of active coordinates */
  Parameters q;              /* Residuals (of the interpolation) [0..1] */
  int grid_point;            /* Grid-point */
  int old_grid_point;        /* Grid-point previous */
  int grid_checks;           /* Nummber of succesive grid checks */
} Flex_Hits;
    

#define MAX_NTR   2
#define MAX_NPAR  5*MAX_NTR

/*** static storage ***/

static  Flex_Hits hits00[MAX_NTR];
static  Hits      hits0[MAX_NTR];
static  float df_dx[MAX_NPAR][NROWS*MAX_NTR];
static  float dertrans[MAX_NTR][MAX_NPAR][5];
static  float a[MAX_NPAR][MAX_NPAR],b[MAX_NPAR];

static  int ntr, npar;

/*** function prototypes ***/

/* These 3 are in this file, for easier profiling they are
   single functions */

static void mult_dermat(int j,float *cder);
static void calc_amat_and_bvect(void);
static int solve_linear_equations(void);

extern void fcsf_get_par(Track *t,float *p);

int TLattice::FastFit(
  Track *t,                /* Pointer to tracks */
  float **sig,             /* Resolutions (dependent on fiber cluster width) */
  int ntrs,                /* Number of tracks */
  float *gpar,             /* Start- and return-values for parameters */
  float *dpar,             /* Deltas for finite differences */
  int npars,               /* Number of parameters */
  //  void get_par_from_gpar(Track*,float*), /* Function to calculate the track-parameters */
  int max_iter )           /* Maximum number of iterations */
{
  Track tp[MAX_NTR],tm[MAX_NTR],t0[MAX_NTR];

  float gp2[MAX_NPAR];
  
  int iter = 0;
  int i, j, itr;

  int ret;

  ntr  = ntrs;     /* assign to 'static' variables */
  npar = npars;

  /********** FOR ALL TRACKS, FILL THE hits00 STRUCTURES **********/

  /* Fill a structure containing only the active coordinates,
     the lattice-index-list, the sigma's for these coordinates,
     the number of active coordinates, the hit-pattern */

  for (itr = 0; itr<ntr; itr++) {
    int k, l, m;
    unsigned int hp;
     
    for (i=k=l = 0, m=32; l<NRW32; l++, m += 32) {
      hp = t->hits.hp[l];
      for (j = 1; i<m; i++, j <<= 1) if (hp&j) {
	hits00[itr].h[k] = t[itr].hits.h[i];
	hits00[itr].index[k] = i;
	hits00[itr].sigma[k] = sig[itr][i];
	k++;
      }
    }
    hp = t->hits.hp[NRW32];
    for (j = 1; i<NROWS; i++, j <<= 1)    if (hp&j) {
      hits00[itr].h[k] = t[itr].hits.h[i];
      hits00[itr].index[k] = i;
      hits00[itr].sigma[k] = sig[itr][i];
      k++;
    }             
    hits00[itr].n_activ = k;
    hits00[itr].grid_checks = 0;
  }

  for(i=0;i<ntr;i++) tp[i] = tm[i] = t0[i] = t[i];
  
  while(1) { /* this loop is left via 'break' from inside */

    /********** LOOP OVER MANY ITERATIONS **********/


    /* first of all, compute the current grid-points for all tracks
       to be able to check for convergence as soon as possible */
  
    //    get_par_from_gpar(t,gpar);
    fcsf_get_par(t,gpar);

    for (itr = 0; itr<ntr; itr++) {
      float* pp = &((t+itr)->par.x);

      /* check parameters against allowed values */
      for (i = 0; i<5; i++) if (fabs(pp[i])>max_par[i]) return -5;

      global_grid_point = hits00[itr].grid_point =
	Lat->get_point_par(&((t+itr)->par),&hits00[itr].q);
                 
      /* (the calling routine relies on 'global_grid_point' in case
          of ntr=1, so we must assign it here ) */
    }    

    /* restore the parameter-fields */
    for (itr = 0; itr<ntr; itr++) t[itr].par = t0[itr].par;
  
    /* the fit is converged if all grid-points are the same as in
       the previous iteration */

    if (iter) {
      if (iter>=max_iter) break;

      for (itr = 0; itr<ntr; itr++ ) {
	int *grid_checks = &hits00[itr].grid_checks;
	if (hits00[itr].old_grid_point==hits00[itr].grid_point) (*grid_checks)++;
	else                                                (*grid_checks) = 0;
	if (*grid_checks<2) break;
      }
      if (itr==ntr) break;  /* this break will leave the for(iter..) loop */
    }

    for (itr = 0; itr<ntr; itr++) hits00[itr].old_grid_point = hits00[itr].grid_point;

    /* we must fill the matrix of derivatives
       first compute the transformation matrix from finite differences */
       
    for (i = 0; i<npar; i++) gp2[i] = gpar[i];

    /* reset the relevant fields of the matrix */

    for (itr = 0; itr<ntr; itr++)
      memset(&(dertrans[itr][0][0]),0,sizeof(int)*npar*5);

    for (i = 0; i<npar; i++) {
      /* restore the parameter-fields */
      for (itr = 0; itr<ntr; itr++) tp[itr].par = tm[itr].par = t0[itr].par;
      
      gp2[i] = gpar[i] + dpar[i];
      //      get_par_from_gpar( tp , gp2 );
      fcsf_get_par(tp,gp2);
      gp2[i] = gpar[i] - dpar[i];
      //      get_par_from_gpar( tm , gp2 );
      fcsf_get_par(tm,gp2);
      gp2[i] = gpar[i];


      for (itr = 0; itr<ntr; itr++)
      {
        register float *tpp = &(tp[itr].par.x);
        register float *tmp = &(tm[itr].par.x);
        float d;

        for (j = 0; j<5; j++)
          if ((d=(*tpp++)-(*tmp++))!=0.)
            dertrans[itr][i][j] = d / dpar[i];
      }
    }
          
    /* now we must calculate the derivatives for all tracks */
    
    for (itr = 0; itr<ntr; itr++) {
      float cder[5*NROWS];

      /* pass the grid-point via a global variable :-( */

      global_grid_point = hits00[itr].grid_point;
    
      Lat->derive_lattice(hits0[itr],cder,
			 hits00[itr].index,hits00[itr].n_activ,
			 &(hits00[itr].q) );
      // if (jter==0) return -6;     /* Working around CoRe NaN bug*/
      /* to get the dfdx, these derivatives must by multiplied
	 with the appropriate row of the transformation matrix */

      mult_dermat(itr,cder);
    }

    /* from this dfdx and the coordinates calculate the a[][] matrix
       and b[] vector that are used to solve the system of linear
       equations to get the parameter changes */

    calc_amat_and_bvect();

    /* solve the linear equations */

    if ( (ret = solve_linear_equations()) != 0 ) return ret;

    /* apply parameter change */

    for(i=0;i<npar;i++) gpar[i] += b[i];

    iter++;
  } /* while(1) (iteration loop) */

  return 0;
}


/* this is a matrix-multiplication which transforms the
   derivatives with respect to the TRACK-parameters
   into the derivatives with respect to the generic
   parameters
   
   the matrix 'dertrans' is now TWICE what it should be
   for speed reasons
*/

static void mult_dermat(int itr,float *cder)
{
  int i,k;
  int jNROWS = NROWS*itr;
  int n = hits00[itr].n_activ;

  for (i = 0; i<npar; i++) {
    int   l = 0;
    float c;

    for (k = 0; k<5; k++) {
      if ((c = dertrans[itr][i][k])) {
	register float *cderp = cder+k*NROWS;
	register float *dfp = &(df_dx[i][jNROWS]);
	if (l)  for (l = 0; l<n; l++) *dfp++ += *cderp++ * c;
	else    for (     ; l<n; l++) *dfp++  = *cderp++ * c;
      }
    }
    if (l==0) { /* just in case there was no entry, delete it */
      register float *dfp = &(df_dx[i][jNROWS]);
      for (; l<n; l++) *dfp++  = 0.;   
      if (ntr==1)
	cout<<"**TLattice::FastFit: Nothing depends on parameter "<<i<<endl;
    }
  }
}
   

static void calc_amat_and_bvect()
{
  int i, j, l, k;

  float hit_dev[NROWS*MAX_NTR];

  /* pre-calculate twice the normalized hit-deviation */

  for (l = 0; l<ntr; l++) {
    int n = hits00[l].n_activ;
    for (k = 0; k<n; k++)
      hit_dev[NROWS*l+k] =
	(hits00[l].h[k]-hits0[l][k])/hits00[l].sigma[k];
  }

  for (i = 0; i<npar; i++) {
    float bb = 0.;

    for (l = 0; l<ntr; l++) {
      register float *sip =  &(hits00[l].sigma[0]);
      register float *dfp =  &(df_dx[i][NROWS*l]);
      register float *hdp =  &(hit_dev[NROWS*l]);

      for (k = hits00[l].n_activ; k; k--) {
	*dfp /= *sip++;
	bb += (*hdp++) * (*dfp++); 
      }
    }
 
    b[i] = bb;
      
    for(j = 0; j<=i; j++) {
      float d2=0.;

      for(l = 0; l<ntr; l++) {
	register float *dfip = &(df_dx[i][NROWS*l]);
	register float *dfjp = &(df_dx[j][NROWS*l]);

	for (k = hits00[l].n_activ; k; k--) d2 += (*dfip++) * (*dfjp++);
      }
      a[i][j] = a[j][i] = d2/2.;
    }
  }
}

static int solve_linear_equations()
{
  int i, j, k;

  /* solve the linear equations */

  for (i = 0; i<npar; i++)
    for (j = i+1; j<npar; j++) {
      float f;

      if (a[i][i]==0.) return -2;

      f=a[i][j]/a[i][i];

      for(k=i+1;k<npar;k++) a[k][j] -= f*a[k][i];
      b[j] -= f*b[i];
    }

  for (i = npar-1; i>0; i--)
    for (j =i-1; j>=0; j--) {
      if (a[i][i]==0.) return -3;
      b[j] -= a[i][j]/a[i][i]*b[i];
    }

  for(i = 0; i<npar; i++) {
    if (a[i][i]==0.) return -4;
    b[i] /= a[i][i];
  }
  return 0;   
}
/*
 *   derive_lattice: compute the coordinates and all 5 derivatives
 *                       for a given set of parameters
 *
 *         INPUT:  - the parameters in preprocessed form via the
 *                   grid-point (in a global variable) and the
 *                   normalized parameter-residuals
 *                 - the number of activ coordinates
 *                 - a vector of lattice-indices
 *                 - hit-pattern masks for left- and right detectors
 *                
 *        OUTPUT:  - the coordinates in 'coords'
 *                 - the derivatives in 'cder'
 */

void TLattice::derive_lattice(
  float *coords,      /* where to write the coordinates */
  float *cder,        /* where to write the derivatives */
  int *index,         /* lattice-indices (0-11 , +13 for left side) */
  int n_activ,        /* number of active coordinates */
  Parameters *qp)     /* parameter residuals (0..1), the grid-point is
                         passed via 'global_grid_point' */
{
  float c_l0, c_l1, c_l2, c_l3, c_l4;
  float c_r0, c_r1, c_r2, c_r3, c_r4;
  float c_0, c_1, c_2, c_3;

  float cp4[4*NROWS], xaver[16*NROWS], xdiff[16*NROWS], *ap, *dp;

  Hits  *glp, *grp;
  int   i, j;
  
  /* The first step is to compute the x-averaged values for the
     2**4 = 16 vectors of NROWS numbers = 16*NROWS values
     We do this for the means and the differences in a single
     step, that is we compute 32*NROWS values
     
     since this is probably the most time consuming step,
     code it carefully !

     Then for the 2 parameters (y,pxzi), compute the average for the
     other 2 parameters (tx,ty)

     Then for the 2 parameters (tx,ty), compute the average for the
     other 2 parameters (y,pxzi)
  */

  /**************** CALCULATE THE x-AVERAGED COORDS ***************/

  c_l0 = 1. - (c_r0 = qp->x);
  c_l1 = 1. - (c_r1 = qp->y);
  c_l2 = 1. - (c_r2 = qp->pxzi);
  c_l3 = 1. - (c_r3 = qp->tx);
  c_l4 = 1. - (c_r4 = qp->ty);

  /********** LOOP OVER 16 GRID-POINTS **********/

  for(i = 0; i<16; i++) {
    float *flp,*frp;

    glp = grid + global_grid_point + multi_jump[2*i];
    grp = glp + jump[0];

    ap = xaver + i*NROWS; dp = xdiff + i*NROWS;
      
    /***** LOOP OVER ONLY ACTIVE COORDS *****/

    for (j = 0; j<n_activ; j++) {
      register float l,r;

      *ap++ = c_l0*(l = (*glp)[index[j]])
	+     c_r0*(r = (*grp)[index[j]]);
      *dp++ = r-l;
    }
  } /* loop over 16 grid-points */

  // if (jter==0) return;     /* Working around CoRe NaN bug*/

  /*************** FOR (y,pxzi) COMPUTE THE AVERAGE FOR (tx,ty) ***************/

  c_0 = c_l3*c_l4; c_1 = c_r3*c_l4; c_2 = c_l3*c_r4; c_3 = c_r3*c_r4;

  for (i = 0; i<4*NROWS; i += NROWS) {
    int je = i+n_activ;
    for (j = i; j<je; j++)
      cp4[j] = c_0*xaver[j        ] + c_1*xaver[j+4*NROWS]
	+      c_2*xaver[j+8*NROWS] + c_3*xaver[j+12*NROWS];
  }

  /*************** (y,pxzi) CENTER VALUES and DERIVATIVES ***************/

  c_0 = c_l1*c_l2; c_1 = c_r1*c_l2; c_2 = c_l1*c_r2; c_3 = c_r1*c_r2;
    
  for (i = 0; i<n_activ; i++) {
    coords[i] =
      c_0*cp4[i] + c_1*cp4[i+NROWS] + c_2*cp4[i+2*NROWS] + c_3*cp4[i+3*NROWS];
    cder[i+  NROWS] = 
      (c_l2*(cp4[i+  NROWS]-cp4[i])+c_r2*(cp4[i+3*NROWS]-cp4[i+2*NROWS]))
      *il_diff.y;
    cder[i+2*NROWS] =
      (c_l1*(cp4[i+2*NROWS]-cp4[i])+c_r1*(cp4[i+3*NROWS]-cp4[i+  NROWS]))
      *il_diff.pxzi;
  }

  /*************** FOR (tx,ty) COMPUTE THE AVERAGE FOR (y,pxzi) ***************/

  for (i = 0; i<16*NROWS; i += 4*NROWS) {
    int je = i+n_activ;
    float *cp4p = cp4 + (i>>2);
      
    for (j = i; j<je; j++)
      *cp4p++ = c_0*xaver[j        ] + c_1*xaver[j+NROWS]
	+       c_2*xaver[j+2*NROWS] + c_3*xaver[j+3*NROWS];
  }

  /**************** "tx" and "ty" DERIVATIVES ***************/
    
  for (i = 0; i<n_activ; i++) {
    cder[i+3*NROWS] = 
      (c_l4*(cp4[i+  NROWS]-cp4[i])+c_r4*(cp4[i+3*NROWS]-cp4[i+2*NROWS]))
      *il_diff.tx;
    cder[i+4*NROWS] =
      (c_l3*(cp4[i+2*NROWS]-cp4[i])+c_r3*(cp4[i+3*NROWS]-cp4[i+  NROWS]))
      *il_diff.ty;
  }

  /* up to here it was very efficient, but the x-derivative is still
     missing and this is cpu-extensive :-( */

  c_0 = c_l3*c_l4; c_1 = c_r3*c_l4; c_2 = c_l3*c_r4; c_3 = c_r3*c_r4;

  for (i = 0; i<4*NROWS; i+= NROWS) {
    int je = i+n_activ;
    for (j = i; j<je; j++)
      cp4[j] = c_0*xdiff[j        ] + c_1*xdiff[j+4*NROWS]
	+      c_2*xdiff[j+8*NROWS] + c_3*xdiff[j+12*NROWS];
    }

  /*************** COMPUTE THE "x" DERIVATIVE ***************/
    
  c_0 = c_l1*c_l2; c_1 = c_r1*c_l2; c_2 = c_l1*c_r2; c_3 = c_r1*c_r2;

  for (i = 0; i<n_activ; i++)
    cder[i] =
      (c_0*cp4[i]+c_1*cp4[i+NROWS]+c_2*cp4[i+2*NROWS]+c_3*cp4[i+3*NROWS])
      *il_diff.x;
}
