// $Id: TLatticeFitInterface.cc,v 1.2 2006/03/13 14:50:11 ybedfer Exp $
/*!
  \file   TLatticeFitInterface.cc
  \brief  Interface to the fast fit method
  \author Yann Bedfer
  */

#include <math.h>
#include "TLattice.h"

extern TLattice *Lat;
// extern int jter;     /* Working around CoRe NaN bug*/
extern float compute_chi_square(Track_Hits *,Track_Hits *);
extern void fcsf_get_par(Track *t,float *p);

/* this is the 1 track fit

   There a some features we have to implement to guarantee compatibilty:
   
   - allow for fixed parameters
   - track back to z=0. if t->z != 0.

   In any case, make sure that on return all coordinates are computed
   correctly

   The requirement for compatibility with the old code isn't good
   for speed-enhancement:
   - a large fraction (30%) of cpu-time accounts for the calculation
     of the coordinates for the fitted parameters. In most cases this
     is not needed
   Therefore 'fast_chi_square_fit' uses a flag to select whether this
   calculation is done. For further speed-up, the tracking code could
   call 'fast_chi_square_fit' directly with mode=0

*/

static float delta[5] = { .25,.5,0.001,0.0001,0.0002 };  /* to give ?mm */


int TLattice::FastChi2Fit(
             Track *t,   /* pointer to track */        
             int mode    /* 1=compute all coords and chi2, 0=dont */)
{

  int    i, np, rval;
  float  paras[5], dparas[5];
  Track_Hits track_hits, track_hits0;
  Parameters return_pars;
  float *sig[1];
   
  track_hits0 = t->hits;

  /***** COUNT THE NUMBER OF PARAMETERS AND FILL THE VECTOR *****/
     
  for (np=i = 0; i<5; i++)
    if (t->par_pattern & (1<<i)) {
      paras[np] = (&(t->par.x))[i];
      dparas[np++] = delta[i];
    }

  /***** CALL THE FIT ROUTINE *****/
  sig[0] = Lat->sigmas;         /* ... fixed weights */
  //  rval = Lat->FastFit(t,sig,1,paras,dparas,np,fcsf_get_par,10);
  rval = Lat->FastFit(t,sig,1,paras,dparas,np,10);
  // if (jter==0) return -6;     /* Working around CoRe NaN bug*/

  if (rval <0) return rval;

  /* write the generalized parameters to the track
     this is for compatibility only: there are already
     the TRACK-parameters. The main difference is that
     for z!=0. we now write the x/y at z!=0 instead of
     x/y at 0.
  */

  return_pars = t->par;    /***** SAVE PARS @ POSSIBLY z!=0 IN "t" *****/

  fcsf_get_par(t,paras);   /***** PARS @ z=0 *****/

  for (np=i=0;i<5;i++)
    if (t->par_pattern & (1<<i))  (&(return_pars.x))[i] = paras[np++];

  if (mode) {

    /***** GET THE EXTRAPOLATED COORDINATES AND THE chi2 *****/

    global_grid_point = Lat->get_point(&t->par);
    Lat->get_coords(t);
    t->chi2 = compute_chi_square(&track_hits0,&t->hits);
  }

  t->par = return_pars; 
  return 0;
}

/********** CALLBACK ROUTINE FOR "fast_chi_square_fit" **********/

void fcsf_get_par(Track *t,float *p)
{
  int i,np;

  /***** EXPAND THE PARAMETER VECTOR TO THE TRACK PARAMETERS *****/

  for (np=i = 0; i<5; i++)
    if (t->par_pattern & (1<<i))  (&(t->par.x))[i] = p[np++];

  /* IF THE Z-COORDINATE IS NOT 0., USE THE TRACK-BACK ROUTINE TO MAKE IT 0 */
  /* not done for the present
  if (t->z !=0. ) {
    float z_keep = t->z;
    
    trafo_track_to_z0(t);
    t->z = z_keep;
  } */
}
