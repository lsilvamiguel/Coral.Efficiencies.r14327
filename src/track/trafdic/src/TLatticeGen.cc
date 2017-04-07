// $Id: TLatticeGen.cc,v 1.9 2010/02/15 21:07:50 ybedfer Exp $

/*!
  \file    TlatticeGen.cc
  \brief   Utils for the generation of the Lattice file == Dico
  \author  Y.B
  \version $Revision: 1.9 $
  \date    $Date: 2010/02/15 21:07:50 $
  */

#include <cmath>
#include <stdio.h>
#include "TSetup.h"
#include "TOpt.h"
#include "TLattice.h"

// TLatticeGen.cc
//  Generate Dico == tracks lattice for fitting:
//  - Generate single-track events == loop on track parameters t in 5-dimension
//  phase space lattice
//  - To be tracked through COMPASS
//  - Write hits(t) to ouput file
//
//  Lattice is defined @ z==0 extrapolated to z==700mm w/ straight line
//  z=700mm corresponds to the (assumed) end of the uniform solenoid field. The
//  choice of this limit is meant to provide for:
//  - avoiding the intricate dependence of the hit positions upon the
//  track paramters that would result from the tracking through several
//  Txm of longitudinal field,
//  - yet easy back-tracking (whether this is indeed needed is not yet
//  clear) using helix trajectories in uniform field (implictly
//  assuming SM1 field is negligible in the target region).
//  The extrapolation to z==0 allows to limit the number of tracks pointing 
//  away from the target while retaining a cubic phase space volume for
//  lattice.

// TLatticeGen.cc: Lattice generator
// void GenInit()                    : Init
// void GenServer(float*,float*,int*): Track server
// void GenOutput(float*)            : Output to file
// void GenClose()                   : Close

// $Log: TLatticeGen.cc,v $
// Revision 1.9  2010/02/15 21:07:50  ybedfer
//  "dZStart", i.e. the distance from the origin of the dico (set @ the
// target center) to its starting point (where extrapolation starts and where
// track parameters are expressed) is now a data member. It's set = 70 cm
// when the target field is active, 0 cm otherwise.
//
// Revision 1.8  2009/03/07 20:40:55  ybedfer
//  The reference plane where the extrapolation starts in the case of an
// active target field is positioned w.r.t. the target center. => In the
// standard case where the target is (approx.) centered @ 0: no change. Only
// change is for those setups where the target magnet is moved away from 0,
// e.g. Drell-Yan.
//
// Revision 1.7  2006/12/13 05:01:25  ybedfer
//  Write newly introduced member array "sieve[6]" (listing which of the
// ``Dico'able'' detectors take part in any given Dico instance) to the Dico
// file.
//
// Revision 1.6  2006/11/10 02:47:20  ybedfer
//  Use "TSetup.TargetCenter[0]" instead of "TOpt::Target[0]". The former is
// now properly defined in all cases, cf. "TSetup::Init".
//
// Revision 1.5  2005/04/10 16:01:37  ybedfer
//  For the muon case: track parameters are specified 70 cm upstream of the
// point where extrapolation starts and transfered to there by straight
// track extrapolation.
//
// Revision 1.4  2004/04/23 16:56:01  ybedfer
//  Switch muons/hadrons based on "MagScale[0]==0". If the latter: no
// straight line extrapolation to ref plane (which is meant to bypass
// target's field). Origin is taken from "TraF Target" (implemented only
// for hadrons (for backward compatibility, would need a versioning of
// the Dicos to extend this to muons)).
//
// Revision 1.3  2002/08/23 13:14:26  compass
// include <cmath> for RH72 compatibility
//
// Revision 1.2  2001/04/06 14:26:49  ybedfer
//  TraFDic,v2.0
//
// Revision 1.1  2000/06/25 18:31:26  ybedfer
// Initial version.
//

static float    step[5];
static FILE  *Dico_File;
static unsigned int full_last_hp;
static float ZStart;

// ******************** INIT THE GENRATION OF THE Dico ********************
void TLattice::GenInit() {

  int i;

  for (i = 4; i>=0; i--)
    step[i] = bound[i].num>1 ? 
      (bound[i].max-bound[i].min)/(bound[i].num-1) : 0;

  //       ***** OPEN OUTPUT FILE
#define DICO_FILE "dicofit.out"
  if ((Dico_File = fopen(DICO_FILE,"w"))==NULL) {
    printf("** TLattice::GenInit:\a Error Opening Output file \"%s\"!\n",
	   DICO_FILE); fflush(stdout); perror(" <="); exit(1);
  }

  //        ***** WRITE BIT PATTERN SIEVE
  fwrite(sieve,sizeof(sieve),1,Dico_File);

  //        ***** WRITE LATTICE BOUNDS
  for (i = 0; i<5; i++)
    fwrite(&bound[i],sizeof(parbound),1,Dico_File);

  // ***** HIT PATTERN

  for (i = 0, full_last_hp = 0; i<NROWS%32; i++)
    full_last_hp |= (1<<i);

  // ***** INIT GRID POINT

  global_grid_point = 0;

  // ***** DERIVE ORIGIN Z and STARTING Z from SETUP
  const TSetup &setup = TSetup::Ref();
  if (setup.MagScale[0]==0)
    dZStart = 0;
  else
    // If the target field is active, we avoid incorporating into the dico the
    // extra intricacy it introduces in the dependence of the track upon the
    // track parameters @ the origin. Therefore, if we still set the origin
    // at the center of the target, we move in straight line from this origin to
    // a plane situated "dZStart" away from it and only then start extrapolating
    // in the magnetic field.
    dZStart = 70;
  ZStart = setup.TargetCenter[0]+dZStart;
}

// ******************** TRACK SERVER ********************
void TLattice::GenServer(float *V,float *P,int *ipart)
{
  // ***** GRID POINT -> Track; INCREMENT GRID POINT

  int i, jdx; Parameters par; float *flp;
  for (i = 4, jdx = global_grid_point++, flp = &par.ty; i>=0; i--, flp--) {
    int num = bound[i].num, kdx = jdx%num;
    *flp  = bound[i].min+kdx*step[i];
    jdx /= num;
  }

  // ***** TKINE TRACK => Convert Track -> TKine
    
  if (fabs(par.pxzi)<0.001) par.pxzi = 0.001;

  // The momentum dimension of the dico is 1/pxz == 1/horiz component of p
  float p = sqrt(1+par.ty*par.ty/(1+par.tx*par.tx))/fabs(par.pxzi);

  float pz = fabs(p)/sqrt(1+par.tx*par.tx+par.ty*par.ty);
  P[0] = pz; P[1] = par.tx*pz; P[2] = par.ty*pz; 
	
  if (par.pxzi>0) *ipart = 8;  // pi+
  else            *ipart = 9;  // pi-
    
  // ***** Fill VERTEX container
  // Extrapolate from point where track parameters are specified to
  // starting point of extrapolation.
  V[1] = par.x+dZStart*par.tx;
  V[2] = par.y+dZStart*par.ty;
  V[0] = ZStart;
}

// ******************** OUTPUT TO FILE ********************
void TLattice::GenOutput(float *coords)
{
  fwrite(coords,NROWS*sizeof(float),1,Dico_File);
}

// ******************** CLOSE FILE ********************
void TLattice::GenClose()
{
  fclose(Dico_File);
}
