// $Id: TEvDumpEvent.cc 13148 2011-12-28 16:55:25Z kbicker $

#include <iostream>
#include <cstdio>
#include <CsSTD.h>
#include "Coral.h"
#include "TEv.h"
#include "TKine.h"
#include "THit.h"
#include "TOpt.h"
#include "CsEvent.h"
#include "CsMCUtils.h"

using namespace std;

  /*!
    \file    TEvDumpEvent.cc
    \brief   Dump for debug purpose
    \author  Y.B
  */

void  TEv::DumpEvent()
{
  if (!CsEvent::Instance()->isAMonteCarloEvent() ) return;

  CsMCUtils utils;
  const TSetup& setup = TSetup::Ref();
  //********** LOOP OVER MC TRACKS

  for (int imc=0; imc<int(vecKine.size()); imc++){
    TKine *mc   = &vecKine[imc];
    int   Q = mc->Q();
    if (Q) {
      float mom =
	sqrt(mc->P(0)*mc->P(0)+mc->P(1)*mc->P(1)+mc->P(2)*mc->P(2))/Q;
      printf("%d %7.2f %6.2f %6.2f\n",imc,mom,mc->P(1),mc->P(2));

      //***** HITS LIST

      int nh, shower_track, jg;
      vector<int>::const_iterator ih;
      for (ih = mc->vHitRef().begin(), nh=shower_track=jg = 0;
	   ih!= mc->vHitRef().end(); ih++) { // loop over MC track hits
	THit& h = vecHit[*ih];
	if (h.IOrig!= 0) {                    // skip hits from showers...
	  if (ih==mc->vHitRef().begin()) shower_track = 1;
	  else if (!shower_track) continue;  // ...but not shower tracks
        }

	// Display hit in zone blocs, numbering them zone by zone
	int ng = int(setup.vIplFirst().size()); ng = ng > 3 ? 3 : ng;
	for (int ig = 0; ig<ng; ig++)                 // loop over det. groups
	  if (h.IPlane>=setup.vIplFirst()[ig] &&
	      h.IPlane<=setup.vIplLast()[ig]) {
	    if (ig!=jg) { nh = 0; printf("\n"); }
	    jg = ig;
	  }
	
	printf(" %3d %3d %7.2f",h.IPlane,h.iHit,h.u);
	if (!(++nh%5)) printf("\n");
      }
      if (nh) printf("\n");
      if (shower_track) printf("S\n");
      if (mc->isPileup()) printf("A\n");
    }  // end of charged MC tracks
  } // end of loop over MC tracks

}













