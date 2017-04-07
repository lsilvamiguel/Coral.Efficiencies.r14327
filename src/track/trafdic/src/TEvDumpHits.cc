// $Id: TEvDumpHits.cc 13148 2011-12-28 16:55:25Z kbicker $

#include <iostream>
#include <cstdio>
#include <CsSTD.h>
#include "Coral.h"
#include "TSetup.h"
#include "TEv.h"

  /*!
    \file    TEvDumpHits.cc
    \brief   Dump for debug purpose
    \author  Y.B
  */

using namespace std;

void  TEv::DumpHits()
{
  //********** UNUSED HITS

  const TSetup& setup = TSetup::Ref();

  for (int imc=0; imc<int(vecKine.size()); imc++) {
    printf("Kine %d\n",imc);
    int nh, igr = -2;
    vector<THit>::iterator ih;
    for (ih = vecHit.begin(), nh = 0; ih!=vecHit.end(); ih++) {
      int ikin = ih->IKine, ikin2 = ih->IKin2, niks, ik;
      if (ikin<0) continue;
      if (ikin2>=0) niks = 2;   // A 2-track hit
      else          niks = 1;
      for (ik = 0; ik<niks; ik++) {
	if (ikin==imc &&
	    ih->sTrackID().empty()) {  // Not appropriate if 2-track hit
	  int igrp = igr, ig;
	  for (ig = 0, igr = -1; ig<(int)setup.vIplFirst().size(); ig++) {
	    if (setup.vIplFirst()[igr]<=ih->IPlane &&
		ih->IPlane<=setup.vIplLast()[igr]) { igr = ig; break; }
	  }
	  if (igrp!=-2 && igr!=igrp) nh = 0;
	  if (!(++nh%8)) printf("\n");
	  if (ih->IOrig)
	    printf("%3d %4d* ",ih->IPlane,ih->IHit);
	  else
	    printf("%3d %4d, ",ih->IPlane,ih->IHit);
	}
	ikin = ikin2;  // Try 2nd track if exists (i.e. niks==2)
      }
    }
    printf("\n");
  }
}













