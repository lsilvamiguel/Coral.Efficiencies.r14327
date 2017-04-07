// $Id: TDisplayDetInfo.cc 13148 2011-12-28 16:55:25Z kbicker $

/*!
  Function (called from TDisplay::DrawMenu) to pickup detector plane and print detector's info. 
*/

#include <iostream>
#include <stdio.h>
#include "CsDetector.h"
#include "TDisplay.h"
#include "TSetup.h"
#include "TPlane.h"
#include "TDetect.h"

using namespace std;

/*
  Version specific to "lattice" alternative.
  Changes w/ respect to "traffic/sources/TDisplayDetInfo.cc".
  i) More info.
*/

void TDisplay::DetInfo(int ip)
{
  const TSetup &setup = TSetup::Ref();
  if (ip<0) {
    cout<<"TDisplay::DetInfo ==> No plane was picked up\n"; return;
  }
  const TPlane  &p = setup.vPlane(ip); 
  const TDetect &d = setup.vDetect(p.IDetRef);  // ref. to detector's vector in Geom
  if (ip!=d.IPlane) {
    cout<<"TDisplay::DetInfo ==> detector <-> plane referencies mismatch\n";
    return;
  }

  printf("Detector plane \"%s\" #%-3d   detID = %4d",
	 d.Name.c_str(),p.IPlane,d.IDet);
  printf("   Center = (%9.4f,%9.4f,%8.2f) cm",d.X(1),d.X(2),d.X(0));
  printf("   #Wires = %4u   Pitch, 1st wire = %7.4f,%9.4f cm\n", 
	 d.Nwires, d.Pitch, d.Uorig);
  printf(" %%X0 = %6.3f   Proj. = #%d,%5.1f deg.  Typical resol. = %5.0f um",
	 100*d.Siz(0)/d.RadLen,p.IProj,setup.vProj()[p.IProj]*.1,d.Resol*10000);
  if (d.TResol>=0) printf(", %.2f ns\n",d.TResol);
  else             printf("\n");
  if (d.DZtype!=NO) {
    printf(" DZ w.r.t. center = (%9.4f,%9.4f) cm   ",d.DZydrs,d.DZzdrs);
    if (d.DZtype==CIRCULAR) printf("Radius = %5.2f",sqrt(d.DZydim));
				      
    else                    printf("Size   = (%5.2f,%5.2f)",d.DZydim,d.DZydim);
    if (d.EZtype==NO) printf(" cm   Empty\n");
    else              printf(" cm\n");
  }
  printf(" #Hits = %3zu\n", p.vHitRef().size());
  if (!p.IFlag) printf("N.B.: \"%s\" is turned OFF\n",d.Name.c_str());
  cout<<endl;
}
