//   $Id: TTrackDumpHits.cc 13148 2011-12-28 16:55:25Z kbicker $

//   $Log: TTrackDumpHits.cc,v $
//   Revision 1.15  2009/12/16 21:35:40  ybedfer
//    Print the v-coordinate of the pixel TPlanes of MPs.
//
//   Revision 1.14  2008/08/21 15:53:38  ybedfer
//    - pixelGEMs: print also v coordinate.
//
//   Revision 1.13  2007/01/15 04:30:34  ybedfer
//    Print momenta, first and last, if available.
//
//   Revision 1.12  2007/01/08 01:48:40  ybedfer
//    PixelGEMs: Reduced chi2 derived from "TTRack::NDFs" instead of "NHits".
//
//   Revision 1.11  2006/12/13 05:11:08  ybedfer
//    - Get mirror info from THit instead of CsCluster.
//    - No longer redefine THit' TKine.
//
//   Revision 1.10  2005/05/23 14:58:02  ybedfer
//    - Changed format.
//    - Print KF chi2, if (supposed to be) available.
//
//   Revision 1.9  2005/05/16 14:50:46  ybedfer
//    Bug fix: consider case where missing planes have been inserted.
//
//   Revision 1.8  2005/04/08 22:18:29  ybedfer
//    Print QN chi2 if iFit&0x8 (and ==0x8).
//
//   Revision 1.7  2004/01/13 01:17:31  ybedfer
//    Bug fix: for when ``inserted planes'' have been added: have to check that
//   "ihp" is >= 0.
//
//   Revision 1.6  2003/03/11 18:34:07  ybedfer
//    QNewton chi2 is now taken from "Chi2aux" only if it is != 0 (which is
//   not the case @ PrePattern time).
//
//   Revision 1.5  2003/02/22 02:29:00  ybedfer
//    Print both KF and Dico chi2's (and bug fix for dico case).
//
//   Revision 1.4  2002/04/22 23:36:40  ybedfer
//    Accessor to TLattice via pointer Ptr().
//
//   Revision 1.3  2002/01/19 21:31:14  ybedfer
//    Determine whether hit is genuine only if event is MC.
//
//   Revision 1.2  2001/04/25 01:09:37  ybedfer
//   TraFDic,v2.1
//
//   Revision 1.1  2001/04/06 14:26:53  ybedfer
//    TraFDic,v2.0
//
//   Revision 1.2  2001/02/10 04:11:35  ybedfer
//   Modifications for:
//    - improved efficiency,
//    - interface to vertex,
//    - setup 2001.
//
//   Revision 1.1  2000/10/23 18:36:20  ybedfer
//   lattice,1.1
//
//----------------------------------------------------------------
//

/*!
  Method for debug purposes: dump hits to screen
  \param mode !=0 <-> Print residuals
*/

#include <stdio.h>
#include "TSetup.h"
#include "TEv.h"
#include "TTrack.h"
#include "THit.h"
#include "TLattice.h"

using namespace std;

void TTrack::DumpHits(int mode = 0)
{
  const TSetup &setup = TSetup::Ref();
  TEv &ev = TEv::Ref();

  printf("TTrack.Id %d\n",Id);

  // ***** PRINT A 1ST LISTING: Plane_#  Hit_#  Hit_Kine                *****
  // *****   OR if NOT GENUINE:                                         *****
  // *****                      Plane_#  Hit_#  Hit_Kine.LRProb ...     *****
  // *****                  ...          U      U_counterpart           *****
  // ***** (BEFOREHAND RECAST as GENUINE a MIRROR CLOSE TO COUNTERPART) *****

  vector<THit>& vHit = const_cast<vector<THit>&>(ev.vHit());
  list<int>::iterator ih; int nh, igr, igrp, nht = 0;
  for (int iter = 0; iter<((mode&0x4) ? 2 : 1); iter++) {
    nh=igr = 0; ih = lHitPat.begin(); while(ih != lHitPat.end()) {
      if (*ih<0) { ih++; continue; } 

      THit &h = vHit[*ih];  // *************** LOOP OVER HITS ***************

      igrp = igr;                                 // ***** NEW ZONE? => NEW LINE
      while (igr<int(setup.vIplFirst().size())) {
	if (setup.vIplFirst()[igr]<=h.IPlane && h.IPlane<=setup.vIplLast()[igr])
	  break;
	igr++;
      }
      if (igr!=igrp) { printf("\n"); nht = 0; }

      if (h.IKine<-1) {                           // ***** NON GENUINE MIRROR...
	printf(" %3d %4d %6.2f",h.IPlane,h.IHit,
	       -2-h.IKine+h.PtrClus()->getLRProb());
	if (!(++nht%5)) printf("\n");  // New line?
	if (h.Mirror>=0) printf(" %7.2f %7.2f",h.U,ev.vHit()[h.Mirror].U);
	else             printf(" ??????? ???????");
      }
      else {                                                    // ***** ...ELSE
	if (iter) printf(" %3d %4d %6.2f",h.IPlane,h.IHit,h.U);
	else      printf(" %3d %4d %6d",h.IPlane,h.IHit,h.IKine);
      }
      nh++;
      if (!(++nht%5)) printf("\n");  // New line?

      const TPlane &p = setup.vPlane(h.IPlane);  // ********** PIXEL **********
      const TDetect &d = setup.vDetect(p.IDetRef);
      if (d.IType==29 || d.IType==32) {
	if (iter) printf(" %3d %4d %6.2f",h.IPlane,h.IHit,h.V);
	else      printf(" %3d %4d %6d",h.IPlane,h.IHit,h.IKine);
	nh++;
	if (!(++nht%5)) printf("\n");  // New line?
      }
      ih++;
    }                                  // End of loop over hits
    if (nh) printf("\nNDFs %d\n",nh);
  }

  if ((mode&0x1) && !vGuests.empty()) {

    // Reference to "TLattice". In view of determining whether given plane
    // has meaningfull guestimates. This is not very satisfying to have a
    // refernce to "TLattice". Best would have been to have "TTrack" (which
    // could had guestimates determined by other mean than "TLattice")
    // return the answer. But for the time being...
    const TLattice *lat = TLattice::Ptr();    

    // ***** RESIDUALS
    // ***** PRINT A 2ND LISTING: Plane_#  Hit_#  Residual             *****
    // *****      IF NOT GENUINE: ...      Hit_#  Residual_counterpart *****

    int ipl0 = -1;
    for (igr = 0; igr<2; igr++)           // First plane where guestimate
      if (GuestsGroups&(1<<igr)) { ipl0 = setup.vIplFirst()[igr]; break; }
    if (ipl0<0) return;

    ih = lHitPat.begin() ; nht=nh=igr = 0; while (ih!=lHitPat.end()) {
      if (*ih<0) { ih++; continue; }

      const THit& h = vHit[*ih]; int jpl = h.IPlane-ipl0;

      if (jpl>=(int)vGuests.size()) break;
      // Guestimates end w/ QNewtonFit range
      if      (lat->dico_idx[h.IPlane]==-1) { ih++; continue; }
      else if (lat->dico_idx[h.IPlane]==-2) break;

      //        *************** LOOP OVER HITS w/in DICO ***************

      igrp = igr;                                 // ***** NEW ZONE? => NEW LINE
      while (igr<int(setup.vIplFirst().size())) {
	if (setup.vIplFirst()[igr]<=h.IPlane && h.IPlane<=setup.vIplLast()[igr])
	  break;
	igr++;
      }
      if (igr!=igrp) { printf("\n"); nht = 0; }

      printf(" %3d %4d %6.2f",h.IPlane,h.IHit,(h.U-vGuests[jpl])/h.SigU);

      if (h.IKine<-1) {                           // ***** NON GENUINE MIRROR...
	if (!(++nh%5)) printf("\n");
	if (h.Mirror>=0)
	  printf("     %4d %6.2f",h.IHit,vGuests[jpl]-ev.vHit()[h.Mirror].U);
	else
	  printf(" ???????? ??????");
      }
      if (!(++nh%5)) printf("\n");

      const TPlane &p = setup.vPlane(h.IPlane);  // ********** PIXEL **********
      const TDetect &d = setup.vDetect(p.IDetRef);
      if (d.IType==29 || d.IType==32) {
	printf(" %3d %4d %6.2f",h.IPlane,h.IHit,(h.V-vGuests[jpl-2])/h.SigV);
	if (!(++nht%5)) printf("\n");  // New line?
      }
      ih++;
    }                     // End of loop over planes
    if (nh) printf("\n");
  }

  // ***** TTrack Kine, CHI2

  FindKine(); printf("Kine = %d, 0x%02x",IKine,IFit);
  if (IFit&0x8) {
    if (NGroups()==1) {
      if (Chi2aux)    printf(" Dico Chi2 %f",Chi2aux/(NDics-4));
      else            printf(" Dico Chi2 %f",Chi2tot/(NDics-4));
    }
    else              printf(" Dico Chi2 %f",Chi2tot/(NDics-5));
  }
  if (IFit&0x66) {
    if (NGroups()==1) {
      if (NDFs>4)    printf(" KF Chi2 %f",Chi2tot/(NDFs-4));
    }
    else             printf(" KF Chi2 %f",Chi2tot/(NDFs-5));
  }
  if (Hfirst(5)) printf(" 1st=%.3f GeV",1/Hfirst(5));
  if (Hlast(5))  printf(" last=%.3f GeV",1/Hlast(5));
  printf("-----------------------------\n");
  return;
}

