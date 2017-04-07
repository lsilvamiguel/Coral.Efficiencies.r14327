// $Id: CsTrafficRefitting.cc 13148 2011-12-28 16:55:25Z kbicker $

/*!
   \file    CsTrafficReffting.cc
   \brief Coral interface to TraFFiC's tracks refitting.
   \author  Yann.Bedfer@cern.ch
*/

#include "CsErrLog.h"
#include "CsTrafficRefitting.h"
#include "Traffic.h"
#include "TEv.h"

//#define CsTrkRF_DEBUG
#ifdef CsTrkRF_DEBUG
#  include "CsGeom.h"
#  include "CsEvent.h"
#  include "CsTrack.h"
#endif

using namespace std;

//Constructor
CsTrafficRefitting::CsTrafficRefitting() {
  if (Traffic::Ptr()==NULL)
    CsErrLog::mes(elFatal,"\"Traffic\" object not instantiated.");
}

//Destructor
CsTrafficRefitting::~CsTrafficRefitting() {}

//Refitting method
bool CsTrafficRefitting::doRefitting(CsVertex *pVertex, list<CsTrack*> &tracksToBeDeleted) {

  if (!TEv::Ptr())  // ********** REQUIRE PREEXISTING TEv **********
    CsErrLog::mes(elFatal,"No pe-existing TEv (\"TraF ReMode[1]!=2\")");
  TEv &ev = TEv::Ref();

#ifndef CsTrkRF_DEBUG
  // "True" returned means actual refit took place
  return ev.TracksRefit(pVertex,tracksToBeDeleted);
#else
  double zTarget = CsGeom::Instance()->getTargetCenter();  // in mm
  list<CsTrack*> eTracks = CsEvent::Instance()->getTracks();
  list<CsTrack*>::iterator it;
  printf("==========\n");
  for (it = eTracks.begin(); it!=eTracks.end(); it++) {
    const CsTrack *t = *it; const vector<CsHelix> &helices = t->getHelices();
    // #parameters in track's fit is expected to be:
    // - Same for all helices.
    // - = 5 for a track w/ momentum, provided it's not a track reco'd in the
    //  sole beam telescope. (These ``beam'' tracks are here singled out by
    //  requiring their starting point to be upstream of target, which is not
    //  100% rigorous.)
    int nPars = (helices[0].getCop() && helices[0].getZ()>zTarget) ? 5 : 4;
    printf("chi2 %4.2f  (mm,GeV)",t->getChi2()/(t->getNDFs()-nPars));
    for (int ih = 0; ih<(int)helices.size(); ih++) {
      double cop = helices[ih].getCop();
      printf(" %6.0f,%6.1f ",helices[ih].getZ(),cop?1/cop:0);
    }
    printf("\n");
  }
  bool ret; if ((ret = ev.TracksRefit(pVertex,tracksToBeDeleted))) {
    printf("========== REFIT w/ event time = %.1f\n",ev.GetEventTime());
    for (it = eTracks.begin(); it!=eTracks.end(); it++) {
      const CsTrack *t = *it; const vector<CsHelix> &helices = t->getHelices();
      int nPars = (helices[0].getCop() && helices[0].getZ()>zTarget) ? 5 : 4;
      printf("chi2 %4.2f  (mm,GeV)",t->getChi2()/(t->getNDFs()-nPars));
      for (int ih = 0; ih<(int)helices.size(); ih++) {
	double cop = helices[ih].getCop();
	printf(" %6.0f,%6.1f ",helices[ih].getZ(),cop?1/cop:0);
      }
      printf("\n");
    }
    printf("==========\n");
  }
  else printf("========== No REFIT\n");
  return ret;
#endif
}










