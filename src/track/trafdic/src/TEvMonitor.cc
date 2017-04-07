// $Id: TEvMonitor.cc,v 1.9 2003/05/19 13:50:54 benigno Exp $

#include <iostream>
#include "TProfile.h"
#include "CsHistograms.h"
#include "CsEvent.h"
#include "Traffic.h"
#include "TEv.h"
#include "TSetup.h"
#include "TOpt.h"

/*! 
  \brief Traffic method for monitoring both MC and RD.
*/

/*
  Version specific to "lattice" alternative.
  Changes w/ respect to "traffic/TEvMCMonitor.cc":
  i) Allow "RDMonitor" w/ MC.
 */

// $Log: TEvMonitor.cc,v $
// Revision 1.9  2003/05/19 13:50:54  benigno
// Againg gcc2 -> gcc3
//
// Revision 1.8  2002/08/15 02:20:46  ybedfer
//  Disable "TraFDIc_HISTO_STAT".
//
// Revision 1.7  2002/07/23 13:37:21  ybedfer
//  Histogram: CPU vs. # hits, # hits vs. event #.
//
// Revision 1.6  2002/07/10 13:11:41  ybedfer
//  Histogram # clusters per event.
//
// Revision 1.5  2002/03/26 01:46:54  ybedfer
//  lattice version of traffic method for monitoring both MC and RD.
//


void TEv::Monitor()
{
  if(TOpt::Hist[1] == 0) return;

  if( IsMC() ) MCMonitor(); // MC specific 
  RDMonitor(); // RD specific 
  
  const TSetup& setup = TSetup::Ref();

  //----------------------------------------

  // Book histograms

  static CsHist1D *mm[2];

#ifdef TraFDIc_HISTO_STAT
  static CsHist2D *cpu_vs_nhits; static TProfile *nhits_vs_evt;
  static double cpu = 0;
#endif

  static bool first = true;
  if(first){
    first = false;
    CsHistograms::SetCurrentPath("/Traffic/Monitor");
    int Nplane = setup.vPlane().size();
    mm[0] =  new CsHist1D( "mm_00","N clusters / plane", Nplane,  0., Nplane);
    mm[1] =  new CsHist1D( "mm_01","N clusters in event", 100,  0., 10000);

#ifdef TraFDIc_HISTO_STAT
    cpu_vs_nhits = new CsHist2D("tScpu"  ,"CPU vs. #Hits" ,5,0,5000,100,0,10);
    nhits_vs_evt = new TProfile("tSnhits","#Hits vs. Evt#",100,0,20000,"S");
#endif

    CsHistograms::SetCurrentPath("/");
  } // end of booking block

  //----------------------------------------
  
  //
  //  Fill general purpose (MC/RD independent) monitoring histograms
  //
  
  for(int ip = 0; ip < int(setup.vPlane().size()); ip++){ //loop over planes
    mm[0]->Fill(ip+0.5, double(setup.vPlane(ip).vHitRef().size())); // Nhits/plane
  }
  
  mm[1]->Fill(vecHit.size()+0.5);

#ifdef TraFDIc_HISTO_STAT
  // CPU per event
  Traffic::Ref().Stopwatch.Stop(8);
  double cpu_prv = cpu; cpu = Traffic::Ref().Stopwatch.SumTime(8);
  cpu_vs_nhits->Fill((double)vecHit.size(),cpu-cpu_prv);
  nhits_vs_evt->Fill((double)CsEvent::Instance()->getEventNumberInBurst(),
		     (double)vecHit.size());
#endif 
}






