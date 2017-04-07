// $Id: TDisplayDrawText.cc,v 1.4 2010/02/03 18:22:23 suhl Exp $

/*!
  Friend function of class TEv
  to draw text in the event display
*/

#include <iostream>
#include <stdio.h>
#include "TDisplay.h"
#include "TOpt.h"
#include "TSetup.h"
#include "TEv.h"
#include "TConstants.h"
#include "higz.h"
#include "CsEvent.h"

/*
  Version specific to "lattice" alternative.
  Changes w/ respect to "traffic/TDisplayDrawText":
   i) Display event number as "event# in run (burst#, event# in burst)"
   ii) Display event time.
*/

void TDisplay::DrawText()
{
  TEv &ev = TEv::Ref();
  bool movie = TOpt::Graph[7]>0;

  ISELNT(0); ISTXCI(6);

  int run = ev.Run(), evNum = ev.Event();
  int evInBurst = ev.EventInBurst(), burst = ev.ptrEvt()->getBurstNumber();
  char str[120];
  sprintf(str,"Run %3u Event %9d (%3d,%5d) %.1f ns Trigger(s) %s  Nhits %4zu",
	  run,evNum,burst,evInBurst,ev.GetEventTime(),
	  (ev.TrigMaskString()).c_str(),ev.vHit().size());

  if (movie) {
    IGTEXT(0.9, 0.01, str, 0.01, 0., "R");
    sprintf(str,"Projection %6.1f deg. ", AngProj);
    IGTEXT(0.15, TConstants_A4-0.01, str, 0.008, 0., "L");
    ISTXCI(4);
    sprintf(str,"TRAFFIC event display");
    IGTEXT(0.01, 0.01, str, 0.008, 0., "L");
    sprintf(str, "(%03u) ", this->ev_count);
    IGTEXT(1-0.01, 0.01, str, 0.006, 0., "R");
  }
  else {
    IGTEXT(0.9, TConstants_A4-0.01, str, 0.01, 0., "R");

    sprintf(str,"Projection %6.1f deg. ", AngProj);
    IGTEXT(0.15, TConstants_A4-0.01, str, 0.008, 0., "L");

    ISTXCI(4);
    sprintf(str,"TRAFFIC (version %4.2f) event display",TConstants_TraFFiC_Version);
    IGTEXT(0.01, 0.01, str, 0.008, 0., "L");

    sprintf(str, "(%03u) ", this->ev_count);
    IGTEXT(1-0.01, 0.01, str, 0.006, 0., "R");
  }
  ISELNT(10);
}
