/*!

  This is the place for TRAFFIC users/developers
  event display extentions.
  This method is called at the end-of-event.

*/

#include <iostream>
#include "TDisplay.h"
#include "TSetup.h"
#include "TEv.h"
#include "TOpt.h"
#include "higz.h"

void TDisplay::UserDraw()
{
  if(TOpt::Graph[0] <= 0) return; // do nothing if the graphics is OFF

  TEv& ev = TEv::Ref();
  const TSetup& setup = TSetup::Ref();

  // ....


  return;
}













