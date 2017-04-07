/*!
  
  Print information about
  hit with index in vHit "ih"
  
*/

#include "TDisplay.h"
#include "TEv.h"
#include <stdio.h>
#include <string.h>

void TDisplay::HitInfo(int ih)
{
  const TEv& ev = TEv::Ref();
  const THit& h = ev.vHit(ih);
  
  h.Print(2);

}

















