
/*!
  Wrapper for TAlgo::Field
  (called from grkuta.F)
*/

#include "TAlgo.h"

extern "C" void gufld_(double vecin[], double f[])
{
  TAlgo::Field(vecin,f);
}
