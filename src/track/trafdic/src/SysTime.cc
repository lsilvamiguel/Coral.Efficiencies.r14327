/*!

  Wrapper for sys/times function. 
  (To avoid name conflict of sys/times.h with STL on old compilers) 
  
  Return real time and cpu time (user+system) in float& parameter
  Times are in seconds. 

*/

#include <sys/times.h>
#include <unistd.h>
#include <iostream>

float SysTime (float& systime)
{
  struct tms cpt; 
  int    tic = sysconf(_SC_CLK_TCK);
  float  rt  = (float)times(&cpt)/tic;
  systime    = (float)(cpt.tms_utime + cpt.tms_stime) / tic;
  return rt;
}

void Wait(float sec)
{
  struct tms cpt; 
  int    tic = sysconf(_SC_CLK_TCK);
  float  rt  = (float)times(&cpt)/tic;
  while( float(times(&cpt)/tic) < (rt+sec) ){};

}
