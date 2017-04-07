/*!
  Produce combination of int indices like
  inside N nested loops.

  \param N - number of nested loops.
  \param ns[k] (k = [0;N-1])  - number of cycles in the loop # k.
  If some ns[k] == 0 no combination produced at all!
  \param i[k] (k = [0;N-1]) - running index pattern  ( 0 <= i[k] < ns[k] ). 
  Returns new one every call.

  \return "false" if all combinations had been done

*/

#include "TAlgo.h"

bool TAlgo::NLoopComb(int N, int ns[], int i[])
{
  int l=0;
  static bool first=true;
  if(first){ // first call
    first=false;
    for(int k = 0; k < N; k++) {
      if(ns[k] <= 0) goto end;
      i[k] = 0;
    }
    return(true);
  }

  l = N-1;
  while(i[l] == ns[l]-1){
    i[l] = 0; l--;
    if(l == -1) goto end;
  }
  i[l]++;
  return(true);

 end:
  first=true;
  return(false);

}
