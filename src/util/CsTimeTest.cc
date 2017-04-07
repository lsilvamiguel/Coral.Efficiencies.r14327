#include <iostream.h>
#include <sys/time.h>
#include <unistd.h>
#include "CsTime.h"

// A test program for the CsTime class

int main() {

  struct timeval tp;
# ifdef __hpux
  void *tzp = 0;
  gettimeofday( &tp, tzp );
# else
  struct timezone tzp;
  gettimeofday( &tp, &tzp );
# endif
  int time = tp.tv_sec;
  int usec = tp.tv_usec;
  
  cout << endl << "Local time, at present: "
       << time << " s, " << usec << " us" << endl; 
   
  cout << endl << "Object defined with the above quantities:" << endl;
  CsTime actualDate( time, usec );
  cout << dmy << actualDate << endl;
  cout << sec << actualDate << endl;
 
  cout << endl << "Object defined with Y, M, D, etc. :" << endl;
  CsTime aDate( 1999, 8, 12, 13, 30, 30 );
  cout << dmy << aDate << endl;
  cout << sec << aDate << endl;

  cout << endl << "Object defined with the default constructor:" << endl;
  CsTime thisDate;
  cout << dmy << thisDate << endl;
  cout << sec << thisDate << endl;

  return 0;

}
