//-*-Mode: C++;-*-
#ifndef _GdbSTD_h_
#define _GdbSTD_h_

#if __GNUG__ >= 2
#  pragma interface
#endif

#include "coral_config.h"

//# ifdef COMPASS_USE_OSPACE_STD
#  include <ospace/std/vector>
#  include <ospace/std/list>
#  include <ospace/std/iostream>
#  include <ospace/std/iomanip>
//# else
//#  ifdef __HP_aCC
//#   include <iostream.h>
//#  include <iomanip.h>
//#  include <vector.h>
//#  else
//#   include <iostream>
//#  include <iomanip>
//#  include <vector>
//#  endif // __HP_aCC
//# endif // COMPASS_USE_OSPACE_STD


#endif
