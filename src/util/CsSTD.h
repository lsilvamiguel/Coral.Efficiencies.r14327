// $Id: CsSTD.h,v 1.6 2010/04/20 10:13:07 tnagel Exp $

/*!
   \file    CsSTD.h
   \brief   Compass STD include file
   Used to include, in a correct way, the STD include files.

   \author  Benigno Gobbo
   \version $Revision: 1.6 $
   \date    $Date: 2010/04/20 10:13:07 $
*/

#ifndef CsSTD_h
#define CsSTD_h

#include "coral_config.h"

#ifdef __HP_aCC
#   include <iostream.h>
#   include <fstream.h>
#   include <iomanip.h>
#   include <strstream.h>
#else
#   include <iostream>
#   include <fstream>
#   include <iomanip>
#   include <sstream>
#endif // __HP_aCC

#include <string>
#include <list>
#include <vector>
#include <map>

#endif // CsSTD_h
