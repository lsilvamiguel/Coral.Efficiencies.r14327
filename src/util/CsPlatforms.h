#ifndef CsPlatform_h
#define CsPlatform_h

// $Id: CsPlatforms.h,v 1.3 2004/10/25 12:56:20 benigno Exp $

/*!
   \file    CsPlatforms.h
   \brief   Compass PLatform definer
   Sets the CS_PLATFORM macro according to the machine OS/C++ compiler.

   \author  Massimo Lamanna
   \version $Revision: 1.3 $
   \date    $Date: 2004/10/25 12:56:20 $
*/

#include "coral_config.h"

#undef CS_PLATFORM

//#if defined (_WIN32) && defined (_M_IX86) && defined (_MSC_VER)
//#define CS_PLATFORM "WindowsNT_VC++";
//#define __wnt
//#endif
#if defined (__sun) && defined (__SUNPRO_CC)
#define CS_PLATFORM "Solaris_CC"
#endif
#if defined (__osf__) && defined (__DECCXX)
#define CS_PLATFORM "OSF_cxx"
#endif
//#if defined (__AIX) && defined (__cplusplus) && !defined (__GNUG__)
//#define CS_PLATFORM "AIX_xlC"
//#endif
#if defined (__linux) && ( defined (__i386) || defined (__x86_64) )&& defined (__GNUG__)
#define CS_PLATFORM "Linux_egcs"
#endif
#if defined (__hpux) && defined (__HP_aCC)
#define CS_PLATFORM "HP-UX_aCC"
#define __unix 1
#endif

#ifndef CS_PLATFORM
#error "Unsupported platform"
#endif

#endif // CsPlatform_h
