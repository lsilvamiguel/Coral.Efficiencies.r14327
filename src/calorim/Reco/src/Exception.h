/*
   $Source: /afs/cern.ch/user/g/gobbo/w0/coral/src/calorim/Reco/src/Exception.h,v $
   $Date: 2010/12/10 14:42:56 $
   $Revision: 1.9 $
   -------------------------------------------------------------------------

   Authors:
     Alexander Zvyagin   ( Alexander.Zviagine@cern.ch, Zvyagin@mx.ihep.su )

   Copyright(C): 2000  A.Zvyagin

     This library is free software; you can redistribute it and/or
     modify it under the terms of the GNU Library General Public
     License as published by the Free Software Foundation; either
     version 2 of the License, or (at your option) any later version.

     This library is distributed in the hope that it will be useful,
     but WITHOUT ANY WARRANTY; without even the implied warranty of
     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
     Library General Public License for more details.

     You should have received a copy of the GNU Library General Public
     License along with this library; if not, write to the Free
     Software Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
*/

#ifndef Reco_Exception__include
#define Reco_Exception__include

#include <cstdarg>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <exception>

namespace Reco {

class Exception : public std::exception
{
 public:
  ~Exception (void) throw() { free(buf); }
  Exception (void) : buf(NULL) {};
  Exception (const Exception &e) : std::exception(e) { buf = strdup(e.buf); }
  Exception (const char* name, ...) {
    va_list ap;
    va_start(ap, name);
    vasprintf(&buf, name, ap);
    va_end(ap);
  }
  const char* name (void) const         { return buf; }
  const char* what (void) const throw() { return buf; } // virtual method of std::exception
 private:
  char           *buf;       // Exception text.
};

}

#endif // Reco_Exception__include
