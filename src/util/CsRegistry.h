// $Id: CsRegistry.h,v 1.2 2000/03/06 15:12:03 benigno Exp $

/*!
   \file    CsRegistry.h
   \brief   Compass Registry Class
   \author  Benigno Gobbo 
   \version $Revision: 1.2 $
   \date    $Date: 2000/03/06 15:12:03 $
*/

#ifndef CsRegistry_h
#define CsRegistry_h

#include "CsEndOfJob.h"
#include "CsEndOfEvent.h"
#include "CsStartOfRun.h"

/*! \class CsRegistry 
    \brief Compass Registry Class

    This class is used the interface to the CsRegistrySing class. The
    CsRegistrySing class is used to keep a list of packages to whom send a 
    particular call to a method in particular cases. At the moment it is used 
    to send an pkg->end() and pkg->eoe() method calls to all packages that 
    registered themselves for. It can be easilly extended to other event calls
    (e.g. end-of-burst call, end-of-run call, etc.).
*/

class CsRegistry {

 public:

  CsRegistry(); //!< Default Constructor

  /*! \fn bool   EOJRegistration( CsEndOfJob* ptr );
    \brief the registration method. It is used by packages to register 
    themselves for the end() method call.
    
    Returns \c true if the package is being correctly registered, \c
    false if the package was already registered for the call.
    \param ptr pointer to the object or singleton to be registered for
    the end() method call.
  */
  bool EOJRegistration( CsEndOfJob* ptr );

  /*! \fn bool   EOERegistration( CsEndOfEvent* ptr );
    \brief the registration method. It is used by packages to register 
    themselves for the eoe() method call.
    
    Returns \c true if the package is being correctly registered, \c
    false if the package was already registered for the call.
    \param ptr pointer to the object or singleton to be registered for
    the eoe() method call.
  */
  bool EOERegistration( CsEndOfEvent* ptr );

  /*! \fn bool   SORRegistration( CsStartOfRun* ptr );
    \brief the registration method. It is used by packages to register 
    themselves for the sor() method call.
    
    Returns \c true if the package is being correctly registered, \c
    false if the package was already registered for the call.
    \param ptr pointer to the object or singleton to be registered for
    the sor() method call.
  */
  bool   SORRegistration( CsStartOfRun* ptr );

};

#endif // CsRegistry_h
