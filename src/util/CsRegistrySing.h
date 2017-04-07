// $Id: CsRegistrySing.h,v 1.4 2003/03/13 16:57:32 benigno Exp $

/*!
   \file    CsRegistrySing.h
   \brief   Compass Registry Singleton Class
   \author  Benigno Gobbo 
   \version $Revision: 1.4 $
   \date    $Date: 2003/03/13 16:57:32 $
*/

#ifndef CsRegistrySing_h
#define CsRegistrySing_h

#include "CsSTD.h"
#include "CsEndOfJob.h"
#include "CsEndOfEvent.h"
#include "CsStartOfRun.h"

/*! \class CsRegistrySing 
    \brief Compass Registry Class

    This class is used to keep a list of packages to whom send a particular
    call to a method in particular cases. At the moment it is used to send
    an pkg->end() method call to all packages that registered themselves for. 
    It can be easilly extended to other event calls (e.g. end-of-burst call,
    start-of-run call, etc.).
*/

class CsRegistrySing {

 public:

  /*! \fn static CsRegistrySing* Instance();
    \brief singleton instantiation (but first).
  */
  static CsRegistrySing* Instance();

  /*! \fn bool   EOJRegistration( CsEndOfJob* ptr );
    \brief the registration method. It is used by packages to register 
    themselves for the end() method call.
    
    Returns \c true if the package is being correctly registered, \c
    false if the package was already registered for the call.
    \param ptr pointer to the object or singleton to be registered for
    the end() method call.
  */
  bool   EOJRegistration( CsEndOfJob* ptr );

  /*! \fn bool   EOERegistration( CsEndOfEvent* ptr );
    \brief the registration method. It is used by packages to register 
    themselves for the eoe() method call.
    
    Returns \c true if the package is being correctly registered, \c
    false if the package was already registered for the call.
    \param ptr pointer to the object or singleton to be registered for
    the eoe() method call.
  */
  bool   EOERegistration( CsEndOfEvent* ptr );

  /*! \fn bool   SORRegistration( CsStartOfRun* ptr );
    \brief the registration method. It is used by packages to register 
    themselves for the sor() method call.
    
    Returns \c true if the package is being correctly registered, \c
    false if the package was already registered for the call.
    \param ptr pointer to the object or singleton to be registered for
    the sor() method call.
  */
  bool   SORRegistration( CsStartOfRun* ptr );

  /*! \fn bool   callEndMethods();
    \brief Sends a pkg->end() call to all registered packages. 

    Returns \c true if all calls returns \c true, returns \c false if at least
    a package returns \c false.
  */
  bool   callEndMethods();

  /*! \fn bool   callEoeMethods();
    \brief Sends a pkg->eoe() call to all registered packages. 

    Returns \c true if all calls returns \c true, returns \c false if at least
    a package returns \c false.
  */
  bool   callEoeMethods();

  /*! \fn bool   callSorMethods();
    \brief Sends a pkg->sor() call to all registered packages. 

    Returns \c true if all calls returns \c true, returns \c false if at least
    a package returns \c false.
  */
  bool   callSorMethods();

 protected:

  CsRegistrySing(); //!< Default Constructor

 private:

  static CsRegistrySing* _instance;    //<! The singleton static pointer
  std::list<CsEndOfJob*>      _EOJRegister; //<! The list of registered for end
  std::list<CsEndOfEvent*>    _EOERegister; //<! The list of registered for eoe
  std::list<CsStartOfRun*>    _SORRegister; //<! The list of registered for sor

};

#endif // CsRegistrySing_h
