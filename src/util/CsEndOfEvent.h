// $Id: CsEndOfEvent.h,v 1.2 2000/04/18 13:44:39 benigno Exp $

/*!
   \file    CsEndOfEvent.h
   \brief   Compass End Of Event Base Class
   \author  Benigno Gobbo 
   \version $Revision: 1.2 $
   \date    $Date: 2000/04/18 13:44:39 $
*/

#ifndef CsEndOfEvent_h
#define CsEndOfEvent_h

/*! \class CsEndOfEvent 
    \brief Compass End Of Job Base Class

    This virtual class is used in association with the CsRegistrySing 
    class. It must be inherited by all packages that want to register
    themselves for eoe() method call just before the end of the job.
*/

class CsEndOfEvent {

 public:

  virtual ~CsEndOfEvent() {}

  /* \fn virtual bool eoe();
     \brief the end() method thet must be implemented by the packages that
     want to be called at end of job. 
  */
  virtual bool eoe() = 0;

};

#endif // CsEndOfEvent_h
