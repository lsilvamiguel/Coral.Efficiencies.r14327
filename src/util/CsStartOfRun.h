// $Id: CsStartOfRun.h,v 1.2 2000/04/18 13:44:39 benigno Exp $

/*!
   \file    CsStartOfRun.h
   \brief   Compass End Of Job Base Class
   \author  Benigno Gobbo 
   \version $Revision: 1.2 $
   \date    $Date: 2000/04/18 13:44:39 $
*/

#ifndef CsStartOfRun_h
#define CsStartOfRun_h

/*! \class CsStartOfRun 
    \brief Compass Start Of Run Base Class

    This virtual class is used in association with the CsRegistrySing 
    class. It must be inherited by all packages that want to register
    themselves for sor() method call at each start of run.
*/

class CsStartOfRun {

 public:

  virtual ~CsStartOfRun() {}

  /* \fn virtual bool sor();
     \brief the sor() method thet must be implemented by the packages that
     want to be called everytime a new run is encountered. 
  */
  virtual bool sor() = 0;

};

#endif // CsStartOfRun_h
