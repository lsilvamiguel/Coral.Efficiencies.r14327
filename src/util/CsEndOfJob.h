// $Id: CsEndOfJob.h,v 1.2 2000/04/18 13:44:39 benigno Exp $

/*!
   \file    CsEndOfJob.h
   \brief   Compass End Of Job Base Class
   \author  Benigno Gobbo 
   \version $Revision: 1.2 $
   \date    $Date: 2000/04/18 13:44:39 $
*/

#ifndef CsEndOfJob_h
#define CsEndOfJob_h

/*! \class CsEndOfJob 
    \brief Compass End Of Job Base Class

    This virtual class is used in association with the CsRegistrySing 
    class. It must be inherited by all packages that want to register
    themselves for end() method call just before the end of the job.
*/

class CsEndOfJob {

 public:

  virtual ~CsEndOfJob() {}

  /* \fn virtual bool end();
     \brief the end() method thet must be implemented by the packages that
     want to be called at end of job. 
  */
  virtual bool end() = 0;

};

#endif // CsEndOfJob_h
