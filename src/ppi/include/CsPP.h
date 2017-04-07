/*!
  \file      CsPP.h
  \brief     Compass PreProcessor Handling Class
  \author    Markus Kraemer
  \version   $Revision: 1.2 $
  \date      $Date: 2010/02/12 19:45:08 $

  \par       History:
  20091127   Design start of this class
*/

#ifndef __CsPP_H__
#define __CsPP_H__

#include "CsPPI.h"
#include "CsPPI_EC02time.h"
#include <list>

class CsPP {
  
  private:
//     CsPP(const string* ppiTableName);
    CsPP();
    virtual ~CsPP() {};

    std::list<CsPPI*> _CsPPI;
  

  public:
    
//     CsPP*               Instance(const string* ppiTableName);
    static CsPP*        Instance();
    int                 ProcessEvent();
    std::list<CsPPI*>&  GetPPI();
    int                 end();

  private:
    static CsPP*        _instance;

};

#endif // __CsPP_H__

