/*!
   \file CsFieldGridSM1.h
   \brief COMPASS Magnetic Filed Database Class
   \author  Takeaki TOEDA
   \version $Revision: 1.1 $
   \date    $Date: 2000/05/23 14:24:42 $
*/

#ifndef CsFieldGridSM1_h
#define CsFieldGridSM1_h
#include "HepODBMS/odbms/HepODBMS.h"                                           

/*! \class CsFieldGridSM1
 *  \brief Grid of magnetic field for SM1
 */
class CsFieldGridSM1{
public:
  CsFieldGridSM1();//!< Default Constructor
  CsFieldGridSM1(float x,float y,float z,float Bx,float By,float Bz);//Constructor

  /* \fn float getX();
    \brief Get X value
  */
  float getX(){return _x;};

  /* \fn float getY();
    \brief Get Y value
  */
  float getY(){return _y;};

  /* \fn float getZ();
    \brief Get Z value
  */

  /* \fn float getBx();
    \brief Get Bx value
  */
  float getBx(){return _Bx;};

  /* \fn float getBy();
    \brief Get By value
  */
  float getBy(){return _By;};

  /* \fn float getBz();
    \brief Get Bz value
  */
  float getBz(){return _Bz;};

  float getZ(){return _z;};
protected:
  float _x,_y,_z,_Bx,_By,_Bz;
};

#endif
