/*!
   \file CsFieldGridSol.h
   \brief COMPASS Magnetic Filed Database Class
   \author  Takeaki TOEDA
   \version $Revision: 1.1 $
   \date    $Date: 2000/05/23 14:24:42 $
*/

#ifndef CsFieldGridSol_h
#define CsFieldGridSol_h

/*! \class CsFieldGridSol
 *  \brief Grid of magnetic Field for Solenoid.
 */

class CsFieldGridSol{
public:
  CsFieldGridSol(); //!< Default Constructor
  CsFieldGridSol(float z,float r,float Bz,float Br); //!< Constructor

  /* \fn float getZ();
    \brief Get Z value
  */
  float getZ(){return _z;};

  /* \fn float getR();
    \brief Get R value
  */
  float getR(){return _r;};

  /* \fn float getBz();
    \brief Get Bz value
  */
  float getBz(){return _Bz;};

  /* \fn float getBr();
    \brief Get Br value
  */
  float getBr(){return _Br;};
protected:
  float _z,_r,_Bz,_Br;
};
#endif


