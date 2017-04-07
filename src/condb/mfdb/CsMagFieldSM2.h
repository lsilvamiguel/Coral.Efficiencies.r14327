/*!
   \file CsMagFieldSM2.h
   \brief COMPASS Magnetic Filed Database Class
   \author  Takeaki TOEDA
   \version $Revision: 1.1 $
   \date    $Date: 2000/05/23 14:24:43 $
*/

#ifndef CsMagFieldSM2_h
#define CsMagFieldSM2_h

/*! \class CsMagFieldSM2
 *  \brief Magnetic Field Class for SM2.
 */

class CsMagFieldSM2{
public:
  /*! \brief Following array contain polypomial coefficients 
    that are result of fitting procedure for field of SM2 magnet.
  */
  float FSMAX1[1770],FSMAY1[2360],FSMAZ1[3540],FSMBX1[472],
    FSMBY1[300],FSMCX1[84],FSMCY1[126],FSMCZ1[168],FSMDX1[312],
    FSMDY1[390],FSMDZ1[468],FSMA01,FSMA11;
};

#endif
