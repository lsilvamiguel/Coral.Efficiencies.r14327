// $Id: CsDecMap.h,v 1.1 2001/02/05 11:25:54 hpereira Exp $

/*!
   \file    CsDecMap.h
   \brief   Compass Decoding Map structure for SaclayMicromegas and Drift Chambers.
   \author  Hugo Pereira
   \version $Revision: 1.1 $
   \date    $Date: 2001/02/05 11:25:54 $
*/

/*! \struct CsDecMap
    \brief decoding map structure. Define a segment of consecutive channels
    and the corresponding wire numbers.
    \param uid0 first uniq ID of the segment
    \param uid1 last uniq ID of the segment
    \param n_w number of wires in segment
    \param first_w wire number corresponding to uid0
    \param step_w wire increment from one channel to the next one.
    \param step_c channel increment from one wire to the next one.
*/

#ifndef CsDecMap_h
#define CsDecMap_h
#include "CsTypes.h"
typedef struct {
  uint32 uid0;
  uint32 uid1;
  int n_w;
  int first_w;
  int step_w;
  int step_c;
} CsDecMap;
#endif
