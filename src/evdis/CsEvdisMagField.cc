// $Id$

/*!
   \file    CsEvdisMagField.cc
   \brief   CORAL Event Display Package.
   \version $Revision$
   \author  Sebastian Uhl
   \date    $Date$
*/

#include "CsEvdisMagField.h"

#include <CsGeom.h>

CsEvdisMagField::CsEvdisMagField() : TEveMagField() {
}

TEveVector CsEvdisMagField::GetField(float x, float y, float z) const {
    float bx, by, bz;

    CsGeom::Instance()->getCsField()->getField(x, y, z, bx, by, bz);

    // TEve's magnetic field orientation is inverted
    return TEveVector(-bx, -by, -bz);
}
