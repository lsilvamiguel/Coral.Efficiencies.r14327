// $Id$

/*!
   \file    CsEvdisMagField.h
   \brief   CORAL Event Display Package.
   \version $Revision$
   \author  Sebastian Uhl
   \date    $Date$
*/

#ifndef CsEvdisMagField_h
#define CsEvdisMagField_h

#include <TEveTrackPropagator.h>

class CsEvdisMagField : public TEveMagField {
    public:
                            CsEvdisMagField();

        virtual TEveVector  GetField(float x, float y, float z) const;

    ClassDef(CsEvdisMagField, 0)
};

#endif // CsEvdisMagField_h
