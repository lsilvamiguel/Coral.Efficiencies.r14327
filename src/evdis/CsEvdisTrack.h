// $Id$

/*!
   \file    CsEvdisTrack.h
   \brief   CORAL Event Display Package.
   \version $Revision$
   \author  Sebastian Uhl
   \date    $Date$
*/

#ifndef CsEvdisTrack_h
#define CsEvdisTrack_h

#include <TEveTrack.h>

class CsTrack;

class CsEvdisTrack : public TEveTrack {
    public:
        CsEvdisTrack();

        void    ImportTrack(const CsTrack* track, TEveTrackPropagator* prop);

    ClassDef(CsEvdisTrack, 0)
};

#endif // CsEvdisTrack_h
