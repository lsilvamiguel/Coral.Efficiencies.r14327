// $Id$

/*!
   \file    CsEvdisTrack.cc
   \brief   CORAL Event Display Package.
   \version $Revision$
   \author  Sebastian Uhl
   \date    $Date$
*/

#include "CsEvdisTrack.h"

#include <TEveVSDStructs.h>

#include <CsTrack.h>

CsEvdisTrack::CsEvdisTrack() : TEveTrack() {
}

void CsEvdisTrack::ImportTrack(const CsTrack* track, TEveTrackPropagator* prop) {
    assert(track->getHelices().size() >= 2);

    TEveRecTrack* eveTrack = new TEveRecTrack();
    const CsHelix& hlx = track->getHelices()[0];
    // charge convention of eve is wrong
    if (hlx.getCop() < 0)
        eveTrack->fSign = -1;
    else
        eveTrack->fSign = +1;
    const double p = eveTrack->fSign / hlx.getCop();
    const double len = sqrt(hlx.getDXDZ()*hlx.getDXDZ() + hlx.getDYDZ()*hlx.getDYDZ() + 1.);
    if (p==0.)
        eveTrack->fP.Set(hlx.getDXDZ(), hlx.getDYDZ(), 1.);
    else
        eveTrack->fP.Set(p * hlx.getDXDZ()/len, p * hlx.getDYDZ()/len, p * 1./len);
    eveTrack->fV.Set(hlx.getX()/10., hlx.getY()/10., hlx.getZ()/10.);

    SetTrackParams(TEveTrack(eveTrack, prop));

    for (size_t i=1; i<track->getHelices().size()-1; i++) {
        const CsHelix& hlx = track->getHelices()[i];
        TEvePathMark* pm = new TEvePathMark(TEvePathMark::kReference);
        const double p = eveTrack->fSign / hlx.getCop();
        const double len = sqrt(hlx.getDXDZ()*hlx.getDXDZ() + hlx.getDYDZ()*hlx.getDYDZ() + 1.);
        if (p==0.)
            pm->fP.Set(hlx.getDXDZ(), hlx.getDYDZ(), 1.);
        else
            pm->fP.Set(p * hlx.getDXDZ()/len, p * hlx.getDYDZ()/len, p * 1./len);
        pm->fV.Set(hlx.getX()/10., hlx.getY()/10., hlx.getZ()/10.);
        AddPathMark(*pm);
    }

    {
        const CsHelix& hlx = track->getHelices()[track->getHelices().size()-1];
        TEvePathMark* pm = new TEvePathMark(TEvePathMark::kReference);
        pm->fP.Set(0., 0., 0.);
        pm->fV.Set(hlx.getX()/10., hlx.getY()/10., hlx.getZ()/10.);
        AddPathMark(*pm);
    }

    MakeTrack();
}
