// $Id$

/*!
   \file    CsEvdisEvent.cc
   \brief   CORAL Event Display Package.
   \version $Revision$
   \author  Sebastian Uhl
   \date    $Date$
*/

#include "CsEvdisEvent.h"

#include <TEveCalo.h>
#include <TEveCaloData.h>
#include <TEveTrackPropagator.h>
#include <TEveTrans.h>

#include <CsCalorimeter.h>
#include <CsEvent.h>
#include <CsGeom.h>
#include <CsTrack.h>

#include "CsEvdis.h"
#include "CsEvdisMagField.h"
#include "CsEvdisTrack.h"

CsEvdisEvent::CsEvdisEvent() : TEveElementList("Event") {
}

void CsEvdisEvent::ImportEvent(const CsEvent* event) {
    { // import tracks
        TEveTrackList* tracks = new TEveTrackList("Tracks");
        AddElement(tracks);

        TEveTrackPropagator* prop = tracks->GetPropagator();
        prop->SetMaxR( 750.);
        prop->SetMaxZ(6500.);
        prop->SetMagFieldObj(CsEvdis::Instance()->GetMagField(), false);

        for (std::list<CsTrack*>::const_iterator it=event->getTracks().begin(); it!=event->getTracks().end(); it++) {
            CsEvdisTrack* track = new CsEvdisTrack();
            track->ImportTrack(*it, prop);
            tracks->AddElement(track);
        }
    }
}
