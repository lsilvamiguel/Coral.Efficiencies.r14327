// $Id$

/*!
   \file    CsEvdisEvent.h
   \brief   CORAL Event Display Package.
   \version $Revision$
   \author  Sebastian Uhl
   \date    $Date$
*/

#ifndef CsEvdisEvent_h
#define CsEvdisEvent_h

#include <TEveElement.h>

class CsEvent;

class CsEvdisEvent : public TEveElementList {
    public:
                CsEvdisEvent();
    
        void    ImportEvent (const CsEvent* event);

    ClassDef(CsEvdisEvent, 0);
};

#endif // CsEvdisEvent_h
