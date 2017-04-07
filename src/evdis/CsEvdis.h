// $Id: CsEvdis.h 13368 2012-04-20 08:53:24Z suhl $

/*!
   \file    CsEvdis.h
   \brief   CORAL Event Display Package.
   \version $Revision: 13368 $
   \author  Sebastian Uhl
   \date    $Date: 2012-04-20 10:53:24 +0200 (Fri, 20 Apr 2012) $
*/

#ifndef CsEvdis_h
#define CsEvdis_h

#include <string>

#include <Rtypes.h>

#include <CsEndOfEvent.h>

class CsEvdisMagField;
class TGTextButton;

class CsEvdis : public CsEndOfEvent {
    public:
                            CsEvdis();

        virtual bool        eoe();

        // static parts of this class
        static void         Init();
        static CsEvdis*     Instance();
        static CsEvdis*     instance;

        // methods called via GUI actions
        void                next();
        void                close();
        void                quit();

        // getters
        CsEvdisMagField*    GetMagField() const;

    private:
        void                InitEve();
        void                InitEveGeom();
        void                InitEveGeomSimple();
        void                InitEveGui();

        bool                fInitializedEve;

        std::string         fGeometryFile;
        bool                fHaveFullGeom;

        TGTextButton*       fBtnNext;
        TGTextButton*       fBtnClose;
        TGTextButton*       fBtnQuit;

        CsEvdisMagField*    fMagField;

        bool                fKeepOn;
        bool                fKeepOnShowing;

    ClassDef(CsEvdis, 0)
};

#endif // CsEvdis_h
