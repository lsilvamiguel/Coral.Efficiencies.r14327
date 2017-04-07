// $Id: CsEvdis.cc 13368 2012-04-20 08:53:24Z suhl $

/*!
   \file    CsEvdis.cc
   \brief   CORAL Event Display Package.
   \version $Revision: 13368 $
   \author  Sebastian Uhl
   \date    $Date: 2012-04-20 10:53:24 +0200 (Fri, 20 Apr 2012) $
*/

#include "CsEvdis.h"

#include <iostream>

#include <TApplication.h>
#include <TEveBrowser.h>
#include <TEveGeoNode.h>
#include <TEveManager.h>
#include <TEveTrans.h>
#include <TEveViewer.h>
#include <TGButton.h>
#include <TGeoBBox.h>
#include <TGeoBoolNode.h>
#include <TGeoCompositeShape.h>
#include <TGeoManager.h>
#include <TGeoNode.h>
#include <TGLViewer.h>
#include <TSystem.h>

#include <CsBMSDetector.h>
#include <CsErrLog.h>
#include <CsEvent.h>
#include <CsGeom.h>
#include <CsOpt.h>
#include <CsRegistrySing.h>

#include "CsEvdisEvent.h"
#include "CsEvdisMagField.h"

CsEvdis* CsEvdis::instance = NULL;

CsEvdis::CsEvdis() : CsEndOfEvent(), fInitializedEve(false), fMagField(NULL), fKeepOn(true), fKeepOnShowing(true) {
    if (instance!=NULL)
        CsErrLog::Instance()->msg(elFatal, __FILE__, __LINE__,
                                  "CsEvdis::CsEvdis: There can only be one event display.");
    instance = this;

    std::string filename;
    if (CsOpt::Instance()->getOpt("event display", "geometry", filename)) {
        fGeometryFile = filename;
    }
}

void CsEvdis::Init() {
    if (instance!=NULL)
        CsErrLog::Instance()->msg(elFatal, __FILE__, __LINE__,
                                  "CsEvdis::Init: Event display package has already been initialised.");

    CsRegistrySing::Instance()->EOERegistration(new CsEvdis());
}

CsEvdis* CsEvdis::Instance() {
    if (instance==NULL)
        CsErrLog::Instance()->msg(elFatal, __FILE__, __LINE__,
                                  "CsEvdis::Instance: Event display has to be initialised before an instance can be retrieved.");

    return instance;
}

void CsEvdis::InitEveGeom() {
    // save pointer to currect gGeoManager
    TGeoManager* oldGeoManager = gGeoManager;

    // special geometry for the event display:
    // using correct axis, not COMGeant definition
    gEve->RegisterGeometryAlias("Default", fGeometryFile);
    fHaveFullGeom = true;
    TGeoManager* newGeoManager(NULL);
    // ignore all ROOT messages while trying to load the geometry
    Int_t oldIgnoreLevel = gErrorIgnoreLevel;
    gErrorIgnoreLevel = kError+1;
    try {
        newGeoManager = gEve->GetDefaultGeometry();
    } catch (const TEveException& e) {
        CsErrLog::Instance()->msg(elError, __FILE__, __LINE__,
                                  "Caught exception while trying to load geometry for event display:\n"
                                  " %s\n"
                                  " The full geometry cannot be shown.", e.what());
        fHaveFullGeom = false;
    }
    // restore previous ROOT error message level
    gErrorIgnoreLevel = oldIgnoreLevel;

    if (fHaveFullGeom) {
        TGeoNode* top_node = newGeoManager->GetTopNode();
        TEveGeoTopNode* eve_top_node = new TEveGeoTopNode(newGeoManager, top_node);
        gEve->AddGlobalElement(eve_top_node);
    }

    // restore old gGeoManager used in tracking, ...
    gGeoManager = oldGeoManager;
}

void CsEvdis::InitEveGeomSimple() {
    TEveElementList* simpleGeom = new TEveElementList("simpleGeometry");
    gEve->AddGlobalElement(simpleGeom);

    // start with the big things:
    // * magnets
    //   magnet 0 -> target magnet (not yet done)
    //   magnet 1 -> SM1
    //   magnet 2 -> SM2
    assert(CsGeom::Instance()->getCsField()->getNumOfMags()==3);
    for (unsigned int i=1; i<3; i++) {
        const CsMagInfo& magInfo = CsGeom::Instance()->getCsField()->getMagInfo()[i];

        TString name;
        name.Form("SM%u", i);

        // TODO: hard-coded magnet sizes below
        const double xHalfSize[3] = { 0., 209.5,  255. };
        const double yHalfSize[3] = { 0., 209.25, 175. };
        const double zHalfSize[3] = { 0.,  55.,   200. };
        const double xHoleSize[3] = { 0., 114.5,  100. };
        const double yHoleSize[3] = { 0., 111.,    50. };

        // create shapes
        TGeoBBox* shBox  = new TGeoBBox(xHalfSize[i], yHalfSize[i], zHalfSize[i]);
        TGeoBBox* shHole = new TGeoBBox(xHoleSize[i], yHoleSize[i], 1.1*zHalfSize[i]);
        TGeoCompositeShape* shMagnet = new TGeoCompositeShape(name, new TGeoSubtraction(shBox, shHole));

        // create object to display
        TEveGeoShape* eveMagnet = new TEveGeoShape(name);
        eveMagnet->SetShape(shMagnet);
        eveMagnet->SetMainColor(kRed);
        eveMagnet->RefMainTrans().SetPos(magInfo.xcm/10., magInfo.ycm/10., magInfo.zcm/10.);

        // add to list
        simpleGeom->AddElement(eveMagnet);
    }

    // * absorbers

    // * individual detector planes
    const std::list<CsDetector*>& dets = CsGeom::Instance()->getDetectors();
    for (std::list<CsDetector*>::const_iterator it=dets.begin(); it!=dets.end(); it++) {
        CsDetector* det = *it;

        // skip BMS detectors, they would (over-)stretch the complete setup
        if (dynamic_cast<CsBMSDetector*>(det))
            continue;

        // create shape
        TGeoBBox* detBox = new TGeoBBox(det->getXsiz()/10., det->getYsiz()/10., det->getZsiz()/10.);

        // create object to display
        TEveGeoShape* detEve = new TEveGeoShape(det->GetTBName().c_str());
        detEve->SetShape(detBox);
        detEve->RefMainTrans().SetPos(det->getXcm()/10., det->getYcm()/10., det->getZcm()/10.);
        const double ang = atan2(det->getRotDRS()(2, 1), det->getRotDRS()(1, 1));
        detEve->RefMainTrans().SetRotByAngles(ang, 0., 0.);

        // add to list
        simpleGeom->AddElement(detEve);
    }
    // getDetectors
    // getRich1Detector
    // getCalorimeters
    // getMiscDetectors
    // CsDet::GetAllDetectors

    // by default do not draw the simple geometry
    //simpleGeom->DisableListElements();
}

void CsEvdis::InitEveGui() {
    TEveBrowser* browser = gEve->GetBrowser();

    // close the event display if the window is closed
    browser->Connect("CloseWindow()", "CsEvdis", this, "close()");

    // create an extra tab for controlling CORAL
    browser->StartEmbedding(TRootBrowser::kLeft);

    TGMainFrame* frmMain = new TGMainFrame(gClient->GetRoot(), 1000, 600);
    frmMain->SetWindowName("XX GUI");
    frmMain->SetCleanup(kDeepCleanup);

    fBtnNext = new TGTextButton(frmMain, "Next");
    frmMain->AddFrame(fBtnNext);
    fBtnNext->Connect("Clicked()", "CsEvdis", this, "next()");

    fBtnClose = new TGTextButton(frmMain, "Close");
    frmMain->AddFrame(fBtnClose);
    fBtnClose->Connect("Clicked()", "CsEvdis", this, "close()");

    fBtnQuit = new TGTextButton(frmMain, "Quit");
    frmMain->AddFrame(fBtnQuit);
    fBtnQuit->Connect("Clicked()", "CsEvdis", this, "quit()");

    frmMain->MapSubwindows();
    frmMain->Resize();
    frmMain->MapWindow();

    browser->StopEmbedding();
    browser->SetTabTitle("Event Control", 0);
}

void CsEvdis::InitEve() {
    // save pointer to current gGeoManager
    TGeoManager* oldGeoManager = gGeoManager;

    // a valid non-batch gApplication object must exist, the default created
    // during reading the ROOT geometry is a batch one
    if ((!gApplication) || (gApplication && gApplication->TestBit(TApplication::kDefaultApplication)))
        new TApplication("eventDisplay", 0, 0);

    // COMPASS is a bit longer than LHC experiments
    TEveTrackPropagator::fgEditorMaxR =  750.;
    TEveTrackPropagator::fgEditorMaxZ = 6500.;

    TEveManager::Create();

    InitEveGeom();
    InitEveGeomSimple();
    InitEveGui();

    gEve->Redraw3D(true);

    fMagField = new CsEvdisMagField();

    fInitializedEve = true;

    assert(oldGeoManager == gGeoManager);
}

void CsEvdis::next() {
    gApplication->Terminate();
}

void CsEvdis::close() {
    fKeepOnShowing = false;

    TEveBrowser* browser = gEve->GetBrowser();
    browser->Disconnect("CloseWindow()", this, "close()");

    gEve->GetBrowser()->UnmapWindow();
    if (gTQSender == fBtnClose)
        TEveManager::Terminate();

    gApplication->Terminate();
}

void CsEvdis::quit() {
    fKeepOn = false;
    fKeepOnShowing = false;
    gApplication->Terminate();
}

bool CsEvdis::eoe() {
    // check whether events should still be displayed
    if (!fKeepOnShowing)
        return fKeepOn;

    if (!fInitializedEve)
        InitEve();

    CsEvdisEvent* event = new CsEvdisEvent();
    event->ImportEvent(CsEvent::Instance());
    gEve->AddElement(event);

    // make a backup of the current gGeoManager pointer
    TGeoManager* oldGeoManager = gGeoManager;

    gApplication->Run(true);

    // gGeoManager might get changed if e.g. the window is closed, so restore it here
    gGeoManager = oldGeoManager;

    return fKeepOn;
}

CsEvdisMagField* CsEvdis::GetMagField() const {
    return fMagField;
}
