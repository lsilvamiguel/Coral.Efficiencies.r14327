/*!
   \file    CsRCDetectors.cc
   \----------------------
   \brief   CsRCDetectors class implementation.
   \author  Paolo Schiavon
   \version 0.01
   \date    11 August 2000, rev. August 2005
*/


//---------------------------
#include <iostream>
#include <fstream>
#include <iomanip>
#include <sstream>
#include <string>
//----------------------------

#  include "CsInit.h"
#  include "CsOpt.h"
#  include "CsGeom.h"

#  include "CsRCDetectors.h"
#  include "CsRCPhotonDet.h"
#  include "CsRCCathode.h"
#  include "CsRCMirrors.h"

#  include "CsRichOne.h"
//----------------------------

#  include <ostream>
#  include <cstdio>

  using namespace std;
  using namespace CLHEP;

  CsRCDetectors* CsRCDetectors::instance_ = 0;

//===========================================================================
  CsRCDetectors* CsRCDetectors::Instance() {
//------------------------------------------
    if( instance_ == 0 ) instance_ = new CsRCDetectors();
    return instance_;
  }

//===========================================================================
  CsRCDetectors::CsRCDetectors() {
//--------------------------------
    lPhoDet_.clear();
    lCathodes_.clear();
    nCatPMT_.clear();
    nCatAPV_.clear();
  }

//===========================================================================
  void CsRCDetectors::setPhotonDet( string name,
//----------------------------------------------
              double xDet0, double  yDet0, double zDet0,
	      HepMatrix rotMatrix ) {

    int kDet = lPhoDet_.size();
    lPhoDet_.push_back( new CsRCPhotonDet( kDet, name, 
				     xDet0, yDet0, zDet0, rotMatrix ) );
  }

//===========================================================================
  void CsRCDetectors::setCathode( int id, string TBname, string name,
//-------------------------------------------------------------------
              double xOffCat0, double yOffCat0, double zOffCat0, 
              int nPadx, int nPady, double padx, double pady,
	      double ddQzW, double ddGap, HepMatrix corrRotMx ) {

    int kCat = lCathodes_.size();
    lCathodes_.push_back( new CsRCCathode( kCat, id, TBname, name,
				     xOffCat0, yOffCat0, zOffCat0,
				     nPadx, nPady, padx, pady, ddQzW, ddGap,
				     corrRotMx ) );

    if( lCathodes_.back()->isPMT() ) {
      CsRichOne::Instance()->setUpRICHJob();
      nCatPMT_.push_back( kCat );
    } else {
      nCatAPV_.push_back( kCat );
    }

  }

//===========================================================================
  CsRCDetectors::CsRCDetectors( const CsRCDetectors &det ) {
//----------------------------------------------------------
    cout << "RICHONE : CsRCDetectors CopyConstructor" << endl;
    instance_ = det.instance_;
    lPhoDet_ = det.lPhoDet_;
    lCathodes_ = det.lCathodes_;
    zEntrWind_ = det.zEntrWind_;
    zExitWind_ = det.zExitWind_;
    nCatPMT_ = det.nCatPMT_;
    nCatAPV_ = det.nCatAPV_;
  }

//===========================================================================
  void CsRCDetectors::printDetAll()  const {
//------------------------------------------
    cout << endl;
    cout << "CsRCDetectors::printDetAll() : Photon Detector Geometry" << endl
         << "-------------------------------------------------------" << endl;
    list<CsRCPhotonDet*>::const_iterator id;
    for( id=lPhoDet_.begin(); id!=lPhoDet_.end(); id++ ) {

      (*id)->print();
    }
  }

//===========================================================================
  void CsRCDetectors::printCatAll()  const {
//------------------------------------------
    cout << endl;
    cout << "CsRCDetectors::printCatAll() : Cathode Geometry" << endl
         << "-----------------------------------------------" << endl;
    list<CsRCCathode*>::const_iterator ic;
    for( ic=lCathodes_.begin(); ic!=lCathodes_.end(); ic++ ) {

      (*ic)->print();
    }
  }

//===========================================================================
  CsRCDetectors::~CsRCDetectors() {
//---------------------------------
    lPhoDet_.clear();
    lCathodes_.clear();
    nCatPMT_.clear();
    nCatAPV_.clear();
  }


//==========================================================================
  void CsRCDetectors::setLocalGeo() {
//-----------------------------------

//--- date  11/02/01


//--- MRS : Main Reference System
//    ---------------------------
//--- MWR : Mirror Working Reference, up and down :
//    ---------------------------------------------
//    ref. origin in mirror centre, z-axis orthogonal to detector
//    surface.


      CsOpt* opt = CsOpt::Instance();
      bool boo;
      int kOpw;
      bool lprint = false;
//    ( print... = 0 : NO print,  print... = 1 : print )
      boo = opt->CsOpt::getOpt( "RICHONE", "PrintDetMirr", kOpw );
      if( boo ) { if( kOpw > 0 ) { lprint = true; } }

      int k = 0;

//--- NOMINAL MIRRORS :
//    -----------------
//--- 'up' = 0;  'down' = 1 :
//    -----------------------
      CsRCMirrors* mirr = CsRCMirrors::Instance();
      if( lprint ) {
	mirr->printNomAll();
	mirr->printEleAll();
      }

      list<CsRCMirrorNom*> lMirrNom = mirr->lMirrNom();
      int nMirror = lMirrNom.size();
      list<CsRCMirrorElem*> lMirrEle = mirr->lMirrEle();
      if( nMirror == 0 || lMirrEle.size() == 0 ) {
        string str = "RICHONE, CsRCDetectors::setLocalGeo() : ";
	string err = "wrong detectors.dat for RICH1";
	str.append( err );
        CsErrLog::Instance()->mes( elFatal, str );
      }
//--- mirror centres and radii
      list<CsRCMirrorNom*>::iterator in;
      Hep3Vector vC0v[nMirror];
      float RRv[nMirror];
      k = 0;
      for ( in=lMirrNom.begin(); in!=lMirrNom.end(); in++ ) {
	vC0v[k] = (*in)->vC0();
      	RRv[k]  = (*in)->RR();
	k++;
      }

//--- DETECTORS :
//    -----------
//--- 'up' = 0;  'down' = 1 :
//    -----------------------
      CsRCDetectors* dets = CsRCDetectors::Instance();

      list<CsRCPhotonDet*> lPhoDet = dets->lPhoDet();
      int nDetector = lPhoDet.size();
      list<CsRCCathode*> lCathodes = dets->lCathodes();
      if( nDetector == 0 || lCathodes.size() == 0 ) {
        string str = "RICHONE, CsRCDetectors::setLocalGeo() : ";
	string err = "wrong detectors.dat for RICH1";
	str.append( err );
        CsErrLog::Instance()->mes( elFatal, str );
      }
//--- detector centres and directions (normal to)
      list<CsRCPhotonDet*>::iterator id;

      Hep3Vector vDet0v[nDetector];
      Hep3Vector vDcDetv[nDetector];
      k = 0;
      for( id=lPhoDet.begin(); id!=lPhoDet.end(); id++ ) {
	vDet0v[k] = (*id)->vDet0();
	vDcDetv[k] = (*id)->vDcDet();
      	k++;
      }

//--- RICH entrance window z-position 
//    ( 6156.3 mm, RICH1Upgraded )
      //zEntrWind_ = lPhoDet.front()->vDet0().z() + 65.;
      zEntrWind_ = lPhoDet.front()->vDet0In().z() + 65.;
//@@---------------------------------------------------

//--- RICH exit window z-position
//    ( 9366.3 mm, RICH1Upgraded )
      //zExitWind_ = lPhoDet.front()->vDet0().z() + 3275.;
      zExitWind_ = lPhoDet.front()->vDet0In().z() + 3275.;
//@@-----------------------------------------------------

      if( nDetector != nMirror ) {
        string str = 
	  "RICHONE, CsRCDetectors::setLocalGeo() : wrong geometry!";
	CsErrLog::Instance()->mes( elFatal, str );
      }

//--- transform detector centres to MWR's, up and down :
//    --------------------------------------------------
      Hep3Vector vDetWw;
      k = 0;
      for( id=lPhoDet.begin(); id!=lPhoDet.end(); id++ ) {
	HepVector vddr( 3, 0 );
	for( int j=0; j<3; j++ ) vddr[j] = (vDet0v[k] - vC0v[k])[j];
	HepVector vDetW = (*id)->rotMatrix() * vddr;
	for( int j=0; j<3; j++ ) vDetWw[j] = vDetW[j];
	(*id)->setDetW( vDetWw );
	//cout << "vDetW   " << vDetWw << endl;
	k++;
      }

//--- transform cathode offsets to MWR's, up and down :
//    -------------------------------------------------
      list<CsRCCathode*>::iterator ic;
      int nCathode = lCathodes.size();
      Hep3Vector vOffCatWw;
      int kc = 0;
      for( ic=lCathodes.begin(); ic!=lCathodes.end(); ic++ ) {
	Hep3Vector vOffCat0 = (*ic)->vCat0();
//@@----------------------------------------
        if( (*ic)->name()[0] == 'U' ) k = 0;
        if( (*ic)->name()[0] == 'D' ) k = 1;
	HepVector vddr( 3, 0 );
	for( int j=0; j<3; j++ ) vddr[j] = (vOffCat0 - vC0v[k])[j];
	HepVector vOffW = (*ic)->rotMatrix() * vddr;
	for( int j=0; j<3; j++ ) vOffCatWw[j] = vOffW[j];
	(*ic)->setOffCatW( vOffCatWw );
	double lim = 0.;
	lim = vOffCatWw.x() - (*ic)->hCatx();
	(*ic)->setXLimMnW( lim );
	lim = vOffCatWw.x() + (*ic)->hCatx();
	(*ic)->setXLimMxW( lim );
	lim = vOffCatWw.y() - (*ic)->hCaty();
	(*ic)->setYLimMnW( lim );
	lim = vOffCatWw.y() + (*ic)->hCaty();
	(*ic)->setYLimMxW( lim );
	kc++;
      }

//--- compute yDRS - yMWR ( Hits <-> Digits ) :
//    -----------------------------------------
      k = 0;
      for( id=lPhoDet.begin(); id!=lPhoDet.end(); id++ ) {
        double ww = (vDet0v[k].z() - vC0v[k].z()) * vDcDetv[k].z();
                  + (vDet0v[k].y() - vC0v[k].y()) * vDcDetv[k].y();
        double zP = ww * vDcDetv[k].z() + vC0v[k].z();
        double yP = ww * vDcDetv[k].y() + vC0v[k].y();
        double zQ = vDet0v[k].z() + vDet0v[k].y() * vDcDetv[k].y() /
                    vDcDetv[k].z();
        double yQ = 0.;
	double dyDRSMWRw = sqrt( (zP-zQ)*(zP-zQ) + (yP-yQ)*(yP-yQ) );
	if( k == 1 ) dyDRSMWRw *= -1.;
	(*id)->setDRSMWR( dyDRSMWRw );
	k++;
      }

      if( lprint ) {
	dets->printDetAll();
	dets->printCatAll();
      }

  }


//==========================================================================
  Hep3Vector CsRCDetectors::vImpDetW( const Hep3Vector vPoIn0,
//------------------------------------------------------------
				      const Hep3Vector vDcIn0 ) const {

//- Paolo  -  September 2005

//- particle or photon impact on detector (MWR) :
//  ---------------------------------------------

    Hep3Vector vPoDetW( 0.);
//- MRS :
//  -----
    CsRCMirrors *mirr = CsRCMirrors::Instance();
    std::list<CsRCMirrorNom*> lMirrNom = mirr->lMirrNom();
    CsRCMirrorNom* mino = NULL;
    CsRCPhotonDet* deno = NULL;
    int kDetPart = 0;
    mino = lMirrNom.front();
    deno = lPhoDet_.front();
    Hep3Vector vC0( mino->vC0() );
    double RR = mino->RR();
//- impact on mirror :
    Hep3Vector DD( vPoIn0 - vC0 );
    double dot = vDcIn0 * DD;
    double delta = dot*dot - ( DD*DD - RR*RR );
    if( delta < 0. ) return  vPoDetW;
    double norm = sqrt( delta ) - dot;
    Hep3Vector vPoMir0 = vPoIn0 + norm * vDcIn0;
    if( vPoMir0.y() < 0. ) kDetPart = 1;
    if( kDetPart == 1 ) {
      mino = lMirrNom.back();
      deno = lPhoDet_.back();
      vC0 = mino->vC0();
      RR = mino->RR();
      Hep3Vector DD( vPoIn0 - vC0 );
      dot = vDcIn0 * DD;
      delta = dot*dot - ( DD*DD - RR*RR );
      if( delta < 0. ) return  vPoDetW;
      norm = sqrt( delta ) - dot;
      vPoMir0 = vPoIn0 + norm * vDcIn0;
    }
    Hep3Vector vDcDe = vDcDetv( kDetPart );
    Hep3Vector vDcNorm0 = (1./RR) * (vPoMir0 - vC0);
    double cosMir = vDcNorm0 * vDcIn0;
    Hep3Vector vDcRefl0 = 2.*cosMir * vDcNorm0 - vDcIn0;
    norm = ( ( vDet0v( kDetPart ) - vPoMir0 ) * vDcDe ) /
      ( vDcRefl0 * vDcDe );
    Hep3Vector vPoDet0 = vPoMir0 + norm * vDcRefl0;
    HepVector vrot( 3, 0 );
    for( int j=0; j<3; j++ ) vrot[j] = (vPoDet0 - vC0)[j];
    HepVector vvW = deno->rotMatrix() * vrot;
    vPoDetW.setX( vvW( 1 ) );
    vPoDetW.setY( vvW( 2 ) );
    vPoDetW.setZ( vvW( 3 ) );

    return  vPoDetW;
  }


//==========================================================================
  CsRCCathode* CsRCDetectors::ptrToCat( int cat ) {
//-------------------------------------------------

//- Paolo  -  August 2006
//  moved from inline

    std::list<CsRCCathode*>::const_iterator ic;
    CsRCCathode* ptr = NULL;
    for( ic=lCathodes_.begin(); ic!=lCathodes_.end(); ic++ ) {
      if( (*ic)->number() == cat ) {
	ptr = (*ic);
	break;
      }
    }
    if( ptr == NULL ) {
      std::cout << std::endl;
      std::cout << " RICHONE, CsRCDetectors : Cathode number " << cat;
      std::cout << "   richEvent  ";
      std::cout << CsRichOne::Instance()->kEvent() << std::endl;
      std::string str =
        "RICHONE, CsRCDetectors::ptrToCat() : wrong cathode number";
      //CsErrLog::Instance()->mes( elFatal, str );
      CsErrLog::Instance()->mes( elError, str );
    }

    return ptr;
  }
