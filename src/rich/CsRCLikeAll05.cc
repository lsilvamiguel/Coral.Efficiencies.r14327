/*!
   \File    CsRCLikeAll05.cc
   \----------------------------
   \brief   CsRCLikeAll05 class implementation.
   \author  Paolo Schiavon
   \version 1.0
   \date    May  2004
*/


  #include <iostream>
  #include <ostream>
  #include <cstdio>

  #include <iomanip>

  #include "CsOpt.h"
  #include "CsHist.h"
  #include "CsHistograms.h"

//------------------------------
  #include "CsRCLikeAll.h"
  #include "CsRCLikeAll05.h"

  #include "CsRCParticle.h"
  #include "CsRCCluster.h"
  #include "CsRCEventPartPhotons.h"
  #include "CsRCPartPhotons.h"
  #include "CsRCPhoton.h"

  #include "CsRCHistos.h"

  #include "CsRCRecConst.h"
  #include "CsRCExeKeys.h"

  #include "CsErrLog.h"
  #include "CsRCMirrors.h"
  #include "CsRCDetectors.h"
  #include "CsRCUty.h"
//------------------------------

  using namespace std;
  using namespace CLHEP;

//===========================================================================
  CsRCLikeAll05::CsRCLikeAll05( const CsRCPartPhotons* papho ):
//-------------------------------------------------------------
    CsRCLikeAll() {

    pPartPhot_ = papho;

  }

//===========================================================================
  CsRCLikeAll05::~CsRCLikeAll05() {}
//----------------------------------


  extern "C" {
    float prob_( const float&, const int& );
    float erf_( const float& );
  }


//===========================================================================
  double CsRCLikeAll05::normSignal( const double theIpo ) {
//---------------------------------------------------------


//- compute Signal normalization for Likelihood
//  -------------------------------------------

//- Paolo  -  May 2004
//  Rev.      January 2005
//  Rev.      May 2006


    CsRCRecConst* cons = CsRCRecConst::Instance();
    static const double TwoPI = cons->TwoPI();
    static const double Drad = 360./ TwoPI;
    int sz = pPartPhot_->lPhotons().size();

    static const bool exeFrac = true;
    //static const bool exeFrac = false;
//@@--------------------------------

    CsRCHistos& hist = CsRCHistos::Ref();
    int level = hist.levelBk();
    bool doHist = false;
    if ( level >= 2 ) doHist = true;
    double xh, yh, wh;

    bool testHist = true;
//@@--------------------
    static int kOffHi = cons->kOffHi();
    static CsHist1D *hRC0011 = NULL;
    static CsHist1D *hRC0012 = NULL;
    static CsHist1D *hRC0013 = NULL;

    static bool firstCall = true;
    if( firstCall ) {
      firstCall = false;

      std::cout << "               Useful Photons ";
      if( exeFrac ) std::cout << "ON";
      else  std::cout << "OFF";
      std::cout << std::endl;

      if ( testHist ) {
        CsHistograms::SetCurrentPath("/RICH");
        stringstream Title;
        stringstream hNtest1;
        stringstream hNtest2;
        stringstream hNtest3;
        Title << "d-nPhotExp";
        hNtest1 << kOffHi + 11;
        hRC0011 = new CsHist1D( hNtest1.str(), Title.str(), 100, -50., 50. );
        Title << "nPhotExpP-K/sqP";
        hNtest2 << kOffHi + 12;
        hRC0012 = new CsHist1D( hNtest2.str(), Title.str(), 200, -0.01, 0.19);
        Title << "nPhotExpP-K/sqP";
        hNtest3 << kOffHi + 13;
        hRC0013 = new CsHist1D( hNtest3.str(), Title.str(), 100, -50., 50. );
        CsHistograms::SetCurrentPath("/");
      }
    }

    //dump_ = true;
    dump_ = false;
    //useOldCode_ = true;
    useOldCode_ = false;
//@@-------------------
    xHole_ = 100.;
    yHole_ =  60.;
    rPipe_ =  50.;
//@@-------------

    bool splitPatt = false;
    bool newPF = false;
    static int kEventC = -1;
    int kEvent = CsRichOne::Instance()->kEvent();
    static const CsRCPartPhotons* pPartPhotC = NULL;
    if( kEvent != kEventC  ||  pPartPhotC != pPartPhot_ ) {
      newPF = true;
      std::list<int> lSelCats = pPartPhot_->lSelCats();
      std::list<int>::iterator isc = lSelCats.begin();
      //for( isc=lSelCats.begin(); isc!=lSelCats.end(); isc++ )
      //std::cout << (*isc) << "  ";
      //std::cout << std::endl;
      int iCast = (*isc);
      if( iCast > 7 ) splitPatt = false;
      else {
	for( isc=lSelCats.begin(); isc!=lSelCats.end(); isc++ ) {
	  if( (*isc) > 7 ) { splitPatt = true; break; }
	}
      }
      splitPatt_ = splitPatt;
//    ----------------------
      if( dump_ )
	std::cout << "CsRCLikeAll05 : === new partPhoton " << sz << std::endl;
    }
    kEventC = kEvent;
    pPartPhotC = pPartPhot_;
    CsRCPartPhotons* papho = const_cast<CsRCPartPhotons*>(pPartPhot_);


    static float nPhotEx = 0.;
//- Compute expected number of Signal photons :
//  -------------------------------------------
    double theIrad = theIpo / 1000.;
    float pathLen = pPartPhot_->pPart()->pathLen();
    //nPhotEx = cons->nZero() * pathLen/10. * theIrad*theIrad;
    nPhotEx = papho->nPhotExpct( theIpo );

    double nPhotExS = 0.;
    nPhotExS = papho->nPhotExpctCat( theIpo );
//-------------------------------------------
    if( testHist ) {
      if( pPartPhot_->likeONLY() ) {
	xh = nPhotExS - nPhotEx;
	if( hRC0011 ) hRC0011->Fill( xh );
//      ---------------------------------
	xh = nPhotExS - sz;
	if( hRC0013 ) hRC0013->Fill( xh );
//      ---------------------------------
	static double nPhotExP;
	if( pPartPhot_->pionONLY() ) nPhotExP = nPhotEx;
	static double nPhotExK;
	if( pPartPhot_->kaonONLY()  &&  nPhotExP > 0. ) {
	  nPhotExK = nPhotEx;
	  //std::cout << "CsRCLikeAll05 " << nPhotExP << "  " << nPhotExK
	  //<< "  " <<  sqrt( nPhotExP ) << "  "
	  //<< papho->pPart()->mom() << std::endl;
	  xh = ( nPhotExP - nPhotExK ) / sqrt( nPhotExP );
	  if( hRC0012 ) hRC0012->Fill( xh );
//        ---------------------------------
	}
      }
    }
//std::cout << "CsRCLikeAll05 " << nPhotEx << "  " << nPhotExS << std::endl;

    nPhotEx_ = nPhotEx;
//  ------------------


    static float fracUse = 0.;
//- Compute average photon survival prob. :
//  ---------------------------------------
    fracUse = 1.;
    fracUsePhiSet_ = false;
    if( exeFrac ) {
      if( dump_ ) 
	std::cout << "CsRCLikeAll05 : == exeFrac " << sz << std::endl;
      int kDetPart = pPartPhot_->kDetPart();
      Hep3Vector vPoPhot0e = pPartPhot_->pPart()->vPosEmP();
      Hep3Vector vDcPart0e = pPartPhot_->pPart()->vDirIn();
      CsRCMirrors *mirr = CsRCMirrors::Instance();
      CsRCDetectors *dets = CsRCDetectors::Instance();
      std::list<CsRCMirrorNom*> lMirrNom = mirr->lMirrNom();
      std::list<CsRCMirrorNom*>::iterator in = lMirrNom.begin();
      std::list<CsRCPhotonDet*> lPhoDet = dets->lPhoDet();
      std::list<CsRCPhotonDet*>::iterator id = lPhoDet.begin();
      //      if( id == NULL  ||  in == NULL ) return  0.;
      if( id == lPhoDet.end()  ||  in == lMirrNom.end() ) return  0.;
//--- inverse rotation
      int kDetClu = kDetPart;
      if( useOldCode_ ) {
        if( kDetClu == 0 ) id++;                  //   ???
        if( kDetClu == 1 ) in++;
      } else {
	if( kDetClu == 1 ) { in++; id++; }
      }
      Hep3Vector vPoC0( (*in)->vC0() );
      vPoC0 += mirr->vCorPoPa0();
      float RR = mirr->RRv( kDetClu );
      Hep3Vector vDcPartW = pPartPhot_->vDcPartW()[kDetClu];
      Hep3Vector vPoPaMirr0 = pPartPhot_->pPart()->vPoPaMir0()[kDetPart];
      //std::cout << vPoPaMirr0.x() << "  " << vPoPaMirr0.y() << std::endl;
      //std::cout << sqrt(pow(vPoPaMirr0.x(),2)+pow(vPoPaMirr0.y(),2))
      //          << std::endl;
      double thePhot = theIpo / 1000.;
      int nStepL = 10;
//@@-----------------
      int nStepL2 = nStepL/2;
      double dPath = pathLen / float( nStepL );
      int nStepP = 36;
//@@-----------------
      double dPhi = 360. / float( nStepP );
      dPhi_ = dPhi;
      int nTot = 0;
      int nGood = 0;
//--- loop on photon angle Phi (to Y axis) :
      for( int stepp=0; stepp<nStepP; stepp++ ) {
	double phiPhot0 = pPartPhot_->lPhotons().front()->phiA() / Drad;
	double phiPhot = 0.;
	int nTotl = 0;
	int nGoodl = 0;
//----- loop on photon position along part path :
	for( int stepl=-nStepL2; stepl<=nStepL2; stepl++ ) {
	  Hep3Vector vPoPhot0 = vPoPhot0e + float( stepl )*dPath * vDcPart0e;
          if( useOldCode_ ) phiPhot = (float( stepp ) * dPhi) / Drad;
	  else  phiPhot = phiPhot0 + (float( stepp ) * dPhi) / Drad;
          double tga, tgb;
          tga = tan( thePhot ) * cos( phiPhot );
          tgb = tan( thePhot ) * sin( phiPhot );
          Hep3Vector vDcPhotWw( tga, tgb, 1.);
          vDcPhotWw = vDcPhotWw.unit();
          double theP, phiP;
          Hep3Vector vDcPhotW = CsRCUty::Instance()->
            rotfbcc( -1., vDcPartW, vDcPhotWw, theP, phiP );
//          ----------------------------------------------
          HepVector vrot( 3, 0 );
          vrot[0] = vDcPhotW.x();
          vrot[1] = vDcPhotW.y();
          vrot[2] = vDcPhotW.z();
          HepVector vPhot = (*id)->rotMatrix() * vrot;
          Hep3Vector vDcPhot0( vPhot[0], vPhot[1], vPhot[2] );
          Hep3Vector vPhotMirr0 =
            mirr->vImpMir( vPoPhot0, vDcPhot0, vPoC0, RR );
//          ----------------------------------------------

          bool Hole = false;
          double elips = pow( vPhotMirr0.x()/xHole_, 2 ) +
                         pow( vPhotMirr0.y()/yHole_, 2 ) - 1.;
          if( elips < 0. ) Hole = true;
          //Hole = false;

  	  bool Pipe = false;
          double xx = vPhotMirr0.x() - vPoPhot0.x();
          double yy = vPhotMirr0.y() - vPoPhot0.y();
          double mm = yy / xx;
          double qq = sqrt( xx*xx + yy*yy );
          //double rr = fabs( (vPoPhot0.y() * xx - vPoPhot0.x() * yy) / qq );
          //double pp = vPhotMirr0.x()*vPhotMirr0.x() +
          //            vPhotMirr0.y()*vPhotMirr0.y();
          double ssp, ssm;
          if( vPoPhot0.x()*vPoPhot0.x() + vPoPhot0.y()*vPoPhot0.y()
	      < rPipe_*rPipe_ ) Pipe = true;
          else if( vPhotMirr0.x()*vPhotMirr0.x()
		   + vPhotMirr0.y()*vPhotMirr0.y()
  	           < rPipe_*rPipe_ ) Pipe = true;
          else if( fabs( (vPoPhot0.y()* xx - vPoPhot0.x()* yy)/ qq )
		   < rPipe_ ) {
  	    //ssp = vPoPhot0.x() + vPoPhot0.y() * mm;
	    //ssm = vPhotMirr0.x() + vPhotMirr0.y() * mm;
	    //if( (ssp < 0.  &&  ssm > 0.)  ||  (ssm < 0.  &&  ssp > 0.) )
  	    //Pipe = true;
	    double aa = tga*tga + tgb*tgb;
	    double bb = tga*vPoPhot0.x() + tgb*vPoPhot0.y();
	    double cc = vPoPhot0.x()*vPoPhot0.x() + vPoPhot0.y()*vPoPhot0.y()
	      - rPipe_*rPipe_;
	    double delta = bb*bb - aa*cc;
	    if( delta > 0. ) {
	      double zp1 = ( -bb + sqrt( delta )) / aa;
	      double zp2 = ( -bb - sqrt( delta )) / aa;
	      if( !(vPoPhot0.z() > zp1  &&  vPoPhot0.z() > zp2) ) Pipe = true;
	    }
	  }
          //Pipe = false;

          nTotl++;
          if( !Hole  &&  !Pipe ) nGoodl++;
          nTot++;
          if( !Hole  &&  !Pipe ) nGood++;
          //std::cout << true << Hole << " " << Pipe << "  " << rr << "  "
          //<< vPoPhot0.x() << "  " << xPhotMirr0 << "  " << mm << "  "
          //<< ssp << "  " << ssm << std::endl;
	}
	//cout << stepp << "  " << nGoodl << "  " << nTotl << endl;
	if( nTotl == 0 ) nTotl = 1;
	fracUsePhi_[stepp] = float( nGoodl ) / float( nTotl );
//      -----------------------------------------------------
	//^fracUsePhiSet_ = true;
//      ---------------------
      }
      fracUsePhiSet_ = true;
//    ---------------------
      //cout << "------  " << nGood << "  " << nTot << endl;
      if( nTot == 0 ) nTot = 1;
      fracUse = float( nGood ) / float( nTot );
//    ----------------------------------------
      //std::cout << "05All-N  " << nPhotEx << "  " << fracUse << std::endl;
    }
    fracUse_ = fracUse;
//  ------------------
    papho->setFracUse( fracUse );
    if( dump_ ) {
      static float fracUseS;
      if( newPF ) {
	fracUseS = fracUse;
	std::cout << setprecision(4);
	std::cout << "CsRCLikeAll05 : fracUseS " << fracUse << std::endl;
	std::cout << setprecision(2);
      } else {
	if( fracUse != fracUseS ) {
	  std::cout << setprecision(4);
	  std::cout << "CsRCLikeAll05 : ... fracUse " << fracUse << std::endl;
	  std::cout << setprecision(2);
	}
      }
    }


//- Histograms :
//  ------------
    if( doHist ) {
      xh = nPhotEx;
      yh = fracUse;
      if( hist.hRC1510 ) hist.hRC1510->Fill( xh, yh );
//    -----------------------------------------------
      Hep3Vector vPoPart0 = pPartPhot_->pPart()->vPosIn();
      double rr0 = sqrt( pow( vPoPart0.x(), 2 ) + pow( vPoPart0.y(), 2 ) );
      Hep3Vector vDcPart0 = pPartPhot_->pPart()->vDirIn();
      double the0 = acos( vDcPart0.z() ) * 1000.;
      xh = the0;
      yh = rr0;
      wh = fracUse;
      if( hist.hRC1514 ) hist.hRC1514->Fill( xh, yh, wh );
//    ---------------------------------------------------
      xh = the0;
      yh = rr0;
      if( hist.hRC1515 ) hist.hRC1515->Fill( xh, yh );
//    -----------------------------------------------
      if( fracUse == 0. ) {
        if( hist.hRC1624 ) hist.hRC1624->Fill( xh, yh );
//      -----------------------------------------------
      }
      if( fracUse == 1. ) {
        if( hist.hRC1625 ) hist.hRC1625->Fill( xh, yh );
//      -----------------------------------------------
      }
    }


    if( exeFrac ) nPhotEx *= fracUse;
//----------------------------------
    //std::cout << theIpo << "  " << pPartPhot_ 
    //          << "  " << fracUse << "  " << nPhotEx << std::endl;

    double normS = nPhotEx;

    if( theIpo <= 0. ) normS = 0.;

    if( doHist ) {
      xh = normS;
      if( hist.hRC1577 ) hist.hRC1577->Fill( xh );
//    -------------------------------------------
    }

    return  normS;

  }


//===========================================================================
  double CsRCLikeAll05::normBackgr( const double theIpo ) {
//---------------------------------------------------------


//- compute Background normalization for Likelihood
//  -----------------------------------------------
//  NOT used - to be CHECKED!

//- Paolo  -  May 2004
//  Rev.      May 2006
//  TO BE CHECKED!


    CsRCRecConst* cons = CsRCRecConst::Instance();
    static const double TwoPI = cons->TwoPI();
    int sz = pPartPhot_->lPhotons().size();

    static const float focLen = pPartPhot_->pPart()->pMirPart()->RR() / 2.;

    //static const bool exeInt = false;
    static const bool exeInt = true;
//@@--------------------------------

    CsRCHistos& hist = CsRCHistos::Ref();
    int level = hist.levelBk();
    bool doHist = false;
    if ( level >= 2 ) doHist = true;
    double xh, yh, wh;

    static bool firstCall = true;
    if( firstCall ) {
      firstCall = false;

      std::cout << "               Back Int ";
      if( exeInt ) std::cout << "ON";
      else  std::cout << "OFF";
    }


//- Compute Background integral :
//-------------------------------
    double normB = 0.;
    static const double theMax = 70./ 1000.;
    static const int nThe = 23;
    static const int nPhi = 24;
//@@--------------------------
    //^long double backInt = 0.;
    //^long double backIntS = 0.;
    static long double backInt = 0.;
    static long double backIntS = 0.;
    static int kEventC = -1;
    int kEvent = CsRichOne::Instance()->kEvent();
    static const CsRCPartPhotons* pPartPhotC = NULL;
    if( exeInt ) {
      if( kEvent != kEventC  ||  pPartPhotC != pPartPhot_ ) {
        if( dump_ )
	  std::cout << "CsRCLikeAll05 : exeInt " << sz << std::endl;
        backInt = 0.;
        double dThe = theMax / float( nThe );
        double dPhi = TwoPI / float( nPhi );
        //std::cout << " " << dThe << " " << dPhi << std::endl;
        std::list<CsRCPhoton*> lPhotons = pPartPhot_->lPhotons();
        CsRCCluster clu( *(lPhotons.front()->ptToClu()) );
        int kDetPart = pPartPhot_->kDetPart();
        double xPade = pPartPhot_->pPart()->vPade()[kDetPart].x();
        double yPade = pPartPhot_->pPart()->vPade()[kDetPart].y();
        //std::cout << " " << xPade << " " << yPade << std::endl;
        int iPhot = 0;
        for( int kThe=0; kThe<nThe; kThe++ ) {
          for( int kPhi=0; kPhi<nPhi; kPhi++ ) {
	    double the = kThe*dThe + dThe/2.;
  	    double phi = kPhi*dPhi + dPhi/2.;
  	    double xClu = xPade + focLen*the * cos( phi );
  	    double yClu = yPade + focLen*the * sin( phi );
	    clu.setX( xClu );
	    clu.setY( yClu );
	    double backWgt = 0.;
	    if( clu.getBackWgt( backWgt ) ) {
	      backInt += backWgt * focLen*the * focLen*dThe * dPhi;
//            ----------------------------------------------------
 	      //std::cout << "BackInt" << "  " << backWgt << "  "
	      //          << iPhot << "  "
	      //          << the << "  " << phi << "  " << xClu << "  "
	      //          << yClu <<std::endl;
	      iPhot++;
	    }
	  }
	}
        //backInt /= (TwoPI/2. * theMax*theMax);
        //std::cout << "BackInt = " << backInt << " " << iPhot << "  "
        //          << lPhotons.size() << std::endl;
	backIntS = 0.;
	if( splitPatt_ ) {
  	  std::vector<double> vySplitLLim;
	  vySplitLLim = pPartPhot_->vySplitLLim();
	  if( vySplitLLim.size() < 2 ) {
	    vySplitLLim.clear();
	    vySplitLLim.push_back( 0. );
	    vySplitLLim.push_back( 0. );
	  }
 	  int jDetPart = 1 - kDetPart;
 	  xPade = pPartPhot_->pPart()->vPade()[jDetPart].x();
	  yPade = pPartPhot_->pPart()->vPade()[jDetPart].y();
	  int sign = 1 - 2 * jDetPart;
	  int jPhot = 0;
	  for( int kThe=0; kThe<nThe; kThe++ ) {
	    for( int kPhi=0; kPhi<nPhi; kPhi++ ) {
	      double the = kThe*dThe + dThe/2.;
	      double phi = kPhi*dPhi + dPhi/2.;
	      double xClu = xPade + focLen*the * cos( phi );
	      double yClu = yPade + focLen*the * sin( phi );
	      if( (yClu - vySplitLLim[jDetPart])* sign > 0. ) {
	        //std::cout << jDetPart << "  " << yClu << "  "
	        //<< vySplitLLim[jDetPart] << std::endl;
	        clu.setX( xClu );
	        clu.setY( yClu );
	        double backWgt = 0;
	        if( clu.getBackWgt( backWgt ) ) {
		  backIntS += backWgt * focLen*the * focLen*dThe * dPhi;
//                -----------------------------------------------------
  	          //std::cout << "BackIntS" << " " << backWgt << " " << iPhot
	          //<< "  " << the << "  " << phi << "  " << xClu << "  "
	          //<< yClu <<std::endl;
		  jPhot++;
		}
	      }
	    }
	  }
	}
	//std::cout << "BackInt  " << backInt << "  " << backIntS
	//  << "  " << lPhotons.size() << "  " << splitPatt_ << std::endl;

//----- Histograms :
//      ------------
	if( doHist ) {
	  xh = backInt;
	  yh = backIntS;
	  if( hist.hRC1511 ) hist.hRC1511->Fill( xh, yh );
//        -----------------------------------------------
	  int kDetPart = pPartPhot_->pPart()->kDetPart();
	  float ddpdet = pPartPhot_->pPart()->ddPaDet()[kDetPart];
	  xh = ddpdet;
	  yh = backInt;
	  if( kDetPart == 0 ) {
	    if( hist.hRC1518 ) hist.hRC1518->Fill( xh, yh );
//          -----------------------------------------------
	  }
	  if( kDetPart == 1 ) {
	    if( hist.hRC1519 ) hist.hRC1519->Fill( xh, yh );
//          -----------------------------------------------
	  }
	}
      }
      kEventC = kEvent;
      pPartPhotC = pPartPhot_;
    
      normB = backInt + backIntS;
//    --------------------------
    }
    else  normB = 1.;
//        ----------

    normB_ = normB;
//  --------------

    if( doHist ) {
      xh = normB;
      if( hist.hRC1578 ) hist.hRC1578->Fill( xh );
//    -------------------------------------------
    }

    return  normB;

  }


//===========================================================================
  double CsRCLikeAll05::likeSignal( const CsRCPhoton* pPhot,
//----------------------------------------------------------
				    const double theIpo ) {


//- compute Signal value for Likelihood
//  -----------------------------------
//- Paolo  -  May 2004
//  Rev.      May 2006


    //if( theIpo < 50. ) return 0.;

    CsRCRecConst* cons = CsRCRecConst::Instance();
    static const double TwoPI = cons->TwoPI();
    static const double sq2pi = sqrt( TwoPI );
    static const double Drad = 360./ TwoPI;
    int sz = pPartPhot_->lPhotons().size();

    CsRCDetectors *dets = CsRCDetectors::Instance();
    //static const float pipeLen = dets->zExitWind() - dets->zEntrWind();

    CsRCHistos& hist = CsRCHistos::Ref();
    int level = hist.levelBk();
    bool doHist = false;
    if ( level >= 2 ) doHist = true;
    double xh, yh, wh;


    static const bool exeSpl = false;
//@@--------------------------------
    static const bool exeSig = true;
    //static const bool exeSig = false;
//@@-------------------------------

    static bool firstCall = true;
    if( firstCall ) {
      firstCall = false;

      std::cout << "   Signal Spl ";
      if( exeSpl ) std::cout << "ON";
      else  std::cout << "OFF";
      std::cout << "   Signal Corr ";
      if( exeSig ) std::cout << "ON";
      else  std::cout << "OFF";
      std::cout << std::endl;
    }


//- Compute Signal term :
//  ---------------------
    double thePhot = pPhot->the();
    double sigPhot = pPhot->sigmaPhoPid( pPartPhot_ );
    double signal = 0.;
    if( sigPhot > 0.) {
      double qSw = (thePhot - fabs( theIpo )) / sigPhot;
      signal = exp( - 0.5* qSw*qSw ) / (sigPhot*sq2pi);
//-----------------------------------------------------------
    }
    //signal *= 1./ (TwoPI*theIpo);             // test int
    //signal *= thePhot / (TwoPI*theIpo);
    //signal *= thePhot / theIpo;
//@@----------------------------------


    static int kDetPart = 0;
    static float nPhotEx = 0.;
    static std::vector<double> vySplitLLim;
    static std::vector<double> cosphi0;
    static const float focLen = pPartPhot_->pPart()->pMirPart()->RR() / 2.;
    int kDetClu = pPhot->kDetClu();

    static CsRCPhoton* pPhotF = NULL;
    pPhotF = pPartPhot_->lPhotons().front();
    double theIrad = theIpo / 1000.;
    if( pPhot == pPhotF ) {

      kDetPart = pPartPhot_->kDetPart();
      vySplitLLim = pPartPhot_->vySplitLLim();
      if( vySplitLLim.size() < 2 ) {
	vySplitLLim.clear();
	vySplitLLim.push_back( 0. );
	vySplitLLim.push_back( 0. );
      }

      nPhotEx = nPhotEx_;
//    ------------------

      if( splitPatt_ ) {
	cosphi0.clear();
	for( int kDet=0; kDet<2; kDet++ ) {
          //double rOnDet = theIrad * ( pPartPhot_->pPart()->vPade()[kDet] -
	  //  	  		        pPartPhot_->vPoPaMirW()[kDet] ).mag();
	  double rOnDet = (theIrad + 0.003) * focLen;
	  double yPade = pPartPhot_->pPart()->vPade()[kDet].y();
	  double dd = fabs( yPade - vySplitLLim[kDet] );
	  double cos = 1.;
	  if( rOnDet > 0. ) cos = fabs( dd / rOnDet );
	  //cout << rOnDet << "  " << cos << endl;
	  if( cos > 1. ) cos = 1.;
	  cosphi0.push_back( cos );
	}
      }
      //std::cout << "-----  " << splitPatt_ 
      //  << "  " << cosphi0[0] << "  " << cosphi0[1]
      //  << "  " << vySplitLLim[0] << "  " << vySplitLLim[1]
      //  << "  " << theIpo << std::endl;

    }


//- Split rings effect (do NOT use) :
//  ---------------------------------
    float normSpl = 1.;
    if( exeSpl ) {
      double cosphi = 0.;
      if( splitPatt_  &&  cosphi0[kDetClu] < 1. ) {
        double xc = pPhot->ptToClu()->xc();
        double yc = pPhot->ptToClu()->yc();
        double xPade = pPartPhot_->pPart()->vPade()[kDetClu].x();
        double yPade = pPartPhot_->pPart()->vPade()[kDetClu].y();
        double tanphi = 0.;
        //std::cout << "--- " << dd << "  " << rOnDet << std::endl;
        //std::cout << xc << "  " << yc << "  " 
        //  << xPade << "  " << yPade << "  " << dd << std::endl;
        int sign = 1 - 2 * kDetPart;
        if( kDetClu == kDetPart ) {
	  if( ( yc - vySplitLLim[kDetClu] )* sign < 0. ) {
	    if((yc - yPade) > 0.) tanphi = (xc - xPade) / (yc - yPade);
  	    cosphi = 1./ sqrt( 1.+ tanphi*tanphi );
	    if( cosphi > cosphi0[kDetClu] ) normSpl = cosphi0[kDetClu]/cosphi;
//                                          ---------------------------------
  	    //else normSpl = 0.;
	  }
    	  //std::cout << kDetClu << "  " << cosphi << "  " << cosphi0[kDetClu]
	  //          << "  " << normSpl << std::endl;
	}
	tanphi = 0.;
        if( kDetClu != kDetPart ) {
	  normSpl = 0.;
//        ------------
	  //if( ( yc - vySplitLLim[kDetClu] )* sign < 0. ) {
	  if( ( yc - vySplitLLim[kDetClu] )* sign > 0. ) {
	    if((yc - yPade) > 0.) tanphi = (xc - xPade) / (yc - yPade);
	    cosphi = 1./ sqrt( 1.+ tanphi*tanphi );
	    if( cosphi > cosphi0[kDetClu] )
	      normSpl = 1.- cosphi0[kDetClu]/cosphi;
//            -------------------------------------
	    else normSpl = 0.;
//               ------------
	  }
	  //std::cout << kDetClu << "  " << cosphi << "  " << cosphi0[kDetClu]
	  //          << "  " << normSpl << std::endl;
	}
      }
    }


//- Photon survival prob. :
//  -----------------------
    float probSig = 1.;
    if( exeSig ) {
      if( useOldCode_ ) fracUsePhiSet_ = false;
      if( fracUsePhiSet_ ) {
	if( dump_ )
	  if( pPhot == pPartPhot_->lPhotons().front() )
	    std::cout << "CsRCLikeAll05 : (Sig) " << sz << std::endl;
	double phiPhot0 = pPartPhot_->lPhotons().front()->phiA();
	double phiPhot = pPhot->phiA();
	double diffPhi = phiPhot - phiPhot0;
	if( diffPhi < 0. ) diffPhi = diffPhi + 360.;
	if( diffPhi > 360. ) diffPhi = diffPhi - 360.;
	//std::cout << "05All-S  " << phiPhot << "  " << phiPhot0 << "  "
	//  << diffPhi << "  " << dPhi_ << std::endl;
	int kFrac = int( diffPhi / dPhi_ );
	if( kFrac < 0 ) kFrac = 0;
	if( kFrac > 35 ) kFrac = 35;
	probSig = fracUsePhi_[kFrac];
//      ----------------------------
        //std::cout << "05All-S  " << nPhotEx << "  " << probSig << std::endl;
      } else {
        if( dump_ )
	  if( pPhot == pPartPhot_->lPhotons().front() )
	    std::cout << "CsRCLikeAll05 : exeSig " << sz << std::endl;
        CsRCMirrors *mirr = CsRCMirrors::Instance();
        float RR = mirr->RRv( kDetClu );
        Hep3Vector vPoPhot0e = pPartPhot_->pPart()->vPosEmP();
        Hep3Vector vDcPart0e = pPartPhot_->pPart()->vDirEmP();
        double phiPhot = pPhot->phiA() / Drad;
        thePhot = theIpo / 1000.;
//------------------------------
        list<CsRCMirrorNom*> lMirrNom = mirr->lMirrNom();
        list<CsRCMirrorNom*>::iterator in = lMirrNom.begin();
        list<CsRCPhotonDet*> lPhoDet = dets->lPhoDet();
        list<CsRCPhotonDet*>::iterator id = lPhoDet.begin();
	//        if( id == NULL  ||  in == NULL ) return  0.;
        if( id == lPhoDet.end()  ||  in == lMirrNom.end() ) return  0.;
//----- inverse rotation
	if( useOldCode_ ) {
          if( kDetClu == 0 ) id++;         //   ???
          if( kDetClu == 1 ) in++;
        } else  if( kDetClu == 1 ) { in++; id++; }
        Hep3Vector vPoC0( (*in)->vC0() );
        vPoC0 += mirr->vCorPoPa0();
        Hep3Vector vDcPartW = pPartPhot_->vDcPartW()[kDetClu];
        double tga, tgb;
        tga = tan( thePhot ) * cos( phiPhot );
        tgb = tan( thePhot ) * sin( phiPhot );
        Hep3Vector vDcPhotWw( tga, tgb, 1.);
        vDcPhotWw = vDcPhotWw.unit();
        double theP, phiP;
        Hep3Vector vDcPhotW = CsRCUty::Instance()->
          rotfbcc( -1., vDcPartW, vDcPhotWw, theP, phiP );
//        -----------------------------------------------
        HepVector vrot( 3, 0 );
        vrot[0] = vDcPhotW.x();
        vrot[1] = vDcPhotW.y();
        vrot[2] = vDcPhotW.z();
        HepVector vPhot = (*id)->rotMatrix() * vrot;
        Hep3Vector vDcPhot0( vPhot[0], vPhot[1], vPhot[2] );
        //std::cout << vDcPhotWw << " " << vDcPhotW << " " << vDcPhot0
        //          <<std::endl;
        double pathLen = pPartPhot_->pPart()->pathLen();
        int nStep = 10;
        double dPath = pathLen / float( nStep );
//----- in MRS
        int nTot = 0;
        int nGood = 0;
//----- loop on photon position along part path :
        for( int step=-nStep/2; step<=nStep/2; step++ ) {
          Hep3Vector vPoPhot0 = vPoPhot0e + float( step )*dPath * vDcPart0e;
          //std::cout << vPoPhot0e << "  " << vPoPhot0 << std::endl;

          Hep3Vector vPhotMirr0 =
	    mirr->vImpMir( vPoPhot0, vDcPhot0, vPoC0, RR );
//          ----------------------------------------------
	  bool Hole = false;
          double elips = pow( vPhotMirr0.x()/xHole_, 2 ) +
                         pow( vPhotMirr0.y()/yHole_, 2 ) - 1.;
          if( elips < 0. ) Hole = true;
          //Hole = false;
          bool Pipe = false;
          double xx = vPhotMirr0.x() - vPoPhot0.x();
	  if( xx == 0.) continue;
          double yy = vPhotMirr0.y() - vPoPhot0.y();
          double mm = yy / xx;
          double qq = sqrt( xx*xx + yy*yy );
          //double rr = fabs( (vPoPhot0.y() * xx - vPoPhot0.x() * yy) / qq );
          //double pp = vPhotMirr0.x()*vPhotMirr0.x() +
          //            vPhotMirr0.y()*vPhotMirr0.y();
          double ssp, ssm;
          if( vPoPhot0.x()*vPoPhot0.x() + vPoPhot0.y()*vPoPhot0.y()
	      < rPipe_*rPipe_ ) Pipe = true;
          else if( vPhotMirr0.x()*vPhotMirr0.x()
		   + vPhotMirr0.y()*vPhotMirr0.y()
		   < rPipe_*rPipe_ ) Pipe = true;
          else if( fabs( (vPoPhot0.y()* xx - vPoPhot0.x()* yy)/ qq )
		   < rPipe_) {
	    //ssp = vPoPhot0.x() + vPoPhot0.y() * mm;
  	    //ssm = vPhotMirr0.x() + vPhotMirr0.y() * mm;
	    //if( (ssp < 0.  &&  ssm > 0.)  ||  (ssm < 0.  &&  ssp > 0.) )
	    //Pipe = true;
	    double aa = tga*tga + tgb*tgb;
	    double bb = tga*vPoPhot0.x() + tgb*vPoPhot0.y();
	    double cc = vPoPhot0.x()*vPoPhot0.x() + vPoPhot0.y()*vPoPhot0.y()
	      - rPipe_*rPipe_;
	    double delta = bb*bb - aa*cc;
	    if( delta > 0. ) {
	      double zp1 = ( -bb + sqrt( delta )) / aa;
	      double zp2 = ( -bb - sqrt( delta )) / aa;
	      if( !(vPoPhot0.z() > zp1  &&  vPoPhot0.z() > zp2) ) Pipe = true;
	    }
	  }
          //Pipe = false;

          nTot++;
          if( !Hole  &&  !Pipe ) nGood++;
          //std::cout << Hole << "  " << Pipe << std::endl;
	}
	if( nTot == 0 ) nTot = 1;
        probSig = float( nGood ) / float( nTot );
//      ----------------------------------------
        //std::cout << "05All-S  " << nPhotEx << "  " << probSig << std::endl;
        //std::cout << nPhotEx << "  " << normSpl << "  " << probSig
        //          << std::endl;
        //std::cout << std::endl;
      }
    }

    signal *= nPhotEx;
//-------------------
    if( exeSpl ) signal *= normSpl;
//--------------------------------
    if( exeSig ) signal *= probSig;
//--------------------------------

//- Histogramms :
//  -------------
    if( doHist ) {
      xh = nPhotEx;
      yh = probSig;
      if( hist.hRC1512 ) hist.hRC1512->Fill( xh, yh );
//    -----------------------------------------------
      Hep3Vector vPoPart0 = pPartPhot_->pPart()->vPosIn();
      double rr0 = sqrt( pow( vPoPart0.x(), 2 ) + pow( vPoPart0.y(), 2 ) );
      Hep3Vector vDcPart0 = pPartPhot_->pPart()->vDirIn();
      double the0 = acos( vDcPart0.z() ) * 1000.;
      xh = the0;
      yh = rr0;
      wh = probSig;
      if( hist.hRC1516 ) hist.hRC1516->Fill( xh, yh, wh );
//    ---------------------------------------------------
      xh = the0;
      yh = rr0;
      if( hist.hRC1517 ) hist.hRC1517->Fill( xh, yh );
//    -----------------------------------------------
      if( probSig == 0. ) {
        if( hist.hRC1626 ) hist.hRC1626->Fill( xh, yh );
//      -----------------------------------------------
      }
      if( probSig == 1. ) {
        if( hist.hRC1627 ) hist.hRC1627->Fill( xh, yh );
//      -----------------------------------------------
      }
      xh = normSpl;
      yh = -1.;
      if( kDetPart == 0 ) {
        if( kDetClu == kDetPart ) yh = 0.5;
        if( kDetClu != kDetPart ) yh = 1.5;
      }
      if( kDetPart == 1 ) {
        if( kDetClu == kDetPart ) yh = 2.5;
        if( kDetClu != kDetPart ) yh = 3.5;
      }
      if( hist.hRC1513 ) hist.hRC1513->Fill( xh, yh );
//    -----------------------------------------------
      if( cosphi0.size() >= 2 ) {
	xh = cosphi0[kDetPart];
	yh = 4.5;
	if( hist.hRC1513 ) hist.hRC1513->Fill( xh, yh );
//      -----------------------------------------------
      }
      if( pPartPhot_->likeONLY() ) {
	xh = signal;
	if( pPhot->isPMT() ) if( hist.hRC1574 ) hist.hRC1574->Fill( xh );
//                           -------------------------------------------
	if( pPhot->isAPV() ) if( hist.hRC1575 ) hist.hRC1575->Fill( xh );
//                           -------------------------------------------
      }
    }


    if( theIpo <= 0. ) signal = 0.;

    return  signal;

  }


//===========================================================================
  double CsRCLikeAll05::likeBackgr( const CsRCPhoton* pPhot,
//----------------------------------------------------------
				    const double theIpo ) {


//- compute Background value for Likelihood
//  ---------------------------------------
//- Paolo  -  May 2004
//  Rev.      May 2006


    //if( theIpo < 50. ) return 0.;

    static const float focLen = pPartPhot_->pPart()->pMirPart()->RR() / 2.;
    static const float factor = focLen/1000. * focLen/1000.;
    int sz = pPartPhot_->lPhotons().size();

    static bool firstCall = true;
    if( firstCall ) {
      firstCall = false;
    }

    CsRCHistos& hist = CsRCHistos::Ref();
    int level = hist.levelBk();
    bool doHist = false;
    if ( level >= 2 ) doHist = true;
    double xh, yh, wh;

//- Warning : kPhot is NOT sequential!
    int kPhot = pPhot->kPhot();
    double thePhot = pPhot->the();

//- from   back-para-... .new-vector74   background map
//         ---------------------------
    //static const double backWgtMin = 0.0002;
    static const double backWgtMin = 0.00001;
//@@----------------------------------------

    static const CsRCPartPhotons* pPartPhotP = NULL;
    if( kPhot >= 1000 ) {
      if( pPartPhotP != pPartPhot_ ) {
        std::cout << " RICHONE, CsRCLikeAll05::likeBackgr : "
	    	  << " BUFFER Overflow!  " << kPhot << "  " 
		  << pPartPhot_->lPhotons().size() << std::endl;
      }
      pPartPhotP = pPartPhot_;
      return  backWgtMin;
//    ------------------
    }

    double backgr = 0.;
    double backWgt = 0.;
    static int kEventC = -1;
    static CsRCPartPhotons* pPartPhotC = NULL;
    static CsRCPhoton* pPhotC = NULL;
    static int kPhotC = -1;
    int kEvent = CsRichOne::Instance()->kEvent();
//- Warning : pPartPhot can be the same for two conseq. events!
//            also pPhot suspicious
    if( kEvent != kEventC  ||  pPartPhot_ != pPartPhotC ) {
      phoBackSet_ = false;
      pPhotC = const_cast<CsRCPhoton*>(pPhot);
      kPhotC = kPhot;
      //std::cout << "CsRCLikeAll05 : --- zeroPF " << std::endl;
    } else {
      if( kPhot == kPhotC ) {
	phoBackSet_ = true;
	//std::cout << "CsRCLikeAll05 : zeroF " << std::endl;
      }
    }
    if( useOldCode_ ) phoBackSet_ = false;

    if( phoBackSet_ ) {
      if( dump_ ) if( kPhot == kPhotC )
	std::cout << "CsRCLikeAll05 : (Back) " << sz << std::endl;
      backgr = phoBack_[kPhot];
//    ------------------------
    } else {
      if( dump_ ) if( kPhot == kPhotC )
	std::cout << "CsRCLikeAll05 : = exeBack " << sz << std::endl;
      backgr = backWgtMin;
      if( pPhot->ptToClu()->getBackWgt( backWgt ) ) {
//        ---------------------------------------
        if( backWgt < backWgtMin ) backWgt = backWgtMin;
        backgr = factor * backWgt * thePhot;             // (3300/1000)**2
        //backgr = factor * backWgt * thePhot/1000.;   
//      -----------------------------------
	if( !useOldCode_ ) phoBack_[kPhot] = backgr;
      }
      else {
	backWgt = backWgtMin;
        backgr = factor * backWgt * thePhot;
	if( !useOldCode_ ) phoBack_[kPhot] = backgr;
	std::cout << " RICHONE, CsRCLikeAll05::likeBackgr() : "
		  << "MINIMUM background assignment! - Ev. "
		  << kEvent << "  Phot. " << kPhot << std::endl;
      }
    }
    kEventC = kEvent;
    pPartPhotC = const_cast<CsRCPartPhotons*>(pPartPhot_);

    //if( backgr == 0.) {
    //if( pPartPhot_->likeONLY() ) {
    //std::cout << setprecision( 12 ) << "LAll05  ";
    //std::cout << kPhot << "  "
    //	  << log(backgr) << "  " << backWgt << "  " << thePhot << std::endl;
    //}
    //}

    float normSpl = 1.;
//  ------------------

//- Histograms :
//  ------------
    if( doHist ) {
      int kDetPart = pPartPhot_->kDetPart();
      int kDetClu = pPhot->kDetClu();
      xh = normSpl;
      yh = 10.5;
      if( kDetPart == 0 ) {
        if( kDetClu == kDetPart ) yh = 5.5;
        if( kDetClu != kDetPart ) yh = 6.5;
      }
      if( kDetPart == 1 ) {
        if( kDetClu == kDetPart ) yh = 7.5;
        if( kDetClu != kDetPart ) yh = 8.5;
      }
      if( hist.hRC1513 ) hist.hRC1513->Fill( xh, yh );
//    -----------------------------------------------

      if( backgr > 0. ) {
        xh = pPhot->ptToClu()->xc();
        yh = pPhot->ptToClu()->yc();
        wh = backgr;
        if( hist.hRC1530 ) hist.hRC1530->Fill( xh, yh, wh );
//      ---------------------------------------------------
        if( hist.hRC1540 ) hist.hRC1540->Fill( xh, yh );
//      -----------------------------------------------
      }
      if( pPartPhot_->likeONLY() ) {
	xh = backgr;
	if( pPhot->isPMT() ) if( hist.hRC1576 ) hist.hRC1576->Fill( xh );
//                           -------------------------------------------
	if( pPhot->isAPV() ) if( hist.hRC1579 ) hist.hRC1579->Fill( xh );
//                           -------------------------------------------
      }
    }

    return  backgr;

  }


//===========================================================================
  double CsRCLikeAll05::getRingBackground( const double theReco ) {
//-----------------------------------------------------------------


//-  Paolo - May 2004

//-  NOT IMPLEMENTED !

    CsRCRecConst * cons = CsRCRecConst::Instance();

    static bool firstCall = true;
    if( firstCall ) {
      firstCall = false;
    }

    double corrBack = 0.;

    return  corrBack;

  }
