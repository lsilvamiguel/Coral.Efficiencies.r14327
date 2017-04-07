
/*!
   \file    CsRCMirrors.cc
   \----------------------
   \brief   CsRCMirrors class implementation.
   \author  Paolo Schiavon
   \version 0.01
   \date    11 August 2000
*/


//---------------------------
#include <iostream>
#include <fstream>
#include <iomanip>
#include  <sstream>
#include <string>
//----------------------------
//----------------------------
  #include "CsInit.h"
  #include "CsErrLog.h"
  #include "CsOpt.h"

  #include "CsRCMirrors.h"

  #include "CsRCCluster.h"

  #include "CsRCEventParticles.h"

  #include "CsRCPartPhotons.h"
  #include "CsRCPhoton.h"

  #include "CsRCEventRings.h"
  #include "CsRCRing.h"

  #include "CsRCEventAnalysis.h"

  #include "CsRCChiSqFit.h"
  #include "CsRCCircleFit.h"

  #include "CsRCRecConst.h"
  #include "CsRCExeKeys.h"
  #include "CsRCHistos.h"
  #include "CsRCUty.h"
//----------------------------

#  include <ostream>
#  include <cstdio>

  using namespace std;
  using namespace CLHEP;

  CsRCMirrors* CsRCMirrors::instance_ = 0;

//===========================================================================
  CsRCMirrors* CsRCMirrors::Instance() {
//--------------------------------------
    if( instance_ == 0 ) instance_ = new CsRCMirrors();
    return instance_;
  }

//===========================================================================
  CsRCMirrors::CsRCMirrors() {
//----------------------------
    lMirrNom_.clear();
    lMirrEle_.clear();
    lMirrEleAlg_.clear();
    lMirrEleK_.clear();

    nMirEleSel_ = 0;

    vCorPoPa0_.set(0., 0., 0.);
    corPaRR_ = 0.; 
    vCorPoPho0_.set(0., 0., 0.);
    corPhoRR_ = 0.;

    hRC1980 = NULL;
    hRC3980 = NULL;

    hRC7450 = NULL;
    hRC7451 = NULL;
    hRC7650 = NULL;
    hRC7651 = NULL;
    hRC7850 = NULL;
  }


//===========================================================================
  void CsRCMirrors::setMirrNom( const string name,
//------------------------------------------------
                                const double xC0, const double yC0,
				const double zC0, const double RR ) {

    int kMir = lMirrNom_.size();
    lMirrNom_.push_back( new CsRCMirrorNom( kMir, name, xC0, yC0, zC0,  RR ));

  }


//===========================================================================
  void CsRCMirrors::setMirrEle( const string name,
//------------------------------------------------
                                const double theta, const double phi,
				const double RR,
                                const double deTheta, const double dePhi,
                                const double delta, const double qfact,
				const int align ) {

    int kMir = lMirrEle_.size();
    lMirrEle_.push_back( new CsRCMirrorElem( kMir, name, theta, phi, RR,
                              deTheta, dePhi, delta, qfact, align ) );

  }

//===========================================================================
  CsRCMirrors::CsRCMirrors( const CsRCMirrors &mirr ) {
//-----------------------------------------------------
    cout << "RICHONE : CsRCMirrors CopyConstructor" << endl;
    instance_ = mirr.instance_;

    lMirrNom_ = mirr.lMirrNom_;
    lMirrEle_ = mirr.lMirrEle_;
    lMirrEleAlg_ = mirr.lMirrEleAlg_;
    lMirrEleK_ = mirr.lMirrEleK_;

    nMirEleSel_ = mirr.nMirEleSel_;

    pMirrElePa_ = mirr.pMirrElePa_;
    vCorPoPa0_ = mirr.vCorPoPa0_;
    corPaRR_ = mirr.corPaRR_;
    vCorPoPho0_ = mirr.vCorPoPho0_;
    corPhoRR_ = mirr.corPhoRR_;

    hRC1980 = mirr.hRC1980;
    hRC3980 = mirr.hRC3980;
    vRC1800 = mirr.vRC1800;
    vRC3800 = mirr.vRC3800;

    vRC7100 = mirr.vRC7100;
    vRC7300 = mirr.vRC7300;
    hRC7450 = mirr.hRC7450;
    hRC7451 = mirr.hRC7451;
    vRC7500 = mirr.vRC7500;
    hRC7650 = mirr.hRC7650;
    hRC7651 = mirr.hRC7651;
    vRC7700 = mirr.vRC7700;
    hRC7850 = mirr.hRC7850;

  }

//===========================================================================
  void CsRCMirrors::printNomAll() const {
//---------------------------------------
    cout << endl;
    cout << "CsRCMirrors::printNomAll() : Nominal Mirror Geometry" << endl
         << "----------------------------------------------------" << endl;
    list<CsRCMirrorNom*>::const_iterator in;
    for( in=lMirrNom_.begin(); in!=lMirrNom_.end(); in++ ) {
      cout << "Mirror  " << (*in)->name()
	   << "   centre  " << (*in)->vC0()
           << ",  Radius  " << (*in)->RR() << endl;
    }
  }

//===========================================================================
  void CsRCMirrors::printEleAll() const {
//---------------------------------------
    double rg = 180. / 3.1415926;
    cout << endl;
    cout << "CsRCMirrors::printEleAll() : Mirror Element Geometry" << endl
         << "----------------------------------------------------" << endl;
    list<CsRCMirrorElem*>::const_iterator ie;
    for( ie=lMirrEle_.begin(); ie!=lMirrEle_.end(); ie++ ) {
      cout << "Mirror  " << (*ie)->name()
           << " :  theta= " << (*ie)->theta()*rg
           << ",  phi= " << (*ie)->phi()*rg
           << ",  RR= " << (*ie)->RR()
           << ",  dTheta= " << (*ie)->deTheta()*rg
           << ",  dPhi= " << (*ie)->dePhi()*rg
           << ",  delta= " << (*ie)->delta()
           << ",  quality =" << (*ie)->qfact()
           << ",  align= " << (*ie)->align()
           << ",  --- pos = " << (*ie)->vpos()
	   << endl;
    }
  }

//===========================================================================
  CsRCMirrors::~CsRCMirrors() {
//-----------------------------
    lMirrNom_.clear();
    lMirrEle_.clear();
    lMirrEleAlg_.clear();
    lMirrEleK_.clear();
    vRC1800.clear();
    vRC3800.clear();
    vRC7100.clear();
    vRC7300.clear();
    vRC7500.clear();
    vRC7700.clear();
  }


//=========================================================================== 
  bool CsRCMirrors::doSelMirrors( CsRCParticle* part ) {
//------------------------------------------------------


//--- version 0.02   july 2000
//    ------------------------
//    corrects nominal values of mirror element position and
//    radius according to the particle impact on mirrror surface.


//--- from "CsRCRecConst"
      CsRCRecConst *cons = CsRCRecConst::Instance();
//--- from "CsRCExeKeys"
      CsRCExeKeys *key = CsRCExeKeys::Instance();

//--- from "CsRCHistos"
      CsRCHistos& hist = CsRCHistos::Ref();
      double xh, yh;

      int k;


//--- particle variables :
//    --------------------
      Hep3Vector vPoPart0 = part->vPosIn();
      Hep3Vector vDcPart0 = part->vDirIn();
      //cout << "vPoPart0   " << vPoPart0;
      //cout << "vDcPart0   " << vDcPart0 << endl;

//--- compute variables in both MWR :
//    -------------------------------
//--- LOOP over up and down MIRROR :
//    ------------------------------
      int kMirPart = -1;
      CsRCMirrorNom* pMirPart = NULL;
      vector<Hep3Vector> vPoPaMir0w;
      list<CsRCMirrorNom*>::iterator in;
      int kMir = 0;
      for( in=lMirrNom_.begin(); in!=lMirrNom_.end(); in++ ) {
	Hep3Vector vPoC0 = (*in)->vC0();
	double RR = (*in)->RR();
	//cout << "doSelMirrors " << kMir << "  " << RR << endl;

//----- particle impact on mirrors (MRS) :
//      ----------------------------------
        Hep3Vector vPoPaMir0 = vImpMir( vPoPart0, vDcPart0, vPoC0, RR );
//                             ----------------------------------------
        //(*in)->setPoPaMir0( vPoPaMir0 );
	//cout << "vPoPaMir0   " << vPoPaMir0 << "  " << RR << endl;

        if ( kMir == 0  &&  vPoPaMir0.y() >= 0. ) {
	  kMirPart = kMir;
	  part->setDetPart( kMirPart );
          pMirPart = (*in); 
	  part->setpMirPart( pMirPart );
	}
        if ( kMir == 1  &&  vPoPaMir0.y() < 0. ) {
          kMirPart = kMir;
	  part->setDetPart( kMirPart );
	  pMirPart = (*in);
	  part->setpMirPart( pMirPart );
	}

	vPoPaMir0w.push_back( vPoPaMir0 );
	kMir++;
      }   /* end of loop on mirrors: kMir */

      part->setPoPaMir0( vPoPaMir0w );

      if ( pMirPart == NULL )  {
        string str = 
	  "RICHONE, CsRCMirrors::doSelMirrors() : wrong particle reflection";
	CsErrLog::Instance()->mes( elError, str );

	//cout << vPoPart0 << "  " << vDcPart0 << "  " << vPoPaMir0w.front()
	//     << "  " << vPoPaMir0w.back() << "  "  << kMirPart << "  "
	//     << part->mom() << endl;

	return false;
      }

      xh = part->vPoPaMir0()[kMirPart].x();
      yh = part->vPoPaMir0()[kMirPart].y();
      if( hist.hRC3505 ) hist.hRC3505->Fill( xh, yh );
//hh                     ----------------------------

//--- corrections for mirror element differences :
//    --------------------------------------------
      pMirrElePa_ = NULL;
      vCorPoPa0_.set(0., 0., 0.);
      corPaRR_ = 0.;
      vCorPoPho0_.set(0., 0., 0.);
      corPhoRR_ = 0.;
      if( key->CorrMirror() ) {

//--- mirror elements hit by the particle photons :
//    ---------------------------------------------
      Hep3Vector vPoPaMir = part->vPoPaMir0()[kMirPart];
      double ddEleQ = 0;
      list<CsRCMirrorElem*>::iterator imirr;
      for( imirr=lMirrEle_.begin(); imirr!=lMirrEle_.end(); imirr++ ) {
        //Hep3Vector vCePos = (*imirr)->vpos() + (*imirr)->mirNo()->vC0();
        Hep3Vector vCePos = (*imirr)->vpos();   //   030122
        ddEleQ = (vPoPaMir - vCePos ).mag2();
        (*imirr)->setddEleQ( ddEleQ );
      }

      nMirEleSel_ = 3;                          // provisional !
      list<CsRCMirrorElem*>::iterator imirrMin[nMirEleSel_];
      list<CsRCMirrorElem*>::iterator kmirrmn = lMirrEle_.begin();
      char mirpow = ' ';
      if( kMirPart == 0 ) mirpow = 'U';
      if( kMirPart == 1 ) mirpow = 'D';
      bool doCorr = false;
      int kEle;
      for( kEle=0; kEle < nMirEleSel_; kEle++) {
        double ddElemn = 1.e+09;
        for( imirr=lMirrEle_.begin(); imirr!=lMirrEle_.end(); imirr++ ) {
          if( (*imirr)->mirpo() == mirpow ) {
            bool comp = true;
            int jEle;
            for( jEle=0; jEle < kEle; jEle++) {
	      if( imirr == imirrMin[jEle] ) { comp = false; break; }
            }
            if( !comp ) continue;
            ddEleQ = (*imirr)->ddEleQ();
            if( ddEleQ < ddElemn ) { 
              ddElemn = ddEleQ;
              kmirrmn = imirr;
            }
            doCorr  = true;
          }
        }
        imirrMin[kEle] = kmirrmn;
      }
      //cout << (*imirrMin[0])->name() << "  " 
      //     << (*imirrMin[1])->name() << "  "
      //     << (*imirrMin[2])->name() << endl;

//--- correction for the mirror element hit by the particle (MAIN RS) :
//    -----------------------------------------------------------------
      if( doCorr ) {
        kmirrmn = imirrMin[0];
	pMirrElePa_ = (*kmirrmn);                            //   030911
        //vCorPoPa0_ = (*kmirrmn)->vC0();
        //corPaRR_ = (*kmirrmn)->RR() - pMirPart->RR();
        vCorPoPa0_ = (*kmirrmn)->vC0() - pMirPart->vC0();    //   030120
        corPaRR_ = (*kmirrmn)->RR() - pMirPart->RR();        //   030120

	part->setVCorPoPa0( vCorPoPa0_ );                    //   101220
	part->setCorPaRR( corPaRR_ );                        //   101220
      }
      //cout << (*kmirrmn)->name() << "  " << vCorPoPa0_
      //     << "  " << (*kmirrmn)->RR() << "  " <<  pMirPart->RR()
      //     << "  " << corPaRR_ << endl;

//--- average correc. for the mirror elements hit by the photons (MAIN RS) :
//    ---------------------------------------------------------------------
      if( doCorr ) {
        Hep3Vector poPho0av(0., 0., 0.);
        double RRav  = 0.;
        double norm = 0;
        for( kEle=0; kEle < nMirEleSel_; kEle++) {
          double ddEleQ = (*imirrMin[kEle])->ddEleQ();
          ddEleQ = 1./ddEleQ;
          poPho0av += (*imirrMin[kEle])->vC0() * ddEleQ;
          RRav += (*imirrMin[kEle])->RR() * ddEleQ;
          norm += ddEleQ;
        }
        if( norm != 0. ) {
          norm = 1./norm;
          poPho0av = poPho0av * norm;
          RRav *= norm;
        }
        //vCorPoPho0_ = poPho0av;
        //corPhoRR_ = RRav - pMirPart->RR();
        vCorPoPho0_ = vCorPoPa0_;                //   030120
        corPhoRR_ = corPaRR_;                    //   030120
      }
      //cout << (*kmirrmn)->name() << "  " << vCorPoPho0_ 
      //     << "  " << corPhoRR_ << endl;
      }

      return true;

  }

/*
//============================================================================
  void CsRCMirrors::doAliMirrors( CsRCPartPhotons* ipf, CsRCRing* ring ) {
//------------------------------------------------------------------------

//--- version 0.02,  July  2000
//    -------------------------
//    OBSOLETE !

//--- from "CsRCRecConst"
      CsRCRecConst *cons = CsRCRecConst::Instance();
//--- from "CsRCExeKeys"
      CsRCExeKeys *key = CsRCExeKeys::Instance();

      CsOpt* opt = CsOpt::Instance();

//--- from "CsRCHistos"
      CsRCHistos& hist = CsRCHistos::Ref();
      double xh, yh, wh;

      static bool firstCall = true;
      if( firstCall ) {
        firstCall = false;
        key->acknoMethod( "CsRCMirrors::MirrorAlignment" );

        list<string> lName;
        bool boo;
        boo = opt->CsOpt::getOpt( "RICHONE", "MirrorElemName", lName );
        if( boo && !lName.empty() ) {
          list<CsRCMirrorElem*>::iterator imirr;
          list<string>::iterator is = lName.begin();
          if( (*is) == "ALL" ) {
            for( imirr=lMirrEle_.begin(); imirr!=lMirrEle_.end(); imirr++ ) {
              (*imirr)->setAlign( true );
	    }
	  }
          else if( (*is) == "UP" ) {
            for( imirr=lMirrEle_.begin(); imirr!=lMirrEle_.end(); imirr++ ) {
              if( (*imirr)->mirpo() == 'U' ) { (*imirr)->setAlign( true ); }
	    }
	  }
          else if( (*is) == "DOWN" ) {
            for( imirr=lMirrEle_.begin(); imirr!=lMirrEle_.end(); imirr++ ) {
              if( (*imirr)->mirpo() == 'D' ) { (*imirr)->setAlign( true ); }
	    }
	  } 
          else {
            for( is=lName.begin(); is!=lName.end(); is++ ) {
              for( imirr=lMirrEle_.begin(); imirr!=lMirrEle_.end(); imirr++ ) {
                if( (*is) == (*imirr)->name() ) { (*imirr)->setAlign( true ); }
	      }
	    }
	  }
          if( lName.front() == "NONE" ) {
            for( imirr=lMirrEle_.begin(); imirr!=lMirrEle_.end(); imirr++ ) {
              (*imirr)->setAlign( false );
	    }
	  }
        }

        int kElm  = 0;
        list<CsRCMirrorNom*>::iterator jmirr = 0;
        list<CsRCMirrorElem*>::iterator imirr = 0;
        for( imirr=lMirrEle_.begin(); imirr!=lMirrEle_.end(); imirr++ ) {
          if( (*imirr)->align() ) { 
            lMirrEleAlg_.push_back( (*imirr) );
            lMirrEleK_.push_back( kElm );
	    //cout << (*imirr)->name() << endl;
	  }

          char mirpo = (*imirr)->mirpo();
          if( mirpo == 'U' ) { jmirr = lMirrNom_.begin(); }
          if( mirpo == 'D' ) { jmirr = lMirrNom_.begin(); jmirr++; }
          //Hep3Vector vCePos = (*imirr)->vpos() + (*jmirr)->vC0();
	  Hep3Vector vCePos = (*imirr)->vpos();   //   030122
	  //cout << (*imirr)->name() << vCePos << endl;
          xh = 6.*float( kElm ) + 0.5;
          wh = vCePos.x();
          if( hist.hRC3750 ) hist.hRC3750->Fill( xh, wh );
//hh                         ----------------------------
          xh = 6.*float( kElm ) + 1.5;
          wh = vCePos.y();
          if( hist.hRC3750 ) hist.hRC3750->Fill( xh, wh );
//hh                         ----------------------------
          xh = 6.*float( kElm ) + 2.5;
          wh = vCePos.z();
          if( hist.hRC3750 ) hist.hRC3750->Fill( xh, wh );
//hh                         ----------------------------
          kElm++;
          xh = vCePos.x();
          yh = vCePos.y();
          if( hist.hRC3751 ) hist.hRC3751->Fill( xh, yh );
//hh                         ----------------------------
        }

        CsHistograms::SetCurrentPath("/RICH");
//@@-----------------------------------------
        static int kOffHi = cons->kOffHi();
        int kmirr = 0;
        list<CsRCMirrorElem*>::iterator amirr;
        for( amirr=lMirrEleAlg_.begin(); amirr!=lMirrEleAlg_.end(); amirr++ )
        {
          string hTitle;
          stringstream hN1800;
          hN1800 << kOffHi + 1800 + kmirr;
          hTitle = "phi photon vs theta photon, align. - in";
          vRC1800.push_back( new CsHist2D( hN1800.str(), hTitle,
       			                   140, 0., 70., 90, 0., 360. ) );
	  hN1800.clear();
	  stringstream nTitle;
          stringstream hN3800;
	  nTitle << (*amirr)->name() << " - fit ";
          hN3800 << kOffHi + 3800 + kmirr;
          //hTitle = "phi photon vs theta photon, align. - fit";
          //vRC3800.push_back( new CsHist2D( hN3800.str(), hTitle,
          vRC3800.push_back( new CsHist2D( hN3800.str(), nTitle.str(),
	     				   140, -70., 70., 140, -70., 70. ) );
	  nTitle.clear();
	  hN3800.clear();
          kmirr++;
        }
        CsHistograms::SetCurrentPath("/");
//@@-------------------------------------
      }

      float CFRefInd = cons->CFRefInd();

      CsRCEventAnalysis* ana = CsRCEventAnalysis::Instance();

      Hep3Vector vPoPaMir;
      vPoPaMir = ipf->pPart()->vPoPaMir0()[ipf->pPart()->kDetPart()];
      double rrElmQ = cons->rrAliMir() * cons->rrAliMir();
      string mirEve = "@@@";
      list<CsRCMirrorElem*>::iterator amirr;
      for( amirr=lMirrEleAlg_.begin(); amirr!=lMirrEleAlg_.end(); amirr++ ) {
        //Hep3Vector vCePos = (*amirr)->vpos() + (*amirr)->mirNo()->vC0();
	Hep3Vector vCePos = (*amirr)->vpos();   //   030122
        if ( (vPoPaMir - vCePos).mag2() < rrElmQ ) {
	  mirEve = (*amirr)->name();  break;
	}
      }

//--- mirror elements alignment :
//    ---------------------------
//--- use aligned points within 6 sigma :
//    -----------------------------------
      double theReco = ring->the();
      //double betaReco = 1./( cos( theReco/1000. ) * CFRefInd );
      //if( betaReco > 1.) { betaReco = 1.; }
      //double recoCut = 6.* ana->sigmaPhoRec( betaReco );
      double momPart = ring->pPartPhot()->pPart()->mom();
      double recoCut = 6.* ring->pPartPhot()->sigmaPhoRec( momPart );
      if( theReco > cons->theAliMir() ) {
        double theLo = theReco - recoCut;
        double theUp = theReco + recoCut;
        int kElmAli = 0;
        list<CsRCMirrorElem*>::iterator amirr;
        for( amirr=lMirrEleAlg_.begin(); amirr!=lMirrEleAlg_.end(); amirr++ ) 
        {
          if( (*amirr)->name() == mirEve ) {
            //xh = (*amirr)->vpos().x() + (*amirr)->mirNo()->vC0().x();
            //yh = (*amirr)->vpos().y() + (*amirr)->mirNo()->vC0().y();
	    xh = (*amirr)->vpos().x();   //   030122
            yh = (*amirr)->vpos().y();   //   030122
	    if( hist.hRC3765 ) hist.hRC3765->Fill( xh, yh );
//hh                           ----------------------------
	    list<CsRCPhoton*> lPhotons = ipf->lPhotons();
	    list<CsRCPhoton*>::iterator ih;
	    for( ih=lPhotons.begin(); ih!=lPhotons.end(); ih++ ) {
              double thePhot = (*ih)->the();
              if ( thePhot > theLo && thePhot < theUp ) {
	        xh = thePhot;
		//yh = (*ih)->phi();   //   !!!
	        yh = (*ih)->phiA();
	        if( vRC1800[kElmAli] ) vRC1800[kElmAli]->Fill( xh, yh );
//hh                                   --------------------------------
	      }
	    }
          }
	  kElmAli++;
	}
      }

  }
*/


//============================================================================
  void CsRCMirrors::doAliMirrors() {
//----------------------------------

//--- Paolo --- November  2001
//    ------------------------


      CsRCRecConst *cons = CsRCRecConst::Instance();
      CsRCExeKeys *key = CsRCExeKeys::Instance();

      CsRCHistos& hist = CsRCHistos::Ref();
      double xh, yh, wh;

      static bool firstCall = true;
      if( firstCall ) {
        firstCall = false;
        key->acknoMethod( "CsRCMirrors::MirrorAlignment" );

	CsOpt* opt = CsOpt::Instance();
        list<string> lName;
        bool boo;
        boo = opt->CsOpt::getOpt( "RICHONE", "MirrorElemName", lName );
        if( boo && !lName.empty() ) {
          list<CsRCMirrorElem*>::iterator imirr;
          list<string>::iterator is = lName.begin();
          if( (*is) == "ALL" ) {
            for( imirr=lMirrEle_.begin(); imirr!=lMirrEle_.end(); imirr++ ) {
              (*imirr)->setAlign( true );
	    }
	  }
          else if( (*is) == "UP" ) {
            for( imirr=lMirrEle_.begin(); imirr!=lMirrEle_.end(); imirr++ ) {
              if( (*imirr)->mirpo() == 'U' ) { (*imirr)->setAlign( true ); }
	    }
	  }
          else if( (*is) == "DOWN" ) {
            for( imirr=lMirrEle_.begin(); imirr!=lMirrEle_.end(); imirr++ ) {
              if( (*imirr)->mirpo() == 'D' ) { (*imirr)->setAlign( true ); }
	    }
	  } 
          else {
            for( is=lName.begin(); is!=lName.end(); is++ ) {
              for( imirr=lMirrEle_.begin(); imirr!=lMirrEle_.end(); imirr++ ) {
                if( (*is) == (*imirr)->name() ) { (*imirr)->setAlign( true ); }
	      }
	    }
	  }
          if( lName.front() == "NONE" ) {
            for( imirr=lMirrEle_.begin(); imirr!=lMirrEle_.end(); imirr++ ) {
              (*imirr)->setAlign( false );
	    }
	  }
        }

        int kElm  = 0;
        list<CsRCMirrorNom*>::iterator jmirr;
        list<CsRCMirrorElem*>::iterator imirr;
        for( imirr=lMirrEle_.begin(); imirr!=lMirrEle_.end(); imirr++ ) {
          if( (*imirr)->align() ) { 
            lMirrEleAlg_.push_back( (*imirr) );
            lMirrEleK_.push_back( kElm );
	    //cout << (*imirr)->name() << endl;
	  }

          char mirpo = (*imirr)->mirpo();
          if( mirpo == 'U' ) { jmirr = lMirrNom_.begin(); }
          if( mirpo == 'D' ) { jmirr = lMirrNom_.begin(); jmirr++; }
          //Hep3Vector vCePos = (*imirr)->vpos() + (*jmirr)->vC0();
	  Hep3Vector vCePos = (*imirr)->vpos();   //   030122
	  //cout << (*imirr)->name() << vCePos << endl;
          xh = 6.*float( kElm ) + 0.5;
          wh = vCePos.x();
          if( hist.hRC3750 ) hist.hRC3750->Fill( xh, wh );
//hh                         ----------------------------
          xh = 6.*float( kElm ) + 1.5;
          wh = vCePos.y();
          if( hist.hRC3750 ) hist.hRC3750->Fill( xh, wh );
//hh                         ----------------------------
          xh = 6.*float( kElm ) + 2.5;
          wh = vCePos.z();
          if( hist.hRC3750 ) hist.hRC3750->Fill( xh, wh );
//hh                         ----------------------------
          kElm++;
          xh = vCePos.x();
          yh = vCePos.y();
          if( hist.hRC3751 ) hist.hRC3751->Fill( xh, yh );
//hh                         ----------------------------
        }

        CsHistograms::SetCurrentPath("/RICH");
//@@-----------------------------------------
        static int kOffHi = cons->kOffHi();
        int kmirr = 0;
        list<CsRCMirrorElem*>::iterator amirr;
        for( amirr=lMirrEleAlg_.begin(); amirr!=lMirrEleAlg_.end(); amirr++ )
        {
          string hTitle;
          stringstream hN1800;
          hN1800 << kOffHi + 1800 + kmirr;
          hTitle = "phi photon vs theta photon, align. - in";
          vRC1800.push_back( new CsHist2D( hN1800.str(), hTitle,
       			                   140, 0., 70., 90, 0., 360. ) );
	  hN1800.clear();
	  stringstream nTitle;
          stringstream hN3800;
	  nTitle << (*amirr)->name() << " - fit ";
          hN3800 << kOffHi + 3800 + kmirr;
          //hTitle = "phi photon vs theta photon, align. - fit";
          //vRC3800.push_back( new CsHist2D( hN3800.str(), hTitle,
          vRC3800.push_back( new CsHist2D( hN3800.str(), nTitle.str(),
	     				   140, -70., 70., 140, -70., 70. ) );
	  nTitle.clear();
	  hN3800.clear();
          kmirr++;
        }
        CsHistograms::SetCurrentPath("/");
//@@-------------------------------------
      }


//--- loop over rings :
//    -----------------
      list<CsRCRing*> lRings = CsRCEventRings::Instance()->lRings();
      list<CsRCRing*>::iterator ir;
      for( ir=lRings.begin(); ir!=lRings.end(); ir++ ) {
        if( !(*ir)->flag() ) continue;

        float CFRefInd = cons->CFRefInd();

        CsRCEventAnalysis* ana = CsRCEventAnalysis::Instance();

        CsRCPartPhotons* ipf = (*ir)->pPartPhot();
        Hep3Vector vPoPaMir;
        vPoPaMir = ipf->pPart()->vPoPaMir0()[ipf->pPart()->kDetPart()];

	xh = vPoPaMir.x();
	yh = vPoPaMir.y();
	if( hist.hRC3507 ) hist.hRC3507->Fill( xh, yh );
//hh                       ----------------------------
        double rrElmQ = cons->rrAliMir() * cons->rrAliMir();
        string mirEve = "@@@";
        list<CsRCMirrorElem*>::iterator amirr;
        for( amirr=lMirrEleAlg_.begin(); amirr!=lMirrEleAlg_.end(); 
             amirr++ ) {
          //Hep3Vector vCePos = (*amirr)->vpos() + (*amirr)->mirNo()->vC0();
	  Hep3Vector vCePos = (*amirr)->vpos();   //   030122
          if ( (vPoPaMir - vCePos).mag2() < rrElmQ ) {
	    mirEve = (*amirr)->name();  break;
	  }
        }

//----- mirror elements alignment :
//      ---------------------------
//----- use aligned points within 6 sigma :
//      -----------------------------------
        double theReco = (*ir)->the();
        double betaReco = 1./( cos( theReco/1000. ) * CFRefInd );
        if( betaReco > 1.) betaReco = 1.;

        if( betaReco > 0.9998 ) {                           //   ???
          //double recoCut = 6.* ana->sigmaPhoRec( betaReco );
	  double momPart = (*ir)->pPartPhot()->pPart()->mom();
	  double recoCut = 6.* (*ir)->pPartPhot()->sigmaPhoRec( momPart );
          if( theReco > cons->theAliMir() ) {
            double theLo = theReco - recoCut;
            double theUp = theReco + recoCut;
            int kElmAli = 0;
            list<CsRCMirrorElem*>::iterator amirr;
            for( amirr=lMirrEleAlg_.begin(); amirr!=lMirrEleAlg_.end(); 
                 amirr++ ) {
              if( (*amirr)->name() == mirEve ) {
                //xh = (*amirr)->vpos().x() + (*amirr)->mirNo()->vC0().x();
                //yh = (*amirr)->vpos().y() + (*amirr)->mirNo()->vC0().y();
		xh = (*amirr)->vpos().x();   //   030122
		yh = (*amirr)->vpos().y();   //   030122
	        if( hist.hRC3765 ) hist.hRC3765->Fill( xh, yh );
//hh                               ----------------------------
	        list<CsRCPhoton*> lPhotons = (*ir)->lPhotons();
	        list<CsRCPhoton*>::iterator ih;
	        for( ih=lPhotons.begin(); ih!=lPhotons.end(); ih++ ) {
                  double thePhot = (*ih)->the();
                  if ( thePhot > theLo && thePhot < theUp ) {
	            xh = thePhot;
		    //yh = (*ih)->phi();   //   !!!
   	            yh = (*ih)->phiA();
	            if( vRC1800[kElmAli] ) vRC1800[kElmAli]->Fill( xh, yh );
//hh                                       --------------------------------
		  }
		}
	      }
	      kElmAli++;
	    }
	  }
	}

      }

  }


//============================================================================
  void CsRCMirrors::doAliMirrAll() {
//----------------------------------

//--- Paolo --- November  2001
//    ------------------------


      CsRCRecConst *cons = CsRCRecConst::Instance();
      CsRCExeKeys *key = CsRCExeKeys::Instance();

      CsRCHistos& hist = CsRCHistos::Ref();
      double xh, yh, wh;

      static bool firstCall = true;
      if( firstCall ) {
        firstCall = false;

        CsHistograms::SetCurrentPath("/RICH");
//@@-----------------------------------------
        static int kOffHi = cons->kOffHi();

        string hTitle;
        stringstream hN1980;
        hN1980 << kOffHi + 1980;
        hTitle = "phi photon vs theta photon, align. - in";
        hRC1980 = new CsHist2D( hN1980.str(), hTitle, 140, 0., 70., 
				90, 0., 360. );
        stringstream hN3980;
        hN3980 << kOffHi + 3980;
        hTitle = "phi photon vs theta photon, align. - fit";
        hRC3980 = new CsHist2D( hN3980.str(), hTitle, 140, -70., 70., 
				140, -70., 70. );

        CsHistograms::SetCurrentPath("/");
//@@-------------------------------------
      }


//--- loop over rings :
//    -----------------
      list<CsRCRing*> lRings = CsRCEventRings::Instance()->lRings();
      list<CsRCRing*>::iterator ir;
      for( ir=lRings.begin(); ir!=lRings.end(); ir++ ) {
        if( !(*ir)->flag() ) continue;

        float CFRefInd = cons->CFRefInd();

        CsRCEventAnalysis* ana = CsRCEventAnalysis::Instance();

        CsRCPartPhotons* ipf = (*ir)->pPartPhot();
        Hep3Vector vPoPaMir;

//----- use aligned points within 6 sigma :
//      -----------------------------------
        double theReco = (*ir)->the();
        double betaReco = 1./( cos( theReco/1000. ) * CFRefInd );
        if( betaReco > 1.) betaReco = 1.;

        if( betaReco > 0.9998 ) {                           //   ???
          //double recoCut = 6.* ana->sigmaPhoRec( betaReco );
	  double momPart = (*ir)->pPartPhot()->pPart()->mom();
	  double recoCut = 6.* (*ir)->pPartPhot()->sigmaPhoRec( momPart );
          double theLo = theReco - recoCut;
          double theUp = theReco + recoCut;
	  list<CsRCPhoton*> lPhotons = (*ir)->lPhotons();
	  list<CsRCPhoton*>::iterator ih;
	  for( ih=lPhotons.begin(); ih!=lPhotons.end(); ih++ ) {
	    double thePhot = (*ih)->the();
            if ( thePhot > theLo && thePhot < theUp ) {
	      xh = thePhot;
	      //yh = (*ih)->phi();  //   !!!
	      yh = (*ih)->phiA();
	      if( hRC1980 ) hRC1980->Fill( xh, yh );
//hh                        -----------------------
	    }
	  }
	}

      }

  }


//===========================================================================
  void CsRCMirrors::fitAliMirrors() {
//-----------------------------------


//--- version 0.02,  July  2000
//    -------------------------


      CsRCRecConst *cons = CsRCRecConst::Instance();
      double radDeg = cons->RadDeg();

      CsRCUty* uty = CsRCUty::Instance();

      CsRCHistos& hist = CsRCHistos::Ref();
      double xh, yh, wh;

      CsHistograms::SetCurrentPath("/RICH");
//@@---------------------------------------

      static double thetaRel = 0;

      static bool firstCall = true;
      if( firstCall ) {
        firstCall = false;
	thetaRel = 1000.* acos( 1./ cons->CFRefInd() );
      }

      cout << endl;
      cout << "RICHONE, CsRCMirrors::fitAliMirrors() : mirror alignment fit"
	   << endl;
      cout << "------------------------------------------------------------"
	   << endl; 

      size_t nBinX = 140;
      size_t nBinY =  90;
      double xMin =  0.;
      double xMax = 70.;
      double dX = (xMax - xMin)/float( nBinX );
      int kElm = 0;
      int kElmAli = 0;
      list<CsRCMirrorElem*>::iterator amirr;
      list<int>::iterator kel;
      if( !lMirrEleK_.empty() ) kel = lMirrEleK_.begin();
      for( amirr=lMirrEleAlg_.begin(); amirr!=lMirrEleAlg_.end(); amirr++ ) {
	  kElm = (*kel);  kel++;
          size_t kBin;
	  size_t kBin0;
          int kPoints = 0;
          double theFit[nBinY];
          double phiFit[nBinY];
          double errFit[nBinY];
	  //- fortran indexes for HBOOK...
	  for ( size_t kBiny=1; kBiny <= nBinY; kBiny++) {
	      kBin0 = kBiny * (nBinX + 2);          //  un.flw & ov.flw  !
              double theAv = 0.;
              double phiAv = 0.;
              double binCto = 0.;
	      for ( size_t kBinx=1; kBinx <= nBinX; kBinx++) {
		kBin = kBin0 + kBinx;
                double binC = 0;
                if( vRC1800[kElmAli] ) 
		  binC = vRC1800[kElmAli]->GetBinContent( kBin );
                if ( binC > 0. ) { 
		  //cout << binC << "   " << kBiny << " " << kBinx << endl;
                  double the = xMin + kBinx * dX;
                  theAv += binC * the;
                  binCto += binC;
		}
              }
              if ( binCto > 0. ) { 
                theAv /= binCto;
                theAv += 0.25;;
                phiAv = float( kBiny )*4. + 2.;
		//cout << theAv << "   " << phiAv << endl;
                theFit[kPoints] = theAv * cos( phiAv/radDeg );
                phiFit[kPoints] = theAv * sin( phiAv/radDeg );
		//errFit[kPoints] = 1./(binCto*binCto);
                errFit[kPoints] = binCto;    //   020414
                double xh = theFit[kPoints];
                double yh = phiFit[kPoints];
		if( vRC3800[kElmAli] ) vRC3800[kElmAli]->Fill( xh, yh );
//hh                                   --------------------------------
                kPoints++;
              }
            }
            int nPoint = kPoints;

//--------- fit ring in S.Y. plane :
//          ------------------------
            bool exeFit = true;
	    if ( nPoint < 8 ) exeFit = false;
//@@----------------------------------------
            if( exeFit ) {
              //uty->printVL( "theFit", nPoint, theFit );
	      //uty->printVL( "phiFit", nPoint, phiFit );
	      //uty->printVL( "errFit", nPoint, errFit );
              int nParam = 3;
              double param[nParam];
              //param[0] = 55.;
              param[0] = thetaRel;
              param[1] = 0.;
              param[2] = 0.;
              int iPaFit[nParam];
              iPaFit[0] = 0;
              iPaFit[1] = 0;
	      iPaFit[2] = 0;
	      //uty->printVL( "param", nParam, param );

              CsRCCircleFit oCirclev( nPoint, theFit, phiFit, errFit,
//            -------------------------------------------------------
                                      nParam, param, iPaFit );
              oCirclev.doChiFit();
//            ------------------- 

              float step = 5.;
	      xh = step*float( kElmAli ) + 0.5;
              wh = kElm;
	      if( hist.hRC3760 ) hist.hRC3760->Fill( xh, wh );
//hh                             ----------------------------

              bool circleFlag = oCirclev.flag();
              if( circleFlag ) {

		xh = step*float( kElmAli ) + 1.5;
                wh = oCirclev.para()[0];
		if( hist.hRC3760 ) hist.hRC3760->Fill( xh, wh );
		xh = step*float( kElmAli ) + 2.5;
                wh = oCirclev.para()[1];
		if( hist.hRC3760 ) hist.hRC3760->Fill( xh, wh );
		xh = step*float( kElmAli ) + 3.5;
                wh = oCirclev.para()[2];
		if( hist.hRC3760 ) hist.hRC3760->Fill( xh, wh );
		xh = step*float( kElmAli ) + 4.5;
                wh = kElm;
		if( hist.hRC3760 ) hist.hRC3760->Fill( xh, wh );
//hh                               ----------------------------

		cout << " mirror element  " << (*amirr)->name() << " :  ";
		//cout << " --------------------  " << endl;
                oCirclev.print();
		cout << endl;
                //oCirclev.doHist();
//              -----------------

              } else {

		xh = step*float( kElmAli ) + 1.5;
                wh = 0.;
		if( hist.hRC3760 ) hist.hRC3760->Fill( xh, wh );
		xh = step*float( kElmAli ) + 2.5;
                wh = 0.;
		if( hist.hRC3760 ) hist.hRC3760->Fill( xh, wh );
		xh = step*float( kElmAli ) + 3.5;
                wh = 0.;
		if( hist.hRC3760 ) hist.hRC3760->Fill( xh, wh );
		xh = step*float( kElmAli ) + 4.5;
                wh = kElm;
		if( hist.hRC3760 ) hist.hRC3760->Fill( xh, wh );
//hh                               ----------------------------
              }

            }      /*   end if exeFit   */

            kElmAli++;

      }

      CsHistograms::SetCurrentPath("/");
//@@-----------------------------------

  }


//===========================================================================
  void CsRCMirrors::fitAliMirrAll() {
//-----------------------------------


//--- Paolo --- November 2001
//    -----------------------


      CsRCRecConst *cons = CsRCRecConst::Instance();
      double radDeg = cons->RadDeg();

      CsRCUty* uty = CsRCUty::Instance();

      CsRCHistos& hist = CsRCHistos::Ref();
      double xh, yh, wh;

      CsHistograms::SetCurrentPath("/RICH");
//@@---------------------------------------

      static double thetaRel = 0;

      static bool firstCall = true;
      if( firstCall ) {
        firstCall = false;
	thetaRel = 1000.* acos( 1./ cons->CFRefInd() );
      }

      if( hRC1980 ) {

        size_t nBinX = 140;
        size_t nBinY =  90;
        double xMin =  0.;
        double xMax = 70.;
        double dX = (xMax - xMin)/float( nBinX );

        size_t kBin;
        size_t kBin0;
        int kPoints = 0;
        double theFit[nBinY];
        double phiFit[nBinY];
        double errFit[nBinY];
        //- fortran indexes for HBOOK...
        for ( size_t kBiny=1; kBiny <= nBinY; kBiny++) {
          kBin0 = kBiny * (nBinX + 2);          //  un.flw & ov.flw  !
          double theAv = 0.;
          double phiAv = 0.;
          double binCto = 0.;
	  for ( size_t kBinx=1; kBinx <= nBinX; kBinx++) {
	    kBin = kBin0 + kBinx;
            double binC = 0;
            binC = hRC1980->GetBinContent( kBin );
            if ( binC > 0. ) { 
	      //cout << binC << "   " << kBiny << " " << kBinx << endl;
              double the = xMin + kBinx * dX;
              theAv += binC * the;
              binCto += binC;
	    }
	  }
          if ( binCto > 0. ) { 
            theAv /= binCto;
            theAv += 0.25;;
            phiAv = float( kBiny )*4. + 2.;
	    //cout << theAv << "   " << phiAv << endl;
            theFit[kPoints] = theAv * cos( phiAv/radDeg );
            phiFit[kPoints] = theAv * sin( phiAv/radDeg );
	    //errFit[kPoints] = 1./(binCto*binCto);
            errFit[kPoints] = binCto;
            double xh = theFit[kPoints];
            double yh = phiFit[kPoints];
            if( hRC3980 ) hRC3980->Fill( xh, yh );
//hh                      -----------------------
            kPoints++;
	  }
	}
        int nPoint = kPoints;

//----- fit ring in S.Y. plane :
//      ------------------------
        bool exeFit = true;
        if ( nPoint < 8 ) exeFit = false;
//@@------------------------------------
        if( exeFit ) {
          //uty->printVL( "theFit", nPoint, theFit );
          //uty->printVL( "phiFit", nPoint, phiFit );
          //uty->printVL( "errFit", nPoint, errFit );
          int nParam = 3;
          double param[nParam];
          param[0] = thetaRel;
          param[1] = 0.;
          param[2] = 0.;
          int iPaFit[nParam];
          iPaFit[0] = 0;
          iPaFit[1] = 0;
          iPaFit[2] = 0;
          //uty->printVL( "param", nParam, param );

          CsRCCircleFit oCirclev( nPoint, theFit, phiFit, errFit,
//        -------------------------------------------------------
                                  nParam, param, iPaFit );
          oCirclev.doChiFit();
//        ------------------- 

          bool circleFlag = oCirclev.flag();
          if( circleFlag ) {
	    float step  = 5.;
	    int kElmAli = 118;
	    xh = step*float( kElmAli ) + 1.5;
            wh = oCirclev.para()[0];
      	    if( hist.hRC3760 ) hist.hRC3760->Fill( xh, wh );
	    xh = step*float( kElmAli ) + 2.5;
            wh = oCirclev.para()[1];
	    if( hist.hRC3760 ) hist.hRC3760->Fill( xh, wh );
	    xh = step*float( kElmAli ) + 3.5;
            wh = oCirclev.para()[2];
      	    if( hist.hRC3760 ) hist.hRC3760->Fill( xh, wh );

            cout << " all mirrors  :  ";
	    //cout << " --------------------  " << endl;
            oCirclev.print();
            cout << endl;
	  }
        }      /*   end if exeFit   */

      }

      CsHistograms::SetCurrentPath("/");
//@@-----------------------------------

  }


//==========================================================================
//Hep3Vector CsRCMirrors::vImpMir( Hep3Vector &vra, Hep3Vector &vda,
//-----------------------------------------------------------------------
//                                 Hep3Vector &vr0, float RR ) {
  Hep3Vector CsRCMirrors::vImpMir( const Hep3Vector vra, const Hep3Vector vda,
//----------------------------------------------------------------------------
                                   const Hep3Vector vr0, const float RR ) {

//--- Paolo  -  May 1999.

//--- particle impact on mirror :
//    ---------------------------
      Hep3Vector DD( vra - vr0 );
      double dot = vda * DD;
      double Delta = dot*dot - ( DD*DD - RR*RR );
      double norm = 0.;
      //if ( Delta > 0. )  {  norm = sqrt( Delta );  }
      //norm = - dot + norm;
      if ( Delta < 0. ) return Hep3Vector(0., 0., 0.);
      norm = sqrt( Delta ) - dot;

      return  vra + norm * vda;
  }


//============================================================================
  void CsRCMirrors::doAliMirrPhi() {
//----------------------------------

//--- Paolo --- March  2004
//    ---------------------


      CsRCRecConst *cons = CsRCRecConst::Instance();
      CsRCExeKeys *key = CsRCExeKeys::Instance();

      CsRCHistos& hist = CsRCHistos::Ref();
      double xh, yh, wh;

      double radDeg = cons->RadDeg();

      static bool firstCall = true;
      if( firstCall ) {
        firstCall = false;
        key->acknoMethod( "CsRCMirrors::doAliMirrPhi()" );

	CsOpt* opt = CsOpt::Instance();
        list<string> lName;
        bool boo;
        boo = opt->CsOpt::getOpt( "RICHONE", "MirrorElemName", lName );
        if( boo && !lName.empty() ) {
          list<CsRCMirrorElem*>::iterator imirr;
          list<string>::iterator is = lName.begin();
          if( (*is) == "ALL" ) {
            for( imirr=lMirrEle_.begin(); imirr!=lMirrEle_.end(); imirr++ ) {
              (*imirr)->setAlign( true );
	    }
	  }
          else if( (*is) == "UP" ) {
            for( imirr=lMirrEle_.begin(); imirr!=lMirrEle_.end(); imirr++ ) {
              if( (*imirr)->mirpo() == 'U' ) { (*imirr)->setAlign( true ); }
	    }
	  }
          else if( (*is) == "DOWN" ) {
            for( imirr=lMirrEle_.begin(); imirr!=lMirrEle_.end(); imirr++ ) {
              if( (*imirr)->mirpo() == 'D' ) { (*imirr)->setAlign( true ); }
	    }
	  } 
          else {
            for( is=lName.begin(); is!=lName.end(); is++ ) {
              for( imirr=lMirrEle_.begin(); imirr!=lMirrEle_.end(); imirr++) {
                if( (*is) == (*imirr)->name() ) (*imirr)->setAlign( true );
	      }
	    }
	  }
          if( lName.front() == "NONE" ) {
            for( imirr=lMirrEle_.begin(); imirr!=lMirrEle_.end(); imirr++ ) {
              (*imirr)->setAlign( false );
	    }
	  }
        }

        int kElm  = 0;
        list<CsRCMirrorNom*>::iterator jmirr;
        list<CsRCMirrorElem*>::iterator imirr;
        for( imirr=lMirrEle_.begin(); imirr!=lMirrEle_.end(); imirr++ ) {
          if( (*imirr)->align() ) { 
            lMirrEleAlg_.push_back( (*imirr) );
            lMirrEleK_.push_back( kElm );
	    //cout << (*imirr)->name() << endl;
	  }

          char mirpo = (*imirr)->mirpo();
          if( mirpo == 'U' ) { jmirr = lMirrNom_.begin(); }
          if( mirpo == 'D' ) { jmirr = lMirrNom_.begin(); jmirr++; }
	  Hep3Vector vCePos = (*imirr)->vpos();
	  //cout << (*imirr)->name() << vCePos << endl;
          xh = 6.*float( kElm ) + 0.5;
          wh = vCePos.x();
          if( hist.hRC3750 ) hist.hRC3750->Fill( xh, wh );
//hh                         ----------------------------
          xh = 6.*float( kElm ) + 1.5;
          wh = vCePos.y();
          if( hist.hRC3750 ) hist.hRC3750->Fill( xh, wh );
//hh                         ----------------------------
          xh = 6.*float( kElm ) + 2.5;
          wh = vCePos.z();
          if( hist.hRC3750 ) hist.hRC3750->Fill( xh, wh );
//hh                         ----------------------------
          kElm++;
          xh = vCePos.x();
          yh = vCePos.y();
          if( hist.hRC3751 ) hist.hRC3751->Fill( xh, yh );
//hh                         ----------------------------
        }

        CsHistograms::SetCurrentPath("/RICH");
//@@-----------------------------------------
        static int kOffHi = cons->kOffHi();
        int kmirr = 0;
        list<CsRCMirrorElem*>::iterator amirr;
        for( amirr=lMirrEleAlg_.begin(); amirr!=lMirrEleAlg_.end(); amirr++ )
        {
          string hTitle1;
          stringstream hN7100;
          hN7100 << kOffHi + 7100 + kmirr;
          hTitle1 = "theta photon  vs  phi photon";
          vRC7100.push_back( new CsHist2D( hN7100.str(), hTitle1,
       			                   90, -90., 270., 2, 0., 2. ) );
	  hN7100.clear();

	  stringstream nTitle3;
          stringstream hN7300;
	  nTitle3 << (*amirr)->name() << " - fit ";
          hN7300 << kOffHi + 7300 + kmirr;
          vRC7300.push_back( new CsHist2D( hN7300.str(), nTitle3.str(),
	     				   150, -15., 15., 150, -15., 15. ) );
	  nTitle3.clear();
	  hN7300.clear();

	  stringstream nTitle5;
          stringstream hN7500;
	  nTitle5 << (*amirr)->name() << " - fit ";
          hN7500 << kOffHi + 7500 + kmirr;
          vRC7500.push_back( new CsHist2D( hN7500.str(), nTitle5.str(),
	     				   150, -15., 15., 150, -15., 15. ) );
	  nTitle5.clear();
	  hN7500.clear();

	  stringstream nTitle7;
          stringstream hN7700;
	  nTitle7 << (*amirr)->name() << " yDet  vs  xDet ";
          hN7700 << kOffHi + 7700 + kmirr;
          vRC7700.push_back( new CsHist2D( hN7700.str(), nTitle7.str(),
	     		     150, -1500., 1500., 150, -1500., 1500. ) );
	  nTitle7.clear();
	  hN7700.clear();

          kmirr++;
        }
        stringstream Title40;
        stringstream hN7450;
        Title40 << " nu vs iter";
        hN7450 << kOffHi + 7450;
        hRC7450 = new CsHist2D( hN7450.str(), Title40.str(),
				20, 0., 20., 50, 0., 50. );
        stringstream Title41;
        stringstream hN7451;
        Title41 << " nu vs chisq";
        hN7451 << kOffHi + 7451;
        hRC7451 = new CsHist2D( hN7451.str(), Title41.str(),
				100, 0., 10., 50, 0., 50. );

        stringstream Title60;
        stringstream hN7650;
        Title60 << " nu vs iter";
        hN7650 << kOffHi + 7650;
        hRC7650 = new CsHist2D( hN7650.str(), Title60.str(),
				20, 0., 20., 50, 0., 50. );
        stringstream Title61;
        stringstream hN7651;
        Title61 << " nu vs chisq";
        hN7651 << kOffHi + 7651;
        hRC7651 = new CsHist2D( hN7651.str(), Title61.str(),
				100, 0., 10., 50, 0., 50. );

        stringstream Title80;
        stringstream hN7850;
        Title80 << " yMirr  vs  xMirr";
        hN7850 << kOffHi + 7850;
        hRC7850 = new CsHist2D( hN7850.str(), Title80.str(),
				150, -1500., 1500., 150, -1500., 1500. );

        CsHistograms::SetCurrentPath("/");
//@@-------------------------------------
      }


//--- loop over rings :
//    -----------------
      list<CsRCRing*> lRings = CsRCEventRings::Instance()->lRings();
      list<CsRCRing*>::iterator ir;
      for( ir=lRings.begin(); ir!=lRings.end(); ir++ ) {
        if( !(*ir)->flag() ) continue;

        double theReco = (*ir)->the();
	if( theReco < cons->theAliMir() ) continue;
//      ------------------------------------------
	if( theReco <= 0. ) continue;

	list<CsRCPhoton*> lPhotons = (*ir)->lPhotons();
	int nPhoRing = lPhotons.size();
	if( nPhoRing < int( cons->nPhoAliMir() ) ) continue;
//      ---------------------------------------------------

        CsRCPartPhotons* papho = (*ir)->pPartPhot();
        Hep3Vector vPoPaMir;
        vPoPaMir = papho->pPart()->vPoPaMir0()[papho->pPart()->kDetPart()];
	xh = vPoPaMir.x();
	yh = vPoPaMir.y();
	if( hist.hRC3507 ) hist.hRC3507->Fill( xh, yh );
//hh                       ----------------------------

        double rrElmQ = cons->rrAliMir() * cons->rrAliMir();
//      ---------------------------------------------------
        string mirEve = "@@@";
        list<CsRCMirrorElem*>::iterator amirr;
	CsRCMirrorElem* pMirEve = NULL;
        int kElmAli = 0;
        for( amirr=lMirrEleAlg_.begin(); amirr!=lMirrEleAlg_.end(); 
             amirr++ ) {
	  Hep3Vector vCePos = (*amirr)->vpos();
          if ( (vPoPaMir - vCePos).mag2() < rrElmQ ) {
	    mirEve = (*amirr)->name();
	    pMirEve = (*amirr);
	    break;
	  }
	  kElmAli++;
        }
	if( mirEve == "@@@" ) continue;

	int kDetPart = papho->kDetPart();
	double xPade = papho->vPoPaDetW()[kDetPart].x();
	double yPade = papho->vPoPaDetW()[kDetPart].y();

	xh = vPoPaMir.x();
	yh = vPoPaMir.y();
	if( hRC7850 ) hRC7850->Fill( xh, yh );
//hh                  -----------------------

	double xPaFit[nPhoRing];
	double yPaFit[nPhoRing];
	double erPaFit[nPhoRing];
	double xDeFit[nPhoRing];
	double yDeFit[nPhoRing];
	double erDeFit[nPhoRing];
	double momPart = papho->pPart()->mom();
        double sigPhoRec = papho->sigmaPhoRec( momPart );
//                                ----------------------
	double errFitw = sigPhoRec;
//      --------------------------

//----- use ring points within n sigma :
//      --------------------------------
        float nSigCut = cons->nSigAliMir();
//      ----------------------------------
        double aliCut = nSigCut * sigPhoRec;
        double theLoL = theReco - aliCut;
        double theUpL = theReco + aliCut;
	//list<int>::iterator kel;
	//if( !lMirrEleK_.empty() ) kel = lMirrEleK_.begin();
        //int kElmAli = 0;
        //for( amirr=lMirrEleAlg_.begin(); amirr!=lMirrEleAlg_.end(); 
	//   amirr++ ) {
	//int kElm = (*kel);  kel++;
	//if( (*amirr)->name() == mirEve ) {

	//xh = (*amirr)->vpos().x();
	//yh = (*amirr)->vpos().y();
	if( pMirEve ) {
  	  xh = pMirEve->vpos().x();
	  yh = pMirEve->vpos().y();
	  if( hist.hRC3765 ) hist.hRC3765->Fill( xh, yh );
//hh                         ----------------------------
	}

	list<CsRCPhoton*>::iterator ih;
	int kPoint = 0;
	for( ih=lPhotons.begin(); ih!=lPhotons.end(); ih++ ) {
	  double thePhot = (*ih)->the();
	  if ( thePhot > theLoL  &&  thePhot < theUpL ) {

	    double phiPhot = (*ih)->phiA();
	    //double phiPhot = (*ih)->phi();                  //   !!!
	    if( phiPhot < 0. ) phiPhot += 360.;
	    if( phiPhot > 360. ) phiPhot -= 360.;
	    if( phiPhot > 270. ) phiPhot -= 360.;
	    xh = phiPhot;
	    yh = 0.5;
	    wh = thePhot / theReco;
	    if( vRC7100[kElmAli] ) vRC7100[kElmAli]->Fill( xh, yh, wh );
//hh                               ------------------------------------
	    yh = 1.5;
	    if( vRC7100[kElmAli] ) vRC7100[kElmAli]->Fill( xh, yh );
//hh                               --------------------------------

	    xPaFit[kPoint] = thePhot * cos( phiPhot/radDeg );
	    yPaFit[kPoint] = thePhot * sin( phiPhot/radDeg );
	    erPaFit[kPoint] = errFitw * errFitw;
	    xDeFit[kPoint] = (*ih)->ptToClu()->xc() - xPade;
	    yDeFit[kPoint] = (*ih)->ptToClu()->yc() - yPade;
	    erDeFit[kPoint] = errFitw * errFitw;

	    kPoint++;
	  }
	}
	int nPoint = kPoint;
	if( nPoint < int( cons->nPhoAliMir() ) ) continue;
//      -------------------------------------------------

//----- fit ring in S.Y. plane :
//      ------------------------
	bool exeFitPa = true;
	if( exeFitPa ) {
	  int nParam = 3;
	  double param[nParam];
	  param[0] = theReco;
	  param[1] = 0.;
	  param[2] = 0.;
	  int iPaFit[nParam];
	  iPaFit[0] = 0;
	  iPaFit[1] = 0;
	  iPaFit[2] = 0;

	  CsRCCircleFit oCirclev( nPoint, xPaFit, yPaFit, erPaFit,
//        --------------------------------------------------------
				  nParam, param, iPaFit );

	  oCirclev.doChiFit();
//        -------------------
	  if( oCirclev.flag() ) {

	    xh = oCirclev.para()[1];
	    yh = oCirclev.para()[2];
	    if( vRC7300[kElmAli] ) vRC7300[kElmAli]->Fill( xh, yh );
//hh                               --------------------------------

	    xh = oCirclev.nIter();
	    yh = oCirclev.degFree();
	    if( hRC7450 ) hRC7450->Fill( xh, yh );
//hh                      -----------------------
	    xh = oCirclev.chiSquare();
	    yh = oCirclev.degFree();
	    if( hRC7451 ) hRC7451->Fill( xh, yh );
//hh                     -----------------------
	  }
	}

//----- fit ring in Detector plane :
//      ----------------------------
	bool exeFitDe = true;
	if( exeFitDe ) {
	  int nParam = 3;
	  double param[nParam];
	  param[0] = theReco * 3.36;
	  param[1] = 0.;
	  param[2] = 0.;
	  int iPaFit[nParam];
	  iPaFit[0] = 0;
	  iPaFit[1] = 0;
	  iPaFit[2] = 0;

	  CsRCCircleFit oCirclev( nPoint, xDeFit, yDeFit, erDeFit,
//        --------------------------------------------------------
				  nParam, param, iPaFit );

	  oCirclev.doChiFit();
//        -------------------
          if( oCirclev.flag() ) {

	    xh = oCirclev.para()[1];
	    yh = oCirclev.para()[2];
	    if( vRC7500[kElmAli] ) vRC7500[kElmAli]->Fill( xh, yh );
//hh                               --------------------------------

	    xh = oCirclev.nIter();
	    yh = oCirclev.degFree();
	    if( hRC7650 ) hRC7650->Fill( xh, yh );
//hh                      -----------------------
	    xh = oCirclev.chiSquare();
	    yh = oCirclev.degFree();
	    if( hRC7651 ) hRC7651->Fill( xh, yh );
//hh                      -----------------------

	    for( int kPho=0; kPho<nPoint; kPho++ ) {
	      xh = xDeFit[kPho] + xPade;
	      yh = yDeFit[kPho] + yPade;
	      if( vRC7700[kElmAli] ) vRC7700[kElmAli]->Fill( xh, yh );
//hh                                 --------------------------------
	    }
	  }
	}

      }
      //kElmAli++;
      //}
      //}

  }

