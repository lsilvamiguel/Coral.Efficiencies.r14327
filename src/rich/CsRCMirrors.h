#ifndef  RCMIRRORS_H
#define  RCMIRRORS_H

/*!
   \file    CsRCMirrors.h
   \-------------------------
   \brief   CsRCMirrors class declaration.
   \author  Paolo Schiavon
   \version 0.02,  rev. 20/6/00
   \date    7 December 1999
*/


  #include "CsHistograms.h"

//---------------------------
  #include "CsRCMirrorNom.h"
  #include "CsRCMirrorElem.h"

  class CsRCParticle;
  class CsRCPartPhotons;
  class CsRCPartRing;
  class CsRCRing;
//---------------------------

  #include <CLHEP/Vector/ThreeVector.h>

//-------------------------------------


    class CsRCMirrors {


    public:


      static CsRCMirrors* Instance();

      void setMirrNom( const std::string,
		       const double, const double, const double,
		       const double );

      void setMirrEle( const std::string, const double, const double,
		       const double, const double, const double,
		       const double, const double, int);

      CsRCMirrors( const CsRCMirrors& );

      bool doSelMirrors( CsRCParticle* );

      void doAliMirrors( CsRCPartPhotons*, CsRCRing* );
      void doAliMirrors();
      void doAliMirrAll();
      void doAliMirrPhi();

      void fitAliMirrors();
      void fitAliMirrAll();
      void fitAliMirrPhi();

      inline  int nMirror() const { return lMirrNom_.size(); };

      inline  const std::list<CsRCMirrorNom*> &lMirrNom() const {
	return lMirrNom_; };
      inline  const std::list<CsRCMirrorElem*> &lMirrEle() const {
	return lMirrEle_; };
      inline  const std::list<CsRCMirrorElem*> &lMirrEleAlg() const {
	return lMirrEleAlg_; };

      inline  int nMiEle() const { return nMirEleSel_; };

      inline  CLHEP::Hep3Vector vC0v( int k ) const {
	std::list<CsRCMirrorNom*>::const_iterator im = lMirrNom_.begin();
	for( int j=0; j<k; j++ ) im++;
	return (*im)->vC0();
      };
      inline  double RRv( int k ) const {
	std::list<CsRCMirrorNom*>::const_iterator im = lMirrNom_.begin();
	for( int j=0; j<k; j++ ) im++;
	return (*im)->RR();
      };

      inline  CsRCMirrorElem* pMirrElePa() { return pMirrElePa_; };
      inline  CLHEP::Hep3Vector vCorPoPa0() const { return vCorPoPa0_; };
      inline  float corPaRR() const { return corPaRR_; };
      inline  CLHEP::Hep3Vector vCorPoPho0() const { return vCorPoPho0_; };
      inline  float corPhoRR() const { return corPhoRR_; };

      void printNomAll() const;
      void printEleAll() const;

      //CLHEP::Hep3Vector vImpMir( CLHEP::Hep3Vector&, CLHEP::Hep3Vector&, CLHEP::Hep3Vector&, float );
      CLHEP::Hep3Vector vImpMir( const CLHEP::Hep3Vector, const CLHEP::Hep3Vector, 
			  const CLHEP::Hep3Vector, const float );


    protected:


      CsRCMirrors();

      ~CsRCMirrors();

  
    private:


      static CsRCMirrors* instance_;

      std::list<CsRCMirrorNom*> lMirrNom_;
      std::list<CsRCMirrorElem*> lMirrEle_;
      std::list<CsRCMirrorElem*> lMirrEleAlg_;
      std::list<int> lMirrEleK_;

      int nMirEleSel_;

      CsRCMirrorElem* pMirrElePa_;
      CLHEP::Hep3Vector vCorPoPa0_;
      float corPaRR_;
      CLHEP::Hep3Vector vCorPoPho0_;
      float corPhoRR_;

      CsHist2D* hRC1980;
      CsHist2D* hRC3980;
      std::vector<CsHist2D*> vRC1800;
      std::vector<CsHist2D*> vRC3800;

      std::vector<CsHist2D*> vRC7100;
      std::vector<CsHist2D*> vRC7300;
      CsHist2D* hRC7450;
      CsHist2D* hRC7451;
      std::vector<CsHist2D*> vRC7500;
      CsHist2D* hRC7650;
      CsHist2D* hRC7651;
      std::vector<CsHist2D*> vRC7700;
      CsHist2D* hRC7850;

  };

#endif

