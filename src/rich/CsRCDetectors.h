
#ifndef  RCDETECTORS_H
#define  RCDETECTORS_H

/*!
   \file    CsRCDetectors.h
   \-------------------------
   \brief   CsRCDetectors class declaration.
   \author  Paolo Schiavon
   \version 0.02,  rev. 20/6/00
   \date    7 December 1999, rev. August 2005
*/

  #include <list>
  #include "CsErrLog.h"

  #include <CLHEP/Vector/ThreeVector.h>
//-------------------------------------

  #include "CLHEP/Matrix/Matrix.h"
  #include "CLHEP/Matrix/SymMatrix.h"
  #include "CLHEP/Matrix/DiagMatrix.h"
  #include "CLHEP/Matrix/Vector.h"
//--------------------------------------

//---------------------------
  #include "CsRCPhotonDet.h"
  #include "CsRCCathode.h"
  #include "CsRichOne.h"
//---------------------------


  class CsRCDetectors {


  public:

    static CsRCDetectors* Instance();

    void setPhotonDet( std::string, double, double, double, CLHEP::HepMatrix );

    void setCathode( int, std::string, std::string, double, double, double,
		     int, int, double, double, double, double, 
		     CLHEP::HepMatrix );

    CsRCDetectors( const CsRCDetectors& );

    void setLocalGeo();

    inline  int nDetector() const { return lPhoDet_.size(); };
    inline  int nCathode() const { return lCathodes_.size(); };

    inline  const std::list<CsRCPhotonDet*> &lPhoDet() const {
      return lPhoDet_; };
    inline  const std::list<CsRCCathode*> &lCathodes() const {
      return lCathodes_; };

    inline  const std::vector<int> &nCatPMT() const { return nCatPMT_; };
    inline  const std::vector<int> &nCatAPV() const { return nCatAPV_; };

    inline  CLHEP::Hep3Vector vDet0v( int k ) const {
      std::list<CsRCPhotonDet*>::const_iterator id = lPhoDet_.begin();
      for( int j=0; j<k; j++ ) id++;
      return (*id)->vDet0();
    };
    inline  CLHEP::Hep3Vector vDcDetv( int k ) const {
      std::list<CsRCPhotonDet*>::const_iterator id = lPhoDet_.begin();
      for( int j=0; j<k; j++ ) id++;
      return (*id)->vDcDet();
    };

    CsRCCathode* ptrToCat( int cat );
    /*
    inline CsRCCathode* ptrToCat( int cat ) const {
      std::list<CsRCCathode*>::const_iterator ic;
      CsRCCathode* ptr = NULL;
      for( ic=lCathodes_.begin(); ic!=lCathodes_.end(); ic++ ) {
	if( (*ic)->number() == cat ) {
	  ptr = (*ic);
	  break;
	}
      }
      if( ptr == NULL ) {
	std::cout << " RICHONE, CsRCDetectors : Cathode number " << cat;
	std::cout << "   richEvent  ";
	std::cout << CsRichOne::Instance()->kEvent() << std::endl;
        std::string str =
	  "RICHONE, CsRCDetectors::ptrToCat() : wrong cathode number";
        CsErrLog::Instance()->mes( elFatal, str );
      }
      return ptr;
    };
    */

    inline int catNfromName( const std::string name ) const {
      std::list<CsRCCathode*>::const_iterator ic;
      int catNum = -1;
      for( ic=lCathodes_.begin(); ic!=lCathodes_.end(); ic++ ) {
	if( (*ic)->name() == name ) {
	  catNum = (*ic)->number();
	  break;
	}
      }
      if( catNum < 0 ) {
        std::string str =
	  "RICHONE, CsRCDetectors::catNfromName() : wrong cathode name";
        CsErrLog::Instance()->mes( elFatal, str );
      }
      return  catNum;
    };

    inline int catNfromTBname( const std::string TBname ) const {
      std::list<CsRCCathode*>::const_iterator ic;
      int catNum = -1;
      for( ic=lCathodes_.begin(); ic!=lCathodes_.end(); ic++ ) {
	if( (*ic)->TBname() == TBname ) {
	  catNum = (*ic)->number();
	  break;
	}
      }
      if( catNum < 0 ) {
        std::string str = 
	  "RICHONE, CsRCDetectors::catNfromTBName() : wrong cathode TBname";
        CsErrLog::Instance()->mes( elFatal, str );
      }
      return  catNum;
    };

    inline CLHEP::Hep3Vector vOffCatW( int cat ) const {
      std::list<CsRCCathode*>::const_iterator ic;
      CLHEP::Hep3Vector vOff( 1000000., 1000000., 1000000. );
      for( ic=lCathodes_.begin(); ic!=lCathodes_.end(); ic++ ) {
	if( (*ic)->number() == cat ) {
//------- detector Up = 0, Down = 1;
	  vOff.setX( (*ic)->vOffCatW().x() );
	  vOff.setY( (*ic)->vOffCatW().y() );
	  vOff.setZ( (*ic)->vOffCatW().z() );
	  break;
	}
      }
      if( vOff.x() == 1000000. ) {
        std::string str =
	  "RICHONE, CsRCDetectors::vOffCatW() : wrong cathode number";
        CsErrLog::Instance()->mes( elFatal, str );
      }
      return vOff;
    };

    inline int cathodePos( int cat ) const {
      std::list<CsRCCathode*>::const_iterator ic;
      int det = -1;
      for( ic=lCathodes_.begin(); ic!=lCathodes_.end(); ic++ ) {
	if( (*ic)->number() == cat ) {
//------- detector Up = 0, Down = 1;
	  if( (*ic)->vOffCatW().y() >= 0. ) det = 0;
	  if( (*ic)->vOffCatW().y() <  0. ) det = 1;
	  break;
	}
      }
      if( det == -1 ) {
        std::string str =
	  "RICHONE, CsRCDetectors::cathodePos() : wrong cathode number";
        CsErrLog::Instance()->mes( elFatal, str );
      }
      return det;
    };

    inline float zEntrWind() const { return zEntrWind_; };
    inline float zExitWind() const { return zExitWind_; };

    void printDetAll() const;
    void printCatAll() const;

    CLHEP::Hep3Vector vImpDetW( const CLHEP::Hep3Vector, const CLHEP::Hep3Vector ) const;


  protected:

    CsRCDetectors();

    ~CsRCDetectors();

  
  private:

    static CsRCDetectors* instance_;

    std::list<CsRCPhotonDet*> lPhoDet_;
    std::list<CsRCCathode*> lCathodes_;

    float zEntrWind_;
    float zExitWind_;

    std::vector<int> nCatPMT_;
    std::vector<int> nCatAPV_;

  };

#endif

