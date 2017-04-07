#ifndef  RCPHOTONDET_H
#define  RCPHOTONDET_H

/*!
   \file    CsRCPhotonDet.h
   \------------------------
   \brief   CsRCPhotonDet class declaration.
   \author  Paolo Schiavon
   \version 1.0
   \date    28 August 2000, rev. 21 February 2001
*/


  #include <string>
  #include <CLHEP/Vector/ThreeVector.h>

  #include "CLHEP/Matrix/Matrix.h"
  #include "CLHEP/Matrix/SymMatrix.h"
  #include "CLHEP/Matrix/DiagMatrix.h"
  #include "CLHEP/Matrix/Vector.h"

//------------------------------------


    class CsRCPhotonDet  {


        public:

    CsRCPhotonDet( int, std::string, double, double, double, CLHEP::HepMatrix );

    inline  int kDet() const { return kDet_; };
    inline  std::string name() const { return name_; };
    inline  CLHEP::Hep3Vector vDet0() const { return vDet0_; };
    inline  CLHEP::Hep3Vector vDet0In() const { return vDet0In_; };

    inline  CLHEP::HepMatrix rotMatrix() const { return rotMatrix_; };
    inline  double detAng() const { return detAng_; };
    inline  CLHEP::Hep3Vector vDcDet() const { return vDcDet_; };
    inline  char detpo() const { return detpo_; };

    inline  CLHEP::Hep3Vector vDetW() const { return vDetW_; };
    inline  void setDetW( CLHEP::Hep3Vector vDetW ) { vDetW_ = vDetW; };
    inline  double dyDRSMWR() const { return dyDRSMWR_; };
    inline  void setDRSMWR( double dyDRSMWR ) { dyDRSMWR_ = dyDRSMWR; };

    void print() const;

    ~CsRCPhotonDet();

  
        private:

    CsRCPhotonDet( const CsRCPhotonDet& );
    CsRCPhotonDet& operator=( const CsRCPhotonDet& );

    int kDet_;
    std::string  name_;
    CLHEP::Hep3Vector vDet0_;
    CLHEP::Hep3Vector vDet0In_;

    CLHEP::HepMatrix rotMatrix_;
    double detAng_;
    CLHEP::Hep3Vector vDcDet_;
    char detpo_;

    CLHEP::Hep3Vector vDetW_;
    double dyDRSMWR_;

  };

  #endif
