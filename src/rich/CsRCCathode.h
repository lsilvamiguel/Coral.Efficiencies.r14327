
#ifndef  RCCATHODE_H
#define  RCCATHODE_H

/*!
   \file    CsRCcathode.h
   \----------------------
   \brief   CsRCCathode class declaration.
   \author  Paolo Schiavon
   \version 4.0
   \date    21 February 2001, rev. August 2005
*/


  #include <string>
  #include <CLHEP/Vector/ThreeVector.h>

  #include "CLHEP/Matrix/Matrix.h"
  #include "CLHEP/Matrix/SymMatrix.h"
  #include "CLHEP/Matrix/DiagMatrix.h"
  #include "CLHEP/Matrix/Vector.h"

//------------------------------------


  class CsRCCathode {


  public:

    CsRCCathode( int, int, std::string, std::string, double, double, double,
                 int, int, double, double, double, double, CLHEP::HepMatrix );

    CsRCCathode( const CsRCCathode& );

    inline  int kCat() const { return kCat_; };
    inline  int ID() const { return ID_; };
    inline  std::string TBname() const { return TBname_; };
    inline  std::string name() const { return name_; };
    inline  int number() const { return number_; };
    inline  CLHEP::Hep3Vector vOffCat0() { return vOffCat0_; };
    inline  CLHEP::Hep3Vector vCat0() const { return vCat0_; };
    inline  int nPadx() const { return nPadx_; };
    inline  int nPady() const { return nPady_; };
    inline  double padx() const { return padx_; };
    inline  double pady() const { return pady_; };
    inline  int nPPadx() const { return nPPadx_; };
    inline  int nPPady() const { return nPPady_; };
    inline  double ppadx() const { return ppadx_; };      // ???
    inline  double ppady() const { return ppady_; };      // ???
    inline  int nPMTx() const { return nPMTx_; };
    inline  int nPMTy() const { return nPMTy_; };
    inline  double PMTx() const { return PMTx_; };
    inline  double PMTy() const { return PMTy_; };
    inline  double hCatx() const { return hCatx_; };
    inline  double hCaty() const { return hCaty_; };
    inline  double xLimMnW() const { return xLimMnW_; };
    inline  double xLimMxW() const { return xLimMxW_; };
    inline  double yLimMnW() const { return yLimMnW_; };
    inline  double yLimMxW() const { return yLimMxW_; };
    inline  void setXLimMnW( double lim ) { xLimMnW_ = lim; };
    inline  void setXLimMxW( double lim ) { xLimMxW_ = lim; };
    inline  void setYLimMnW( double lim ) { yLimMnW_ = lim; };
    inline  void setYLimMxW( double lim ) { yLimMxW_ = lim; };
    inline  double ddQzW() const { return ddQzW_; };
    inline  double ddGap() const { return ddGap_; };
    inline  CLHEP::HepMatrix rotMatrix() const { return rotMatrix_; };

    inline  CLHEP::Hep3Vector vOffCatW() const { return vOffCatW_; };
    inline  void setOffCatW( CLHEP::Hep3Vector vOffCatW ) { vOffCatW_ = vOffCatW; };

    inline  bool isPMT() const { return isPMT_; };
    inline  bool isAPV() const { return isAPV_; };

    inline  double* xcePMTWC() { return xcePMTWC_; };
    inline  double* ycePMTWC() { return ycePMTWC_; };
    inline  double* xCoPPadPMT() { return xCoPPadPMT_; };
    inline  double* yCoPPadPMT() { return yCoPPadPMT_; };
    inline  double* xcePPadPMT() { return xcePPadPMT_; };
    inline  double* ycePPadPMT() { return ycePPadPMT_; };

    bool ccePMTWC( const int, double&, double& );
    bool ccePPadWC( const int, const int, double&, double& );

    void print() const;
    void plotGeo() const;

    ~CsRCCathode();

  
  private:


    int kCat_;
    int ID_;
    std::string  TBname_;
    std::string  name_;
    int number_;
    CLHEP::Hep3Vector vOffCat0_;
    CLHEP::Hep3Vector vCat0_;
    CLHEP::Hep3Vector vOffCatW_;
    int nPadx_, nPady_;
    double padx_, pady_;
    int nPPadx_, nPPady_;
    double ppadx_, ppady_;
    int nPMTx_, nPMTy_;
    double PMTx_, PMTy_;
    double hCatx_, hCaty_;
    double xLimMnW_, xLimMxW_;
    double yLimMnW_, yLimMxW_;
    double ddQzW_;
    double ddGap_;
    CLHEP::HepMatrix rotMatrix_;

    bool isPMT_, isAPV_;

    double xcePMTWC_[12];
    double ycePMTWC_[12];
    double xCoPPadPMT_[25];
    double yCoPPadPMT_[25];
    double xcePPadPMT_[16];
    double ycePPadPMT_[16];

  };

#endif
