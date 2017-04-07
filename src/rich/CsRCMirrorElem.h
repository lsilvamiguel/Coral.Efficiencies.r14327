#ifndef  MIRRORELEM_H
#define  MIRRORELEM_H

/*!
   \file    CsRCMirrorElem.h
   \------------------------
   \brief   CsRCMirrorElem class declaration.
   \author  Paolo Schiavon
   \version 1.0
   \date    10 August 2000
*/

  #include "CsOpt.h"

  #include <CLHEP/Vector/ThreeVector.h>

  class CsRCMirrorNom;


  class CsRCMirrorElem  {


  public:


    CsRCMirrorElem();

    CsRCMirrorElem( const int, const std::string, 
		    const double, const double,
		    const double, const double, const double,
                    const double, const double, const int );

    CsRCMirrorElem( const CsRCMirrorElem& );
    CsRCMirrorElem& operator=( const CsRCMirrorElem& );

    inline  int kMir() const { return kMir_; };
    inline  std::string name() const { return name_; };
    inline  double theta() const { return theta_; };
    inline  double phi() const { return phi_; };
    inline  double RR() const { return RR_; };
    inline  double deTheta() const { return deTheta_; };
    inline  double dePhi() const { return dePhi_; };
    inline  double delta() const { return delta_; };
    inline  double qfact() const { return qfact_; };
    inline  char mirpo() const { return mirpo_; };
    inline  CsRCMirrorNom* mirNo() const { return mirNo_; }
    inline  bool align() const { return align_; };
    inline  void setAlign( bool a )  { align_ = a; };

    inline  double ddEleQ() const { return ddEleQ_; };
    inline  void setddEleQ( double dd )  { ddEleQ_ = dd; };

    inline CLHEP::Hep3Vector vpos() const { return vpos_; };
    inline CLHEP::Hep3Vector vC0() const { return vC0_; };

    inline  CLHEP::Hep3Vector vtan() const { return vtan_; };

    void print() const;

    ~CsRCMirrorElem();

  
  private:


    int kMir_;
    std::string name_;
    double theta_;
    double phi_;
    double RR_;
    double deTheta_;
    double dePhi_;
    double delta_;
    double qfact_;
    char mirpo_;
    CsRCMirrorNom* mirNo_;
    bool align_;

    CLHEP::Hep3Vector vpos_;
    CLHEP::Hep3Vector vC0_;

    CLHEP::Hep3Vector vtan_;
    //??    CLHEP::Hep3Vector vdeTan_;

    double ddEleQ_;

  };

  #endif
