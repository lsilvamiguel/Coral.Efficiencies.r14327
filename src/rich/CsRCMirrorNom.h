#ifndef  RCMIRRORNOM_H
#define  RCMIRRORNOM_H

/*!
   \file    CsRCMirrorNom.h
   \------------------------
   \brief   CsRCMirrorNom class declaration.
   \author  Paolo Schiavon
   \version 1.0
   \date    25 August 2000
*/


  #include <string>

  #include <CLHEP/Vector/ThreeVector.h>

  class CsRCMirrorNom {


  public:


    CsRCMirrorNom();

    CsRCMirrorNom( const int, const std::string,
		   const double, const double, const double, const double );

    CsRCMirrorNom( const CsRCMirrorNom& );
    CsRCMirrorNom& operator=( const CsRCMirrorNom& );

    inline  int kMir() const { return kMir_; };
    inline  std::string name() const { return name_; };
    inline  CLHEP::Hep3Vector vC0nom() const { return vC0nom_; };
    inline  CLHEP::Hep3Vector vC0() const { return vC0_; };
    inline  double RRnom() const { return RRnom_; };
    inline  double RR() const { return RR_; };
    inline  char mirpo() const { return mirpo_; };

    void print() const;

    ~CsRCMirrorNom();

  
  private:


    int kMir_;
    std::string name_;
    CLHEP::Hep3Vector vC0nom_;
    CLHEP::Hep3Vector vC0_;
    double RRnom_;
    double RR_;
    char mirpo_;

  };

#endif
