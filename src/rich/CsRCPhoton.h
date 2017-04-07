
#ifndef  RCPHOTON_H
#define  RCPHOTON_H

/*!
   \file    CsRCPhoton.h
   \-----------------------
   \brief   CsRCPhoton class declaration.
   \author  Paolo Schiavon
   \version 1.0
   \date    1 October 2000, rev. August 2005
*/


  class CsRCPhotonDet;
  class CsRCPartPhotons;

  #include "CsRCCluster.h"

  class CsRCPhoton {


  public:


    CsRCPhoton();

    CsRCPhoton( int, double, double, double, double, double, double,
		CsRCCluster*, CsRCPhotonDet*, CsRCPartPhotons*, int );

    CsRCPhoton( const CsRCPhoton& );

    inline  int kPhot() const { return kPhot_; };
    inline  double the0() const { return the0_; };
    inline  double the() const { return the_; };
    inline  double theNorm() const { return theNorm_; };
    inline  double phi() const { return phi_; };
    inline  double phiA() const { return phiA_; };
    inline  double theB() const { return theB_; };
    inline  double theM() const { return theM_; };
    inline  double PH() const { return PH_; };
    inline  CsRCCluster* ptToClu() const { return ptToClu_; };
    inline  void setPtToClu( CsRCCluster* ptToClu ) { ptToClu_ = ptToClu; };

    inline  void setThe( const double the ) { the_ = the; }

    inline  CsRCPhotonDet* ptToDet() const { return ptToDet_; };
    inline  void setPtToDet( CsRCPhotonDet* ptToDet ) { ptToDet_ = ptToDet; };

    inline  CsRCPartPhotons* pPartPhot() const { return pPartPhot_; };

    double sigmaPhoPid( const CsRCPartPhotons* ) const;
    double sigmaPhot6() const;
    double sigmaPhot6APV() const;
    double sigmaPhot6PMT() const;

    inline  int kDetClu() const { return kDetClu_; };

    inline  bool flag() const { return flag_; };
    inline  void setFlag( bool flag ) { flag_ = flag; };

    inline  bool isPMT() const { return ptToClu_->isPMT(); };
    inline  bool isAPV() const { return ptToClu_->isAPV(); };

    inline  double CFRefInd() const { return CFRefInd_; };
    inline  double thetaIpo( int kPa ) const {
      if( kPa <= 0 ) return 0.;
      return thetaIpo_[(kPa-2)/3+1]; };

    double facVStoUV() const;
    double thetaVStoUV( const double ) const;
    double thetaUVtoVS( const double ) const;
    double getThetaIpo( const double ) const;
    inline  void setTheToNorm() { the_ = theNorm_; };
    inline  void setTheTo0() { the_ = the0_; };

    inline  bool likeFirst() const { return likeFirst_; };
    inline  void setLikeFirst( const bool first ) { likeFirst_ = first; };

    void  flagAllPhotons( const CsRCPhoton* ) const;

    double getPMTOptCorr();

    void print() const;

    ~CsRCPhoton();


  private:


    int kPhot_;
    double the0_;
    double the_;
    double theNorm_;
    double phi_;
    double phiA_;
    double theB_;
    double theM_;
    double PH_;
    CsRCCluster* ptToClu_;
    CsRCPhotonDet* ptToDet_;
    CsRCPartPhotons* pPartPhot_;
    int kDetClu_;
    bool flag_;

    double CFRefInd_;
    double thetaIpo_[5];
    double facVStoUV_;

    bool likeFirst_;

  };

  #endif
