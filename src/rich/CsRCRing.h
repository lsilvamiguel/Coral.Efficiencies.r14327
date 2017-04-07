
#ifndef  RCRING_H
#define  RCRING_H

/*!
   \file    CsRCRing.h
   \------------------------
   \brief   CsRCRing class declaration.
   \author  Paolo Schiavon
   \version 0.01
   \date    October 2000
*/


#include <list>
#include <vector>


  #include "CsHist.h"
  #include "CsRCRecConst.h"

//----------------------------
  class CsRCPartPhotons;
  class CsRCPhoton;
  class CsRCCircleFit;
  class CsRCEllipseFit;
//----------------------------


  class CsRCRing {


  public:


    CsRCRing();

    void getRingPk( const int, CsRCPartPhotons* );

    CsRCRing( const CsRCRing& );

    //^void getRingBf();
    //^void getRingWa();
    //^void getRing3S();
    void setTime();

    void getRingMaxLike( const int, CsRCPartPhotons* );

    void histCntRatio();
    void peakScan();
    bool rejectWrongPeak( const float );
    bool peakCountSearch( const CsRCPartPhotons*, int& );
    //^bool peakCountSearchTEST( const CsRCPartPhotons*, int& );
    //^bool peakCountSearchVB( const CsRCPartPhotons*, int& );
    //^bool peakCountSearchMC( const CsRCPartPhotons*, int& );
    //^bool peakChiSearch( const CsRCPartPhotons*, int& );
    //^bool peakMassSearch( const CsRCPartPhotons*, int& );
    //^bool peakLikeSearch( const CsRCPartPhotons*, int& );
    void histPeakSearch( const CsRCPartPhotons* );
    int getVarMcan( const CsRCPartPhotons*, const int );
    int getVarNMore( const CsRCPartPhotons*, const int );

    void setThetaType();
    double getThetaWgav();
    double getThetaRFit();

    long double getLikelihood( const double, const double, int& );
    long double getLikeRing( const double, const double, int& );

    long double ringPoissProb( const double );
    double ringGaussProb( const double );
    double getNPhotExpected( const double );
    double getNPhotCorrected( const double );
    bool getRingness( const std::list<CsRCPhoton*>&, double& );

    void checkChiSquare( const double );
    bool getThetaLikeMax( long double&, double& );
    void checkLikelihood( const double );
    void checkLikeVsIndex( const double );
    double getQsquare( const double, int& );
    double getRingQsq( const double, int& );
    void checkChiReso();

    float recoCut( const double );
    float sigmaRing( const double );
    float sigmaRing3( const double );
    float sigmaRing4( const double );

    bool getDetRFit( const std::vector<CsHist2D*>&, CsRCCircleFit& );
    bool getDetEFit( const std::vector<CsHist2D*>&, CsRCEllipseFit& );

    bool ringPMTonly();
    bool ringAPVonly();

    inline  int kRing() const { return kRing_; };
    inline  CsRCPartPhotons* pPartPhot() const { return pPartPhot_; };

    inline  double the() const { return the_; };
    inline  const std::list<CsRCPhoton*> &lPhotons() const { return lPhotons_; };
    inline  int kRingLoc() const { return kRingLoc_; } ;

    inline  bool flagReco() const { return flagReco_; };
    inline  bool flagOverThrs() const { return flagOverThrs_; };
    inline  bool flag() const { return flag_; };

    inline  double theLoL() const { return theLoL_; };
    inline  double theUpL() const { return theUpL_; };

    inline  double theReco() const { return theReco_; };
    inline  double thetaWgav() const { return thetaWgav_; };
    inline  double thetaRFit() const { return thetaRFit_; };
    inline  double thetaLike() const { return thetaLike_; };
    inline  bool thetaLikeSet() const { return thetaLikeSet_; };
    inline void setThetaLike( const double thetaLike ) {
      thetaLike_ = thetaLike;
      thetaLikeSet_ = true;
    }
    inline  double mTime() const { return mTime_; };
    inline  double mT0() const { return mT0_; };
    inline  bool mTimeSet() const { return mTimeSet_; };
    inline  int nPhotPMT() { return nPhotPMT_; };

    inline  float* binCont() { return binCont_; };

    inline  double* partProbs() { return partProbs_; };
    inline  double partProb( const int k ) { return partProbs_[k]; };
    inline  bool partProbsSet() const { return partProbsSet_; };
    inline  void setPartProbs( const double* paPro ) {
      int nProb = CsRCRecConst::Instance()->outBufferSize();
      for( int j=0; j<nProb; j++ ) partProbs_[j] = paPro[j];
      partProbsSet_ = true;
    };
    inline  void setPartProb( const int k, const double prob ) {
      partProbs_[k] = prob; }
    inline  void setPartProbsSet() { partProbsSet_ = true; };

    inline  double* probaLKRing() { return probaLKRing_; };
    inline  double* derivLKRing() { return derivLKRing_; };
    inline  double probaLKBgRing() { return probaLKBgRing_; };
    inline  void setProbaLKRing( const double* proba ) {
      for( int j=0; j<31; j++ ) probaLKRing_[j] = proba[j]; };
    inline  void setDerivLKRing( const double* deriv ) {
      for( int j=0; j<31; j++ ) derivLKRing_[j] = deriv[j]; };
    inline  void setProbaLKBgRing( const double proba ) { 
      probaLKBgRing_ = proba; };

    inline  int nPhotQsQ() { return nPhotQsQ_; };
    inline  double ringQsQ() { return ringQsQ_; };
    inline  double* qSquareRing() { return qSquareRing_; };

    inline  void setnPhotQsQ( const int npho ) { nPhotQsQ_ = npho; };
    inline  void setRingQsQ( const double qsq ) { ringQsQ_ = qsq; };
    inline  void setQSquareRing( double* square ) {
      for( int j=0; j<31; j++ ) qSquareRing_[j] = square[j]; };

    inline  void setFlagReco( const bool flag ) { flagReco_ = flag; };
    inline  void setFlagOverThrs( const bool flag ) { flagOverThrs_ = flag; };
    inline  void setFlag( const bool flag ) { flag_ = flag; };

    inline  void setFlagBack( const bool flag ) { flagBack_ = flag; };
    inline  bool flagBack() const { return flagBack_; };
    inline  double ringBack() const { 
      if( flagBack_ ) return ringBack_;
      else  return  0.;
    };
    inline void setRingBack( const double back ) {
      ringBack_ = back;
      flagBack_ = true;
    };

    void print() const;
    void printPhotons() const;
    void printPhoFlag() const;

    ~CsRCRing();


  private:


    int kRing_;
    CsRCPartPhotons* pPartPhot_;

    double the_;
    std::list<CsRCPhoton*> lPhotons_;
    int kRingLoc_;

    std::list<CsRCPhoton*> lDoublePhotons_;

    double theLoL_;
    double theUpL_;

    double theReco_;
    double thetaWgav_;
    double thetaRFit_;
    double thetaLike_;
    bool thetaLikeSet_;

    double mTime_;
    double mT0_;
    bool mTimeSet_;
    int nPhotPMT_;

    float binCont_[200];

    double partProbs_[50];
    bool partProbsSet_;

    double probaLKRing_[31];
    double derivLKRing_[31];
    double probaLKBgRing_;

    int nPhotQsQ_;
    double ringQsQ_;
    double qSquareRing_[31];

    bool flagReco_;
    bool flagOverThrs_;
    bool flag_;

    bool flagBack_;
    double ringBack_;

  };

#endif
