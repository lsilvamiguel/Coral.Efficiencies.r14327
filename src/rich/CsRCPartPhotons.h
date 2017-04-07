
#ifndef  PARTPHOTONS_H
#define  PARTPHOTONS_H

/*!
   \file    CsRCPartPhotons.h
   \-------------------------
   \brief   CsRCPartPhotons class declaration.
   \author  Paolo Schiavon
   \version 0.01
   \date    1 October 2000, rev. August 2005
*/

  #include <CLHEP/Vector/ThreeVector.h>

//-------------------------------------
  #include "CsHistograms.h"

  #include <vector>

  #include "CsRCRecConst.h"

//-------------------------
  class CsRCPhotonDet;

  class CsRCEventClusters;
  class CsRCCluster;

  class CsRCParticle;
  class CsRCRing;

  class CsRCPhoton;
//-------------------------


  class CsRCPartPhotons {


  public:


    CsRCPartPhotons();

    CsRCPartPhotons( int, CsRCParticle* );

    CsRCPartPhotons( const CsRCPartPhotons& );

    void clearPartPhotons();

    void doPartToMWR( CsRCParticle* );

    void doSelCathodes( CsRCParticle* );

    void getPartPhotons( CsRCParticle* );

    bool mcsPartCorr( CsRCParticle*, CLHEP::Hep3Vector&, CLHEP::Hep3Vector& );
    bool MIPReject( CsRCCluster* );
    bool killHalo( CsRCCluster* );

    bool getSplitCondition( CsRCParticle* );
    std::vector<double> getSplitLimits( CsRCParticle* );
    float yDetLLimit( const int );

    bool getSplitCondition();
    std::vector<double> getSplitLimits();
    inline bool bSplitLLimSet() const { return bSplitLLimSet_; };
    std::vector<double> xDetLimits( const int );
    std::vector<double> yDetLimits( const int );
    bool getHitCathodes( const double, int*, int*, double* );
    double nPhotExpctCat( double );

    double sigmaPhoRec( const double ) const;
    double sigmaPhot5( const double ) const;
    double sigmaPhot7( const double ) const;
    double sigmaPhot8( const double ) const;
    double getCorrMom( const double ) const;
    double getCorrFactor() const;

    long double getLikelihood( const double, int& );
    long double getLikeAll( const double, int& );
    long double getLikeAllMap( const double, int& );   // map 050204 OBSOLETE
    long double getLikelihoodAPV( const double, int& );
    long double getLikelihoodPMT( const double, int& );
    void getLikeDeriv( const int, const double, const double,
		       double&, double&,
		       double&, double&, int&, int& );

    double getQsquare( const double, int& );

    void compThetaIpo( CsRCParticle* );

    bool getThetaLikeMax( long double&, double& );
    bool getThetaLikeMax();
    void checkLikelihood( const double );
    void checkLikeRatios( const CsRCPhoton*,  const double,  const double );
    void checkLikeReso();
    void getLikeProb();
    int getRelMaxima( int, double*, int&, int* );
    bool getLikeProfile();
    bool histLikeProfile( std::vector<CsHist1D*>&, const int, const int );
    bool histQsqProfile( std::vector<CsHist1D*>&, const int, const int );
    double nPhotExpct( double );

    double getThetaWgav();

    inline  int kPaPhot() const { return kPaPhot_; };
    inline  CsRCParticle* pPart() const { return pPart_; };
    inline  CsRCRing* pRing() const { return pRing_; };
    void  setpRing( CsRCRing* ptr ) { pRing_ = ptr; };

    inline  const std::list<CsRCPhoton*> &lPhotons() const { 
      return lPhotons_; };
    inline  const std::list<CsRCCluster*> &lCluSignal() const {
      return lCluSignal_; };
    void setCluSignal( CsRCCluster* );

    inline  std::vector<CLHEP::Hep3Vector> vPoPartW() const { return vPoPartW_; };
    inline  std::vector<CLHEP::Hep3Vector> vDcPartW() const { return vDcPartW_; };
    inline  std::vector<CLHEP::Hep3Vector> vPoPhotW() const { return vPoPhotW_; };
    inline  std::vector<CLHEP::Hep3Vector> vPoPaMirW() const { return vPoPaMirW_; };
    inline  std::vector<CLHEP::Hep3Vector> vCorPoPhoW() const { return vCorPoPhoW_; };
    inline  std::vector<CLHEP::Hep3Vector> vDcPaReflW() const { return vDcPaReflW_; };
    inline  std::vector<CLHEP::Hep3Vector> vPoPaDetW() const { return vPoPaDetW_; };
    inline  std::vector<double> vySplitLLim() const { return vySplitLLim_; };

    inline  int kDetPart() const { return kDetPart_; };
    inline  CsRCPhotonDet* pDetPart() const { return pDetPart_; };
    inline  int iCaPa() const { return iCaPa_; };
    inline  int iXpPa() const { return iXpPa_; };
    inline  int iYpPa() const { return iYpPa_; };

    inline  const std::list<int> &lSelCats() const { return lSelCats_; };

    inline  double thetaLike() const { return thetaLike_; };
    inline  double pLikeMax() const { return pLikeMax_; };
    inline  bool thetaLikeSet() const { return thetaLikeSet_; };
    inline  void setThetaLike( const double thetaLike ) {
      thetaLike_ = thetaLike;
      thetaLikeSet_ = true;
    }
    inline  float fracUse() const { return fracUse_; };
    inline  void setFracUse( float frac ) { fracUse_ = frac; };
    inline  double normSg() const { return normSg_; };
    inline  double normBk() const { return normBk_; };

    inline  double* partProbs() { return partProbs_; };
    inline  bool partProbsSet() const { return partProbsSet_; };
    inline  void setPartProbs( double* paPro ) {
      int nProb = CsRCRecConst::Instance()->outBufferSize();
      for( int j=0; j<nProb; j++ ) partProbs_[j] = paPro[j];
      partProbsSet_ = true;
    };
    inline  void setPartProb( const int k, const double prob ) {
      partProbs_[k] = prob; }
    inline  void setPartProbsSet() { partProbsSet_ = true; };

    inline  double* probaLKAll() { return probaLKAll_; };
    inline  double* derivLKAll() { return derivLKAll_; };
    inline  double* probaLKAllAPV() { return probaLKAllAPV_; };
    inline  double* probaLKAllPMT() { return probaLKAllPMT_; };
    inline  double* derivLKAllAPV() { return derivLKAllAPV_; };
    inline  double* derivLKAllPMT() { return derivLKAllPMT_; };
    inline  void setProbaLKAll( double* proba ) {
      for( int j=0; j<31; j++ ) probaLKAll_[j] = proba[j]; };
    inline  void setDerivLKAll( double* deriv ) {
      for( int j=0; j<31; j++ ) derivLKAll_[j] = deriv[j]; };
    inline  void setProbaLKAllAPV( double* proba ) {
      for( int j=0; j<31; j++ ) probaLKAllAPV_[j] = proba[j]; };
    inline  void setProbaLKAllPMT( double* proba ) {
      for( int j=0; j<31; j++ ) probaLKAllPMT_[j] = proba[j]; };
    inline  void setDerivLKAllAPV( double* deriv ) {
      for( int j=0; j<31; j++ ) derivLKAllAPV_[j] = deriv[j]; };
    inline  void setDerivLKAllPMT( double* deriv ) {
      for( int j=0; j<31; j++ ) derivLKAllPMT_[j] = deriv[j]; };
    inline  double probaLKBgAll() { return probaLKBgAll_; };
    inline  double nPhotAll() { return nPhotAll_; };
    inline  double nPhotAllAPV() { return nPhotAllAPV_; };
    inline  double nPhotAllPMT() { return nPhotAllPMT_; };
    inline  void setProbaLKBgAll( const double proba ) { 
      probaLKBgAll_ = proba; };
    inline  void setNPhotAll( const int npho ) { nPhotAll_ = npho; };
    inline  void setNPhotAllAPV( const int npho ) { nPhotAllAPV_ = npho; };
    inline  void setNPhotAllPMT( const int npho ) { nPhotAllPMT_ = npho; };
    inline  double* probaChiAll() { return probaChiAll_; };
    inline  void setProbaChiAll( double* proba ) {
      for( int j=0; j<31; j++ ) probaChiAll_[j] = proba[j]; };

    inline  bool flag() const { return flag_; };
    inline  void setFlag( bool flag ) { flag_ = flag; };

    inline  int nLikePro() const { return nLikePro_; };
    inline  double theLikeProMn() const { return theLikeProMn_; };
    inline  double theLikeProMx() const { return theLikeProMx_; };
    inline  double* theLikePro() { return theLikePro_; };
    inline  double* pLikePro() { return pLikePro_; };
    inline  bool likeProSet() const { return likeProSet_; };

    CLHEP::Hep3Vector getCluAngle( const CLHEP::Hep3Vector&, const CLHEP::Hep3Vector&,
			    const float&, bool& );

    CLHEP::Hep3Vector doQzWCorr( CLHEP::Hep3Vector, CLHEP::Hep3Vector, CLHEP::Hep3Vector,
			  float, float );
    CLHEP::Hep3Vector doQzWCorr( const CLHEP::Hep3Vector, const CLHEP::Hep3Vector, CLHEP::Hep3Vector&,
			  const float, const float, bool& );
    CLHEP::Hep3Vector doQzWCorr( CsRCCluster*, const CLHEP::Hep3Vector, const CLHEP::Hep3Vector,
			  CLHEP::Hep3Vector&, const float, const float, bool& );

    CLHEP::Hep3Vector photImpDet( CLHEP::Hep3Vector, CLHEP::Hep3Vector, const double,
			   const float );

    CLHEP::Hep3Vector doPMTOptCorr( CsRCCluster*, const CLHEP::Hep3Vector, const CLHEP::Hep3Vector,
		             CLHEP::Hep3Vector&, const float, bool&, bool& );
    //double getPMTOptCorr( CsRCCluster*, const double, const double );
    
    inline  double thetaIpo( const int kPa ) const { return thetaIpo_[kPa]; };
    inline  double thetaIpoUV( const int kPa ) const {
      return thetaIpoUV_[kPa]; };
    inline  double thetaIpoVS( const int kPa ) const {
      return thetaIpoVS_[kPa]; };
    double howPMT();

    void setPhoTheNorm();

    double thetaUVtoVS( const double ) const;
    double thetaVStoUV( const double ) const;

    inline  bool likeONLY() const { return likeONLY_; };
    inline  void setLikeONLY( const bool on ) { likeONLY_ = on; };
    inline  bool likeFirst() const { return likeFirst_; };
    inline  void setLikeFirst( const bool first ) { likeFirst_ = first; };
    inline  bool pionONLY() const { return pionONLY_; };
    inline  void setPionONLY( const bool on ) { pionONLY_ = on; };
    inline  bool kaonONLY() const { return kaonONLY_; };
    inline  void setKaonONLY( const bool on ) { kaonONLY_ = on; };

    inline  bool likeAPVONLY() const { return likeAPVONLY_; };
    inline  void setLikeAPVONLY( const bool on ) { likeAPVONLY_ = on; };
    inline  bool likePMTONLY() const { return likePMTONLY_; };
    inline  void setLikePMTONLY( const bool on ) { likePMTONLY_ = on; };

    bool GetPMTPhotonPosition( int, int, double, double, double&, double& );
    bool GetPadCentre( int, double&, double& );

    double getThetaIpo( double, double, double );

    void print() const;

    ~CsRCPartPhotons();


    private:


    int kPaPhot_;
    CsRCParticle* pPart_;
    CsRCRing* pRing_;
    std::list<CsRCPhoton*> lPhotons_;
    std::list<CsRCCluster*> lCluSignal_;

    std::vector<CLHEP::Hep3Vector> vPoPartW_;
    std::vector<CLHEP::Hep3Vector> vDcPartW_;
    std::vector<CLHEP::Hep3Vector> vPoPhotW_;
    std::vector<CLHEP::Hep3Vector> vPoPaMirW_;
    std::vector<CLHEP::Hep3Vector> vCorPoPhoW_;
    std::vector<CLHEP::Hep3Vector> vDcPaReflW_;
    std::vector<CLHEP::Hep3Vector> vPoPaDetW_;

    std::vector<double> vySplitLLim_;
    bool bSplitLLimSet_;

    double thetaIpo_[31];
    double thetaIpoUV_[31];
    double thetaIpoVS_[31];

    int kDetPart_;
    CsRCPhotonDet* pDetPart_;
    int iCaPa_;
    int iXpPa_;
    int iYpPa_;
    std::list<int> lSelCats_;

    double thetaLike_;
    double pLikeMax_;
    bool thetaLikeSet_;
    float fracUse_;
    double normSg_;
    double normBk_;

    double partProbs_[50];
    bool partProbsSet_;
    double probaLKAll_[31];
    double derivLKAll_[31];
    double probaLKBgAll_;
    int nPhotAll_;
    double probaChiAll_[31];
    double probaLKAllAPV_[31];
    double probaLKAllPMT_[31];
    double derivLKAllAPV_[31];
    double derivLKAllPMT_[31];
    int nPhotAllAPV_;
    int nPhotAllPMT_;

    int nLikePro_;
    double theLikeProMn_;
    double theLikeProMx_;
    double theLikePro_[400];
    double pLikePro_[400];
    bool likeProSet_;

    bool likeONLY_;
    bool likeFirst_;
    bool pionONLY_;
    bool kaonONLY_;

    bool likeAPVONLY_;
    bool likePMTONLY_;

    bool flag_;

  };

#endif
