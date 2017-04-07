  #ifndef  RCPARTICLE_H
  #define  RCPARTICLE_H

/*!
   \file    CsRCParticle.h
   \------------------------
   \brief   CsRCParticle class declaration.
   \author  Paolo Schiavon
   \version 1.0
   \date    1 October 2000
*/

# include <sstream>
# include <string>
# include <list>
# include <vector>

# include <CLHEP/Vector/ThreeVector.h>

//---------------------------------------

    class CsRCMirrorNom;
    class CsTrack;
    class CsMCTrack;
    class CsRCRing;


    class CsRCParticle  {


        public:

    CsRCParticle();

    CsRCParticle( int, CLHEP::Hep3Vector, CLHEP::Hep3Vector, CLHEP::Hep3Vector, double,
		  int, int, int );

    CsRCParticle( int, CLHEP::Hep3Vector, CLHEP::Hep3Vector, CLHEP::Hep3Vector,
		  CLHEP::Hep3Vector, CLHEP::Hep3Vector, CLHEP::Hep3Vector, CLHEP::Hep3Vector,
		  double, int, double, int, int, int );

    CsRCParticle( int, CsMCTrack*, CLHEP::Hep3Vector, CLHEP::Hep3Vector, CLHEP::Hep3Vector,
		  CLHEP::Hep3Vector, CLHEP::Hep3Vector, CLHEP::Hep3Vector, CLHEP::Hep3Vector,
		  double, int, double );

    CsRCParticle( int, CsTrack*, CsMCTrack*,
		  CLHEP::Hep3Vector, CLHEP::Hep3Vector, CLHEP::Hep3Vector,
		  CLHEP::Hep3Vector, CLHEP::Hep3Vector, CLHEP::Hep3Vector, CLHEP::Hep3Vector,
		  double, int, double );

    CsRCParticle( int, CsTrack*, CLHEP::Hep3Vector, CLHEP::Hep3Vector, CLHEP::Hep3Vector,
		  CLHEP::Hep3Vector, CLHEP::Hep3Vector, CLHEP::Hep3Vector, CLHEP::Hep3Vector,
		  double, int );

    CsRCParticle( const CsRCParticle& );

    inline  int kPart() const { return kPart_; };
    inline  CsTrack* pTrack() const { return pTrack_; };
    inline  CLHEP::Hep3Vector vPosIn() const { return vPosIn_; };
    inline  CLHEP::Hep3Vector vDirIn() const { return vDirIn_; };
    inline  CLHEP::Hep3Vector vPosOut() const { return vPosOut_; };
    inline  CLHEP::Hep3Vector vPosEmP() const { return vPosEmP_; };
    inline  CLHEP::Hep3Vector vDirEmP() const { return vDirEmP_; };
    inline  CLHEP::Hep3Vector vPosExW() const { return vPosExW_; };
    inline  CLHEP::Hep3Vector vDirExW() const { return vDirExW_; };
    inline  double xa() const { return vPosIn_.x(); };
    inline  double ya() const { return vPosIn_.y(); };
    inline  double za() const { return vPosIn_.z(); };
    inline  double ld() const { return vDirIn_.x(); };
    inline  double md() const { return vDirIn_.y(); };
    inline  double nd() const { return vDirIn_.z(); };
    inline  double xo() const { return vPosOut_.x(); };
    inline  double yo() const { return vPosOut_.y(); };
    inline  double zo() const { return vPosOut_.z(); };
    inline  double mom() const { return mom_; };
    inline  int charge() const { return charge_; };
    inline  bool flag() const { return flag_; };
    inline  CsRCRing* pRing() const { return pRing_; };
    void setpRing( CsRCRing* ptr ) { pRing_ = ptr; } ;

    inline  CsMCTrack* pMCTrack() const { return pMCTrack_; };
    inline  int iTrack() const { return iTrack_; };
    inline  int iPartT() const { return iPartTy_; };
    inline  int ioExWd() const { return ioExWd_; };
    inline  float MCmTime() const { return MCmTime_; };
    inline  void setMCmTime( float time )  { MCmTime_ = time; };
    inline  bool flagS() const { return flagS_; };
    inline  int nPhoCer() const { return nPhoCer_; };
    inline  double thetaCer() const { return thetaCer_; };

    inline  float zHelix0() const { return zHelix0_; };
    inline  float zHelix1() const { return zHelix1_; };
    inline  int nClus() const { return nClus_; };
    inline  float chiSq() const { return chiSq_; };

    inline  double pathLen() const { return pathLen_; };
    inline  double thPamir() const { return thPamir_; };
    inline  std::vector<CLHEP::Hep3Vector> vPade() const { return vPade_; };
    inline  std::vector<double> ddPaDet() const { return ddPaDet_; };
    inline  std::vector<double> thPade() const { return thPade_; };
    inline  int kDetPart() const { return kDetPart_; };
    inline  CsRCMirrorNom* pMirPart() const { return pMirPart_; };
    double mass( int );
    double paMass();

    inline  void setFlag( bool flag ) { flag_ = flag; };
    inline  void setFlagS( bool flag ) { flagS_ = flag; };

    inline  void setPathLen( double pathLen ) { pathLen_ = pathLen; };
    inline  void setThPamir( double thPamir ) { thPamir_ = thPamir; };
    inline  void setDirIn( CLHEP::Hep3Vector vDir ) { vDirIn_ = vDir; }
    inline  void setvPade( std::vector<CLHEP::Hep3Vector> vPade )  { vPade_ = vPade; };
    inline  void setDdPaDet( const std::vector<double> ddPaDet )
      { ddPaDet_ = ddPaDet; };
    inline  void setThPade( std::vector<double> thPade ) { thPade_ = thPade; };
    inline  void setDetPart( int kDetPart )  { kDetPart_ = kDetPart; };
    inline  void setpMirPart( CsRCMirrorNom* pMirPart ) {
      pMirPart_ = pMirPart; }

    inline  std::vector<CLHEP::Hep3Vector> vPoPaMir0() const { return vPoPaMir0_; };
    inline  void setPoPaMir0( std::vector<CLHEP::Hep3Vector> vPoPaMir0 ) { 
      vPoPaMir0_ = vPoPaMir0;
    };

    inline void setQualityPars( const float& zHlx0,  const float& zHlx1, 
        			const int& nClus,  const float& chi ) {
      zHelix0_ = zHlx0;
      zHelix1_ = zHlx1;
      nClus_ = nClus;
      chiSq_ = chi;
    }

    int getClosestCat() const;

    //as from 101220
    inline  void setVCorPoPa0( const CLHEP::Hep3Vector vCorPoPa0 )
      { vCorPoPa0_ = vCorPoPa0; };
    inline  CLHEP::Hep3Vector vCorPoPa0() const { return vCorPoPa0_; };
    inline  void setCorPaRR( const double corPaRR ) { corPaRR_ = corPaRR; }
    inline  double corPaRR() const { return corPaRR_; }

    void print() const;
    void printMC() const;

    ~CsRCParticle();
  

        private:

    int kPart_;
    CsTrack* pTrack_;
    CLHEP::Hep3Vector vPosIn_;
    CLHEP::Hep3Vector vDirIn_;
    CLHEP::Hep3Vector vPosOut_;
    CLHEP::Hep3Vector vPosEmP_;
    CLHEP::Hep3Vector vDirEmP_;
    CLHEP::Hep3Vector vPosExW_;
    CLHEP::Hep3Vector vDirExW_;
    double mom_;
    int charge_;
    bool flag_;

    double pathLen_;
    double thPamir_;
    std::vector<CLHEP::Hep3Vector> vPade_;
    std::vector<double> ddPaDet_;
    std::vector<double> thPade_;
    int kDetPart_;
    CsRCMirrorNom* pMirPart_;
    CsRCRing* pRing_;

    std::vector<CLHEP::Hep3Vector> vPoPaMir0_;

    //as from 101220
    CLHEP::Hep3Vector vCorPoPa0_;
    double corPaRR_;

    CsMCTrack* pMCTrack_;
    int iTrack_;
    int iPartTy_;
    int ioExWd_;
    float MCmTime_;
    bool flagS_;

    int nPhoCer_;
    double thetaCer_;

    //as from 020322
    float zHelix0_;
    float zHelix1_;
    int nClus_;
    float chiSq_;

    };

  #endif
