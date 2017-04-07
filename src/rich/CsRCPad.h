
#ifndef  RCPAD_H
#define  RCPAD_H

/*!
   \file    CsRCPad.h
   \-----------------
   \brief   CsRCPad class declaration.
   \author  Paolo Schiavon
   \version 0.01
   \date    October 2000, rev. September 2005
*/


  class CsDigit;


  class CsRCPad {


  public:

    CsRCPad();

    CsRCPad( int, CsDigit*, int, int, int, double );

    CsRCPad( int, int, int, int, double );

    CsRCPad( const CsRCPad& );

    inline  int kPad() const { return kPad_; };
    inline  CsDigit* pDigit() const { return pDigit_; };
    inline  int ic() const { return ic_; };
    inline  int ix() const { return ix_; };
    inline  int iy() const { return iy_; };
    inline  double PH() const { return PH_; };
    inline  double PHPack() const { return PHPack_; };
    inline  void setPHPack( const double PH ) { PHPack_ = PH; };
    inline  double PMTT0() const { return PMTT0_; };
    inline  void setPMTT0( float time ) { PMTT0_ = time; };
    inline  int PMTnum() const { return PMTnum_; };
    inline  int PMTcha() const { return PMTcha_; };
    bool setPMTpars();
    inline  float MCTime() const { return MCTime_; };
    inline  void setMCTime( float time ) { MCTime_ = time; };
    inline  bool flag() const { return flag_; };
    inline  void setFlag( bool flag ) { flag_ = flag; };
    inline  bool flagS() const { return flagS_; };
    inline  void setFlagS( bool flag ) { flagS_ = flag; };

    void print() const;

    ~CsRCPad();

  
  private:

    int kPad_;
    CsDigit*  pDigit_;
    int ic_, ix_, iy_;
    double PH_;
    double PHPack_;
    double PMTT0_;
    int PMTnum_;
    int PMTcha_;
    float MCTime_;
    bool flag_;
    bool flagS_;

  };

#endif
