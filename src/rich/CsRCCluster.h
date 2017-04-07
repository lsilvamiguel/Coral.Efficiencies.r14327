
#ifndef  RCCLUSTER_H
#define  RCCLUSTER_H

/*!
   \file    CsRCCluster.h
   \----------------------
   \brief   CsRCCluster class declaration.
   \author  Paolo Schiavon
   \version 0.01
   \date    October 2000, rev. August 2005
*/


//---------------------------
  class CsRCPad;
//---------------------------

  #include <list>


  class CsRCCluster {


  public:


    CsRCCluster();

    CsRCCluster( int, std::list<CsRCPad*>, int );

    CsRCCluster( int, std::list<CsRCPad*>, int, double,  double, double );

    CsRCCluster( const CsRCCluster& );

    inline  int kClu() const { return kClu_; };
    inline  const std::list<CsRCPad*> &lPads() const { return lPads_; };
    inline  int ic() const { return ic_; };
    inline  double xc() const { return xc_; };
    inline  double yc() const { return yc_; };
    inline  double PH() const { return PH_; };
    inline  bool flag() const { return flag_; };

    inline  void setFlag( bool flag ) { flag_ = flag; };

    inline void setPadList( std::list<CsRCPad*> lPads ) { lPads_ = lPads; };
    inline void pushPadList( CsRCPad* ppad ) {
      lPads_.push_back( ppad );
    }
    inline void setX( double xc ) { xc_ = xc; };
    inline void setY( double yc ) { yc_ = yc; };
    inline void setPH( double PH ) { PH_ = PH; };

    inline  bool isPMT() const { return isPMT_; };
    inline  bool isAPV() const { return isAPV_; };

    bool getBackWgt( double &back );
    bool getBackWgtX( double &back );

    void print() const;

    ~CsRCCluster();

  
  private:


    int kClu_;
    std::list<CsRCPad*> lPads_;
    int ic_;
    double xc_, yc_;
    double PH_;

    bool flag_;

    bool isPMT_, isAPV_;

  };

#endif
