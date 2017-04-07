#ifndef  LIKEALL05_H
#define  LIKEALL05_H

/*!
   \file    CsRCLikeAll05.h
   \--------------------------
   \brief   CsRCLikeAll05 class declaration.
   \author  Paolo Schiavon
   \version 1.0
   \date    May 2004
*/


//-----------------------------
  #include "CsRCLikeAll.h"

  class CsRCPartPhotons;
  class CsRCPhoton;
//-----------------------------



  class CsRCLikeAll05 : public CsRCLikeAll {


  public:


    CsRCLikeAll05( const CsRCPartPhotons* );

    double normSignal( const double );
    double normBackgr( const double );
    double likeSignal( const CsRCPhoton*, const double );
    double likeBackgr( const CsRCPhoton*, const double );
    double getRingBackground( const double );

    inline bool splitPatt() const { return splitPatt_; };
    inline float fracUse() const { return fracUse_; };

    virtual ~CsRCLikeAll05();


  private:

    double xHole_;
    double yHole_;
    double rPipe_;

    bool splitPatt_;
    float nPhotEx_;
    float fracUse_;
    float fracUsePhi_[100];
    bool fracUsePhiSet_;
    double dPhi_;
    float normB_;
    double phoBack_[1000];
    bool phoBackSet_;

    bool dump_;
    bool useOldCode_;


  };

#endif
