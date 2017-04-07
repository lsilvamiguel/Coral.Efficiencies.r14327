
#ifndef  LIKERING04_H
#define  LIKERING04_H

/*!
   \file    CsRCLikeRing04.h
   \--------------------------
   \brief   CsRCLikeRing03 class declaration.
   \author  Paolo Schiavon
   \version 1.0
   \date    16 October 2003
*/


//-----------------------------
#include "CsRCLikeRing.h"
  class CsRCPhoton;
  class CsRCRing;
//-----------------------------



  class CsRCLikeRing04 : public CsRCLikeRing {


  public:


    CsRCLikeRing04( const CsRCRing* );

    double normSignal( const double );
    double normBackgr( const double, const double );
    double likeSignal( const double, const double, const double );
    double likeBackgr( CsRCPhoton*, const double );
    double getRingBackground( const double );

    virtual ~CsRCLikeRing04();


  private:


  };

#endif
