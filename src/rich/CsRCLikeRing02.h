
#ifndef  LIKERING02_H
#define  LIKERING02_H

/*!
   \file    CsRCLikeRing02.h
   \--------------------------
   \brief   CsRCLikeRing02 class declaration.
   \author  Paolo Schiavon
   \version 1.0
   \date    9 December 2002
*/


//-----------------------------
#include "CsRCLikeRing.h"
  class CsRCPhoton;
  class CsRCRing;
//-----------------------------



  class CsRCLikeRing02 : public CsRCLikeRing {


  public:


    CsRCLikeRing02( const CsRCRing* );

    double normSignal( const double );
    double normBackgr( const double, const double );
    double likeSignal( const double, const double, const double );
    double likeBackgr( CsRCPhoton*, const double );
    double getRingBackground( const double );

    virtual ~CsRCLikeRing02();


  private:


  };

#endif
