
#ifndef  LIKERING05_H
#define  LIKERING05_H

/*!
   \file    CsRCLikeRing05.h
   \--------------------------
   \brief   CsRCLikeRing05 class declaration.
   \author  Paolo Schiavon
   \version 1.0
   \date    May 2004
*/


//-----------------------------
#include "CsRCLikeRing.h"

  class CsRCPhoton;
  class CsRCRing;
//-----------------------------



  class CsRCLikeRing05 : public CsRCLikeRing {


  public:


    CsRCLikeRing05( const CsRCRing* );

    double normSignal( const double );
    double normBackgr( const double, const double );
    double likeSignal( const double, const double, const double );
    double likeBackgr( CsRCPhoton*, const double );
    double getRingBackground( const double );

    virtual ~CsRCLikeRing05();


  private:


  };

#endif
