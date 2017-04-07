
#ifndef  LIKERING03_H
#define  LIKERING03_H

/*!
   \file    CsRCLikeRing03.h
   \--------------------------
   \brief   CsRCLikeRing03 class declaration.
   \author  Paolo Schiavon
   \version 1.0
   \date    23 December 2002
*/


//-----------------------------
#include "CsRCLikeRing.h"
  class CsRCPhoton;
  class CsRCRing;
//-----------------------------



  class CsRCLikeRing03 : public CsRCLikeRing {


  public:


    CsRCLikeRing03( const CsRCRing* );

    double normSignal( const double );
    double normBackgr( const double, const double );
    double likeSignal( const double, const double, const double );
    double likeBackgr( CsRCPhoton*, const double );
    double getRingBackground( const double );

    virtual ~CsRCLikeRing03();


  private:


  };

#endif
