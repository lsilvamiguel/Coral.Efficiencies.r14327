
#ifndef  LIKEALL04_H
#define  LIKEALL04_H

/*!
   \file    CsRCLikeAll04.h
   \--------------------------
   \brief   CsRCLikeAll04 class declaration.
   \author  Paolo Schiavon
   \version 1.0
   \date    May 2004
*/


//-----------------------------
  #include "CsRCLikeAll.h"

  class CsRCPartPhotons;
  class CsRCPhoton;
//-----------------------------



  class CsRCLikeAll04 : public CsRCLikeAll {


  public:


    CsRCLikeAll04( const CsRCPartPhotons* );

    double normSignal( const double );
    double normBackgr( const double );
    double likeSignal( const CsRCPhoton*, const double );
    double likeBackgr( const CsRCPhoton*, const double );

    virtual ~CsRCLikeAll04();


  private:


  };

#endif
