#ifndef  LIKEALL03_H
#define  LIKEALL03_H

/*!
   \file    CsRCLikeAll03.h
   \--------------------------
   \brief   CsRCLikeAll03 class declaration.
   \author  Paolo Schiavon
   \version 1.0
   \date    23 December 2002
*/


//-----------------------------
  #include "CsRCLikeAll.h"

  class CsRCPartPhotons;
  class CsRCPhoton;
//-----------------------------



  class CsRCLikeAll03 : public CsRCLikeAll {


  public:


    CsRCLikeAll03( const CsRCPartPhotons* );


    double normSignal( const double );
    double normBackgr( const double );
    double likeSignal( const CsRCPhoton*, const double );
    double likeBackgr( const CsRCPhoton*, const double );
    double getRingBackground( const double );

    virtual ~CsRCLikeAll03();


  private:


  };

#endif
