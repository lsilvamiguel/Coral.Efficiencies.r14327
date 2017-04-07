
#ifndef  LIKEALL02_H
#define  LIKEALL02_H

/*!
   \file    CsRCLikeAll02.h
   \--------------------------
   \brief   CsRCLikeAll02 class declaration.
   \author  Paolo Schiavon
   \version 1.0
   \date    23 December 2002
*/


//-----------------------------
  #include "CsRCLikeAll.h"

  class CsRCPartPhotons;
  class CsRCPhoton;
//-----------------------------



  class CsRCLikeAll02 : public CsRCLikeAll {


  public:


    CsRCLikeAll02( const CsRCPartPhotons* );

    double normSignal( const double );
    double normBackgr( const double );
    double likeSignal( const CsRCPhoton*, const double );
    double likeBackgr( const CsRCPhoton*, const double );

    virtual ~CsRCLikeAll02();


  private:


  };

#endif
