
#ifndef  LIKEALL_H
#define  LIKEALL_H

/*!
   \file    CsRCLikeAll.h
   \---------------------
   \brief   CsRCLikeAll class declaration.
   \author  Paolo Schiavon
   \version 1.0
   \date    12 January 2003
*/


//-----------------------------
  class CsRCPartPhotons;
  class CsRCPhoton;
//-----------------------------



  class CsRCLikeAll {


  public:


    virtual double normSignal( const double );
    virtual double normBackgr( const double );
    virtual double likeSignal( const CsRCPhoton*, const double );
    virtual double likeBackgr( const CsRCPhoton*, const double );
    virtual double getRingBackground( const double );

    inline void setPPartPhot( CsRCPartPhotons* pp ) { pPartPhot_ = pp; };

    virtual ~CsRCLikeAll();


  protected:


    CsRCLikeAll();

    const CsRCPartPhotons* pPartPhot_;

  };

#endif
