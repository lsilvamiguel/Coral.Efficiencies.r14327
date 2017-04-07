
#ifndef  LIKERING_H
#define  LIKERING_H

/*!
   \file    CsRCLikeRing.h
   \---------------------
   \brief   CsRCLikeRing class declaration.
   \author  Paolo Schiavon
   \version 1.0
   \date    12 January 2003
*/


//-----------------------------
  class CsRCPartPhotons;
  class CsRCPhoton;
//-----------------------------



  class CsRCLikeRing {


  public:


    virtual double normSignal( const double );
    virtual double normBackgr( const double, const double );
    virtual double likeSignal( const double, const double, const double );
    virtual double likeBackgr( CsRCPhoton*, const double );
    virtual double getRingBackground( const double );

    virtual ~CsRCLikeRing();


  protected:


    CsRCLikeRing();

    const CsRCPartPhotons* pPartPhot_;
    double theUpL_;
    double theLoL_;

  };

#endif
