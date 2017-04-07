
#ifndef  LIKEALLMAP_H
#define  LIKEALLMAP_H

/*!
   \file    CsRCLikeAllMap.h
   \--------------------------
   \brief   CsRCLikeAllMap class declaration.
   \author  Stefano Panebianco
   \version 1.0
   \date    19 March 2004
*/


//-----------------------------
  #include "CsRCLikeAll.h"

  class CsRCPartPhotons;
  class CsRCCluster;
  
//-----------------------------


  class CsRCLikeAllMap : public CsRCLikeAll {


  public:

    CsRCLikeAllMap( const CsRCPartPhotons* );

    double normSignal( const double, double );
    double normSignal( const double var)  { return CsRCLikeAll::normSignal(var); }
    double normBackgr(std::list<CsRCPhoton*, std::allocator<CsRCPhoton*> >&);
    double normBackgr( const double var)  { return CsRCLikeAll::normBackgr(var); }
    double likeSignal( const double, const double, const double );
    double likeSignal( const CsRCPhoton* var1, const double var2)  { return CsRCLikeAll::likeSignal(var1,var2); }
    double likeBackgr( const double, const double );
    double likeBackgr( const CsRCPhoton* var1, const double var2)  { return CsRCLikeAll::likeBackgr(var1,var2); }
    virtual double normSignalCorr(const double x0, const double y0, const double tx0, const double ty0, 
const double theta, const double xph, const double yph); 
    virtual double likeSignalCorr( double x0,  double y0,  double tx0,  double ty0, 
     double theta,  double phi, double xph, double yph, const int top);
    virtual double thetaphitoxy(const double x0, const double y0, const double tx0, const double ty0, 
 double z, double theta, double phi, double &xph, double &yph);


    virtual ~CsRCLikeAllMap();


  private:


  };

#endif
