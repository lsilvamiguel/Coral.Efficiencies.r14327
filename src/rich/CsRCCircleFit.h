#ifndef  CIRCLEFIT_H
#define  CIRCLEFIT_H

/*!
   \file    CsRCCircleFit.h
   \-----------------------
   \brief   CsRCCircleFit class declaration.
   \author  Paolo Schiavon
   \version 1.0
   \date    May 2000
*/



  class CsRCCircleFit : public CsRCChiSqFit  {


  public:


    CsRCCircleFit();

    CsRCCircleFit( const int, const double*, const double*,
		   const double, const int, const double*, const int* );

    CsRCCircleFit( const int, const double*, const double*, const double*,
		   const int, const double*, const int* );

    bool doChiFit();
    bool fitStart();
    bool getDeriv();
    bool fitIteration();
    bool fitCheckConv();
    bool fitEnd();
    bool fitPulls();

    bool doHist( const std::vector<CsHist2D*> );

    void print() const;

    virtual ~CsRCCircleFit();


  private:


    CLHEP::HepVector convTol_;


  };

#endif
