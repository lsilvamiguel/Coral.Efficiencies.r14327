#ifndef  ELLIPSEFITTEST_H
#define  ELLIPSEFITTEST_H

/*!
   \file    CsRCEllipseFitTest.h
   \-----------------------
   \brief   CsRCEllipseFit class declaration.
   \author  Paolo Schiavon
   \version 1.0
   \date    January  2003
*/



  class CsRCEllipseFitTest : public CsRCChiSqFit {


  public:


    CsRCEllipseFitTest();

    CsRCEllipseFitTest( const int, const double*, const double*,
			const double, const double,
			const int, const double*, const int* );

    CsRCEllipseFitTest( const int, const double*, const double*,
			const double*, const double*,
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

    virtual ~CsRCEllipseFitTest();


  private:


    CLHEP::HepVector convTolY_;
    CLHEP::HepVector convTolP_;


  };

#endif
