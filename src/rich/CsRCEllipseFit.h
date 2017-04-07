#ifndef  ELLIPSEFIT_H
#define  ELLIPSEFIT_H

/*!
   \file    CsRCEllipseFit.h
   \-----------------------
   \brief   CsRCEllipseFit class declaration.
   \author  Paolo Schiavon
   \version 1.0
   \date    January  2003
*/



  class CsRCEllipseFit : public CsRCChiSqFit {


  public:


    CsRCEllipseFit();

    CsRCEllipseFit( const int, const double*, const double*,
		    const double, const double,
		    const int, const double*, const int* );

    CsRCEllipseFit( const int, const double*, const double*,
		    const double*, const double*,
		    const int, const double*, const int* );

    //CsRCEllipseFit( const CsRCEllipseFit& fit );

    bool doChiFit();
    bool fitStart();
    bool getDeriv();
    bool fitIteration();
    bool fitCheckConv();
    bool fitEnd();
    bool fitPulls();

    bool doHist( const std::vector<CsHist2D*> );

    void print() const;

    virtual ~CsRCEllipseFit();


  private:


    CLHEP::HepVector convTolY_;
    CLHEP::HepVector convTolP_;


  };

#endif
