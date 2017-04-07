#ifndef  GAUPOLFIT_H
#define  GAUPOLFIT_H

/*!
   \file    CsRCGauPolFit.h
   \-----------------------
   \brief   CsRCGauPolFit class declaration.
   \author  Paolo Schiavon
   \version 1.0
   \date    December 2003
*/



  class CsRCGauPolFit : public CsRCChiSqFit  {


  public:


    CsRCGauPolFit();

    CsRCGauPolFit( const int, const double*, const double*,
		   const int, const double*, const int* );

    CsRCGauPolFit( const int, const double*, const double*, const double*,
		   const int, const double*, const int* );

    CsRCGauPolFit( const int, const double*, const double*, const double*,
		   const int, const double*, const double*, const int* );

    void getConvTol( const double* );

    bool doChiFit();
    bool fitStart();
    bool getDeriv();
    bool fitIteration();
    bool fitCheckConv();
    bool fitEnd();
    bool fitPulls();

    bool doHist( const std::vector<CsHist2D*> );

    void print() const;

    virtual ~CsRCGauPolFit();


  private:


    CLHEP::HepVector convTol_;
    bool convTolSet_;


  };

#endif
