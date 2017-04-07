#ifndef  DOCHIFIT_H
#define  DOCHIFIT_H

/*!
   \file    CsRCChiSqFit.h
   \----------------------
   \brief   CsRCChiSqFit class declaration.
   \author  Paolo Schiavon
   \version 1.0
   \date    May 2000
*/


  #include "CLHEP/Matrix/Matrix.h"
  #include "CLHEP/Matrix/SymMatrix.h"
  #include "CLHEP/Matrix/DiagMatrix.h"
  #include "CLHEP/Matrix/Vector.h"

//------------------------------------

  #include "CsHistograms.h"


  class CsRCChiSqFit {


  public:


    CsRCChiSqFit();
    CsRCChiSqFit( const int, const double*, const double*, const int,
		  const double*, const int* );
    CsRCChiSqFit( const int, const double*, const int );
    CsRCChiSqFit( const int, const double*, const double*, const int,
		  const int, const double*, const int* );

    virtual void getCovMax( const double );
    virtual void getCovMax( const double* );
    virtual void getCovMax( const double, const double );
    virtual void getCovMax( const double*, const double* );
    virtual void getCovMax( const CLHEP::HepMatrix& );

    virtual bool doChiFit();

    virtual bool fitStart();
    virtual bool getDeriv();
    virtual bool fitIteration();
    virtual bool fitCheckConv();
    virtual bool fitEnd();
    virtual bool fitPulls();

    virtual bool doHist( const std::vector<CsHist2D*> );

    inline double chiSquare() const { return chiSquare_; };
    inline int nIter() const { return nIter_; };
    inline bool flag() const { return flagFit_; };
    inline int degFree() const { return nDegFree_; };
    inline CLHEP::HepVector fitted() const { return Yf_; };
    inline CLHEP::HepMatrix covYf() const { return Erf_; }
    inline CLHEP::HepVector para() const { return Pf_; };
    inline CLHEP::HepMatrix covPf() const { return Epf_; }
    inline CLHEP::HepVector pull() const { return pully_; }
    inline void setFlag( const bool boo ) { flagFit_ = boo; };

    void print() const;

    virtual ~CsRCChiSqFit();


  protected:


    double chiSquare_;
    int nDegFree_;
    int nIter_;

    int nPoint_;
    CLHEP::HepVector Ym_, Ya_, Yf_;
    CLHEP::HepVector dYm_, dYw_, dYd_;
    CLHEP::HepVector Xm_;
    CLHEP::HepVector pully_;

    CLHEP::HepMatrix Erm_, Wgm_, Erf_;

    int nConst_;
    CLHEP::HepVector Co_, Ro_, Lk_;
    CLHEP::HepMatrix Eco_;

    CLHEP::HepMatrix mA_;
    CLHEP::HepMatrix mB_;

    int nParam_;
    CLHEP::HepVector PaIn_, Pa_, iPaF_, Pf_;
    CLHEP::HepMatrix Epf_;
    int nPaF_;
    CLHEP::HepVector Xp_;
    CLHEP::HepVector dXp_;

    int kFail_;
    bool flagFit_;

    bool printKey_ ;

  };

#endif
