#ifndef  CIRCTRACKMOMFIT_H
#define  CIRCTRACKMOMFIT_H

/*!
   \file    CsRCTrackMomFit.h
   \-----------------------
   \brief   CsRCTrackMomFit class declaration.
   \author  Paolo Schiavon
   \version 1.0
   \date    June 2005
*/



  class CsRCTrackMomFit {


  public:


    static CsRCTrackMomFit* Instance();

    CsRCTrackMomFit();

    bool test( CLHEP::Hep3Vector&, CLHEP::Hep3Vector&, double, double );

    bool CsRCTrackThruField( CLHEP::Hep3Vector, CLHEP::Hep3Vector,
			     double, float&, int&, int,
			     std::vector<CLHEP::Hep3Vector>&,
			     std::vector<CLHEP::Hep3Vector>&,
			     CLHEP::Hep3Vector&,  CLHEP::Hep3Vector& );

    bool CsRCTrackThruField( CLHEP::Hep3Vector, CLHEP::Hep3Vector,
			     double, float&,
			     CLHEP::Hep3Vector&,  CLHEP::Hep3Vector& );

    ~CsRCTrackMomFit();




  private:


    static CsRCTrackMomFit* instance_;


  };

#endif
