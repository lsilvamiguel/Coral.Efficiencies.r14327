// $Id: CsMCParticle.h,v 1.2 2000/04/07 12:40:18 benigno Exp $

/*!
   \file    CsMCParticle.h
   \brief   Compass Monte Carlo particle Class.
   \author  Benigno Gobbo
   \version $Revision: 1.2 $
   \date    $Date: 2000/04/07 12:40:18 $
*/

#ifndef CsMCParticle_h
#define CsMCParticle_h

/*! \class CsMCParticle 
    \brief   Compass Monte Carlo particle Class.
*/

class CsMCParticle {

 public:

  /*! \fn CsMCParticle(); 
    \brief Default Constructor.
    It creates a photon.
  */
  CsMCParticle();  

  /*! \fn CsMCParticle( int gnumber, int isGeantNumber=1 );
    \brief Constructor
    \param gnumber Geant particle number
    \param isGeantNumber (1) if number is Geant ParticleId or (0) if PDG-ID
  */
  CsMCParticle( int gnumber,  int isGeantNumber=1 );

  /*! \fn CsMCParticle( char pname[], int pcharge );
    \brief Constructor
    \param pname PDG particle name
    \param pcharge particle charge
  */
  CsMCParticle( const char* pname, int pcharge );

  /*! \fn CsMCParticle( int pnumber, bool antiparticle );
    \brief Constructor
    \param pnumber PDG particle number
    \param antiparticle if true create the antiparticle
  */
  CsMCParticle( int pnumber, bool antiparticle );

  CsMCParticle( const CsMCParticle& ); //!< Copy Costructor

  bool operator==( const CsMCParticle& ) const; //!< Assigment operator

  /*! \fn   char*  getName();
    \brief Returns PDG the particle name
  */
  char*  getName();

  /*! \fn   int    getCharge();
    \brief Returns the particle charge
  */
  int    getCharge();

  /*! \fn   double getMass();
    \brief Returns the particle mass (GeV/c^2)
  */
  double getMass();

  /*! \fn   double getMassErrP();
    \brief Returns the positive error on particle mass
  */
  double getMassErrP();

  /*! \fn  double getMassErrN();
    \brief Returns the negative error on particle mass
  */
  double getMassErrN();

  /*! \fn  double getWidth();
    \brief Returns the particle width (GeV/c^2)
  */
  double getWidth();

  /*! \fn  double getWidthErrP();
    \brief Returns the negative error on particle width
  */
  double getWidthErrP();

  /*! \fn  double getWidthErrN();
    \brief Returns the negative error on particle width
  */
  double getWidthErrN();

  /*! \fn  int getNumber();
    \brief Returns the PDG particle number
  */
  int getNumber();

  /*! \fn int getGeantNumber();
    \brief Returns the Geant particle number
  */
  int getGeantNumber();

 private:
  char   name_[20];      //!< PDG particle name
  int    charge_;        //!< particle charge
  double mass_;          //!< particle mass (GeV/c^2)
  double massErr_[2];    //!< +/- errors on particle mass 
  double width_;         //!< particle width 9GeV/c^2)
  double widthErr_[2];   //!< +/- errors on particle mass 
  int    number_;        //!< PDG particle number
  int    GeantNumber_;   //!< Geant particle number
  bool   antiparticle_;  //!< true if antiparticle

  /* \fn int Gpar2PDGpar( int Gpar );
     \brief Converts from Geant Particles numbering to PDG particles 
     numbering convention.

     NOTE: (Hope) a temporary method.
  */
  int Gpar2PDGpar( int Gpar );
  
  /* this is even worse then the above, but we need to keep it compatible */
  int PDGpar2GPar2( int PDGpar);
};

#endif // CsMCParticle_h
