  #ifndef  EVENTPARTICLES_H
  #define  EVENTPARTICLES_H

/*!
   \file    CsRCEventParticles.h
   \----------------------------
   \brief   CsRCEventParticles class declaration.
   \author  Paolo Schiavon
   \version 0.01
   \date    1 October 2000
*/


    #include <CLHEP/Vector/ThreeVector.h>
//---------------------------------------

    class CsTrack;
    class CsMCTrack;

    class CsRCEventSimul;

    #include "CsRCParticle.h"


    class CsRCEventParticles {


        public:

    static CsRCEventParticles* Instance();

    CsRCEventParticles( const CsRCEventParticles& );

    CsRCEventParticles& operator=( const CsRCEventParticles& );

    void clearEventParticles();

    void getEventParticles();

    void setMyParticle();
    bool initMyFile( int& );
    int openMyFile( char* );

    void setMCParticle();
    void setMCRecParticle();

    void setDataParticle();
    void setGeoCorrections( CLHEP::Hep3Vector&, CLHEP::Hep3Vector& );
    void setRunIndex( char * );

    inline  int nPart() const { return lParticles_.size(); };
    inline  const std::list<CsRCParticle*> &lParticles() const 
      { return lParticles_; }
    inline  bool flag() const { return flag_; };

    inline  bool flagSimul() const { return flagSimul_; };
    inline  int myFile() const { return myFile_; };
    inline  int myFileRun() const { return myFileRun_; };

    inline  void setFlag( bool flag ) { flag_ = flag; }

    inline  void printMem() const {
      std::cout << "  EventParticles pt  " << &lParticles_ << ",";
    }

    void partAnalysis();

    bool partFilter( CsRCParticle* );

    bool partSelect( CsRCParticle* );

    bool exitWindow( CsRCParticle* );
    double getCutXq( const double );
    double getCutYq( const double );
    double cirCut( double );
    double cirCu1( double );

    void partCorr( CsRCParticle* );

    void partResol( CsRCParticle* );

    void print() const;

    inline  std::vector<double> getPartVar() const { return partVar_; };

        protected:

      CsRCEventParticles();

      ~CsRCEventParticles();
  

        private:

    static CsRCEventParticles* instance_;

    std::list<CsRCParticle*> lParticles_;
    bool flag_;

    bool flagSimul_;

    int myFile_;
    int myFileRun_;

    std::vector<double> partVar_;

    };

  #endif
