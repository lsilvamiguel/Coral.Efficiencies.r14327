  #ifndef  RICHONE_H
  #define  RICHONE_H

/*!
   \file    CsRichOne.h
   \--------------------------
   \brief   CsRichOne class declaration.
   \author  Paolo Schiavon
   \version 0.02
   \date    December 1999, rev. June 2000
*/


  #include <iostream>
  #include <ostream>

  #include "CsRegistry.h"
  #include "CsEndOfJob.h"

  class CsTrack;
  class CsMCTrack;


  class CsRichOne : public CsEndOfJob  {


      public:

  static CsRichOne* Instance();

  CsRichOne( const CsRichOne& );


  inline  int kEvent() const { return kEvent_; };
  inline  bool skipMyEvent() const { return skipMyEvent_; };
  inline  bool readMyEvent();
  inline  int kMyWrEvent() const { return kMyWrEvent_; };
  inline  bool flag() const { return flag_; };
  inline  bool stopMyFileJob() const { return stopMyFileJob_; };
  inline  void setStopMyFileJob( const std::string name ) {
    stopMyFileJob_ = true;
    std::cout << std::endl;
    std::cout << "Stop of MyFileJob requested by " << name << std::endl;
    std::cout << "------------------------------" << std::endl;
  };
  inline  bool UpRICHJob() const { return UpRICHJob_; };
  inline  void setUpRICHJob() { UpRICHJob_ = true; };
  inline  bool endMyFile() const { return endMyFile_; };
  inline  void setEndMyFile( const bool flag ) { endMyFile_ = flag; };

  inline  double richTime() const { return richTime_; };
  inline  double partsTime() const { return partsTime_; };
  inline  double padsTime() const { return padsTime_; };
  inline  double clusTime() const { return clusTime_; };
  inline  double paphsTime() const { return paphsTime_; };
  inline  double ringsTime() const { return ringsTime_; };
  inline  double anasyTime() const { return anasyTime_; };
  inline  double dispTime() const { return dispTime_; };

  inline  int evGood() const { return evGood_; };
  inline  int kMuonPart() const { return kMuonPart_; };
  inline  float zVertex() const { return zVertex_; };
  inline  float phiMass() const { return phiMass_; };

  void doRichOne();

  void MCRecEffic();
  void dataRecEffic();
  void execTime();

  bool getPartProbs( CsTrack*, double* );
  bool getPartProbs( CsMCTrack*, double* );

  void print() const;
  void printMemory( int ) const;

  bool end();

  void testMCTracks();

  void writeMyFile();
  void closeMyFile();

  void setFullPrintOn( const int );

  void checkOutput();
 
  void physSelec();
  void physAppend();


      protected:

  CsRichOne();

  ~CsRichOne();
  

      private:

  static CsRichOne* instance_;

  int kEvent_;
  bool skipMyEvent_;
  int kMyWrEvent_;
  bool stopMyFileJob_;
  bool UpRICHJob_;
  bool endMyFile_;

  float *probPart_;
  bool flag_;

  double richTime_;
  double partsTime_;
  double padsTime_;
  double clusTime_;
  double paphsTime_;
  double ringsTime_;
  double anasyTime_;
  double dispTime_;

//Physics
  int evGood_;
  int kMuonPart_;
  float zVertex_;
  float phiMass_;

  };

  #endif
