#ifndef  EVENTRINGS_H
#define  EVENTRINGS_H

/*!
   \file    CsRCEventRings.h
   \------------------------
   \brief   CsRCEventRings class declaration.
   \author  Paolo Schiavon
   \version 0.03,  rev. 1/10/00
   \date    23 June 1999
*/


  #include <list>

// ----------------------------
  class CsRCPartPhotons;
  class CsRCRing;
// ----------------------------


  class CsRCEventRings {


  public:


    static CsRCEventRings* Instance();

    CsRCEventRings( const CsRCEventRings& );

    void clearEventRings();

    void getEventRings();

    void flagPhoSignal();

    void ringSelection();

    void singlePhoton();
    void singlePhotonPMT();
    void checkSignal();
    void checkAMPSCorr();

    void checkNoise();
    void checkRingPH();

    void checkPhotonAngle();
    void checkPMTTimeDif();

    void ringSignal();
    void ringSignal( CsRCRing* );

    void setBackgrType();
    void partChiSqIdent();
    void partRingAllIdent();
    void partRingIdent();

    void checkEventSignal();

    void checkRichMom();
    void checkPhiTheta();
    void checkRingFit();

    void testPhiKzero();

    inline  const std::list<CsRCRing*> &lRings() const { return lRings_; };

    inline  bool flag() const { return flag_; }
    inline  bool flagRingSel() const { return flagRingSel_; };

    inline  void setFlag( bool flag ) { flag_ = flag; }

    void print() const;
    inline  void printMem() const {
      std::cout << "  EventRings pt  " << &lRings_ << ",";
    }


  protected:


    CsRCEventRings();

    ~CsRCEventRings();


  private:


    static CsRCEventRings* instance_;

    std::list<CsRCRing*> lRings_;

    bool flag_;

    bool flagRingSel_;

  };

#endif
