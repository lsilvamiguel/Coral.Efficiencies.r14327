#ifndef  EVENTPARTPHOTONS_H
#define  EVENTPARTPHOTONS_H

/*!
   \file    CsRCEventPartPhotons.h
   \------------------------------
   \brief   CsRCEventPartPhotons class declaration.
   \author  Paolo Schiavon
   \version 0.01
   \date    1 October 2000
*/


//----------------------------
  class CsRCPartPhotons;
//----------------------------



  class CsRCEventPartPhotons {


  public:


    CsRCEventPartPhotons();

    static CsRCEventPartPhotons* Instance();

    CsRCEventPartPhotons( const CsRCEventPartPhotons& );

    void clearEventPartPhotons();

    void getEventPartPhotons();

    void moniCFRefInd();

    void bkgrPhotons();

    void rawSignal();

    void partPhoSelection();

    void checkSignal();
    void checkAMPSCorr();
    void checkLikeCorr();

    void setBackgrType();
    bool getBackgrParam02();
    bool getBackgrParam03();
    bool getBackgrParam04();
    void partAllIdent();
    void setProbToTrack();

    inline  const std::list<CsRCPartPhotons*> &lPartPhotons() const
      { return lPartPhotons_; };

    inline  bool flag() const { return flag_; };
    inline  bool flagPartPhoSel() const { return flagPartPhoSel_; };

    inline std::vector<float> vPara() const { return vPara_; };

    inline  void setFlag( bool flag ) { flag_ = flag; };

    void print() const;


    ~CsRCEventPartPhotons();


  private:


    static CsRCEventPartPhotons* instance_;

    std::list<CsRCPartPhotons*> lPartPhotons_;

    bool flag_;

    bool flagPartPhoSel_;

    std::vector<float> vPara_;

  };

#endif
