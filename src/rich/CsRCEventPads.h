
#ifndef  EVENTPADS_H
#define  EVENTPADS_H

/*!
   \file    CsRCEventPads.h
   \-----------------------
   \brief   CsRCEventPads class declaration.
   \author  Paolo Schiavon
   \version 0.03
   \date    October 2000, rev, August 2005
*/

// All active pads of the RICH photon chambers.
// --------------------------------------------


  #include <list>

  #include "CsHist.h"


  class CsRCPad;
  class CsDigit;


  class CsRCEventPads {


  public:

    static CsRCEventPads* Instance();

    CsRCEventPads( const CsRCEventPads& );
    CsRCEventPads& operator=( const CsRCEventPads& );

    void clearEventPads();
    void clearEventPads( std::list<CsRCPad*> &lPads );
    inline  void clearPadList() { lPads_.clear(); };
    inline  void copyPadList( const std::list<CsRCPad*> &lPads ) { 
      lPads_ = lPads; };

    bool getEventPads();

    bool setMyPads();
    bool setMCRecPads();
    bool setDataPads();

    bool checkPMTPPad( const double );
    bool checkAPVPad( const double, const double );

    void checkPPads();
    void swapPMTChs( int& x, int& y );

    bool setThreshold( int&, int&, int&, double& );

    bool killHalo( int&, int&, int&, double& );

    inline  const std::list<CsRCPad*> &lPads() const { return lPads_; };

    inline  bool flag() const { return flag_; };

    void print() const;
    void print( const int& ) const;

    inline  void printMem() const {
      std::cout << "  EventPads pt  " << &lPads_ << ",";
    }

    void genMCRings();

    inline std::vector<double> getPadVar() const { return padVar_; };
    inline void setPadVar( std::vector<double> padVar ) { padVar_ = padVar; };
 

  protected:


    CsRCEventPads();

    ~CsRCEventPads();

  
  private:

    static CsRCEventPads* instance_;

    std::list<CsRCPad*> lPads_;

    bool flag_;

    std::vector<double> padVar_;

  };

#endif
