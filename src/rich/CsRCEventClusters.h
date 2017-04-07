
#ifndef  EVENTCLUSTERS_H
#define  EVENTCLUSTERS_H

/*!
   \file    CsRCEventClusters.h
   \---------------------------
   \brief   CsRCEventClusters class declaration.
   \author  Paolo Schiavon
   \version 0.01
   \date    October 2000, rev. August 2005
*/

// All found clusters of the RICH photon chambers.
// -----------------------------------------------


//--------------------------
  class CsRCEventPads;
  class CsRCCluster;
//--------------------------


  class CsRCEventClusters {


  public:


    static CsRCEventClusters* Instance();

    CsRCEventClusters( const CsRCEventClusters& );

    void clearEventClusters();
    void getEventClusters();

    void getClusFromPads();
    void getClusFromPads( int );

    void doClustering();
    void doClustering( int );
    void doClusteringAPV( int );

    void killHaloClus();
    void checkClusters();
    void checkPMTClus();

    inline  const std::list<CsRCCluster*> &lClusters() const {
      return lClusters_;};

    inline  bool flag() const { return flag_; };

    void print() const;
    void printFull() const;
    void printFull( const int& ) const;

    inline  void printMem() const {
      std::cout << "  EventClusters pt  " << &lClusters_ << ",";
    }


  protected:


    CsRCEventClusters();

    ~CsRCEventClusters();

  
  private:


    static CsRCEventClusters* instance_;

    std::list<CsRCCluster*> lClusters_;

    bool flag_;


  };

#endif
