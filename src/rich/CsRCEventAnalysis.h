
#ifndef  EVENTANALYSIS_H
#define  EVENTANALYSIS_H

/*!
   \file    CsRCEventAnalysis.h
   \--------------------------
   \brief   CsRCEventAnalysis class declaration.
   \author  Paolo Schiavon
   \version 1.0
   \date    6 December 1999
*/


  class CsMCHit;
  class CsHist;

//-----------------------------
  class CsRCEventParticles;
  class CsRCEventClusters;
  class CsRCPartPhotons;
  class CsRCPhoton;
  class CsRCEventRings;
  class CsRCRing;
  class CsRCCircleFit;
//-----------------------------


  class CsRCEventAnalysis {


  public:


    CsRCEventAnalysis();

    CsRCEventAnalysis( const CsRCEventAnalysis& );

    void doEveAnalysis();

    void partPhoSelection();
    void ringSelection();

    void MCMonitor();

    void dataMonitor();

    void cluStructure( CsRCRing* );

    void PIDMonitor();

    void setProbsToTrack();

    void chkHitLists( CsRCRing*, int&, int&, std::list<CsMCHit*>&, 
		      std::list<CsMCHit*>& );
    void chkHitPrint( CsRCRing*, int&, int&, std::list<CsMCHit*>&,
		      std::list<CsMCHit*>& );

    void hitDisplay( int, std::list<CsMCHit*> );
    void BMDisplay( CsRCRing* );

    void doMirrDetFit();

    //^double getCorrTheta( const double );
    //^float offTheRing( const double );
    //^float offTheRing3( const double );

    void moniCFRefInf( CsRCRing* );
    bool fitCFRefInf( double& ,double& );
    bool fitGCFRefInd( double& ,double& );
    void moniCFRefIndR();

    void checkRingEff();
    void checkLikeDistr();
    void checkThetaLikeMax();

    void checkMCPMTs();
    void checkDataPMTs();
    void checkOptCorr();
    void singlePhotonCAT();
    void checkLikeDisTheta();

    std::string getVarStr( double, int, int );

    inline int nMCProcPart() const { return nMCProcPart_; };
    inline int nMCRecPart() const { return nMCRecPart_; };
    inline int nDaProcPart() const { return nDaProcPart_; };
    inline int nDaRecPart() const { return nDaRecPart_; };

    inline  bool flag() const { return flag_; };
    inline  bool flagMCMon() const { return flagMCMon_; };
    inline  bool flagDataMon() const { return flagDataMon_; };
    inline  bool flagPID() const { return flagPID_; };
    inline  bool flagAliM() const { return flagAliM_; };

    void print();

    ~CsRCEventAnalysis();

    static CsRCEventAnalysis* Instance();
  

  private:


    static CsRCEventAnalysis* instance_;

    int nMCProcPart_;
    int nDaProcPart_;
    int nMCRecPart_;
    int nDaRecPart_;

    bool flag_;
    bool flagMCMon_;
    bool flagDataMon_;
    bool flagPID_;
    bool flagAliM_;

  };

#endif
