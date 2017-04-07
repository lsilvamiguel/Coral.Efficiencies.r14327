// $Id: CsBeamRecons.h,v 1.59 2010/11/22 00:19:51 ybedfer Exp $

/*!
   \file CsBeamRecons.h
   \brief Compass Beam  reconstruction Class.
   \author  G. Khaustov
   \version $Revision: 1.59 $
   \date    $Date: 2010/11/22 00:19:51 $
*/
#ifndef CsBeamRecons_h
#define CsBeamRecons_h
#include<iostream>
#include "CsBMSrecons.h"
#include "CsBmFiHod.h"
#include "CsHistograms.h"
#include "CsBeam.h"

typedef struct {
  CsTrack* track;
  bool     Trig;
  double   time;
  bool     inComb;
} SFtrack;

typedef struct {
  int    SFtrackNr;
  int    BMStrackNr;
  double dt;
  double time;
  bool   Trig;
  bool   rescue;
  int    spec; //special flag: 0: nothing special, 1: BM05 used
  double Ex_bms_Y[6]; //backpropagated values of y at BMS postions
  double Ex_bms_Y_err[4][4]; //backpropagated values of y at BMS postions
  double Ex_bms_DY[2]; //backpropagated values of ThetaY at BMS postions
  double BackProp_Chi2; //Chi2 of backpropagated track (vs measured one)
  double BackProp_LH; //P(Chi2) of backpropagated track
  double UseIt; //flag to resolve ambiguous events (0: ignore it, 1: use it)
  double BmResol;
} BeamComb;

typedef struct {
    int BMSi;
    int SFi;
    unsigned int BMCombI;
    double Chi2;
} CombCand;

//This class is used to perform back propagation from SciFi/Si to BMS
class BackPropagated {
    public:
	BackPropagated();
	bool Extrapolate(double BMS_p, double BmResol, const CsHelix &Hin);
	void GetX(double (&X)[6]);
	double GetX(int plane) {return _X[plane];}
	void GetXerr(double (&X_err)[4][3]);
	double GetXerr(int plane) {return _X_err[plane][2];}
	void GetDX(double (&DX)[2]);
	double GetDX(int plane) {return _DX[plane/2];}
	void GetY(double (&Y)[6]);
	double GetY(int plane) {return _Y[plane];}
	void GetYerr(double (&Y_err)[4][4]);
	double GetYerr(int plane) {return _Y_err[plane][3];}
	void GetDY(double (&DY)[2]);
	double GetDY(int plane) {return _DY[plane/2];}

    private:
	std::vector<double> _X;
	std::vector<std::vector<double> > _X_err;
	std::vector<double> _DX;
	std::vector<double> _Y;
	std::vector<std::vector<double> > _Y_err;
	std::vector<double> _DY;
};

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
class CsBeamRecons {
public:
  static CsBeamRecons* Instance();
//
         CsBeamRecons();
         CsBMSrecons *BMS;
         CsBmFiHod BmHod;
//
         void beamprint(void);
	 void bmrecons(void);
	 void bmreconsTR(std::list<CsTrack*> tracks);
         void reconsMcExect(void);
	 bool IsGoodEvent(void){
	   //if(goodRandomTrigger) cout<<"is good event!"<<endl;
	   //else cout<<"is not a good event!"<<endl;
	   return goodRandomTrigger;
	 }; // gives true if eg in random trigger all SF have fired
         inline std::list<CsBeam*> getBeam() { return(bmtracks); }
         ~CsBeamRecons();
//
private:
//
//     coral.options file parameters
//
        int BMSrec;                    // flag: 0 - BMS reconstruction is OFF, 1 - is ON  
	int histoLevel;
	double useTRAFFIC;             //flag: 0 - the beam package is used as standalone, 1- Traffic is used for tracks 
        double BmChrg;                 // beam charge  
        double BmMoment;               // beam moment
        double BmResol;                // beam momentum resolution: dP = BmResol*P
        double MxTmdiff;               // BMS/FiHod max. track time difference
	double MxTmdiff_min;
	double MxTmdiff_max;
	double MxTmdiff_cntr;
	double MaxSFHitTime;
	double RecMethod;
	double DebugMode;
	double ThetaCorrCut4,ThetaCorrCut3;
	double thetaXcut,thetaYcut;
        int Printlev;                  // print level
        bool MCExect,MC;
        long int nrunold;
	CsDetector*  SFpnt[27];
	CsDetector*  UpTpnt[14];
	std::string detname[27];
	std::string detnameUpT[14];
	bool SFfired[10];
	bool goodRandomTrigger;
	int  SFMult[4];
	bool IT,MT,LT,OT,IMT,VT,VTO,VTI,BT,CT;
	bool Trigger[6];
	int  Filter_Flag;
	bool doRescue;
	bool useMultiples;
	float F1bin;
	double TotMeanCut;
	double TotMeanCut_min,TotMeanCut_max;
	double goodBMStrack_min,goodBMStrack_max;
	int year;
	int doRandom;
	int killBoS;
	int randomSFnumber;
	int randomSFnumber_SAS;
	int max_randomSFplanes;
	int verb;
	double GoodTrackMomCut;
	double RandomTargetCut;
	double EffTimeCut15;
	double EffTimeCutJ;
	double EffTimeCutJ_min,EffTimeCutJ_max;
	double TimeCut_VT_cntr[4];
	double TimeCutVt[4];
	double RandomTrackTimeCut;
	CsZone *BMS_zone;

	// ***** BACKTRACKING BMS <- BeamTelescope
#define bTNPLANES 4      //!< # of BMS planes considered in backTracking
	double bTResBMS[bTNPLANES]; //!< Resolution of considered BMS planes
	double bTCorr[bTNPLANES];   //!< Offset corrections
	bool   bTDefault;           //!< True if default bT coeffs requested

//
        std::list<CsBeam*>   bmtracks;      // beam track pointer list - main results of this staff
	std::list<CsDetector*> BMSzoneDets;
	std::list<CsHelix> goodSFhelixes; //X,Y position of bactracked tracks to BMS (this is only for testing BMS spatial acceptance)
//
	double *cov;
//
        void clear(void);
        void bookhist(void);
	void bookhist1(void);
	int  FiredSF();
	int SFEfficiency(int mode, double RefTime);
	int SFTrigger(void);
	int CheckVeto(int mode, double RefTime);
        void tmcorrel(int nbmhod, int nbms, int & ntrack);
        static CsBeamRecons* instance_;       // The Singleton Static pointer
        void FillCovMx(const double Dp);
	bool checktracktiming(CsTrack* track,int &nCl,int ntracks[]);
	bool checktrackorientation(CsTrack* track);
	bool findTriggerMask(int year);
	bool doRandomTrAna(double &GoodTrackTime,double &GoodTrackMom, CsHelix &GoodTrackStartHelix);
	void makeWiremap();
	void FillChi2(double BMSY[],double Ex_BMSY[],double &Chi2,double &LH);
	void HandleAmbig(int nBMSpassed, int nSFpassed,bool goodTrigger);
	void BmsGeometryCheck();//procedure used to check BMS geometrical acceptance
	std::vector<SFtrack>  SFtracks;
	std::vector<BeamComb> BeamCombV;
//
//       histograms
//
        CsHist1D *HiNtrack, *HiBmom,   *HiXprof,   *HiYprof;
        CsHist1D *HiXincl,  *HiYincl;
        CsHist1D *HiTime,   *HiTMchiq, *HiSPchiq;
        CsHist2D *HiNBNhod, *HiTBThod, *HiTBTh1;
        CsHist1D *HiTmdiff, *HiTmdif1, *HiNBMStrack, *HiNScFbtrack;
	std::map<int, CsHist1D*>   mH1;  //! 1D histogram pointers
	std::map<int, CsHist1D*>   mD1;  //! 1D histogram pointers
	std::map<int, CsHist2D*>   mH2;  //! 2D histogram pointers
	std::map<int, CsHist2D*>   mD2;  //! 2D histogram pointers
	std::map<int, CsHist1D*>   mHisto1;  //! 1D histogram pointers
	std::map<int, CsHist2D*>   mHisto2;  //! 2D histogram pointers
	std::map<int, CsHist1D*>   Htrig1;
	std::map<int, CsHist2D*>   Htrig2;
	std::map<int, CsHist1D*>   Ratrig1;
	std::map<int, CsHist2D*>   Ratrig2;
	std::vector<CsHist2D*> HiExBmsXY;
	std::vector<CsHist2D*> HiExBmsCorelX;
	std::vector<CsHist2D*> HiExBmsCorelY;
	CsHist1D *HiExBmsFired;
	CsHist1D *HiExBmsOldFired;
	CsHist1D *HiExBmsStuff;
	CsHist1D *HiExBmsEff;
	CsHist1D *HiBeamHitTimes;//Hit time distribution of hits that form a beam candidate
};
#endif // CsBeamRecons
