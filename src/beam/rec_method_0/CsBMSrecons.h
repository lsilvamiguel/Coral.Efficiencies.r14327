// $Id: CsBMSrecons.h,v 1.70 2010/06/26 01:23:24 ybedfer Exp $

/*!
   \file CsBMSrecons.h
   \brief Compass BMS reconstruction Class.
   \author  G. Khaustov
   \version $Revision: 1.70 $
   \date    $Date: 2010/06/26 01:23:24 $
*/
#ifndef CsBMSrecons_h
#define CsBMSrecons_h
//#include "CsBMSconstants.h"

#include "CsHistograms.h"
#include "CsDetector.h"
#include "CsBeam.h"

#include <string>
#include <vector>

//
#define NTRACKMX 50                              // Max. number of BMS tracks
#define BMSBUFSIZE 100                            // Input buffer size
#define NBMSHODS    6                            // Number of Hodoscope planes
#define BMS_IDEAL_MC
//
typedef struct {
  int nplanes;
  int indx[NBMSHODS];
  int kfl[NBMSHODS];
  double tmean;
  double tchiq;
  double chiq;
  double p;
  double z[NBMSHODS];
} BMScomb;

//
class CsBMSrecons {
public:
        std::ifstream input_file;      // for TEST only !!!!!!!!!!!
         CsBMSrecons();
	 void test(void);
         void rawevprint(std::ostream &out) const;
         void prntclb(std::ostream &out, double calib[][NBMSHODS]) const;
         void bmsprnt(std::ostream &out) const;
         void getBMStrack(const int n, double &pmom, double &pchiq, 
                          double &mntime, double &tmchiq, 
                          int &nhit, int& nhttot,int &spec ) const;
	 void getBMStrack_resc(const int n, double &pmom, double &pchiq, 
                          double &mntime, double &tmchiq, 
                          int &nhit, int& nhttot, double& BmResol_resc ) const;
         inline int getBMSnmb(void) {return ntrack;};
         inline double getBMStm(const int n) {
                      if((n<0)||(n>=ntrack)) return(0);
                      return (TRKtime[n]); };
         void getBMSmoments(const int n, double& m0, double& m1) const;
         double getNewChi2(const int n, const double nmtime) const;
	 void bmsrec(void);
	 int getBMSMult(int plane);
	 double getBMSmom(const int n) const;
	 int getBMSThetaY(const int n, double &thetaY) const;
	 int getBMSY(const int n, double Y[]) const;
	 int getBMS_HitInfo(const int n, double Time[], int Channel[]) const;
	 int SFTrigger(double RefTime,int Mode,int trackID);
	 void TriggerChecks(bool Trigger[],int mode);
	 void CommonHits(std::vector<int> &vec);
	 void testRescueMethod(int TrackId, double slopeTg);
	 int  track_resc(void);
	 void moment_resc(int , double , double *, double *);
	 void MultBMStracks(std::vector<int> trackID);
	 //get raw hits info (time & y pos) for all BMS planes
	 bool getBMSRawHits(int plane, std::vector<double>& time, std::vector<double>& Y);
	 //get hit ids for track n
	 bool getBMStrackHits(int n, std::vector<int>& hits);
	 //get number of hodoscopes used
	 int getBMSnhod(void) {return nhod;}
         // Set BMS bits in argument CsBeam's hit patterns.
         bool setHitPattern(CsBeam* beam, int iBMStrack, bool rescue);
private:
//
//      constants and calibrations
//
        int Printlev;
        int itrg16;
        int mxnmbhit;                           //  mx.number of hits/plane
	int histoLevel;
	int DebugMode;
	int clusterIt;
	int allowedCommon;
	int year;
	double Cl_ChannelLimit;
	double Cl_TimeLimit;
	double BMS05TmDiff_min,BMS05TmDiff_max;
        double TMt0[64][NBMSHODS];
	bool   Fired[NBMSHODS];
	bool   FiredHit[NBMSHODS];
        double TDCSLP;                          // TDC time slope  
        double BMSMN1[10]; 
        double BMSpar[10];                      // Selection parameters
        double tmresol[NBMSHODS];                      // BMS time resolution
        double disp[NBMSHODS];
	double useTRAFFIC;
	std::string coeff_path[9];
	double p0;
	double coeff_avsg16[3][3][4];
	double coeff_evec16[3][3][4][4];
	double coeff_coff16[3][3][4][3][5];
	double method;
	bool method3;
	bool useBM05;				//Flag that indicates usage of BM05 plane
	bool useBM06;				//Flag that indicates usage of BM06 plane (implies useBM05)
	double SpecTrTmdiff_cntr;
	double SpecTrTmdiff_max;
	double SpecTrTmdiff_min;
	double SpecTrTmdiff;
	double EffTimeCut_min;
	double EffTimeCut_max;
	double BmResol_1_2,BmResol_1_3,BmResol_1_4,BmResol_2_3,BmResol_2_4;
//         
//       BMS description          
//
        int nhod;                               // Number of the BMS hodoscopes 
	int nhod_old;
        double BMStable[128][NBMSHODS];		// BMS hodoscope decoding table (only BM06 has 128 channels)
        CsDetector*  Id[NBMSHODS];                  // IDs for BMS Hods 0-3
        std::string hodnames[NBMSHODS];                  // names for BMS Hods 0-3
        int   hodsize[NBMSHODS];                    // Number of elements in the BMS hodoscopes
        int  Idpln0, Idpln1, Idpln2, Idpln3;    // IDs for BMS Hods 0-3
        std::string nmhod0, nmhod1, nmhod2, nmhod3;  // names for BMS Hods 0-3
//
//      beam track parameters
//
        int ntrackmx;
        int ntrack;
	int ntrack_resc;
        int nhittotal; 
        double TRKp[NTRACKMX]; 
        double TMchiq[NTRACKMX];
        double SPchiq[NTRACKMX];
        int    TRKhit[NTRACKMX];
        double TRKtime[NTRACKMX];
        double TRKy[NTRACKMX][NBMSHODS];
	int    HitUsedCl[BMSBUFSIZE][NBMSHODS];
	double TRKp_resc[NTRACKMX];
	double TMchiq_resc[NTRACKMX];
	double SPchiq_resc[NTRACKMX];
        int    TRKhit_resc[NTRACKMX];
        double TRKtime_resc[NTRACKMX];
        double TRKy_resc[NTRACKMX][NBMSHODS];
	double TRK_resol[NTRACKMX];
 
//
//       common working arrays
//
        int nfired;                            // Number of BMS planes fired 
        int khit16[NBMSHODS];                         // Number hits in hodoskope
        int ibmz16[BMSBUFSIZE][NBMSHODS];                   // hit cell numbers 
        int iflc16[BMSBUFSIZE][NBMSHODS];
        double tbmz16[BMSBUFSIZE][NBMSHODS];
        double zhit16[BMSBUFSIZE][NBMSHODS];                   // decoded hit coordinates
	double rawHtime[BMSBUFSIZE][NBMSHODS];
	int rawHchnl[BMSBUFSIZE][NBMSHODS];
	int    nRawH[NBMSHODS];
	int TRKindx[NTRACKMX][NBMSHODS]; 
	int TRKindx_resc[NTRACKMX][NBMSHODS]; 
        long int nrunold;
//
        void bookhist(void); 
	void bookhist1(void);
	void inithist(void);
	int  moment_resc(double *, double , double *, double *);
	void moment3(double *, double *, double *);
	void moment4(double *, double *, double *);
	int decode(void);
	int decodeCl(void);
	void tstdecode(void);
	void track(void);
	void SaveThisTrack(int nplanes, int ktrack, int indx[],int kfl[], double tmean,double tchiq,double z[]);
	void SaveThisTrack(int nplanes, int ktrack, int indx[],int kfl[],double tmean,double tchiq,double chiq,double p,double z[]);
	void ReplaceThisTrack(int ToBeReplaced, int Replacing);
	void StoreThisTrack(BMScomb &thisComb,int TrIndex);
	void BMStbl(void);
        bool readcalib(void);
        void init(void);
	void Efficiency(int i, double refTime);
	void TimeRes(void);
	void checkClustering(void);
	void checkBMStbl(void);
	void readBMSconst(std::string filename, int plane1, int plane2);
	void prntBMSconst(void);

//
        CsHist1D *HiNhit[NBMSHODS], *HiFired[NBMSHODS], *HiProf[NBMSHODS];
        CsHist1D *HiHtime[NBMSHODS],*HiHtimeTCut[NBMSHODS],*HiHtime1[NBMSHODS],*HiHtime2[NBMSHODS];
        CsHist1D *HiHpro4[NBMSHODS],*HiHpro3[NBMSHODS], *HiHpro1[NBMSHODS];
        CsHist1D *HiYpro4[NBMSHODS],*HiNHitafter[NBMSHODS];
	CsHist1D *Delta1H1, *Delta1H2, *Delta1H3, *Delta1H4, *DeltaG, *Delta, *DeltaT, *DeltaTclean,*DeltaTcleanBMSnew, *DeltaCut;
        CsHist1D *HiNTrack,  *HiBmom4,    *HiBmom3;
        CsHist1D *HiTMchi4,  *HiTMchi3,   *HiSPchi4, *HiSPchi3;      
        CsHist1D *HiTime4,   *HiTime3,    *HiNhtev;
        CsHist1D *HiBmom1,   *HiTMchi1,   *HiSPchi1, *HiTime1, *BMSPlanesfired;
        CsHist1D *Hit0t1dif, *Hit0t2dif,  *Hit0t3dif, *NuTrpEv,*NuTrpEvClean, *BMSPlanesfiredNT,*BMSPlanesfiredNTcl,*BMSPlanesfiredTr, *BMSdtTRX ;
	CsHist1D *BMSeff0,*BMSeff1,*BMSeff2,*NeventsCl,*TimeRes1_0102,*TimeRes2_0102,*TimeRes1_0304,*TimeRes2_0304,*TimeRes1_0104,*TimeRes2_0104;
        CsHist2D *Hit0t1,    *Hit0t2,     *Hit0t3, *EvSiVStri, *RescueDeltaMom, *RescueMom2H;
        CsHist2D *HiTmVsHt[NBMSHODS],  *HiN4h3h, *NHitP2vsP1, *NHitP2vsP3, *NHitP4vsP3, *P1Clus,*P2Clus,*P3Clus, *P4Clus,  *ThetaYVsMom;
	CsHist2D *HitBMS01VsHitBMS02,*HitBMS03VsHitBMS04,*HitBMS01VsHitBMS04,*CommonHitsHist3,*HitBMS02VsHitBMS01y,*HitBMS04VsHitBMS03y,*HitBMS03VsHitBMS02y,*HitBMS02VsHitBMS01y3,*HitBMS04VsHitBMS03y3,*HitBMS03VsHitBMS02y3;
	CsHist1D *Cl_chDiff,*Cl_tmDiff,*Cl_clSize,*Cl_HiProf[NBMSHODS],*Cl_HiHtime[NBMSHODS],*Cl_HiNhit[NBMSHODS]; 
	CsHist1D *TrackDec,*TchiqDiff,*BadDec4,*TchiqDiff3,*BadDec3,*RecPerfVStrig,*CommonHitsHist1,*CommonHitsHist2;
	std::map<int, CsHist1D*>   mH1;  //! 1D histogram pointers for level 1
	std::map<int, int, CsHist1D*>   mH2;
	std::map<int, CsHist1D*>   trigH;  //! 1D histogram pointers for trigger studies
	std::map<int, CsHist2D*>   trig2H;  //! 2D histogram pointers for trigger studies
	std::map<int, CsHist1D*>   Resc1D;  //! 1D histogram pointers for test of rescue algo
	std::map<int, CsHist1D*>   BeamH1D;
	std::map<int, CsHist2D*>   BeamH2D;
	CsHist1D *HiPlanesCheck[NBMSHODS];
	std::vector<CsHist1D*> HiCleanFired;
};


//
#endif // CsBMSrecons





