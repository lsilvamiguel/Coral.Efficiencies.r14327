// $Id: CsBmFiHod.h,v 1.8 2005/07/29 12:35:04 conrad Exp $

/*!
   \file CsBmFiHod.h
   \brief Compass Beam SciFb reconstruction Class.
   \author  G. Khaustov
   \version $Revision: 1.8 $
   \date    $Date: 2005/07/29 12:35:04 $
*/
#ifndef CsBmFiHod_h
#define CsBmFiHod_h
#define NCOMBMX 1024                               // Working space
#define BMHODTRKMX 200                             // Max. number of BmHod tracks
#define INPBUFSIZE 100                              // Max. inp buff size
#define NBMHODS    4                            // Number of Hodoscope planes
//
#include<iostream>
#include<string>
#include "CsHistograms.h"
#include "CsDetector.h"
#include "CsStartOfRun.h"
#include "CsCluster.h"

class CsBmFiHod : public CsStartOfRun  {
public:
         CsBmFiHod();
	 void bhodrecons(void);
         void rawevprint(std::ostream &out) const;
	 void bmhodprnt(std::ostream &out) const;
	 void test(void);
         inline int getBmHodnmb(void) const {return ntrack;};
         void getBmFiHodtrack(const int n, double &x, 
              double & y, double& z, double& dxdz, double& dydz, 
              int& nhits, int& tothits, bool& trig, double &mntime, double &tmchiq) const;
         inline double getBmFiHodtm(const int n) const {
                if((n<0)||(n>=ntrack)) return(0);
                return (HDtime[n]);}
         void getBmHodmoments(const int n, double& m0, double& m1) const;
         double getNewChi2(const int n, const double mntime) const;
         void getBmErrors(double& dx, double& dy, double& ddxdz, double& ddydz) const;
         inline std::list<CsZone*> getZones(void) const { return myzones;};
         std::list<CsCluster*> getClusters(const int n) const;
         bool sor(void);
         ~CsBmFiHod();
//
private:
//
//      constants and calibrations
//
        int  Printlev;
	double useTRAFFIC;             //flag: 0 - the beam package is used as standalone, 1- Traffic is used for tracks 
        double BHODt0[NBMHODS][96];               // Hodoscope time calibration
        double TDCslope;                          // TDC slope    
        double selparam[10];                      // selection constants for rec. procedure 
        int MxNmbHit;                             // Max number of hits/plane
        static const int ntrackmx;
//         
//       Beam hodoscope description          
//
        static const int nhod;                     // Number of the fiber hodoscope planes 
        static const std::string hodnames[NBMHODS];     // names for BMS Hods 0-3
        CsDetector*  Id[NBMHODS];                  // IDs for BMS Hods 0-3
        int   hodsize[NBMHODS];                     // Number of elements in the hodoscope
        double   tmresol[NBMHODS];                   //  hodoscope time resolution
        double disp[NBMHODS];
        std::list<CsZone*> myzones;
//
//      output beam track parameters 
//
        int ntrack;
        int ntothits;                         // tot. number of hits in 4 planes
        double HDtchiq[BMHODTRKMX];
        double HDtime[BMHODTRKMX];
        double HDhittm[2][BMHODTRKMX][NBMHODS];
        int HDhit[3][BMHODTRKMX][NBMHODS];
        int HDlabel[BMHODTRKMX];
        CsCluster* outclust[2][BMHODTRKMX][NBMHODS];
//
//       inputs/common working arrays
//
        int nfired;                            // Number of BMS planes fired 
        int khit[NBMHODS];                     // Number hits in hodoscope plane
        int ihit[INPBUFSIZE][NBMHODS];         // hit channels 
        double thit[INPBUFSIZE][NBMHODS];      // hittimes
        CsCluster* inpclust[INPBUFSIZE][NBMHODS];
//
        int ihitsv[INPBUFSIZE][NBMHODS];       // copy of working arrays 
        double thitsv[INPBUFSIZE][NBMHODS];     
//
        void bookhist(void);   
	int decode(void);
	void onehit(void);
	void multyhit(void);
        bool readcalib(void);
        bool init(void);
        bool selonehit(int ipl, double trigwin, bool& onehit, int& nh, double& t);
        void eff(double t1, double t2, int ipl3, int nh, bool goodhod);
        bool readT0 ( int nchanread, std::string filename, double* coeff); 
        void prntclb(std::ostream &out, double calib[][96]) const;
//
        CsHist1D *hist1[200];
        CsHist2D *hist2[200];
};
//
#endif // CsBmFiHod
