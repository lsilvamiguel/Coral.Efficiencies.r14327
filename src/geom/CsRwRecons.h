// $Id: CsRwRecons.h,v 1.4 2008/10/14 12:28:34 alekseev Exp $

/*!
   \file CsBeamRecons.h
   \brief Compass Beam  reconstruction Class.
   \author  G. Khaustov
   \version $Revision: 1.4 $
   \date    $Date: 2008/10/14 12:28:34 $
*/
#ifndef CsRwRecons_h
#define CsRwRecons_h
#include<iostream>
#include "CsHistograms.h"
#include "CsECAL1.h"

typedef struct  ec1rec{
       int  zone;
       double  X;
       double  Y;
       double  Zx;
       double  Zy;
       double  E;
       double nhits;
       double  Mom;
       CsHelix hel;
       int    nfirst;
       int    npart;
       int    calobj;
} ec1rec;
typedef struct  ec1out{
       int    npart;
       Reco::CalorimeterParticle* calobj;
} ec1out;
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
class CsRwRecons {
public:
  static CsRwRecons* Instance(void);
  CsRwRecons();
  ~CsRwRecons();
  void RwRecons(void);
  std::vector<ec1out*>  dataout;
//
private:
       static CsRwRecons* instance_;       // The Singleton Static pointer
//
       void corcoeff(double Eg,double wtot, int nz, double& a,double& b);
       void PRrecons_consts(void);
       void coordw(CsCluster* mycl, double &weight);
       void CrCorrection(int npl, int lb, double Xg, double Zg, double w4[],double xmean4[],
                        double z0cl[], double & Xg4,double & Zg4x, double &wxmn);
       void  gm_gm_distance(std::vector<ec1rec> & ec1,  int ng, 
                           double & Dmin, double& Dx, double& Dy, double & x, double& y,
                           double& zx, double & x1, double & y1,double& zy1);
       void  kill_ghnost(std::list<ec1rec>&  ec1, bool& yes, int mod);
       void  kill_all_ghnosts(std::list<ec1rec>&  ec1, bool& yes1);
       void  sum_digits(int ibg, double X,double Y, double Z,double deltaE,
                       double deltaCR,double &wtot,double wpl[],double xmean[]);
       void  sum_CRdigits( int ibg,double delta, double Xg, double Yg, 
                          double xmean[], double w4[], double xmean4[]);
       void  ec1energy_cor(double xcut, double ycut,  std::list<CsDigit*>& ec1dig, double X, double Y,
                           const  std::vector<Reco::Cell> & cells, double &Enew);
       void find_best_MCvertex(double x, double y, double zx,  double zy, 
                              double xtarg,double ytarg, double ztarg, double &alfa,double & beta);
       double  GetLike( int tag, CsTrack *track );
       int LikePid(CsTrack *track);
       void Clear(void);
       bool gmchk(double Xg, double Yg, double delta);
       int wrchk(int prj, CsCluster* hit);
       bool htchk(int ipl, double Xg, double Yg, int nz, double delta);
       void m2g(std::vector<ec1rec>&  ec1,int nh1,int nh2, int nh3, int nh4);
//
//     coral.options file parameters
//
        bool   RW_rec, RW_MC;
        int RW_histoLevel;
        int EnergyCorrection,CoordCorrection,PhotonSelection;
        bool readcal;
        int Ncorpoints;
        double Epoints[20];
        double Ecorr_points[3][20], hits_slp[3][20];
        double hits_slope[3],Ecorr[3];
        double deltaEn[3], deltaCRd[3];
        double deltaCR1;
        double deltaPR_overlap[3];
        double deltaSpace_overlap[3];
//
        static const int nplanes=8;
        CsDetector*  RWdet[nplanes];
        int plmap[nplanes];
        int goodmap[nplanes];
        double z0cl[nplanes];
        int nwires[nplanes];
        CsECAL1 *calorimeter;
        CsDet *idet1;
//
//       histograms
//
        CsHist1D  *hi1D[200];
        CsHist2D  *hi2D[100];
};
#endif // CsBeamRecons
