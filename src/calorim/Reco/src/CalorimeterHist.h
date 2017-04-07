/*
   $Source: /local/reps/coral/coral/src/calorim/Reco/src/CalorimeterHist.h,v $
   $Date: 2010/06/18 10:44:19 $
   $Revision: 1.21 $
   -------------------------------------------------------------------------

   Set of test histograms for calorimeter's reconstruction.

   Authors:
     Vladimir  Kolosov   ( Vladimir.Kolosov@cern.ch,Kolosov@mx.ihep.su )
     Alexander Zvyagin   ( Alexander.Zviagine@cern.ch, Zvyagin@mx.ihep.su )
     Denis     Murashev  ( Denis.Mourachev@cern.ch, Murashev@sirius.ihep.su )

   Copyright(C): 1999-2001  V.Kolosov, A.Zvyagin, D.Murashev

     This library is free software; you can redistribute it and/or
     modify it under the terms of the GNU Library General Public
     License as published by the Free Software Foundation; either
     version 2 of the License, or (at your option) any later version.

     This library is distributed in the hope that it will be useful,
     but WITHOUT ANY WARRANTY; without even the implied warranty of
     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
     Library General Public License for more details.

     You should have received a copy of the GNU Library General Public
     License along with this library; if not, write to the Free
     Software Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
*/

#ifndef CalorimeterHist___include
#define CalorimeterHist___include
#include <vector>

#include "Reco_config.h"

class TH1D;
class TH2D;
class TProfile;
class TDirectory;
class TCanvas;

////////////////////////////////////////////////////////////////////////////////

namespace Reco {

class Calorimeter;

/*! \brief Set of histograms for reconstruction code testing.

*/
class CalorimeterHist
{
  //============================================================================
  // Constructors and destructor
  //============================================================================

  public:

    /// Destructor
    virtual            ~CalorimeterHist         (void);

    /// Constructor initialize all pointers to NULL
                        CalorimeterHist         (Calorimeter *calo);
  private:

    /// You can not use copy constructor.
                        CalorimeterHist         (const CalorimeterHist &);

  //============================================================================
  // Operators
  //============================================================================

  private:

    /// You can not use assignment operator.
    CalorimeterHist&    operator=               (const CalorimeterHist &);

  public:
    /// .Reset histograms
    virtual  void        Reset(void);

  //============================================================================
  // Attributes
  //============================================================================

  public:

    bool  raw_histos_booked;
    bool  led_histos_booked;
    bool  ped_histos_booked;
    bool  test_histos_booked;

    /// ROOT directory corresponed to the calorimeter
    TDirectory                         *root_dir;
    TDirectory                         *stat_dir;
    TDirectory                         *test_dir;
    TDirectory                         *cells_dir;
    TDirectory                         *cells_time_dir;
    TDirectory                         *cells_raw_dir;
    TDirectory                         *cells_led_dir;
    TDirectory                         *cells_ped_dir;
    TDirectory                         *noise_dir;
    TDirectory                         *cells_calib_dir;
    TDirectory                         *cells_XY_dir;

    TH1D                               *h1_ngen;
    TH1D                               *h1_ngen_Out;
    TH1D                               *h1_nrec;
    TH2D                               *h2_ngen_nrec;

    TH1D                               *h1_Xgen_Out;
    TH1D                               *h1_Ygen_Out;
    TH1D                               *h1_Egen_Out;
    TH2D                               *h2_XgYg_Out;

    TH1D                               *h1_Xgen;
    TH1D                               *h1_Ygen;
    TH1D                               *h1_Rgen;
    TH1D                               *h1_Egen;
    TH2D                               *h2_XgYg;

    TH1D                               *h1_Xgen_Norec;
    TH1D                               *h1_Ygen_Norec;
    TH1D                               *h1_Egen_Norec;
    TH2D                               *h2_XgYg_Norec;

    TH1D                               *h1_Xrec;
    TH1D                               *h1_Yrec;
    TH1D                               *h1_Rrec;
    TH1D                               *h1_Erec;
    TH1D                               *h1_Sumrec;
    TH1D                               *h1_SumrecCut;
    TH2D                               *h2_XrYr;

    TH1D                               *h1_Xdelta;
    TH1D                               *h1_Ydelta;
    TH2D                               *h2_EgEr;
    TH2D                               *h2_TgTr;
    TH2D                               *h2_XgXr;
    TH2D                               *h2_YgYr;
    TH2D                               *h2_XYdelta;
    TH1D                               *h1_Edelta;
    TH1D                               *h1_EdeltaN;

    TH1D                               *h1_XdeltaHadron;
    TH1D                               *h1_YdeltaHadron;
    TH1D                               *h1_XdeltaHadron40;
    TH1D                               *h1_YdeltaHadron40;
    TH2D                               *h2_EgErHadron;
    TH1D                               *h1_XdeltaMuon;
    TH1D                               *h1_YdeltaMuon;
    TH2D                               *h2_EgErMuon;

    TH2D                               *h2_EdeltaEg;
    TH2D                               *h2_EdeltaN_Eg;
    TH2D                               *h2_XdeltaXg;
    TH2D                               *h2_YdeltaYg;
    TH2D                               *h2_XdeltaYdelta;

    // For Cells Test
    TH2D                               *h2_XgYgCellTest;
    TH2D                               *h2_XrYrCellTest;
    TH2D                               *h2_XrYrMissedCell;
    TH2D                               *h2_XgXrCellTest;
    TH2D                               *h2_YgYrCellTest;
    TProfile                           *p1_DxXgCell;
    TProfile                           *p1_DyYgCell;

    // For Near Gamma Test
    TProfile                           *p1_Hi2Ndf;
    TProfile                           *p1_Hi2Prob;
    TH2D                               *h2_Hi2E;

    TH2D                               *h2_ProbE;
    TH2D                               *h2_ProbE_4;
    TH2D                               *h2_ProbE_9m4;
    TH2D                               *h2_ProbE_25m9;

    TH1D                               *h1_Prob;
    TH1D                               *h1_Prob_4;
    TH1D                               *h1_Prob_9m4;
    TH1D                               *h1_Prob_25m9;

    TH1D                               *h1_Pool_4;
    TH1D                               *h1_Pool_9m4;
    TH1D                               *h1_Pool_25m9;

    TH2D                               *h2_DENorm_4;
    TH2D                               *h2_DENorm_9m4;
    TH2D                               *h2_DENorm_25m9;

    TProfile                           *p1_Hi2Cells;
    TProfile                           *p1_DENorm;
    TProfile                           *p1_DEabsNorm;
    TProfile                           *p1_DEcalc;
    TProfile                           *p1_HiMaxOK;
    TH2D                               *h2_HiHiMaxOK;
    TProfile                           *p1_HiMaxBad;
    TH2D                               *h2_HiHiMaxBad;
    TH2D                               *h2_Dist2;
    TH2D                               *h2_Dist2emin;
    TH2D                               *h2_Dist2Cell;

    // For total hited cells and total energy at a time
    TH1D                               *h1_ETotal;
    TH1D                               *h1_ESumADC;
    TH2D                               *h2_NHitedEtotal;
    TH2D                               *h2_XhitYhit;
    TH2D                               *h2_XhitYhitE;
    // Collect general statistic on calibrations
    TH1D                               *h1_CalibFactorBase;
    TH1D                               *h1_CalibFactorLED;
    TProfile                           *p1_CalibFactorEdep;
    TProfile                           *p1_CalibFactorTis;

    std::vector <TH1D *>               h1_ErCell;
    std::vector <TH1D *>               h1_ADCmCell;

    TH1D                               *h1_Time;
    TH2D                               *h2_TimeE;

    std::vector <TH1D *>               h1_TimeCell;

    std::vector <TH1D *>               h1_RAWCell;
    std::vector <TH1D *>               h1_LEDCell;
    std::vector <TH1D *>               h1_PEDCell;

    std::vector <TH1D *>               h1_NOISECell;
    std::vector <TH1D *>               h1_NOISEGamCell;

    TH2D                               *h2_XYAll;
    std::vector <TH2D *>               h2_XdeltaYdeltaCell;
    std::vector <TH1D *>               h1_Time_subset;

    // Some more calibration histo
    std::vector <TH1D *>               h1_CalibCell;
    std::vector <TH2D *>               h2_CalibInSpill;

    // Some general statistic histograms
    TH1D                               *h1_calib_new;
    TH1D                               *h1_calib_old;
    TH1D                               *h1_calib_prim;
    TH1D                               *h1_calib_new_over_calib_prim;
    TH1D                               *h1_calib_new_over_calib_old;
    TH1D                               *h1_calib_old_over_calib_prim;
    TH1D                               *h1_calib_monitor_over_calib_old;
    TH2D                               *h2_XYout_calib_monitor_over_old;
    TH1D                               *h1_led_monitor_over_led_old;
    TH1D                               *h1_led_new;
    TH1D                               *h1_led_old;
    TH1D                               *h1_ped_new;
    TH1D                               *h1_ped_old;
    TH1D                               *h1_time_new;
    TH1D                               *h1_time_old;
    TH1D                               *h1_ecut_low_old;
    TH1D                               *h1_ecut_high_old;
    TH1D                               *h1_ecut_low_new;
    TH1D                               *h1_ecut_high_new;
    TH1D                               *h1_ecut_sparse_mode;
    TH1D                               *h1_ecut_min;

    std::vector <TH1D *>               h1_calib_old_subset;
    std::vector <TH1D *>               h1_calib_monitor_over_old_subset;
    std::vector <TH2D *>               h2_XYout_calib_monitor_over_old_subset;

    TCanvas                            *c;

    protected:
    Calorimeter                        *calorimeter_;
};

////////////////////////////////////////////////////////////////////////////////

class InternalHisto
{

  public:
    InternalHisto(Calorimeter *calo) : calorimeter_(calo) {}

    virtual ~InternalHisto(void) {}

    TProfile *p1_Hi2Ndf;
    TProfile *p1_Hi2Prob;
    TProfile *p1_Hi2Cells;
    TProfile *p1_DENorm;
    TProfile *p1_DEabsNorm;
    TProfile *p1_DEcalc;
    TProfile *p1_HiMaxOK;
    TH2D *h2_HiHiMaxOK;
    TProfile *p1_HiMaxBad;
    TH2D *h2_HiHiMaxBad;
    TProfile *p1_ClusterSizeOfInterest;
    TProfile *p1_ClusterSize_absEcut;
    TProfile *p1_ClusterSize_relEcut;

    protected:
    Calorimeter                        *calorimeter_;

};

////////////////////////////////////////////////////////////////////////////////

class ProfilesHisto
{

  public:
    const static unsigned BIG_ANGLE_BINS_SIZE = 10;
    const static unsigned SMALL_ANGLE_BINS_SIZE = 10;
    const static unsigned ENERGY_BINS_SIZE = 10;

    ProfilesHisto( void ) {
                             double e_range[ENERGY_BINS_SIZE+1] =      { 0., 0.5, 1., 2., 4., 10., 20., 50., 100., 200., 1000.};
                             double b_angle_range[BIG_ANGLE_BINS_SIZE+1] =   { -0.25,-0.20,-0.15,-0.10,-0.05,0.00,0.05,0.10,0.15,0.20,0.25 };
                             double  s_angle_range[SMALL_ANGLE_BINS_SIZE+1] = { -0.10,-0.08,-0.06,-0.04,-0.02,0.00,0.02,0.04,0.06,0.08,0.10};
                             for ( unsigned i=0; i< ENERGY_BINS_SIZE+1; i++) energy_range[i]=e_range[i];
                             for ( unsigned i=0; i< BIG_ANGLE_BINS_SIZE+1; i++) big_angle_range[i]=b_angle_range[i];
                             for ( unsigned i=0; i< SMALL_ANGLE_BINS_SIZE+1; i++) small_angle_range[i]=s_angle_range[i];
                           }

    virtual ~ProfilesHisto(void) {}

    double  energy_range[ENERGY_BINS_SIZE+1];
    double  big_angle_range[BIG_ANGLE_BINS_SIZE+1];
    double  small_angle_range[SMALL_ANGLE_BINS_SIZE+1];

    int IndexEnergy ( double e ) const {  for ( unsigned i=0; i< ENERGY_BINS_SIZE; i++)
                                            {
                                               if( e >= energy_range[i] && e < energy_range[i+1] ) return i;
                                            }
                                            return ENERGY_BINS_SIZE-1;
                                         }
    int IndexBigAngle ( double a ) const {  for ( unsigned i=0; i< BIG_ANGLE_BINS_SIZE; i++)
                                             {
                                               if( a >= big_angle_range[i] && a < big_angle_range[i+1] ) return i;
                                             }
                                             return BIG_ANGLE_BINS_SIZE-1;
                                           }
    int IndexSmallAngle ( double a ) const {  for ( unsigned i=0; i< SMALL_ANGLE_BINS_SIZE; i++)
                                              {
                                                if( a >= small_angle_range[i] && a < small_angle_range[i+1] ) return i;
                                              }
                                              return SMALL_ANGLE_BINS_SIZE-1;
                                            }

    TProfile *ProfileX_Electrons;
    TProfile *ProfileY_Electrons;

    TProfile *ProfileX_Muons;
    TProfile *ProfileY_Muons;

    TProfile *ProfileX_Pions;
    TProfile *ProfileY_Pions;

    TProfile              *ProfileX_BigAngle_Energy_Electrons[BIG_ANGLE_BINS_SIZE][ENERGY_BINS_SIZE];
    TProfile              *ProfileX_SmallAngle_Energy_Electrons[SMALL_ANGLE_BINS_SIZE][ENERGY_BINS_SIZE];
    TProfile              *ProfileY_BigAngle_Energy_Electrons[BIG_ANGLE_BINS_SIZE][ENERGY_BINS_SIZE];
    TProfile              *ProfileY_SmallAngle_Energy_Electrons[SMALL_ANGLE_BINS_SIZE][ENERGY_BINS_SIZE];

    TProfile              *ProfileX_BigAngle_Energy_Muons[BIG_ANGLE_BINS_SIZE][ENERGY_BINS_SIZE];
    TProfile              *ProfileX_SmallAngle_Energy_Muons[SMALL_ANGLE_BINS_SIZE][ENERGY_BINS_SIZE];
    TProfile              *ProfileY_BigAngle_Energy_Muons[BIG_ANGLE_BINS_SIZE][ENERGY_BINS_SIZE];
    TProfile              *ProfileY_SmallAngle_Energy_Muons[SMALL_ANGLE_BINS_SIZE][ENERGY_BINS_SIZE];

    TProfile              *ProfileX_BigAngle_Energy_Pions[BIG_ANGLE_BINS_SIZE][ENERGY_BINS_SIZE];
    TProfile              *ProfileX_SmallAngle_Energy_Pions[SMALL_ANGLE_BINS_SIZE][ENERGY_BINS_SIZE];
    TProfile              *ProfileY_BigAngle_Energy_Pions[BIG_ANGLE_BINS_SIZE][ENERGY_BINS_SIZE];
    TProfile              *ProfileY_SmallAngle_Energy_Pions[SMALL_ANGLE_BINS_SIZE][ENERGY_BINS_SIZE];
};

////////////////////////////////////////////////////////////////////////////////

class CorrelationsHisto
{

  public:
    CorrelationsHisto(Calorimeter *calo) : calorimeter_(calo) {}

    virtual ~CorrelationsHisto(void) {}

    TH2D                               *h2_XYAll;
    std::vector <TH2D *>               h2_XdeltaYdeltaCell;
    std::vector <TH2D *>               h2_XYCellCorrelations;

    protected:

    Calorimeter                        *calorimeter_;

};

////////////////////////////////////////////////////////////////////////////////

class TrigGroupHist
{
  //============================================================================
  // Constructors and destructor
  //============================================================================

  public:

    /// Constructor initialize all pointers to NULL
                        TrigGroupHist         (void);
  private:

    /// You can not use copy constructor.
                        TrigGroupHist         (const TrigGroupHist &);

  //============================================================================
  // Operators
  //============================================================================

  private:

    /// You can not use assignment operator.
    TrigGroupHist&    operator=               (const TrigGroupHist &);

  public:
    /// .Reset hstograms
    void        Reset(void) {};

  //============================================================================
  // Attributes
  //============================================================================

  public:

    bool  test_histos_booked;

    /// ROOT directory corresponed to the calorimeter
    TDirectory                         *root_dir;

    std::vector <TH1D *>                h1_energy_;
    std::vector <TH1D *>                h1_max_energy_;
    std::vector <TH2D *>                h2_XYe_;
    std::vector <TH2D *>                h2_Ee_;
    std::vector <TH1D *>                h1_E_;
    std::vector <TProfile *>            p1_FindMe_;

    TCanvas                            *c;
};

//class MakeParticlesHisto :  public HistoCollection
class MakeParticlesHisto
{
  public:
  /// Base constructor
    MakeParticlesHisto   (const std::string &histo_path ) : root_dir(NULL), booked_ok_(false) {}

    //============================================================================
    // Attributes, data
    //============================================================================
  public:
    TDirectory                         *root_dir;
    bool                booked_ok_;
    TH1D              * caloE;            //  Energy of Calorimeter clusters
    TH1D              * caloSumE;         //  Summ Energy of Calorimeter clusters
    TH1D              * deltaX;           //  XCalorimeter - XTrack distribution for all clusters
    TH1D              * deltaXelectron;
    TH2D              * deltaXelectronX;
    TProfile             * MeanDeltaXelectronX;
    TH2D              * deltaXelectronX_selected;
    TProfile             * MeanDeltaXelectronX_selected;

    TH1D              * deltaXpion;
    TH2D              * deltaXpionX;
    TH1D              * deltaXmuon;
    TH2D              * deltaXmuonX;
    TH2D              * deltaXelectronE;
    TH2D              * deltaXpionE;
    TH2D              * deltaXmuonE;
    TH1D              * deltaXwEas;       //  XCalorimeter - XTrack distribution for completely associated clusters
    TH1D              * deltaXS;          //  (XCalorimeter - XTrack)/SigmaXCalo distribution for all clusters
    TH1D              * deltaY;           //  YCalorimeter - YTrack distribution for all clusters
    TH1D              * deltaYelectron;
    TH2D              * deltaYelectronY;
    TProfile             * MeanDeltaYelectronY;
    TH2D              * deltaYelectronY_selected;
    TProfile             * MeanDeltaYelectronY_selected;

    TH1D              * deltaYpion;
    TH2D              * deltaYpionY;
    TH1D              * deltaYmuon;
    TH2D              * deltaYmuonY;
    TH2D              * deltaYelectronE;
    TH2D              * deltaYpionE;
    TH2D              * deltaYmuonE;
    TH1D              * deltaYwEas;       //  YCalorimeter - YTrack distribution for completely associated clusters
    TH1D              * deltaYS;          //  (YCalorimeter - YTrack)/SigmaYCalo distribution for all clusters
    TH2D              * caloXtrackX;      //  XCalorimeter vs XTrack distribution for all clusters
    TH2D              * caloYtrackY;      //  YCalorimeter vs YTrack distribution for all clusters
    TH2D              * deltaXdeltaY;      //
    TH2D              * deltaXdeltaYall;      //

    TH1D              * deltaXelectron_NearBadCell;
    TH1D              * deltaYelectron_NearBadCell;
    TH1D              * deltaEsElectron_NearBadCell;

    TH1D              * deltaXelectron_InGoodZone;
    TH1D              * deltaYelectron_InGoodZone;
    TH1D              * deltaEsElectron_InGoodZone;

    TH2D              * caloXtrackXcell;      //  XCalorimeter vs XTrack in cell distribution for all clusters
    TH2D              * caloYtrackYcell;      //  YCalorimeter vs YTrack in cell distribution for all clusters
    TH2D              * caloXmuonXcell;      //  XCalorimeter vs XTrack in cell distribution for all clusters
    TH2D              * caloYmuonYcell;      //  YCalorimeter vs YTrack in cell distribution for all clusters
    TH2D              * caloXpionXcell;      //  XCalorimeter vs XTrack in cell distribution for all clusters
    TH2D              * caloYpionYcell;      //  YCalorimeter vs YTrack in cell distribution for all clusters
    TH2D              * caloXelectronXcell;      //  XCalorimeter vs XTrack in cell distribution for all clusters
    TH2D              * caloYelectronYcell;      //  YCalorimeter vs YTrack in cell distribution for all clusters

    TH2D              * deltaXtrackXcell;      //  XCalorimeter - XTrack vs XTrack in cell distribution for all clusters
    TH2D              * deltaYtrackYcell;      //  YCalorimeter - YTrack vs YTrack in cell distribution for all clusters
    TH2D              * deltaXmuonXcell;      //  XCalorimeter - XTrack vs XTrack in cell distribution for all clusters
    TH2D              * deltaYmuonYcell;      //  YCalorimeter - YTrack vs YTrack in cell distribution for all clusters
    TH2D              * deltaXpionXcell;      //  XCalorimeter - XTrack vs XTrack in cell distribution for all clusters
    TH2D              * deltaYpionYcell;      //  YCalorimeter - YTrack vs YTrack in cell distribution for all clusters
    TH2D              * deltaXelectronXcell;      //  XCalorimeter - XTrack vs XTrack in cell distribution for all clusters
    TH2D              * deltaYelectronYcell;      //  YCalorimeter - YTrack vs YTrack in cell distribution for all clusters
    TH2D              * deltaXelectronXcell_selected;      //  XCalorimeter - XTrack vs XTrack in cell distribution for all clusters
    TH2D              * deltaYelectronYcell_selected;      //  YCalorimeter - YTrack vs YTrack in cell distribution for all clusters

    TH2D              * caloXY;           //  YCalorimeter vs XCalorimeter distribution for all clusters
    TH2D              * trackXY;          //  YTrack vs XTrack distribution for all tracks which hit calorimeter
    TH2D              * deltaXtrackX;     //  XTrack vs XCalorimeter - XTrack distribution for all clusters
    TH2D              * deltaYtrackY;     //  YTrack vs YCalorimeter - YTrack distribution for all clusters
    TH2D              * deltaSqXtrackSqX; //  covTrackX  vs (XCalorimeter - XTrack)^2
    TH2D              * deltaSqYtrackSqY; //  covTrackY  vs (YCalorimeter - YTrack)^2
    TH2D              * deltaSqXelectronSqX; //  covTrackX  vs (XCalorimeter - XTrack)^2 for electrons
    TH2D              * deltaSqYelectronSqY; //  covTrackY  vs (YCalorimeter - YTrack)^2 for electrons
    TH2D              * deltaSqXpionSqX; //  covTrackX  vs (XCalorimeter - XTrack)^2 for pions
    TH2D              * deltaSqYpionSqY; //  covTrackY  vs (YCalorimeter - YTrack)^2 for pions
    TH2D              * deltaSqXmuonSqX; //  covTrackX  vs (XCalorimeter - XTrack)^2 for muons
    TH2D              * deltaSqYmuonSqY; //  covTrackY  vs (YCalorimeter - YTrack)^2 for muons
    TH2D              * deltaXtrackXwEas; //  XTrack vs XCalorimeter - XTrack distribution for associated(including E) clusters
    TH2D              * deltaYtrackYwEas; //  YTrack vs YCalorimeter - YTrack distribution for associated(including E) clusters
    TH2D              * trackXYas;        //  YTrack vs XTrack distribution for all tracks associated with calorimeter cluster
    TH2D              * trackXYnotas;     //  YTrack vs XTrack distribution for all tracks which hit but not associated with calorimeter cluster

    TH2D              * electronXYas;        //  YTrack vs XTrack distribution electron associated with calorimeter cluster
    TH2D              * electronXYnotas;     //  YTrack vs XTrack distribution electron which hit but not associated with calorimeter cluster
    TH2D              * pionXYas;         //  YTrack vs XTrack distribution pion associated with calorimeter cluster
    TH2D              * pionXYnotas;      //  YTrack vs XTrack distribution pion which hit but not associated with calorimeter cluster
    TH2D              * muonXYas;         //  YTrack vs XTrack distribution muon associated with calorimeter cluster
    TH2D              * muonXYnotas;      //  YTrack vs XTrack distribution muon which hit but not associated with calorimeter cluster


    TH2D              * trackXYas_E;      //  YTrack vs XTrack distribution for all tracks with momentum associated with calorimeter cluster
    TH2D              * trackXYnotas_E;   //  YTrack vs XTrack distribution for all tracks with momentum not associated with calorimeter cluster
    TH2D              * trackXYas_noE;    //  YTrack vs XTrack distribution for all tracks without momentum associated with calorimeter cluster
    TH2D              * trackXYnotas_noE; //  YTrack vs XTrack distribution for all tracks without momentum not associated with calorimeter cluster
    TH2D              * caloEtrackE;      //  ETrack vs ECalorimeter for all tracks associated with calorimeter cluster
    TH2D              * deltaEtrackE;     //  ECalorimeter-ETrack vs ETracks for all tracks associated with calorimeter cluster
    TH1D              * deltaEs;          //  Ecalo/Etrack
    TH1D              * deltaEselectron;          //  Ecalo/Etrack
    TProfile             * MeanDeltaEsElectronX;
    TProfile             * MeanDeltaEsElectronY;

    TH1D              * deltaEspion;          //  Ecalo/Etrack
    TH1D              * deltaEsmuon;          //  Ecalo/Etrack
    TH2D              * deltaEsEtrack;    //  ETrack vs Ecalo/Etrack
    TH2D              * deltaEsEelectron; //  Eelectron vs Ecalo/Eelectron
    TH2D              * deltaEsEpion;     //  Epion vs Ecalo/Epion
    TH2D              * deltaEsEmuon;     //  Emuon vs Ecalo/Emuon

    TH2D              * caloTtrackT;      //  Time Track vs Time Calorimeter for all tracks associated with calorimeter cluster
    TH1D              * neutralE;         //  Eneutral
    TH2D              * neutralXY;        //  Y vs X distribution of neutral particles produced from not associated calorimeter clusters
    TH1D              * gammaE;            //  Egamma
    TH2D              * gammaXY;          //  Y vs X distribution of gamma particles produced from not associated calorimeter clusters
    TH1D              * electronE;        //  Egamma
    TH2D              * electronXY;       //  Y vs X distribution of gamma particles produced from not associated calorimeter clusters
    TH1D              * pionE;            //  Egamma
    TH2D              * pionXY;           //  Y vs X distribution of gamma particles produced from not associated calorimeter clusters
    TH1D              * muonE;            //  Egamma
    TH2D              * muonXY;           //  Y vs X distribution of gamma particles produced from not associated calorimeter clusters
    TH2D              * caloEelectronE;    //  E electron vs ECalorimeter
    TH2D              * caloEpionE;       //  E pion vs ECalorimeter
    TH2D              * caloEmuonE;       //  E muon vs ECalorimeter
    TH2D              * distGamPiE;         //  E deposit in calorimeter vs pion-gamma distance
    TH2D              * muonEmE4low;
    TH2D              * muonEmE4high;
    TH2D              * pionEmE4low;
    TH2D              * pionEmE4high;
    TH2D              * electronEmE4low;
    TH2D              * electronEmE4high;
    TH2D              * gammaEmE4low;
    TH2D              * gammaEmE4high;
    TH2D              * distGamGamE;
    TH2D              * distXYGamGam;
    TH1D              * MGamGam;
    TH1D              * MGamGam_NearBadCell;

    TH1D              * Ngamma;            // Number gamma particles associated with calorimeter clusters
    TH1D              * NKl;               // Number Kl particles associated with calorimeter clusters
    TH1D              * NPiplus;           // Number Piplus particles associated with calorimeter clusters
    TH1D              * NPiminus;          // Number Piminus particles associated with calorimeter clusters
    TH1D              * NMuplus;           // Number Muplus particles associated with calorimeter clusters
    TH1D              * NMuminus;          // Number Muminus particles associated with calorimeter clusters
    TH1D              * NElectrons;        // Number electrons associated with calorimeter clusters
    TH1D              * NPositrons;        // Number electrons associated with calorimeter clusters

    TH2D              * distXYElectronGam;
    TH2D              * distXYPiGam;

    TProfile             * effCaloAss;    //  Probability of Calorimeter cluster association with track (What is the sense? - Just fraction of clusters from tracks)

    TProfile             * effTrackAss;        //  Probability of Track Calorimeter association
    TProfile             * effPiplusAss;       //  Probability of Piplus Calorimeter association
    TProfile             * effPiminusAss;      //  Probability of Piminus Calorimeter association
    TProfile             * effMuplusAss;       //  Probability of Muplus Calorimeter association
    TProfile             * effMuminusAss;      //  Probability of Muminus Calorimeter association
    TProfile             * effElectronAss;    //  Probability of electrons Calorimeter association
    TProfile             * effPositronAss;    //  Probability of positrons Calorimeter association

    TH1D              * ZlastTrackAss;        // Zlast of Track in Calorimeter associated
    TH1D              * ZlastPiplusAss;       // Zlast of Piplus in Calorimeter associated
    TH1D              * ZlastPiminusAss;      // Zlast of Piminus in Calorimeter associated
    TH1D              * ZlastMuplusAss;       // Zlast of Muplus in Calorimeter associated
    TH1D              * ZlastMuminusAss;      // Zlast of Muminus in Calorimeter associated
    TH1D              * ZlastElectronAss;     // Zlast of electrons in Calorimeter associated
    TH1D              * ZlastPositronAss;     // Zlast of positrons in Calorimeter associated

    TH1D              * ZlastTrackNotAss;        // Zlast of Track in Calorimeter not associated
    TH1D              * ZlastPiplusNotAss;       // Zlast of Piplus in Calorimeter not associated
    TH1D              * ZlastPiminusNotAss;      // Zlast of Piminus in Calorimeter not associated
    TH1D              * ZlastMuplusNotAss;       // Zlast of Muplus in Calorimeter not associated
    TH1D              * ZlastMuminusNotAss;      // Zlast of Muminus in Calorimeter not associated
    TH1D              * ZlastElectronNotAss;     // Zlast of electrons in Calorimeter not associated
    TH1D              * ZlastPositronNotAss;     // Zlast of positrons in Calorimeter not associated

    // For Near Gamma Test
    TH2D                               *electronHi2OK;
    TH2D                               *electronHi2Bad;
    TProfile                           *electronERNorm;
    TProfile                           *electronETNorm;
    TProfile                           *electronDENorm;
    TProfile                           *electronDEabsNorm;
    TProfile                           *electronDEcalc;
    TH2D                               *electronHiMaxOK;
    TH2D                               *electronHiHiMaxOK;
    TH2D                               *electronHiMaxBad;
    TH2D                               *electronHiHiMaxBad;

    TH2D                               *h2_Dist2;
    TH2D                               *h2_Dist2emin;
    TH2D                               *h2_Dist2Cell;

    TH1D                               *gammaHi2OK;
    TH1D                               *gammaHi2Bad;
    TH2D                               *gammaHi2EOK;
    TH2D                               *gammaHi2EBad;

    TH1D                               *gammaHi2OK_1gEvt;
    TH1D                               *gammaHi2Bad_1gEvt;
    TH2D                               *gammaHi2EOK_1gEvt;
    TH2D                               *gammaHi2EBad_1gEvt;

    TH1D                               *gammaHi2OK_2gEvt;
    TH1D                               *gammaHi2Bad_2gEvt;
    TH2D                               *gammaHi2EOK_2gEvt;
    TH2D                               *gammaHi2EBad_2gEvt;


    TProfile                           *gammaERNorm;
    TProfile                           *gammaETNorm;
    TProfile                           *gammaDENorm;
    TProfile                           *gammaDEabsNorm;
    TProfile                           *gammaDEcalc;

    TH2D                               *gammaHiMaxOK;
    TH2D                               *gammaHiHiMaxOK;
    TH2D                               *gammaHiMaxBad;
    TH2D                               *gammaHiHiMaxBad;

    TH2D                               *gammaHiSignOK;
    TH2D                               *gammaHiHiSignOK;
    TH2D                               *gammaHiSignBad;
    TH2D                               *gammaHiHiSignBad;

    TH2D                               *gammaHiMaxOK_1gEvt;
    TH2D                               *gammaHiHiMaxOK_1gEvt;
    TH2D                               *gammaHiMaxBad_1gEvt;
    TH2D                               *gammaHiHiMaxBad_1gEvt;

    TH2D                               *gammaHiSignOK_1gEvt;
    TH2D                               *gammaHiHiSignOK_1gEvt;
    TH2D                               *gammaHiSignBad_1gEvt;
    TH2D                               *gammaHiHiSignBad_1gEvt;

    TH2D                               *gammaHiMaxOK_2gEvt;
    TH2D                               *gammaHiHiMaxOK_2gEvt;
    TH2D                               *gammaHiMaxBad_2gEvt;
    TH2D                               *gammaHiHiMaxBad_2gEvt;

    TH2D                               *gammaHiSignOK_2gEvt;
    TH2D                               *gammaHiHiSignOK_2gEvt;
    TH2D                               *gammaHiSignBad_2gEvt;
    TH2D                               *gammaHiHiSignBad_2gEvt;
};

class FitInfoHisto
{
  public:
  /// Base constructor
//    FitInfoHisto   (const std::string &histo_path ) : root_dir(NULL), booked_ok_(false) {}
    FitInfoHisto   ( void ) : root_dir(NULL), booked_ok_(false) {}

    //============================================================================
    // Attributes, data
    //============================================================================
  public:
    TDirectory        *root_dir;
    bool               booked_ok_;
    TH1D              * MeanCALIB;            //  Mean of calibration
    TH1D              * SigmaCALIB;           //  Sigma of calibration
    TH1D              * MeanNCALIB;           //  Mean Norm of calibration
    TH1D              * SigmaNCALIB;          //  Sigma Norm of calibration
    TH1D              * StatCALIB;            //  Stat of calibration
};

} // namespace Reco

#endif // CalorimeterHist___include
