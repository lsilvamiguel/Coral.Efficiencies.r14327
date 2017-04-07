/*
   $Source: /afs/cern.ch/user/g/gobbo/w0/coral/src/calorim/Reco/src/Calorimeter.h,v $
   $Date: 2011/02/01 22:05:52 $
   $Revision: 1.125 $
   -------------------------------------------------------------------------

   This file is part of cellular calorimeter reconstruction program.

   Authors:
     Vladimir  Kolosov   ( Vladimir.Kolosov@cern.ch, Kolosov@mx.ihep.su )
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

#ifndef Calorimeter___include
#define Calorimeter___include
#define inside_Calorimeter_h

#include <cassert>
#include <ctime>
#include <list>
#include <map>
#include <string>
#include <utility>   // for std::pair
#include <vector>

#include "Reco_config.h"
#include "CalorimeterParticle.h"
#include "Cell.h"
#include "Exception.h"
#include "StatInfo.h"

class TH1D;
class TH2D;
class TDirectory;
class TLorentzVector;
class StoreLED;
class DevMonitorFEM ;

////////////////////////////////////////////////////////////////////////////////

namespace Reco {

class CellDataRaw;
class CellType;
class Cluster;
class DataBase;
class GUICalorimeter;
class MCConstruction;
class OneParticleResponse;
class Reconstruction;

class CalorimeterHist;
class CorrelationsHisto;
class FitInfoHisto;
class InternalHisto;
class MakeParticlesHisto;
class ProfilesHisto;
class TrigGroupHist;

/// Cellular calorimeter.
class Calorimeter
{
  //============================================================================
  // Types, constants
  //============================================================================

  public:

    enum                DataType                {DataReal,DataMC};

    enum                CellInfoType            {CALIB=0,LED=1,PED=2,TIME=3,NOISE=4,CLUSTER=5,RAW=6};
    static const unsigned CellInfoTypeSize=7;

    enum                CalibTimeType           {OLD=0,NEW=1,MONITOR=2,PRIM=3,MC=4};
    static const unsigned CalibTimeTypeSize=5;

    enum                CalorimeterType         {ElectroMagnetic,Hadronic};
    // Recommended meaning for misc(used outside Reco:: code) bool options
    enum                MiscBool                { USED_IN_MP=0,
						  USED_TO_MAKE_GAMMA=1,
						  SET_ELECTRON_ID=2,
						  EXCLUDE_CENTRAL_EP_REGION_FROM_GAMMA_SEARCH=3,
						  SELECT_IN_MP_BY_TIME=4,
                                                  STORE_POSITION_INFO=5 };
    // Recommended meaning for misc(used outside Reco:: code) double options
    enum                MiscDouble              {EMIN_MAKE_GAMMA=0, RELAX_GATE_FOR_ASS=1, TIME_MIN=2, TIME_MAX=3,
                                                  MONITOR_TIME_BIN=4 };
   // Recommended meaning for misc(used outside Reco:: code) time options
    enum                MiscTime                {MONITOR_TIME_MIN=0, MONITOR_TIME_MAX=1};

    static const char  *CellInfoNames[];

    enum                FitMethod             { Simple=0, Normal=1, NoFit=3};

#include "CalorimeterSubSet.h"
#include "CalorimeterFitInfo.h"
#include "CalorimeterOptions.h"
#include "CalorimeterSummaryInspectLED.h"
#include "CalorimeterCellsMatrix.h"
#include "CalorimeterAlignmentInTime.h"

// class CaloEvent
// {
//   public:
//      CaloEvent ( size_t ncells, const EventID &event_id ) : event_id_(event_id)  { for( size_t i=0; i<ncells; i++ ) data_.push_back(0.);}
//   public:
//     EventID  event_id_;
//     std::vector <double>   data_;
// };
//
// class StatInSpill
// {
//   public:
//      StatInSpill ( void ) : ev_cnt(0) {}
//   public:
//    int    ev_cnt;
//    size_t    evmin;
//    size_t    evmax;
//    time_t    tmin;
//    time_t    tmax;
// };
//
// class StoreCaloEvents
// {
//   public:
// //    static std::map < size_t,  std::vector<double> > time_start_run_map;
//     const static double spill_gate =   22.7;
//     const static double spill_tcycle = 46.792;
//   public:
//      StoreCaloEvents ( Calorimeter *c );
//      void InitStatic ( void );
//      void AddEvent ( CaloEvent  * ev);
//      void Process ( void );
//      void ProcessInSpill ( void );
//      bool Processed ( void ) const { return   ((vev_id_.size() >0)&& (led_no_cuts.size() > 0)); }
//      size_t  NCells ( void ) const { return ncells_;}
//      void SetSpillsInfo ( size_t run, size_t spill );
//      const std::pair<  int, int> & GetSpillsInfo( size_t run ) const;
//
//      const StatInfo &  GetValue( size_t run, size_t spill, size_t cell ) const;
//      const StatInfo &  GetStrictAverageValue ( size_t cell ) const;
//
//      static std::pair< bool,std::pair< double, double>  >  GetSpillTime ( size_t run, size_t spill );
//
//   public:
//     Calorimeter   * c_;
//     size_t  run_min_;
//     size_t  run_max_;
//     size_t  spill_min_;
//     size_t  spill_max_;
//     std::map <  std::pair< size_t, size_t >,  std::vector< CaloEvent  * > >  events_;
//     std::map <  size_t, std::pair<  int, int >  > spill_min_max_;
//     std::map <  std::pair< size_t, size_t >,  StatInSpill >  statinspill_;
//     std::vector < EventID > vev_id_;
//
//     std::vector<StatInfo> led_no_cuts;
//     std::vector<StatInfo> led_thr;
//     std::vector<StatInfo> led_strict;
//     std::map <  std::pair< size_t, size_t >,  std::vector<StatInfo> > leds_in_spill;
//     StatInfo zero;
//     std::pair<  int, int >  dummy_spill_info_;
//     size_t   ncells_;
// };

  // ============================================================================
  // Constructors, destructor
  // ============================================================================

  public:

    /// The destructor.
    virtual             ~Calorimeter             (void);

    /// Construct calorimeter with given name from the file with geometry description
                         Calorimeter             (const std::string &the_name,
						  const std::string &geom_file = "");


  // ==========================================
  // General Methods
  // ==========================================

  public:

     /*! \brief Calorimeter initialization after Calorimeter object creation.
                Should not be called from constructor.
                Some initializations depends on the options setting and it is not really clear
                during construction whether the options have final settings.
     */
    /// Initialization used outside Coral
    virtual void        InitializeWithOptions   ( const std::list<std::string>& reco_opt_list );
    virtual void        Initialize              (void);

    /// Print properties
    void Print (std::ostream &o=std::cout, const std::string &prefix="") const;

  protected:
    /// Calorimeter initialization.
    virtual void        Init                    (void);
    /// Allocate memory during Calorimeter initialization.
    void                InitMemory              (void);

    /// Allocate memory during Calorimeter creation.
  public:
    void                InitMemoryPublic        (void) { InitMemory(); }
//     /// Options initialization (virtual method), to be overridden in derived classes.
//     virtual void        InitOptions             (void);
//     virtual void        UpdateAfterOptionsSettings (void);

  public:

    /// \return Calorimeter name.
    const std::string       &GetName                 (void) const {return name;}

    /// \return Calorimeter type.
    CalorimeterType          GetType            (void) const {return options.calo_type;}

    ///  Print General Calorimeter Info
    void       PrintGeneralInfo                 (void) const;

    /// Clear Calorimeter cell amplitudes.
    virtual void        Clear                   (void);

    /// \return list of calorimeter cells.
    const std::vector<Cell> &GetCells                (void) const {return cells;}

    /// \return false if calibrated data out of reasonable size /energy .
    virtual bool CheckCalibData     (const std::vector<CellDataRaw>& data, double threshold,
				     double sum_energy_cut, int nmax_cells_cut) const;

    /// \return Calorimeter SubSets.
    const std::vector<SubSet>      &GetSubSets       (void) const {return calorimeter_sub_sets;}
    std::vector<SubSet>            &GetSubSetsToModify (void) {return calorimeter_sub_sets;}
    /// \return Calorimeter SubSet index(if any) by SubSets Name, \return -1 if there is no such SubSet.
    int                            GetSubSetIndex     (const std::string &name) const;

    /// Add Calorimeter SubSet
    void                            AddSubSet       ( SubSet s ) {calorimeter_sub_sets.push_back(s);}
    void                            AddSubSet       ( SubSet s, bool is_trigger_group )
                                                        {calorimeter_sub_sets.push_back(s);
                                                         if( is_trigger_group )
                                                           trigger_groups_index_.push_back(calorimeter_sub_sets.size()-1);}

    void                            ClearSubSet      ( void ) { calorimeter_sub_sets.clear(); }

  // ==========================================
  // Methods for Calorimeter construction anyway
  // ==========================================

    /// Insert new CellType (without checks)
   void                             InsertCellType  ( CellType cell_type);
    /// Insert Cell ( CellTypes Shoulde be there already)
   void                             InsertCell  ( Cell cell );

  // ==========================================
  // Methods for dealing with cell's information
  // ==========================================

  public:

    /// Clear OLD/NEW info for CALIB,LED,PED
    virtual void        CellsInfoClear          (int what=-1,int when=-1);

    /// Clear OLD/NEW info for CALIB,LED,PED
    virtual void        CellsInfoPrint          (int what=-1,int when=-1);

    /// Use NEW valuest instead of OLD one.
    virtual void        CellsInfoUpdate         (CellInfoType what);

    /// \return cell's information
    const StatInfo     &GetCellInfo             (CellInfoType what, CalibTimeType when,size_t cell) const;
    double              GetEnergyCutInCell        (size_t cell) const {return energy_cut_bad_cells_old[cell];} // VeryBad!
    double              GetEnergyGammaCutInCell        (size_t cell) const {return energy_gamma_cut_bad_cells_old[cell];} // VeryBad!
    double              GetMinEnergyCutInCell     (size_t cell) const {return energy_cut_[0][cell];} // VeryBad also!
    /// \return energy depending on signal and time in spill
    /// \param signal signal in ADC counts
    /// \param tis    Time in Spill, if set to -1 tis-dependent calib is ignored
    double              SignalToEnergy          (size_t cell, double signal, double tis) const;
    double              GetCorrTimeInSpillLED ( size_t cell, double tis ) const;
    void                 SetCorrTimeInSpillLED ( size_t cell, double m0, double m1, double m2, double stat );

    /// Set cell's information
    void                SetCellInfo             (size_t what, size_t when,size_t cell,const StatInfo &ci);

    /// Copy cell's information from one place to another
    void                CopyCellsInfo             (size_t what_source, size_t when_source,
                                                     size_t what_target, size_t when_target);

    /// \return cell's fit information
    const FitInfo      &GetFitInfo             (size_t what,size_t cell) const { return fit_info[what][cell]; }

    /// \return cell's fit information
    void                SaveFit2CellInfo       (int what=-1);

    /// Set LED amplitude depending on run and spill
    void                SetInSpillMeanLED ( int icell, int run, int ispill, double v);

    /// \return LED amplitude depending on run and spill
    double              GetInSpillMeanLED ( int icell, int run, int ispill ) const;

    /// \return LED amplitude depending on run and spill for current spill
    double              GetInCurrentSpillMeanLED (int icell) const;
    /// \return base calibration constant
    double              GetCalibrationFactorBase(size_t cell ) const;
    /// \return LED/LASER correction factor to calibration
    double              GetCalibrationFactorLED(size_t cell ) const;
    /// \return energy \param energy correction factor to calibration
    double              GetCalibrationFactorEdep(size_t cell, double energy ) const;
    /// \return Time in spill \param tis correction factor to calibration
    double              GetCalibrationFactorSdep(size_t cell, double tis ) const;

    double              GetECal0_nonlinear_calib_correction(size_t cell, double signal ) const;

  // ==========================================
  // Methods for Calorimeter calibration.
  // ==========================================

  public:

    /// Set calibration procedure and arguments for it.
    void                SetCalibProc            (bool flag,double energy) const {options.calibration_flag=flag; options.calibration_energy=energy;}

    /// This is some prototype for electron beam calibration.
    bool                CalibrateSummEnergy     (double esumm,
                                                 const  std::vector <OneParticleResponse> &hop);

    /// Perform time calibration if it is possible
    void                CalibrateTime           (void);
    /// Perform time calibration using reconstructed particles
    void                CalibrateTime  (const std::vector<CalorimeterParticle> &particles );

    /// Scale calibration by the options.scale_calibration. Information recovery function.
    void                ScaleCalibration        ( size_t when );

    /// Scale calibration by the cells corrections. Information recovery function.
    void                ScaleCalibrationByCellCorrections        ( size_t when );

    /// Scale calibration by the options.scale_calibration
    double              GetScaleCalibration     (void) const {return options.scale_calibration;}

    /// Scale calibration by using LED info, practical importance, if the options.scale_calibration_by_led is set true.
    void                ScaleCalibrationByLED    (void);

    /// Correct calorimeter position using options.correct_position. Information recovery function.
    void                CorrectPosition          (void);

    ///  \return correction constant to Calorimeter X position
    double              GetPositionCorrectionX   (void) const {return position_correction_[0];}

    ///  \return correction constant to Calorimeter Y position
    double              GetPositionCorrectionY   (void) const {return position_correction_[1];}

    ///  \return correction constant to Calorimeter Z position
    double              GetPositionCorrectionZ   (void) const {return position_correction_[2];}

    ///  \return success flag and vector with X,Y position of the  Calorimeter depending on time stamp
    std::pair<bool,std::vector< double > > GetPosition ( const time_t  &time ) const;


    /// This is a prototype for physic data calibration.
    void                CalibrateAtMesonMass    (size_t ipart,size_t jpart,double mass,double mass_sigma,double mass_tab,
                                                const std::vector<CalorimeterParticle> &particles);
    /// This is a prototype for physic data calibration.
    void                CalibrateParticle      (size_t ipart,double energy,
                                                const std::vector<CalorimeterParticle> &particles,
                                                 bool fill_default_calib_histo);
    /// Calibrate using calibration options.
    bool                CalibrateParticle      (size_t ipart,
                                                const std::vector<CalorimeterParticle> &particles);
    /// Calibrate particular cell
    bool                CalibrateCell          (size_t mcell,double factor,double weight,
                                                                    bool fill_default_calib_histo);

    /// Make inspection of calibration results and some other similar methods
    virtual int        MakeInspections         (void);
    int                InspectCalibration      (void) const;

    int                MakeInspectionsDB       (void) const;
    int                InspectCalibrationDB    (void) const;
    /// Display some
    int                InspectLEDDB            (void) const;
    int                InspectPorogiDB         (void) const;


    virtual void       MonitorCalibration      (void);
    StatInfoStore      *GetStatInfoStore        (size_t what, size_t when,size_t cell);
    virtual void       FillMonitorHisto        ( void );
    virtual void       InitMonitoring          ( int run_start, int run_finish );
    virtual void       StoreMonitorInfo        ( int run );
    virtual void       CalculateMonitorCalibration ( int run );
    virtual void       MonitorLEDInBurst       ( int spill );
     void              ResetMonitorLEDInBurstHisto ( void );

    /// Fill histogramms for visual calibration control
    /// Default calibration histo
    void                FillCalibrationHisto(size_t icell,double factor, double weight);
    void                FillCalibrationMesonHisto(size_t ipart,size_t jpart,double mass,double mas_tab,
                                     const std::vector<CalorimeterParticle> &particles);
    void                FillCalibrationHisto(size_t ipart,double energy,
                                     const std::vector<CalorimeterParticle> &particles);
    // Fill calibrations using options
    void                FillCalibrationHisto(size_t ipart,
                                     const std::vector<CalorimeterParticle> &particles);

  // ==========================================
  // Fast Monte-Carlo Methods to generate calorimeter response (mostly under development)
  // ==========================================

  public:

    /*! \brief  Get a flat distributed random number
     *
     * This method can be overwritten in, e.g., CORAL to use the same random
     * number generator as for all detector classes, so that no special care
     * about seeding has to be taken.
     * The random number generator used is not properly seeded.
    */
    virtual double      GetRandomFlat(void) const;

    /*! \brief  Get a gaussian distributed random number
     *
     * This method can be overwritten in, e.g., CORAL to use the same random
     * number generator as for all detector classes, so that no special care
     * about seeding has to be taken.
     * The random number generator used is not properly seeded.
    */
    virtual double      GetRandomGaus(void) const;

    /// Fast Monte-Carlo method.  Generate response in calorimeter from set of Monte-Carlo particles.
    virtual void        MonteConstruction      (const std::vector<CalorimeterParticle> &particles,
                                                std::vector<CellDataRaw>& data);

    /*! \brief Fast Monte-Carlo method.  Generate and store response in calorimeter from one Monte-Carlo particle.
               MakeMCResponce() function should be called afterwards once per event.
    */
    void                MakeMCHits(const CalorimeterParticle &particle_in,CalorimeterParticle &particle_out,
				   std::vector<CalorimeterParticle> &other_particles_out);

    /// Fast Monte-Carlo method.  Generate response in calorimeter from stored hits.
    void                MakeMCResponse(void);

    /// Fast Monte-Carlo method.  Generate response in calorimeter from stored hits.
    void                MakeMCDigitization(void);

    /// Fast Monte-Carlo method. Add constant term fluctuations to Monte-Carlo amplitudes in Calorimeter.
    void                MonteAddConstantFluctuations (void);

    /// Fast Monte-Carlo method. Add photo-electron statistic fluctuations to Monte-Carlo amplitudes in Calorimeter.
    void                MonteAddStochasticFluctuations (void);

    /// Fast Monte-Carlo method. Add readout fluctuations to Monte-Carlo amplitudes in Calorimeter.
    void                MonteAddReadoutFluctuations (void);

    /// Add smearing to current calibration.
    void                AddCalibrationSmearing ( double s );

  public:
    /// Add COMGEANT MC Hit
    virtual bool        AddMCHit                (int id_hit, const void *data);

  // ==========================================
  // Reconstruction Methods
  // ==========================================

  public:

    /*! \brief Reconstruction.
    */
    virtual bool        Reconstruction          (void);

    /*! \brief Process monitoring and self-research functions.void Calorimeter::
    */
  protected:
    void                CheckPrintManyParticlesResponse( void );
  public:

    virtual int         ProcessOnDuty               (void);
    virtual int         ProcessOnDutyLEDEvent       (void){ return 0;}
    virtual int         ProcessOnDutyPedestalsEvent (void){ return 0;}
    /*! \brief Sometimes we need to repeat Reconstruction.
	     if  option.repeat_reconstruction is FALSE RepeatReconstruction function is DUMMY
    */
    virtual bool        RepeatReconstruction    (void);

    /// Some recovery(if any) of raw data.
    virtual void        RawDataRecovery         (void){}

    /// Processing LED event. Function is used outside Coral to obtain LED calibrations.
    virtual void        ProcLED                 (void);
    const StoreLED *GetStoreLED ( void ) const { return store_led_;}
    void                  SetStoreLED (StoreLED * s ) {store_led_=s;}
     StoreLED       * GetStoreLED2Modify ( void ) { return store_led_;}
    const DevMonitorFEM *GetFEM ( void ) const { return fem_;}
    void                  SetFEM( DevMonitorFEM * fem ) {fem_=fem;}
    DevMonitorFEM   * GetFEM2Modify ( void ) { return fem_;}
//     virtual void        ProcInSpillLED                 (void);
//     virtual void        ProcessStoredLED           (void);

    /// Inspection LED info
    int                 InspectLED                (void);
//     void                SetLedEventsID ( const std::vector < EventID > &led_events_id ) {led_events_id_=led_events_id;}
//   std::vector < EventID > & GetLedEventsIDs (void) const { return led_events_id_; }

//   const std::vector < EventID > & GetLedEventsIDs (void) const { if( store_leds_all_== NULL) assert(false); return store_leds_all_->vev_id_;}
//    std::vector < EventID > v_dummy_;
//    const std::vector < EventID > & GetLedEventsIDs (void) const {   std::cout << " GetLedEventsIDs:: Dont use this function please !!! " << std::endl; return v_dummy_; }

    /// Special functions for PEDESTAL event. Implemented for GAMS ADC data.
    virtual bool        CheckPedestalIsOK       ( const CellDataRaw& data ) const;
    virtual void        RecoverBadPed           (void);
    /// Processing PEDESTAL event.
    virtual void        ProcPedestal            (void);
    int                 InspectPedestal         (void) const;

    /// Filling RAW/LED/PED histograms for each cell
    virtual void        FillHistoRAW            (void);
    virtual void        FillHistoLED            (void);
    virtual void        FillHistoPED            (void);

    bool CellIsBad(  size_t icell, size_t when) const;
    int GetBadCellStatus(  size_t icell, CalibTimeType when) const;
    void FilterOutBadCells(const std::vector<CellDataRaw>& data_in, std::vector<CellDataRaw>& data_out, size_t when) const;
    /// \return vector of Bad Cells around cell
    std::vector<size_t>  BadCellsAround          (size_t cell) const;
    const std::list<size_t>&  GetBadCells             (CalibTimeType when) const { return bad_cells_[when]; }

    const double* GetEnergyCutSparseMode() const { return energy_cut_sparse_mode; }
 public:
    /// Options initialization (virtual method), to be overridden in derived classes.
    virtual void        InitOptions             (void);
    virtual void        UpdateAfterOptionsSettings (void);
    virtual void        SetRecoOptions          (const std::string &opt) { options.SetRecoOptions(opt); }
    /// TODO Clarify usage. It seems function is only used outside coral.
    void                 SetOptions ( const std::list<std::string>& reco_opt_list);
    /// TODO ToBe removed ???
    virtual void        SetOptions              (const std::string &opt) { return SetRecoOptions(opt); }

    /// Get reconstruction options
    virtual const Options& GetOptions           (void) const {return options;}

    /// Get reconstruction options to modify
    Options&                GetOptions2Modify    (void) {return options;}


  public:

    /// Allow to manipulate options reading externally. Not used in Coral.
    bool                 OptionsSettingIsClosed  (void ) { return options.Ok();}
    void                 SetOptionsOK            (void ) { options.SetOK();}
    void                 AllowToSetMoreOptions   (void ) { options.More();}
    /*! \brief  NEW!!! GetCalorimeterParticles has replaced old GetRecoParticles.
                Really public function return the result of reconstruction in the
                form of CalorimeterParticle (something different from Calorimeter::Particle) objects
                Warning:: CalorimeterParticles NOW in Mother Reference System  !!!
    \arg particles Output - the result of reconstruction. CalorimeterParticles in MRS.
    */
    virtual std::vector<CalorimeterParticle>  GetCalorimeterParticles    (void) const;

    /*! \brief  Set the result of reconstruction in the Calorimeter externaly without running Reconstruction()
                getting the result from the CalorimeterParticle  objects
       WARNINGS: 1.CalorimeterParticles in Mother Reference System  !!!
                 2.As Reconstruction() was not running Calorimeter state may not be self-consistent!!
    \arg particles Intput - the result of reconstruction. CalorimeterParticles in MRS.
    */
    void        SetCalorimeterParticles    (const std::vector<CalorimeterParticle> &particles);

    /*! \brief   Recalculate particles position according to GetPositionCorrectionX(YZ) values
    \arg particles Intput - particles to be corrected
    */
    void           RecalculateParticlesPosition    ( std::vector<CalorimeterParticle> &particles) const;
    void           RecalculateActualToStartOfRunPosition    ( std::vector<CalorimeterParticle> &particles) const;
    /*! \brief   Recalculate particles energy and time according to new calibartion constants
    \arg particles Intput - particles to be corrected
    \arg particles Output - particles in the result of corrections.
    */
    std::vector <Reco::CalorimeterParticle>   RecalculateParticles( const std::vector<CalorimeterParticle> &particles) const;

    /*! \brief   Recalculate particles energy according to new coefficients
    \arg particles Intput - particles to be corrected
    */
    void           RecalculateParticlesEnergy      ( std::vector<CalorimeterParticle> &particles) const;

    /*! \brief  Set some hint particles.
        aimed to use in reconstuction algorithms, calorimeter alignment and reconstruction tests
    \arg particles Intput - may be useful in reconstruction
    */
    void                SetHintParticles        (const std::vector<CalorimeterParticle> &particles);

    /*! \brief   Add new hint particle
                WARNINGS: in contrrary to SetHintParticles here particle in Mother Reference System  !!!
    */
    void                AddHintParticle         (const CalorimeterParticle &particle);

    /*! \brief  GetHintParticles.
    \arg particles Output -  return back hint particles
    */
    const std::vector<CalorimeterParticle> &GetHintParticles        (void) const {return hint_particles;}

    /*! \brief  Convert point from MRS to DRS.
    \arg  Intput -  vector <double> &x_in
          Output -  vector<double> &x_out .
    */
     virtual void          ConvertMRS2DRS   (const std::vector <double> &x_in, std::vector<double> &x_out ) const;

    /*! \brief  Convert point from DRS to MRS.
    \arg  Intput -  vector <double> &x_in
          Output -  vector<double> &x_out .
    */
     virtual void          ConvertDRS2MRS   (const std::vector <double> &x_in, std::vector<double> &x_out ) const;

    /*! \brief  Convert Particle from DRS to MRS.
    \arg particle Intput - some particle to be converted
                  Output - converted particle .
    */
     virtual CalorimeterParticle      ConvertParticleMRS2DRS   (const CalorimeterParticle &particle) const;

    /*! \brief  Convert Particle from DRS to MRS.
    \arg particle Intput - some particle to be converted
                  Output - converted particle .
    */
     virtual CalorimeterParticle      ConvertParticleDRS2MRS   (const CalorimeterParticle &particle) const;

    /*! \brief  Calorimeter Test Reconstruction Functions
        This method aimed to compare generated set of particles (hint_particles)
                            with reconstructed set of particles (reco_particles),
                            plus      some other useful histograms.
    */

    virtual void           ReconstructionTest      (void);
    virtual void           TestHistoRawInfo        (void);
    virtual void           HistInternalInfo        (void);
    virtual void           TrigGroupTest           (void);
    void                   MakeProfiles            (void);
    void                   FillParticlesCaloAssociationHist(void);
    void                   FillFitInfoHist        (void);

  protected:

    /*! Test
    */
    void                TestEnergyDelivery   (void);

  public:

    double              EnergyDispersionInCell ( double energy_in_cell, size_t jcell ) const;
    double              EnergyDispersionInCellBasic ( double energy_in_cell, size_t jcell ) const;
    double              EnergyDispersionInCellReadout ( double energy_in_cell, size_t jcell ) const;
    // well, well, well develoopmeeent
    double              EnergyDispersionInCellDigitization ( double energy_in_cell, size_t jcell ) const;
    double              ScaleDispersionInSparceMode(double e, double dispe, size_t jcell) const;


    double              EnergySigmaInCell ( double energy_in_cell, size_t jcell ) const;
    double              EnergySigmaInCellBasic ( double energy_in_cell, size_t jcell ) const;
    double              EnergySigmaInCellReadout ( double energy_in_cell, size_t jcell ) const;

  protected:
    virtual void       InitPreshowerCorrections( void );
  public:
    virtual double     MakePreshowerCorrections( double eg, double xpart, double ypart, CalorimeterParticle::ParticleID pid) const;

  // ==========================================
  //  Data Manipulation Functions
  // ==========================================

 protected:

    ///  Store RAW data statistic
    void                StoreRawDataStatistic              (void);
 public:

    /*! \brief Store data to calorimeter.
        \arg data This data will be inserted.

    */
    void                StoreData              (const std::vector<CellDataRaw>& data);
    /*! \brief Store data to calorimeter in some not standard case. Mainly for Re-Reconstruction.
        \arg p Raw data from p will be extracted and inserted in Calorimeter.
    */
    void                StorebackData              (const std::vector<CalorimeterParticle>& p);

    /// Set map. If map is incorrect, nothing will be changed.
    void                SetMap                  (const std::map<size_t,size_t> &chan_to_cell);

    void                ReadMap                 (const std::string &file_name);
    /// X-check Mapping info
    virtual bool        XCheckMapping           ( void ) { return true;}

    const std::map<size_t,size_t> &
                        GetMap        (void) const {return map___chan_to_cell;}

    virtual void       Remap          (const std::vector<CellDataRaw> &data_in,std::vector<CellDataRaw> &data_out) const;

    /// Data recovery functions
    void               RemapCalibInternal  ( void );
    void               RecalculateCalibData  ( void );
    void               RecalculateTimeCalibData  ( void );
    void               StoreRawFromCalibData  ( void );
    void               RemapCalibCellInfo     ( CellInfoType what, CalibTimeType when );

    void                SubstractPED            (std::vector<CellDataRaw>& data) const;

    /// \return Total calibrated energy of all calorimeter cells - just usual user request
    double              GetEtotalCalib          (void) const;

    /// \return Total raw energy of all calorimeter cells - also usual user request
    double              GetEtotalRaw            (void) const;

  public:
    const std::vector<CellDataRaw>& GetSignals() const { return signals_; }
    std::vector<CellDataRaw>& GetSignals2Modify() { return signals_; }

  // ==========================================
  //  Noise Info Methods
  // ==========================================

  public:
    /*! \brief Noise Info Methods are in the development array
        Sorry, not well documented for the moment.
    */
    void                ClearCellsNoise         (int icell);
    void                StoreCellsNoise         (void);
    void                StoreCellsNoise4RandomTrigger         (void);
    void                StoreGammaNoise         (void);
    int                 InspectCellsNoise4RandomTrigger       (void) const;
    virtual int         InspectCellsNoise       (void);
    virtual int         InspectCellsGamNoise    (void);

  private:

    void                DebugChecksCellNoise      ( int icell ) const;
    double              GetCellNoiseStatEnties    ( int icell ) const;
    double              GetCellNoiseStatOverflow  ( int icell ) const;
    double              GetCellNoiseStatUnderflow ( int icell ) const;
    double              GetCellNoiseStatTotal     ( int icell ) const;
    double              GetCellNoiseStat          ( int icell, double emin, double emax) const;
    std::vector <double> GetCellNoiseStatAround   ( int icell, double emin, double emax) const;
    int                 InspectCellsNoiseContinue ( void );

  protected:
    const StatInfo     &GetCellsNoiseInfo       (size_t icell) const;
    void                PrintCellsNoiseInfo     (void) const;
    void                BookNoiseHisto          (void);

  // ==========================================
  //  More special Methods
  // ==========================================

    void                FillInternalCorrelations        (void);
    int                 InspectInternalCorrelations     (void) const;
    void                FillExternalCorrelations        (void);
    int                 InspectExternalCorrelations     (void) const;

  // ==========================================
  //  Call to some finalization methods before exit
  // ==========================================

  public:

    virtual void       EndOfJob                        (void);
    virtual void       ProcessPedestalsAtEndOfJob      (void);
    virtual void       FilterBadCellsOutOfStatisticLED (void);
    void                AddNewBadCell                   ( int icell, unsigned status );

  // ==========================================
  //  Manipulations with data base
  // ==========================================

  public:

    /// Reading from Data Base StatInfo  about calibration, leds and pedestals to OLD and NEW cells_info
    void                GetInfoFromDataBase       (const DataBase &db);

    /// Saving to Data Base information about NEW cells_info for calibration, leds and pedestals
    void                PutInfoToDataBase         (DataBase &db) const;

    /*! \brief Reading information from Data Base about NEW cells_info for calibration, leds and pedestals
               using time points
        \arg db  calibration data base
        \arg t   time point
    */
    void                ReadFromDataBase(DataBase &db,tm &t);

    /*! \brief Saving to Data Base information about NEW cells_info for calibration, leds and pedestals
               using time points
        \arg db  calibration data base
        \arg t_start   starting time for validity interval
        \arg t_finish  finished time for validity interval
    */
    void                WriteToDataBase(DataBase &db,tm &t_start,tm &t_finish);

  // ==========================================
  //  I/O functions
  // ==========================================

  public:
    /// I/O calibration info specifed by Calorimeter options
    /// Formatting information requested by tag

    virtual int InputAnyCalibInfo(const std::string &tag, const  std::string &s);

    /// \return available output info tags
    const std::vector< std::string > &GetOutputInfoNames( void ) const;
    const std::vector< std::string > &GetInputInfoNames( void ) const;
    void                              SetOutputInfoName( const std::string &tag );
    void                              SetInputInfoName( const std::string &tag );
    /// Formatting information requested by tag
    virtual int OutputAnyCalibInfo(const std::string &tag, std::string &s,const std::string &comment ="") const;

    /// Formatting information about cells_info statistic
    virtual int         InputStatInfo(const std::string &s, size_t what, size_t when);
    virtual int         OutputStatInfo(std::string& s, size_t what, size_t when) const;
    int                 InputStatInfo(const std::string &s,std::vector<StatInfo> &stat_info );
    int                 OutputStatInfo(std::string& s,const std::vector<StatInfo> &stat_info ) const;

    ///  Clear statistic information
    void                   ClearStatInfo(size_t what, size_t when);

    // I/O functions
    /// Formatting information about NEW cells_info for calibration
    virtual int            InputCalibInfo( size_t when, const std::string &s );
    virtual int            OutputCalibInfo( std::string &s) const;
    int                    InputCalibInfoXY(const std::string &s);
    int                    InputCalibInfoXY( size_t when, const std::string &s);
    int                    OutputCalibInfoXY( std::string &s) const;

    /// Formatting information about NEW cells_info for leds
    virtual int         InputLEDInfo( size_t when, const std::string &s);
    virtual int         OutputLEDInfo( std::string &s, const std::string &comment = "") const;
    virtual int         InputLEDInfoXY( size_t when, const std::string &s);
    virtual int         OutputLEDInfoXY( std::string &s, const std::string &comment = "") const;

    virtual int         InputLEDInfoInSpills( size_t when, const std::string &s);
//     virtual int         OutputLEDInfoInSpills( std::string &s, const std::string &comment = "") const;

    virtual int          InputLEDInfoInSpillsXY( size_t when, const std::string &s);
//     virtual int          OutputLEDInfoInSpillsXY( std::string &s, const std::string &comment = "") const;

//     /// read time in spill dependent LED correction factors
//     int                    InputTimeInSpillLED( const std::string &s );
//     /// write time in spill dependent LED correction factors
//     int                    OutputTimeInSpillLED( std::string &s ) const;

    /// read energy dependent correction factors
    void                InputEdepCorr(const std::string &s);
    /// read time in spill dependent correction factors
    void                InputTiSdepCorr(const std::string &s);


    /// Formatting information about NEW cells_info for pedestals
    virtual int         InputPEDInfo(const std::string &s);
    virtual int         OutputPEDInfo( std::string &s) const;
    virtual int         InputPEDInfoXY(const std::string &s);
    virtual int         OutputPEDInfoXY( std::string &s) const;

    /// Formatting information about NEW cells_info for time calibration
    virtual int         InputTimeCalibInfo( size_t when, const std::string &s);
    virtual int         OutputTimeCalibInfo( std::string &s) const;
    int                 InputTimeCalibInfoXY( size_t when, const std::string &s);
    int                 OutputTimeCalibInfoXY( std::string &s) const;

    /// Formatting information about bad cells
    virtual int         InputBadCellsInfo( size_t when, const std::string &s);
    virtual int         OutputBadCellsInfo( std::string &s) const;

    /// Formatting information about bad cells
    virtual int         InputEcutCellsInfo(const std::string &s);
    virtual int         OutputEcutCellsInfo( std::string &s ) const;
    int                 InputEcutCellsInfoXY(const std::string &s);
    int                 OutputEcutCellsInfoXY( std::string &s ) const;
    int                 InputProbEnT(const std::string &s);
    int                 OutputProbEnT( std::string &s ) const;

    /// Formatting information about bad cells
    virtual int         InputCellsInfo(const std::string &s,size_t what, size_t when);
    virtual int         OutputCellsInfo( std::string &s, size_t what, size_t when ) const;
    int                 InputCellsInfoXY(const std::string &s,size_t what, size_t when);
    int                 OutputCellsInfoXY( std::string &s, size_t what, size_t when ) const;

    int                  InputPositionInBurst(const std::string &s);
    int                  OutputPositionInBurst( std::string &s) const;
    ///  Accept generic format to get detector position in time
    int                  InputPositionInTimeGeneric(const std::string &s);
    int                  InputPositionInTime(const std::string &s);
    int                  InputPositionInTimeRedundant(const std::string &s);
    virtual int         InputJOUJOUInfo(const std::string &s);

    void                OperationsWithCalibInfo( void );
    virtual void        SpecialOperationsWithCalibInfo( void );

    virtual void        UpdateInternalAfterNewSettings( void );
    virtual void        UpdateFrontEndDependentSettings( void );

  // ==========================================
  //  Functions for Clusterisation and Gamma Search
  // ==========================================

  public:

    /// \return Cell number (in list GetCells()) with point (x,y,z) inside or return -1 if there is no such cell.
    virtual int         FindCell                (const CalorimeterParticle &p) const;

  protected:

    virtual int         FindCell                (const std::vector <double> &x) const;
    int                 FindCellInternal        (double x,double y,double z) const;
    virtual int         FindCellUltimate        (const std::vector <double> &x) const;

    /*! \brief Find all cells around point (x,y,z) on distance less or equal \b radius.

         Search for cells with distance from cell's center to point (x,y,z) less or equal \b radius.
         \b ATTENTION! Call FindCells(x,y,z,0) does NOT equal to FindCell(x,y,z).
    */
    virtual std::vector<size_t> FindCells            (double x, double y, double z, double radius) const;

  public:

    /// \brief  User suppled function with energy cut for reconstructed object (depends on variety of physics conditions)
    virtual double      UserEnergyCut            (double x, double y) const {return 0;}

  // ==========================================
  //  Calorimeter Geometry Functions
  // ==========================================

  public:

    /*! \brief Set Calorimeter coordinates to point (x,y,z) */
    void                SetPosition             (double x, double y, double z);
    /*! \brief Set Calorimeter coordinates at start of run. Used to monitor calorimeter movements. */
    void                SetPositionAtStartOfRun  (double x,double y,double z)
                                                { position_atstartofrun_[0]=x,
                                                   position_atstartofrun_[1]=y,
                                                    position_atstartofrun_[2]=z;}
    void                SetPositionAtStartOfRun  (void)
                                                { position_atstartofrun_[0]=position[0],
                                                   position_atstartofrun_[1]=position[1],
                                                    position_atstartofrun_[2]=position[2];}

    /*! \brief Set Calorimeter default vertex coordinates at point (x,y,z) */
    void                SetVertexPosition       (double x, double y, double z);

    /// \return Calorimeter X coordinate.
    double              GetPositionX            (void) const {return position[0];}

    /// \return Calorimeter Y coordinate.
    double              GetPositionY            (void) const {return position[1];}

    /// \return Calorimeter Z coordinate.
    double              GetPositionZ            (void) const {return position[2];}

    /// \return Calorimeter X coordinate.
    double              GetVertexPositionX      (void) const {return vertex_position[0];}

    /// \return Calorimeter Y coordinate.
    double              GetVertexPositionY      (void) const {return vertex_position[1];}

    /// \return Calorimeter Z coordinate.
    double              GetVertexPositionZ      (void) const {return vertex_position[2];}

    /// Left-most edge of all blocks (including padding)
    virtual double      GetXmin                 (void)  const;

    /// Left-most edge of all blocks (not including padding)
    virtual double      GetTrueXmin             (void)  const;

    /// Right-most edge of all blocks (including padding)
    virtual double      GetXmax                 (void)  const;

    /// Lower-most edge of all blocks (including padding)
    virtual double      GetYmin                 (void)  const;

    /// Lower-most edge of all blocks (not including padding)
    virtual double      GetTrueYmin             (void)  const;

    /// Top-most edge of all blocks (including padding)
    virtual double      GetYmax                 (void)  const;

    /// Most upstream edge (smallest z) of all blocks (not including padding)
    virtual double      GetZmin                 (void)  const;

    /// Most downstream edge (largest z) of all blocks (not including padding)
    virtual double      GetZmax                 (void)  const;

    /// \return Minimal Calorimeter Cell size in X
    double              GetMinCellSizeX         (void) const;
    double              GetMaxCellSizeX         (void) const;
    double              GetMinCellSizeY         (void) const;
    double              GetMaxCellSizeY         (void) const;
    double              GetMinCellSizeZ         (void) const;
    double              GetMaxCellSizeZ         (void) const;

    /// \return true of Calorimeter position was updated
    virtual bool        UpdatePositionInBurst    ( int busrst_number, const time_t  &time );
    bool                UpdatePositionInTime     (  const time_t  &time );
    void                StorePositionInfo        ( int burst_number, const time_t  &time, double x, double y, double w );
    void                StorePositionInfo( const std::vector< std::pair< StatInfo, StatInfo> > &diff_in_bursts );
    void                UpdateOnNewRun           ( int previous_run, int new_run );

 public:
    /// \return true if point MRS(x,y,z) is in active area of
    /// Calorimeter. (Nota bene: this routine is slow! But fairly general.)
    virtual std::pair<bool,double>  InActiveAreaExternal    ( double x, double y, double z, double dxdz, double dydz )  const;
    virtual bool        InActiveArea            ( const std::vector<double> &x )  const;
    virtual std::pair<bool,double>  InActiveAreaExternal    ( const std::vector<double>& vpar, int syst )  const;
    virtual bool        Associated    (const CalorimeterParticle &pe, const std::vector<double>& vpar, int syst, double gate )  const;

    /// \return true if point MRS(x,y,z) is in active area of Calorimeter and not on boundary.
//    virtual bool        WellInActiveArea        (double x,double y,double z)  const;
     virtual bool        WellInActiveArea        (const std::vector<double> &x)  const;

    /// \return true if point MRS(x,y,z) is on the Calorimeter boundary.
//     virtual bool        OnBoundary              (double x,double y,double z)  const;
    virtual bool        OnBoundary              (const std::vector<double> &x)  const;

    virtual bool        GoodZone4Electrons       ( double x, double y, double e)  const {return true;}

  // ==========================================
  //  GUI
  // ==========================================

  public:

    /*! \brief Create GUI interface
    */
    virtual void         CreateGUI               ( void );
    virtual void         UpdateGUI               ( void );

    const
    GUICalorimeter     &GetGUI                  (void) const {if(gui==NULL) throw "Reco::Calorimeter::GetGUI() no gui"; return *gui;}
    GUICalorimeter     &GetGUI                  (void)       {if(gui==NULL) throw "Reco::Calorimeter::GetGUI() no gui"; return *gui;}

  // ==========================================
  //  Visualisation and other GUI methods
  // ==========================================

  public:

    /// Clear statistic in histograms
    virtual void       ResetHisto              (void);

    /// Clear statistic
    virtual void       ResetStatistic          (void);

    /// Print canvas to PS file
    virtual void       PrintPS                 (void);

    /*! \brief Show some test histograms
    */
    virtual void       ShowEvent               (void);

    /// Draw histograms for one cell
    virtual void       DrawCellHisto           (int icell);

    /// Draw monitoring histograms for one cell
    virtual void       DrawMonitorCellHisto    (int icell);

    /// Draw in spills monitoring histograms for one cell
    virtual void       DrawMonitorCellHistoInSpill     (int icell);

    /// Draw LED histograms for several cells
    virtual void       DrawCellsHistoLED       (int from_cell,int how_much,int step);

    /// Draw Energy histograms for several cells
    virtual void       DrawCellsHistoGAM       (int from_cell,int how_much,int step);

    /// Fit LED histograms for several cells
    virtual void       FitCellsHistoLED       (int from_cell,int how_much,int step);

    /// Fit LED histograms for several cells
    virtual void       FitCellsHistoGAM       (int from_cell,int how_much,int step);

    /// Fit histograms for several cells
    virtual void       DrawFitCellsHisto      (CellInfoType info_type,const char *opt,int from_cell,
                                                  int how_much,int step,double fmin=0.,double fmax=0.);

    /// New Histo Fitting for all cells
    void       FitAllCellsHisto      (CellInfoType info_type,const char *opt );

 private:
    /// New Histo Fitting for one cell, private not to check all parameters
    void       FitCellHisto      (int cell, CellInfoType info_type,const char *opt );
  public:

    /// GUI with XY access
    virtual void       DrawFitCellsHistoXY   (CellInfoType info_type,const char *opt,
                                                                   int xmin,int xhow_much,int xstep,
                                                                   int ymin,int yhow_much,int ystep,
                                                                         double fmin=0.,double fmax=0.);
    virtual void      DrawFitCellsHistoXYSubset   (int sub_set, CellInfoType info_type,const char *opt,
                                                                   int xmin,int xhow_much,int xstep,
                                                                   int ymin,int yhow_much,int ystep,
                                                                         double fmin=0.,double fmax=0.);
    /// Draw monitoring histograms for several cells
    virtual void       DrawMonitorHisto      (CellInfoType info_type,const char *opt,int from_cell,
                                                  int how_much,int step,double fmin=0.,double fmax=0.);
  // ==========================================
  //  Methods ...
  // ==========================================

  public:

    /// \return Number of Calorimeter cells.
    size_t              NCells                  (void) const {return cells.size();}
    /// List of all cell's types
    const std::list<CellType>  & GetCellsTypeList (void) const {return  cells_type;}

    /// \return cell's name.
    virtual  std::string    GetCellName        (int icell) const;
    /// \return -1 in case cell name is not found otherwice \return cell index
    int      FindCellByName        (const std::string &cell_name);

    /// map XY to cell number, asserts XYRegularGrid()
    /// \return -1 if XY is in the hole or outside the calorimeter
    int                     GetCellOfColumnRow     (int x, int y) const;
    /// map cell number to X, asserts XYRegularGrid()
    int                     GetColumnOfCell        (int icell) const;
    /// map cell number to Y, asserts XYRegularGrid()
    int                     GetRowOfCell           (int icell) const;
    /// \return number of columns, asserts XYRegularGrid()
    int                     GetNColumns            ( void ) const { assert( XYRegularGrid() ); return fNcols_; }
    /// \return number of rows, asserts XYRegularGrid()
    int                     GetNRows               ( void ) const { assert( XYRegularGrid() ); return fNrows_; }

    /*! \brief  Set hited cell info to particles. This method sets nothing to Calorimeter.
        Aimed to recover this information from calorimeter and set it to particles
    \arg particles - vector of CalorimeterParticle
    */
    void                    SetHitedCells           (std::vector<CalorimeterParticle> &particles);
    bool                    SetHitedCells           (CalorimeterParticle &particle);
    void                    SetHitedCellsUltimate   (CalorimeterParticle &particle);

   // Set EventID info
    void UpdateEventID(int ev_type, int run_num, int spill_num, int ev_in_spill, int ev_in_run, const unsigned int& trigger_mask,
                       const time_t &time, const double time_in_spill, bool is_led_event );
    void UpdateEventID( const EventID &evid );
    const EventID           &GetEventID(void ) const;
    virtual double           GetTimeInSpill ( void ) const { return GetEventID().GetTimeInSpill();}

    /*! \brief \return cell's response. Interface DEVELOPMENT method !!! Not any checks applied !!!
                       This is first implementation of cell's response function. Not so much fit OO design :(
        \arg jcell - cell's index
        \arg ehit  - hit's energy
        \arg xhit  - hit's x-coordinate (in Calorimeter RS)
        \arg yhit  - hit's y-coordinate
        \arg zhit  - hit's z-coordinate
    */
    virtual  double     CellResponse            (int jcell,double ehit,double xhit,double yhit,double zhit) const
                                                                                                {return ehit;}
    /// \return true if all calorimeter blocks have the same size
    bool                  XYRegularGrid        (void) const { return !options.mixed_blocks; }

    /// Does the same as XYRegularGrid(), for backward compatibility (PHAST) only.
    bool                  XYInitialized        (void) const { return !options.mixed_blocks; }

 protected:
    /// Try to initialise XY structure in calorimeter
    void                  InitXY                  (void);


  // ==========================================
  //  Other Methods ...
  // ==========================================

  private:

    /// Determine the cells neighbors
    void                DetectCellNeighbors     (void);

    /// Determine the cells on boundary and boundary regions in selected cells
    void                DetectCellsOnBoundary   (void);

  public:

    /// Check cell index
    bool                CellIsValid        (int cell ) const { return ( cell >= 0 && (unsigned)cell < NCells() ); }

  private:

    /// \return vector of gamma indexes around the ibcell
    std::vector<size_t > GetGammasAroundCell  ( size_t ibcell ) const;

    bool  CheckBoundaryX4MakeProfiles( double xb, const std::vector < double > &boundary,
                                                        const std::vector<size_t> &long_list ) const;
    bool  CheckBoundaryY4MakeProfiles( double xb, const std::vector < double > &boundary,
                                                        const std::vector<size_t> &long_list ) const;
  public:

    TDirectory* GetHistogramsBaseDir()                const { return histograms_basedir_; }
    void        SetHistogramsBaseDir(TDirectory* dir)       { histograms_basedir_ = dir;  }

  // ==========================================
  //  Attributes, data
  // ==========================================

  protected:

    /// The calorimeter's name.
    std::string                              name;

    /// The calorimeter's position in MRS.
    double                              position[3];
    /// The calorimeter's position in MRS at the begining of the run... Well if it moves and soft is not adopted..
    double                              position_atstartofrun_[3];
    /// The calorimeter's position in MRS at the begining of the run... Well if it moves and soft is not adopted..
    double                              position_correction_[3];

    /// Default vertex position in MRS.
    double                              vertex_position[3];

    /// Channel is key (first element in the pair) and cell number is data (second element).
    std::map<size_t,size_t>             map___chan_to_cell;

    /// List of all cell's types
    std::list<CellType>                     cells_type;

    /// Calorimeter's cells.
    std::vector<Cell>                        cells;

    /// vector of signals used during reconstruction
    std::vector<CellDataRaw>                 signals_;

    /// vector of MC truth signals
    std::vector<CellDataRaw>                 mc_input_;

    /// Calorimeter's XY structure
    std::vector<int>                        map_cell_xy_;
    std::vector<int>                        map_xy_cell_[2];
    int                                     fNcols_;
    int                                     fNrows_;
    double                                  fXStep_;
    double                                  fYStep_;

    /*!  \brief cells information

	 CellInfoType            {CALIB=0,LED=1,PED=2,TIME=3,NOISE=4,CLUSTER=5,RAW=6};
	 CalibTimeType           {OLD=0,NEW=1,MONITOR=2,PRIM=3,MC=4};

	 cells_info[CALIB][OLD]     - calibrations which we read from DB
	 cells_info[CALIB][NEW]     - calibrations which we collecting on-line
	 cells_info[CALIB][MC]      - calibrations which we apply to MC amplitudes
	 cells_info[CALIB][PRIM]    - calibrations references
	 cells_info       [MONITOR] - special usage for monitoring

	 cells_info[LED][OLD]  - leds which we read from DB
	 cells_info[LED][NEW]  - leds which we collecting on-line
	 cells_info[LED][PRIM] - leds references

    */
    std::vector<StatInfo>                     cells_info[CellInfoTypeSize][CalibTimeTypeSize];

    // Some stuff for investigations of internal correlations
    //
    // Is used in Calorimeter::FillInternalCorrelations file
    // ReconstructionTest to fill some cell to cell correlations. For examle
    // in the search of bad mapping.
    // Used: in the test/check purpose
    // Option: options.fill_internal_correlations

    double                                    stat_total_internal_correlations_;
    std::vector<StatInfo>                      internal_correlations_all_info_;
    std::vector<StatInfo>                     *internal_correlations_info_;

    // Some stuff for investigations of external correlations
    //
    // Is used in Calorimeter::FillInternalCorrelations file
    // ReconstructionTest to fill some cell correlations with external
    // particles.  To store external particles(MC for example) hint particles
    // storage is used.
    // Used: in the test/check purpose
    // Option: options.fill_external_correlations

    double                                    stat_total_external_correlations_;
    std::vector<StatInfo>                      external_correlations_all_info_;
    std::vector<StatInfo>                     *external_correlations_info_;

    /// Vector of fit parameters, fit_OK flags etc used to store fit results
    std::vector<FitInfo>                      fit_info[CellInfoTypeSize];

    /// Noise stuff
    bool                                noise_histo_booked;

    /// Hint particles
    std::vector<CalorimeterParticle> hint_particles;

    /// Reconstructed particles
    std::vector<CalorimeterParticle> reco_particles;

    /// calorimeter's SubSets
    std::vector<SubSet> calorimeter_sub_sets;
    /// calorimeter's SubSets indexes of trigger groups
    std::vector<size_t> trigger_groups_index_;

    /// GUI calorimeter object
    mutable GUICalorimeter*     gui;

    /// Options for calorimeter
    mutable Options             options;

    /// Object for histograms code testing.
    CalorimeterHist*            calo_hist;
    InternalHisto*              internal_hist_;
    ProfilesHisto*              profiles_hist_;
    CorrelationsHisto*          internal_correlations_hist_;
    CorrelationsHisto*          external_correlations_hist_;
    TrigGroupHist*              trig_group_hist_;
    MakeParticlesHisto*         make_particles_histo_;
    FitInfoHisto*               fit_info_histo_;

    std::list <size_t>          bad_cells_[CalibTimeTypeSize];

    /// Array of bad cells
    std::vector<bool>           cell_is_bad_[CalibTimeTypeSize];

    /// cells and gamma noise statistic
    double                      cells_noise_statistic;
    double                      cells_gamma_noise_statistic;
    double                      led_statistic_;
// ??
    std::map<long unsigned int,double>     *leds_in_spills_old_;

// leds store
    StoreLED                              *store_led_;
    DevMonitorFEM                    *fem_;

    ///  Array of Bad pedestals
    int*                        bad_ped_;

    /// List of cells with bad pedestals
    std::list <size_t>          bad_ped_cells_;

    /// Array of bad cells status
    std::vector<int>              status_bad_cells_[CalibTimeTypeSize];

    /// Array of individual cells thresholds which is in use
    double                     *energy_cut_sparse_mode;

    /// Array of individual cells thresholds which is in use
    double                     *energy_cut_bad_cells_old;

    /// Array of individual cells thresholds, new one
    double                     *energy_cut_bad_cells_new;

    /// Array of individual gamma thresholds which is in use
    double                     *energy_gamma_cut_bad_cells_old;

    /// Array of individual gamma thresholds, new one
    double                     *energy_gamma_cut_bad_cells_new;

    /*!  \brief Alternative way to store cells information (calibration,LED,PED,TIME,NOISE)
         In the development array.
    */

    //  Summ of statistic
    StatInfo                   *cells_stat_info[CellInfoTypeSize][CalibTimeTypeSize];

    StatInfoStore              *cells_store[CellInfoTypeSize][CalibTimeTypeSize];
    StatInfoStore               LEDTotal;
    StatInfoStore               CALIBTotal;
    std::vector<double>        prob_cut_;
    std::vector<double>        ebin_stat_;
    /// ProbEnT
    double                    *energy_cut_[11];
    // Factor for individual cell characteristics like (e/mju ratio)
    std::vector<double>        individual_calib_factor_;

    std::vector < size_t > subsets_list2store_calib_histo_;

    SummaryInspectLED  summary_led_;
    SummaryInspectLED  prev_summary_led_;


    std::map< int, int>         map2remap_;
    std::map< int, std::vector< double> >             map_burst2position_;
//    std::vector< std::pair< time_t, std::vector< double> >  >          map_time2position_;
    AlignmentInTime              alignment_in_time_;


    /// X Position of moving detector during calibration run
    std::vector<StatInfo>                   positionX_in_burst_new_;
    /// Y Position of moving detector during calibration run
    std::vector<StatInfo>                   positionY_in_burst_new_;

    bool                                    final_initialization_;

  protected:
    EventID                                 evid_;
    bool                                    evid_updated_;

  protected:
    std::vector<double>                     preshower_corrections_;

    Reco::MCConstruction*                   monteconstruction_;
    bool                                    monteconstruction_cleared_;

    Reco::Reconstruction*                   reconstruction_;

    TDirectory*                             histograms_basedir_;

    // from CalorimeterMatrixSet
 protected:
    std::list <CellsMatrix> 	matrixes;

    /*! \brief First cell number from detectors.dat file for given cell matrix

        The size of this list is equal to size of Calorimeter::matrixes.
    */
    std::vector<size_t>      comgeant_first_cell_n_;

  private:
    /*! \brief map to speed up association of cell name to cell number
     */
    std::map<std::string, int> name_to_id_;

 protected:
    std::vector<TProfile *>      p1_led_in_spills_;
    TCanvas              *       c_led_in_spills_;

};

////////////////////////////////////////////////////////////////////////////////

} // using namespace Reco

#undef inside_Calorimeter_h
#endif // Calorimeter___include
