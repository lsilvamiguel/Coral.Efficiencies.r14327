/*
   $Source: /afs/cern.ch/user/g/gobbo/w0/coral/src/calorim/Reco/src/StatInfo.h,v $
   $Date: 2011/01/31 20:35:45 $
   $Revision: 1.11 $
   -------------------------------------------------------------------------

   This file is part of cellular calorimeter reconstruction program.

   Authors:
     Vladimir  Kolosov   ( Vladimir.Kolosov@cern.ch, Kolosov@mx.ihep.su )
     Alexander Zvyagin   ( Alexander.Zviagine@cern.ch, Zvyagin@mx.ihep.su )

   Copyright(C): 2002  V.Kolosov,A.Zvyagin

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

#ifndef RecoStatInfo___include
#define RecoStatInfo___include

#include <string>
#include <vector>
#include <iostream>
#include <list>
#include <cstdlib>   // for exit()

// --- ROOT include files ----
#include "TROOT.h"
#include "TH1.h"
#include "TH2.h"
#include "TDirectory.h"
#include "TCanvas.h"
#include "TF1.h"

#include "myROOT_utils.h"

namespace Reco {

/// Generic class to hold statistical information (count, mean, sigma) of a
/// physical quantity.  The data is stored as 0th, 1st and 2nd moment so that
/// it is easy to add more and more statistcs.  The class is commonly used to
/// during calibration, but also to read in calibration data from disk.
class StatInfo
{
  //============================================================================
  // Types, constants
  //============================================================================

 private:

  static const size_t momenta_size=3;

  //============================================================================
  // Constructors and destructor
  //============================================================================

 public:

  /// Default constructor
  StatInfo                (void);

  /// Constructor with given values
  StatInfo                (double weights_sum,double mean,double sigma);

  /// Copy constructor
  StatInfo                (const StatInfo &c) {*this=c;}

  //============================================================================
  // Operators
  //============================================================================

 public:

  /// Assignment operator
  StatInfo           &operator =              (const StatInfo &c);

  StatInfo           &operator +=             (const StatInfo  &);

  friend std::ostream     &operator <<       (std::ostream &o, const Reco::StatInfo &c);
  friend std::istream     &operator >>       (std::istream &in, Reco::StatInfo &c);

  //============================================================================
  // Methods
  //============================================================================

 public:

  /// Clear all data.
  void                Clear                   (void);

  /// Add an amplitude
  void                Add                     (double a, double weight=1);

  /// Set Cell Info
  void                Set                     (double weights_sum, double mean, double sigma);

  /// \return entries amount
  double              GetEntries              (void) const {return momenta[0];}

  /// \return average
  double              GetMean                 (void) const;

  /// \return sigma
  double              GetSigma                (void) const;

  /// return momentum of \b n-th order
  double              GetMomenta              (size_t n) const { if(n>=momenta_size) throw "StatInfo::GetMomenta()  bad number"; return momenta[n]; }

  void                Result                  (double &n,double &mean,double &sigma) const;

  double*             GetMomenta              (void)       {return momenta;}

  const double*       GetMomenta              (void) const {return momenta;}

  /// Add constant shift
  void                Translate               (double c);

  /// Scale on constant fraction
  void                Scale                   (double c);

  //============================================================================
  // Attribute
  //============================================================================

 private:

  /// momenta
  double       momenta[momenta_size];                // 0,1,2,3...
};

////////////////////////////////////////////////////////////////////////////////

class StatInfoStore
{
  //============================================================================
  // Types, constants
  //============================================================================

   const static int                      maximum_intervals = 100;

  public:

    /// Object for StatInfoStore histograms
    class Monitor_Histo
    {
      public:

      /// ROOT directory
      TDirectory                         *root_dir;

      /// Monitoring histogram
      TH1D                               *h1_monitor;

      /// Monitoring entries histogram
      TH1D                               *h1_monitor_entries;

      /// Monitoring sigma histogram
      TH1D                               *h1_monitor_sigma;

      /// Monitoring histogram with relative normalisation to average
      TH1D                               *h1_monitor_norm;
    };

  //============================================================================
  // Constructors and destructor
  //============================================================================

  public:

    /// Default constructor
                        StatInfoStore           (void);

    /// Base constructor
                        StatInfoStore           ( int run_start, int run_finish );

    /// Constructor with data setting
                        StatInfoStore           (int run_start, std::vector<StatInfo> &stat_info);

  //============================================================================
  // Operators
  //============================================================================

  public:

    /// Assignment operator
    StatInfoStore       &operator =             (const StatInfoStore &c);

    StatInfoStore       &operator +=            (const StatInfoStore &c);

  //============================================================================
  // Methods
  //============================================================================

  public:

    /// Initialization
    void                Init                    ( int run_start, int run_finish );

    /// Store range
    int                 Size                    ( void ) const { return  run_finish_ - run_start_ +1;}

    /// Clear all data.
    void                Clear                   (void);

    /// Add data.
    void                Add                     (int run, StatInfo &stat_info);

    /// Add data.
    void                Add                     (int run, double data);

    /// Check is run in range
    bool                RunInRange               (int run) const { return ( run >= run_start_ &&
                                                                           run <= run_finish_);}

    /// Set all data in one go
    void                SetInfo                 (int run_start, std::vector<StatInfo> &stat_info);

    /// Set all data for particular run
    void                SetInfo                 (int run, StatInfo &stat_info);

    /// \return StatInfo vector
    std::vector<StatInfo>   &GetStatInfo        (void) {return stat_info_;}

    /// \return list of points
    std::list<int>         &GetIntervals        (void) {return points_;}

    /// Set list of points
    void                SetIntervals            (std::list<int> &points);

    /// \return array of parameters
    double             *GetParameters           (void) {return parameters_;}

    /// \return number of intervals
    int                 NIntervals              (void) const { return  points_.size()-1;}

    /// Split in intervals
    void                FindIntervals           (int max_intervals ,double hi_break, double hi_penalty );

    /// Fit in intervals
    double             FitInIntervals           (void) { return FitInIntervals( points_, parameters_ ); }

    /// Calculate average StatInfo
    StatInfo            CalculateAverage        (void);

    /// Calculate average StatInfo
    StatInfo            CalculateAverage        ( int run_min, int run_max );

    /// Calculate average in interval
    StatInfo            CalculateAverageInInterval ( int run );

    /// Calculate average in interval
    StatInfo            GetBinValue             ( int run ) const;

    /// Re-calculate StatInfo vector with relative normalization to vector average
    StatInfo            MakeRelativeNormalizationToAverage        (void);

    /// Book monitoring histograms
    void                BookHisto               (const std::string &name,
                                                 const std::string &histo_path, int histo_level);

    /// Fill monitoring histograms
    void                FillHisto               (void);

    /// Draw monitoring histograms
    void                DrawHisto               (void);
    Monitor_Histo*      GetMonitorHisto          (void)  {return  monitor_hist_;}

    /// Add fitted curves in intervals to monitoring histograms
    void                AddFitToHisto           (void);

    /// Set minimal range
    void                SetMinimalRange         (double r) { min_range_ = r; }

    /// Set minimal statistic
    void                SetMinimalStat          (double stat) { min_stat_ = stat; }

    /// Set minimal sigma (additive)
    void                SetMinimalSigmaA        (double s_a) { min_disp_add_ = s_a * s_a; }

    /// Set minimal sigma (fraction)
    void                SetMinimalSigmaF        (double s_f) { min_disp_frac_ = s_f * s_f; }

    /// Set minimal range, minimal statistic, minimal sigma (additive), minimal sigma (fraction)
    void                SetTogether             (double r, double stat, double s_a , double s_f )
                                                  { SetMinimalRange(r); SetMinimalStat(stat);
                                                    SetMinimalSigmaA( s_a ); SetMinimalSigmaF( s_f );}

    /// Set minimal range, minimal statistic, minimal sigma (additive), minimal sigma (fraction)
    double              CalculateFitValue       ( int run );

    /// Print stored and fitted info
    void                Print                   ( void );

  private:

    /// Allocate memory for data
    void                Init                    ( void );

    /// Fit in interval
    double              FitInInterval           (int from_bin, int to_bin, double *parameters);

    /// Fit in interval
    double              FitInIntervals          (std::list<int> &points, double *parameters);

    /// Add new point
    double              AddNewPoint             (std::list<int> &points, double *parameters, double penalty);

    /// Add new interval
    double              AddNewInterval          (std::list<int> &points, double *parameters, double penalty);

    /// Calculate statistic average StatInfo
    double              EstimateSigma           (int from_bin, int to_bin );

  //============================================================================
  // Attribute
  //============================================================================

  private:

    /// First run
    int                                         run_start_;

    /// Last run
    int                                         run_finish_;

    /// StatInfo vector
    std::vector<StatInfo>                       stat_info_;

    /// List with intervals, each element points at the begining of interval
    std::list<int>                              points_;

    /// Array to store fit parameters
    double                                      parameters_[200];

    /// Object for StatInfoStore histograms
    Monitor_Histo*                              monitor_hist_;

    /// Minimal range width
    double                                      min_range_;

    /// Minimal statistic to consider
    double                                      min_stat_;

    /// Minimal additional dispersion term
    double                                      min_disp_add_;

    /// Minimal additional relative dispersion term
    double                                      min_disp_frac_;

    ///
    bool                                        data_updated_;

};

////////////////////////////////////////////////////////////////////////////////

} /// using namespace Reco

std::ostream &operator << (std::ostream &o ,const std::vector<Reco::StatInfo> &v);
std::istream &operator >> (std::istream &in,      std::vector<Reco::StatInfo> &v);

#endif //RecoStatInfo___include

#ifndef RecoEventID___include
#define RecoEventID___include
namespace Reco {
/*! \brief Event ID Info - pilot implementation for general usage??? we shall see..
                            Initialization at decoding stage.
                   Ported from MN::EventID with very slight modifications
*/
class EventID
{
  public:

    /// Destructor
     virtual           ~EventID (void){}

       EventID ( void ) : initialized_(false),
                          event_type_(0),
                          run_number_(-1),
                          spill_number_in_run_(-1),
                          event_number_in_spill_(-1),
                          event_number_in_run_(-1),
                          trigger_mask_(0),
                          time_in_spill_(-1.),
//                           processed_events_(0),
//                           processed_events_in_spill_(0),
//                           processed_events_in_run_(0),
			   is_led_event_(false)  {}
    /// Copy constructor
                       EventID   (const EventID &id) {*this=id;}

  //============================================================================
  // Operators
  //============================================================================

  public:

    /// Assignment operator
     EventID          &operator =              (const EventID &id) { if( &id!=this )
                                                                       {
                                                                         initialized_ = id.initialized_;
                                                                         event_type_ = id.event_type_;
                                                                         setup_event_type_ = id.setup_event_type_;
                                                                         run_number_ = id.run_number_;
                                                                         spill_number_in_run_ = id.spill_number_in_run_;
                                                                         event_number_in_spill_ = id.event_number_in_spill_;
                                                                         event_number_in_run_ = id.event_number_in_run_ ;
                                                                         trigger_mask_ = id.trigger_mask_;
//                                                                          processed_events_ = id.processed_events_;
//                                                                          processed_events_in_spill_ = id.processed_events_in_spill_;
//                                                                          processed_events_in_run_ = id.processed_events_in_run_;
                                                                         time_ = id.time_;
                                                                         time_in_spill_ = id.time_in_spill_;
                                                                         is_led_event_ = id.is_led_event_;
                                                                       }
                                                                       return *this;
                                                                     }


  //============================================================================
  // Methods
  //============================================================================
       void Update ( int ev_type, int run_num, int spill_num,int ev_in_spill,
                             int ev_in_run, const unsigned int& trigger_mask,
                             const time_t &time, const double time_in_spill,
                             bool is_led_event )
                                                                             {
//                                                                                  std::cout << " EventID::Update " << std::endl;
                                                                                 initialized_=true;
                                                                                 event_type_= ev_type;
										 bool new_run=false;
										 if( run_number_ != run_num ) new_run=true;
                                                                                 run_number_=run_num;
										 bool new_spill=false;
										 if( spill_number_in_run_ != spill_num ) new_spill=true;
                                                                                 spill_number_in_run_=spill_num;
                                                                                 event_number_in_spill_= ev_in_spill;
                                                                                 event_number_in_run_= ev_in_run;
                                                                                 trigger_mask_ = trigger_mask;
                                                                                 time_ = time;
                                                                                 time_in_spill_ = time_in_spill;
                                                                                 is_led_event_ = is_led_event;
//                                                                                  processed_events_++;
// 										 if( new_spill )
// 										   processed_events_in_spill_ = 1;
// 										 else
//                                                                                    processed_events_in_spill_++;
//
// 										 if( new_run )
// 										   processed_events_in_run_ = 1;
// 										 else
//                                                                                    processed_events_in_run_++;
                                                                             }

      bool  Initialized ( void ) const { return initialized_;}

      const int  &GetType ( void ) const { if( !initialized_ )
                                               {std::cerr << " GetType EventID::Not Init " << std::endl; exit(1);}
                                             return  event_type_;
                                           }

      const int  &GetUniType ( void ) const { if( !initialized_ )
                                                  {std::cerr << " GetUniType EventID::Not Init " << std::endl; exit(1);}
                                                return  setup_event_type_;
                                              }

      const int  &GetRunNumber ( void ) const { if( !initialized_ )
                                                    {std::cerr << "GetRunNumber EventID::Not Init " << std::endl; exit(1);}
                                                  return  run_number_;
                                                }

      const int  &GetEventNumberInRun ( void ) const { if( !initialized_ )
                                                           {std::cerr << " GetEventNumberInRun EventID::Not Init " << std::endl; exit(1);}
                                                         return  event_number_in_run_;
                                                       }

      const int  &GetEventNumberInBurst ( void ) const { if( !initialized_ )
                                                             {std::cerr << " GetEventNumberInBurs EventID::Not Init " << std::endl; exit(1);}
                                                           return  event_number_in_spill_;
                                                         }

      const int  &GetBurstNumber ( void ) const { if( !initialized_ )
                                                      {std::cout << "GetBurstNumber EventID::Not Init " << std::endl; exit(1);}
                                                    return  spill_number_in_run_;
                                                  }

      const unsigned int  &GetTriggerMask ( void ) const { if( !initialized_ )
                                                             {std::cerr << " GetTriggerMask EventID::Not Init " << std::endl; exit(1);}
                                                           return  trigger_mask_;
                                                         }

//       const int  &GetProcessedEventNumber ( void ) const { if( !initialized_ )
//                                                                {std::cout << " GetProcessedEventNumber EventID::Not Init " << std::endl; exit(1);}
//                                                              return  processed_events_;
//                                                            }

      const time_t  &GetTime ( void ) const { if( !initialized_ )
                                                  {std::cout << " GetTime EventID::Not Init " << std::endl; exit(1);}
                                                return time_;
                                              }

      const double  &GetTimeInSpill ( void ) const { if( !initialized_ )
                                                  {std::cout << " GetTime EventID::Not Init " << std::endl; exit(1);}
                                                return time_in_spill_;
                                              }

      bool           IsLEDEvent ( void ) const { if( !initialized_ )
                                                     {std::cout << " IsLEDEvent EventID::Not Init " << std::endl; exit(1);}
                                                   return is_led_event_;
                                                 }

      void        Print ( void ) const { if( !initialized_ )
                                             {std::cout << " EventID::Not Init " << std::endl;}
                                           {std::cout << " EventType " << event_type_ << std::endl;}
                                           {std::cout << " SetupEventType " << setup_event_type_ << std::endl;}
                                           {std::cout << " RunNumber " << run_number_ << std::endl;}
                                           {std::cout << " EventNumberInRun " << event_number_in_run_ << std::endl;}
                                           {std::cout << " EventNumberInSpill " << event_number_in_spill_ << std::endl;}
                                           {std::cout << " SpillNumber " << spill_number_in_run_ << std::endl;}
                                           {std::cout << " TimeInSpill " << time_in_spill_ << std::endl;}
//                                            {std::cout << " ProcessedEvents " << processed_events_ << std::endl;}
//                                            {std::cout << " ProcessedEvents in Run" << processed_events_in_run_ << std::endl;}
//                                            {std::cout << " ProcessedEvents in Spill" << processed_events_in_spill_ << std::endl;}
                                         }

      void        SetSetupEventType ( int type ) { setup_event_type_ = type; }
      int         GetSetupEventType ( void ) const { return setup_event_type_; }

  private:
    bool        initialized_;
    int    event_type_;
    int    setup_event_type_;
    int    run_number_;
    int    spill_number_in_run_;
    int    event_number_in_spill_;
    int    event_number_in_run_;
    unsigned int trigger_mask_;
    time_t time_;
    double time_in_spill_;
    bool   is_led_event_;
//    std::pair < time_t, size_t> time_;   // from now on we missed usec
//   public:
//     int    processed_events_;
//     int    processed_events_in_run_;
//     int    processed_events_in_spill_;
};
} /// using namespace Reco
#endif //RecoEventID___include
