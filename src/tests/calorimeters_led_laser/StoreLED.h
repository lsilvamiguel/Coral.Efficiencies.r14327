#ifndef StoreLED___include
#define StoreLED___include

#include "Reco/Calorimeter.h"
#include "Reco/StatInfo.h"

class HistoInSpillLED;

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////

class CaloEvent
{
  public:
     CaloEvent ( size_t ncells, const Reco::EventID &event_id ) : event_id_(event_id)  { for( size_t i=0; i<ncells; i++ ) data_.push_back(0.);}
  public:
    Reco::EventID  event_id_;
    std::vector <double>   data_;
};

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////

class StatInSpill
{
  public:
     StatInSpill ( void ) : ev_cnt(0) {}
  public:
   int    ev_cnt; 
   size_t    evmin; 
   size_t    evmax; 
   time_t    tmin;
   time_t    tmax;
};

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////

class StoreCaloEvents
{
  public:
//    static std::map < size_t,  std::vector<double> > time_start_run_map;
    const static double spill_gate;
    const static double spill_tcycle;
  public:
     StoreCaloEvents ( Reco::Calorimeter *c );
//      void InitStatic ( void );   
     void AddEvent ( CaloEvent  * ev);
     void Process ( void );   
     bool Processed ( void ) const { return   ((vev_id_.size() >0)&& (led_no_cuts.size() > 0)); }
     size_t  NCells ( void ) const { return ncells_;}
     void SetSpillsInfo ( size_t run, size_t spill ); 
     const std::pair<  int, int> & GetSpillsInfo( size_t run ) const;

     const Reco::StatInfo &  GetValue( size_t run, size_t spill, size_t cell ) const;
     const Reco::StatInfo &  GetStrictAverageValue ( size_t cell ) const;
     
//      static std::pair< bool,std::pair< double, double>  >  GetSpillTime ( size_t run, size_t spill );

  public:
    Reco::Calorimeter   * c_;
    size_t  run_min_;
    size_t  run_max_;
    size_t  spill_min_;
    size_t  spill_max_;
    std::map <  std::pair< size_t, size_t >,  std::vector< CaloEvent  * > >  events_;

    std::map <  size_t, std::pair<  int, int >  > spill_min_max_;
    std::map <  std::pair< size_t, size_t >,  StatInSpill >  statinspill_;
    std::vector < Reco::EventID > vev_id_;

    std::vector<Reco::StatInfo> led_no_cuts;
    std::vector<Reco::StatInfo> led_thr;
    std::vector<Reco::StatInfo> led_strict;
    std::map <  std::pair< size_t, size_t >,  std::vector<Reco::StatInfo> > leds_in_spill;
    Reco::StatInfo zero;
    std::pair<  int, int >  dummy_spill_info_;
    size_t   ncells_;
};

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////

class StoreLED
{
  public:
              StoreLED ( Reco::Calorimeter *c  );
              void  AddEvent ( bool off_spill );
              void  Process ( void );
              size_t NCells ( void ) const { return c_->NCells(); } 
              int     InputTimeInSpillLED(const std::string &s);
              int    OutputTimeInSpillLED( std::string &s ) const;
              int    OutputLEDInfoInSpills( std::string &s, const std::string &comment ) const;
              int    OutputLEDInfoInSpillsXY( std::string &s, const std::string &comment ) const;
             void   FillHistoInSpillLED( size_t icell, double framp, double frspill, bool allok );
             void   FitHistoInSpillLED( void );
  public:
    Reco::Calorimeter * c_;
// leds store
    StoreCaloEvents                 *store_leds_all_;
    StoreCaloEvents                 *store_leds_in_spill_all_;
    std::vector<Reco::StatInfo> ledinspill0_;
    std::vector<Reco::StatInfo> ledinspill1_;
    std::vector<Reco::StatInfo> ledinspill2_;

    std::vector <double> ledtis0_;
    std::vector <double> ledtis1_;
    std::vector <double> ledtis2_;
    
    HistoInSpillLED       *histo_in_spill_led_;   
};

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#endif // StoreLED___include
