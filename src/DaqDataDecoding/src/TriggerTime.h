#ifndef CS__TriggerTime____include
#define CS__TriggerTime____include

#include <map>
#include <set>
#include <string>

#include "config.h"
#include "Exception.h"
#include "DetID.h"
#include "ChipF1.h"
#include "Stat.h"

/*! \file TriggerTime.h
    \author Alexander Zvyagin
*/

namespace CS {

class DaqErrors;
class DaqOption;
class Trigger;
class TriggerTimeShift;
class TriggerTimeStat;


/*! @brief Event independent data.

    This object should be kept between events.
    A good place for its storage is the DaqOption class.
*/
class TriggerTimeConfig
{
  public:

    class           TTC
    {
      public:

      TTC (const DetID &det     = DetID(""),
	   double _time_unit    = 0,
	   int TT_index         = -1,
	   int TT_index_recover = -1,
	   int _channels        = -1,
	   int TT_overolling    = 0,
	   double _sigma        = 0) :
	id(det),
	srcID(0),
	port(-1),
	index(TT_index), 
	index_recover(TT_index_recover),
	channels(_channels),
	time_unit(_time_unit),
	overolling(TT_overolling),
	sigma(_sigma),
	MT_shift(0) {}

        bool        operator ==         (const TTC &d) const;
        bool        operator !=         (const TTC &d) const {return !(*this==d);}
        void        Print               (const char *prefix="") const;

        DetID       id;
        uint16      srcID; // expected sourceID
        int         port;
        int         index;
        int         index_recover;
        int         channels;
        double      time_unit;
        int         overolling;
        double      sigma;
        double      MT_shift;
    };  // class TTC

    TriggerTimeConfig (void) :
      TCS_phase_DetID(""),
      TCS_phase_srcID(0),
      TCS_phase_port(-1),
      TCS_phase_channel(-1),
      trigger_mask_DetID(""),
      trigger_mask_srcID(0),
      trigger_mask_port(-1),
      trigger_mask_sources_min(0),
      trigger_mask_sources_max(0),
      trigger_mask_TT_index(-1),
      time_jitter(0.2),
      time_precise(1),
      triggers_diff_sigma(2) {}
		    
    void            Print               (const std::string &prefix="") const;

    void            Add                 (const TTC &ttc);

    const TTC *     Find                (double time_unit) const {return const_cast<TTC*>(Find(time_unit));}
          TTC *     Find                (double time_unit);

    void            SetTriggerMaskDataLimit (int min,int max);

    std::map<int,TTC>  tt;
    DetID              TCS_phase_DetID;
    uint16             TCS_phase_srcID;
    int                TCS_phase_port;
    int                TCS_phase_channel;
    DetID              trigger_mask_DetID;
    uint16             trigger_mask_srcID;
    int                trigger_mask_port;
    unsigned           trigger_mask_sources_min;
    unsigned           trigger_mask_sources_max;
    int                trigger_mask_TT_index;
    float              time_jitter;
    float              time_precise;
    float              triggers_diff_sigma;
    std::set<std::string> tis_tbnames;  //< scalers from which time in spill is read
};

class TriggerTimeShift
{
  public:
                        TriggerTimeShift(uint16 _bit,uint32 tr, float _mt,float _tt_mt_diff )
                                        : bit(_bit),trigger(tr),mt(_mt),tt_mt_diff(_tt_mt_diff) {}
    uint16  bit;                // from the above pattern: 0 or 2
    uint32  trigger;            // Event trigger pattern:  000101
    float   mt;                 // master time (in counts, not [ns])
    float   tt_mt_diff;         // trigger_time - master_time
};

/*! \brief COMPASS Trigger Time measurements

    TriggerTime class holds measurement of trigger times (TT) of COMPASS experiment.
    There are several of them. They are indexed by an integer number, have different
    time resolution, different names and different time overollings.

    \author Alexander Zvyagin 
*/
class TriggerTime
{
  public:
    
    /*! @brief Data of a single TT.
    */
    struct Data
    {
                        Data                (const TriggerTimeConfig::TTC *_config=NULL) : config(_config) {Clear();}
        void            Clear               (void) { data.clear(); decoded=false; time=0; }
        const TriggerTimeConfig::TTC & GetConfig (void) const {if(config==NULL) throw Exception("TriggerTime::GetConfig()==NULL"); return *config;}

        double          time;
        bool            decoded;
        std::map<Chip::DataID,uint16> data;
        const TriggerTimeConfig::TTC *config;
    };
  
  public:
                        TriggerTime         (void) :
			  config(NULL),
			  ltt(-1),
			  event_trigger(0),
			  trigger_of_master_time(NULL),
			  TCS_phase(-1e6),
			  TimeInSpill(-1.)
			  {}

  public:
  
    /// Decode trigger time information
    void                Decode              (const DaqOption &options,DaqErrors &errors);

    /// Decode trigger time information
    void                Decode              (const Chip::Digits& digits,const DaqOption &options,DaqErrors &errors);
  public:

    /// Clear TT data, be ready for next data entries
    void                InitAndClear        (const TriggerTimeConfig &ttc);
    
    /// Decode trigger time information and subtract it from TT.
    void                DecodeAndSubtract   (Chip::Digits& digits,const DaqOption &options,DaqErrors &errors);

    /// Find the Data which holds given time unit \c tu.
    const TriggerTime::Data *Find           (double tu) const;

    /// Find the Data which holds given time index
    const TriggerTime::Data *Find           (int index) const;

    /// Get time for given index.
    double              GetTime             (int index) const;
    
    bool                HasPhaseTCS         (void) const {return !TCS_phase_data.empty();}
    float               GetPhaseTCS         (void) const;

    bool                HasTimeInSpill      (void) const { return TimeInSpill >= 0.; }
    double              GetTimeInSpill      (void) const;
    
    /*!  \return -1 if it was not decoded (no measuremnts)
         \return -2 if the measurements were different

        From Bernhard Ketzer' mail (2003-09-19):
        
        ...a useful feature was 
        implemented in the ADC headers sent by GEM and Silicon: The time since 
        the last trigger in units of a 19.44 MHz clock (half the standard 
        COMPASS clock of 38.88 MHz) is stored in the first 24 bits of the second 
        word of the ADC header (which was already labelled 'time', but used to 
        be the time since the last reset).

        Could you make this time available in the event loop? It is NOT a GEM or 
        Si feature, but rather a TCS value, and should be treated as such in 
        CORAL or COOOL. In principle it has to be the same for all ADC headers 
        of a given event, but it may be useful to check it. It is very useful to 
        verify the dead time of the experiment, and other things.
    */
    int                 GetLastTriggerTicks (void) const {return ltt;}
    
    bool                SetLastTriggerTicks (int t);
    
    const std::vector<TriggerTimeShift> & GetShifts (void) const {return shifts;}

    void                GetFiredTrigger     (const std::multimap<float,const Trigger*> &triggers,
                                             std::multimap<float,const Trigger*>::const_iterator &trig,
                                             float &mt_correction);
    
    /*! \return trigger which is responsable for MT calculation */
    const Trigger *     GetTriggerMT        (void) const {return trigger_of_master_time;}
    
  private:

    void                CorrectMasterTime   (const DaqOption &options,DaqErrors &errors);

    bool                CheckShiftTT        (const TriggerTimeStat &tt_stats,const DaqOption &options,DaqErrors &errors) const;
    
  private:
   
    const TriggerTimeConfig *config;
    std::map<int,Data>  data;       // indexed by time_unit
    int                 ltt;        // last triger ticks
    std::vector<const ChipF1::Digit*> trigger_mask_data;
    std::vector<const ChipF1::Digit*> TCS_phase_data;
    std::vector<TriggerTimeShift> shifts;
    uint32              event_trigger;
    const Trigger *     trigger_of_master_time;
    float               TCS_phase;
    double              TimeInSpill;
};

class TriggerTimeStat
{
  public:
    
                        TriggerTimeStat         (void);
    void                Add                     (unsigned bit,float tt_mt_diff);
    void                Add                     (unsigned bit1,unsigned bit2,float diff);
    void                Print                   (void) const;
    unsigned            GetFillsCounter         (void) const {return fills_counter;}
    
    const
    std::map<uint32,Stat> &  GetTriggersStat    (void) const {return triggers_stat;}
    std::map<uint32,Stat> &  GetTriggersStat    (void)       {return triggers_stat;}
    
    void                SetBufferSize           (unsigned s);
    unsigned            GetBufferSize           (void) const {return buffer_size;}
    
    void                SetBufferCutRaitio      (float r);
    float               GetBufferCutRaitio      (void) const {return buffer_cut_ratio;}
    
    unsigned            GetCheckOnFill          (void) const {return check_on_fill;}
    void                SetCheckOnFill          (unsigned a) {check_on_fill=a;}

    unsigned            GetCheckEntriesMin      (void) const {return check_entries_min;}
    void                SetCheckEntriesMin      (unsigned a) {check_entries_min=a;}
    
    float               GetAllowedShift         (void) const {return allowed_shift;}
    void                SetAllowedShift         (float v) {allowed_shift=v;}
    
  private:
  
    unsigned            buffer_size;
    float               buffer_cut_ratio;
    unsigned            check_entries_min;
    unsigned            check_on_fill;
  
    float               allowed_shift;

    unsigned            fills_counter;
    std::map<uint32,Stat>    triggers_stat;
};

}

#endif // CS__TriggerTime____include
