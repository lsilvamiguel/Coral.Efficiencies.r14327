#ifndef CompassSoft_DaqOptions__include
#define CompassSoft_DaqOptions__include

#include <map>
#include <set>
#include <vector>

#include "DaqMap.h"
#include "Chip.h"
#include "Trigger.h"
#include "TriggerTime.h"

namespace CS {

struct hufftree_s;
typedef struct hufftree_s hufftree_t;

////////////////////////////////////////////////////////////////////////////////

/*! \brief Options used in the decoding and mapping process.
    \author Alexander Zvyagin
*/
class DaqOption: public DaqMap
{
  //============================================================================
  // Types, constants
  //============================================================================

  public:
  
  //============================================================================
  // Constructors and destructor
  //============================================================================

  public:
    
                       ~DaqOption               (void);

    /*! \brief Base constructor.
    */
                        DaqOption               (void);

                        DaqOption               (const DaqOption &o) : DaqMap(ObjectXML("")),chip_SADC_ht(NULL) {*this=o;}

                        DaqOption               (const ObjectXML &o);

  //============================================================================
  // Operators
  //============================================================================

  public:
  
    DaqOption &         operator =              (const DaqOption &);
   
  //============================================================================
  // Methods
  //============================================================================

  public:
  
    void                ReadXML                 (const ObjectXML &o);
  
    /// Print options to output stream
    virtual void        Print                   (std::ostream &o=std::cout,const std::string &prefix="") const;
    
    void                Clear                   (void);
    
    const std::vector<Chip::CreationRule>& GetChipCreationRules(void) const {return chip_creation_rules;}
          std::vector<Chip::CreationRule>& GetChipCreationRules(void)       {return chip_creation_rules;}
    
    const std::set<uint16>&  GetSrcIDScan       (void) const {return srcID_scan;}
          std::set<uint16>&  GetSrcIDScan       (void)       {return srcID_scan;}
    
    /*! @brief @return true if a given srcID should be scanned.*/
    bool                NeedScanSrcID           (uint16 srcID) const;
    
    bool                Use_ERROR_ECC_checksum  (void) const {return use_ERROR_ECC_checksum;}

    bool                Use_TDC_WRPLL           (void) const {return use_TDC_WRPLL;}
    
    bool                Use_SLinkFormat         (void) const {return use_SLinkFormat;}

    bool                FixSLinkMultiplexerSize (void) const {return fixSLinkMultiplexerSize;}
    
    int                 GetOnlineFilterSrcID    (void) const {return online_filter_srcID;}
    void                SetOnlineFilterSrcID    (int v) {online_filter_srcID=v;}
    
    const std::set<int>& GetLastTriggerTicksSrcIgnore(void) const {return last_trigger_ticks_srcID_ignore;}
    const std::set<int>& GetLocalTriggerTicksSrcIgnore(void) const {return local_trigger_ticks_srcID_ignore;}
    
    void                SetCalibration          (const Chip::DataID &data_id,const Chip::Calibration &calib);
    const Chip::Calibration*  FindCalibration   (const Chip::DataID &data_id) const;

    bool                IsMasterTimeSrcID       (uint32 n) const {return master_time_srcID.count(n);}
    bool                SkipMasterTimeCheckSrcID(uint32 n) const {return master_time_srcID_dont_check.count(n);}

    void                ReadChipSADC_HT         (const CS::ObjectXML&);
    const hufftree_t *  GetChipSADC_HT          (void) const {return chip_SADC_ht;}
    void                SetChipSADC_HT          (hufftree_t *t) {chip_SADC_ht=(hufftree_t*)t;}

    const
    std::set<uint16> &  GetAlwaysCMC4APV        (void) const {return cmc_apv_always;}
    std::set<uint16> &  GetAlwaysCMC4APV        (void)       {return cmc_apv_always;}

    const
    std::set<uint16> &  GetSADCDataVer2         (void) const {return sadc_data_ver2;}
    std::set<uint16> &  GetSADCDataVer2         (void)       {return sadc_data_ver2;}

    const
    std::set<uint16> &  GetSADCDataVer3         (void) const {return sadc_data_ver3;}
    std::set<uint16> &  GetSADCDataVer3         (void)       {return sadc_data_ver3;}
    
    const
    std::set<uint16> &  GetEventSrcIDs          (void) const {return event_srcIDs;}
    std::set<uint16> &  GetEventSrcIDs          (void)       {return event_srcIDs;}

    const
    std::set<uint32> &  GetDataF1CMC            (void) const {return f1_cmc_data;}
    std::set<uint32> &  GetDataF1CMC            (void)       {return f1_cmc_data;}
    void                AddDataF1CMC            (uint32 srcID,uint32 port) {f1_cmc_data.insert((srcID<<16)+port);}
    bool                IsDataF1CMC             (uint32 srcID,uint32 port) const {return f1_cmc_data.count((srcID<<16)+port);}
    
    const
    TriggerTimeConfig & GetTTConfig             (void) const {return tt_config;}
    TriggerTimeConfig & GetTTConfig             (void)       {return tt_config;}
    
    unsigned            GetTriggerBit           (const std::string &name) const;
    const Trigger &     GetTrigger              (uint16 bit) const {return const_cast<DaqOption*>(this)->GetTrigger(bit);}
    Trigger &           GetTrigger              (uint16 bit);

    const
    std::set<uint16> &  GetTriggerBitsIgnore    (void) const {return trigger_bits_ignore;}
    std::set<uint16> &  GetTriggerBitsIgnore    (void)       {return trigger_bits_ignore;}
    
    const
    std::map<uint16,Trigger> & GetTriggers      (void) const {return triggers;}
    std::map<uint16,Trigger> & GetTriggers      (void)       {return triggers;}
    
    const
    TriggerTimeStat &   GetStatTT               (void) const {return tt_stats;}
    TriggerTimeStat &   GetStatTT               (void)       {return tt_stats;}

    const
    std::map<uint16,std::set<uint32> > & GetSrcIDPresentPorts (void) const {return srcID_present_ports;}
    void                       AddSrcIDPresentPorts (uint16 srcID, uint32 port);

  //============================================================================
  // Attributes
  //============================================================================

  private:

    bool                        use_ERROR_ECC_checksum;

    bool                        use_TDC_WRPLL;
    
    bool                        use_SLinkFormat;
    
    bool                        fixSLinkMultiplexerSize;

    /// Which chip type should be created (ChipF1,ChipADC,...)
    std::vector<Chip::CreationRule>  chip_creation_rules;

    // List of srcIDs to be present in a event
    std::set<uint16>            event_srcIDs;
    
    /*! @brief  List of ((srcID<<16)+geoID) which are known to have problems (we should not add data from those srcIDs). */
    std::set<uint32>            srcID_geoID_bad;

    /// List of srcID which shoud be scanned. If the list is empty, all srcID are scaned.
    std::set<uint16>            srcID_scan;

    /// Triggers
    std::map<uint16,Trigger>    triggers;

    int                         online_filter_srcID;
    
    /// List of cathes to be ignore in last trigger time calculations.
    std::set<int>               last_trigger_ticks_srcID_ignore;
    
    /// List of cathes to be ignore in last trigger time calculations across a single GeSiCA.
    std::set<int>               local_trigger_ticks_srcID_ignore;
    
    std::map<Chip::DataID,const Chip::Calibration*> calibrations;
    
    std::set<uint16>            master_time_srcID;

    /*! List of srcID with should not be checked for TT. */
    std::set<uint16>            master_time_srcID_dont_check;

    /// SADC Huffman Decoding data.
    hufftree_t *                chip_SADC_ht;
    
    /*! List of srcID (ChipAPV) for which common mode noise correction words always present in data.*/
    std::set<uint16>            cmc_apv_always;
    
    /*! List of srcID (ChipSADC) with format2 (2006) of the data. */
    std::set<uint16>            sadc_data_ver2;

    /*! List of srcID (ChipSADC) with format3 (2008) of the data. */
    std::set<uint16>            sadc_data_ver3;
    
    /*! Data ((catch<<16)+port) with set F1-CMC flag. */
    std::set<uint32>            f1_cmc_data;
    
    TriggerTimeConfig           tt_config;      // indexed by time_unit
    
    TriggerTimeStat             tt_stats;
    
    std::set<uint16>            trigger_bits_ignore;

    /*! Ports (or ADC IDs) that must be present (currently only used for APV equipment). */
    std::map<uint16,std::set<uint32> >    srcID_present_ports;
};

////////////////////////////////////////////////////////////////////////////////

} // namespace CS

#endif // CompassSoft_DaqOptions__include
