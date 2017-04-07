#ifndef CompassSoft_Chip__include
#define CompassSoft_Chip__include

#include <iostream>
#include <string>
#include <map>
#include <vector>
#include <list>
#include <set>
#include <cstdlib>
#include "config.h"
#include "DaqError.h"
#include "SLink.h"
#include "DetID.h"
#include "DaqMap.h"
#include "version.h"

// Old versions of ROOT define this.
#ifdef Check
#undef Check
#endif

namespace CS {

class DaqOption;
class DaqEvent;

////////////////////////////////////////////////////////////////////////////////

/*! \brief Class with description of basic properties of a DAQ chip data.

    This class gives you access to structure of data from a COMPASS
    hardware chip. Each COMPASS chip is constructed from memory buffer
    with a DATE equipment. So every chip has corresponding DATE equipment
    header. See Chip::Header.
    
    Some documentation was taken from "Compass-Note 2000-8" and "Compass-Note 2001-8".
    See http://wwwcompass.cern.ch/compass/notes

    \author Alexander Zvyagin
    \author Wolfgang Kastaun
*/
class Chip
{
  //============================================================================
  // Types
  //============================================================================

  public:

    enum {SourceID_MAX=1024};

    typedef uint64 DataID;

    class Maps;
    
    /*! \brief Chip map for DAQ connection to a real detector
        
        This map describes how the chip's channals are connected to detector's channels.
    */
    class Map: public DaqMap
    {
      public:
        
        /// The destructor.
        virtual        ~Map                     (void) {}

      protected:

        /// Base constructor
                        Map                     (const ObjectXML &o);

      public:

        /*! Global Chip::Map class method for a Map creation.
        */
        static
        Chip::Map*      Create                  (const ObjectXML &o);

      public:

        /// Print the map definition
        virtual void    Print                   (std::ostream &o=std::cout,const std::string &prefix="") const;

      public:
        
        /// \return source ID for the chip.
        uint16          GetSourceID             (void) const {return source_id;}

        /// \return Detector identification correspoinding to this map
        const DetID&    GetDetID                (void) const {return id;}

        /// \return Detector identification correspoinding to this map
              DetID&    GetDetID                (void)       {return id;}
              
        /// \return first chip channel number that will be used for mapping
        int32           GetChanF                (void) const {return chanF;}

        /// \return last chip channel number that will be used for mapping
        int32           GetChanL                (void) const {return chanL;}

        /// \return chip channels step
        int32           GetChanS                (void) const {return chanS;}

        /// \return first detector channel number that will be used for mapping
        int32           GetWireF                (void) const {return wireF;}

        /// \return last detector channel number that will be used for mapping
        int32           GetWireL                (void) const {return wireL;}

        /// \return detector channels step
        int32           GetWireS                (void) const {return wireS;}

        /// \return number of channels
        uint32          GetChanN                (void) const {return 1+(chanL-chanF)/chanS;}
        
        /// Calculate wire number for this channel or throw exception (int)0 if it is not possible.
        int32           CalcWire                (int32 chan) const;
        
        /// Check the map
        virtual void    Check                   (void);

        /// Add the current map to the list of maps.
        virtual void    AddToMaps               (Maps &maps,DaqOption &opt) const {}
        
        bool            IsMultiDigit            (void) const {std::string s; GetAttribute("multi_digit",s); return s=="yes";}

      protected:

        std::string     map_line;
        std::string     dec_line;

        DetID           id;
        uint16          source_id;

        // channel_XXX  is input

        int32           chanF;
        int32           chanL;
        int32           chanS;
        
        // wire_XXX   is output
        
        int32           wireF;
        int32           wireL;
        int32           wireS;

      friend class Chip::Maps;
    };

    ////////////////////////////////////////////////////////////////////////////

    /// Chip creation
    class CreationRule
    {
      public:
                        CreationRule            (const std::string &s);

        void            Print                   (std::ostream &o=std::cout,const std::string &prefix="") const;

        bool            Accept                  (int32 _Dtype,int32 _Did,int32 _Stype,int32 _Sid,int32 _format) const;
        bool            Accept                  (int32 _Stype,int32 _Sid,int32 _format) const;

        int32           Dtype;                  ///< DATE equipment type (use -1 to accept all types)
        int32           Did;                    ///< DATE equipment id (use -1 to accept all IDes)
        int32           Stype;                  ///< S-Link type (use -1 to accept all types)
        int32           Sid;                    ///< S-Link ID (use -1 to accept all IDes)
        int32           format;                 ///< Format mask (use -1 to accept all formats)
        char            srcID;                  ///< 's' for Slink and 'e' for DATE equipment
        
        std::string     chip_name;              ///< Name of the chip (ChipF1, ChipAPV, ChipADC, ...)
    };

    ////////////////////////////////////////////////////////////////////////////
    
    /// Chip Digit.
    class Digit
    {
      public:
  
        /// Destructor
        virtual        ~Digit                   (void) {}

        /// Base constructor
                        Digit                   (const DataID &d,const DetID &id) : data_id(d), det_id(id) {}

        /// Print Digit data to output stream.
        virtual void    Print                   (std::ostream &o=std::cout,const std::string &prefix="") const;

        /// \return Data format for ntuple.
        virtual const char*   GetNtupleFormat   (void) const {return "";}
        
        /// \return Vector of data.
        virtual std::vector<float> GetNtupleData (void) const {return std::vector<float>();}

        /*! \return detector ID corresponding to this digit

            If the detector ID is unknown for this digit, then the detector name
            will be an empty string "".
        */
        const DetID&    GetDetID                (void) const {return det_id;}
        void            SetDetID                (const DetID &d) {det_id=d;}
        
        const DataID&   GetDataID               (void) const {return data_id;}
        void            SetDataID               (const DataID &id) {data_id=id;}

      protected:
    
        /// The data identification.
        DataID          data_id;

        /// Detector identification corresponding to this digit.
        DetID           det_id;
    };

    ////////////////////////////////////////////////////////////////////////////

    /// List of digits    
    class Digits : public std::multimap<DetID,Digit*>
    {
      public:
      
                       Digits                   (void) {}

        /// Destructor.
                      ~Digits                   (void) {Clear();}


      private:

                       Digits                   (const Digits &d);

        Digits&        operator =               (const Digits&);

      public:

        /// Clear the digits list.
        void            Clear                   (void) {for(std::multimap<DetID,Digit*>::iterator it=begin(); it!=end(); it++) delete it->second; std::multimap<DetID,Digit*>::clear();}
        
        /// Clear the digits list (the same as Clear()).
        void            clear                   (void) {Clear();}

        /// Print all digits.
        void            Print                   (std::ostream &o=std::cout,const std::string &prefix="") const {for(std::multimap<DetID,Digit*>::const_iterator it=begin(); it!=end(); it++) it->second->Print(o,prefix);}
    };

    ////////////////////////////////////////////////////////////////////////////

    class Maps : public std::multimap<DataID,Digit*>
    {
      public:
                        Maps                    (void) {}
        virtual        ~Maps                    (void);
        void            Print                   (std::ostream &o=std::cout,const std::string &prefix="") const;

        void            Clear                   (void);
        void            clear                   (void) {Clear();}

        /// \return Number of wires for given dtector name \b det_name in a run number \b run_number
        uint32          GetWires                (const DetID &id) const;
        
        void            SetWires                (const DetID &id,uint32 wires) {channels[id]=wires;}
        
        /// Find all SourceIDs for given DetID
        void            GetSrcIDs               (const DetID& id,std::set<uint16> &srcIDs) const;

        /// Find all SourceIDs for given regular expression pattern
        void            GetSrcIDs               (const std::string &pattern,std::set<uint16> &srcIDs) const;
        
        void            PrintSrcIDs             (const DetID& id) const;

        const std::map<DetID,std::set<uint16> >& GetDetSrcIDes(void) const {return det_srcIDes;}
              std::map<DetID,std::set<uint16> >& GetDetSrcIDes(void)       {return det_srcIDes;}
        
      private:
                        Maps                    (const Maps&);
        Maps&           operator =              (const Maps&) const;
    
      private:

        std::map<DetID,std::set<uint16> > det_srcIDes;
        std::map<DetID,uint16>       channels;
    };

    ////////////////////////////////////////////////////////////////////////////
    
    class ErrorWord
    {
      public:
        
                        ErrorWord               (uint32 word,uint32 error_marker,uint32 source_id,Chip &chip,DaqOption &o);

        ErrorWord&      operator =              (uint32 d) {data.all=d; return *this;}

        uint16          GetMarker               (void) const {return data.ew.error_marker;}
        uint16          GetGeographicID         (void) const {return data.ew.geo_id;}
        uint16          GetPort                 (void) const {return data.ew.port;}
        uint16          GetErrorCode            (void) const {return data.ew.error_code;}
        uint16          GetSpecialErrorCode     (void) const {return data.ew.serror_code;}
        bool            IsError                 (void) const {return is_error;}

      private:

        union
        {
          struct
          {
            uint32      error_code   :  4;
            uint32      serror_code  :  2;
            uint32      port         :  4;
            uint32      geo_id       : 10;
            uint32      error_marker : 12;
          } ew;
          uint32 all;
        } data;

        bool is_error;
    };

    class Calibration
    {
      public:
        virtual ~Calibration(void) {}
    };
    
  //============================================================================
  // Static methods.
  //============================================================================

  public:
    
    /*! \brief Create Chip from a memory buffer.

        \param buffer pointer to memory buffer with a DATE equipment
    
        There are many different chips. For example we
        can have chips A,B,C (these are all classes derived
        from Chip). And there is a memory buffer with
        a DATE equipment of one of the above types. The method
        Chip::Create() will detect the chip type
        corresponding to the DATE equipment on given address of memory and
        will create appropriate C++ object. Because it is assumed
        that all chips will be inherited from Chip, it
        is possible to return a Chip* pointer to the new object.

        Chip::Create() method never returns NULL pointer.
        If it is not possible to create
        a chip (format of buffer is not recognized) the exception 
        DaqError will be thrown.

        Example:
        \code
          #include "DaqError.h"
          // ...
          try  // Start try-catch block
          {
            Chip *eq = Chip::Create(buf,false,rules);
            eq->Print();
            // ...
            delete eq; // Release the memory.
          }
          catch( const DaqError &e )
          {
            // An error! Print message.
            cerr << "Exception: " << e.what() << "\n";
          }
          catch( ... )
          {
            cerr << "Unknown exception.\n";
          }
        \endcode


        \param copy_buf if \b true then data will be copied
                        from \c buffer. Otherwise only reference to
                        \c buffer will be copied.

        With \c copy_buf = \b false the code will be fast, with
        \c copy_buf = \b true the code will be more robust.
        
        Example:
        \code
          CreationRules *rules=....
          char *buf=get_next_equip();  // 'buf' contains next event
          Chip *c1 = Chip::Create(buf,false,*rules);
          Chip *c2 = Chip::Create(buf,true,*rules);
            buf[1]=buf[1]+10; // Damage the memory!
            // c1 is broken
            // c2 is fine
            // ...
            delete c1;   // free memory
            delete c2;   // free memory
        \endcode
        
        \attention Do not forget to \b delete this \b new Chip object!
    */
    static
    Chip               *Create                  (void const * const buffer,bool copy_buf,DaqOption &opt,DaqEvent &ev);

    /*! \brief Read chips maps from the input stream to the map argument.
        \attention New maps will be added to the arguments.
        
        \arg path file or directory with xml files.
        
        \arg dets For every detector name it will be assigned number. This number will be set in DetID object.
        The last argument vector<string> will have list of all names in accordance to thus numbering.
    */
    static
    bool                ReadMaps                (uint32 run,const std::string &path,Maps &maps,
                                                 DaqOption &option,std::vector<std::string> &dets);
    
    /// Construct a data ID from 16-bits + 8-bits + 8-bits + 8-bits + ...
    static DataID       CreateDataID7           (uint16 v1,uint8 v2,uint8 v3, uint8 v4,uint8 v5,uint8 v6,uint8 v7)
                                                {
                                                  return (uint64(v1)<<48)|
                                                         (uint64(v2)<<40)|
                                                         (uint64(v3)<<32)|
                                                         (uint64(v4)<<24)|
                                                         (uint64(v5)<<16)|
                                                         (uint64(v6)<< 8)|
                                                                 v7;
                                                }

    /// Construct a data ID from 16-bits + 16-bits + 8-bits + 8-bits + ...
    static DataID       CreateDataID5           (uint16 v1,uint16 v2,uint16 v3, uint8 v4,uint8 v5)
                                                {
                                                  return (uint64(v1)<<48)|
                                                         (uint64(v2)<<32)|
                                                         (uint64(v3)<<16)|
                                                         (uint64(v4)<< 8)|
                                                                 v5;
                                                }

  private:

    static void         RegisterMap             (Map *m,Maps &maps,std::vector<std::string> &dets,DaqOption &options);

  //============================================================================
  // Constructors and destructor
  //============================================================================

  public:
    
    /// Destructor
    virtual            ~Chip                    (void) {if(flag_buffer_delete) free(buf);}
    
    /*! \brief Construct an chip from memory buffer.

        \param buffer is non-NULL pointer to buffer starting with SLink

        \param copy_buf if \b true then chip data will be copied
                        from \c buffer. Otherwise only reference to
                        \c buffer will be copied.

        \exception DaqError if a chip can not be created.
    */
                        Chip                    (void const * const buffer,bool copy_buf,DaqOption &opt,DaqEvent &ev);
    
//    /// Create Chip.
//                        Chip                    (void);

  private:

    /// Copy constructor is not implemented.
                        Chip                    (const Chip &e);
  
  //============================================================================
  // Operators
  //============================================================================

  private:

    /// Assignment operator is not implemented.
    Chip               &operator =              (const Chip &e);

  //============================================================================
  // Static methods
  //============================================================================
  
  public:
    
    /// Trigger Time Contol system Frequency
    static double       GetTCSFrequency         (void) {return 38.88e6;}

    /*! \brief Set/clear flag to collect some information.
        
        When Chip base constructor will be called, some statistics will be
        saved in static attributes of class Chip. Use method
        Chip::PrintStatistics() to print it. Default condition is to collect such
        information.
        \sa PrintStatistics()
        
        Example:
        \code
          ...
          try
          {
            Chip::SetFillStatistics(true);
            char *buf=...;
            while(some_cond)
            {
              Chip chip(buf);
              ....
            }
          }
          catch(...)
          {
            cerr << "An error!\n";
          }
        \endcode
          
        Chip::PrintStatistics(cout,"Chips statistics: ");
    */
    static void         SetFillStatistics       (bool f) {fill_stat_flag=f;}

    /*! \brief Print equipment statistics.
        \sa SetFillStatistics().
    */
    static void         PrintStatistics         (std::ostream &o=std::cout,const std::string &prefix="");
    
  //============================================================================
  // Methods
  //============================================================================

  public:
  
    /// Print chip information
    virtual void        Print                   (std::ostream &o=std::cout,const std::string &prefix="") const;

    /// Print chip all data
    virtual void        Dump                    (std::ostream &o=std::cout,const std::string &prefix="") const;
    
    /// Clear the chip
    virtual void        Clear                   (void) {geoIDs.clear();is_scaned=false;}

    /*! \return Chip source identification.
    
        If there is S-Link header for this chip the method return S-Link source id, othervise
        it returns DATE equipment id from DATE header.
    */
    uint16              GetSourceID             (void) const {return GetSLink().GetSourceID();}
    
    /// Set sourceID.
    void                SetSourceID             (uint16 n) {GetSLink().SetSourceID(n);}

    /// \return equipment link
    const SLink        &GetSLink                (void) const {return *reinterpret_cast<SLink*>(buf);}

    /// \return equipment link
          SLink        &GetSLink                (void)       {return *reinterpret_cast<SLink*>(buf);}

    /// \return Full chip data length in bytes.
    unsigned            GetLength               (void) const {return GetSLink().GetEventSize()*4;}
    
    virtual
    const uint32*       GetDataStart            (void) const {return buf + sizeof(SLink)/4;}
    
    virtual
    const uint32*       GetDataEnd              (void) const {return buf+GetLength()/4;}
    
    /// Decode data and \b add new digits to \b digits_list.
    virtual void        Decode                  (const Maps &maps,Digits &digits_list,DaqOption &opt) {}
    virtual void        DecodeCalibEvent        (const Maps &maps,DaqOption &opt,...) const {}
    
    /// \return chip name that created this equipment
    virtual std::string GetName                 (void) const {return "ChipUnknown";}
    
    /// \return number 4
    virtual uint16      DataSize                (void) const {return 4;}
    
    /// \return true if you can modify the chip data
    bool                CanModify               (void) const {return flag_can_modify;}
    
    /*! \brief Add 'n' of 32-bits words of data from buffer to the chip
    
        \exception DaqError If you do not have rights to modify the buffer.
    */
    void                Add                     (const void *buffer,size_t n);
    
    virtual void        AddMap                  (const ObjectXML &o,Maps &maps) {}

    /// Add error to the list and take an appropriate action (if needed).
    void                AddError                (const DaqError &e,DaqOption &o);

    /// Set the owner of the chip
    void                SetDaqEvent             (DaqEvent &e) {event=&e;}
    
    bool                IsScaned                (void) const {return is_scaned;}

  protected:

    /*! \brief Scan chip data, detect errors.

        This method should be reimplemented by all derived chip classes.
    */
    virtual void        Scan                    (DaqOption &opt);
    
    void                RegisterGeoID           (uint16 geoID) { if( geoIDs.size()==0 || geoIDs.back()!=geoID )
                                                                 geoIDs.push_back(geoID); }
    
  //============================================================================
  // Attributes, data
  //============================================================================

  private:
  
    /// It is used only for initialisation.
//    static const bool   init;

  private:

    /// Do we need to delete the buffer in destructor?
    bool                flag_buffer_delete;
    
    /// True if user con modify the buffer. If flag_can_modify==true then flag_buffer_delete==true also.
    bool                flag_can_modify;
    
  protected:

    /// Pointer to real data. 
    uint32             *buf;
    
    bool                is_scaned;
    
    /// The owner of the chip.
    DaqEvent*           event;
    
    /// Chip's geoIDs
    std::vector<uint16> geoIDs;
    
    /// Set it to \b true to fill statistic information.
    static bool         fill_stat_flag;
    
    /*! \brief Global equipment information

        First key is equipment ID, second one is amount of this equipment.
    */
    static std::map<int,int> map__slink_source_id;

    static std::map<int,int> map__slink_type;

    static std::map<int,int> map__slink_format;
    
    static std::map<std::pair<int,int>,int> map__src_fmt;
};

////////////////////////////////////////////////////////////////////////////////

} // namespace CS

#endif // CompassSoft_Chip__include
