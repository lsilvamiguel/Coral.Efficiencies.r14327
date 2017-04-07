#ifndef CompassSoft_ChipGandalf__include
#define CompassSoft_ChipGandalf__include

#include "Chip.h"

namespace CS {

using std::ostream;
using std::string;
using std::vector;
using std::cout;

class TriggerTimeConfig;

////////////////////////////////////////////////////////////////////////////////

/*! \brief This is ChipGandalf class.
    

    \author Boris Iven
*/
class ChipGandalf: public Chip
{
  //============================================================================
  // Types, constants
  //============================================================================

  private:

    enum {ErrorMarker=0xfff};

  public:
  
    /// Operation modes
  	enum OpMode
  	{
  		GADC_UNKNOWN			=		-1,
  		GADC_PROC,				//!Gandalf-ADC, normal mode
  		GADC_PROC_IL,			//!Gandalf-ADC, normal mode, interleave
  		GADC_FRAME,				//!Gandalf-ADC, frame mode
  		GADC_FRAME_IL,			//!Gandalf-ADC, frame mode, interleave
  		GADC_DEBUG,				//!Gandalf-ADC, debug mode
  		GADC_DEBUG_IL,			//!Gandalf-ADC, debug mode, interleave
  		GSCALER,				//!Gandalf Scaler
  		GTDC					//!Gandalf TDC
  		//TODO: all the other modes...?
  	};
  	


    /*! \brief DataID class for Gandalf
    	
    	\author Boris Iven
    */


    class DataID
    {
      public:
      	DataID(){ u.data_id = 0; }
        DataID(const Chip::DataID &d) {u.data_id=d;}
        DataID(uint16 srcID, uint16 geoIdorPort, uint16 chan) {u.s.src_id=srcID; u.s.port=0; u.s.geo_id=geoIdorPort; u.s.chan=chan; u.s.none=0;}
        //DataID(uint16 srcID, uint8 port, uint16 geo_id, uint16 chan) {u.s.src_id=srcID; u.s.port=port; u.s.geo_id=geo_id; u.s.chan=chan; u.s.none=0;}        
        operator Chip::DataID (void) const {return u.data_id;}
        union
        {
            Chip::DataID data_id;
            struct
            {
                Chip::DataID  src_id: 16, port: 8, geo_id: 16, chan: 16, none:8;
            } s;
        } u;
    };

    /*! \brief base class for all Gandalf Digit classes
    	
    	\author Boris Iven
    */
	

    class Digit: public Chip::Digit
    {
      public:
        virtual        ~Digit                   (void) {}

            			Digit       (const DataID &d,const DetID &id):Chip::Digit(d, id), m_mapMode(-1),m_channel(-1),
            							m_x(-1), m_y(-1), m_z(-1), m_r(-1), m_phi(-1), time_unit(-1), time_decoded(-1){};
						Digit		(const DataID &d,const DetID &id, int32 channel):Chip::Digit(d, id), m_mapMode(Map::MM_NORMAL), m_channel(channel), m_x(-1), m_y(-1), m_z(-1), \
										m_r(-1), m_phi(-1), time_unit(-1), time_decoded(-1){} //!CTOR for normal mapping
						Digit		(const DataID &d,const DetID &id, int32 channel, int32 x, int32 y, int32 z):Chip::Digit(d, id), m_mapMode(Map::MM_XYZ), m_channel(channel), \
										m_x(x), m_y(y), m_z(z), m_r(-1), m_phi(-1), time_unit(-1), time_decoded(-1){} //!CTOR for XYZ-mapping
						Digit		(const DataID &d,const DetID &id, int32 channel, int32 r, float phi, int32 z):Chip::Digit(d, id), m_mapMode(Map::MM_RPHIZ), m_channel(channel), \
										m_x(-1), m_y(-1), m_z(z), m_r(r), m_phi(phi), time_unit(-1), time_decoded(-1){} //!CTOR for R-Phi-Z-mapping

        virtual void    Print                   (ostream &o=cout,const string &prefix="") const;

        virtual const char*     GetNtupleFormat         (void) const {return "";}
        virtual vector<float>   GetNtupleData           (void) const{vector<float> v; return v;};
        
        int32			getMapMode				(void){ return m_mapMode; }
        
        void			setChannel				(int32 channel){ m_channel = channel; }
        int32			getChannel				(void) const { return m_channel; }

        void			setX					(int32 x){ m_x = x; }
        int32			getX					(void) const { return m_x; }
        void			setY					(int32 y){ m_y = y; }
        int32			getY					(void) const { return m_y; }
        void			setZ					(int32 z){ m_z = z; }
        int32			getZ					(void) const { return m_z; }
        
        void			setR					(int32 r){ m_r = r; }
        int32			getR					(void) const { return m_r; }
        void			setPhi					(float phi){ m_phi = phi; }
        float			getPhi					(void) const { return m_phi; }

		double          GetTimeUnit             (void) const {return time_unit;}

		/// The cost of one time bin
		void            SetTimeUnit             (double t) {time_unit=t;}

		double          GetTimeDecoded          (void) const {return time_decoded;}
		void            SetTimeDecoded          (double t) {time_decoded=t;}

  
      private:
      
      	//mapping attributes
      	int32			m_mapMode;
      	int32     		m_channel;				//the channel AFTER the mapping (wire)!
		
        int32			m_x;
        int32			m_y;
        int32			m_z;
        
        int32			m_r;
        float			m_phi;  
        
		double          time_unit;              ///< unit of time measurement
		double          time_decoded;           ///< time in ns with the respect to trigger time

    };
    
    
    
    /*! \brief Digit class for thr Gandalf ADC mode
    	
    	\author Boris Iven
    */
    
    class DigitGADC: public Digit
    {
      public:
      
      					DigitGADC	(const DataID &d,const DetID &id, int32 opMode):Digit(d, id),
      							m_opMode(opMode),m_pulseDet(0),m_baseLErr(-1),m_baseLine(-1),m_coarseTimeLSB(-1),m_coarseTimeMSB(-1),
      							m_frameTime(-1),m_hiResTime(-1),m_integral(-1),m_maxAmplitude(-1),m_threshold(-1),m_fracfact(-1),m_delay(-1) {}
      					DigitGADC	(const DataID &d,const DetID &id, int32 opMode, int32 channel):Digit(d, id, channel),
      							m_opMode(opMode),m_pulseDet(0),m_baseLErr(-1),m_baseLine(-1),m_coarseTimeLSB(-1),m_coarseTimeMSB(-1),
      							m_frameTime(-1),m_hiResTime(-1),m_integral(-1),m_maxAmplitude(-1),m_threshold(-1),m_fracfact(-1),m_delay(-1) {}
      					DigitGADC	(const DataID &d,const DetID &id, int32 opMode, int32 channel, int32 x, int32 y, int32 z):Digit(d, id, channel, x, y, z),
      							m_opMode(opMode),m_pulseDet(0),m_baseLErr(-1),m_baseLine(-1),m_coarseTimeLSB(-1),m_coarseTimeMSB(-1),
      							m_frameTime(-1),m_hiResTime(-1),m_integral(-1),m_maxAmplitude(-1),m_threshold(-1),m_fracfact(-1),m_delay(-1) {}
      					DigitGADC	(const DataID &d,const DetID &id, int32 opMode, int32 channel, int32 r, float phi, int32 z):Digit(d, id, channel, r, phi, z),
      							m_opMode(opMode),m_pulseDet(0),m_baseLErr(-1),m_baseLine(-1),m_coarseTimeLSB(-1),m_coarseTimeMSB(-1),
      							m_frameTime(-1),m_hiResTime(-1),m_integral(-1),m_maxAmplitude(-1),m_threshold(-1),m_fracfact(-1),m_delay(-1) {}
      					
        virtual void    Print                   (ostream &o=cout,const string &prefix="") const;

        
        virtual const char*     GetNtupleFormat (void) const {return "chan:amp:int:time_dec:msb:lsb:hr:time_un(:frame_time(:N_samples:Sample_1:..:Sample_N))";}
        
        virtual vector<float>   GetNtupleData   (void) const;
        
        int32			getOpMode				(void) const{ return m_opMode; }
        void			setOpMode				(int32 opMode) { m_opMode = opMode; }
        
        int32			getPulseDet				(void) const{ return m_pulseDet; }
        void			setPulseDet				(int32 pulseDet) { m_pulseDet = pulseDet; }
        
        int32			getBaseLErr				(void) const{ return m_baseLErr; }
        void			setBaseLErr				(int32 baseLErr) { m_baseLErr = baseLErr; }
        
        int32			getBaseLine				(void) const{ return m_baseLine; }
        void			setBaseLine				(int32 baseLine) { m_baseLine = baseLine; }

        void			addFrameData			(uint16 data){ m_frameData.push_back( data ); }
        int32			getNumSamples			(void) const{ return m_frameData.size(); }
        uint16			getSample				(int32 idx)const{ return m_frameData[idx]; }

		void			setMaxAmplitude		(uint32 maxAmp){ m_maxAmplitude = maxAmp; }
		uint32			getMaxAmplitude		(void) const { return m_maxAmplitude; }
		void			setIntegral				(uint32 integral){ m_integral = integral; }
		uint32			getIntegral				(void) const { return m_integral; }
		void			setFrameTime			(uint32 ft){ m_frameTime = ft; }
		uint32			getFrameTime			(void) const { return m_frameTime; }
		void			setCoarseTimeLSB	(uint32 ctLSB){ m_coarseTimeLSB = ctLSB; }
		uint32			getCoarseTimeLSB	(void) const { return m_coarseTimeLSB; }
		void			setCoarseTimeMSB	(uint32 ctMSB){ m_coarseTimeMSB = ctMSB; }
		uint32			getCoarseTimeMSB	(void) const { return m_coarseTimeMSB; }
		void			setHiResTime			(uint32 hrt){ m_hiResTime = hrt; }
		uint32			getHiResTime			(void) const { return m_hiResTime; }
		void			setFrameData			(vector<uint16> frameData){ m_frameData = frameData; }
		vector<uint16>	getFrameData		(void) const {return m_frameData; }
		void			setThreshold		(uint32 threshold){ m_threshold = threshold; }
		uint32			getThreshold		(void) const { return m_threshold; }
		void			setFracfact 		(uint32 ff){ m_fracfact = ff; }
		uint32			getFracfact			(void) const { return m_fracfact; }
		void			setDelay 			(uint32 delay){ m_delay = delay; }
		uint32			getDelay			(void) const { return m_delay; }
		
		double 			getTime (void) const {
			return ( ( ((uint64)getCoarseTimeMSB())<<21 ) +getCoarseTimeLSB()+(getHiResTime()/(double)1024));}


      private:      	

        
        int32			m_opMode;				//!Operation mode (normal, frame, debug, interleave..., see ChipGandalf::OpMode)
        int32			m_pulseDet;
        int32			m_baseLErr;
        int32			m_baseLine;
        uint32          m_coarseTimeLSB;
        uint32          m_coarseTimeMSB;
        uint32          m_frameTime;
        uint32          m_hiResTime;
        uint32          m_integral;
        uint32          m_maxAmplitude;
        uint32			m_threshold;
        uint32			m_fracfact;
        uint32			m_delay;
        vector<uint16>	m_frameData;

    };

    
	//Structs used for data decoding
	//Gandalf ADC:
	
	
	//Gandalf ADC Header/Trailer
	struct GADCHeadTrail
	{
		uint32		RDM				:		 4; //!Reaout mode
		uint32		WindowLength	:		11; //!number of frame words
		uint32		SysMon			:		 5; //!System monitor
		uint32		ChID			:		 4; //!Channel ID
		uint32		EventNum		:		 6; //!Event number
		uint32		Head_Trail		:		 2; //!00 for head;01 for trail

		bool operator== (const GADCHeadTrail &a) const
		{
		    return a.RDM==this->RDM && a.WindowLength==this->WindowLength && a.SysMon==this->SysMon && a.ChID == this->ChID && a.EventNum == this->EventNum;
		}
		bool operator!= (const GADCHeadTrail &a) const
		{
			return !( *this==a );
		}
	};

	
	//Gandalf ADC frame word
	struct GADCFrameWord
	{
		uint32		DataLow			:		14; //!Data word bit 0-13
		uint32		Spacing			:		 2;
		uint32		DataHigh		:		14; //!Data word bit 16-29
		uint32		Zero			:		 1; //!This bit must be zero
		uint32		DataFlag		:		 1; //!This bit must be one
	};
	
	//Gandalf ADC word
	struct	GADCProcWord
	{
		uint32		Data			:		31; //!data
		uint32		DataFlag	:		 1; //!This bit must be one
	};

//	//Gandalf ADC word for frametime ctimeMSB data debug mode(FrameTime, CTime)
//	struct	GADCNormWordFTimeCTime
//	{
//		uint32		CTdataMSB	:		17; //!Integral data
//		uint32		FTime			:		12; //!Type identifier
//		uint32		DataID		:		 2; //!Type identifier
//		uint32		DataFlag	:		 1; //!This bit must be one
//	};

	//Gandalf ADC word for normal data normal mode(Channel,Baseline,Integral)
	struct	GADCProcChIntWord
	{
		uint32		Integral		:		16; //!Integral data
		uint32		Baseline		:   11; //!Baseline data
		uint32		ChID				:		4; //!channel number
		uint32		DataFlag	:		 1; //!This bit must be one
	};
	
		//Gandalf ADC word for normal data normal mode(CoarseTimeMSB, Amplitude)
	struct	GADCProcCTimeAmpWord
	{
		uint32		Amplitude		:		14; //!max Amplitude
		uint32		CTdataMSB		:   17; //!CoarseTime MSB
		uint32		DataFlag	:		 1; //!This bit must be one
	};
	
		//Gandalf ADC word for normal data normal mode(Channel,Baseline,Integral)
	struct	GADCProcCTimeHResWord
	{
		uint32		HResTime		:		10; //!HighRes data
		uint32		CTdataLSB		:   21; //!CoarseTime LSB
		uint32		DataFlag	:		 1; //!This bit must be one
	};	
	//Gandalf ADC word for dbg mode; contains constants used to calculate CFD and the FrameTime
	struct	GADCDbgFTimeWord
	{
		uint32		FrameTime		:	12; //!FrameTime data
		uint32		Delay			:   5; //!Delay for cfd
		uint32		FracFact		:	6; //!Fraction Factor for cfd
		uint32		Threshold		:	8; //!Threshold for cfd
		uint32		DataFlag		:	1; //!This bit must be one
	};

    class Map : public Chip::Map
    {
      public:
      
        Map                     (const ObjectXML &o);
                        
                        
	  	//Mapping modes
	  	enum MapMode
	  	{
	  		MM_NORMAL,
	  		MM_XYZ,
	  		MM_RPHIZ
	  	};
                        

      public:

        /// Print the map definition
        virtual void            Print                   (ostream &o=cout,const string &prefix="") const;

        virtual void            AddToMaps               (Maps &maps,DaqOption &option) const;
        
        int32					getMapMode				(void){ return m_mapMode; }
        
        /*! @brief Read trigger time configuration. */
        void            		ReadTTConfig            (TriggerTimeConfig &tt_conf) const;

      private:
    
      	///Mapping mode
  		int32					m_mapMode;

      	
      	int16           m_port;
      	int32			m_geoID;
      	
      	int32			m_x;
      	int32			m_y;
      	int32			m_z;
      	int32			m_r;      	      	      	
      	int32			m_phi;      	
      
        double    time_unit;              ///< unit of time measurement
	

    };

  //============================================================================
  // Constructors and destructor
  //============================================================================

  public:
    
    /// Destructor
                       ~ChipGandalf                  (void) {}
    
    /*! \brief Base constructor.
        \sa DaqEvent::DaqEvent(void const * const buf,bool copy_buf)
    */
                        ChipGandalf                  (void const * const buf,bool copy_buf,DaqOption &opt,DaqEvent &ev);

  private:

     /// Copy constructor
                        ChipGandalf                  (const ChipGandalf &e);
  
  //============================================================================
  // Operators
  //============================================================================

  private:

    /// Assignment operator
    ChipGandalf             &operator =              (const ChipGandalf &e);

  private:
  
  	/// Operation mode
  	int32					m_opMode;
  	std::list<DigitGADC*>	um_digits; //!list of unmapped digits for Gandalf ADC mode (to be filled in 'scan' methods)
  	
  

  //============================================================================
  // Methods
  //============================================================================

  public:
  
    /// Print equipment
    virtual void                Print                   (ostream &o=cout,const string &prefix="") const;

    /// Clear the chip status (do \b not modify the real data buffer).
    virtual void                Clear                   (void) {}

    /// \return Chip name. This is "ChipGandalf"
    virtual string              GetName                 (void) const {return "ChipGandalf";}

    /// \return \b true if data is TDC readout.
    bool                isUndefined              (void) const {return (GetSLink().GetFormat()==0);}

    /// \return \b true if data is TDC readout.
    bool                isTDCReadout             (void) const {return (GetSLink().GetFormat()&(1));}
    
    /// \return \b true if TDC data is in \c debug mode, otherwise the data are in sparsified mode.
    bool                isTDCDebugMode             (void) const {return !(GetSLink().GetFormat()&(1<<1));}

    /// \return \b true if TDC data is in \c latch mode.
    bool                isTDCLatchMode             (void) const {return GetSLink().GetFormat()&(1<<3);}

    /// \return \b true if TDC data is in high resolution mode, otherwise data is in normal resolution.
    bool                isTDCHighResolutionMode    (void) const {return GetSLink().GetFormat()&(1<<4);}

    /// \return \b true if TDC data is with leading AND trailing edge, otherwise data is with leading OR trailing edge.
    bool                isTDCLeadTrail             (void) const {return GetSLink().GetFormat()&(1<<5);}

    /// \return \b true if TDC data is in HOTL-CMC mode, otherwise data is in TDC-CMC/GANDALF-TDC mode.
    bool                isHOTLCMC               (void) const {return GetSLink().GetFormat()&(1<<5);}

    /// \return \b true if TDC data is in TDC-CMC/GANDALF-TDC mode.
    bool                isGandalfTDC             (void) const {return !isHOTLCMC();}

    /// \return \b true if this the the first event of a run    
    bool				isFEOR				(void) const {return ((GetSLink().GetEventType() & 0x1F)) == 0x1C;}

    /// \return \b true if this the the first event of a spill    
    bool				isFEOS				(void) const {return ((GetSLink().GetEventType() & 0x1F)) == 0x1E;}

    
    /// \return \b true if data is ADC readout.
    bool                isADCReadout             (void) const {return !(isTDCReadout());}
    
    /// \return \b true if data is in GANDALF-ADC mode.    
    bool				isGandalfADC			(void) const {return ((GetSLink().GetFormat() & 0x7F) >> 4) == 1;}
    
    /// \return \b true if data is in GANDALF-Scaler mode.    
    bool				isGandalfScaler			(void) const {return ((GetSLink().GetFormat() & 0x7F) >> 4) == 2;}
    
    /// \return \b true if ADC data is in normal mode.    
    bool				isADCProcMode			(void) const {return ((GetSLink().GetFormat() & 0x07) >> 1) == 0;}
        
    /// \return \b true if ADC data is in frame mode.    
    bool				isADCFrameMode			(void) const {return ((GetSLink().GetFormat() & 0x07) >> 1) == 1;}

    /// \return \b true if ADC data is in debug mode.    
    bool				isADCDbgMode			(void) const {return ((GetSLink().GetFormat() & 0x07) >> 1) == 2;}

    /// \return \b true if ADC data (normal, frame, debug) is interleaved mode.        
    bool				isADCIlm				(void) const {return ((GetSLink().GetFormat() & 0x0F) >> 3) == 1;}
    
    /// Decode data and \b add new digits to \b digits_list.
    virtual void                Decode                  (const Maps &maps,Digits &digits_list,DaqOption &opt);

    /// \return the time difference between the time point and trigger time, taking into account the overolling
    static double       TimeDifference          (int time,double trigger_time);

  protected:

    virtual void                Scan                    (DaqOption &opt);
    void						scanADC 				(DaqOption &opt, bool interleaved);

  private:


};

////////////////////////////////////////////////////////////////////////////////

inline ostream &operator << (ostream &o,const ChipGandalf &e)
{
  e.Print(o);
  return o;
}

////////////////////////////////////////////////////////////////////////////////

} // namespace CS



#endif // CompassSoft_ChipGandalf__include
