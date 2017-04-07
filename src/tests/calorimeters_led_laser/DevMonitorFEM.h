#ifndef DevMonitorFEM_h
#define DevMonitorFEM_h

#include <string>

class CsCalorimeter;

// Basic functionality of  SourceLaserLED object is to keep map of cells flashed by this source.
// If we have only one source in detector( usual situation) we'll avoid
// creation of SourceLaserLED object  not to make complification of elementary case.
class SourceLaserLED 
{
  public:
    SourceLaserLED       ( CsCalorimeter *c, const std::string &name );
    void                         Initialize( void );
    void                         SetName               ( const std::string &name ) { name_= name;}
    const std::string       &GetName( void ) const { return name_;}
    bool                          IsMyCell(  int cell_id ) const;
  private:
    void                          InitCellsMap  ( void );
  // ==========================================
  //  Attributes, data
  // ==========================================
  private:
    std::string             name_;
    CsCalorimeter      *c_;
    std::map<int,bool>   map_cell2source_;
};  

class DevMonitorFEM 
{
  // ============================================================================
  // Types, constants
  // ============================================================================
    public:      
      const static int NFEMs_MAX = 20;
    // =========================================================================
    // Constructor
    // =========================================================================
      DevMonitorFEM( void );

    public:

      void          Init                  ( CsCalorimeter *c, const std::string &name = "" );
  //  We can think about some virtual Configure functions later
      void          ConfigureECAL1FEM ( void );
      void          ConfigureECAL0FEM ( void );

      void          InitReadOut           ( void );

      void          SetName               ( const std::string &name ) { name_= name;}
      const std::string &GetName          ( void ) const { return name_;}
      void          Clear                 ( void );

      int             GetFEMIndexForCellNormalization(  int cell_id ) const;

      std::pair< bool,double>        GetCorrectionFactor   ( int cell_id ) const;
      void          MonitorInSpill            ( void );
//   private:    
//       void          PutDummyInSpillValues ( void );
  public:    
      void          DecodeChipDigits    (const CS::Chip::Digits& digits);
      void          DecodeChipDigitsECAL0    (const CS::Chip::Digits& digits);
      void          DecodeChipDigitsECAL1    (const CS::Chip::Digits& digits);
      void          DecodeChipDigitsHCAL1    (const CS::Chip::Digits& digits);
      void          TurnedOutGood4Corrections ( void ) {  fem_is_good4correction_=true;}

  /// FEMs monitoring functions which uses FEMs defined for LED/LASER signals normalization also for 
  /// other FEMs cross normalization to see relative fluctuations
      void          CalculateInternalCrossNormalization( void );

      void          SetSignal             ( int fem_index, double signal);

      int           InputFEMInfo          ( size_t when, const std::string &s);
      int           OutputFEMInfo         ( std::string &s ) const;

//       int           InputFEMInfoInSpills( size_t when, const std::string &s);
      int           OutputFEMInfoInSpills( std::string &s ) const;
      
      std::string   GetSummaryComment     ( void ) const;

    private:
      void          InitMemory            ( double fem_ref_default );
      void          InitMapCell2FEM  ( void );

    public:

    // ==========================================
    //  Attributes, data
    // ==========================================
      std::string        name_;
//    OptionsFEM                      options_fem_;
      bool              fem_is_configured_;
      unsigned          nFEMs_;
      
      std::vector<int>    fem_norm_id_;

      double            fem_range_min_;
      double            fem_range_max_;

//      double            fem_ref_default_;
//
      CsCalorimeter      *c_;
      std::vector < SourceLaserLED * > sources_;
      std::map<int,int>             map_cell2fem_;

    private:
// Data flags ?? Not used ??
    bool fem_is_good4correction_;
    public:
// Data 
    std::vector <CsDigitizerSADC *> dig_max_fem_;
    std::vector <  std::pair<bool,double> > fem_signal_;
    std::vector <  int > map_fem2source_;
//    bool               fem_measured_;
    
    std::vector < Reco::StatInfo >  fem_norm_new_;
    std::vector < Reco::StatInfo >  fem_new_;
    std::vector < Reco::StatInfo >  fem_old_;
    std::vector < Reco::StatInfo >  fem_ref_;

// Data store
    std::vector < double >         store_fem_all_[NFEMs_MAX];
    std::vector < double >         store_femn_all_[NFEMs_MAX];
   std::vector < Reco::EventID > led_events_id_;
};  // end class DevMonitorFEM

////////////////////////////////////////////////////////////////////////////////

#endif // DevMonitorFEM_h
