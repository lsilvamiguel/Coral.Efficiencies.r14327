/*!
   \file    
   \brief   
   \author  Kolosov Vladimir            Vladimir.Kolosov@cern.ch  
   \version 
   \date    
*/

#ifndef CathodeAPV_h
#define CathodeAPV_h

#include "CsRICH1UpGrade.h"

class ElectronicAPV
{
    public:
  class Shaper
  {
    public:
      Shaper ( const std::vector < double > &shape, double ballistic_deficit, 
                                                             double time_tik,unsigned zero_time_ch);
      double  GetSignal( double amp_signal, double time_signal, double time_measuring ) const;
    private:
      std::vector < double > shape_;
      double ballistic_deficit_;
      double norm_shape_max_;
      double time_tik_;
      unsigned zero_time_ch_;   
  };


  class DataMCAPV
  {
    public:
       DataMCAPV        ( void ):id_(0){data_.clear();}
       DataMCAPV        ( int id, const std::vector <int> d ) :id_(id),data_(d){}
    /// Destructor
       virtual         ~DataMCAPV                    (void) {}
    /// Copy constructor
                         DataMCAPV                   (const DataMCAPV &d) {*this=d;}
    public:
    /// Assignment operator
       DataMCAPV         &operator =              (const DataMCAPV &d) {  if( &d!=this ) {
                                                                              id_  =d.id_;
                                                                              data_=d.data_;
                                                                            }                                   
                                                                            return *this;
                                                                         } 
    /// Data memebers
    public:
    /// Channel Id
       int id_;
    /// data
       std::vector<int> data_;
  };

  class ChannelData
  {
    public:
      ChannelData                ( void ){} 
      void               Set    ( double amp_unit, double noise, double signal_cut, double  ped ); 
    public:
       double amp_unit_;
       double noise_;
       double signal_cut_in_sigmas_;
       double pedestal_;

       int int_pedestal_;
       int signal_cut_in_units_;
  };

  public:
    ElectronicAPV        ( void );
   
  /// Cleaning
    void               Clear              ( void );
    void               AddSignal          ( double amp, double time );
    std::vector <int>  Read               ( const ChannelData &ch );
    std::vector <double>  ReadTestSignal  ( void );
    void               PrintSettings      ( void );
    static void       Encode             ( int a0, int a1, int a2, std::vector <int> &out );
    static double     ShaperSignal0      ( double tns );
    static void        TestSignal         ( void );

  private:
    void               InitCsOpt          ( void );

  public:
    double noise_;   
    double pedestal_;   
    double amp_unit_;   
    double signal_cut_in_sigmas_;   
    double noise_scale_factor_;   
  private:
    const static unsigned SIZE_=2000;
    const static  double ballistic_deficit_;

//     double noise_;   
//     double pedestal_;   
//     double amp_unit_;   
    double time_tik_;
     
    unsigned zero_time_ch_;   
    unsigned read_ch0_;   
    unsigned read_ch1_;   
    unsigned read_ch2_;   
    double data_[SIZE_];
    
    Shaper  *shaper_;  
//    Options  
    bool add_noise_;
    bool calculate_signal_only_in_measuring_poins_; 

};


class CathodeAPV : public CsRICH1UpGrade::CathodePlane 
{
  public:
    class Avalanche 
    {
      public:
       Avalanche ( double x, double y, double e, double time, double multiplication );
      public:
     ///  Avalanche position and amplitude
       double     x_;  
       double     y_;  
       double     amp_;  
       double     time_;  
    };

    class Signal 
    {
      public:
       Signal ( double amp, double time, const CsMCRICH1Hit *hit ) : amp_(amp),time_(time),hit_(hit){}
      public:
       double             amp_;  
       double            time_;
       const CsMCRICH1Hit *hit_;  
    };
  private:

  class Flags
  {
    public:
    Flags ( void ) {
                      histo_mc_decoding_booked_ = false;
                    } 
    bool            histo_mc_decoding_booked_;
  };

  class Histo
  { 
    public:
                  Histo (void) {
                                  h1_A2 = NULL;
                                  h1_rA1A2 = NULL;
                                  h1_rA1A2cut = NULL;
                                  h2_A2A1 = NULL;
                                  h2_A2A0 = NULL;
                                  h2_A2A1_A0corr = NULL;
                                  h2_Chisq_Afit = NULL;
                                } 
    public:
     CsHist1D                   *h1_A2;
     CsHist1D                   *h1_rA1A2;
     CsHist1D                   *h1_rA1A2cut;
     CsHist2D                   *h2_A2A1;
     CsHist2D                   *h2_A2A0;
     CsHist2D                   *h2_A2A1_A0corr;
     CsHist2D                   *h2_Chisq_Afit;
  };

  public:
    CathodeAPV ( int id );

    void    MakeMCResponse( void );
    void    Clear( void );
    void    FillMCDecodingHisto( void );

    int     GetRDreadoutDirX( void ) const { return  rd_read_out_dirx_;}                                                
    int     GetRDreadoutDirY( void ) const { return  rd_read_out_diry_;}                                                

  protected:
    void    InitCsOpt( void );


  private:
    void    SetRDreadoutDirs( int dirx, int diry ) { if(dirx < 0) rd_read_out_dirx_=-1;
                                                      if(diry < 0) rd_read_out_diry_=-1;
                                                    }
     void         SetMCReadoutParmeters( void );
     bool         MakeAvalanche ( const CsMCRICH1Hit *hit );
     void         MakeElectronicMCDigits( void );
     double CathodeQuantumEfficiency ( double e, double x, double y );
     static double CathodeQuantumEfficiencyTable ( double e );
     static double PadResponse( double dxpad, double dypad, int ixs, int iys );
    
  private:
     ///  Flags    
     CathodeAPV::Flags  flags_;
     /// Options    
     bool  use_quantum_efficiency_corrections_;
  
     double peak_eff_ ;
     /// Factor for quantum efficiency
     double electron_drift_velocity_;
     double multiplication_gain_;
     /// Internal Chamber geometry     
     unsigned nwires_;
     double wire0_position_;
     double wire_pitch_size_;
     
     int rd_read_out_dirx_; 
     int rd_read_out_diry_; 

     ElectronicAPV apv_;
     ElectronicAPV::ChannelData  ch_apv_[72][72];

     std::vector <CathodeAPV::Avalanche> mc_avalanches_; 
     std::vector< CathodeAPV::Signal >    signals_[72][72];
     std::list <ElectronicAPV::DataMCAPV> data_apv_;
     
     CathodeAPV::Histo    histo_;
};
#endif // CathodeAPV_h

