/*!
   \file    
   \brief   
   \author  Kolosov Vladimir            Vladimir.Kolosov@cern.ch  
   \version 
   \date    
*/

#ifndef CsRICH1UpGrade_h
#define CsRICH1UpGrade_h

#include "coral_config.h"
#include <list>
#include <vector>
#include <string>

#include "DaqDataDecoding/Chip.h"

#include "CsRICH1Detector.h"
class CsMCHit;
class CsDigit;
class CsMCRICH1Hit;
#include "CsHistograms.h"
#include "CsHist.h"

class RICH1Analyser;

class CsRICH1UpGrade : public CsRICH1Detector
{
 public:
  /// Some usefull structures
  class Options
  {
    public:
    Options ( void ) {
                         printinfo_ = false; 
                         printdebug_ = false;
                         fill_histo_  = false;
                         clean_digits_ = false;
                         ignore_mc_hits_ = false;
                         mchits_time_min_=-2000.;
                         mchits_time_max_= 2000.;
                         invert_mc_time_gate_=false;
                     } 
    bool               printinfo_;
    bool               printdebug_;
    bool               fill_histo_;
    bool               clean_digits_;
    bool               ignore_mc_hits_;
    double             mchits_time_min_;
    double             mchits_time_max_;
    bool               invert_mc_time_gate_;
  };

  class Flags
  {
    public:
    Flags ( void ) {
                         histo_booked_ = false;
                         initialization_done_ = false;
                     } 
    bool               histo_booked_;
    bool               initialization_done_;
  };
  
  class CathodePAD 
  {
    public:
     CathodePAD ( void ) : cathode_id_(-1),ix_(-1),iy_(-1),amp_(0.),time_(0.) {}
     CathodePAD ( int ic, int ix, int iy, double e, double time ) : cathode_id_(ic),
                                             ix_(ix),iy_(iy),amp_(e),time_(time) {}
    /// Copy constructor
                          CathodePAD                 (const CathodePAD &d) {*this=d;}
    public:
    /// Assignment operator
       CathodePAD         &operator =              (const CathodePAD &d) {  if( &d!=this ) {
                                                                               cathode_id_ = d.cathode_id_;
                                                                               ix_ = d.ix_;  
                                                                               iy_ = d.iy_;  
                                                                               amp_ = d.amp_;  
                                                                               time_ = d.time_;  
                                                                               el_digits_ = d.el_digits_;  
                                                                            }                                   
                                                                            return *this;
                                                                         } 
    public:
   ///  To identify RICH1 cahode pad we need 3 integers 
     int     cathode_id_;
     int     ix_;  
     int     iy_;  
   ///  pad amplitude 
     double  amp_;  
   ///  time 
     double  time_;  
   ///  electronic digits 
     std::vector <int> el_digits_;  
  };

  class Cluster 
  {
    public:
     Cluster ( double x, double y, double amp, double time ) : x_(x),y_(y),amp_(amp),time_(time){}
    public:
     double  x_;
     double  y_;
     double  amp_;
     double  time_;
  };

  class CathodePlane  
  {
    public:
     enum                ReadOutType  { IDEAL=0, GASSIPLEX=1, APV=2, MAPMT=3 };
    
     CathodePlane ( int id );
     CathodePlane ( int id, double center_position[], double dir_normal[] );
      /// Destructor
     virtual            ~CathodePlane   (void) {}

     virtual std::list < CsMCRICH1Hit* > &GetPadMCHits( int ixpad, int iypad );

     void                  SetReadOutType   ( ReadOutType t ) {read_out_type_=t;}
     ReadOutType           GetReadOutType   ( void )  const { return read_out_type_;}
     virtual void         InitFromCoral    ( void );                             
     virtual void         Clear            ( void );
     virtual void         MakeMCResponse   ( void );
     virtual void         SortMCHits       ( void );
     virtual void         MakeClusters     ( void );
     virtual void         FillMCDecodingHisto ( void ) {}

     unsigned     NMCHits    ( void ) const { return mchits_.size(); }
     unsigned     NPadsFired    ( void ) const { return pads_.size(); }
     const std::vector< CathodePAD > &GetPadsFired    ( void ) const { return pads_; }
     /// Temporary??
     std::vector< CathodePAD > &GetPads    ( void ) { return pads_; }

     void         AddMCHit    ( CsMCRICH1Hit* hit );
     
     std::pair < int, int>        GetPadIndexes( double x, double y );
     std::pair < double, double>  GetPadPosition( int ixpad, int iypad );

     static int DecodeIndxCathode( int iadr );
     static int DecodeIndxPadX( int iadr );
     static int DecodeIndxPadY( int iadr );

     /// Set efficiency reduction factor   
     void                 SetReductionFactor( double fact ) {efficiency_reduction_factor_ = fact;}

     ///  \return efficiency reduction factor   
     const double       &GetReductionFactor( void ) const {return efficiency_reduction_factor_;}

     void                 SetTBname        ( const std::string &name ) { tbname_ = name;}

     const std::string   &GetTBname        ( void ) const { return tbname_;}

      protected:
     virtual void        InitCsOpt         ( void );

     virtual int          SetPadAdr        ( int ixpad, int iypad ) const;

    public:

     /// Cathode ID    
     int    id_;
     double x0_image_;
     double y0_image_;

      private:

     /// PD  TBname  
      std::string              tbname_;

      protected:
      
     /// ReadOut Type    
     ReadOutType                   read_out_type_;
     /// MC hits    
     std::vector < CsMCRICH1Hit* > mchits_; 
     /// pads    
     std::vector< CathodePAD >            pads_;
     /// clusters    
     std::vector< CsRICH1UpGrade::Cluster >  padclusters_;
     /// MC hits in pads    
     std::list< CsMCRICH1Hit* >  padmc_hits_[6000];
     /// Efficiency reduction factor    
     double efficiency_reduction_factor_; 
     /// Geometry    
     double size_x_; 
     double size_y_; 
     /// Position in MRS    
     double center_position_[3];
     double dir_x_[3];
     double dir_y_[3];
     double dir_normal_[3];

     ///  Old pads structure geometry    
     unsigned nxpads_;
     unsigned nypads_;
     unsigned npads_;

     double pad_size_;
     double pad0_x_position_;
     double pad0_y_position_;

  };

 private:

  class Histo
  { 
    public:
                  Histo (void) {
                                  h1_MC_Chambers = NULL;
                                  h1_MC_Time = NULL;
                                  h1_MC_PreciseTime = NULL;
                                  h1_MC_PhotonE = NULL;
                                  h2_MC_ChambersXY = NULL;
                                  h2_MC_ChambersXYProblems = NULL;
                                  h1_RD_Amp = NULL;
                                  h2_RD_ChambersXY = NULL;
                                  h1_Digit_Amp = NULL;
                                  h2_Digit_ChambersXY = NULL;
                                  for( int ich=0; ich< 16; ich++ )
                                  {
                                    h2_MC_ChamberXY[ich] = NULL;
                                    h1_MC_ChamberTime[ich] = NULL;
                                    h1_MC_ChamberPreciseTime[ich] = NULL;
                                  } 
                                } 
    public:
     CsHist1D                               *h1_MC_Chambers;
     CsHist1D                               *h1_MC_Time;
     CsHist1D                               *h1_MC_PreciseTime;
     CsHist1D                               *h1_MC_PhotonE;
     CsHist2D                               *h2_MC_ChambersXY;
     CsHist2D                               *h2_MC_ChambersXYProblems;
     CsHist2D                               *h2_MC_ChamberXY[16];
     CsHist1D                               *h1_MC_ChamberTime[16];
     CsHist1D                               *h1_MC_ChamberPreciseTime[16];
     CsHist1D                               *h1_RD_Amp;
     CsHist2D                               *h2_RD_ChambersXY;
     CsHist1D                               *h1_Digit_Amp;
     CsHist2D                               *h2_Digit_ChambersXY;
  };

 public:

  CsRICH1UpGrade ( const int id, const std::string &TBname );

 /// Really public
 public:
  virtual void getMCDigits( std::list<CsMCHit*>& );


 public:

  /// Decode raw data
  void       DecodeChipDigits        (const CS::Chip::Digits &digits);
  void       DecodeChipDigitPD       (int icath, const CS::Chip::Digit &digit);

  void    MakeMCResponse                     ( std::list<CsMCHit*>& richHits );
//   const std::vector< CathodePAD >   &GetPads ( void ) const { return pads_;}

  unsigned                   NMCHits        ( void ) const {return mchits_.size();}
  std::vector< CsMCRICH1Hit* >       &GetMCHits      ( void ) {return mchits_;}
  unsigned     NCathodes                    ( void ) const {return cathodes_.size();}
  /// Development function
  void          Analysis                     ( void );

 private:
  void    InitCsOpt ( void );

  void    Clear ( void );                             

  void    InitFromCoral         ( void );                             

  void    StoreMCHits                       ( std::list<CsMCHit*>& richHits, double tmin, double tmax );

  void    MakeMCResponse                     ( void );

  void    FillHisto                          ( void );
  void    BookHisto                          ( void );
  void    BookMCHisto                        ( void );
  void    BookRDHisto                        ( void );
  void    FillMCDecodingHisto                ( void );
  void    FillDecodingHisto                  ( void );
  void    FillMCAnalysisHisto                ( void );
  void    Reconstruction                     ( void );

 public:
  /// Check that RICH1 initialization in Coral is done
  static  bool CoralIsReady                  ( void );
  /// Utilities temporary
  static std::pair<double,double> RanGau    ( void );
  static double                    Randm     ( void );
  static double                    RanExp    ( double b );


 private:

 /// MC hits
  std::vector< CsMCRICH1Hit* >           mchits_;
  
 /// Geometry
  std::vector<CathodePlane*>             cathodes_;

 /// Options
  Options                                options_;
  
 /// Flags
  Flags                                  flags_;

 /// Histo
  Histo                                  histo_;

 /// Analyser object
//  RICH1Analyser*                         analyser_;
};
#endif // CsRICH1UpGrade_h

