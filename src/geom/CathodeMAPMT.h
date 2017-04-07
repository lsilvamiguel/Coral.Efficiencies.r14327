
#ifndef CathodeMAPMT_h
#define CathodeMAPMT_h

#include "CsRICH1UpGrade.h"

class CathodeMAPMT : public CsRICH1UpGrade::CathodePlane 
{
 private:
  // to set values from options file
  bool use_quantum_efficiency_corrections_;
  bool use_obsolete_lens_numbering_;
  bool use_old_zebra_format_;
  bool rescale_wavelength_;
  double peak_eff_ ;
  void               InitCsOpt          ( void );
  
  
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
      h1_Time = NULL;
      h1_diffX = NULL;
      h1_diffY = NULL;
      h2_padXY = NULL;
      h2_pseudoXY = NULL;
    } 
  public:
    CsHist1D                   *h1_Time;
    CsHist1D                   *h1_diffX;
    CsHist1D                   *h1_diffY;
    CsHist2D                   *h2_padXY;
    CsHist2D                   *h2_pseudoXY;
  };
  
 public:
  CathodeMAPMT ( int id );
  void    MakeMCResponse( void );
  void    Clear         ( void );
  void    FillMCDecodingHisto( void );
  
 private:
  ///  Flags    
   CathodeMAPMT::Flags  flags_;
   CathodeMAPMT::Histo    histo_;
   
 private:
   std::vector< CsRICH1UpGrade::CathodePAD >            pseudopads_;
};

#endif // CathodeMAPMT_h
