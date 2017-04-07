#ifndef HistoInSpillLED___include
#define HistoInSpillLED___include

class TH1D;
class TH2D;
class TProfile;
class TDirectory;
class TCanvas;

////////////////////////////////////////////////////////////////////////////////

class Calorimeter;

/*! \brief Set of histograms for reconstruction code testing.

*/
class HistoInSpillLED
{
  //============================================================================
  // Constructors and destructor
  //============================================================================

  public:

    /// Destructor
    virtual            ~HistoInSpillLED         (void) {}

    /// Constructor initialize all pointers to NULL
//                         HistoInSpillLED         (Calorimeter *calo); 
                        HistoInSpillLED         (void); 

  //============================================================================
  // Operators
  //============================================================================

  private:
  
//     /// You can not use assignment operator.
//     HistoInSpillLED&    operator=               (const HistoInSpillLED &);

//   public:
//     /// .Reset histograms
//     virtual  void        Reset(void);

  //============================================================================
  // Attributes
  //============================================================================

  public:

//     bool  led_histos_booked;
    TH2D                                     *h2fitpar_good4monitoring;
    TH2D                                     *h2fitpar_bad4monitoring;
    std::vector <TProfile *>               p1_InSpillLEDCell;

//     protected:
//     Calorimeter                        *calorimeter_;
};

#endif // HistoInSpillLED___include
