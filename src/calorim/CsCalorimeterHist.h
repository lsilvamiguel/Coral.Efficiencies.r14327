////////////////////////////////////////////////////////////////////////////////

/*! \brief CsCalorimeterHist class in COMPASS experiment
   Set of test histograms for COMPASS calorimeters.
   Authors:
     Vladimir  Kolosov   ( Vladimir.Kolosov@cern.ch,Kolosov@mx.ihep.su )
     Alexander Zvyagin   ( Alexander.Zviagine@cern.ch, Zvyagin@mx.ihep.su )
     Denis     Murashev  ( Denis.Mourachev@cern.ch, Murashev@sirius.ihep.su )


*/
//class CsCalorimeterHist;
#ifndef CsCalorimeterHist___include
#define CSCalorimeterHist___include

#include "coral_config.h"

class TH1D;
class TH2D;
class TDirectory;

////////////////////////////////////////////////////////////////////////////////

/*! \brief Set of histograms .

*/
class CsCalorimeterHist
{
  //============================================================================
  // Constructors and destructor
  //============================================================================

  public:

    /// Empty default constructor
                        CsCalorimeterHist         (void) {}

  private:

    /// You can not use copy constructor.
                        CsCalorimeterHist         (const CsCalorimeterHist &);

  //============================================================================
  // Operators
  //============================================================================

  private:

    /// You can not use assignment operator.
    CsCalorimeterHist&    operator=               (const CsCalorimeterHist &);

  //============================================================================
  // Attributes
  //============================================================================

  public:
    const static   int  trig_groups_max = 100;

    /// ROOT directory corresponed to the calorimeter
    TDirectory                         *root_dir;
    TDirectory                         *all_groups_dir;
    TDirectory                         *all_gcorr_dir;
    TDirectory                         *timing_dir;
    TDirectory                         *geom_test_dir;
    TDirectory                         *cells_dir;

    TH1D                               *h1_HitedGroups;
    TH2D                               *h2_GroupVsGroup;
    TH2D                               *h2_GroupVsCellX;
    TH2D                               *h2_GroupVsCellY;
    TH2D                               *h2_T1GroupVsEGroup;
    TH2D                               *h2_T2GroupVsEGroup;
    TH2D                               *h2_TGroupVsNCellE;
    TH1D                               *h1_NHits;
    TH1D                               *h1_Egroups[trig_groups_max];
    TH2D                               *h2_EtgEcm[trig_groups_max];
    TH2D                               *h2_EtgEcmT1[trig_groups_max];
    TH2D                               *h2_EtgEcmT2[trig_groups_max];

    TH2D                               *h2_Geom[trig_groups_max];
// Search for bad mapping
    TH2D                               *h2_SearchADC[trig_groups_max];
    TH2D                               *h2_SearchT1[trig_groups_max];
    TH2D                               *h2_SearchT2[trig_groups_max];

    TH2D                               *h2_Timing1[trig_groups_max];
    TH2D                               *h2_Timing2[trig_groups_max];

    TH1D                               *h1_ECell[5000];
};

class TestHistoSADC
{
  public:
   TestHistoSADC   ( void ) {
                               h2_TimevsTime=NULL;
                               h2_AmpvsAmp=NULL;
                               h2_DeltaTimevsTime=NULL;
                               h2_DeltaAmpvsAmp=NULL;

                               h2_TimevsTimeIsNoise0=NULL;
                               h2_AmpvsAmpIsNoise0=NULL;
                               h2_DeltaTimevsTimeIsNoise0=NULL;
                               h2_DeltaAmpvsAmpIsNoise0=NULL;

                               h2_TimevsTimeIsNoise1=NULL;
                               h2_AmpvsAmpIsNoise1=NULL;
                               h2_DeltaTimevsTimeIsNoise1=NULL;
                               h2_DeltaAmpvsAmpIsNoise1=NULL;

                             }

  std::vector<TH1D*> h1_NSamplesSADC;
  std::vector<TH1D*> h1_SizeSamplesSADC;
  std::vector<TH2D*> h2_SummMax;
  std::vector<TH2D*> h2_ChMaxAmp;
  std::vector<TH2D*> h2_SignalMax;
  std::vector<TH1D*> h1_Ped;
  std::vector<TH1D*> h1_Summ;
  std::vector<TH1D*> h1_Signal;
  std::vector<TProfile*> p1_Shape10SADC;
  std::vector<TProfile*> p1_Shape100SADC;
  std::vector<TProfile*> p1_Shape500SADC;
  std::vector<TProfile*> p1_Shape1000SADC;
  std::vector<TProfile*> p1_ShapeNM10SADC;
  std::vector<TProfile*> p1_ShapeNM100SADC;
  std::vector<TProfile*> p1_ShapeNM500SADC;
  std::vector<TProfile*> p1_ShapeNM1000SADC;
  std::vector<TH1D*> h1_TimeSADC;
  std::vector<TH2D*> h2_TimeSADCvsE;
  std::vector<TH2D*> h2_TimeSADCvsTCS;

  TH2D*              h2_TimevsTime;
  TH2D*              h2_AmpvsAmp;
  TH2D*              h2_DeltaTimevsTime;
  TH2D*              h2_DeltaAmpvsAmp;

  TH2D*              h2_TimevsTimeIsNoise0;
  TH2D*              h2_AmpvsAmpIsNoise0;
  TH2D*              h2_DeltaTimevsTimeIsNoise0;
  TH2D*              h2_DeltaAmpvsAmpIsNoise0;

  TH2D*              h2_TimevsTimeIsNoise1;
  TH2D*              h2_AmpvsAmpIsNoise1;
  TH2D*              h2_DeltaTimevsTimeIsNoise1;
  TH2D*              h2_DeltaAmpvsAmpIsNoise1;
};

class TestHistoSADCMore
{
  public:
   TestHistoSADCMore   ( void ) {}

  std::vector<TH1D*> h1_NSamplesSADC;
  std::vector<TH2D*> h2_AmpSummVsMax;
  std::vector<TH2D*> h2_AmpMaxAdVsMax;
  std::vector<TH2D*> h2_TimeSummVsMax;
  std::vector<TH2D*> h2_TimeMaxAdVsMax;
  std::vector<TH2D*> h2_DAmpSummVsMax;
  std::vector<TH2D*> h2_DAmpMaxAdVsMax;
  std::vector<TH2D*> h2_DTimeSummVsMax;
  std::vector<TH2D*> h2_DTimeMaxAdVsMax;
  std::vector<TH2D*> h2_TimeSADCvsAmpMax;
  std::vector<TH1D*> h1_TimeSADC;
  std::vector<TH2D*> h2_TimeSADCvsTCS;
  std::vector<TH2D*> h2_SampleSADC;

  std::vector<TProfile*> p1_Shape10SADC;
  std::vector<TProfile*> p1_Shape100SADC;
  std::vector<TProfile*> p1_Shape500SADC;
  std::vector<TProfile*> p1_Shape1000SADC;
  std::vector<TProfile*> p1_ShapeNM10SADC;
  std::vector<TProfile*> p1_ShapeNM100SADC;
  std::vector<TProfile*> p1_ShapeNM500SADC;
  std::vector<TProfile*> p1_ShapeNM1000SADC;

  std::vector<TProfile*> p1_SADCMaxad2Ref;
  std::vector<TH2D*> h2_SampleSADCcells;
};

class TestHistoSADC_LED
{
  public:
   TestHistoSADC_LED   ( void ) {}

  std::vector<TH1D*> h1_NSamplesSADC;
  std::vector<TH2D*> h2_AmpSummVsMax;
  std::vector<TH2D*> h2_AmpMaxAdVsMax;
  std::vector<TH2D*> h2_TimeSummVsMax;
  std::vector<TH2D*> h2_TimeMaxAdVsMax;
  std::vector<TH2D*> h2_DAmpSummVsMax;
  std::vector<TH2D*> h2_DAmpMaxAdVsMax;
  std::vector<TH2D*> h2_DTimeSummVsMax;
  std::vector<TH2D*> h2_DTimeMaxAdVsMax;
  std::vector<TH2D*> h2_TimeSADCvsAmpMax;
  std::vector<TH1D*> h1_TimeSADC;
  std::vector<TH2D*> h2_TimeSADCvsTCS;
  std::vector<TProfile*> p1_Shape10SADC;
  std::vector<TProfile*> p1_Shape100SADC;
  std::vector<TProfile*> p1_Shape500SADC;
  std::vector<TProfile*> p1_Shape1000SADC;
  std::vector<TProfile*> p1_ShapeNM10SADC;
  std::vector<TProfile*> p1_ShapeNM100SADC;
  std::vector<TProfile*> p1_ShapeNM500SADC;
  std::vector<TProfile*> p1_ShapeNM1000SADC;
};

#endif // CsCalorimeterHist___include
