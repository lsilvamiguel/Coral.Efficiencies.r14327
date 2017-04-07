#ifndef DigitizerSADCN_h
#define DigitizerSADCN_h

#include <string>
#include <iostream>
#include <cstdio>
#include <list>
#include <vector>

#include <math.h>
#include <string>


#include "CsDigitizerSADC.h"

//namespace MN {

namespace Reco{ class StatInfo;}
class TabulatorShapeSADC;

class MyPulse;
class CsCalorimeter;

class ManagerShapeTableSADC
{
  private:
    const static size_t  NTABS;
    const static double  HFWmin;
    const static double  HFWmax;
    const static double  Xmin;
    const static double  Xmax;
    const static double  NXbins;

    class ProfileMaps
     {
        public:
           ProfileMaps(  std::map <double,double> * pmap,  std::map <double,std::pair<double,double > > *  wmap ) : pmap_(pmap), wmap_(wmap) { }

        public:
          std::map <double,double> * pmap_;
          std::map <double,std::pair<double,double > > * wmap_;
     };

  public:
     ManagerShapeTableSADC( CsCalorimeter *calo ) : calorimeter_(calo) {  stepHFW_ = (HFWmax-HFWmin)/(double)NTABS; stepX_=(Xmax-Xmin)/NXbins;}
     int           InputShapeTableSADC ( const std::string &s );
     bool        CheckTables( void ) const;
     double    Value( double hfw, double x ) const;
     std::pair<double, double>    WValue( double hfw, double x ) const;
  private:
     bool        AddTable( double hfwmin, const std::vector <int> &profile );
     const std::map<double,double> *  GetTab( double hfw ) const;
     const std::map <double,std::pair<double,double > > *  GetWTab( double hfw ) const;
//     const  ProfileMaps &GetTabs( double hfw ) const;

     double     ValueTab( double x, const std::map<double,double> * tab ) const;
     std::pair<double,double> WValueTab( double v, const std::map<double, std::pair<double,double> > * tab ) const;

  private:
     double      stepHFW_;
     double      stepX_;

//  Complicated code, not clean structure all pair,map,pair, double!! Must be improved or at least deeply encapsulated
     std::map <double,  ProfileMaps   > vtab_;
     CsCalorimeter                        *calorimeter_;
};

////////////////////////////////////////////////////////////////////////////////
//  TODO Implement Development version

class DigitizerSADCN : public DigitizerSADCBase
{
  public:
    virtual            ~DigitizerSADCN() ;

    DigitizerSADCN ( CsCalorimeter* parent, Type type ) : DigitizerSADCBase(parent,  type ),mtab_(NULL) { SetDefaultCalibrations(); }

    virtual bool FitBase(  const std::vector<uint16> &sample );

    virtual bool Check(  void ) const;

  private:
    void Init( void );
    void FitMaxAdvancedN(  const std::vector<uint16> &sample );
    void SetDefaultCalibrations ( void );

  public:
    /// \return FWHM (ns)
    double   GetFWHMTime        ( void ) const { return result_.timeFWHM; }
// More like development functions
    /// \return FWHM (ns)
    double   GetHFWTime        ( void ) const { return result_.timeHFW; }
    double   GetTime1        ( void ) const { return result_.time1; }
    double   GetTime2        ( void ) const { return result_.time2; }
    double   GetFWHMTimeNorm ( void ) const;
    double   GetSigmaFWHMTime ( void ) const;
    double   GetHFWTimeNorm    ( void ) const ;

    double   GetAmpOver( void ) const {return result_.ampl_overflow;}
    double   GetTimeOver( void ) const {return result_.time_overflow;}

    ///
    bool   BadShape       ( void ) const;
    bool   Fit2dcutFWHMnHFW   ( double r ) const;
    double   GetRadiusFWHMnHFW   ( void ) const;

  private:
   std::pair<double,double> CalcTime2_halfMax(const std::vector<short unsigned int> &sample, const unsigned int max_position);
// Interface function to extarct value from the shape table, aimed to be used in the Fit
   double                           GetTabShapeValue( double time ) const;
   std::pair<double,double>  GetTabShapeWValue( double time ) const;

// Calibrations
  public:
   // Calibrations, shape parametrization parameters
    double                      calib_fwhm_tab_;
    double                      calib_fwhm_sigma_tab_;  // not used yet
    double                      calib_hfwn_tab_;
    double                      calib_hfwn_sigma_tab_;
// Noise selection Calibrations
    double                      parFWHMn_HFW_fwhmn_c_;
    double                      parFWHMn_HFW_fwhmn_r_;
    double                      parFWHMn_HFW_hfw_c_;
    double                      parFWHMn_HFW_hfw_r_;

// Add this object temporary not to break code transition
    const ManagerShapeTableSADC *mtab_;
};

////////////////////////////////////////////////////////////////////////////////
//  Implementation of Development version
class DigitizerSADCNdevel : public DigitizerSADCN
{
  public:
    virtual            ~DigitizerSADCNdevel() {}

    DigitizerSADCNdevel ( CsCalorimeter* parent, Type type ) : DigitizerSADCN(parent,  type ) { }
  public:
   void AddStatCalibTime2AndHFWn ( double time2, double hfwn_sadc ) { stat_calibTime2_.Add(time2); stat_calibhfwn_.Add(hfwn_sadc);}
//    void AddStatCalibTime2 ( double time2 ) { stat_calibTime2_.Add(time2);}
//    void AddStatCalibHFWn ( double hfwn_sadc ) { stat_calibhfwn_.Add(hfwn_sadc);}
// StatInfo Store Development
  void SetTagOverflowValuesToArtificial ( void )
  {
// we are changing overflow_amplitude_  parameter
    if( GetParentName() == "EC02P1__" || GetParentName() == "EC00P1__" )
        overflow_amplitude_ = 1023;
    else if( GetParentName() == "HC01P1__" || GetParentName() ==  "HC02P1__" )
        overflow_amplitude_ = 511;
    else
        overflow_amplitude_ = 253;
  }
  public:
   Reco::StatInfo                   stat_calibTime2_;
   Reco::StatInfo                   stat_calibhfwn_;
};

//} // namespace MN

#endif  //         DigitizerSADCN_h
