// $Id: RTManager.h,v 1.1 2002/12/09 08:50:39 hpereira Exp $
 
/*!
   \file    RTManager.h
   \brief   RT Relation Managment Interface Class.
   \author  Hugo Pereira
   \version $Revision: 1.1 $
   \date    $Date: 2002/12/09 08:50:39 $
*/

#ifndef RTManager_h
#define RTManager_h

#include "HManager.h"

#include <string>

#include <TROOT.h>
#include <TObject.h>
#include <vector>

class DetectorInfo;
class RTGridPoint;
class RTInfo;
class TH2;

 
/*!
   \class   RTManager
   \brief   RT Relation Managment Interface Class.
*/

class RTManager: public HManager {
  public:
  RTManager( const HManager h ): HManager( h ), 
    H_single_RT_(0),
    H_single_RT_abs_(0),
    single_rtInfo_( 0 ) {}
  
  /*! \fn  RTManager( 
    const char* trackfileselection="", 
    const char* detectorfile="", 
    bool magnets_on=false );
    \brief constructor.
    \param trackfileselection input tree(s) obtained from coral. Wildcards are accepted.
    \param detectorfile the name of the det.dat file to be loaded
  */
  RTManager( 
    const char* fileselection="", 
    const char* detectorfile=""  ): 
    HManager( fileselection, detectorfile ), single_rtInfo_( 0 ) {} 
  
  //! dump selected RT relation infos to screen, if any
  void DumpRTInfo( const char* selection = "*" );
  
  //! Fit RT Relation with straight lines for selected detectors to get T0
  void FitRTStraight( const char* selection = "*", const char* cut = "T_fnd", double rMin = 1, double rMax = -1, int nEnt = 0 );
  
  //! Fit RT Relation with polynom for selected detectors
  void FitRTRelation( const char* selection = "*", const char* cut = "T_fnd", int nEnt = 0 );

  //! Fit RT Relation slice by slice with 2 gaussians for selected detectors
  void FitRTGrid( const char* selection = "*", const char* cut = "T_fnd", unsigned int nBin = 0, int nEnt = 0 );

  //! write RT relation fits to calibration files
  void WriteRTParsToDB ( const char* selection = "*", const bool useMyT0 = false, const double myT0 = 0 );

  //! write RT relation grids to calibration files
  void WriteRTGridsToDB( const char* selection = "*", const bool useMyT0 = false, const double myT0 = 0 );

  //! Adds RT Relations for a set of detectors
  void BuildSingleRT( const char* selection = "*", const char* cut = "T_fnd", int nEnt = 0 );

  //! fits sum of RT Relations corresponding to a set of detectors using straight lines
  void FitSingleRTStraight( double rMin = 1, double rMax = -1  ) {_FitRTStraight( H_single_RT_, single_rtInfo_, rMin, rMax ); }

  //! fits sum of RT Relations corresponding to a set of detectors slice by slice using two gaussians
  void FitSingleRTGrid( unsigned int nBins = 0 ) {_FitRTGrid( H_single_RT_, single_rtInfo_, nBins ); }

  //! fits sum of RT Relations corresponding to a set of detectors slice by slice using two gaussians
  void FitSingleRTRelation( void ) {_FitRTRelation( H_single_RT_abs_, single_rtInfo_ ); }

  //! write summed RT relation fits to calibration files
  void WriteSingleRTParsToDB ( const char* selection = "*", const bool useMyT0 = false, const double myT0 = 0 );

  //! write summed RT relation grids to calibration files
  void WriteSingleRTGridsToDB( const char* selection = "*", const bool useMyT0 = false, const double myT0 = 0 );
  
  //! set validity for calibration files associated to RT relation
  void SetValidity ( const char* selection = "*", const char* start = "YYYY-MM-DD-hh:mm:ss", const char* stop = "YYYY-MM-DD-hh:mm:ss" );
 
  private:  
  
  //! fit 2D histogram using straight line. Store the result in RTInfo object
  void _FitRTStraight( TH2* h, RTInfo* rt, double rMin = 1, double rMax = -1  );

  //! Fit 2D histogram using polynom. Store the result in RTInfo object
  void _FitRTRelation( TH2* h, RTInfo* rt );
  
  //! Fit 2D histogram slice by slice with 2 gaussians. Store the result in RTInfo object
  void _FitRTGrid( TH2* h, RTInfo* rt, unsigned int nBins = 0 );
  
  TH2*    H_single_RT_;    //!< summed RT relation 2D histogram for a given set of detector
  TH2*    H_single_RT_abs_; //!< summed RT relation 2D histogram for a given set of detector
  RTInfo* single_rtInfo_;  //!< summed RT relation RTInfo object for a given set of detector
  
  ClassDef(RTManager,0)
  
};

#endif
