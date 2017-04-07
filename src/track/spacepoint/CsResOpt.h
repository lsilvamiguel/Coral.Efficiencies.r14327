// $Id: CsResOpt.h,v 1.6 2003/04/24 07:23:26 benigno Exp $

/*!
   \file    CsResOpt.h
   \brief   Compass Resolution optimiser from Space point Residuals
   \author  Hugo Pereira
   \version $Revision: 1.6 $
   \date    $Date: 2003/04/24 07:23:26 $
*/


#ifndef CsResOpt_h
#define CsResOpt_h

#include "CsSTD.h"
#include "CsDetFamily.h"
class CsDetector;

//___________________________________________________________________________
class CsResOpt{
 public:
  CsResOpt( const CsDetFamily& df ); //!< The constructor
  
  bool setResidsFromHists( void );
  bool setResidsFrom2DHists( void );
  bool setResidsByHand( std::vector< double > resid );
  bool setResolutionsByHand( std::vector< double > res );
  void testResCalc( void );

  inline bool resInitOK( void )               const { return resIOK_; }
  inline double getResAve_Init( void )        const { return resIAve_; }
  inline std::vector< double > getRes_Init( void ) const { return resI_; }

  inline bool resFitOK( void )               const { return resFOK_; }
  inline double getResAve_Fit( void )        const { return resFAve_; }
  inline std::vector< double > getRes_Fit( void ) const { return resF_; }
  
  bool fitResAve( double& resAve, double& resAveErr );
  inline bool fitRes( void ) { return fitRes( stdout ); }
  bool fitRes( FILE* resOut );
  
  double getZoptX( double res );  // optimal z position for x resolution measurement
  double getZoptY( double res );  // optimal z position for y resolution measurement
  double getZoptX( std::vector< double > res );  // optimal z position for x resolution measurement
  double getZoptY( std::vector< double > res );  // optimal z position for y resolution measurement
  bool dumpResAtZ( double z, double res );           // Det family resolution at z
  bool dumpResAtZ( double z, std::vector< double > res ); // Det family resolution at z

  double resCalc( unsigned int it, std::vector< double > res );
  double resCalc( unsigned int it, double res );
  inline const CsDetFamily* getDetFamily( void ) const { return df_; }
  inline std::vector<CsDetector*> getDetectors( void ) const { return df_->getDetectors(); }
  inline std::vector<double> getResids( void ) const { return resid_; } 
 private:
  const CsDetFamily* df_;

  std::vector< double > resid_;
  double resIAve_;         // initial average resolution ( from detectors info )
  double resFAve_;         // final average resolution   ( from fit )
  std::vector< double > resI_;  // initial resolutions        ( from detectors info )
  std::vector< double > resF_;  // final resolutions          ( from fit )

  bool residOK_;
  bool resIOK_;
  bool resFOK_;
};

//___________________________________________________________________________
extern void fcnAve( int &npar, double *gin, double &f, double *X, int flag );
extern void fcn( int &npar, double *gin, double &f, double *x, int flag );

#endif
 
