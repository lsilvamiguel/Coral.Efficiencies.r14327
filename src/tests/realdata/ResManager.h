// $Id: ResManager.h,v 1.3 2008/11/03 08:49:48 aaust Exp $
 
/*!
   \file    ResManager.h
   \brief   Resolution Managment Interface Class.
   \author  Hugo Pereira
   \version $Revision: 1.3 $
   \date    $Date: 2008/11/03 08:49:48 $
*/

#ifndef ResManager_h
#define ResManager_h

#include <vector>

#include <TROOT.h>
#include <TObject.h>

#include "HManager.h"

/*!
   \class   ResManager
   \brief   Resolution Managment Interface Class.
*/

class ResManager: public HManager {

  public:
  ResManager( const HManager h ): HManager( h ) {} 
  
  /*! \fn  ResManager( 
    const char* trackfileselection="", 
    const char* detectorfile="", 
    bool magnets_on=false );
    \brief constructor.
    \param trackfileselection input tree(s) obtained from coral. Wildcards are accepted.
    \param detectorfile the name of the det.dat file to be loaded
  */
  ResManager( 
    const char* fileselection="", 
    const char* detectorfile=""  ): 
    HManager( fileselection, detectorfile ) {} 
  
  //! fit residual with 1 or 2 gaussian depending on OneGaus_ value
  void FitDU(      const char* selection = "*", const char* cut = "T_fnd", int nEnt = 0 );
  void FitDV(      const char* selection = "*", const char* cut = "T_fnd", int nEnt = 0 );
  void FitDU_MWPC( const char* selection = "*", const char* cut = "T_fnd", int nEnt = 0 );
  void FitDUvsU(   const char* selection = "*", const char* cut = "T_fnd", int nEnt = 0 );
  void FitDUvsV(   const char* selection = "*", const char* cut = "T_fnd", int nEnt = 0 );
  void DrawDU_LR(  const char* selection = "*", const char*  cut = "T_fnd", 
    double OffsetN = 0, 
    double OffsetP = 0, 
    int nEnt = 0 );
   bool FitResolutions( const char* selection, const char* cut = "T_fnd", int nEnt = 0 );
  
  protected:
  double _ResCalc( unsigned int it, std::vector< double > res );
  friend class Fit;
   
  ClassDef(ResManager,0)
  
};

#endif
