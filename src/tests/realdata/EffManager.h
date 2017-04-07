// $Id: EffManager.h,v 1.2 2007/06/01 09:13:51 conrad Exp $

/*!
   \file    EffManager.h
   \brief   Efficiency Managment Interface Class.
   \author  Hugo Pereira
   \version $Revision: 1.2 $
   \date    $Date: 2007/06/01 09:13:51 $
*/

#ifndef EffManager_h
#define EffManager_h

#include <TROOT.h>
#include <TObject.h>
#include <TCut.h>

#include "HManager.h"

#include <iostream>

/*!
   \class   EffManager
   \brief   Efficiency Managment Interface Class.
*/
//_______________________________________________________________________________
class EffManager: public HManager {
  public:
  EffManager( const HManager h ): HManager( h ) {} //!< 'Extend' constructor
  
  /*! \fn  EffManager( 
    const char* trackfileselection="", 
    const char* detectorfile="", 
    bool magnets_on=false );
    \brief constructor.
    \param trackfileselection input tree(s) obtained from coral. Wildcards are accepted.
    \param detectorfile the name of the det.dat file to be loaded
  */
  EffManager( 
    const char* fileselection="", 
    const char* detectorfile=""  ): 
    HManager( fileselection, detectorfile ), 
    cut_("T_fnd") {}

  /*! \fn DrawEff( const char* var, const char* selection = "*", TCut cut="1", const char* opt="", int ent=0  )
    \brief draw efficiency plot wrt any tree variable for a set of detectors
    \param var the variable wrt which the eff is drawn
    \param selection the detectors TBNames for wich efficiency is drawn. example: "St03X1*b", "GM01X1__"
    \param cut track selection cut
    \param opt drawing options
    \param ent number of selected tree entries.
  */
  void DrawEff( const char* var, const char* selection = "*", TCut cut="1", const char* opt="", int ent=0 );

  /*! \fn DrawEffUProf( const char* var, const char* selection = "*", TCut cut="1", const char* opt="", int ent=0  )
    \brief draw profiles for all tracks matching cut and tracks for which the det is found efficient on the same plot
    \param var the variable wrt which the eff is drawn
    \param selection the detectors TBNames for wich efficiency is drawn. example: "St03X1*b", "GM01X1__"
    \param cut track selection cut
    \param ent number of selected tree entries.
  */
  void DrawEffUProf(  const char* selection = "*", TCut cut="T_inActive", const char* opt="", int nEnt = 0 );
  
  //! set the cut used to tell if a detector is efficient or not  
  inline void SetEffCut( TCut cut="T_fnd" ) { cut_ = cut; }
  
  private:
  TCut cut_; //!< cut used to tel if a detector is efficient or not
  ClassDef(EffManager,0)
  
};

#endif
