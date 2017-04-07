// $Id: CsSPMaker.cc,v 1.16 2010/01/28 12:51:26 tnagel Exp $

/*!
   \file    CsSPMaker.cc
   \brief   Compass event/event Space Point Maker
   \author  Hugo Pereira
   \version $Revision: 1.16 $
   \date    $Date: 2010/01/28 12:51:26 $
*/

#include "CsSPMaker.h"
#include "CsSTD.h"
#include "CsErrLog.h"
#include "CsOpt.h"
#include "CsGeom.h"

#include "CsSpacePoint.h"
#include "CsDetector.h"
#include "CsDetFamily.h"
#include "CsSPUtils.h"
#include "CsHistograms.h"

#include <cstring>

using namespace std;

//_____________________________________________________________________________
CsSPMaker* CsSPMaker::instance_ = 0;

CsSPMaker* CsSPMaker::Instance( void ) {
  if( instance_ == 0 ) instance_ = new CsSPMaker( );
  return( instance_ );
}

//_____________________________________________________________________________
CsSPMaker::CsSPMaker( void )
{
  // init families
  df_.clear();  
  if( !readDetFamilies() ) 
  CsErrLog::Instance()->mes(elError,"Error reading Detector Families.");
}

//_____________________________________________________________________________
bool CsSPMaker::readDetFamilies( void )
{
  list<string> fMember;          // configuration string family member. must have at least three items
  list<string>::iterator IfM;    // corresponding iterator
  
  list<CsDetector*> dets = CsGeom::Instance()->getDetectors(); 
  list<CsDetector*>::iterator Id;
  
  CsDetFamily *df = NULL;
  CsDetector *det = NULL;
  
  CsOpt* opt = CsOpt::Instance();

  string tag = "SP";
  string key = "Family";
  int oldfid = -1;
  int fid = -1;
  double correl = -1;
  
  cout << endl <<"======= detector families for space points ======"<<endl;

  //___ FIRST LOOK AFTER DETECTORS FOR EACH FAMILY ___
  while( opt->getOptRec( tag, key, fMember ) ) {
    
    if( fMember.size()<2 ) {
      CsErrLog::Instance()->mes(elWarning,"wrong format for family member.");
      continue;
    }
    
    string TBName;    
    int maxMlt = -1;
    
    int i = 0;
    for( IfM=fMember.begin(); IfM!=fMember.end(); IfM++, i++ ) 
    switch ( i ) {
    case 0: istringstream( (*IfM) ) >> fid;    break;
    case 1: TBName = *IfM;   break;
    case 2: istringstream( (*IfM) ) >> maxMlt; break;
    default: break;
    }
        
    // check fid
    if( fid < 0 ) {
      CsErrLog::Instance()->mes(elFatal,"negative family ID.");
      return false;
    }
    
    // look for detector matching TBName
    for (Id=dets.begin();Id!=dets.end();Id++)
    if ( (*Id)->GetTBName() == TBName ) det = (*Id);  
        
    // check det
    if( det == NULL ) {
      CsErrLog::Instance()->mes(elError,"no detector matching TBName.");
      continue;
    }
    
    // fill families
    
    if ( oldfid < 0 ) {
      
      // create first family
      df = new CsDetFamily( fid );    
      oldfid = fid;

    } else if( fid!=oldfid ) {
      
      // save old family
      if( df->detSize() < 4 ) {
        CsErrLog::Instance()->mes( elError , "Two small detector family" );  
        delete df;
      } else df_.push_back( df ); 

      // start a new one
      df = new CsDetFamily( fid );    
      oldfid = fid;

    } 
    
    if( det != NULL && df != NULL ) df->addDetector( *det, maxMlt );
  }  // next CsOpt Line
  
  // add last registered family to df_
  if( df == NULL || df->detSize() < 2 ) {
    CsErrLog::Instance()->mes( elError , "Troubles with last detector family" );  
    delete df;
  } else if( df->detSize() < 4 ) {
    CsErrLog::Instance()->mes( elError , "Two small detector family" );  
    delete df;
  } else df_.push_back( df ); 
  
  if( df_.empty() ){
    CsErrLog::Instance()->mes( elError , "No detector families booked." );
    return false;
  }

  //___ SECOND SET PARAMETERS FOR EACH FAMILY ___
  
  for( unsigned int i = 0; i < df_.size(); i++ ) {
    double cut;
    string word;
    
    if( CsSPUtils::Instance()->getOptForFamily( df_[i], "SP", "Name", word ) )
    df_[i]->setName( word );  
    
    if( !CsSPUtils::Instance()->getOptForFamily( df_[i], "SP", "zrec", cut ) ) { 
      cout << "family " << df_[i]->getID() << " zrec not set." << endl;
      return false;
    } else df_[i]->setZRec( cut );
    
    if( !CsSPUtils::Instance()->getOptForFamily( df_[i], "SP", "NClCut", cut ) ) {
      cout << "family " << df_[i]->getID() << " NClCut not set." << endl;
      return false;
    } else df_[i]->setNClCut( (int) cut );
    
    if( !CsSPUtils::Instance()->getOptForFamily( df_[i], "SP", "chi2Cut", cut ) ) 
    cout << "family " << df_[i]->getID() << " chi2cut not set." << endl;
    else df_[i]->setChi2Cut( cut );
     
    if( !CsSPUtils::Instance()->getOptForFamily( df_[i], "SP", "chi2Cut_Fast", cut ) ) 
    cout << "family " << df_[i]->getID() << " chi2cut_Fast not set." << endl;
    else df_[i]->setChi2Cut_Fast( cut );
     
    if( !CsSPUtils::Instance()->getOptForFamily( df_[i], "SP", "mode", word ) )
    cout << "family " << df_[i]->getID() << " mode not set." << endl;
    else {
      char smode[9];
      strncpy( smode, word.c_str(), 9 ); 
      if( strncmp( smode, "TPOINT", 6 ) == 0 ) df_[i]->setMode( CsDetFamily::TPOINT );
      else if( strncmp( smode, "STRAIGHT", 8 ) == 0 ) df_[i]->setMode( CsDetFamily::STRAIGHT );
      else { 
        cout << "family " << df_[i]->getID() << " wrong mode: " << word << "." << endl;
        return false;
      }
    }
    
    // detFamily geometry
    if( CsSPUtils::Instance()->getOptForFamily( df_[i], "SP", "geometry", fMember )  && 
      fMember.size() >= 4) {
      double xMin, xMax, yMin, yMax;
      int j = 0;
      for( IfM = fMember.begin(); IfM != fMember.end(); IfM++, j++ ) 
      switch (j) {
        case 0: istringstream( (*IfM) ) >> xMin; break;
        case 1: istringstream( (*IfM) ) >> xMax; break;
        case 2: istringstream( (*IfM) ) >> yMin; break;
        case 3: istringstream( (*IfM) ) >> yMax; break;
        default: break;
      }
      df_[i]->setGeometry( xMin, xMax, yMin, yMax );
    }
    
  }  
      
  //___FOURTH WRITE CONFIGURATION TO SCREEN___
  for( unsigned int i = 0; i < df_.size(); i++ ) df_[i]->dumpConfig();
  return true;
} 

//_____________________________________________________________________________
void CsSPMaker::cleanEvent( ) 
{
  for( unsigned int i=0; i< df_.size(); i++ ) df_[i]->cleanEvent();
  return;
}
	
