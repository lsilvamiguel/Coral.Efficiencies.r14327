// $Id: AddOffset.cc,v 1.5 2011/01/18 14:59:06 fsozzi Exp $

/*!
   \file    AddOffset.cc
   \brief   some macro to add offset to gems, micromegas by hand on alignment file when switching from alignment run to physics run
   \author  Hugo Pereira
   \version $Revision: 1.5 $
   \date    $Date: 2011/01/18 14:59:06 $
*/

#include "Macro.h"
#include "DetFileManager.h"
#include "DetectorInfo.h"
#include "Utils.h"

using namespace std;

//________________________________________
void Macro::AddGemOffset( DetFileManager &df )
{
  DetectorInfo *det = 0;
  
  //=== Values taken from fieldoff 102486 to fieldon 102479
  cout << "Macro::AddGemOffset - \"" << df.detFileName_ << "\". GM01.\n";
  if( ( det = df.GetDetInfo( "GM01U1" ) ) ) det->DUtoCenter( -0.277 );
  if( ( det = df.GetDetInfo( "GM01V1" ) ) ) det->DUtoCenter( -0.347 );
  if( ( det = df.GetDetInfo( "GM01X1" ) ) ) det->DUtoCenter(  0.380 );
  df.MatchCenters( "GM01X","GM01Y");
  df.MatchCenters( "GM01U","GM01V");
  df.MoveDeadZone( "GM01" );

  cout << "Macro::AddGemOffset - \"" << df.detFileName_ << "\". GM02.\n";
  if( ( det = df.GetDetInfo( "GM02U1" ) ) ) det->DUtoCenter( -0.218 );
  if( ( det = df.GetDetInfo( "GM02V1" ) ) ) det->DUtoCenter( -0.195 );
  if( ( det = df.GetDetInfo( "GM02X1" ) ) ) det->DUtoCenter(  0.083 );
  df.MatchCenters( "GM02X","GM02Y");
  df.MatchCenters( "GM02U","GM02V");
  df.MoveDeadZone( "GM02" );

  cout << "Macro::AddGemOffset - \"" << df.detFileName_ << "\". GM03.\n";
  if( ( det = df.GetDetInfo( "GM03U1" ) ) ) det->DUtoCenter( -0.145 );
  if( ( det = df.GetDetInfo( "GM03V1" ) ) ) det->DUtoCenter( -0.125 );
  if( ( det = df.GetDetInfo( "GM03X1" ) ) ) det->DUtoCenter( -0.010 );
  df.MatchCenters( "GM03X","GM03Y");
  df.MatchCenters( "GM03U","GM03V");
  df.MoveDeadZone( "GM03" );

  return;
}

//________________________________________
void Macro::AddMMOffset( DetFileManager &df )
{
  DetectorInfo *det = 0;
  //=== Values taken from fieldoff 102486 to fieldon 102479

  cout << "Macro::AddMMOffset - \"" << df.detFileName_ << "\". MM03.\n";
  if( ( det = df.GetDetInfo( "MM03U1" ) ) ) for( unsigned int i=0; i< det->sub_.size(); i++ ) det->sub_[i]->DUtoCenter(  0.393 );
  if( ( det = df.GetDetInfo( "MM03V1" ) ) ) for( unsigned int i=0; i< det->sub_.size(); i++ ) det->sub_[i]->DUtoCenter( -0.447 );
  if( ( det = df.GetDetInfo( "MM03X1" ) ) ) for( unsigned int i=0; i< det->sub_.size(); i++ ) det->sub_[i]->DUtoCenter( -0.631 );
  df.MatchCenters( "MM03X","MM03Y");
  df.MatchCenters( "MM03U","MM03V");
  df.MoveDeadZone( "MM03" );
  

  return;
}
 
