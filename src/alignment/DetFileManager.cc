// $Id: DetFileManager.cc,v 1.34 2008/06/05 13:08:55 rgazda Exp $
/*!
  \file    DetFileManager.cc
  \brief   read, manipulate, format detectors and dead zones in det.dat file
  \author  Hugo Pereira
  \version $Revision: 1.34 $
  \date    $Date: 2008/06/05 13:08:55 $
*/

#include "DetFileManager.h"
#include "Utils.h"
#include "Defs.h"
#include "Obj3D.h"
#include "Point.h"
#include "Utils.h"

#include <map>
#include <stdlib.h>
#include <unistd.h>
#include <stdio.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <list>
#include <functional>
#include <TRotMatrix.h>
#include <TCanvas.h>

#include <TBRIK.h>
#include <TNode.h>

using namespace std;

char* operator+( std::streampos&, char* );
//_______________________________________________________________________________
ClassImp( DetFileManager )

  DetFileManager::DetFileManager( void ): TObject() {
    // not useless !
    // DetFileManager::Streamer will be called through TObject
  }

DetFileManager::DetFileManager( const char* detFileName ): TObject() {
  vector<string> files = Utils::GetFiles( detFileName );
  if( !files.size() ) {
    cout << "DetFileManager::DetFileManager - ERROR: cannot read file \"" 
         << detFileName_  << "\"." 
         << endl;
    return;
  }
  
  detFileName_ = files[0];  

  // set detFileContents_ to NULL to be able later to check if it points somewhere (jj)
  detFileContents_=NULL;
  Init();
}

void DetFileManager::Init() {
  
  vector< int > nDetProj; // number of detectors / projection
  
  // Check detFileName_
  if(!detFileName_.size() ) return;
  DetInfo_.clear();
  DeadZInfo_.clear();
  
  std::ifstream in( detFileName_.c_str(), ios::in );
  if( !in ) {
    cout << "DetFileManager::DetFileManager - ERROR: cannot read file \"" 
         << detFileName_  << "\"." 
         << endl;
    return;
  }
  // std::ostringstream fileconts;
  ostringstream fileconts;
  cout << "DetFileManager::DetFileManager - INFO: reading file \"" << detFileName_ << "\".\n";   

  const int linesize = 512;
  char* line = new char[linesize];
  unsigned int lineNumber = 0;
  int pos=0;    
  while(!in.eof()) {
    in.seekg(pos);
    in.getline(line,linesize, '\n');
    if(!in) break;
    // pos = in.tellg() + 1;
    pos = in.tellg();
    // pos = max(0,pos);
    // pos++;
    lineNumber++; 
 
    fileconts<<line<<endl;
    
    if( line[0] == '\0' ) { continue; }
    
    istringstream s(line);
    string opt;  s >> opt;

    // magnet info
    if( opt == "mag" ){ MagnetInfo_.push_back( new MagnetInfo( detFileName_, lineNumber, line ) );
    }

    // detector info
    else if( opt == "det" ) {       
      DetectorInfo *di =  new DetectorInfo( detFileName_, lineNumber, line );
      vector< DetectorInfo* > tmp;
      
      // special case for straws (last TBName digit is ignored)
      // not true right now (jj)
      // if( di->TBName_.substr(0,2) == "ST" ) {
      //   string TBTmp = di->TBName_.substr(0,7);
      //   tmp = _GetDetSelection( TBTmp );
      // }
   
      // special case for muon wall2 (last TBName digit is ignored)
      if( di->TBName_.substr(0,2) == "MB" ) { // changed else if to if (jj)
        string TBTmp = di->TBName_.substr(0,7);
        tmp = _GetDetSelection( TBTmp );
      }  
      
      // standard case
      else tmp = _GetDetSelection( di->TBName_ );
      switch( tmp.size() ) {
      case 0: {
        DetInfo_.push_back( di );

	// see if detector corresponds to new projection
        double ang_d = double( atan2(di->sinTheta_,di->cosTheta_) );
        int ang_i = int(ang_d*180/PI + 0.5 * (ang_d < 0 ? -1 : 1));
        unsigned int iproj = 0;
        for(; iproj<projections_.size(); iproj++ )
	  if( projections_[iproj]==ang_i ) {
	    di->projection_ = iproj;
	    nDetProj[iproj]++;
	    break;
	  }     
        if( iproj == projections_.size() ) {
          di->projection_ = iproj;
          projections_.push_back( ang_i );
          nDetProj.push_back( 1 );
        }
      }
	break;

      case 1:
        // copy projection id from parent
        di->projection_ = tmp.front()->projection_;
        
        // add to parent list of childs
        tmp.front()->Add( di );
        break;
      
      default:
        printf( "DetFileManager::DetFileManager - FATAL - line %4i several dets with TBName \"%s\" saved independent.\n",
		lineNumber, 
		tmp.front()->TBName_.c_str() );        
        exit( 0 );
      }   
      // dead zone info
    } else if( opt == "dead" ) {
      DeadZoneInfo* dead = new DeadZoneInfo( detFileName_, lineNumber, line );
      vector< DetectorInfo* > tmp = _GetDetSelection( dead->TBName_ );
      switch( tmp.size() ) {
      case 0:  
        printf( "DetFileManager::DetFileManager - FATAL - line %4i no detector associated to deadZone \"%s\".\n", lineNumber, dead->TBName_.c_str() );
        exit( 0 );
        break;
       
      case 1:
        // link detectorInfo and DeadZone
        tmp.front()->GetMain()->dead_ = dead;
        dead->det_ = tmp.front()->GetMain();
        
        // store dead zone into list     
        DeadZInfo_.push_back( dead );          
        break;
     
      default:  
        printf( "DetFileManager::DetFileManager - FATAL - line %4i several independant detectors associated to deadZone \"%s\".\n", 
		lineNumber, 
		dead->TBName_.c_str() );
        exit( 0 );
      }
    }       
  }
  // delete line
  delete[] line;
  _SortVects();

  // first check if data on detFileContents_ wasn't already created (jj)
  if(detFileContents_) delete[] detFileContents_;
  detFileContents_ = new char[ strlen(fileconts.str().c_str())+1];
  strcpy(detFileContents_ , fileconts.str().c_str());

  cout << "DetFileManager::DetFileManager - INFO: " << DetInfo_.size() << " detectors loaded.\n";
  cout << "DetFileManager::DetFileManager - INFO: " << DeadZInfo_.size() << " dead zones loaded.\n";
  cout << "DetFileManager::DetFileManager - WARNING: all length have been converted to [mm].\n";
  return;
}

#ifndef __CINT__  
//_______________________________________________________________________________
vector< DetectorInfo* > DetFileManager::GetDetInfo( void )
{
  vector< DetectorInfo* > out;
  for( unsigned int id=0; id< DetInfo_.size(); id++ ) 
    out.push_back( DetInfo_[id]->GetMain() );
  return out;
} 

//_______________________________________________________________________________
vector< DetectorInfo* > DetFileManager::GetAllDetInfo( void )
{
  vector< DetectorInfo* > out;
  for( unsigned int id=0; id< DetInfo_.size(); id++ ) 
    for( unsigned int is=0; is< DetInfo_[id]->sub_.size(); is++ )
      out.push_back( DetInfo_[id]->sub_[is] );
  return out;
} 
#endif

//_______________________________________________________________________________
void DetFileManager::Dump( const char* detselection )
{
  vector< DetectorInfo* > sel = _GetDetSelection( detselection );
  for( unsigned int i = 0; i < sel.size(); i++ ) 
    cout << "DetFileManager::Dump - " << i << " - Got \"" << sel[i]->GetMain()->TBName_ << "\".\n";

}

//_______________________________________________________________________________
void DetFileManager::DumpAll( const char* detselection )
{
  unsigned int n=0;
  vector< DetectorInfo* > sel = _GetDetSelection( detselection );
  for( unsigned int i = 0; i < sel.size(); i++ ) 
    for( unsigned int j = 0; j < sel[i]->sub_.size(); j++, n++ )
      cout << "DetFileManager::Dump - " << n << " - Got \"" << sel[i]->sub_[j]->TBName_ << "\".\n";
}


//_______________________________________________________________________________
bool DetFileManager::FixWirD( const char* detselection, double value )
{
  vector< DetectorInfo* > sel = _GetDetSelection( detselection );
  for( unsigned int i = 0; i < sel.size(); i++ ) {
    DetectorInfo *di = sel[i]->GetMain();
    double du = value-di->wirD_;
    for( unsigned int j=0; j < di->sub_.size(); j++ ) {
      printf("DetFileManager::FixWirD - \"%s\" [%9.3f %9.3f %9.3f]",
	     di->sub_[j]->TBName_.c_str(),
	     di->sub_[j]->wirD_,
	     di->sub_[j]->xcm_,
	     di->sub_[j]->ycm_ );
      di->sub_[j]->xcm_-=di->sub_[j]->irotM_(0,0)*du;
      di->sub_[j]->ycm_-=di->sub_[j]->irotM_(0,1)*du;
      di->sub_[j]->wirD_ += du;
      printf( " ->  [%9.3f %9.3f %9.3f].\n",
	      di->sub_[j]->wirD_,
	      di->sub_[j]->xcm_,
	      di->sub_[j]->ycm_ );
    }
  }
  return true;
}    


//_______________________________________________________________________________
bool DetFileManager::DumpToFile( const char* file )
{

  //! get output file name
  string outFile = (file) ? string( file ): detFileName_;
  cout << "DetFileManager::DumpToFile - File: \"" << outFile << "\".\n";
  
  //! make output backup, if needed
  string svFile = Utils::MakeBackup( outFile );
      
  // update input file 
  ifstream in( (outFile==detFileName_) ? svFile.c_str():detFileName_.c_str(), ios::in );
  ofstream out( outFile.c_str(), ios::out );
  cout<<" OUTPUT = "<<outFile.c_str()<<endl;

  char* line = new char[512];
  unsigned int lineNumber = 0;  
  while( !in.eof() ) { 
    in.getline( line, 512, '\n');
    if( in.eof() || !in.good() ) continue;
    lineNumber++;
    // cout<<" OUT0 = "<<lineNumber<<endl;    
    if( line[0] == '\0' ) {
      out << line << endl;
      continue;
    }
    // cout<<" OUT1 = "<<lineNumber<<endl;

    bool found = false;
    for( unsigned int i=0; i < DetInfo_.size(); i++ )
      for( unsigned int j=0; j < DetInfo_[i]->sub_.size(); j++ )
	if( lineNumber == DetInfo_[i]->sub_[j]->lineNumber_ ) {
	  found = true;
	  out << DetInfo_[i]->sub_[j]->Dump() << endl;
	}
    // cout<<" OUT2 = "<<lineNumber<<endl;
    if( !found )
      for( unsigned int i=0; i < DeadZInfo_.size(); i++ )
	if( lineNumber == DeadZInfo_[i]->lineNumber_ ) {
	  found=true;
	  out << DeadZInfo_[i]->Dump() << endl;
	}
    // cout<<" OUT3 = "<<lineNumber<<endl;
    if(!found) out << line << endl;
  }
  // cout<<" OUT4 = "<<lineNumber<<endl;
  
  out.close();
  return true;
}  

//_______________________________________________________________________________
bool DetFileManager::Sort( const char* selection )
{
  // Sort detectors
  vector< DetectorInfo* > sel = _GetDetSelection( selection );
  if( !sel.size() ) return false;
  
  list<DetectorInfo* > selList; 
  list<unsigned int> lines;
  
  // store DetectorInfos in a list and corresponding lines in a vector
  for( unsigned int i=0; i< sel.size(); i++ )
    for( unsigned int j=0; j< sel[i]->sub_.size(); j++ ) {
      selList.push_back( sel[i]->sub_[j] );
      lines.push_back( sel[i]->sub_[j]->lineNumber_ );
    }
  
  // Sort the list
  selList.sort( DetFileManager::sortDet_() );
  lines.sort();
  
  // Reatribute Line Number
  list<unsigned int>::iterator Il=lines.begin();
  for( list< DetectorInfo* >::iterator Id=selList.begin(); Id != selList.end(); Id++, Il++ ) {
    printf( "DetFileManager::SortDetInfo - INFO: %s (%10.4f) %4i ->", (*Id)->TBName_.c_str(), (*Id)->zcm_, (*Id)->lineNumber_ );
    (*Id)->lineNumber_ = (*Il);
    printf( "%4i.\n", (*Id)->lineNumber_ );
  }
  
  // Sort dead zones
  vector< DeadZoneInfo* > seld = _GetDeadSelection( selection );
  if( !seld.size() ) return false;
  
  list<DeadZoneInfo* > seldList; 
  lines.clear();
  
  // store deadZoneInfos and corresponding lines in lists
  for( unsigned int i=0; i< seld.size(); i++ ) {
    seldList.push_back( seld[i] );
    lines.push_back( seld[i]->lineNumber_ );
  }
  
  // Sort the list
  seldList.sort( DetFileManager::sortDeadZ_() );
  lines.sort();
  
  // Reatribute Line Number
  Il=lines.begin();
  for( list< DeadZoneInfo* >::iterator Id=seldList.begin(); Id != seldList.end(); Id++, Il++ ) {
    printf( "DetFileManager::SortDeadZInfo - INFO: %s (%10.4f) %4i ->", (*Id)->TBName_.c_str(), (*Id)->zcm_, (*Id)->lineNumber_ );
    (*Id)->lineNumber_ = (*Il);
    printf( "%4i.\n", (*Id)->lineNumber_ );
  }
  
  return true;
}

//_______________________________________________________________________________
bool DetFileManager::MatchCenters( const char* det0, const char* det1 )
{
  map<DetectorInfo*, DetectorInfo* > sel;

  // Get Detector infos
  vector< DetectorInfo* > sel0;
  vector< DetectorInfo* > sel1;
  for( unsigned int i0=0; i0 < DetInfo_.size(); i0++ ) {
    char key0[8]; 
    strcpy( key0,"");
    bool accept = true;
    if( strlen(det0) > DetInfo_[i0]->TBName_.size() ) accept = false;
    else for( unsigned int j=0; j < strlen(det0); j++ ) {
      if( det0[j]=='*' ) sprintf(key0,"%s%c",key0,DetInfo_[i0]->TBName_[j] );
      else if( det0[j] != DetInfo_[i0]->TBName_[j] ) accept = false;
    }
      
    if(!accept) continue;
    for( unsigned int i1=0; i1< DetInfo_.size(); i1++ ) {
      accept = true;
      char key1[8]; 
      strcpy( key1,"");
      if( strlen(det1) > DetInfo_[i1]->TBName_.size() ) accept = false;
      else for( unsigned int j=0; j < strlen(det1); j++ ) {
        if( det1[j]=='*' ) sprintf(key1,"%s%c",key1,DetInfo_[i1]->TBName_[j] );
        else if( det1[j]-DetInfo_[i1]->TBName_[j] ) { accept = false; break;}
      }
      if( !(accept && ( (strlen(key0)==0 && strlen(key1)==0) || strncmp( key0, key1, strlen( key1 ) ) == 0 ) ) ) continue;
      sel.insert( pair<DetectorInfo*, DetectorInfo* >( DetInfo_[i0], DetInfo_[i1] ) );
    }
  }
  
  // Matrix equations
  for( map<DetectorInfo*,DetectorInfo*>::iterator I = sel.begin(); I != sel.end(); I++ ) {
    DetectorInfo *di0 = (*I).first;  DetectorInfo *di0_M = di0->GetMain();
    DetectorInfo *di1 = (*I).second; DetectorInfo *di1_M = di1->GetMain();
    //printf("\nDetFileManager::MatchCenters - [\"%s\",\"%s\"].\n",
    //	   di0->TBName_.c_str(),
    //	   di1->TBName_.c_str());
    
    // To calculate the new center positions, parameters are taken from _Main_ detectors corresponding to selected detectors  
    // This to handle properly the case of detectors booked with subdetectors
    TMatrix A(2,2);
    A(0,0) = di0_M->irotM_(1,0);  A(0,1) = -di1_M->irotM_(1,0);
    A(1,0) = di0_M->irotM_(1,1);  A(1,1) = -di1_M->irotM_(1,1);
  
    TMatrix B(2,1);
    B(0,0) = di1_M->xcm_ - di0_M->xcm_; 
    B(1,0) = di1_M->ycm_ - di0_M->ycm_;
    //     TMatrix AInv = TMatrix( TMatrix::kInverted, A );
    TMatrix AInv = Utils::InvertMatrix( A );
    TMatrix AInvB = TMatrix( AInv, TMatrix::kMult, B );
   
    // calculate x, y displacements
    double dx0 = di0_M->irotM_(1,0)*AInvB(0,0);
    double dy0 = di0_M->irotM_(1,1)*AInvB(0,0);
    double dx1 = di1_M->irotM_(1,0)*AInvB(1,0);
    double dy1 = di1_M->irotM_(1,1)*AInvB(1,0);
    
    // update di0 subdetectors
    for( unsigned int i=0; i<di0->sub_.size(); i++ ) {
      DetectorInfo* sub = di0->sub_[i];
      if( !sub ) continue;
      //      printf("DetFileManager::MatchCenters - \"%s\" [(%9.3f %9.3f)]",
      //	     sub->TBName_.c_str(),
      //	     sub->xcm_,
      //	     sub->ycm_ );
      sub->xcm_ += dx0;
      sub->ycm_ += dy0;
      //      printf(" -> [(%9.3f %9.3f)]\n",
      //	     sub->xcm_,
      //	     sub->ycm_ );
    } 
  
    // update di1 subdetectors
    for( unsigned int i=0; i<di1->sub_.size(); i++ ) {
      DetectorInfo* sub = di1->sub_[i];
      if( !sub ) continue;
      //printf("DetFileManager::MatchCenters - \"%s\" [(%9.3f %9.3f)]",
      //	     sub->TBName_.c_str(),
      //	     sub->xcm_,
      //	     sub->ycm_ );
      sub->xcm_ += dx1;
      sub->ycm_ += dy1;
      //printf(" -> [(%9.3f %9.3f)]\n",
      //	     sub->xcm_,
      //	     sub->ycm_ );
    }
  
    // update di0_M deadzone, if any
    DeadZoneInfo* dead0 = di0_M->dead_;
    if( dead0 ) {
      dead0->xcm_ += dx0;
      dead0->ycm_ += dy0;
    }
  
    // update di1_M deadzone, if any
    DeadZoneInfo* dead1 = di1_M->dead_;
    if( dead1 ) {
      dead1->xcm_ += dx1;
      dead1->ycm_ += dy1;
    }
    
  }
  return true;  
}
  
//______________________________________________________________________________
bool DetFileManager::MoveCenter( const char* detSelection, double dX, double dY )
{
  // Detectors
  vector< DetectorInfo* > dets = _GetDetSelection( detSelection );
  for( unsigned int i=0; i < dets.size(); i++ ) {
    DetectorInfo* det = dets[i];
    
    // move subdetectors
    for( unsigned int j=0; j< det->sub_.size(); j++ ) {
      DetectorInfo* sub = det->sub_[j];
      if( !sub ) continue;
      printf("DetFileManager::MoveCenter (det)  - \"%s\" [(%10.4f %10.4f)]",
	     sub->TBName_.c_str(),
	     sub->xcm_,
	     sub->ycm_ );
      sub->xcm_+=dX;
      sub->ycm_+=dY;
      printf(" -> [(%10.4f %10.4f)]\n",
	     sub->xcm_,
	     sub->ycm_ );
    } 
    
    // move deadzone, if any
    DeadZoneInfo *dead = det->GetMain()->dead_;
    if( !dead ) continue;
    printf("DetFileManager::MoveCenter (dead) - \"%s\" [(%10.4f %10.4f)]",
	   dead->TBName_.c_str(),
	   dead->xcm_,
	   dead->ycm_ );
    dead->xcm_+=dX;
    dead->ycm_+=dY;
    printf(" -> [(%10.4f %10.4f)]\n",
	   dead->xcm_,
	   dead->ycm_ );
    
  }
    
  return true;
}
  
//______________________________________________________________________________
bool DetFileManager::MoveDeadZone( const char* detSelection, 
				   double xOffset, 
				   double yOffset ) 
{ 
  vector< DetectorInfo* > dets = _GetDetSelection( detSelection );
  for( unsigned int idet=0; idet < dets.size(); idet++ ) {
    DetectorInfo* det  = dets[idet]->GetMain();
    DeadZoneInfo* dead = dets[idet]->GetMain()->dead_;
    if( !( det && dead ) ) continue; 
    printf( "DetFileManager::MoveDeadZone - \"%s\" [(%9.3f %9.3f),(%9.3f %9.3f)]",
	    det->TBName_.c_str(),
	    det->xcm_,
	    det->ycm_,
	    dead->xcm_,
	    dead->ycm_ );
    dead->xcm_ =   det->xcm_ + xOffset;
    dead->ycm_ =   det->ycm_ + yOffset;
    printf( " -> [(%9.3f %9.3f)]\n", dead->xcm_, dead->ycm_ );
  }     
  return true;
}  
  
//______________________________________________________________________________
bool DetFileManager::ChangePitch( const char* detSelection, double pitch, bool changeWirD )
{
  vector< DetectorInfo* > dets = _GetDetSelection( detSelection );
  for( unsigned int idet=0; idet < dets.size(); idet++ ) {
    DetectorInfo *det = dets[idet]->GetMain();
    printf( "DetFileManager::ChangePitch - \"%s\" [(%9.3f %9.3f)]",
	    det->TBName_.c_str(), 
	    det->wirP_, 
	    det->wirD_ );
    if( changeWirD ) det->wirD_ *= pitch/det->wirP_;
    det->wirP_ = pitch;
    printf( " -> [(%9.3f %9.3f)] \n", det->wirP_, det->wirD_ );
  } 
  
  return true;
}

//______________________________________________________________________________
TNode* DetFileManager::Draw3D( const char* detSelection, bool drawSub, bool drawWir ) 
{ 
  vector< DetectorInfo* > dets = _GetDetSelection( detSelection );
  if( !dets.size() ) {
    cout << "DetFileManager::Draw3D - no detectors selected.\n";
    return 0;
  }
  
  TCanvas *cv = new TCanvas("cv","Compass spectrometer",200,10,1000,500);
  
  TBRIK* mainShape = new TBRIK( "mainS", "mainS", "void", 0.1, 0.1, 0.1 );
  TNode* mainNode  = new TNode( "mainN", "mainN", mainShape );
  
  mainNode->cd();
  
  TBRIK* rotShape = new TBRIK( "rotS", "rotS", "void", 0.1, 0.1, 0.1 );
  TNode* rotNode  = new TNode( "rotN", "rotN", rotShape );
  
  double rot[] = {
    0, 1, 0, 
    0, 0, 1, 
    1, 0, 0 };
  TRotMatrix *rotM = new TRotMatrix( "rotM", "rotM", rot );
  rotNode->SetMatrix( rotM );
  rotNode->cd();
  
  double zmin = 0;
  double zmax = 0;
  
  for( unsigned int i=0; i< dets.size(); i++ ) {
    DetectorInfo *det = dets[i]->GetMain();
    if( det->zcm_ < zmin || i == 0 ) zmin = det->zcm_;
    if( det->zcm_ > zmax || i == 0 ) zmax = det->zcm_;

    Obj3D* active;
    
    // Draw Active area
    if( ( active = det->GetArea3D() ) )
      active->MakeNodes( rotNode, Utils::GetColorCode( det->TBName_.c_str() ) );
    
    // Draw Wires
    if( drawWir )
      for( int iw=0; iw<det->nWir_; iw++ ) {
	Obj3D* wire = det->GetWire3D(iw);
	if( wire ) wire->MakeNodes( rotNode, Utils::GetColorCode( det->TBName_.c_str() ) );
      }
    
    // Draw DeadZone
    DeadZoneInfo *dead = det->dead_;
    if( dead && (active = dead->GetArea3D() ) )
      active->MakeNodes( rotNode, Utils::GetColorCode( det->TBName_.c_str() ) );

    // draw subdetectors   
    if( drawSub )
      for( unsigned int j=0; j<dets[i]->sub_.size(); j++ )
	if( dets[i]->sub_[j] != det && ( active = dets[i]->sub_[j]->GetArea3D() ) ) {
	  DetectorInfo *sub = dets[i]->sub_[j];
	  active->MakeNodes( rotNode, Utils::GetColorCode( sub->TBName_.c_str() ) );
    
	  // Draw Wires
	  if( drawWir )
	    for( int iw=0; iw<sub->nWir_; iw++ ) {
	      Obj3D* wire = sub->GetWire3D(iw);
	      if( wire ) wire->MakeNodes( rotNode, Utils::GetColorCode( sub->TBName_.c_str() ) );
	    }
	}
  }
      
  mainNode->Draw();
  cv->Update();
  return mainNode;
}  
  
//========================================
// Accessing DetectorInfo and DeadZonesInfo     
//========================================

//_______________________________________________________________________________
DetectorInfo* DetFileManager::GetDetInfo( const char* selection, bool quiet )
{
  string detselection( selection );
  if( detselection.size() == 0 ) return DetInfo_.front();
  for( unsigned int iInf = 0; iInf < DetInfo_.size(); iInf++ ) {
    if( Utils::Match( DetInfo_[iInf]->TBName_, detselection ) ) return DetInfo_[iInf];
    
    // check subdetectors
    for( unsigned int iSub=0; iSub < DetInfo_[iInf]->sub_.size(); iSub++ )
      if( Utils::Match( DetInfo_[iInf]->sub_[iSub]->TBName_, detselection ) ) return DetInfo_[iInf]->sub_[iSub];
  
  }
  if( !quiet) cout << "DetFileManager::GetDetInfo - ERROR: no match for \"" << detselection << "\".\n";
  return 0;

}

//_______________________________________________________________________________
DetectorInfo* DetFileManager::GetDetInfo( int id )
{
  for( unsigned int iInf = 0; iInf < DetInfo_.size(); iInf++ ) {
    
    if( DetInfo_[iInf]->id_ == id ) return DetInfo_[iInf];
    
    // Check subdetectors
    for( unsigned int iSub=0; iSub < DetInfo_[iInf]->sub_.size(); iSub++ )
      if( DetInfo_[iInf]->sub_[iSub]->id_ == id ) return DetInfo_[iInf]->sub_[iSub];

  }

  cout << "DetFileManager::GetDetInfo - ERROR: no match for id " << id << ".\n";
  return 0;
}

//================
// Private methods
//================
    
//_______________________________________________________________________________
vector< DetectorInfo* >  DetFileManager::_GetDetSelection( string detselection )
{
 
  if( detselection.size() == 0 ) return DetInfo_;
  vector< string > select;
  istringstream in( detselection.c_str() );
  while( in ) { 
    string name; 
    in >> name; 
    if( name.size() ) select.push_back( name );
  }
  
  vector< DetectorInfo* > selection;
  for( unsigned int iInf = 0; iInf < DetInfo_.size(); iInf++ )
    for( unsigned int isel = 0; isel < select.size(); isel++ ) {
      if( Utils::Match( DetInfo_[iInf]->TBName_, select[isel] ) ) {
	selection.push_back( DetInfo_[iInf] );
	break;
      }

      // Check if any subdetector match
      for( unsigned int iSub=0; iSub < DetInfo_[iInf]->sub_.size(); iSub++ )
	if( Utils::Match( DetInfo_[iInf]->sub_[iSub]->TBName_, select[isel] ) ) {
	  selection.push_back( DetInfo_[iInf] );
	  break;
	}
      
    }   // loop oved detectors
  return selection;

}
 
//_______________________________________________________________________________
vector< DeadZoneInfo* > DetFileManager::_GetDeadSelection( string deadselection )
{
 
  if( deadselection.size() == 0 ) return DeadZInfo_;
  vector< string > select;
  istringstream in( deadselection.c_str() );
  while( in ) { 
    string name; 
    in >> name; 
    if( name.size() ) select.push_back( name );
  }
  vector< DeadZoneInfo* > selection;
  for( unsigned int iInf = 0; iInf < DeadZInfo_.size(); iInf++ )
    for( unsigned int isel = 0; isel < select.size(); isel++ ) {
      if( Utils::Match( DeadZInfo_[iInf]->TBName_, select[isel] ) )
	selection.push_back( DeadZInfo_[iInf] );
      break;
    }
  return selection;

}

//_______________________________________________
struct DetFileManager::sortDet_Line_ : public binary_function< DetectorInfo*, DetectorInfo*, bool > 
{ 
  bool operator() ( DetectorInfo* di0, DetectorInfo* di1 ) 
  { return ( di0->lineNumber_ < di1->lineNumber_ ); } 
};

//_______________________________________________
struct DetFileManager::sortDeadZ_Line_ : public binary_function< DeadZoneInfo*, DeadZoneInfo*, bool > 
{ 
  bool operator() ( DeadZoneInfo* dzi0, DeadZoneInfo* dzi1 ) 
  { return ( dzi0->lineNumber_ < dzi1->lineNumber_ ); } 
};

//_______________________________________________
void DetFileManager::_SortVects( void )
{

  // DetectorInfos
  list< DetectorInfo* > tmpDet;
  for( unsigned int i=0; i < DetInfo_.size(); i++ ) tmpDet.push_back( DetInfo_[i] );
  tmpDet.sort( DetFileManager::sortDet_Line_() );
  DetInfo_.clear();  
  for( list< DetectorInfo* >::iterator I = tmpDet.begin(); I != tmpDet.end(); I++ ) DetInfo_.push_back(*I);

  // DeadZoneInfos
  list< DeadZoneInfo* > tmpDead;
  for( unsigned int i=0; i < DeadZInfo_.size(); i++ ) tmpDead.push_back( DeadZInfo_[i] );
  tmpDead.sort( DetFileManager::sortDeadZ_Line_() );
  DeadZInfo_.clear();  
  for( list< DeadZoneInfo* >::iterator I = tmpDead.begin(); I != tmpDead.end(); I++ ) DeadZInfo_.push_back(*I);
  
  return;
}
  

  
//========================================
// Comparison with concurent detector file
//========================================

//_______________________________________________
bool DetFileManager::SameZ( const char* detFileName, const char* selection )
{
  DetFileManager DF( detFileName );
  
  // Process detfile infos
  vector< DetectorInfo* > dets  =    _GetDetSelection( selection );
  vector< DetectorInfo* > rdets = DF._GetDetSelection( selection );
  for( unsigned int i=0;      i < dets.size();   i++ )
    for( unsigned int i_s=0;    i_s  < dets[i]->sub_.size();   i_s++ ) 
      for( unsigned int iR = 0;   iR < rdets.size(); iR++ )
	for( unsigned int iR_s = 0; iR_s < rdets[iR]->sub_.size(); iR_s++ )
	  if( rdets[iR]->sub_[iR_s]->TBName_ == dets[i]->sub_[i_s]->TBName_ ) {
	    printf("DetFileManager::SameZ - Info: Det %s (%s) - %10.4f -> %10.4f.\n",
		   dets[i]->sub_[i_s]->TBName_.c_str(),
		   rdets[iR]->sub_[iR_s]->TBName_.c_str(),
		   dets[i]->sub_[i_s]->zcm_,
		   rdets[iR]->sub_[iR_s]->zcm_ );
	    dets[i]->sub_[i_s]->zcm_ = rdets[iR]->sub_[iR_s]->zcm_;
	    break;
	  }
  
  // process deadzone infos
  vector< DeadZoneInfo* > deads  =    _GetDeadSelection( selection );
  vector< DeadZoneInfo* > rdeads = DF._GetDeadSelection( selection );
  for( unsigned int i=0;      i <  deads.size();  i++ )
    for( unsigned int iR = 0;   iR < rdeads.size(); iR++ )
      if( rdeads[iR]->TBName_ == deads[i]->TBName_ ) {
	printf("DetFileManager::SameZ - Info: Dead %s (%s) - %10.4f -> %10.4f.\n",
	       deads[i]->TBName_.c_str(),
	       rdeads[iR]->TBName_.c_str(),
	       deads[i]->zcm_,
	       rdeads[iR]->zcm_ );
	deads[i]->zcm_ = rdeads[iR]->zcm_;
	break;
      }
  
  return true;
}  

//_______________________________________________
bool DetFileManager::SameSpSli( const char* detFileName, const char* selection )
{
  DetFileManager DF( detFileName );
  
  // Process detfile infos
  vector< DetectorInfo* > dets  =    _GetDetSelection( selection );
  vector< DetectorInfo* > rdets = DF._GetDetSelection( selection );
  for( unsigned int i=0;      i < dets.size();   i++ )
    for( unsigned int i_s=0;    i_s  < dets[i]->sub_.size();   i_s++ ) 
      for( unsigned int iR = 0;   iR < rdets.size(); iR++ )
	for( unsigned int iR_s = 0; iR_s < rdets[iR]->sub_.size(); iR_s++ )
	  if( rdets[iR]->sub_[iR_s]->TBName_ == dets[i]->sub_[i_s]->TBName_ ) {
	    printf("DetFileManager::SameSpSli - Info: Det %s (%s) - %10.4f -> %10.4f.\n",
		   dets[i]->sub_[i_s]->TBName_.c_str(),
		   rdets[iR]->sub_[iR_s]->TBName_.c_str(),
		   dets[i]->sub_[i_s]->spSli_,
		   rdets[iR]->sub_[iR_s]->spSli_ );
	    dets[i]->sub_[i_s]->spSli_ = rdets[iR]->sub_[iR_s]->spSli_;
	    break;
	  }
    
  return true;
}  
 
//_______________________________________________
bool DetFileManager::SameVel( const char* detFileName, const char* selection )
{
  DetFileManager DF( detFileName );
  
  // Process detfile infos
  vector< DetectorInfo* > dets  =    _GetDetSelection( selection );
  vector< DetectorInfo* > rdets = DF._GetDetSelection( selection );
  for( unsigned int i=0;      i < dets.size();   i++ )
    for( unsigned int i_s=0;    i_s  < dets[i]->sub_.size();   i_s++ ) 
      for( unsigned int iR = 0;   iR < rdets.size(); iR++ )
	for( unsigned int iR_s = 0; iR_s < rdets[iR]->sub_.size(); iR_s++ )
	  if( rdets[iR]->sub_[iR_s]->TBName_ == dets[i]->sub_[i_s]->TBName_ ) {
	    printf("DetFileManager::SameVel - Info: Det %s (%s) - %10.4f -> %10.4f.\n",
		   dets[i]->sub_[i_s]->TBName_.c_str(),
		   rdets[iR]->sub_[iR_s]->TBName_.c_str(),
		   dets[i]->sub_[i_s]->vel_,
		   rdets[iR]->sub_[iR_s]->vel_ );
	    dets[i]->sub_[i_s]->vel_ = rdets[iR]->sub_[iR_s]->vel_;
	    break;
	  }
    
  return true;
}

//_______________________________________________
bool DetFileManager::SameNames( const char* detFileName, const char* selection )
{
  DetFileManager df( detFileName );
  
  // change detector names 
  vector< DetectorInfo* > dets  =    _GetDetSelection( selection );
  vector< DetectorInfo* > rdets = df._GetDetSelection( selection );
  for( unsigned int i=0;      i < dets.size();   i++ )
    for( unsigned int iR = 0;   iR < rdets.size(); iR++ ){
    
      // Get Main detectors, compare TBNames
      DetectorInfo *rdet = rdets[iR]->GetMain();
      DetectorInfo *det  = dets[i]->GetMain();
      if( rdet->TBName_ != det->TBName_ ) continue;
    
      // Change names for all subdetectors
      for( unsigned int i_s=0;    i_s  < dets[i]->sub_.size();   i_s++ ) {
	det  = dets[i]->sub_[i_s];
	rdet = rdets[iR]->sub_[0];
	double d = pow(det->xcm_-rdet->xcm_, 2 ) + pow(det->ycm_-rdet->ycm_, 2 ); 
	for( unsigned int iR_s=0; iR_s < rdets[iR]->sub_.size(); iR_s++ ) {
        
	  // special case for straws
	  if( det->TBName_.substr(0,2) == "ST" && det->TBName_ == rdets[iR]->sub_[iR_s]->TBName_ )  { 
	    rdet = rdets[iR]->sub_[iR_s];
	    break;
          
	    // standard case
	  } else if( pow( det->xcm_ - rdets[iR]->sub_[iR_s]->xcm_, 2 ) + pow( det->ycm_ - rdets[iR]->sub_[iR_s]->ycm_, 2 ) < d ) {
	    rdet = rdets[iR]->sub_[iR_s];
	    d = pow( det->xcm_ - rdet->xcm_, 2 ) + pow( det->ycm_ - rdet->ycm_, 2 );     
	  }
	} 
	printf( "DetFileManager::SameNames - INFO: det %s (%s) - (%4i,%s,%2i) -> (%4i,%s,%2i).\n",
		det->TBName_.c_str(),
		rdet->TBName_.c_str(),
		det->id_,  det->name_.c_str(),  det->unit_,
		rdet->id_, rdet->name_.c_str(), rdet->unit_);
         
	det->id_ = rdet->id_;
	det->name_ = rdet->name_;
	det->unit_ = rdet->unit_;
      }
    
      // Change name for DeadZones
      DeadZoneInfo*  dead =  dets[i]->GetMain()->dead_;
      DeadZoneInfo* rdead = rdets[i]->GetMain()->dead_;
      if( !( dead && rdead ) ) continue;
      printf( "DetFileManager::SameNames - INFO: dead %s (%s) - (%4i,%s,%2i) -> (%4i,%s,%2i).\n",
	      dead->TBName_.c_str(),
	      rdead->TBName_.c_str(),
	      dead->id_,  dead->name_.c_str(),  dead->unit_,
	      rdead->id_, rdead->name_.c_str(), rdead->unit_);
       
      dead->id_   = rdead->id_;
      dead->name_ = rdead->name_;
      dead->unit_ = rdead->unit_;
    
    }
  return true;
} 

void DetFileManager::Streamer(TBuffer &b) {
  
  if (b.IsReading()) {
    // Version_t v = b.ReadVersion(); // unused variable (jj)
    
    int length;
    b >> length;
    //length++;

    // first check if data on detFileContents_ wasn't already created (jj)
    if(detFileContents_) delete[] detFileContents_;
    detFileContents_ = new char[length];
    b.ReadString(detFileContents_,length);
    
    // this class needs to modify the contents of a file
    // giving a tmp file
    // changed tmpnam to mkstemp (jj)
    char tmpfile[]="dfmanagerXXXXXX";
    // int IsOk = mkstemp(tmpfile); // unused variable (jj)
    if(tmpfile != NULL) {
      ofstream out(tmpfile);
      out<<detFileContents_<<endl;
      detFileName_ = tmpfile;
      out.close();
    } else {
      cerr<<"cannot create a tmp file /tmp/dfmanagerXXXXXX. exit";
      delete this;
      return;
    }

    // ok, now can follow standard build process on the tmp file
    Init();
  } else {
    b.WriteVersion(DetFileManager::IsA());
    b << strlen(detFileContents_);
    b.WriteString(detFileContents_);    
  }
}

