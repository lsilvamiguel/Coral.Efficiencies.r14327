// $Id: UpdateFromAlign.cc,v 1.14 2009/08/31 01:38:57 ybedfer Exp $
/*!
  \file    UpdateFromAlign.cc
  \brief   DetFileManagerFunction to update detector file from alignment output
  \author  Hugo Pereira
  \version $Revision: 1.14 $
  \date    $Date: 2009/08/31 01:38:57 $
*/

#include "DetFileManager.h"
#include <iostream>
#include <fstream>
#include <sstream>

char* operator+( std::streampos&, char* );
using namespace std;

//_______________________________________________________________________________
bool DetFileManager::UpdateFromAlign( const char* file )
{
  ifstream in( file, ios::in );
  if( !in ) {
    cout << "DetFileManager::UpdateFromAlign - ERROR: cannot read file \"" 
         << file  << "\"." 
         << endl;
    return false;
  }
  
     
  char* line = new char[512];
  unsigned int lineNumber = 0;
  
  while( !in.eof() ) { 
    in.getline( line, 512, '\n');
    if( in.eof() || !in.good() ) continue;
    lineNumber++;      
    
    if( line[0] == '\0' ) { continue; }
    
    istringstream s(line);

    // Read TBName
    string TBName; s >> TBName;
   
    // test the line validity
    if( !s.good() ) continue;
    
    // avoid reading comments (jj)
    if(TBName=="TBName") continue;
    if(TBName=="Options") break;

    // special treatment of pixel (2 projection) detectors
    bool isPixel=false;
    bool projU=true;
    if ((TBName.substr(0,2)=="GP" ||
	 TBName.substr(0,2)=="MP") &&
	(TBName.substr(4,1)=="P" ||
	 TBName.substr(4,1)=="M") ) { // PixelGM/MM
      isPixel=true;
      if (TBName[7]=='U') projU=true;   // 8th character in name defines the projection which is to be changed
      else                projU=false;  // U proj. is used for U, pitch, angle, Z
      TBName[7]='_';                 // V proj. is only used for V
    }

    // Get and Check DetectorInfos
    vector< DetectorInfo* > dets = _GetDetSelection( TBName ); 
    if( !dets.size() ) {
      cout << "DetFileManager::UpdateFromAlign - WARNING: no match for \"" << TBName << "\". Skipping this line...\n";
      continue;
    }
    
    if( dets.size() > 1 ) {
      cout << "DetFileManager::UpdateFromAlign - ERROR: more than 1 match for \"" << TBName << "\".\n";
      continue;
    }
    
    DetectorInfo *det = dets.front();
    
    // Put all strings into a vector
    vector< string > values;
    string buf;
    while( s.good() ) { 
      s >> buf;
      if( buf.size() ) values.push_back( buf );
    }
    
    // Parse valid lines
    bool foundU = false;
    bool foundZ = false;
    bool foundT = false;
    bool foundP = false;
    
    double du=0,     du_err = -999.9;          
    double dz=0,     dz_err = -999.9;          
    double dteta=0,  dteta_err = -999.9;   
    double dpitch=0, dpitch_err = -999.9; 
    for( unsigned int i=0; i < values.size(); i++ ) {
      istringstream val( values[i].c_str() );
      switch( i ) {
      case 0: val >> du;     foundU = true; break; //!< [mm]
      case 2: val >> dz;     foundZ = true; break; //!< [mm] 
      case 4: val >> dteta;  foundT = true; break; //!< [mm]
      case 6: val >> dpitch; foundP = true; break; //!< [deg]
      case 1: val >> du_err;     break;
      case 3: val >> dz_err;     break;
      case 5: val >> dteta_err;  break;   
      case 7: val >> dpitch_err; break;   
      default: break;
      }
    }

    if (isPixel && !projU) { // second pixel projection, only use for V alignment
      foundZ = false;
      foundT = false;
      foundP = false;
    }
    
    // store detectorInfo 'main' detector
    DetectorInfo* dMain = det->GetMain();
    
    // Update dead zone, if any
    DeadZoneInfo *dead = dMain->dead_;
    if( dead ) {
      
      // Update center Position (mm)
      if( foundU && du && du_err != -999.9 ) {
        dead->xcm_ += dMain->cosTheta_*du;
        dead->ycm_ += dMain->sinTheta_*du;
      }
            
      // Update Z
      if( foundZ && dz && dz_err != -999.9 ) dead->zcm_+= dz;
      
      // update angle. 
      // WARNING: only the center position is updated. 
      // The angle given by rotation matrix remains unchanged, as it have no meaning in case of circular dead zone
      if( foundT && dteta && dteta_err != -999.9 ) {
        double dtRad = dteta*TMath::Pi()/180;
        double x = cos(dtRad)*dead->xcm_ - sin(dtRad)*dead->ycm_; 
        double y = sin(dtRad)*dead->xcm_ + cos(dtRad)*dead->ycm_;
        dead->xcm_ = x;
        dead->ycm_ = y;
      }
      
      // Update pitch (no unit)        
      // WARNING: only the center position is updated. The deadzone size remains unchanged
      if( foundP && dpitch && dpitch_err != -999.9 ) {
        dead->xcm_ *= (1-dpitch );
        dead->ycm_ *= (1-dpitch );
      }
    }
    
    // Update det, and det subdetectors, if any
    for( unsigned int isub=0; isub < det->sub_.size(); isub++ ) { 
      DetectorInfo* dSub = det->sub_[isub];
      cout << dSub->TBName_ << " -> ";

    
      // Update center Position (mm)
      if( foundU && du && du_err != -999.9 ) {
        if ( isPixel ) {
          if ( dSub==dMain ) {
            if ( projU ) { dSub->DUtoCenter(  du ); cout << "U"; }
            else         { dSub->DVtoCenter( -du ); cout << "V"; }
          }  else {
            if ( projU ) { dSub->DVtoCenter(  du ); cout << "V"; }
            else         { dSub->DUtoCenter(  du ); cout << "U"; }
          }
	    } else {
    	  dSub->DUtoCenter( du );
	      cout << "U";         
    	}
      }

	      
      // Update Z, position along the beam (mm)
      if( foundZ && dz && dz_err != -999.9 ) {
        dSub->zcm_ += dz;
        cout << "Z";
      }
 
      
      // Update angle (deg)
      if( foundT && dteta && dteta_err != -999.9 ) {
        // update angles, rotation matrices and center position
        dSub->ang_ += dteta;
        dSub->UpdateRotMatrices();
        double dtRad = dteta*TMath::Pi()/180;
        double x = cos(dtRad)*dSub->xcm_ - sin(dtRad)*dSub->ycm_; 
        double y = sin(dtRad)*dSub->xcm_ + cos(dtRad)*dSub->ycm_;
        dSub->xcm_ = x;
        dSub->ycm_ = y;
        cout << "T";
      }
      
      // Update pitch (no unit)         
      // For Straws, only main detector is updated
      // above statement is not true right now (jj)
      if(foundP && dpitch && dpitch_err != -999.9 /* && (det->TBName_.substr(0,2) != "ST" || dSub == dMain ) */ ){
        dSub->wirP_ *= (1-dpitch );
        dSub->wirD_ *= (1-dpitch );
        dSub->xcm_ *= (1-dpitch );
        dSub->ycm_ *= (1-dpitch );
        cout << "P";
      }
      cout << "\n"; 
    }
  
  }   // loop over lines in file
    
  delete[] line;
  return true;
}  
