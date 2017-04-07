// $Id: Align.cc,v 1.22 2010/02/03 18:22:23 suhl Exp $
 
/*!
  \file    Align.cc
  \brief   Alignment Interface Class.
  \author  Hugo Pereira
  \version $Revision: 1.22 $
  \date    $Date: 2010/02/03 18:22:23 $
*/

#include "millepede.h"
#include "Align.h"
#include "Tracks.h"
#include "DetFileManager.h"
#include "DetectorInfo.h"
#include "Utils.h"

#include <iostream>
#include <sstream>
#include <fstream>
#include <list>
#include <iomanip>
char* operator+( std::streampos&, char* );
using namespace std;

//!________________________________________
ClassImp(Align)
  Align::Align(const char* trackfileselection, const char* detectorfile, bool magnets_on ):
    TObject(),
    isBatch_( false ),
    df_(0),
    tracks_(0),
    detectorFile_( "" ),
    cut_(""),
    magnets_on_( magnets_on ),

    //! Alignment switches
    alignU_( false ),
    alignZ_( false ),
    alignT_( false ),
    alignP_( false ),
    alignR_( false ),
    alignL_( false ),

    //! configuration switches
    parInit_( false ),
    iterate_( false ),
    nStdDev_( 1 ),
    nTracks_( 0 )
{
  if( strlen( detectorfile ) )       LoadDetectorFile( detectorfile );
  if( strlen( trackfileselection ) ) LoadTracks( trackfileselection, magnets_on ); 
}


//!________________________________________
DetFileManager* Align::LoadDetectorFile( const char* name )
{
  if( df_ ) SafeDelete(df_);
  df_ = new DetFileManager( name );
  detectorFile_ = string(name);
  LoadAllDetectors();
  return df_;
}

//!________________________________________
Tracks* Align::LoadTracks( const char* fileselection, bool magnets_on )
{
  if( tracks_ ) SafeDelete(tracks_);

  magnets_on_ = magnets_on;
  cout << "Align::LoadTracks - magnets: " << ((magnets_on_)?"on":"off") << ".\n";

  tracks_ = new Tracks( magnets_on_ );  
  tracks_->AddToChain( fileselection, isBatch_ );
  trackFileSelection_.clear();
  trackFileSelection_.push_back( string(fileselection) );
  return tracks_;
}

//!________________________________________
Tracks* Align::AddToTracks( const char* fileselection )
{
  if( !tracks_ ) {
    cout << "Align::AddToTracks - Tracks object is not defined.\n";
    return 0;
  }
  tracks_->AddToChain( fileselection, isBatch_ );
  trackFileSelection_.push_back( string(fileselection) );
  return tracks_;
}

//!________________________________________
bool Align::LoadAllDetectors( void )
{
  if( !df_ ) {
    cout << "Align::LoadAllDetectors - ERROR: no detector File.\n";
    return false;
  }
  
  dets_.clear();
  nDetTracks_.clear();
  dets_ = df_->GetDetInfo();
  return true;
}   

//!________________________________________
int Align::UseDetectors( const char* detselection )
{
  nDetTracks_.clear();
  if( !( df_ && dets_.size() ) ) {
    cout << "Align::UseDetectors - ERROR: no detectors loaded.\n";
    return 0;
  }
  
  vector< string > select;
  istringstream in( detselection );
  while( in ) { 
    string name; 
    in >> name; 
    if( name.size() ) select.push_back( name );
  }
  
  vector< DetectorInfo* > tmp = dets_;
  dets_.clear();
  for( unsigned int iinf=0; iinf < tmp.size(); iinf++ )
    for( unsigned int isel=0; isel < select.size(); isel++ )
      if( Utils::Match( tmp[iinf]->TBName_, select[isel] ) ) {
	dets_.push_back( tmp[iinf] );
	printf("Align::UseDetectors - Adding \"%s\".\n", tmp[iinf]->TBName_.c_str() );
	break;
      }
  parInit_=false;
  return dets_.size();
}

//!_____________________________________________________
int Align::ExcludeDetectors( const char* detselection )
{
  nDetTracks_.clear();
  if( !( df_ && dets_.size() ) ) {
    cout << "Align::ExcludeDetectors - ERROR: no detectors loaded.\n";
    return 0;
  }
  
  vector< string > select;
  istringstream in( detselection );
  while( in ) { 
    string name; 
    in >> name; 
    if( name.size() ) select.push_back( name );
  }
  
  vector< DetectorInfo* > tmp = dets_;
  dets_.clear();
  for( unsigned int iinf=0; iinf < tmp.size(); iinf++ ) {
    bool accepted = true;
    for( unsigned int isel=0; isel < select.size(); isel++ )
      if( Utils::Match( tmp[iinf]->TBName_, select[isel] ) ) {
	printf("Align::ExcludeDetectors - Removing \"%s\".\n", tmp[iinf]->TBName_.c_str() );
	accepted = false;
	break;
      }
    if( accepted ) dets_.push_back( tmp[iinf] );
  }
  parInit_=false;
  return dets_.size();
}


//!_____________________________________________________
void Align::ChangeResolution( const char* detselection, double res )
{
  if( !( df_ && dets_.size() ) ) {
    cout << "Align::ChangeResolution - ERROR: no detectors loaded.\n";
    return;
  }
  
  vector< string > select;
  istringstream in( detselection );
  while( in ) { 
    string name; 
    in >> name; 
    if( name.size() ) select.push_back( name );
  }
  
  for( unsigned int iinf=0; iinf < dets_.size(); iinf++ ) {
    for( unsigned int isel=0; isel < select.size(); isel++ )
      if( Utils::Match( dets_[iinf]->TBName_, select[isel] ) ) {
	dets_[iinf]->res_ = res;     
	printf("Align::ChangeResolution - \"%s\" -> %f mm.\n", dets_[iinf]->TBName_.c_str(), dets_[iinf]->res_ );
	break;
      }
  }
  return;
}


//!________________________________________
bool Align::SortDetectors( void )
{
  if( !( df_ && dets_.size() ) ) {
    cout << "Align::SortDetectors - ERROR: no detectors loaded.\n";
    return false;
  }
  
  list< DetectorInfo *> detList;
  for( unsigned int i=0; i<dets_.size(); i++ ) detList.push_back( dets_[i] );
  detList.sort( DetFileManager::sortDet_() );

  dets_.clear();
  for( list< DetectorInfo* >::iterator I = detList.begin(); I!=detList.end(); I++ )
    dets_.push_back( *I );
  return true;
}

//!________________________________________
bool Align::SetBias( const char* TBName, const char* biasL )
{
  bool found = false;
  for( unsigned int i=0; i<dets_.size(); i++ )
    if( dets_[i]->TBName_ == string( TBName ) ) {
      found = true;
      dets_[i]->SetBias( biasL );
      break;
    }
  
  return found;
}
  
//!________________________________________
bool Align::AddCut( const char* cut )
{
  if( !tracks_ ) {
    cout << "Align::AddCut - ERROR: no tracks loaded.\n";
    return false;
  }
  
  tracks_->AddCut( cut );
  cut_ = string( cut );
  return true;
}

//!_______________________________________________________________
void Align::DumpDetectors( const char* selection, std::ostream &out )
{
  for( unsigned int i=0; i < dets_.size(); i++ ) 
    out << *dets_[i];
}  

//!_______________________________________________________________
bool Align::InitParameters( int nStdDev, bool dumpMille )
{
  if( !( df_ && dets_.size() ) ) {
    cout << "Align::InitParameters - ERROR: no detectors loaded.\n";
    parInit_ = false;
    return false;
  }

  // cout<<"jj "<<dets_.size()<<endl;
  // for (unsigned int i=0;i<dets_.size();i++) cout<<dets_[i]->TBName_.c_str()<<endl ;  

  if( dets_.size() >= NPLAN ) {
    cout << "Align::InitParameters - ERROR: too many detectors selected (>" << NPLAN << ").\n";
    return false;
  }
  
  //! initialize millepede 
  int nGlb=dets_.size()*NPARPLAN;
  nStdDev_ = nStdDev;
  C_INITGL(nGlb,NPARTRCK,nStdDev_, int(dumpMille)-1);
  parInit_ = true;
  return true;
}
   
//!_______________________________________________________________
bool Align::FixU( const char* detselection )
{
  if( !parInit_ ){
    cout << "Align::Fix - ERROR: parameters not initialized.\n";
    return false;
  }
  
  vector< string > select;
  istringstream in( detselection );
  while( in ) { 
    string name; 
    in >> name; 
    if( name.size() ) select.push_back( name );
  }

  for( unsigned int isel=0; isel < select.size(); isel++ ) {
    bool found = false;
    for( unsigned int iinf=0; iinf < dets_.size(); iinf++ )
      if( Utils::Match( dets_[iinf]->TBName_, select[isel] ) ) {
        C_PARSIG( iinf*NPARPLAN+1,0.0 );

        // for pixel detectors use entry 5 for V
        if ( (dets_[iinf]->TBName_.substr(0,2)=="GP" ||
	      dets_[iinf]->TBName_.substr(0,2)=="MP") &&
	     dets_[iinf]->TBName_.substr(4,1)=="P" ) // PixelGEM/MMs are the only pixel detectors
          C_PARSIG( iinf*NPARPLAN+5, 0.0 );

        found = true;
      }
    if( !found ) cout << "Align::FixU - no match for \"" << select[isel] << "\".\n";
  }
  return true;
}
   
//!_______________________________________________________________
bool Align::FixZ( const char* detselection )
{
  if( !parInit_ ){
    cout << "Align::FixZ - ERROR: parameters not initialized.\n";
    return false;
  }
  
  vector< string > select;
  istringstream in( detselection );
  while( in ) { 
    string name; 
    in >> name; 
    if( name.size() ) select.push_back( name );
  }

  for( unsigned int isel=0; isel < select.size(); isel++ ) {
    bool found = false;
    for( unsigned int iinf=0; iinf < dets_.size(); iinf++ )
      if( Utils::Match( dets_[iinf]->TBName_, select[isel] ) ) {
	C_PARSIG( iinf*NPARPLAN+2,0.0 );
	found = true;
      }
    if( !found ) cout << "Align::FixZ - no match for \"" << select[isel] << "\".\n";
  }
  return true;
}
   
//!_______________________________________________________________
bool Align::FixT( const char* detselection )
{
  if( !parInit_ ){
    cout << "Align::FixT - ERROR: parameters not initialized.\n";
    return false;
  }
  
  vector< string > select;
  istringstream in( detselection );
  while( in ) { 
    string name; 
    in >> name; 
    if( name.size() ) select.push_back( name );
  }

  for( unsigned int isel=0; isel < select.size(); isel++ ) {
    bool found = false;
    for( unsigned int iinf=0; iinf < dets_.size(); iinf++ )
      if( Utils::Match( dets_[iinf]->TBName_, select[isel] ) ) {
	C_PARSIG( iinf*NPARPLAN+3,0.0 );
	found = true;
      }
    if( !found ) cout << "Align::FixT - no match for \"" << select[isel] << "\".\n";
  }
  return true;
}
  
//!_______________________________________________________________
bool Align::FixP( const char* detselection )
{
  if( !parInit_ ){
    cout << "Align::FixP - ERROR: parameters not initialized.\n";
    return false;
  }
  
  vector< string > select;
  istringstream in( detselection );
  while( in ) { 
    string name; 
    in >> name; 
    if( name.size() ) select.push_back( name );
  }

  for( unsigned int isel=0; isel < select.size(); isel++ ) {
    bool found = false;
    for( unsigned int iinf=0; iinf < dets_.size(); iinf++ )
      if( Utils::Match( dets_[iinf]->TBName_, select[isel] ) ) {
	C_PARSIG( iinf*NPARPLAN+4,0.0 );
	found = true;
      }
    if( !found ) cout << "Align::FixP - no match for \"" << select[isel] << "\".\n";
  }
  return true;
}
  
//!_______________________________________________________________
bool Align::FixR( const char* detselection )
{
  if( !parInit_ ){
    cout << "Align::FixR - ERROR: parameters not initialized.\n";
    return false;
  }
  
  vector< string > select;
  istringstream in( detselection );
  while( in ) { 
    string name; 
    in >> name; 
    if( name.size() ) select.push_back( name );
  }

  for( unsigned int isel=0; isel < select.size(); isel++ ) {
    bool found = false;
    for( unsigned int iinf=0; iinf < dets_.size(); iinf++ )
      if( Utils::Match( dets_[iinf]->TBName_, select[isel] ) ) {

        // for pixel detectors use entry 5 (R) for V, so do not set it to zero
        if ( !((dets_[iinf]->TBName_.substr(0,2)=="GP" ||
		dets_[iinf]->TBName_.substr(0,2)=="MP") &&
	       dets_[iinf]->TBName_.substr(4,1)=="P") ) // PixelGEM/MMs are the only pixel detectors
          C_PARSIG( iinf*NPARPLAN+5, 0.0 );

	    found = true;
      }
    if( !found ) cout << "Align::FixR - no match for \"" << select[isel] << "\".\n";
  }
  return true;
}
  
//!_______________________________________________________________
bool Align::FixL( const char* detselection )
{
  if( !parInit_ ){
    cout << "Align::FixL - ERROR: parameters not initialized.\n";
    return false;
  }
  
  vector< string > select;
  istringstream in( detselection );
  while( in ) { 
    string name; 
    in >> name; 
    if( name.size() ) select.push_back( name );
  }

  for( unsigned int isel=0; isel < select.size(); isel++ ) {
    bool found = false;
    for( unsigned int iinf=0; iinf < dets_.size(); iinf++ )
      if( Utils::Match( dets_[iinf]->TBName_, select[isel] ) ) {
	C_PARSIG( iinf*NPARPLAN+6,0.0 );
	found = true;
      }
    if( !found ) cout << "Align::FixL - no match for \"" << select[isel] << "\".\n";
  }
  return true;
}
    
//!_______________________________________________________________
bool Align::Minimize( unsigned int nTracks, unsigned int refresh, const char* reparamdetselection )
{
  //! Check detectors
  if( !( df_ && dets_.size() ) ) {
    cout << "Align::Minimize - ERROR: no detectors loaded.\n";
    return false;
  } else 
    printf("Align::Minimize - %zu detectors selected.\n", dets_.size() ); 
  
  //! Check parameters
  if( !( parInit_ ) ) {
    cout << "Align::Minimize - ERROR: parameters not initialized.\n";
    return false;
  }
    
  //! Check which parameters to align
  if( !( alignU_ || alignZ_ || alignT_ || alignP_ || alignR_ ) ) {
    cout << "Align::Minimize - ERROR: nothing to minimize.\n";
    return false;
  } else {
    cout << "Align::Minimize - AlignU: " << ((alignU_)?"yes":"no") << endl; 
    cout << "Align::Minimize - AlignZ: " << ((alignZ_)?"yes":"no") << endl; 
    cout << "Align::Minimize - AlignT: " << ((alignT_)?"yes":"no") << endl; 
    cout << "Align::Minimize - AlignP: " << ((alignP_)?"yes":"no") << endl; 
    cout << "Align::Minimize - AlignR: " << ((alignR_)?"yes - drift-like detectors only (type==11)":"no") << endl; 
    cout << "Align::Minimize - AlignL: " << ((alignL_)?"yes - drift-like detectors only (type==11)":"no") << endl; 
  }
  
  //! Check tracks
  if( ! (tracks_ && tracks_->GetEntries() ) ) {
    cout << "Align::Minimize - ERROR: no track loaded.\n";
    return false;
  } else if( !nTracks ) {
    cout << "Align::Minimize - All tracks required.\n";
    nTracks = tracks_->GetEntries();
  } else {
    if( nTracks > tracks_->GetEntries() ) nTracks = tracks_->GetEntries();
    cout << "Align::Minimize - " << nTracks << " tracks required.\n";
  }
  
  nTracks_ = nTracks;

  //! fix all parameters if alignment not required
  for( unsigned int i=0; i<dets_.size(); i++) {
    // for pixel detectors use entry 5 (R) for V
    bool isPixel=false;
    if ( (dets_[i]->TBName_.substr(0,2)=="GP" ||
	  dets_[i]->TBName_.substr(0,2)=="MP") &&
	 dets_[i]->TBName_.substr(4,1)=="P" ) // PixelGEM/MMs are the only pixel detectors
      isPixel=true;

    if (!alignU_) C_PARSIG(i*NPARPLAN+1,0.0)	        //!< fix all u
    if (!alignU_ && isPixel) C_PARSIG(i*NPARPLAN+5,0.0)	//!< fix all v
    if (!alignZ_) C_PARSIG(i*NPARPLAN+2,0.0)	        //!< fix all z 
    if (!alignT_) C_PARSIG(i*NPARPLAN+3,0.0)	        //!< fix all theta 
    if (!alignP_) C_PARSIG(i*NPARPLAN+4,0.0)	        //!< fix all pitch 
    if (!(alignR_ && dets_[i]->type_ == 11 ) && !isPixel) C_PARSIG(i*NPARPLAN+5,0.0)	//!< fix all R0 
    if (!(alignL_ && dets_[i]->type_ == 11 )            ) C_PARSIG(i*NPARPLAN+6,0.0)	//!< fix all lorentz angle scaling 
  }
  
  //! Initialise nTracks
  nDetTracks_.clear();
  for( unsigned int i=0; i<dets_.size(); i++ ) nDetTracks_.push_back( 0 );
  
  //! Start iterations
  zerloc_(dergb,derlc); 
  if (iterate_) C_INITUN(11,10000.);
  for ( unsigned int itrack = 0; itrack < nTracks_;) {	
    
    if( ! tracks_->GetNextEntry( ) ) {
      cout << "Align::Minimize - Cannot load entry " << itrack << ".\n";
      break;
    }
    
    //! Check number of entries
    if( !tracks_->AcceptEntry( dets_, NPARTRCK ) ) continue;
    
    //! Dump if required
    if( refresh && (itrack%refresh)==0 )
      cout << "Event: " << itrack << endl;

    //! Get Tracks static data into local vars
    int nFiredDet = (int) tracks_->T_nDets;
    for (int id=0;id< nFiredDet;id++) {		// loop on det

      //! look for detector in det list
      bool found = false;
      unsigned int iP = 0;
      for( iP=0; iP<dets_.size() && !found; iP++ ) { 
        // make a copy of _main_ detector
        DetectorInfo *dMain = dets_[iP]->GetMain();
      
        // scan subdetectors
        for( unsigned int j=0; j<dets_[iP]->sub_.size() && !found; j++ ) {
          DetectorInfo *dSub = dets_[iP]->sub_[j];
          
          // for straws, only tracks passing through the 'main' (central) detector are accepted
          // Right now true only if option "align OuterST YES" is not present (jj)
          if( (dMain->TBName_.substr(0,2) != "ST" || dSub->TBName_[7]=='b' || alignOuterST_ ) && dSub->id_== tracks_->T_detVect[id]) found = true;

        }
      }
      iP--;
      
      if( !found ) continue;
      DetectorInfo *det  = dets_[iP];
      if( iP < nDetTracks_.size() ) nDetTracks_[iP]++; // update number of tracks passing through detector iP
      
      if( !magnets_on_ ) {
        //! first case - No magnetic field, straight tracks
        //! Track coordinates at the detector
        float x = float( tracks_->T_xLx ); float tx = float( tracks_->T_txLx );
        float y = float( tracks_->T_yLx ); float ty = float( tracks_->T_tyLx );
        float z = float( tracks_->T_zLx );

        float r_ = tracks_->T_rVect[id];
        float u_ = (1.0+det->biasP_)*tracks_->T_uVect[id] 
                   + det->biasU_ 
                   + float(Utils::GetSign(r_))*det->biasR_
                   + r_*det->biasL_; 
	
        float z_ = det->zcm_ + det->biasZ_;  
        float dT = det->biasT_*PI/180.0;

        float xTrkDet = x+tx*(z_-z);
        float yTrkDet = y+ty*(z_-z);

        float cosT_  = det->cosTheta_*cos(dT) - det->sinTheta_*sin(dT);
        float sinT_  = det->sinTheta_*cos(dT) + det->cosTheta_*sin(dT);
        float sigma_ = det->res_;
                  
     	//***********************
    	//! Z alignment "move full station")
	    int samedet=1;
    	vector< string > select;
	    istringstream in( reparamdetselection );
    	while( in ) { 
	      string name; 
    	  in >> name; 
	      if( name.size() ) select.push_back( name );
    	}
    	int reparam = 0;
	    for( unsigned int isel=0; isel < select.size(); isel++ ) {
    	  if((det->TBName_.substr(0,4))== select[isel]) {
	        reparam = reparam+1;
    	  }
	    }    
    	if(reparam>0) {
  
	      for(samedet=1;iP-samedet<1000;samedet++){
    	    DetectorInfo *detold  = dets_[iP-samedet];
	        if((!(detold->TBName_.substr(0,4) == det->TBName_.substr(0,4)))/*||(det->zcm_<255000)*/) break;
    	  }
	    }
    	//*********************   

        // if the detector is a pixel detector some things are changed
        if ( (det->TBName_.substr(0,2)=="GP" ||
	      det->TBName_.substr(0,2)=="MP") &&
	     det->TBName_.substr(4,1)=="P" ) { // PixelGEMs are the only pixel detectors
          float resU =  xTrkDet*cosT_+yTrkDet*sinT_ - tracks_->T_uVect[id];
          float resV = -xTrkDet*sinT_+yTrkDet*cosT_ - tracks_->T_vVect[id];
          float res  = sqrt(resU*resU + resV*resV);

          u_     = ((1.0+det->biasP_)*tracks_->T_uVect[id] + det->biasU_) * ((1.0+det->biasP_)*tracks_->T_uVect[id] + det->biasU_);
 	      u_    += ((1.0+det->biasP_)*tracks_->T_vVect[id] + det->biasV_) * ((1.0+det->biasP_)*tracks_->T_vVect[id] + det->biasV_);
          u_     = sqrt(u_);
          sigma_ = sigma_;

          //! local derivatives
          derlc[0] = (resU*cosT_        - resV*sinT_       ) / res;
          derlc[1] = (resU*cosT_*(z_-z) - resV*sinT_*(z_-z)) / res;
          if (NPARTRCK>2) {
            derlc[2] = (resV*cosT_        + resU*sinT_       ) / res;
            derlc[3] = (resV*cosT_*(z_-z) + resU*sinT_*(z_-z)) / res;
          }

          //! global derivatives 
          dergb[NPARPLAN*iP  ]= -resU/res;                                                                       // u
          dergb[NPARPLAN*iP+1]= (resU*(tx     *cosT_+ty     *sinT_) + resV*(     ty*cosT_-tx*sinT_))      / res; // z
          dergb[NPARPLAN*iP+2]= (resU*(yTrkDet*cosT_-xTrkDet*sinT_) - resV*(xTrkDet*cosT_+yTrkDet*sinT_)) / res; // theta 
          dergb[NPARPLAN*iP+3]= (resU*(xTrkDet*cosT_+yTrkDet*sinT_) + resV*(yTrkDet*cosT_-xTrkDet*sinT_)) / res; // pitch
          dergb[NPARPLAN*iP+4]= -resV/res;                                                                       // v
          dergb[NPARPLAN*iP+5]= 0.;                                                                              // ...
       } else {
          //! store std local derivatives
          derlc[0]= cosT_;
          derlc[1]= cosT_*(z_-z);   
        
          if (NPARTRCK>2) { 
            derlc[2]= sinT_;        
            derlc[3]= sinT_*(z_-z); 
          }

         //! calculate/store global derivatives   
          dergb[NPARPLAN*iP]=-1.;                                                         //!< /d u
          //dergb[NPARPLAN*iP+1]=  cosT_*tx+sinT_*ty;                                       //!< /d z
          dergb[NPARPLAN*(iP-samedet+1)+1]= cosT_*tx+sinT_*ty;                            //!< /d z
          dergb[NPARPLAN*iP+2]= -sinT_*( x+tx*(z_-z) ) + cosT_*( y+ty*(z_-z) );           //!< /d theta
          dergb[NPARPLAN*iP+3]=  cosT_*( x+tx*(z_-z) ) + sinT_*( y+ty*(z_-z) );           //!< /d Pitch
          dergb[NPARPLAN*iP+4]= (det->type_ == 11 ) ? -float( Utils::GetSign( r_ ) ) : 0; //!< /d R0
          dergb[NPARPLAN*iP+5]= (det->type_ == 11 ) ? -r_ : 0;                            //!< /d lorentz angle scaling
        }
        
        equloc_(dergb,derlc,&u_,&sigma_);  //!< book local/global derivatives, measurement, error   

      } else {
        //! second case - magnetic field, bended tracks (not fitted. use staight correction to track instead)
        //! Track coordinates at the detector
        float x = float( tracks_->T_xVect[id] ); float tx = float( tracks_->T_txVect[id] );
        float y = float( tracks_->T_yVect[id] ); float ty = float( tracks_->T_tyVect[id] );
        float z = float( tracks_->T_zVect[id] );

        float r_  = tracks_->T_rVect[id];
        float du_ = (1.0+det->biasP_)*tracks_->T_duVect[id] 
                    + det->biasU_ 
                    + float(Utils::GetSign(r_))*det->biasR_
                    + r_*det->biasL_;

    	float z_  = det->zcm_ + det->biasZ_;  
        float dT  = det->biasT_*PI/180.0;
        
        float xTrkDet = x+tx*(z_-z);
        float yTrkDet = y+ty*(z_-z);

        float cosT_  = det->cosTheta_*cos(dT) - det->sinTheta_*sin(dT);
        float sinT_  = det->sinTheta_*cos(dT) + det->cosTheta_*sin(dT);
        float sigma_ = det->res_;
  	
    	//*************************************
	    int samedet=1;
    	vector< string > select;
	    istringstream in( reparamdetselection );
    	while( in ) { 
	      string name; 
    	  in >> name; 
	      if( name.size() ) select.push_back( name );
    	}
	    int reparam = 0;
    	for( unsigned int isel=0; isel < select.size(); isel++ ) {
	      if((det->TBName_.substr(0,4))== select[isel]) {
	        reparam = reparam+1;
    	  }
	    }    
    	if(reparam>0) {	
	      for(samedet=1;iP-samedet<1000;samedet++){
	        DetectorInfo *detold  = dets_[iP-samedet];
    	    if((!(detold->TBName_.substr(0,4) == det->TBName_.substr(0,4))) /*||(det->zcm_<255000)*/ ) break;
    	  }
	    }
    	//*******************************
	      
        // if the detector is a pixel detector some things are changed
        if ( (det->TBName_.substr(0,2)=="GP" ||
	      det->TBName_.substr(0,2)=="MP") &&
	     det->TBName_.substr(4,1)=="P" ) { // PixelGEMs are the only pixel detectors
          float resU = tracks_->T_duVect[id];
          float resV = tracks_->T_dvVect[id];
          float res  = sqrt(resU*resU + resV*resV);

          du_    = ((1.0+det->biasP_)*resU + det->biasU_) * ((1.0+det->biasP_)*resU + det->biasU_);
 	      du_   += ((1.0+det->biasP_)*resV + det->biasV_) * ((1.0+det->biasP_)*resV + det->biasV_);
          du_    = sqrt(du_);
          sigma_ = sigma_;

          //! local derivatives
          derlc[0] = (resU*cosT_        - resV*sinT_       ) / res;
          derlc[1] = (resU*cosT_*(z_-z) - resV*sinT_*(z_-z)) / res;
          if (NPARTRCK>2) {
            derlc[2] = (resV*cosT_        + resU*sinT_       ) / res;
            derlc[3] = (resV*cosT_*(z_-z) + resU*sinT_*(z_-z)) / res;
          }

          //! global derivatives 
          dergb[NPARPLAN*iP  ]= -resU/res;                                                                       // u
          dergb[NPARPLAN*iP+1]= (resU*(tx     *cosT_+ty     *sinT_) + resV*(     ty*cosT_-tx*sinT_))      / res; // z
          dergb[NPARPLAN*iP+2]= (resU*(yTrkDet*cosT_-xTrkDet*sinT_) - resV*(xTrkDet*cosT_+yTrkDet*sinT_)) / res; // theta 
          dergb[NPARPLAN*iP+3]= (resU*(xTrkDet*cosT_+yTrkDet*sinT_) + resV*(yTrkDet*cosT_-xTrkDet*sinT_)) / res; // pitch
          dergb[NPARPLAN*iP+4]= -resV/res;                                                                       // v
          dergb[NPARPLAN*iP+5]= 0.;                                                                              // ...
        } else {
          //! local derivatives
          derlc[0]= cosT_;
          derlc[1]= cosT_*(z_-z);
          if (NPARTRCK>2) {
            derlc[2]= sinT_;
            derlc[3]= sinT_*(z_-z);
          }

          //! global derivatives 
          dergb[NPARPLAN*iP]=-1.;                                                         //!< /d u
          //        dergb[NPARPLAN*iP+1]=  cosT_*tx+sinT_*ty;                                       //!< /d z
          dergb[NPARPLAN*(iP-samedet+1)+1]=  cosT_*tx+sinT_*ty;                                       //!< /d z
          dergb[NPARPLAN*iP+2]= -sinT_*( x+tx*(z_-z) ) + cosT_*( y+ty*(z_-z) );           //!< /d theta z_-z ==0!!!
          dergb[NPARPLAN*iP+3]=  cosT_*( x+tx*(z_-z) ) + sinT_*( y+ty*(z_-z) );           //!< /d Pitch
          dergb[NPARPLAN*iP+4]= (det->type_ == 11 ) ? -float( Utils::GetSign( r_ ) ) : 0; //!< /d R0
          dergb[NPARPLAN*iP+5]= (det->type_ == 11 ) ? -r_ : 0;                            //!< /d lorentz angle scaling
        }
        
        equloc_(dergb,derlc,&du_,&sigma_);  //!< book local/global derivatives, measurement, error 

      }
    
    }  //!< loop over dets
    
    fitloc_(); //!< Minimize track parameters (local minimisation)
    itrack++;
  }
  
  fitglo_(par); //!< minimize alignment parameters 
  C_PRTGLO(20); //!< Dump to screen
  return true;
}    

//!________________________________________
bool Align::DumpToFile( const char* filename, const char* reparamdetselection ) 
//bool Align::DumpToFile( string filename )
{
  Utils::MakeBackup( filename );
 
  //  std::ofstream out2( filename, ios::out );
  ofstream out( filename );
  if( !out ) {
    cout << "Align::DumpToFile - ERROR: cannot write to file \"" << filename << "\".\n";
    return false;
  }
  _Dump( out, reparamdetselection );

  out.close();
  return true;
}

//!________________________________________
void Align::_Dump( std::ostream &out, const char* reparamdetselection )
{  
  if (NPARPLAN>=4) {
    // printf("%-12s %-22s %-22s %-22s %-22s %-22s\n", 
    // "TBName",
    // "U(mm)", "Z(mm)", "T(deg)", "P",/* "R0(mm)", "L",*/
    // "nTracks");
    out << "TBName   " << setfill (' ') << setw (12) << left << "U(mm)" << setfill (' ') << setw (15) << left << "dU(mm)" << setfill (' ') << setw (12) << left << "Z(mm)" << setfill (' ') << setw (15) << left << "dZ(mm)" << setfill (' ') << setw (12) << left << "T(deg)" << setfill (' ') << setw (15) << left << "dT(deg)" << setfill (' ') << setw (12) << left << "P" << setfill (' ') << setw (15) << left << "dP" << setfill (' ') << setw (7) << left << "nTracks" << endl;; 

    bool proj2nd=false;
    for (unsigned int i=0;i<dets_.size();i++) {
      
      //! Dump TBName
      // printf("%8s ", dets_[i]->TBName_.c_str() );
      if ( (dets_[i]->TBName_.substr(0,2)=="GP" ||
	    dets_[i]->TBName_.substr(0,2)=="MP") &&
	   (dets_[i]->TBName_.substr(4,1)=="P" ||
	    dets_[i]->TBName_.substr(4,1)=="M") ) {
        if (proj2nd) 
          out << (dets_[i]->TBName_).substr(0,7).c_str() << "V " ;
        else
          out << (dets_[i]->TBName_).substr(0,7).c_str() << "U " ;
      } else
        out << dets_[i]->TBName_.c_str() << " " ;

      //! Dump U
      int j;
      if (proj2nd) j=NPARPLAN*i+5;
      else         j=NPARPLAN*i+1;
      float err=errpar_(&j);
      if( err==0 && alignU_ ) err=-999.9;
      else if( alignU_ && proj2nd ) dets_[i]->biasV_ += par[NPARPLAN*i+4];
      else if( alignU_ ) dets_[i]->biasU_ += par[NPARPLAN*i];
      // printf("%10.4f %10.4f ", par[NPARPLAN*i],err);
      out << setfill (' ') << setw (12) << setprecision (5) << left << showpoint << showpos << par[NPARPLAN*i+(proj2nd?4:0)] << setfill (' ') << setw (15) << setprecision (5) << left << err;
      //******************************
	  //! Dump Z

	  int samedet=1;
      vector< string > select;
      istringstream in( reparamdetselection );
      while( in ) { 
	string name; 
	in >> name; 
	if( name.size() ) select.push_back( name );
      }
      int reparam = 0;
      for( unsigned int isel=0; isel < select.size(); isel++ ) {
	if((dets_[i]->TBName_.substr(0,4))== select[isel]) {
	  reparam = reparam+1;
	}
      }    
      if(reparam>0) {
	for(samedet=1;i-samedet<1000;samedet++){
	  if((!(dets_[i]->TBName_.substr(0,4) == dets_[i-samedet]->TBName_.substr(0,4)))/*||(dets_[i]->zcm_<255000)*/) break;
	}
      }
       

      // if DV coordinate is aligned, z/pitch/angle are not changed !! 
      // ^ not true for the moment, let's see what happens (aa)
      
      
      j=NPARPLAN*(i-samedet+1)+2;
      err=errpar_(&j);
      if( (err==0 && alignZ_) /*|| (alignDV_ && (dets_[i]->TBName_.substr(0,2)=="GP" || dets_[i]->TBName_.substr(0,2)=="MP") && dets_[i]->TBName_.substr(4,1)=="P")*/) err=-999.9;
      else if( alignZ_ ) dets_[i]->biasZ_ += par[NPARPLAN*i-(samedet+1)+1];
      // out.form("%10.4f %10.4f ", par[NPARPLAN*(i-samedet+1)+1],err);
      out << setfill (' ') << setw (12) << setprecision (5) << left << showpoint <<  par[NPARPLAN*(i-samedet+1)+1] << setfill (' ') << setw (15) << setprecision (5) << left << err;
      //******************************************      
      
      //! Dump Z
      //    j=NPARPLAN*i+2;
      // err=errpar_(&j);
      // if( err==0 && alignZ_ ) err=-999.9;
      // else if( alignZ_ ) dets_[i]->biasZ_ += par[NPARPLAN*i+1];
      // printf("%10.4f %10.4f ", par[NPARPLAN*i+1],err);
      // out << setfill (' ') << setw (12) << setprecision (5) << left << showpoint <<  par[NPARPLAN*i+1] << setfill (' ') << setw (15) << setprecision (5) << left << err;
      
      
      //! Dump T
      j=NPARPLAN*i+3;
      err=errpar_(&j);
      if ( (err==0 && alignT_) /*|| (alignDV_ && (dets_[i]->TBName_.substr(0,2)=="GP" || dets_[i]->TBName_.substr(0,2)=="MP") && dets_[i]->TBName_.substr(4,1)=="P")*/) err=-999.9;
      else if( alignT_ ){ dets_[i]->biasT_ += par[NPARPLAN*i+2];
	err*=180.0/PI;}
      // printf("%10.4f %10.4f  ", par[NPARPLAN*i+2]*180.0/PI,err);
      out << setfill (' ') << setw (12) << setprecision (5) << left << showpoint <<  par[NPARPLAN*i+2]*180.0/PI <<  setfill (' ') << setw (15) << setprecision (5) << left << err; 

      //! Dump P
      j=NPARPLAN*i+4;
      err=errpar_(&j);
      if( (err==0 && alignP_) /* || (alignDV_ && (dets_[i]->TBName_.substr(0,2)=="GP" || dets_[i]->TBName_.substr(0,2)=="MP") && dets_[i]->TBName_.substr(4,1)=="P")*/) err=-999.9;
      else if( alignP_ ) dets_[i]->biasP_ += par[NPARPLAN*i+3];
      // printf("%10.5f %10.5f  ", par[NPARPLAN*i+3],err);
      out << setfill (' ') << setw (12) << setprecision (5) << left << showpoint << par[NPARPLAN*i+3] <<  setfill (' ') << setw (15) << setprecision (5) << left << err;

      //! Dump R
      j=NPARPLAN*i+5;
      err=errpar_(&j);
      if( err==0 && alignR_ ) err=-999.9;
      else if( alignR_ ) dets_[i]->biasR_ += par[NPARPLAN*i+4];
      // printf("%10.4f %10.4f  ", par[NPARPLAN*i+4],err);
      //      out << par[NPARPLAN*i+4] << " " << err << "  ";

      //! Dump L
      j=NPARPLAN*i+6;
      err=errpar_(&j);
      if( err==0 && alignL_ ) err=-999.9;
      else if( alignL_ ) dets_[i]->biasL_ += par[NPARPLAN*i+5];
      //      printf("%10.4f %10.4f  ", par[NPARPLAN*i+5],err);
      //      out << par[NPARPLAN*i+5] << " " << err << "  ";

      //! Dump number of tracks 
      if( i < nDetTracks_.size() ) 
	//printf( "%10i ", nDetTracks_[i] );
	out << setfill (' ') << setw (7) << left << noshowpos << nDetTracks_[i];
      // std::cout.form( "%10i ", nDetTracks_[i] );
      
      out << endl ;

      if (proj2nd) {
        proj2nd = false;
      } else if ( (dets_[i]->TBName_.substr(0,2)=="GP" ||
		   dets_[i]->TBName_.substr(0,2)=="MP") &&
		  dets_[i]->TBName_.substr(4,1)=="P" ) {
        proj2nd = true;
        i--;
      }
    }
  }
  
  out << endl;
  out << "Magnet        :" << ((magnets_on_)?"on":"off")  << endl;
  out << "Align U       :" << ((alignU_)?"yes":"no")  << endl;
  out << "Align Z       :" << ((alignZ_)?"yes":"no")  << endl;
  out << "Align angle   :" << ((alignT_)?"yes":"no")  << endl;
  out << "Align pitch   :" << ((alignP_)?"yes":"no")  << endl;
  out << "Align R0      :" << ((alignR_)?"yes":"no")  << endl;
  out << "Align L       :" << ((alignL_)?"yes":"no")  << endl;
  out << "Iterations    :" << ((iterate_)?"yes":"no") << endl;
  out << "nStdDev       :" << nStdDev_ << endl;
  out << "nTracks       :" << nTracks_ << endl;
  out << "Selection cut :" << (( cut_.size() ) ? cut_.c_str():"none") << endl;
  out << "detector table:" << detectorFile_ << endl; 
  out << endl;
  
  for( unsigned int i = 0; i < trackFileSelection_.size(); i ++ ) 
    out << "input tracks   :" << trackFileSelection_[i] << endl;
  out << endl;
  return;
}  
