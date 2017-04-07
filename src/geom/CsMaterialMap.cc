// $Id: CsMaterialMap.cc 14069 2015-09-17 20:44:46Z lsilva $

/*!
   \file CsMaterialMap.h
   \brief Material Map Class.
   \authors  Alexandre Korzenev, Jan P. Nassalski
   \version  $Revision: 14069 $
   \date     $Date: 2015-09-17 22:44:46 +0200 (Thu, 17 Sep 2015) $
*/

# include "CsOpt.h"
# include <string.h>
# include <stdlib.h>
# include <stdio.h>
# include "CsMaterialMap.h"

#include "TROOT.h"
#include "TMacro.h"
#include "TGeoManager.h"

#include "coral_config.h"

# define DEG2RAD 0.017453292519943295769    //   Pi/180

using namespace std;

CsMaterialMap::MatMap::MatMap(){

  CoordSyst = 0;   
  nx = ny = nz = 0; 
  minX = minY = minZ = 0;
  maxX = maxY = maxZ = 0;
  x.clear(); y.clear(); z.clear();
  Cos.clear(); Sin.clear();
  Val.clear();

  return;
}



inline double CsMaterialMap::MatMap::StepToBorder( double Z, bool direc ) {
  double step;
  if( CoordSyst > 0 ){
    if     ( Z >= maxZ ) step = direc ?     9999999    : maxZ - Z ;
    else if( Z <= minZ ) step = direc ? minZ - Z :     -9999999    ;
    else                 step = direc ? maxZ - Z : minZ - Z ;}
  else{
    if     ( Z >= maxY ) step = direc ?     9999999    : maxY - Z ;
    else if( Z <= minY ) step = direc ? minY - Z :     -9999999    ;
    else                 step = direc ? maxY - Z : minY - Z ;
  }
  
  return step;
}; 



inline bool CsMaterialMap::MatMap::insideZ( double Z ) { 
  if( CoordSyst > 0 ){
    if( Z>minZ && Z<maxZ ) return true; 
    else                   return false;}
  else{
    if( Z>minY && Z<maxY ) return true; 
    else                   return false;}
}



CsMaterialMap::CsMaterialMap() : 
  ROOTGeometry_ (NULL)
{
  
#if USE_TGEANT
  usingROOTGeometryTGEANT_ = false;
#endif
  
  string path;
  if( CsOpt::Instance()->getOpt( "CsROOTGeometry", "file", path ) ){
    CsErrLog::Instance()->mes( elInfo, "ROOT geometry will be used" );

    ROOTGeometry_  = new TMacro("detectors");
    if (!ROOTGeometry_->ReadFile(path.c_str())
	|| (ROOTGeometry_->Exec() && !gGeoManager))
      {
	CsErrLog::Instance()->mes( elFatal, "ROOT geometry couldn't be loaded" );
	abort();
      }

    usingROOTGeometry_ = true;

    massDefault_ = 0.139579018;  // Default to Pion mass
    CsOpt::Instance()->getOpt( "CsROOTGeometry", "massDefault", massDefault_ );

    simpleELoss_ = 1;
    CsOpt::Instance()->getOpt( "CsROOTGeometry", "simpleELoss", simpleELoss_ );

    return;
  }
  
  #if USE_TGEANT
  if( CsOpt::Instance()->getOpt( "CsGDMLGeometry", "file", path ) ){
    CsErrLog::Instance()->mes( elInfo, "GDML-imported ROOT geometry will be used" );

    TGeoManager* myGeo = new TGeoManager;
    myGeo->Import(path.c_str());

    
    
    

    if (!gGeoManager)
    {
      CsErrLog::Instance()->mes( elFatal, "ROOT-GDML-imported geometry couldn't be loaded" );
      abort();
    }

    
    //copy the whole GDML-File into the vector of string to make it exportable to PHAST
    std::ifstream file(path.c_str(), std::ios::binary);
    file.seekg(0, std::ios::end);
    std::streamsize size = file.tellg();
    file.seekg(0, std::ios::beg);
    GDMLFile.resize(size);
    file.read(GDMLFile.data(), size);
    
    
    
    usingROOTGeometry_ = false;
    usingROOTGeometryTGEANT_ = true;

    massDefault_ = 0.139579018;  // Default to Pion mass
    CsOpt::Instance()->getOpt( "CsROOTGeometry", "massDefault", massDefault_ );

    simpleELoss_ = 1;
    CsOpt::Instance()->getOpt( "CsROOTGeometry", "simpleELoss", simpleELoss_ );

    return;
  }
  #endif

  usingROOTGeometry_ = false;

  list<string> ZONES;
  CsOpt* opt = CsOpt::Instance();

  ifstream in;
  string msg;

  /********************* RL *******************/
  int i,itmp;
  string strZone;
  char cN[6];
  vector<MatMap*>::iterator iz;  
  for( i=0; i<200; i++ ) {

    sprintf(cN,"%i",i);
    strZone = "Zone_";
    strZone += cN;

    if( opt->getOpt( "CsMaterialMap", strZone, path) ){
      in.open( path.c_str(), ios::in);
      if( !in ){
	msg =  "Map file for " + strZone + ": " + path + " was not found. If you don't like to use Material Map you should comment all strings in option file with tag 'CsMaterialMap  Zone_*'. Be sure that there are not other packages (Traffic, Vertex reconstruction) which require Material Map.";
	CsErrLog::Instance()->mes( elFatal, msg );
      }
      else {
	in.close();
	RL_.push_back( new MatMap );
      }
    }

  }

  /********************* EL *******************/
 
  if( opt->getOpt( "CsMaterialMap", "ELossTarget", path) ) {
    in.open( path.c_str(), ios::in);
    if( !in ){
      msg =  "Map file for 'ELossTarget': " + path + " was not found.";
      CsErrLog::Instance()->mes( elFatal, msg );
    }
    else {
      in.close();
      EL_ = new MatMap();
    } 
  } 
  else EL_=0;

  return;
}



CsMaterialMap::~CsMaterialMap(){
}



bool CsMaterialMap::ReadMaps( string VERS )
{
 if (usingROOTGeometry_ 
    #if USE_TGEANT
    || usingROOTGeometryTGEANT_
    #endif 
     )
   return true;

  bool bin = false;
  char ctmp[50];
  int i,j,k,I,cosys;
  float ftmp;
  int itmp;

  ifstream in;
  string path,msg;

  CsOpt* opt = CsOpt::Instance();

  /********************* RL *******************/

  string strZone = "Zone_0";
  char cN[6];
  vector<MatMap*>::iterator iz;
  
  for( I=0, iz=RL_.begin(); iz!=RL_.end(); iz++, I++ ) {

    sprintf(cN,"%i",I);
    strZone = "Zone_";
    strZone += cN;
    
    while( !opt->getOpt( "CsMaterialMap", strZone, path) ){
      I++;
      sprintf(cN,"%i",I);
      strZone = "Zone_";
      strZone += cN;
      //cout<<" !!!!!! "<<strZone<<", "<<I<<endl;
    };
    
    if(      path.rfind(".dat") == path.size()-4 ) bin = false;
    else if( path.rfind(".map") == path.size()-4 ) bin = true;
    else {
      msg = "Unclear format of file: Nuzhen *.dat or *.map file.";
      CsErrLog::Instance()->mes( elFatal, msg );
    }
    
    in.open( path.c_str(), ios::in);
    if( !in ){
      msg =  "Map file for " + strZone + ": " + path + " was not found.";
      CsErrLog::Instance()->mes( elFatal, msg );
    }
    else{

      msg = "+++ Reading Material map for `" + strZone + "` ... ";
      CsErrLog::Instance()->mes( elInfo, msg );
      
      if( atoi( &VERS[3] ) < 6 && VERS[0] != 'y' ) {
	msg =  "No crosscheck of map for version < 6 .";
	CsErrLog::Instance()->mes( elWarning, msg );
      } 
      else if( atoi( &VERS[3] ) > 6 || 
	       ( atoi( &VERS[3] ) == 6 && atoi( &VERS[7] ) > 2 ) ||
	         VERS[0] == 'y' ) {
	string sline;
	char chead[5];
	char cline[20];
	if( bin ) {
	  in.read(chead,5 );
	  in.read(cline,20);
	  in.read(chead,5 );
	  in.read(cline,20);
	} else {
	  in>>chead>>cline;
	  in>>chead>>cline;
	}
	sline = cline;
	
	if (VERS!=sline) // ***** GEOMETRY VERSION: MATMAP'S != DETECTORS.DAT'S
	  //  We want to be able to reconstruct MC data generated w/ a later,
	  // but accurate, COMGeant geometry version w/ older matmaps in order
	  // to simulate the discrepancy between the real geometry and its
	  // description available to coral at mass production time.
	  //  All files, zebra file of data, detectors.dat and matmaps, have
	  // their geometry pedigree encoded in a tag. Full compatibility is
	  // achieved when all 3 tags agree. But the compatibility requirement
	  // is strong only between zebra and detectors.dat, because it controls
	  // the decoding by coral, based on detectors.dat, of the hit info
	  // stored in the zebra file: any discrepancy there is punished by a
	  // fatal error, cf. "CsGeant3::readGeantHead". On the contrary,
	  // discrepancies between matmaps and the rest only rarely pose serious
	  // threats. There is, yet, the case where the, typically newer,
	  // geometry encoded in detectors.dat, moves some detector outside the
	  // boundary of the, older, matmap, which leads coral to double count
	  // its material.
	  // => Therefore let's emit a warning, but not a fatal error.
	  CsErrLog::msg(elBasicInfo,__FILE__,__LINE__,
	    "%s: The geometry versions of material map(\"%s\") and "
	    "detectors.dat(\"%s\") differ",
			strZone.c_str(),sline.c_str(),VERS.c_str());
      }

      if( bin ) in.read( (char *) & ( (*iz)->CoordSyst ), 4 );
      else      in >> (*iz)->CoordSyst;

      if( abs( (*iz)->CoordSyst ) != 1 ){
	msg = "Radiative length map is needed for " + strZone;
	CsErrLog::Instance()->mes(elFatal, msg );
      }

      if( bin ) {
	in.read( (char*)( &((*iz)->nz)   ), 4 );
	in.read( (char*)( &ftmp ), 4 ); (*iz)->minZ = ftmp;
	in.read( (char*)( &ftmp ), 4 ); (*iz)->maxZ = ftmp;
	in.read( (char*)( &((*iz)->nx)   ), 4 );
	in.read( (char*)( &ftmp ), 4 ); (*iz)->minX = ftmp;
	in.read( (char*)( &ftmp ), 4 ); (*iz)->maxX = ftmp;
	in.read( (char*)( &((*iz)->ny) ),   4 );
	in.read( (char*)( &ftmp ), 4 ); (*iz)->minY = ftmp;
	in.read( (char*)( &ftmp ), 4 ); (*iz)->maxY = ftmp;
	
      } else {
	in>>ctmp>>ctmp>>ctmp>>ctmp>>ctmp>>ctmp>>ctmp>>ctmp>>ctmp; 
	in>>(*iz)->nz>>(*iz)->minZ>>(*iz)->maxZ
	  >>(*iz)->nx>>(*iz)->minX>>(*iz)->maxX
	  >>(*iz)->ny>>(*iz)->minY>>(*iz)->maxY;
      }

      msg = strZone + ": ";
      if ( (*iz)->CoordSyst < 0 ) { // cylindrical coordinates
        if( (*iz)->nx < 0 ) msg += "Non-equidistant step in R; ";
        else                msg += "Equidistant step in R; ";
        if( (*iz)->nz < 0 ) msg += "Non-equidistant step in PHI; "; 
        else                msg += "Equidistant step in PHI; ";
        if( (*iz)->ny < 0 ) msg += "Non-equidistant step in Z."; 
        else                msg += "Equidistant step in Z.";
      } else {                      // rectangular coodinates
        if( (*iz)->nx < 0 ) msg += "Non-equidistant step in X; "; 
        else                msg += "Equidistant step in X; ";
        if( (*iz)->ny < 0 ) msg += "Non-equidistant step in Y; ";
        else                msg += "Equidistant step in Y; ";
        if( (*iz)->nz < 0 ) msg += "Non-equidistant step in Z."; 
        else                msg += "Equidistant step in Z.";
      }
      CsErrLog::Instance()->mes( elInfo, msg );

      if( !bin ) 
	in>>ctmp>>ctmp>>ctmp>>ctmp>>ctmp>>ctmp>>ctmp; // 7 words: IAX IEQD ...
      
      /************  READ grid in Z  ***************/
      if( bin ) {
	in.read( (char *)& itmp , 4 );
	in.read( (char *)& itmp , 4 );
	in.read( (char *)& itmp , 4 );
      } else
	in>>ftmp>>ftmp>>ftmp;    //skip: IAX IEQD Ncells
      
      (*iz)->nz = abs( (*iz)->nz );
      for( i = 0; i < (*iz)->nz; i++ ){  // Low edges (grid) in Z
	if( bin ) in.read( (char *)& ftmp , 4 );
	else      in>>ftmp;
	(*iz)->z.push_back( ftmp );
      }
      
      /************  READ grid in X  ***************/
      if( bin ) {
	in.read( (char *)& ftmp , 4 );
	in.read( (char *)& ftmp , 4 );
	in.read( (char *)& ftmp , 4 );
      } else
	in>>ftmp>>ftmp>>ftmp;    // IAX IEQD Ncells:
      
      (*iz)->nx = abs( (*iz)->nx );
      for( i = 0; i < (*iz)->nx; i++ ){  // Low edges (grid) in X
	if( bin ) in.read( (char *)& ftmp , 4 );
	else      in>>ftmp;
	(*iz)->x.push_back( ftmp );
      }
      
      /************  READ grid in Y  ***************/
      if( bin ) {
	in.read( (char *)& ftmp , 4 );
	in.read( (char *)& ftmp , 4 );
	in.read( (char *)& ftmp , 4 );
      } else
	in>>ftmp>>ftmp>>ftmp;    // IAX IEQD Ncells:
      
      (*iz)->ny = abs( (*iz)->ny );
      for( i = 0; i < (*iz)->ny; i++ ){  // Low edges (grid) in Y
	if( bin ) in.read( (char *)& ftmp , 4 );
	else      in>>ftmp;
	(*iz)->y.push_back( ftmp );
      } 
      
      /*************  READ the table  **************/
      
      if( !bin ) in>>ctmp>>ctmp>>ctmp>>ctmp;   // read 4 words

      int count = 0;
      for( k = 0; k < (*iz)->nz; k++ ){             // PHI points
	for( j = 0; j < (*iz)->nx; j++ ){           // R   points
	  for( i = 0; i < (*iz)->ny; i++ ){         // Z   points
	    if( in.eof() ) CsErrLog::Instance()->mes(elFatal,"File for Zone_? is too short.");
	    if( bin ) 
	      in.read( (char *)& ftmp , 4 );
	    else
	      in>>ftmp>>ftmp>>ftmp>>ftmp;   // first 3 not needed
	    (*iz)->Val.push_back( ftmp );
	  }
	}
	
      }
      in.close();

      if( (*iz)->CoordSyst < 0 ) {
	msg = "Coordinate system is CYLINDRICAL for " + strZone;
	CsErrLog::Instance()->mes(elInfo, msg );
	(*iz)->minZ *= DEG2RAD;
	(*iz)->maxZ *= DEG2RAD;
	for( i = 0; i < (*iz)->nz; i++ ){  // Low edges (grid) in PHI
	  (*iz)->z[i] *= DEG2RAD;
	  (*iz)->Cos.push_back( cos( (*iz)->z[i] ) );
	  (*iz)->Sin.push_back( sin( (*iz)->z[i] ) );
	}
      }
      else if ( (*iz)->CoordSyst > 0 ) {
	msg = "Coordinate system is RECTANGULAR for " + strZone;
	CsErrLog::Instance()->mes( elInfo, msg );
      } 
      else {
	msg = "Coordinate system is indefinite for " + strZone;
	CsErrLog::Instance()->mes( elFatal, msg );
      }

    }

  }
  
  /********************* EL *******************/

  if( opt->getOpt( "CsMaterialMap", "ELossTarget", path) ) {
    
    if(      path.rfind(".dat") == path.size()-4 ) bin = false;
    else if( path.rfind(".map") == path.size()-4 ) bin = true;
    else {
      msg = "Unclear format of file: Nuzhen *.dat or *.map file.";
      CsErrLog::Instance()->mes( elFatal, msg );
    }
    
    in.open( path.c_str(), ios::in);
    if( !in ){
      msg =  "Map file for ELossTarget: " + path + " was not found.";
      CsErrLog::Instance()->mes( elFatal, msg );
    }
    else{

      msg = "+++ Reading map for `ELossTarget` ... ";
      CsErrLog::Instance()->mes( elInfo, msg );

      if( atoi( &VERS[3] ) < 6 && VERS[0] != 'y' ) {
	msg =  "No crosscheck of map for version < 6 .";
	CsErrLog::Instance()->mes( elWarning, msg );
      } 
      else if( atoi( &VERS[3] ) > 6 || 
	       ( atoi( &VERS[3] ) == 6 && atoi( &VERS[7] ) > 2 ) ||
	        VERS[0] == 'y' ) {
	string sline;
	char chead[5];
	char cline[20];
	if( bin ) {
	  in.read(chead,5 );
	  in.read(cline,20);
	  in.read(chead,5 );
	  in.read(cline,20);
	} else {
	  in>>chead>>cline;
	  in>>chead>>cline;
	}
	sline = cline;

	if (VERS!=sline) // ***** GEOMETRY VERSION: MATMAP'S != DETECTORS.DAT'S
	  // Cf. supra for explanations.
	  CsErrLog::msg(elBasicInfo,__FILE__,__LINE__,
	    "%s: The geometry versions of material map(\"%s\") and "
	    "detectors.dat(\"%s\") differ",
			strZone.c_str(),sline.c_str(),VERS.c_str());
      }
     
      if( bin ) in.read( (char *)& ( EL_->CoordSyst ), 4 );
      else      in >> EL_->CoordSyst;

      if( abs( EL_->CoordSyst ) != 3 ){
	msg = "Map of dE/dX is needed for 'ELossTarget'";
	CsErrLog::Instance()->mes(elFatal, msg );
      }

      if( bin ) {
	in.read( (char*)& ( EL_->nz   ), 4 );
	in.read( (char*)& ( ftmp ), 4 ); EL_->minZ = ftmp;
	in.read( (char*)& ( ftmp ), 4 ); EL_->maxZ = ftmp;
	in.read( (char*)& ( EL_->nx   ), 4 );
	in.read( (char*)& ( ftmp ), 4 ); EL_->minX = ftmp;
	in.read( (char*)& ( ftmp ), 4 ); EL_->maxX = ftmp;
	in.read( (char*)& ( EL_->ny   ), 4 );
	in.read( (char*)& ( ftmp ), 4 ); EL_->minY = ftmp;
	in.read( (char*)& ( ftmp ), 4 ); EL_->maxY = ftmp;
      } else {
	in>>ctmp>>ctmp>>ctmp>>ctmp>>ctmp>>ctmp>>ctmp>>ctmp>>ctmp; 
	in>>EL_->nz>>EL_->minZ>>EL_->maxZ
	  >>EL_->nx>>EL_->minX>>EL_->maxX
	  >>EL_->ny>>EL_->minY>>EL_->maxY;
      }
      
      msg = "ELossTarget: ";
      if ( EL_->CoordSyst < 0) { // cyclindrical coordinates
        if( EL_->nx < 0 ) msg += "Non-equidistant step in R; "; 
        else              msg += "Equidistant step in R; ";
        if( EL_->nz < 0 ) msg += "Non-equidistant step in PHI; ";
        else              msg += "Equidistant step in PHI; ";
        if( EL_->ny < 0 ) msg += "Non-equidistant step in Z."; 
        else              msg += "Equidistant step in Z.";
      } else {                   // rectangular coordinates
        if( EL_->nx < 0 ) msg += "Non-equidistant step in X; "; 
        else              msg += "Equidistant step in X; ";
        if( EL_->ny < 0 ) msg += "Non-equidistant step in Y; ";
        else              msg += "Equidistant step in Y; ";
        if( EL_->nz < 0 ) msg += "Non-equidistant step in Z."; 
        else              msg += "Equidistant step in Z.";
      }
      CsErrLog::Instance()->mes( elInfo, msg );

      if( !bin ) 
	in>>ctmp>>ctmp>>ctmp>>ctmp>>ctmp>>ctmp>>ctmp;    // 7 words: IAX IEQD ...
     
      /************  READ grid in Z  ***************/
      if( bin ) {
	in.read( (char *)& itmp , 4 );
	in.read( (char *)& itmp , 4 );
	in.read( (char *)& itmp , 4 );
      } else
	in>>ftmp>>ftmp>>ftmp;    //skip: IAX IEQD Ncells
      
      EL_->nz = abs( EL_->nz );
      for( i = 0; i < EL_->nz; i++ ){  // Low edges (grid) in Z
	if( bin ) in.read( (char *)& ftmp , 4 );
	else      in>>ftmp;
	EL_->z.push_back( ftmp );
      }
      
      /************  READ grid in X  ***************/
      if( bin ) {
	in.read( (char *)& itmp , 4 );
	in.read( (char *)& itmp , 4 );
	in.read( (char *)& itmp , 4 );
      } else
	in>>ftmp>>ftmp>>ftmp;    // IAX IEQD Ncells:
      
      EL_->nx = abs( EL_->nx );
      for( i = 0; i < EL_->nx; i++ ){  // Low edges (grid) in X
	if( bin ) in.read( (char *)& ftmp , 4 );
	else      in>>ftmp;
	EL_->x.push_back( ftmp );
      }
      
      /************  READ grid in Y  ***************/
      if( bin ) {
	in.read( (char *)& itmp , 4 );
	in.read( (char *)& itmp , 4 );
	in.read( (char *)& itmp , 4 );
      } else
	in>>ftmp>>ftmp>>ftmp;    // IAX IEQD Ncells:
      
      EL_->ny = abs( EL_->ny );
      for( i = 0; i < EL_->ny; i++ ){  // Low edges (grid) in Y
	if( bin ) in.read( (char *)& ftmp , 4 );
	else      in>>ftmp;
	EL_->y.push_back( ftmp );
      } 

      /*************  READ the table  **************/
      
      if( !bin ) in>>ctmp>>ctmp>>ctmp>>ctmp;   // read 4 words
      
      int count = 0;
      for( k = 0; k < EL_->nz; k++ ){             // PHI points
	for( j = 0; j < EL_->nx; j++ ){           // R   points
	  for( i = 0; i < EL_->ny; i++ ){         // Z   points
	    if( in.eof() ) CsErrLog::Instance()->mes(elFatal,"File for 'ELossTarget' is too short.");
	    if( bin ) 
	      in.read( (char *)& ftmp , 4 );
	    else 
	      in>>ftmp>>ftmp>>ftmp>>ftmp;   // first 3 not needed
	    EL_->Val.push_back( ftmp );
	  }
	}
	
      }
      in.close();

      if( EL_->CoordSyst < 0 ) {
	msg = "Coordinate system is CYLINDRICAL for 'ELossTarget'";
	CsErrLog::Instance()->mes(elInfo, msg );
	EL_->minZ *= DEG2RAD;
	EL_->maxZ *= DEG2RAD;
	for( i = 0; i < EL_->nz; i++ ){  // Low edges (grid) in PHI
	  EL_->z[i] *= DEG2RAD;
	  EL_->Cos.push_back( cos( EL_->z[i] ) );
	  EL_->Sin.push_back( sin( EL_->z[i] ) );
	}
      }
      else if ( EL_->CoordSyst > 0 ) {
	msg = "Coordinate system is RECTANGULAR for 'ELossTarget'";
	CsErrLog::Instance()->mes( elInfo, msg );
      } 
      else {
	msg = "Coordinate system is indefinite for 'ELossTarget'";
	CsErrLog::Instance()->mes( elFatal, msg );
      }

    }

  }


  return true;
} 



void CsMaterialMap::getRadLength( const CsHelix &hel, bool direc, float &RadLen, float &StepZ)
{
  THlx Thel;
  Thel.ImportHelix(hel);

  getRadLength( Thel, direc, RadLen, StepZ);

  RadLen *= 10;
  StepZ  *= 10;

  return;
}



void CsMaterialMap::getRadLength( const THlx &hel, bool direc, float &RadLen, float &StepZ)
{
  if (usingROOTGeometry_
#if USE_TGEANT
    && !usingROOTGeometryTGEANT_
#endif    
  )
    {
      int sign = direc ? 1 : -1;
      // The geometry uses the COMGEANT coordinate system.
      gGeoManager->InitTrack(hel(0), hel(1), hel(2),
                             sign*hel.DirCos(1), sign*hel.DirCos(2), sign*hel.DirCos(3));
      // Find current material.
      TGeoMaterial *cmat = gGeoManager->GetCurrentVolume()->GetMedium()->GetMaterial();

      // Normal case.
      RadLen = cmat->GetRadLen();

      // Special case for the vacuum which is the predefined GEANT3 material with
      // IMATE == 16.  g2root makes this the material's UniqueID in the ROOT geometry:
      if (cmat->GetUniqueID() == 16)
	RadLen = 1.e16; // Value used in GEANT3.

      // Find the next boundary within 5 cms along the current
      // direction.  In gases, we allow longer steps.
      if (RadLen < 3244) // i.e. shorter RadLen than the RICH gas.
	gGeoManager->FindNextBoundaryAndStep(5);
      else
	gGeoManager->FindNextBoundaryAndStep();

      StepZ = sign * gGeoManager->GetStep() / sqrt(1 + hel(3)*hel(3) + hel(4)*hel(4));
      // Make sure the boundary is crossed if tracking actually
      // extrapolates across the whole calculated step length.
      StepZ += direc ? 0.01 : -0.01;

      return;
    }
    
    #if USE_TGEANT
    //for tgeant root geometry from gdml we use coral coordinate system
    if (usingROOTGeometryTGEANT_)
    {
      int sign = direc ? 1 : -1;
      // The geometry uses the CORAL coordinate system.
      // therefore init the track with transformed coordinates
      gGeoManager->InitTrack(hel(1), hel(2), hel(0),
                             sign*hel.DirCos(2), sign*hel.DirCos(3), sign*hel.DirCos(1));
      // Find current material.
      TGeoMaterial *cmat = gGeoManager->GetCurrentVolume()->GetMedium()->GetMaterial();

      // Normal case.
      RadLen = cmat->GetRadLen();


      // Find the next boundary within 5 cms along the current
      // direction.  In gases, we allow longer steps.
      if (RadLen < 3244) // i.e. shorter RadLen than the RICH gas.
  gGeoManager->FindNextBoundaryAndStep(5);
      else
  gGeoManager->FindNextBoundaryAndStep();

      StepZ = sign * gGeoManager->GetStep() / sqrt(1 + hel(3)*hel(3) + hel(4)*hel(4));
      // Make sure the boundary is crossed if tracking actually
      // extrapolates across the whole calculated step length.
      StepZ += direc ? 0.01 : -0.01;
     return;      
    }
    
#endif

  const float correction = 0.01;
  float pos_z = hel(0);

  vector<MatMap*>::iterator iz;
  for( iz=RL_.begin(); iz!=RL_.end(); iz++ ) {
    if( (*iz)->insideZ( pos_z ) ) { // input point inside of map
      getRadLengthZone( (*iz), hel, direc, RadLen, StepZ);
      StepZ += ( direc ? correction : -correction );
      break;
    }
  }

  if( iz==RL_.end() ) {             // input point outside of maps
    StepZ = ( direc ? 9999999 : -9999999 );
    float Step;
    for( iz=RL_.begin(); iz!=RL_.end(); iz++ ) {
      Step = (*iz)->StepToBorder( pos_z, direc );
      if( fabs( Step ) < fabs( StepZ ) ) StepZ = Step;
    }
    RadLen = 30423;
  }
  StepZ += ( direc ? correction : -correction );

  //cout << "Coral " << RadLen << " " << StepZ << endl;

  return;
}




void CsMaterialMap::getRadLengthZone( MatMap *matmap, THlx hel, bool direc, float &RadLen, float &StepZ)
{
  assert(!usingROOTGeometry_);

  double pos_x,pos_y,pos_z;
  int i,j,k;
  
  pos_x=hel(1);
  pos_y=hel(2);
  pos_z=hel(0);

  if( matmap->CoordSyst < 0 ) {
    double PHI,RAD,SIN,COS;
    if( pos_x==0 && pos_y==0 ){
      RAD=0.000000001;
      COS=1;
      SIN=0;
    }else{
      RAD = sqrt(pos_x*pos_x + pos_y*pos_y);
      if(RAD > matmap->maxX) {RadLen=30423; StepZ=( direc ? 100 : -100 ); return; }
      
      COS = pos_x/RAD;
      SIN = pos_y/RAD;
    }
    
    if( matmap->nz > 1 ) {
     i=1;
     if( SIN < 0){
       while ( i<matmap->nz && matmap->Sin[i] < 0 && COS > matmap->Cos[i] ) {i++;}
     }else {
       while ( i<matmap->nz && matmap->Sin[i] <= 0 ) {i++;}
       while ( i<matmap->nz && COS <= matmap->Cos[i] ) {i++;}    
     }
     i=i-1;
    }
    else i=0;
    
    j=1;
    while ( j < matmap->nx && RAD > matmap->x[j] ) {j++;} 
    j=j-1;
    
    k=1;
    while ( k < matmap->ny && pos_z > matmap->y[k] ) {k++;}
    k=k-1;
    
    RadLen = matmap->Val[ i * (matmap->nx)*(matmap->ny) + j * (matmap->ny) + k ];
    if( hel(5)==0 ) hel(5)=1.;
    hel(5) = (direc ? fabs(hel(5)) : -fabs(hel(5)) );
    StepZ  = getStepSize(hel,matmap,i,j,k,RadLen);
  }
  else {
    if( pos_x <= matmap->minX || pos_x >= matmap->maxX ) {
      RadLen=30423; StepZ=( direc ? 100 : -100 ); return;}
    if( pos_y <= matmap->minY || pos_y >= matmap->maxY ) {
      RadLen=30423; StepZ=( direc ? 100 : -100 ); return;}

    // par 5 shows direction along Z
    hel(5) = (direc ? 1 : -1 );
    
    i=1;
    while ( i < matmap->nx && pos_x > matmap->x[i] ) {i++;} 
    if( i < matmap->nx && pos_x == matmap->x[i] ) {
      if(       direc && hel.DirCos(2) < 0 ) i=i-1;
      else if( !direc && hel.DirCos(2) > 0 ) i=i-1;
    } else i=i-1;
    
    j=1;
    while ( j < matmap->ny && pos_y > matmap->y[j] ) {j++;} 
    if( j < matmap->ny && pos_y == matmap->y[j] ) {
      if(       direc && hel.DirCos(3) < 0 ) j=j-1;
      else if( !direc && hel.DirCos(3) > 0 ) j=j-1;
    } else j=j-1;
    
    k=1;
    while ( k < matmap->nz && pos_z > matmap->z[k] ) {k++;} 
    k=k-1;
    
    RadLen = matmap->Val[ k * (matmap->nx)*(matmap->ny) + i * (matmap->ny) + j ];
    
    StepZ   = getStepSize(hel,matmap,i,j,k,RadLen);
    
  }

  return;
}



///////////////////////////////////////////////////////////////////////////



float CsMaterialMap::getStepSize(const THlx &hel, MatMap *matmap,
				 int i, int j, int k , float RadLen)
{
  assert(!usingROOTGeometry_);

  int iZone;
  float Step,StepZ;

  double pos_x,pos_y,pos_z;
  double vct_x,vct_y,vct_z;

  double intr_x,  intr_y,  intr_z;
  double intrYZ_x,intrYZ_y,intrYZ_z;
  double intrXZ_x,intrXZ_y,intrXZ_z;
  double intrXY_x,intrXY_y,intrXY_z;

  float t1_y;
  float b1_x,b1_y,b1_z;
  float b2_z;
  float b4_x;
  
  pos_x = hel(1);
  pos_y = hel(2);
  pos_z = hel(0);
  
  if( hel(5)<0 ) {
    vct_z = -hel.DirCos(1);
    vct_x = -hel.DirCos(2);
    vct_y = -hel.DirCos(3);
  }else {
    vct_z = hel.DirCos(1);
    vct_x = hel.DirCos(2);
    vct_y = hel.DirCos(3);
  }

  if( matmap->CoordSyst < 0 ) {
    
    if( k < matmap->ny - 1 ) {
      if( vct_z > 0 ) StepZ = matmap->y[k+1] - pos_z;
      else            StepZ = matmap->y[k]   - pos_z;}
    else if( k == matmap->ny - 1 ) {
      if( vct_z > 0 ) StepZ = matmap->maxY  - pos_z;
      else            StepZ = matmap->y[k]  - pos_z;}
    else {
      StepZ = 0.;  // just to avoid compiler warning
      CsErrLog::Instance()->mes( elFatal, "Can't determine step length" );
    }

    if( StepZ > 2 ) StepZ = 2;
    if( StepZ < -2 ) StepZ = -2;
    
    //cout<<"CMM: StepZ="<<StepZ<<""<<endl;
    return StepZ;

  }
  else {

    if( i < matmap->nx - 1 ){
      b1_x = matmap->x[i];
      b4_x = matmap->x[i+1];}
    else if( i == matmap->nx - 1 ){
      b1_x = matmap->x[i];
      b4_x = matmap->maxX;}
    
    if( j < matmap->ny - 1 ){
      b1_y = matmap->y[j];
      t1_y = matmap->y[j+1];}
    else if( j == matmap->ny - 1 ){
      b1_y = matmap->y[j];
      t1_y = matmap->maxY;}

    if( k < matmap->nz - 1 ){
      b1_z = matmap->z[k];
      b2_z = matmap->z[k+1];}
    else if( k == matmap->nz - 1 ){
      b1_z = matmap->z[k];
      b2_z = matmap->maxZ;}

  }

  //********* YZ plane **********
  
  if(vct_z != 0){
    if( vct_y != 0 ){
      if( vct_z > 0 ) intrYZ_y = pos_y + (b2_z-pos_z)*vct_y/vct_z;
      else            intrYZ_y = pos_y + (b1_z-pos_z)*vct_y/vct_z;
      
      if( intrYZ_y > t1_y ){
        intrYZ_y=t1_y;
        intrYZ_z=pos_z+(t1_y-pos_y)*vct_z/vct_y;
      }
      else if( intrYZ_y < b1_y ){
        intrYZ_y=b1_y;
        intrYZ_z=pos_z+(b1_y-pos_y)*vct_z/vct_y;
      }
      else{
        intrYZ_z = (vct_z > 0 ? b2_z : b1_z);
      }
    }
    else{
      intrYZ_y=pos_y;
      intrYZ_z = (vct_z > 0 ? b2_z : b1_z);
    }
  }
  else{
    if( vct_y != 0 ){
      intrYZ_z=pos_z;
      intrYZ_y = (vct_y >= 0 ? t1_y : b1_y);
    }
    else{
      intrYZ_z=pos_z;
      intrYZ_y=pos_y;
      //cout<<"NO MOTION in YZ plane"<<endl;
    }
  }
   
  //********* XZ plane **********  
  if(vct_z != 0){
    if( vct_x != 0 ){
      if( vct_z > 0 ) intrXZ_x = pos_x + (b2_z-pos_z)*vct_x/vct_z;
      else            intrXZ_x = pos_x + (b1_z-pos_z)*vct_x/vct_z;
      
      if( intrXZ_x > b4_x ){
        intrXZ_x=b4_x;
        intrXZ_z=pos_z+(b4_x-pos_x)*vct_z/vct_x;
      }
      else if( intrXZ_x < b1_x ){
        intrXZ_x=b1_x;
        intrXZ_z=pos_z+(b1_x-pos_x)*vct_z/vct_x;
      }
      else{
        intrXZ_z = (vct_z > 0 ? b2_z : b1_z);
      }
    }
    else{
      intrXZ_x=pos_x;
      intrXZ_z = (vct_z > 0 ? b2_z : b1_z);
    }
  }
  else{
    if( vct_x != 0 ){
      intrXZ_z=pos_z;
      intrXZ_x = (vct_x >= 0 ? b4_x : b1_x);
    }
    else{
      intrXZ_z=pos_z;
      intrXZ_x=pos_x;
      //cout<<"NO MOTION in XZ plane"<<endl;
    }
  }
   
  //********* XY plane **********
  
  if(vct_x != 0){
    if( vct_y != 0 ){
      if( vct_x > 0 ) intrXY_y = pos_y + (b4_x-pos_x)*vct_y/vct_x;
      else            intrXY_y = pos_y + (b1_x-pos_x)*vct_y/vct_x;
      
      if( intrXY_y > t1_y ){
        intrXY_y=t1_y;
        intrXY_x=pos_x+(t1_y-pos_y)*vct_x/vct_y;
      }
      else if( intrXY_y < b1_y ){
        intrXY_y=b1_y;
        intrXY_x=pos_x+(b1_y-pos_y)*vct_x/vct_y;
      }
      else{
        intrXY_x = (vct_x > 0 ? b4_x : b1_x);
      }
    }
    else{
      intrXY_y=pos_y;
      intrXY_x = (vct_x > 0 ? b4_x : b1_x);
    }
  }
  else{
    if( vct_y != 0 ){
      intrXY_x=pos_x;
      intrXY_y = (vct_y >= 0 ? t1_y : b1_y);
    }
    else{
      intrXY_x=pos_x;
      intrXY_y=pos_y;
      //cout<<"NO MOTION in XY plane"<<endl;
    }
  }
   
  //******* CHECK FOR INTERSECTIONS: YZ, XZ, XY *******  

  if( vct_z > 0 )
    intr_z = (intrYZ_z < intrXZ_z ? intrYZ_z : intrXZ_z);
  else if( vct_z < 0 )
    intr_z = (intrYZ_z < intrXZ_z ? intrXZ_z : intrYZ_z);
  else 
    intr_z = pos_z;
  
  if( vct_x > 0 )
    intr_x = (intrXZ_x < intrXY_x ? intrXZ_x : intrXY_x);
  else if( vct_x < 0 )
    intr_x = (intrXZ_x < intrXY_x ? intrXY_x : intrXZ_x);
  else 
    intr_x = pos_x;
  
  if( vct_y > 0 )
    intr_y = (intrYZ_y < intrXY_y ? intrYZ_y : intrXY_y);
  else if( vct_y < 0 )
    intr_y = (intrYZ_y < intrXY_y ? intrXY_y : intrYZ_y);
  else
    intr_y = pos_y;
  
  
  if(intr_x==b4_x)      i++;
  else if(intr_x==b1_x) i--;
  
  if(intr_y==t1_y)      j++;
  else if(intr_y==b1_y) j--;

  if(intr_z==b2_z)      k++;
  else if(intr_z==b1_z) k--;

  Step  = sqrt((intr_x-pos_x)*(intr_x-pos_x)+
	       (intr_y-pos_y)*(intr_y-pos_y)+
	       (intr_z-pos_z)*(intr_z-pos_z) );
  StepZ = (intr_z-pos_z);

  //cout<<"   ME: pos_x = "<<pos_x<<", pos_y = "<<pos_y<<", pos_z = "<<pos_z<<endl;
  /*
  if( Step/RadLen < 0.1 ){
    if( i<matmap->nx && i>=0  &&  j<matmap->ny && j>=0  &&  k<matmap->nz && k>=0 ){
      if( fabs( matmap->Val[ k * (matmap->nx)*(matmap->ny) + i * (matmap->ny) + j ] - RadLen ) < 0.05 * RadLen ){
	hel(1)=intr_x;
	hel(2)=intr_y;
	hel(0)=intr_z;
	StepZ += getStepSize(hel,matmap,i,j,k,RadLen);
      }
    }
  }
  else
    StepZ *= 0.1 * RadLen / Step;
  */
  
  return StepZ;
}



///////////////////////////////////////////////////////////////////////////



float CsMaterialMap::getRadLength( float pos_x, float pos_y, float pos_z)
{
  if (usingROOTGeometry_
#if USE_TGEANT
    && !usingROOTGeometryTGEANT_
#endif
    
  )
    {
      // The geometry uses the COMGEANT coordinate system.
      gGeoManager->IsSameLocation(pos_z, pos_x, pos_y, kTRUE);
      // Find current material.
      TGeoMaterial *cmat = gGeoManager->GetCurrentVolume()->GetMedium()->GetMaterial();

      // Normal material.
      double radLen = cmat->GetRadLen()*10;

      // Special case for the vacuum which is the predefined GEANT3 material with
      // IMATE == 16.  g2root makes this the material's UniqueID in the ROOT geometry:
      if (cmat->GetUniqueID() == 16)
	radLen = 1.e16; // Value used in GEANT3.

      return radLen;
    }
    
#if USE_TGEANT    
      if (usingROOTGeometryTGEANT_)
    {
      // The geometry uses the COMGEANT coordinate system.
      gGeoManager->IsSameLocation(pos_x, pos_y, pos_z, kTRUE);
      // Find current material.
      TGeoMaterial *cmat = gGeoManager->GetCurrentVolume()->GetMedium()->GetMaterial();

      // Normal material.
      double radLen = cmat->GetRadLen()*10;
      return radLen;
    }
#endif    

  float RadLen = 0;
  int i,j,k;

  pos_x /= 10;
  pos_y /= 10;
  pos_z /= 10;
  
  vector<MatMap*>::iterator iz;
  for( iz=RL_.begin(); iz!=RL_.end(); iz++ ) {
    
    // go to the next map, if the z position is not in the current
    if( !((*iz)->insideZ(pos_z)) ) continue;
    
    if( (*iz)->CoordSyst < 0 ) {
      double PHI,RAD,SIN,COS;
      if( pos_x==0 && pos_y==0 ){
	RAD=0.000000001;
	COS=1;
	SIN=0;
      }else{
	RAD = sqrt(pos_x*pos_x + pos_y*pos_y);
	if(RAD > (*iz)->maxX) {RadLen=30423;  break;; }
	
	COS = pos_x/RAD;
	SIN = pos_y/RAD;
      }
      
      if( (*iz)->nz > 1 ) {
       i=1;
       if( SIN < 0){
 	 while ( i < (*iz)->nz && (*iz)->Sin[i] < 0 && COS > (*iz)->Cos[i] ) {i++;}
       }else {
 	 while ( i < (*iz)->nz && (*iz)->Sin[i] <= 0 ) {i++;}
	 while ( i < (*iz)->nz && COS <= (*iz)->Cos[i]  ) {i++;}    
       }
       i=i-1;
      }
      else i=0;
      
      j=1;
      while ( j < (*iz)->nx && RAD > (*iz)->x[j] ) {j++;} 
      j=j-1;
      
      k=1;
      while ( k < (*iz)->ny && pos_z > (*iz)->y[k] ) {k++;}
      k=k-1;
      
      RadLen = (*iz)->Val[ i * ((*iz)->nx)*((*iz)->ny) + j * ((*iz)->ny) + k ];
    
    }
    else {
      if( pos_x <= (*iz)->minX || pos_x >= (*iz)->maxX ) {
	RadLen=30423; break;}
      if( pos_y <= (*iz)->minY || pos_y >= (*iz)->maxY ) {
	RadLen=30423; break;}
      
      i=1;
      while ( i < (*iz)->nx && pos_x > (*iz)->x[i]  ) {i++;} 
      i=i-1;
      
      j=1;
      while ( j < (*iz)->ny && pos_y > (*iz)->y[j] ) {j++;} 
      j=j-1;
      
      k=1;
      while ( k < (*iz)->nz && pos_z > (*iz)->z[k] ) {k++;} 
      k=k-1;
      
      RadLen = (*iz)->Val[ k * ((*iz)->nx)*((*iz)->ny) + i * ((*iz)->ny) + j ];
      
    }
  }

  return RadLen;
}

 

float CsMaterialMap::getdE( const THlx& hel, float Len )
{
  if (usingROOTGeometry_ 
#if USE_TGEANT    
    || usingROOTGeometryTGEANT_
#endif   
     )
    {
      
      // The geometry uses the COMGEANT coordinate system.
      if (usingROOTGeometry_)
        gGeoManager->InitTrack(hel(0), hel(1), hel(2),
                             hel.DirCos(1), hel.DirCos(2), hel.DirCos(3));
#if USE_TGEANT
      else if (usingROOTGeometryTGEANT_)
        gGeoManager->InitTrack(hel(1), hel(2), hel(0),
                             hel.DirCos(2), hel.DirCos(3), hel.DirCos(1));
#endif
      // Find current material.
      const TGeoMaterial *cmat = gGeoManager->GetCurrentVolume()->GetMedium()->GetMaterial();

      // Special case for the vacuum which has IMATE == 16.  This ends
      // up being the material's UniqueID after the translation to a
      // ROOT geometry by g2root.
      if (cmat->GetUniqueID() == 16)
	return 0;

      double A = cmat->GetA();
      double Z = cmat->GetZ();
      double rho = cmat->GetDensity();
      double p = fabs(1/hel(5));
      double E = hypot(p, massDefault_);
      double beta = p / E;
      double betagamma = p / massDefault_;

      const double me = .51099891e-3;         // electron mass

      double res;
      if ( simpleELoss_ ) {
	// simple most-probable energy loss based on PDG 2010 (27.10)

	// Use most probable Eloss:
	// Ksi is the collisional thickness.
	double Ksi = 0.153537e-3 * Z/A * rho / beta / beta * Len;
	// ion is the average ionization energy: 
	double ion = 16.e-9 * pow(Z,0.9);
	
	// in the formula below, sigma is neglected (density correction):
	res = log(2*me*betagamma*betagamma*Ksi/ion/ion);
	res = Ksi * (res - beta*beta + 0.198); 

	//cout << "root " << res << endl;

      } else {
	// Catarina Quintans' code:

	double gamma = E / massDefault_;

	const double avog =6.02214199e+23;      // Avogadro constant
	const double re = 2.817940325e-13;
	const double vlight=29979245800.;
	const double h=6.58211915e-25;
	const double alpha=1./137.0359991;

	// new calculations including all terms, and Bethe-Bloch
	// Bethe-Bloch according to D.E.Groom et al (2001) (see PDG)
        double ksix = 0.15353747e-3 * Z/A * rho/beta/beta;
	//        double ion = 16.e-9*pow(Z,0.9);
	// mean excitation energy according to Sternheimer and Peierls PhysRev B3 (1971) 3681
        double ion = (9.76*Z + 58.8*pow(Z,-0.19))*1.e-6;
        double Tmax = 2.*me*betagamma*betagamma/(1.+ 2. * gamma * me/massDefault_ + me*me / massDefault_ / massDefault_);
        res = log(2.*me*Tmax*betagamma*betagamma/ion/ion) - 2.*beta*beta;

        double X = log10(betagamma);
        double C = 2. * log(h*sqrt(avog*rho*Z/A*4.*M_PI*re*vlight*vlight)/ion) - 1.;
        double X0, X1, m;
        double delta;

	X0=0.;
        if (ion < 1.e-7) {
                X1 = 2.;
                m = 3.;
                if (C > -3.681) X0 = 0.2;
                if (C <= -3.681) X0 = -0.326 * C - 1.;
        }
        if (ion >= 1.e-7) {
                X1 = 3.0;
                m = 3.0;
                if (C > -5.215) X0 = 0.2;
                if (C <= -5.215) X0 = -0.326 * C - 1.5;
        }

	// density term delta:
        if (X < X0) delta = 0.;
        if (X >= X0 && X < X1) {
                double Xa = -C / 4.606;
                double a = 4.606 * (Xa - X0)/pow(X1 - X0,m);
                delta = 4.606 * X + C + a * pow(X1 - X,m);
        }
        if (X  > X1) delta = 4.606 * X + C;

	// average energy loss from ionization:
        double extra = avog * re*re * me * Z/A * alpha *(log(2.*E/massDefault_) - 1./3.*log(2.*Tmax/me)) * log(2.*Tmax/me) * log(2.*Tmax/me);
	res = ksix * (res - delta + Tmax*Tmax/4./gamma/gamma / massDefault_ / massDefault_) + extra*rho;

	// contribution from Bremmstrahlung:
	// from Richard-Serre, CERN 71-18:
        double brem = 4.* alpha * avog * re*re * Z*Z/A * me*me / massDefault_ / massDefault_ * E;
        brem *= (log(12.*gamma*pow(Z,-1./3.)/5.) - 1./3.) * rho;

	res += brem;

	// contribution from electron-positron pairs creation:
	// from Richard-Serre, CERN 71-18:
        double Ez = (28. - 0.36 *Z + 0.002 *Z*Z);

        double f=0.;
        if (E < Ez) f = 1.;
        else if (E >=Ez) f = (16./9. * log(183. * pow(Z,-1./3.)) + 1.)/(16./9. * log(gamma) - 14./9. + log(2.));
        double pair = avog/A * me/massDefault_ * (alpha*Z*re)*(alpha*Z*re)/M_PI * E * (19.3 * log(gamma) - 53.7) *f;
        pair = pair * rho;

	res += pair;

	// contribution from nuclear interactions:
	// my parametrization from PDG tabulated values:
	//	double nuc = 4.3195e-7 * E - 8.312e-10 * E*E + 2.79164e-12 * E*E*E;
	// From Richard-Serre CERN 71-18: http://cdsweb.cern.ch/record/190197/files/p1.pdf
	double nuc = 3.357e-7*E;
	nuc = nuc*rho;
	
        res += nuc;
        res *= Len;
      }
      return res; 
    }


  float MomLoss = 0;
  double pos_x,pos_y,pos_z;
  int i,j,k;

  if( !EL_ ) {
    goto end;
  }

  pos_x=hel(1);
  pos_y=hel(2);
  pos_z=hel(0);

  if( EL_->CoordSyst < 0 ) {
    
    if( pos_z <= EL_->minY || pos_z >= EL_->maxY || !EL_ ) goto end;

    double PHI,RAD,SIN,COS;
    if( pos_x==0 && pos_y==0 ){
      RAD=0.000000001;
      COS=1;
      SIN=0;
    }else{
      RAD = sqrt(pos_x*pos_x + pos_y*pos_y);
      if(RAD > EL_->maxX) goto end;
      
      COS = pos_x/RAD;
      SIN = pos_y/RAD;
    }

    if( EL_->nz > 1 ) {
     i=1;
     if( SIN < 0){
       while ( i<EL_->nz && EL_->Sin[i] < 0 && COS > EL_->Cos[i] ) {i++;}
     }else {
       while ( i<EL_->nz && EL_->Sin[i] <= 0 ) {i++;}
       while ( i<EL_->nz && COS <= EL_->Cos[i] ) {i++;}
     }
     i=i-1;
    }
    else i = 0;
    
    j=1;
    while ( j < EL_->nx && RAD > EL_->x[j] ) {j++;}
    j=j-1;

    k=1;
    while ( k < EL_->ny && pos_z > EL_->y[k] ) {k++;}
    k=k-1;

    MomLoss = EL_->Val[ i * (EL_->nx)*(EL_->ny) + j * (EL_->ny) + k ];
  }
  else {
    if( pos_z <= EL_->minY || pos_z >= EL_->maxY || !EL_ ) goto end;
    if( pos_x <= EL_->minX || pos_x >= EL_->maxX ) goto end;
    if( pos_y <= EL_->minY || pos_y >= EL_->maxY ) goto end;

    i=1;
    while ( i < EL_->nx && pos_x > EL_->x[i] ) {i++;}
    i=i-1;

    j=1;
    while ( j < EL_->ny && pos_y > EL_->y[j] ) {j++;}
    j=j-1;

    k=1;
    while ( k < EL_->nz && pos_z > EL_->z[k] ) {k++;}
    k=k-1;

    MomLoss = EL_->Val[ k * (EL_->nx)*(EL_->ny) + i * (EL_->ny) + j ];
    
  }

 end:
  //cout << "coral " << MomLoss*Len << endl;

  return MomLoss*Len;
}

float CsMaterialMap::getdEStraggling(const THlx& hel, float Len )
{
  if (usingROOTGeometry_ 
#if USE_TGEANT    
    || usingROOTGeometryTGEANT_
#endif  
     )
    
    {
      // The geometry uses the COMGEANT coordinate system.
      if (usingROOTGeometry_)
      gGeoManager->InitTrack(hel(0), hel(1), hel(2),
                             hel.DirCos(1), hel.DirCos(2), hel.DirCos(3));
#if USE_TGEANT
      else if (usingROOTGeometryTGEANT_)
      gGeoManager->InitTrack(hel(2), hel(0), hel(1),
                             hel.DirCos(3), hel.DirCos(1), hel.DirCos(2));
#endif       
      // Find current material.
      const TGeoMaterial *cmat = gGeoManager->GetCurrentVolume()->GetMedium()->GetMaterial();

      // Special case for the vacuum which has IMATE == 16.  This ends
      // up being the material's UniqueID after the translation to a
      // ROOT geometry by g2root.
      if (cmat->GetUniqueID() == 16)
        return 0;

      double A = cmat->GetA();
      double Z = cmat->GetZ();
      double rho = cmat->GetDensity();
      double p = fabs(1/hel(5));
      double E = hypot(p, massDefault_);
      double beta = p / E;
      double gamma = E / massDefault_;
      double betagamma = p / massDefault_;

      const double me = .51099891e-3;

      double ksix = 0.15353747e-3 * Z/A * rho/beta/beta;
      double Tmax = 2.*me*betagamma*betagamma/(1.+ 2. * gamma * me/massDefault_ + me*me / massDefault_ / massDefault_);

// use gaussian (as approx to Landau :-() to get the sigma of the energy loss distrib
//        double straggling = ksix * Len * Tmax * (1. - beta*beta/2.);
// use Landau
	double straggling = ksix * Len * 4.018;


//      cout << "straggling " << straggling << " p " << p << endl;

      return straggling*straggling/4.;
      }

      return 0;
}

vector<float> CsMaterialMap::getZoneBorders()
{
  assert (!usingROOTGeometry_);
  vector<float> Zones;

  vector<MatMap*>::iterator iz;
  for( iz=RL_.begin(); iz!=RL_.end(); iz++ ) {
    if( (*iz)->CoordSyst < 0 ) {
      Zones.push_back( 10 * (*iz)->minY );
      Zones.push_back( 10 * (*iz)->maxY );
    }
    else {
      Zones.push_back( 10 * (*iz)->minZ );
      Zones.push_back( 10 * (*iz)->maxZ );
    }
  }

  return Zones;
}

