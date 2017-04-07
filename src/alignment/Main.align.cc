// $Id: Main.align.cc,v 1.14 2009/07/21 13:30:19 aaust Exp $
/*!
   \file    Main.align.cc
   \brief   Instanciate a Align object, performs alignment minimisation using option file
   \author  Hugo Pereira
   \version $Revision: 1.14 $
   \date    $Date: 2009/07/21 13:30:19 $
*/
#include "Tracks.h"
#include "Align.h"
#include "Defs.h"
#include "Opt.h"


#include <TROOT.h>

#include <stdio.h>
#include <unistd.h>
#include <iostream>
#include <fstream>
#include <list>
#include <string>
using namespace std;
//! TROOT singleton
TROOT root("align", "global alignment using Millepede");

/*! 
  \fn int main( int argc, char *argv[] )
  \brief arguments are [-b] (batch mode) <option file>
*/
//!___________________________________________________
int main( int argc, char *argv[] )
{
  
  //! Parse command line arguments - Initialise Opt
  if( argc < 2 ) {
    cout << "Usage: align <OptionFile>\n";
    return 0;
  }

  Opt *opt = Opt::Instance(argc, argv);

  list<DetectorInfo*> detList;
  vector<DetectorInfo*> dets;
  
  //! Create alignment object
  Align al;
  
  //!= Check if job is batch job
  for( int i=1; i<argc; i++ )
  if( string(argv[i])=="-b" ) {
    cout << "align - INFO: this is a batch job.\n";
    al.isBatch_=true;
  }
  
  //! Get Alignment options 
  bool magnets_on = opt->getOpt( "main", "magnets on" ); 
  al.AlignU(opt->getOpt( "align", "U" ));   //!< Offset perp to the wire
  al.AlignZ(opt->getOpt( "align", "Z" ));   //!< Offset along the beam
  al.AlignT(opt->getOpt( "align", "T" ));   //!< angular offset perp to the beam
  al.AlignP(opt->getOpt( "align", "P" ));   //!< pitch modification
  al.AlignR(opt->getOpt( "align", "R" ));   //!< distance to the wire (Offset on t0 for drift-like detectors)
  al.AlignL(opt->getOpt( "align", "L" ));   //!< Lorentz angle R scaling (for drift-like detectors)

  //! millepede configuration
  al.Iterate(opt->getOpt( "make", "iterations" ));
  
  //! get detector table
  string detectorTable;

  if( !opt->getOpt("detector","table", detectorTable ) ) {
    cout << "align - FATAL - could not read \"detector\" \"table\" in option file.\n";
    return 0;
  } 
  cout << "detector table is " << detectorTable << endl;
  opt->expand( detectorTable );
  al.LoadDetectorFile( detectorTable.c_str() );
  //! Set Excluded detectors
  list< string > excludeDets;
  string exclude("");
  opt->getOpt("align","excludeDets",excludeDets );
  cout << "exclude some dets... " << endl;
  for( list<string>::iterator I=excludeDets.begin(); I!=excludeDets.end(); I++ )
  exclude+=(*I)+" ";
  al.ExcludeDetectors( exclude.c_str() );
  
  //! Set Used Detectors
  list< string > useDets;
  string use("");
  opt->getOpt("align","useDets", useDets );
  cout << "include some dets... " << endl;
  for( list<string>::iterator I=useDets.begin(); I!=useDets.end(); I++ )
  use+=(*I)+" ";
  al.UseDetectors( use.c_str() );
  al.SortDetectors( );

  //! Set OuterST option
  string aost;
  bool alignOuterST = opt->getOpt("align", "OuterST", aost); 
  if( alignOuterST && aost != "YES" ) alignOuterST=false;
  else if (alignOuterST) cout << "Align::AlignOuterST YES" << endl;
  al.AlignOuterST(alignOuterST);

  //! Change resolutions  
  list< string > resL;
  while( opt->getOptRec( "align","change resolution", resL ) ) {
    if( resL.size()<2 ) continue;
    string selection;
    double res;
    int i=0;
    for( list< string >::iterator I=resL.begin(); I!=resL.end(); I++, i++ ) {
      istringstream s( (*I).c_str() );
      switch (i) {
        case 0:  selection = (*I);  break;
        case 1:  s >> res; break;
        default: break;
      }
    }
    
    al.ChangeResolution( selection.c_str(), res );
  }

  //! Dump detectors if needed
  if( opt->getOpt( "debug", "dump dets" ) ) al.DumpDetectors("*");
  
//   //! Get Solenoid info
//   vector<double> solInfo;  //!< solenoid information ( zmin[mm], zmax[mm], bz[T], [scale] ) 
//   bool hasSolInfo = false;
//   if( ( hasSolInfo = opt->getOpt( "solenoid", "info", solInfo ) ) && solInfo.size() > 2 ) {
//     if( !magnets_on ) al.SetSolenoidInfo( solInfo[0], solInfo[1], solInfo[2], (solInfo.size() > 3) ? solInfo[3]:1 );
//     else { cout << "align - WARNING: magnets on. Solenoid info ignored.\n"; hasSolInfo = false; }
//   } else cout << "align - WARNING: no solenoid info available.\n";
//   
//   //! book transfer matrix
//   double poc=0;  //!< momentum over charge [GeV] for ref particles
//   if( hasSolInfo && !( opt->getOpt( "particle", "info", poc ) && al.BuildTransferMatrix(poc) ) ) 
//   cout << "align - WARNING: troubles building solenoid transfer matrix. solenoid info ignored.\n";
    
  //! initialize millepede 
  int stdDev;
  if( !opt->getOpt( "align", "nStdDev", stdDev ) ) stdDev=1;

  al.InitParameters( stdDev, opt->getOpt("debug", "dump millepede" ) ); 
  
  //! Fix U
  list< string > fixDets;
  string fix("");
  opt->getOpt("fix","U", fixDets );
  for( list<string>::iterator I=fixDets.begin(); I!=fixDets.end(); I++ )
  fix+=(*I)+" ";
  al.FixU( fix.c_str() );
  
  //! Fix Z
  fixDets.clear();
  fix="";
  opt->getOpt("fix","Z", fixDets );
  for( list<string>::iterator I=fixDets.begin(); I!=fixDets.end(); I++ )
  fix+=(*I)+" ";
  al.FixZ( fix.c_str() );
  
  //! Fix T
  fixDets.clear();
  fix="";
  opt->getOpt("fix","T", fixDets );
  for( list<string>::iterator I=fixDets.begin(); I!=fixDets.end(); I++ )
  fix+=(*I)+" ";
  al.FixT( fix.c_str() );
  
  //! Fix P
  fixDets.clear();
  fix="";
  opt->getOpt("fix","P", fixDets );
  for( list<string>::iterator I=fixDets.begin(); I!=fixDets.end(); I++ )
  fix+=(*I)+" ";
  al.FixP( fix.c_str() );
   
  //! Fix R
  fixDets.clear();
  fix="";
  opt->getOpt("fix","R", fixDets );
  for( list<string>::iterator I=fixDets.begin(); I!=fixDets.end(); I++ )
  fix+=(*I)+" ";
  al.FixR( fix.c_str() );
    
  //! Fix L
  fixDets.clear();
  fix="";
  opt->getOpt("fix","L", fixDets );
  for( list<string>::iterator I=fixDets.begin(); I!=fixDets.end(); I++ )
  fix+=(*I)+" ";
  al.FixL( fix.c_str() );

  //! Add Bias
  list< string > biasL;
  while( opt->getOptRec( "align", "add bias", biasL ) ) {
    string bias("");
    string TBName("");
    int i=0;
    for( list< string >::iterator I=biasL.begin(); I!=biasL.end(); I++, i++ )
    switch (i) {
      case 0:  TBName = (*I);  break;
      default: bias+=(*I)+" "; break;
    }
    al.SetBias( TBName.c_str(), bias.c_str() );
  }
     
  //! Load Tracks
  string file;
  while( opt->getOptRec( "input","tree", file ) )
  if( !al.GetTracks() ) al.LoadTracks( file.c_str(), magnets_on );
  else al.GetTracks()->AddToChain( file.c_str() );
    
  //! Load cut
  string selectionCut;
  if( opt->getOpt("selection", "cut", selectionCut ) ) al.AddCut( selectionCut.c_str() );
  
  //! Perform minimisation
  int nEnt;    if( !opt->getOpt("align", "nTracks", nEnt ) )         nEnt = 0;
  int refresh; if( !opt->getOpt("debug", "counter rate", refresh ) ) refresh=0;
  // new part - seb
  list< string > reparamDets;
  string reparam("");
  opt->getOpt("reparam","Z",reparamDets);
  for( list<string>::iterator I=reparamDets.begin(); I!=reparamDets.end(); I++ )
    reparam+=(*I)+" ";
  // end of new part -seb
  al.Minimize( (unsigned int) nEnt, (unsigned int) refresh );
  
  //! Dump Result to files
  string outFile;
  if( opt->getOpt("output","file", outFile ) ) opt->expand(outFile);
  else outFile = "align.out";
  cout << "align - INFO: results dump to \"" << outFile << "\".\n";
  
  al.DumpToFile( outFile.c_str(), reparam.c_str() );

}
