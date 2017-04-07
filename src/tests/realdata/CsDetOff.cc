// $Id: CsDetOff.cc,v 1.18 2010/03/16 08:33:04 suhl Exp $
#include <string>
#include <vector>
#include <list>
#include <iostream>
#include <sstream>
#include <fstream>
#include <cstdlib>

#include "Opt.h"

#include <sys/time.h>
#include <sys/resource.h>
#include <unistd.h>
#include <stdio.h>

using namespace std;

//______________________________________________________________________
int main( int argc, char *argv[] )
{
  //=== First load Arguments
  if( argc < 4 ){
    cout << "Usage: DetOff <ExecProg> <optionFile> <DetSelection1> <DetSelection2> <...>" << endl;
    return 0;
  }
    
  list< string > TBNSels;
  list< string > TBNames;
  string execFile( argv[1] );
  string optFile( argv[2] );
  list< string > optOutFiles;

  cout << "DetOff - INFO: Option file      : \"" << optFile << "\".\n";
  for( unsigned int i = 3; i< (unsigned int) argc; i++ ) {
    TBNSels.push_back(string( argv[i] )); 
    cout << "DetOff - INFO: Detector name    : \"" << TBNSels.back() << "\".\n";
  }
  
  //=== Scan Option File to look for Detector.dat File
  Opt* opt = Opt::Instance( argc, argv );
  string detFile;
  if( !opt->getOpt("detector","table", detFile ) ) {
    cout << "DetOff - FATAL: No detector table found.\n";
    return 0;
  } else cout << "DetOff - INFO: Detector Table   : \"" << detFile << "\".\n\n";

  //=== Scan Detector File to look for detectors to be desactivated
  ifstream f( detFile.c_str(), ios::in );
  if( !f ) {
    cout << "DetOff - FATAL: Can't open detector table.\n";
    return 0;
  }
   

  const int lineSize = 512; 
  char line[lineSize];
  do {
    f.getline( line, lineSize, '\n' );
    if( f.eof() ) continue;
    if( line[0] == '\0' || line[0] != ' ' ) continue;
    if( line[1] == ' ' || line[1] == '-' || line[1] == '=' ) continue;

    // something good in this line...
    istringstream s(line);
    string opt; s >> opt;
    if( opt != "det" ) continue;
    
    int iDet;   s >> iDet;
    string TBName; s >> TBName;
    for( list<string>::iterator Is = TBNSels.begin(); Is != TBNSels.end(); Is++ ) {
      if( TBName.size() < (*Is).size() ) continue;
      bool accept = true;
      for( unsigned int j=0; j < (*Is).size(); j++ ) {
        if( (*Is)[j]!='_' && (*Is)[j] != TBName[j] ) {  
          accept = false;
          break;
        }
      }
      if( accept ){
        
        TBNames.push_back( TBName );
        cout << "DetOff - INFO: \"" << TBName << "\" matches.\n"; 
      }
    }
  } while( !f.eof() );
  f.close();
  cout << "\n";
  
  TBNames.unique();
  if( TBNames.empty() ) {
    cout << "DetOff - FATAL: no match found.\n"; 
    return 0;
  }
    
  //=== Get Traffic Dead Detectors
  list< string > TrafDetOff;
  bool DetOffFound = opt->getOpt( "TraF","DetNameOff", TrafDetOff );
  //=== Remove Redoundant DetNames
  for( list<string>::iterator It = TrafDetOff.begin(); It != TrafDetOff.end(); It++ ) 
  for( list<string>::iterator Is = TBNames.begin(); Is != TBNames.end(); Is++ ) {
    string TB((*Is),0, (*It).size() );
    if( TB == (*It) ) {
      cout << "DetOff - INFO: \"" << (*Is) << "\" removed from the list, as already off in option file.\n"; 
      TBNames.erase( Is );
      Is--;
      continue;
    }
  }
  
  //=== Get Histogram File
  string hFile;
  if( ! opt->getOpt("histograms","home", hFile ) ){
    cout << "DetOff - FATAL: Could not find \"histograms home\" key in option file.\n";
    return 0;
  } else cout << "\nDetOff - INFO: Histogram file : \"" << hFile << "\".\n";
  
  //=== Get DataFile and Run Number
  string DataFile;
  if( ! opt->getOpt("Data","file", DataFile ) ) cout << "DetOff - WARNING: Could not find \"Data File\" key in option file.\n";
  int iStop1 = DataFile.find( ".dat" ); cout <<"DetOff - INFO: iStop1 " << iStop1 << endl;
  int iStop2 = DataFile.find( ".dmp" ); cout <<"DetOff - INFO: iStop2 " << iStop2 << endl;
  int iStop = (iStop1 >= 0 && iStop1<iStop2) ? iStop1:iStop2;
  int iStart = DataFile.rfind( "-", iStop );
  string optRNumber("");
  if( iStop != 0 && iStart != 0 ){
    optRNumber = string(DataFile,iStart+1, iStop-iStart-1);
    cout << "DetOff - INFO: Last run number: " << optRNumber << "\n";
  } else cout << "DetOff - WARNING: Trouble looking into Data File name.\n";

  string Dir = "";
  if( optRNumber.size() ) {
    Dir = optRNumber+"/";
    if( access( Dir.c_str(), R_OK ) != 0 ) {
      string Task = string("mkdir ")+Dir;
      system( Task.c_str() );
      cout << "DetOff - INFO - Directory : \"" << Dir << "\".\n";
    }  
  }  

  //=== Get Number of Events to read
  string optNEvt;
  if( ! opt->getOpt("events","to read", optNEvt ) ) cout << "DetOff - WARNING: Could not find \"number of events\" key in option file.\n";
  cout << "DetOff - INFO: #events        : " << optNEvt << "\n";
      
  //=== Generate Option files
  for( list<string>::iterator Is = TBNames.begin(); Is != TBNames.end(); Is++ ) {
    string optOutFile( 
      Dir +
      ( (optRNumber.size() ) ? 
      optFile+"."+"run"+optRNumber+"."+(*Is):
      optFile+"."+(*Is) ) );
    cout << "DetOff - INFO: building \"" << optOutFile << "\"."; 
    if( access( optOutFile.c_str(), R_OK ) == 0 ) {
      cout << " ERROR -> file already exists.\n";
      optOutFiles.push_back( optOutFile );
      continue;
    }

    //=== Scan Detector File to look for detectors to be desactivated
    ifstream in(  optFile.c_str(), ios::in );
    ofstream out( optOutFile.c_str(), ios::out );
    if( !in ) return 0;
    if( !out ) {
      cout << " ERROR -> can't open.\n";
      continue;
    }

    //=== Print Header
    out << "// Automaticaly generated File using command <";
    for( unsigned int i = 0; i < (unsigned int) argc; i++ ) out << " " << argv[i];
    out << ">\n";
    
    bool hFileFound = false;
    bool detOffFound = false;
      
    char l[lineSize];
    do {
    
      in.getline( l, lineSize, '\n' );
      
      istringstream s(l);
      
      //=== Check if line is "histograms","home"
      string tag; s >> tag;
      if( s.good() && tag == "histograms" ) {
        string key; s >> key; 
        if( s.good() && key == "home" ){

          //=== Generate Histogram output file name
          string hOutFile( hFile );
          if( optRNumber.size() ) hOutFile+=(".run"+optRNumber);
          if( optNEvt.size() )    hOutFile+=("."+optNEvt);
          hOutFile+=("."+(*Is));
          
          cout << " HFile is \"" << hOutFile << "\".";
          out << "histograms home " << hOutFile << endl;
          hFileFound = true;
        } else out << l << endl;
        
      //=== Check if line is "TraF","DetNameOff"
      } else if( s.good() && tag == "TraF" ) {
        string key; s >> key; 
        if( s.good() && key == "DetNameOff" ){
          out << "TraF DetNameOff";
          for( list<string>::iterator It = TrafDetOff.begin(); It != TrafDetOff.end(); It++ )  out << " " << (*It);
          
          out << " ";
          
          //=== remove underscores from the TBName <- needed for Traffic
          for( unsigned int i = 0; i< ( (*Is).size() ); i++ ) 
          if( (*Is)[i] != '_' ) out << (*Is)[i];
          detOffFound = true;
          out << endl;  
        } else out << l << endl;
        
      //=== copy line as is  
      } else out << l << endl; 
    } while( !in.eof() );
    in.close();
    
    if( !hFileFound ) {
      cout << " WARNING -> could not find hFile.\n";

      string hOutFile( hFile+"."+(*Is) );
      cout << " HFile is \"" << hOutFile << "\".";
      out << "histograms home " << hOutFile << endl;
    }
    
    if( !detOffFound ){
      cout << " WARNING -> could not find DetNameOff.\n";
      out << "TraF DetNameOff ";
      for( unsigned int i = 0; i< ( (*Is).size() ); i++ ) 
      if( (*Is)[i] != '_' ) out << (*Is)[i];
      detOffFound = true;
      out << endl;  
    }
    
    out.close();
    optOutFiles.push_back( optOutFile );
    cout << " Done\n";
  }
    
  //=== Generate the script
  string outScript;
  outScript = "DetOff";
  for( list<string>::iterator It = TBNSels.begin(); It != TBNSels.end(); It++ ) outScript+="_"+(*It);
  outScript+=( (optRNumber.size() ) ? ".run"+optRNumber+".csh":".csh" );
  
  cout << "\nDetOff - INFO: building script file: \"" << outScript << "\".";
  if( access( outScript.c_str(), R_OK ) == 0 ) {
    cout << " FATAL -> file already exists.\n";
    return 0;
  }
  
  ofstream out( outScript.c_str(), ios::out );
  if( !out ) {
    cout << " FATAL -> can't open.\n";
    return 0;
  }
  

  //=== Print Header
  out << "#Automaticaly generated file using command <";
  for( unsigned int i = 0; i < (unsigned int) argc; i++ ) out << " " << argv[i];
  out << ">\n";
  
  //=== scan optOut
  unsigned int i = 1;
  for( list<string>::iterator Io = optOutFiles.begin(); Io != optOutFiles.end(); Io++ ) {
    out << "set Log=$HOME/coral.`date +%d_%m_%Y.%H_%M_%S`.log\n";
    out << "bsub -R \"select[cpuf>=1.5]\" -q8nh -J " 
        << (*Io) 
        << ".`date +%d_%m_%Y.%H_%M_%S` -o $Log `pwd`/" 
        << execFile 
        << " `pwd`/" 
        << (*Io) 
        << "; sleep 5\n"; 
    i++;
  }
  out.close(); 
  cout << " Done.\n";
  return 0;
}
   
