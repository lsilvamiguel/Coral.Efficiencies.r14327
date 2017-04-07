// $Id: CsManyChunk.cc,v 1.8 2010/03/16 08:33:04 suhl Exp $
#include <string>
#include <vector>
#include <list>
#include <iostream>
#include <sstream>
#include <fstream>
#include <cstdlib>

#include "Opt.h"
#include "fdstream.hpp"

#include <sys/time.h>
#include <sys/resource.h>
#include <unistd.h>
#include <stdio.h>
#include <pwd.h>

using namespace std;

#define year "2002"

string Path;
string Run;
string RunNumber;
string Chunk;

//______________________________________________________________________
void ParseDataFile( string DataFile )
{
  //=== Get Path
  int iStop = DataFile.rfind("/");
  Path = DataFile.substr( 0, iStop );
  
  //=== Get chunk ID
  iStop  = DataFile.rfind("-");
  int iStart = DataFile.rfind("cdr");
  Chunk = DataFile.substr( iStart+strlen("cdr"), iStop-iStart-strlen("cdr") );
  
  //=== Get endof string
  Run = DataFile.substr( iStop+1, DataFile.size()-iStop );
  RunNumber = Run.substr( 0, Run.find( "." ) );
}


//______________________________________________________________________
int main( int argc, char *argv[] )
{
  //=== load current path
  char* buf = new char[256];
  string pwd("");
  if (!getcwd(buf, 256)) {
    cout << "ManyChunk - FATAL: getcwd failed\n";
    return 0;
  } else {
    pwd = string( buf );
    cout << "ManyChunk - INFO: pwd is \"" << pwd << "\".\n";
  }
  
  delete buf; 	
  
  //=== load coral path
  string coral;
  if( getenv("CORAL") ) {
    coral = string( getenv("CORAL") );
    cout << "ManyChunk - INFO: coral path is \"" << coral << "\".\n";
  } else {
    cout << "ManyChunk - FATAL: getenv(\"CORAL\") failed\n";
    return 0;
  }
  
  //=== load home path
  string home;
  if( getenv("HOME") ) {
    home = string( getenv("HOME") );
    cout << "ManyChunk - INFO: home path is \"" << home << "\".\n";
  } else {
    cout << "ManyChunk - FATAL: getenv(\"HOME\") failed\n";
    return 0;
  }
  
  //=== load Arguments
  if( argc < 3 ){
    cout << "Usage: ManyChunk [-r] <ExecProg> <optionFile> [<ChunkID1> <ChunkID2> <...>]" << endl;
    return 0;
  }
    
  unsigned int iArg=1;
  bool doReload = false;
  while( iArg < (unsigned int) argc && argv[iArg][0] == '-') {
    if( string( argv[iArg] ) == "-r" ) doReload = true;
    iArg++;
  }
  
  string execFile("");
  string optFile("");
  if( iArg < (unsigned int) argc ){ execFile = argv[iArg]; iArg++; }
  if( iArg < (unsigned int) argc ){ optFile  = argv[iArg]; iArg++; }
  list< string > optOutFiles;
  
  cout << "ManyChunk - INFO: Exec file        : \"" << execFile << "\".\n";
  cout << "ManyChunk - INFO: Option file      : \"" << optFile << "\".\n";
  if( doReload ) 
  cout << "ManyChunk - INFO: Reloading chunks required.\n";
   
  vector< string > Chunks;
  
  for( unsigned int i = iArg; i< (unsigned int) argc; i++ ) {
    Chunks.push_back(string( argv[i] )); 
    cout << "ManyChunk - INFO: add chunk : " << Chunks.back() << endl;
  }
    
  Opt *opt = Opt::Instance(argc, argv);
  
  //=== Get and parse Input Files
  unsigned int iDF = 0;
  string DataFile;   
  
  bool useDB = (!opt->getOptRec("Data", "file", DataFile ));
  if( useDB && !opt->getOptRec("Data run", "select", RunNumber )) {
    cout << "ManyChunk - FATAL: Could not find \"Data run select\" in option file.\n";
    return 0;
  } 
  
  string SavePath( "" );
  if( !useDB ) { 
    ParseDataFile( DataFile ); 
    SavePath = Path;
    cout << "ManyChunk - INFO - Data file : \"" << DataFile << "\".\n";
    cout << "ManyChunk - INFO - Path      : \"" << SavePath << "\".\n";
    cout << "ManyChunk - INFO - Chunk     : \"" << Chunk << "\".\n";
    cout << "ManyChunk - INFO - Run       : \"" << Run << "\".\n";
  }
  cout << "ManyChunk - INFO - RunNumber : \"" << RunNumber << "\".\n";
    
  //=== Do CastorReload if doReload
  vector< string > DBNames;
  string DBFile;
  if( doReload ) {
    string Task( string("castorFile ") + RunNumber + " " + year );
    FILE* fp;
    fp=popen(Task.c_str(),"r");
    boost::fdistream in(fileno(fp));
        
    const int lineSize = 512; 
    char line[lineSize];
    do {
      in.getline( line, lineSize, '\n' );
      if( line[0] == '\0' );
      istringstream s(line);
      string file ;  s >> file;
      string status; s >> status;
      if( status == "onTAPE" || status == "onDISK" ) {
        ParseDataFile( file );
        if( Chunk[Chunk.size()-1] == '0' ) continue;
        cout << "ManyChunk - INFO - got \"" << file << "\".\n";
        DBNames.push_back( file );
      } 
    } while( in );
  }

  //=== Create directory
  string Dir( RunNumber );
  if( access( Dir.c_str(), R_OK ) != 0 ) {
    string Task = string("mkdir ")+Dir;
    system( Task.c_str() );
    cout << "ManyChunk - INFO - Directory : \"" << Dir << "\".\n";
  }  
  const int lineSize = 512; 

  //===Check if all Chunks are Required - Get all chunks if yes
  if( Chunks.empty() ) {
    cout << "ManyChunk - INFO - All reloaded chunks required.\n";
    string Task = string( "nsls ") + Path + " | grep " + RunNumber;
    FILE* fp;
    fp=popen(Task.c_str(),"r");
    boost::fdistream in(fileno(fp));
    
    char* buf = new char[lineSize];
  	while( in.getline( buf, lineSize, '\n' ) ) {
      ParseDataFile( buf );
      
      //=== Reject X0000 chunks
      if( Chunk[Chunk.size()-1] == '0' ) continue;
      
      Chunks.push_back( Chunk );
    }
    delete buf;
  }
  
  //=== Get Histogram File
  string hFile;
  if( ! opt->getOpt("histograms","home", hFile ) ){
    cout << "DetOff - FATAL: Could not find \"histograms home\" key in option file.\n";
    return 0;
  }  
  
  //=== Create DataFiles accoring to chunks in Chunk
  vector< string > OptFiles;
  vector< string > ChunksOK;
  for( unsigned int i=0; i< Chunks.size(); i++ ) {
    string newFile("");
    string container("");
    if( useDB ){
      container = string("cdr")+Chunks[i]+"-"+RunNumber;
      cout << "ManyChunk - \"" << container << "\" -> ";
    } else {
      newFile = SavePath+"/"+"cdr"+Chunks[i]+"-"+RunNumber+".dmp";
      cout << "ManyChunk - \"" << newFile << "\" -> ";
    }
    
    string optOutFile( Dir+"/"+optFile + "." + Chunks[i] );
    cout << "\"" << optOutFile << "\" - ";

    //=== Check if chunk was found
    if( doReload ) {
      bool found = false;
      for( unsigned int j=0; j<DBNames.size(); j++ )
      if( DBNames[j].find( Chunks[i] ) != string::npos ) {
        found = true;
        break;
      }
    
      if( !found ) {
        cout << " ERROR -> chunk not found" << endl;
        continue;
      }
    }
    
    ifstream in(  optFile.c_str(), ios::in );
    ofstream out( optOutFile.c_str(), ios::out );
    if( !in ) return 0;
    if( !out ) {
      cout << " ERROR -> can't open";
      continue;
    }
    
    //=== Print Header
    out << "// Automaticaly generated File using command <";
    for( unsigned int j = 0; j < (unsigned int) argc; j++ ) out << " " << argv[j];
    out << ">\n";
    
    bool hFileFound = false;
    bool dataFileFound = false;
   
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
          string hOutFile = hFile;
          if( hFile.find( RunNumber ) == string::npos ) hOutFile += "."+RunNumber;
          hOutFile+="."+Chunks[i];
          out << "histograms home " << hOutFile << endl;
          hFileFound = true;
        } else out << l << endl;
        
      //=== Check if line is "Data","file"
      } else if( s.good() && tag == "Data" ) {
        string key; s >> key; 
        
        // CastorReload
        if( s.good() && (!useDB) && key == "file" ){
          out << "Data file " << newFile << endl;
          out << endl;  
          dataFileFound = true;
        
        // objectivity
        } else if( s.good() && useDB && key == "container" ) {
          out << "Data container " << container << endl;
          dataFileFound = true;
  
        //=== copy line as is  
        } else out << l << endl; 
      
      //=== copy line as is  
      } else out << l << endl;
        
    } while( !in.eof() );
    
    if( !dataFileFound ) {
      if( useDB ) out << "Data container " << container << endl;
      else out << "Data file " << newFile << endl;
    } 
    in.close();
    out.close();
    OptFiles.push_back( optOutFile );
    ChunksOK.push_back( Chunks[i] );
    cout << " Done.\n";
  }
    
  //=== Generate the script
  string outScript;
  outScript = "ManyChunk";
  outScript+= "."+Dir+".csh";
  
  cout << "\nDetOff - INFO: building script file: \"" << outScript << "\".\n";  
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
  for( unsigned int i=0; i<OptFiles.size(); i++ ) {
    
    cout << "ManyChunk - Adding \"" << OptFiles[i] << "\" - ";
        
    //=== if doReload, try getting ChunkOK in DBNames
    string DBName;
    if( doReload ) {
      bool found = false;
      for( unsigned int j=0; j<DBNames.size(); j++ )
      if( DBNames[j].find( ChunksOK[i] ) != string::npos ) {
        found = true;
        DBName = DBNames[j];
        break;
      }
  
      if( !found ) {
        cout << " ERROR -> chunk not found.\n"; 
        continue;
      }
    }  
    string outBSubScript;
    outBSubScript = Dir+"/bsub";
    outBSubScript+= "."+optFile + "."+ChunksOK[i];
    
    cout << "\"" << outBSubScript << "\" - ";  
    ofstream subout( outBSubScript.c_str(), ios::out );
    if( !subout ) {
      cout << " FATAL -> can't open.\n";
      continue;
    }
    
    //=== Print Header
    subout << "#Automaticaly generated file using command <";
    for( unsigned int j = 0; j < (unsigned int) argc; j++ ) subout << " " << argv[j];
    subout << ">\n";
    
    //=== Print BSub options
    string log = home + "/log." + optFile + "." + ChunksOK[i];
    subout << "#BSUB -q 8nh" << endl;
    subout << "#BSUB -L /bin/csh" << endl;
    subout << "#BSUB -R \"select[cpuf>=1.5]\" "  << endl;     
    subout << "#BSUB -J " << OptFiles[i] << endl;
    subout << "#BSUB -o " << log << endl;
    subout << endl;
    
    //=== Print execs
    subout << "cd " << coral << "; source setup.sh" << endl;
    subout << "cd " << pwd << "; " << execFile << " "  << OptFiles[i] << endl;
    subout.close();
    
    if( doReload && !useDB ) out << endl << "/afs/cern.ch/user/o/objsrvvy/castorFile/castorReload " << DBName << endl;
    out << "rm -f " << log << "; bsub < " << outBSubScript << "; sleep 5" << endl;
    cout << " Done.\n";
      
  }
  out.close(); 
  
  //=== Make output file executable
  string Task = string("chmod u+x ")+outScript;
  system( Task.c_str() );
  
  return 0;
}
