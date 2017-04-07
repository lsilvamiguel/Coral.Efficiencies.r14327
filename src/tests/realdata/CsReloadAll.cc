#include <string>
#include <vector>
#include <iostream>
#include <sstream>

#include <stdio.h>

#include "fdstream.hpp"

using namespace std;

string Path;
string Run;
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
}

//______________________________________________________________________
int main( int argc, char *argv[] )
{
  //=== First load Arguments
  if( argc < 3 ){
    cout << "Usage: ReloadAll <runNumber> <year> [Chunk1] [Chunk2] ..." << endl;
    return 0;
  }
  
  string runNumber( argv[1] );
  string year( argv[2] );
  vector< string > Chunks;
  for( int iArg=3; argc > iArg; iArg++ ) Chunks.push_back( argv[iArg] ); 
  
  vector< string > rawNames;
  string Task( string("castorFile ") + runNumber + " " + year );
  FILE* fp;
  fp=popen(Task.c_str(),"r");
  boost::fdistream in(fileno(fp));

  const int lineSize = 512; 
  char line[lineSize];
  do {
    in.getline( line, lineSize, '\n' );
    if( line[0] == '\0' );
    istringstream s(line);
    string file; s >> file;
    string status; s >> status;
    if( status == "onTAPE" || status == "onDISK" ) {
      ParseDataFile( file );
      if( Chunk[Chunk.size()-1] == '0' ) continue;
      
      bool accepted = true;
      if( Chunks.size() ) {
        accepted=false;
        for( unsigned int i=0; i< Chunks.size(); i++ )
        if( Chunk == Chunks[i] ) { accepted = true; break; }
      }
      if( !accepted ) continue;
      
      rawNames.push_back( file );
      cout << "ReloadAll - INFO - got \"" << file << "\".\n";
    } 
  } while( in );
  
  for( unsigned int i = 0; i < rawNames.size(); i++ ) {
    string Task( string("castorReload ") + rawNames[i] );
    FILE* fp;
    fp=popen(Task.c_str(),"r");
    boost::fdistream in(fileno(fp));

    const int lineSize = 512; 
    char line[lineSize];
    do {
      in.getline( line, lineSize, '\n' );
      cout << line << endl;
    } while( in );
  
  }
  
  return 0;
}
