// $Id: DBShift.cc,v 1.2 2002/08/17 20:54:26 hpereira Exp $
#include <string>
#include <iostream>
#include <strstream>
#include <fstream>
#include <vector>
#include <stdio.h>

//______________________________________________________________________
int main( int argc, char *argv[] )
{
  //=== First load Arguments
  cout << "Usage: DBShift <new_t0 (ns)> <inDBFileSelection>\n";
  if( argc < 3 ) return 0;
    
  float newT; istrstream( (const char*) argv[1] ) >>  newT;
  vector< string > inFile;
  for( int i=2; i< argc; i++ ) inFile.push_back( string( argv[i] ) );
 
  //=== Process all input files ===
  for( unsigned int i = 0; i< inFile.size(); i++ ) {
    
    //=== open input file
    ifstream in( inFile[i].c_str(), ios::in );
    if( !in ) {
      cout << "DBShift - ERROR: Can't open input DB file \"" << inFile[i] << "\".\n";
      continue;
    }
   
    //=== open output file
    string outFile = inFile[i]+".out";
    ofstream out( outFile.c_str(), ios::out );
    if( !out ) {
      cout << "DBShift - ERROR: Can't open output DB file \"" << outFile << "\".\n";
      in.close();
      continue;
    }
  
    cout << "DBShift - INFO: Updating \"" << inFile[i] << "\" [NewT0 = " << newT << " (ns)].\n";
  
    //=== Initialize oldTime
    double oldT;
    bool oldTSet = false;
    bool RTTypeFound = false;
    string RTType; 
    
    //=== Parse input file
    const int lineSize = 512; 
    char line[lineSize];
    do {
      in.getline( line, lineSize, '\n' );
      if( in.eof() ) continue;
      if( line[0] == '\0' ) {
        out << line;
        continue;
      }
      
      //=== Read line
      istrstream s(line);
      
      //=== Try getting old T0 value
      if( !oldTSet ) {
        s >> oldT;
        if( !s.good() ) {
          cout << "DBShift - FATAL: Can't read old T0.\n";
          in.close(); out.close();
          return 0;
        }
          
        cout << "DBShift - INFO: OldT0 = " << oldT << " (ns).\n";
        if( newT > oldT ) cout << "DBShift - WARNING: old<new.\n";
        oldTSet = true;
        out.form("%9.1f\n", newT );
        continue;
      }
      
      //=== Try looking for RTKeyWord
      if( !RTTypeFound ) {
        s >> RTType;
        if( ( RTType == "RTFit1" || RTType == "RTFit2" || RTType == "RT" ) && oldTSet ) {
          cout << "DBShift - INFO: found RT Relation of type \"" << RTType << "\".\n";
          out << RTType;
          double p0; s >> p0;
          if( ! s.good() ) {
            cout << "DBShift - FATAL: Wrong format for RTRelation.\n";
            in.close(); out.close();
            return 0;
          }
          out.form( "%12.4e", p0+oldT-newT);  
          while( s.good() ) {
            double p; s >> p;
            if( ! s.good() ) continue;
            out.form( "%12.4e", p);  
          }
          out << endl;
          
          RTTypeFound = true;
          continue;
        } else if( RTType == "RTGrid" && oldTSet ) {
          cout << "DBShift - INFO: found RT Relation of type \"" << RTType << "\".\n";
          out << line << endl;
          RTTypeFound = true;
        } else out << line << endl;
      } else if( RTType == "RTGrid" ) {
        double t;  s >> t;
        double r;  s >> r;
        double dr; s >> dr;
        if( !s.good() ) {
          out << line << endl;
          RTTypeFound = false;
        } else out.form( "%10.4e %10.4e %10.4e\n", t+oldT-newT, r, dr );
      }
      
    } while( !in.eof() );
    
    in.close();
    out.close();
    string command = string("mv ")+outFile+" "+inFile[i];
    system( command.c_str() );
  }
    
  cout << "DBShift - Done.\n";
  
  return 0;
}
      
