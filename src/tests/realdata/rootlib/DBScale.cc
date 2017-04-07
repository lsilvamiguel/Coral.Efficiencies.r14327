// $Id: DBScale.cc,v 1.3 2002/05/23 11:49:46 hpereira Exp $
#include <string>
#include <iostream>
#include <strstream>
#include <fstream>
#include <vector>
#include <stdio.h>
#include <math.h>

//______________________________________________________________________
int main( int argc, char *argv[] )
{
  //=== First load Arguments
  cout << "Usage: DBScale <scale> <inDBFileSelection>\n";
  if( argc < 3 ) return 0;
    
  double scale; istrstream( (const char*) argv[1] ) >>  scale;
  vector< string > inFile;
  for( int i=2; i< argc; i++ ) inFile.push_back( string( argv[i] ) );
 
  //=== Process all input files ===
  for( unsigned int i = 0; i< inFile.size(); i++ ) {
    
    //=== open input file
    ifstream in( inFile[i].c_str(), ios::in );
    if( !in ) {
      cout << "DBScale - ERROR: Can't open input DB file \"" << inFile[i] << "\".\n";
      continue;
    }
   
    //=== open output file
    string outFile = inFile[i]+".out";
    ofstream out( outFile.c_str(), ios::out );
    if( !out ) {
      cout << "DBScale - ERROR: Can't open output DB file \"" << outFile << "\".\n";
      in.close();
      continue;
    }
  
    cout << "DBScale - INFO: Updating \"" << inFile[i] << "\"\n[";
  
    //=== Initialization
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
          cout << "ERROR: Can't read old T0].\n";
          in.close(); out.close();
          continue;
        }
          
        cout << "INFO: OldT0 = " << oldT << "][";
        oldTSet = true;
        double newT = oldT*scale;
        out.form("%9.1f\n", newT );
        cout << "INFO: NewT0 = " << newT << "][";
        continue;
      }
      
      //=== Try looking for RTKeyWord
      if( !RTTypeFound ) {
        s >> RTType;
        if( ( RTType == "RTFit1" || RTType == "RTFit2" || RTType == "RT" ) && oldTSet ) {
          cout << "INFO: found \"" << RTType << "\"][";
          out << RTType;
          double p0; s >> p0;
          if( ! s.good() ) {
            cout << "ERROR: Wrong format].\n";
            in.close(); out.close();
            continue;
          }
          out.form( "%12.4e", p0*scale);  
          
          unsigned int pw = 0;
          while( s.good() ) {
            double p; s >> p;
            pw++;
            if( ! s.good() ) continue;
            out.form( "%12.4e", p/pow( scale, double( pw ) ) );  
          }
          cout << "OK].\n";
          out << endl;
          
          RTTypeFound = true;
          continue;
        } else if( RTType == "RTGrid" && oldTSet ) {
          cout << "INFO: found \"" << RTType << "\"][";
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
        } else out.form( "%10.4e %10.4e %10.4e\n", t*scale, r, dr );
      }
    } while( !in.eof() );
    
    if(  RTType == "RTGrid" ) cout << "OK].\n";
    
    in.close();
    out.close();
    string command = string("mv ")+outFile+" "+inFile[i];
    system( command.c_str() );
  }
    
  cout << "DBScale - Done.\n";
  
  return 0;
}
      
