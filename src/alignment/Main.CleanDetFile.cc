// $Id: Main.CleanDetFile.cc,v 1.4 2006/06/16 15:21:47 conrad Exp $
/*!
   \file    Main.CleanDetFile.cc
   \brief   Instanciate a DetFileManager object, parse detectors, executes all macro in CleanDetFile.cc
   \author  Hugo Pereira
   \version $Revision: 1.4 $
   \date    $Date: 2006/06/16 15:21:47 $
*/

#include <iostream>
#include "DetFileManager.h"

using namespace std;

/*! 
  \fn int main( int argc, char *argv[] )
  \brief arguments are <detector.dat> [<new detector.dat>]
*/
//_______________________________________________________________________________
int main( int argc, char *argv[] )
{
  if( argc < 2 ) {
    cout << "Usage: cleanDetFile <detector.dat> [<new detector.dat>]" << endl;
    return 0;
  }
  
  string detFile(argv[1]);
  
  DetFileManager df( detFile.c_str() );
  df.CleanDetFile();
  
  if( argc > 2 ) df.DumpToFile( argv[2] );
  else df.DumpToFile( );
  
  return 0;
}
