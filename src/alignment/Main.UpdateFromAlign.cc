// $Id: Main.UpdateFromAlign.cc,v 1.5 2006/06/16 15:21:47 conrad Exp $
/*!
   \file    Main.UpdateFromAlign.cc
   \brief   Instanciate a DetFileManager object, update and save det.dat file from alignment output file
   \author  Hugo Pereira
   \version $Revision: 1.5 $
   \date    $Date: 2006/06/16 15:21:47 $
*/

#include <iostream>
#include "DetFileManager.h"

using namespace std;

/*! 
  \fn int main( int argc, char *argv[] )
  \brief arguments are <alignment output file> <detector.dat> [<new detector.dat>]
*/
//_______________________________________________________________________________
int main( int argc, char *argv[] )
{
  if( argc < 3 ) {
    cout << "Usage: updateFromAlign <alignment output file> <detector.dat> [<new detector.dat>]" << endl;
    return 0;
  }  
  
  string alignFile(argv[1]);
  string detFile(argv[2]);
  
  DetFileManager df( detFile.c_str() );
  df.UpdateFromAlign( alignFile.c_str() );
  df.CleanDetFile();
  
  if( argc > 3 ) df.DumpToFile( argv[3] );
  else df.DumpToFile( );
  
  return 0;
}
