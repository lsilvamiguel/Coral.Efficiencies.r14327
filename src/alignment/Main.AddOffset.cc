// $Id: Main.AddOffset.cc,v 1.3 2006/06/16 15:21:47 conrad Exp $
/*!
   \file    Main.AddOffset.cc
   \brief   Instanciate a DetFileManager object, parse detectors, adds offset to Gems and Micromegas by hand, needed to go from alignment run to physics run
   \author  Hugo Pereira
   \version $Revision: 1.3 $
   \date    $Date: 2006/06/16 15:21:47 $
*/

#include <iostream>
#include <string>
#include "Macro.h"
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
    cout << "Usage: AddOffset <detector.dat> <new detector.dat>" << endl;
    return 0;
  }
  
  string detFile(argv[1]);
  DetFileManager df( detFile.c_str() );
  
  Macro::AddGemOffset( df );
  Macro::AddMMOffset( df );
  
  if( argc > 2 ) df.DumpToFile( argv[2] );
  else df.DumpToFile( );
  
  return 0;
}
