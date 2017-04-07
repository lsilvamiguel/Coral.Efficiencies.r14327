/*!
   \file    CsRolandPattern.cc
   \brief   Compass Vertex Pattern Class.
   \author  Alexandre Korzenev
   \version $Revision: 1.8 $
   \date    $Date: 2010/01/28 12:51:26 $ 

*/
#include "CsRolandPattern.h"
#include "CsEvent.h"
#include "CsGeant3.h"
#include "CsOpt.h"

using namespace std;

CsRolandPattern::CsRolandPattern() 
{
  CsOpt* opt = CsOpt::Instance();

  int NUMB;

  hist_ = false;
  opt->getOpt( "CsRolandPattern", "Hist", NUMB );
  hist_ = (bool)NUMB;
  
  NSpec_ = 0;
  opt->getOpt( "CsRolandPattern", "NSpec", NSpec_ );

  list<string> Cuts;
  c_[0] = 40;
  c_[1] = 10;
  c_[2] = 5;
  c_[3] = 0.75;
  
  if( opt->getOpt("CsRolandPattern","CUTS",Cuts) ){
    list<string>::iterator Is;
    int i = 0;
    for( Is=Cuts.begin(); Is!=Cuts.end() && i<4; Is++, i++ ) {
      istringstream( (*Is) ) >> c_[i];
    }
  }
}

