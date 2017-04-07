/*!
   \file    getPatterns.cc
   \brief   Compass Vertex Pattern Class.
   \author  Alexandre Korzenev
   \version $Revision: 1.2 $
   \date    $Date: 2003/04/24 08:24:20 $ 

*/
#include "CsRolandPattern.h"
#include "CsEvent.h"
#include "CsGeant3.h"


bool CsRolandPattern::getPatterns( std::list<CsVertex*>& vrts, std::map<CsTrack*,bool>& specials )
{
  if( vrts_.size() ) {
    vrts = vrts_;
    return true;
  } else
    return false;
}
