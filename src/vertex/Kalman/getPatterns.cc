/*!
   \file    getPatterns.cc
   \brief   Compass Vertex Pattern Class.
   \author  Alexandre Korzenev
   \version $Revision: 1.2 $
   \date    $Date: 2003/04/24 08:21:31 $ 

*/

#include "CsAverPattern.h"

bool CsAverPattern::getPatterns( std::list<CsVertex*>& vrts, std::map<CsTrack*,bool>& specials )
{
  
  if( vrts_.size() ) {
    vrts = vrts_;
    specials = specials_;
    return true;
  } else
    return false;

}
