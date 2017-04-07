/*
   $Source: /afs/cern.ch/user/g/gobbo/w0/coral/src/calorim/Reco/src/ShowerProfile.cc,v $
   $Date: 2010/11/29 19:16:07 $
   $Revision: 1.2 $
*/

#include "ShowerProfile.h"

#include "Calorimeter.h"

namespace Reco {

////////////////////////////////////////////////////////////////////////////////

ShowerProfile::ShowerProfile(const Calorimeter* c, const CellType& ct) :
    calorimeter_(c), cellType_(ct) {
}

////////////////////////////////////////////////////////////////////////////////

ShowerProfile::~ShowerProfile() {
}

////////////////////////////////////////////////////////////////////////////////

} // namespace Reco

