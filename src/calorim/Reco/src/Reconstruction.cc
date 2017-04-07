/*
   $Source: /afs/cern.ch/user/g/gobbo/w0/coral/src/calorim/Reco/src/Reconstruction.cc,v $
   $Date: 2010/03/31 14:50:08 $
   $Revision: 1.60 $
*/

// --- Standard C/C++ library ---
#include <cassert>
#include <cmath>
#include <iostream>

// --- Internal files ---
#include "Reconstruction.h"

namespace Reco {

////////////////////////////////////////////////////////////////////////////////

Reconstruction::Reconstruction(const Calorimeter* c) : calorimeter_(c) {
}

////////////////////////////////////////////////////////////////////////////////

} // namespace Reco
