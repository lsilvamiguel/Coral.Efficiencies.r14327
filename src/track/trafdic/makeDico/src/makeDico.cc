// $Id: makeDico.cc,v 1.1 2000/06/25 18:31:26 ybedfer Exp $

#include "Coral.h"
#include "CsGeom.h"
// tracking stuff
#include "CsTrkUtils.h"
#include "CsTrafficUtils.h"

/*!
  \fn int main( int argc, char *argv[] ) {
  \brief Generate Dico == Lattice file (Lattice is a tracks look-up table
    used for track fitting. Generation uses Traffic track propagation
  \author Y.B
  \version $Revision: 1.1 $
  \date $Date: 2000/06/25 18:31:26 $
*/

int main( int argc, char *argv[] ) {

  // Package Initialization 
  Coral* coral        = Coral::init( argc, argv );

  std::cout << "Comgeant geometry version :"
       << CsGeom::Instance()->getGeomVers() << std::endl;

  // Define tracking procedures (by hand for the moment)...
  CsTrafficUtils DicoUtil;

  CsTrkUtils*    trkutil = &DicoUtil;

  return (trkutil->genLattice());
}


