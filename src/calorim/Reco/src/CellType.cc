/*
   $Source: /afs/cern.ch/user/g/gobbo/w0/coral/src/calorim/Reco/src/CellType.cc,v $
   $Date: 2010/11/29 19:54:49 $
   $Revision: 1.21 $
   -------------------------------------------------------------------------

   This file is part of cellular calorimeter reconstruction program.

   Authors:
     Vladimir  Kolosov   ( Kolosov@mx.ihep.su )
     Alexander Zvyagin   ( Alexander.Zviagine@cern.ch, Zvyagin@mx.ihep.su )
     Denis     Murashev  ( Denis.Mourachev@cern.ch, Murashev@sirius.ihep.su )

   Copyright(C): 1999-2000  V.Kolosov,A.Zvyagin

     This library is free software; you can redistribute it and/or
     modify it under the terms of the GNU Library General Public
     License as published by the Free Software Foundation; either
     version 2 of the License, or (at your option) any later version.

     This library is distributed in the hope that it will be useful,
     but WITHOUT ANY WARRANTY; without even the implied warranty of
     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
     Library General Public License for more details.

     You should have received a copy of the GNU Library General Public
     License along with this library; if not, write to the Free
     Software Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
*/

#include <algorithm>  // for transform()
#include <cctype>     // for tolower()
#include <iostream>
#include <string>

#include "CellType.h"
#include "ShowerProfile.h"

#include "Exception.h"

using namespace std;

namespace Reco {

////////////////////////////////////////////////////////////////////////////////

CellType::CellType(const string &name,
                   double sizeX, double sizeY, double sizeZ,
                   double stepX, double stepY,
                   double RadiationLength, double NuclearLength,
                   double CriticalEnergy,
                   double StochasticTerm, double ConstantTerm, double ReadoutTerm,
                   double active_material, double density,
                   const std::string &hwtype ) :
  fName                 (name),
  fSizeX                (sizeX),
  fSizeY                (sizeY),
  fSizeZ                (sizeZ),
  fStepX                (stepX),
  fStepY                (stepY),
  fRadiationLength      (RadiationLength),
  fNuclearLength        (NuclearLength),
  fCriticalEnergy       (CriticalEnergy),
  fConstantTerm         (ConstantTerm),
  fStochasticTerm       (StochasticTerm),
  fReadoutTerm          (ReadoutTerm),
  fActiveMaterial       (active_material),
  fDensity              (density),
  fShowerProfile        (0)
{
  fHwType = GetHwTypeFromString(hwtype);

  if( fSizeX<=0. || fSizeY<=0. || fSizeZ<=0. )
    throw Exception("Reco::CellType::CellType():  SizeXYZ must be >0.  SizeX=%g SizeY=%g SizeZ=%g",
                      fSizeX, fSizeY, fSizeZ);

  if( fRadiationLength <= 0. )
    throw Exception("Reco::CellType::CellType():  Radiation length must be >0. fRadiationLength=%g",
                    fRadiationLength);

  if( fNuclearLength <= 0. )
    throw Exception("Reco::CellType::CellType():  Nuclear interaction length must be >0. fNuclearLength=%g",
                    fNuclearLength);

  if( fCriticalEnergy <= 0. )
    throw Exception("Reco::CellType::CellType():  Critical energy length must be >0. fCriticalEnergy=%g",
                    fCriticalEnergy);

  if( fConstantTerm <= 0.   || fConstantTerm > 1. )
    throw Exception("Reco::CellType::CellType():  Constant term must be >0 and <=1 but not %g",
                    fConstantTerm);

  if( fStochasticTerm <= 0. || fStochasticTerm > 1. )
    throw Exception("Reco::CellType::CellType():  Sigma at 1 GeV must be >0 and <=1 but not %g",
                    fStochasticTerm);

  if( fReadoutTerm <= 0.    || fReadoutTerm > 1. )
    throw Exception("Reco::CellType::CellType():  Readout noise must be >0 and <=1 but not %g",
                    fReadoutTerm);

  if( fActiveMaterial <= 0. || fActiveMaterial > 1. )
    throw Exception("Reco::CellType::CellType():  active material must be >0 and <=1 but not %g",
                    fActiveMaterial);
}

////////////////////////////////////////////////////////////////////////////////

void CellType::Print(ostream &o,const string &prefix) const
{
  o << prefix.c_str()<<fName.c_str()<<"  SizeXYZ=(" << fSizeX <<","<<fSizeY<<","<<fSizeZ<<
         ")  RadLen=" << fRadiationLength <<"  NuclearLen=" << fNuclearLength <<
         " ConstantTerm=" << fConstantTerm <<" StochasticTerm="<<fStochasticTerm<<
         " ReadoutTerm="<<fReadoutTerm<<" ActiveMat="<<fActiveMaterial<<" Density="<<fDensity<<endl;

}

////////////////////////////////////////////////////////////////////////////////

CellType::HwType CellType::GetHwTypeFromString (const std::string& hwtype) {
  string hwt(hwtype);
  transform(hwt.begin(), hwt.end(), hwt.begin(), ::tolower);
  if ( hwt == "olga" )
    return hwt_olga;
  else if ( hwt == "mainz" )
    return hwt_mainz;
  else if ( hwt == "gams" )
    return hwt_gams;
  else if ( hwt == "rhgams" )
    return hwt_rhgams;
  else if ( hwt == "shashlik" )
    return hwt_shashlik;
  else if ( hwt == "hcal1" )
    return hwt_hcal1;
  else if ( hwt == "hcal2" )
    return hwt_hcal2;
  else if ( hwt == "unknown" || hwt == "" )
    return hwt_unknown;
  else throw Exception( ("Bad HwType: "+hwt).c_str() );
}

////////////////////////////////////////////////////////////////////////////////

std::string CellType::GetHwTypeString (HwType hwtype) {
  if ( hwtype == hwt_olga )
    return "olga";
  if ( hwtype == hwt_mainz )
    return "mainz";
  if ( hwtype == hwt_gams )
    return "gams";
  if ( hwtype == hwt_rhgams )
    return "rhgams";
  if ( hwtype == hwt_shashlik )
    return "shashlik";
  if ( hwtype == hwt_hcal1 )
    return "hcal1";
  if ( hwtype == hwt_hcal2 )
    return "hcal2";

  return "unknown";
}

////////////////////////////////////////////////////////////////////////////////

} // namespace Reco
