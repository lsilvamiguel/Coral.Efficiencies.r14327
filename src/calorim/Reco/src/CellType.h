/*
   $Source: /afs/cern.ch/user/g/gobbo/w0/coral/src/calorim/Reco/src/CellType.h,v $
   $Date: 2010/11/29 19:54:49 $
   $Revision: 1.10 $
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

#ifndef RecoCellType___include
#define RecoCellType___include

#include <string>
#include <vector>
#include <iostream>

namespace Reco {

////////////////////////////////////////////////////////////////////////////////

// forward declaration
class ShowerProfile;

/*! \brief Properties of calorimeter cell.

    This is the base cell properties: size,density,radiation length,...

    Do not put CellType in std::vector<CellType>. Use std::list<CellType> or std::slist<CellType>
    instead of this.
*/
class CellType
{
  //============================================================================
  // Constructors and destructor
  //============================================================================

  public:

  enum HwType { hwt_olga, hwt_mainz, hwt_gams, hwt_rhgams, hwt_shashlik,
                hwt_hcal1, hwt_hcal2, hwt_unknown };

    /// Destructor
    virtual            ~CellType                (void) {}

    /*! \brief Construct CellType with given properties.

        \param name is the name of this cell type (do not use spaces)
        \param sizeX
        \param sizeY size of cell (cell shape is box)
        \param sizeZ
        \param RadiationLength radiation length
        \param CriticalEnergy energy at which electron ionisation losses are equal to radiative losses
        \param NuclearLength nuclear interaction length
        \param StochasticTerm dispertion of photoelectrons

        All values are in GeV:
        \fn
          \delta E = \sigma_{Readout} + \sigma_{Stochastic} * \sqrt{E} + \sigma_{Constant}*E
        \endfn
    */
                            CellType            ( const std::string &name,
                                                  double sizeX, double sizeY, double sizeZ,
                                                  double stepX, double stepY,
                                                  double RadiationLength, double NuclearLength,
                                                  double CriticalEnergy,
                                                  double StochasticTerm, double ConstantTerm, double ReadoutTerm,
                                                  double active_material, double density,
                                                  const std::string &hwtype );
                            CellType() {};

  //============================================================================
  // Operators
  //============================================================================

  public:

    /// Print CellType info to output stream
    friend std::ostream&    operator <<      (std::ostream &o,const CellType &cell_type)
                                                            {cell_type.Print(o);return o;}

  //============================================================================
  // Static methods
  //============================================================================

  public:

    /// \return cell type extracted from input string
    static HwType           GetHwTypeFromString     (const std::string& hwtype);

    /// \return cell type name
    static std::string      GetHwTypeString         (HwType hwtype);

  //============================================================================
  // Methods
  //============================================================================

  public:

    /// \return Name of this cell type.
    const std::string       GetName                 (void) const {return fName;}

    /// \return cell X step in mm
    double                  GetSizeX                (void) const {return fStepX;}

    /// \return cell Y step in mm
    double                  GetSizeY                (void) const {return fStepY;}

    /// \return cell Z size in mm
    double                  GetSizeZ                (void) const {return fSizeZ;}

    /// \return cell X size in mm
    double                  GetTrueSizeX            (void) const {return fSizeX;}

    /// \return cell Y size in mm
    double                  GetTrueSizeY            (void) const {return fSizeY;}

    /// \return cell Z size in mm
    double                  GetTrueSizeZ            (void) const {return fSizeZ;}

    /// X-position of cell center in mm.
    double                  GetStepX                (void) const {return fStepX;}

    /// Y-position of cell center in mm.
    double                  GetStepY                (void) const {return fStepY;}

    /// \return cell's material radiation length in mm
    double                  GetRadiationLength      (void) const {return fRadiationLength;}

    /// \return cell's material radiation length in mm
    void                    SetRadiationLength      (double v) {fRadiationLength=v;}

    /// \return cell's material nuclear length in mm
    double                  GetNuclearLength        (void) const {return fNuclearLength;}

    /// \return cell's material nuclear length in mm
    void                    SetNuclearLength        (double v)  {fNuclearLength=v;}

    /// \return cell's material critical energy in GeV
    double                  GetCriticalEnergy       (void) const {return fCriticalEnergy;}

    /// \return SigmaE/E constant term.
    double                  GetConstantTerm         (void) const {return fConstantTerm;}

    /// \return SigmaE/E at 1 GeV. Statistic of photo-electrons.
    double                  GetStochasticTerm       (void) const {return fStochasticTerm;}

    /// \return SigmaE Readout term.
    double                  GetReadoutTerm          (void) const {return fReadoutTerm;}
    void                    SetReadoutTerm          (double v) {fReadoutTerm=v;}

    /// \return ratio of active material
    double                  GetActiveMaterial       (void) const {return fActiveMaterial;}

    /// Set ratio of active material.
    void                    SetActiveMaterial       (double v) {fActiveMaterial=v;}

    /// \return hardware type of the calorimeter cell (eg. "gams" or "mainz")
    HwType                  GetHwType               (void) const {return fHwType;}

    /// \return pointer to shower profile
    const ShowerProfile*    GetShowerProfile        (void) const { return fShowerProfile; }

    /// set the shower profile for this cell type
    void                    SetShowerProfile        (ShowerProfile* showerProfile) { fShowerProfile = showerProfile; }

    /// Print properties
    void                    Print                   (std::ostream &o=std::cout,const std::string &prefix="") const;

  //============================================================================
  // Attributes, data
  //============================================================================

  private:

    /// Name of the cell type
    std::string         fName;

    /// X cell size [mm]
    double              fSizeX;
    /// Y cell size [mm]
    double              fSizeY;
    /// Z cell size [mm]
    double              fSizeZ;

    /// X distance to next cell (center-to-center) [mm]
    double              fStepX;
    /// Y distance to next cell (center-to-center) [mm]
    double              fStepY;

    /// Material radiation length.
    double              fRadiationLength;

    /// Material nuclear length.
    double              fNuclearLength;

    /// Critical energy
    double              fCriticalEnergy;

    /// SigmaE/E  constant term.
    double              fConstantTerm;

    /*! SigmaE/E at 1 GeV. Statistic of photo-electrons.
        Depends from read-out that is why and seems the same for the definite cell type.
    */
    double              fStochasticTerm;

    /// SigmaE read-out noise in GeV
    double              fReadoutTerm;

    /// Ratio of active material. (0,1]
    double              fActiveMaterial;

    /// Average density
    double              fDensity;

    /// hardware type of the calorimeter cell (eg. "gams" or "mainz")
    HwType              fHwType;

    /// pointer to used shower profile
    ShowerProfile*      fShowerProfile;
};

////////////////////////////////////////////////////////////////////////////////

} /// using namespace Reco

#endif //RecoCellType___include
