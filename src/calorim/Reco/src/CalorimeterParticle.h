/*
   $Source: /afs/cern.ch/user/g/gobbo/w0/coral/src/calorim/Reco/src/CalorimeterParticle.h,v $
   $Date: 2010/06/29 15:47:03 $
   $Revision: 1.24 $
   -------------------------------------------------------------------------

   This file is part of cellular calorimeter reconstruction program.

   Authors:
     Vladimir  Kolosov   ( Vladimir.Kolosov@cern.ch, Kolosov@mx.ihep.su )

   Copyright(C): 1999-2002  V.Kolosov, A.Zvyagin, D.Murashev

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

#ifndef CalorimeterParticle___include
#define CalorimeterParticle___include

#include <iostream>
#include <string>
#include <vector>
#include <map>

#include "Reco_config.h"

#include <TLorentzVector.h>

#include "Exception.h"

////////////////////////////////////////////////////////////////////////////////

namespace Reco {

  class Calorimeter;

 /*! \brief This is internal and output interface class.

     Particle properties from a calorimeter's point of view.
     This is not general particle description. A calorimeter can not detect
     directly particle's spin, charge, lifetime... The only things it can detect
     are particle coordinates on calorimeter plane (X,Y,Z) and  sometimes
     particle energy (E) and time.
 */
  class CalorimeterParticle
  {
    //============================================================================
    // Types, constants
    //============================================================================

    public:

      /// Particle identification. The same numbers as in GEANT 3.
      enum ParticleID
      {
        UNKNOWN       =  0,
        GAMMA         =  1,
        POSITRON      =  2,
        ELECTRON      =  3,
        NEUTRINO      =  4,
        MUON_PLUS     =  5,
        MUON_MINUS    =  6,
        PION_0        =  7,
        PION_PLUS     =  8,
        PION_MINUS    =  9,
        KAON_0_LONG   = 10,
        KAON_PLUS     = 11,
        KAON_MINUS    = 12
     };
      /// Misc data identification.
      enum MiscDataID
      {
        CHISQ_SADC_SIGNAL,
        NDF_SADC_SIGNAL,
        CHISQ_SADC_NOISE,
        NDF_SADC_NOISE,
        CHISQ_CLUSTER_DATA,
        NDF_CLUSTER_DATA,
        RW_CHRG_PROB,
        SHAPE_SADC_OK,
        VALUE_R4,
        VALUE_R9,
        CLUSTER_SIZE,
        CLUSTER_SHOWERS,        // number of showers in the same cluster as this
                                // shower (including the current shower)
        CHI_GAM,
        NDF_GAM,

        UNCORR_X,               // uncorrected positions and energy
        UNCORR_Y,
        UNCORR_E
      };

    //============================================================================
    // Constructors and destructor
    //============================================================================

    public:

      /// Destructor
      ~CalorimeterParticle            (void) {}

      /*! \brief Particle constructor

          \param id identification
          \param prob
          \param e energy
          \param x X coordinate
          \param y Y coordinate
          \param z Z coordinate
          \param calo pointer to calorimeter the particle belongs to
          \param e_err energy error
          \param x_err X coordinate error
          \param y_err Y coordinate error
          \param z_err Z coordinate error
          \param angle_x angle in ZY plane
          \param angle_y angle in ZX plane
      */
      CalorimeterParticle ( ParticleID id_=UNKNOWN, double prob_=0., double e_=0., double x=0., double y=0., double z=0.,
                            const Calorimeter *calo=NULL,
			    double e_err_=0., double x_err=0., double y_err=0., double z_err=0.,
			    double angle_x=0., double angle_y=0. ) :
        id( id_ ), prob( prob_ ), e( e_ ), e_err( e_err_ ),
        has_time( false ), time( -100000000. ), time_err( 100000000. ),
        cluster_size_in_cells( 0 ),
        calorimeter_name( "" ), calorimeter( calo )
	{
          xyz[0]      = x;
          xyz[1]      = y;
          xyz[2]      = z;
          xyz_err[0]  = x_err;
          xyz_err[1]  = y_err;
          xyz_err[2]  = z_err;
          xyz_gate[0] = -1.;
          xyz_gate[1] = -1.;
          xyz_gate[2] = -1.;
          angle[0]    = angle_x;
          angle[1]    = angle_y;
          if ( e <= 0 ) throw Exception( "CalorimeterParticle has negative energy e=%f", e);
	}

    //============================================================================
    // Operators
    //============================================================================

    public:

      /// Print particle properties
      friend std::ostream &operator <<             (std::ostream &o,const CalorimeterParticle &particle) {particle.Print(o);return o;}

    //============================================================================
    // Methods
    //============================================================================

      /// Print particle properties
      void            Print                   (std::ostream &o=std::cout) const;

      /// \return Particle identification (GEANT 3 conversions)
      ParticleID      GetID                   (void) const {return id;}

      /// \return Particle energy
      double          GetE                    (void) const {return e;}

      /// \return Particle X coordinate
      double          GetX                    (void) const {return xyz[0];}

      /// \return Particle Y coordinate
      double          GetY                    (void) const {return xyz[1];}

      /// \return Particle Z coordinate
      double          GetZ                    (void) const {return xyz[2];}

      /// \return Particle energy error
      double          GetEerr                 (void) const {return e_err;}

      /// \return Particle X coordinate error
      double          GetXerr                 (void) const {return xyz_err[0];}

      /// \return Particle Y coordinate error
      double          GetYerr                 (void) const {return xyz_err[1];}

      /// \return Particle Z coordinate error
      double          GetZerr                 (void) const {return xyz_err[2];}

      /// \return Set cell size related X gate
      void             SetGateX                 ( double v ) { xyz_gate[0]= v;}

      /// \return Set cell size related Y gate
      void             SetGateY                 ( double v ) { xyz_gate[1]= v;}

      /// \return Set cell size related Z gate
      void             SetGateZ                 ( double v ) { xyz_gate[2]= v;}

      /// \return cell size related X gate
      double          GetGateX                 (void) const {
            if( xyz_gate[0] < 0.) std::cerr << " WARNING! Usage of not initialized CalorimeterParticle::GetGateX " << std::endl;
                                                              return xyz_gate[0];}

      /// \return cell size related Y gate
      double          GetGateY                 (void) const {
            if( xyz_gate[1] < 0. ) std::cerr << " WARNING! Usage of not initialized CalorimeterParticle::GetGateY " << std::endl;
                                                              return xyz_gate[1];}

      /// \return cell size related Z gate
      double          GetGateZ                 (void) const {
            if( xyz_gate[2] < 0. ) std::cerr << " WARNING! Usage of not initialized CalorimeterParticle::GetGateZ " << std::endl;
                                                              return xyz_gate[2];}

      /// \return Angle in ZY plane
      double          GetAngleX               (void) const {return angle[0];}

      /// \return Angle in ZX plane
      double          GetAngleY               (void) const {return angle[1];}

      /// \return Angle error in ZY plane
      double          GetAngleXerr            (void) const {return angle_err[0];}

      /// \return Angle error in ZY plane
      double          GetAngleYerr            (void) const {return angle_err[1];}

      /// \return Probability
      double          GetProb                 (void) const {return prob;}

      /// Get main cells (for simple calorimeters there is only one main cell,
      /// but in general there might be more, like eg. for a calorimeter that
      /// is segmented in z)
      const std::vector<size_t> &GetMainCells     (void) const {return main_hited_cells;}

      /// Deprecated function for backwards compatibility with PHAST, use
      /// GetMainCells() instead!
      /// \todo remove after 2011-06-30
      const std::vector<size_t> &GetHitedCells    (void) const {return main_hited_cells;}

      /// \return energy in hitted(main) cells. If  vector is empty this information is not provided.
      const std::vector<double> &GetEnergyInMainCells     (void) const {return energy_in_main_cells;}

      /// \return complete information about calorimeter cluster
      const std::vector < std::pair< size_t,double> >  &GetClusterData     (void) const {return cluster_data;}

      /// Set complete information about calorimeter cluster
      void            SetClusterData     (const std::vector < std::pair< size_t,double> > &cluster);
      // The same as SetClusterData but no warnings in case of data replacement
      void            ReplaceClusterData     (const std::vector < std::pair< size_t,double> > &cluster);

      /// \return complete information about calorimeter cluster
      size_t          GetClusterSize     (void) const {return cluster_size_in_cells;}

      /// Set size of cluster with >90% of total energy
      void            SetClusterSize        (size_t ncells);

      /// Set particle identification
      void            SetID                   (ParticleID v)            {id=v;}

      /// Set particle energy and energy error
      void            SetE                    (double v,double v_err=0)
                                                {
                                                  if( v<=0 || v_err<0 )
                                                    throw Exception("Calorimeter::Particle::SetE: negative energy e=%g e_error=%g",v,v_err);
                                                  e=v;
                                                  e_err=v_err;
                                                }

      /// Set particle X coordinate and coordinate error
      void            SetX                    (double v,double v_err=0) {xyz[0]=v; xyz_err[0]=v_err;}

      /// Set particle Y coordinate and coordinate error
      void            SetY                    (double v,double v_err=0) {xyz[1]=v; xyz_err[1]=v_err;}

      /// Set particle Z coordinate and coordinate error
      void            SetZ                    (double v,double v_err=0) {xyz[2]=v; xyz_err[2]=v_err;}

      /// Set particle angle in ZY plane (with error)
      void            SetAngleX               (double v,double v_err=0) {angle[0]=v; angle_err[0]=v_err;}

      /// Set particle angle in ZX plane (with error)
      void            SetAngleY               (double v,double v_err=0) {angle[1]=v; angle_err[1]=v_err;}

      /// Set probability.
      void            SetProb                 (double v)                {prob=v;}

      ///  Set hitted cell number to the list of main hited_cells
      void            SetHitedCell            (size_t m) {main_hited_cells.push_back(m);}

      /// Set main cells (for simple calorimeters there is only one main cell,
      /// but in general there might be more, like eg. for a calorimeter that
      /// is segmented in z)
      void            SetMainCells            (const std::vector<size_t> &main_cells) {main_hited_cells=main_cells;}

      /// \return validity of timing information
      bool            HasTime                 (void) const {return has_time;}

      /// \return time
      double          GetTime                 (void) const {return time;}

      /// \return time error
      double          GetTimeErr              (void) const {return time_err;}

      ///  Set timing information
      void            SetTime                 (double v,double v_err=0) {has_time=true,time=v,time_err=v_err;}

      /// \return Calorimeter name.
      const std::string   &GetCalorimeterName (void) const;

      /// \return Pointer to calorimeter.
      const Calorimeter   *GetCalorimeter     (void) const { return calorimeter; }

      /// Set fictional Calorimeter name, for RW/ECAL1 reconstruction "RW_ECAL1" is used.
      void            SetFictionalCalorimeterName(const std::string &name) {calorimeter_name=name;}

      /// Set misc info.
      void            SetMiscInfo         (MiscDataID data_id, double data);

      /// Get misc info.
      std::pair <bool, double >           GetMiscInfo         (MiscDataID data_id) const;

      /// Get misc data.
      std::map< MiscDataID, double > &    GetMiscData(void) {return misc_data_;}

      TLorentzVector         GetMomentum ( void ) const;
      /// Set mass
      void                   SetMass     ( double mass ) { mass_ = mass; }
      /// Get mass
      double                 GetMass     ( void ) const { return mass_; }

      /// Calculate Z position under the assumption, that the particle belongs
      /// to the charged track with the given parameters.
      /// \param mip track is from MIP (i.e. pion/muon for ECAL or muon for HCAL)
      double                 CalcZforTrack( double dxdz, double dydz, bool mip );

    //============================================================================
    // Attributes, data
    //============================================================================

      private:

      /// Particle identification
      ParticleID      id;

      /// Probability that particle has these properties.
      float           prob;

      /// Particle mass
      double          mass_;

      /// Particle energy
      double          e;

      /// Particle energy error
      double          e_err;

      /// Particle coordinates
      double          xyz[3];

      /// Particle coordinates error
      double          xyz_err[3];

      /// Relaxed gate for coordinates related to cell size
      double          xyz_gate[3];

      /// Particle track projection angle X,Y.
      double          angle[2];

      /// Particle track projection angle error X,Y.
      double          angle_err[2];

      /*! \brief  List of MAIN Calorimeter cells hited by the Particle.
                  If the Particle missed list is empty.
          IMPORTANT WARNING!!! THIS IS NOT LIST OF HITED CELLS !!!
          ********************
          NORMALY THE SIZE OF   main_hited_cells.size() = 1 !!!!
          ONLY FOR SPECIAL (SEVERAL LAYERS OR(AND) PROJECTION TYPE CALORIMETERS) IT COULD BE > 1 !!!!
                   *******
      */
      std::vector<size_t>  main_hited_cells;

      /// Energy in main cells. This information may not be provided.
      std::vector<double>  energy_in_main_cells;

      /// Particle might have some timing information in Calorimeter
      bool            has_time;

      /// time
      double          time;

      /// time error
      double          time_err;

      /*! \brief  Detailed information about calorimeter object.
          Each vector element corresponds to the cell ( first element is a cell number in Calorimeter's Cells list)
          and the energy in the cell (second element).
          If  cluster_data is empty this information is not provided.
      */
      std::vector < std::pair< size_t,double> >     cluster_data;

      /// This is not the same as cluster_data.size(), this is amount of cells with >90% of total energy deposit.
      size_t                              cluster_size_in_cells;

      /// The calorimeter's name which produce this particle. This field for sure may be empty.
      std::string                              calorimeter_name;

      /// Pointer to calorimeter this particle belongs to.
      const Calorimeter                       *calorimeter;

      ///  Misc data
      std::map< MiscDataID, double >          misc_data_;

  };

} // using namespace Reco

#endif // CalorimeterParticle___include
