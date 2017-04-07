/*
   $Source: /afs/cern.ch/user/g/gobbo/w0/coral/src/calorim/Reco/src/OneParticleResponse.h,v $
   $Date: 2010/03/31 14:50:08 $
   $Revision: 1.1 $
   -------------------------------------------------------------------------

   This file is part of cellular calorimeter reconstruction program.

   Authors:
     Vladimir  Kolosov   ( Vladimir.Kolosov@cern.ch, Kolosov@mx.ihep.su )
     Alexander Zvyagin   ( Alexander.Zviagine@cern.ch, Zvyagin@mx.ihep.su )
     Denis     Murashev  ( Denis.Mourachev@cern.ch, Murashev@sirius.ihep.su )

   Copyright(C): 1999-2001  V.Kolosov, A.Zvyagin, D.Murashev

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

#ifndef OneParticleResponse___include
#define OneParticleResponse___include

#include <vector>

#include "Calorimeter.h"

////////////////////////////////////////////////////////////////////////////////

namespace Reco {

/*! \brief Hypothesis for reconstruction: one particle
 */
class OneParticleResponse {
    // =========================================================================
    // Constructors and destructor
    // =========================================================================

    public:


        /// Default constructor
        OneParticleResponse   (void) {}

        OneParticleResponse   (const CalorimeterParticle &p, const Calorimeter* c,bool for_development=false);

        OneParticleResponse   (const CalorimeterParticle &p, const Calorimeter* c,
                               const Cell &cell, const std::vector<std::pair<double, double> >& realAmplitudes);
        ///  Make   OneParticleResponse from cluster data
        OneParticleResponse   (const CalorimeterParticle &p, const Calorimeter* c,
                std::vector<CellDataRaw>& cluster_data);
        /// Copy constructor
        OneParticleResponse   (const OneParticleResponse &h);

        /// Destructor
        virtual            ~OneParticleResponse   (void) {}

    // =========================================================================
    // Operators
    // =========================================================================

    public:

        /// Assignment operator
        virtual  OneParticleResponse &operator =    (const OneParticleResponse &h);

        /// Return Data Member   fParticle
        bool                operator ==             (OneParticleResponse &h);

        /// Print hypothesis info to output stream
        friend std::ostream     &operator <<             (std::ostream &o,OneParticleResponse &h);

    // =========================================================================
    // Methods
    // =========================================================================

    public:

        /// \return  amount of hited cells
        size_t              NCells                      (void)  const {return hited_cells.size();}

        /// \return Data Member   fParticle
        const CalorimeterParticle     &GetParticle                 (void)  const {return fParticle;}

        /// \return Data Member   fParticle
        CalorimeterParticle     &GetParticleToModify                 (void)  {return fParticle;}

        /// \return vector of cell numbers (list hitted_cells ).
        const std::vector<size_t> &GetCells                  (void)  const {return hited_cells;}

        /// \return vector of energy deposit in hitted_cells .
        const std::vector<double >   &GetCellsEnergy   (void) const {return e_real_;}
        const std::vector<double >   &GetCellsDispEnergy   (void) const {return de_real_;}

        /// \return vector of expected energy deposit in hitted_cells .
        const std::vector< double >   &GetExpectedCellsEnergy   (void) const {return e_expected_;}
        const std::vector< double >   &GetExpectedCellsDispEnergy   (void) const {return de_expected_;}

        /// \return mask with classification of hitted_cells .
        const std::vector<std::pair<int,int> >& GetMask   (void) const {return fMaskCells;}

        /// Calculate Expected One Particle Response in the Cells
        void                CalculateExpectedResponse   (void);
        void                CalculateExpectedResponse   (const std::vector<std::pair<double, double> >& realAmplitudes);

        /// Calculate Expected One Particle Response in cells used for development FMC
        void                CalculateExpectedResponseForFMC   (void);

        /// \return Sum of expected energy in the Cells
        double              GetEnergyDeposit            (void)  const {return fEnergy_deposit;}

        /// Calculate Derivatives Expected One Particle Response in the Cells. Used in Fit() method.
        void                CalculateDerivativesExpectedResponse(void);
        /// \return  Derivatives in the Cell.
        double              GetDDE( int cell);
        double              GetDDX( int cell);
        double              GetDDY( int cell);
        /// Update Particle Parametes.
        void                UpdateParticleParameters( double e, double se, double x, double sx, double y, double sy);

        /*! Find List of Cells (hited_cells[] vector) were One Particle Response from Calorimeter::Particle &p is Expected.
          Here we use cell search : int p_main_cell =  calorimeter->FindCell(p.GetX(),p.GetY(),p.GetZ());
          which can be rather expensive. Foreseen for some external usage.??
          */
        virtual void        FindList                    (void);

        /*! The best to use in gamma search when the main cell is known.
          For this moment already not all FindList methods return only cell's neighbours
          */
        void                FindList                    (const Cell& cell);

        ///  ExtendList method add neighbors of the hitted_cells to this list

        void                ExtendList                  (void);

        ///  ExtendList method add neighbors located not far distmax from the particle's trace
        void                ExtendList                  (double distmax);

        /// Fast Monte-Carlo method. Add shower fluctuations to One Particle Response in the Cells.
        void                MonteRandomise              (void);

        /*! Share Cell's energy between several Particles
          This method is used if we don't use simultanious fit for all(at least in the same cluster) Particles.
          Return vector of cells wich need recovery procedure
          */
        bool       CalculateRealResponse       ( std::vector < size_t > &problems,
                                                const std::vector<std::pair<double, double> >& realAmplitudes,
                                                const std::vector<std::pair<double, double> >& expectedAmplitudes);

        // Tuning reconstruction, fixing imperfections, creating new
        void                ChangeMain                  (void);

        /// Development method to CalculateParticleResponse more and more ...
        double              TraceInCell                 (const Cell &pcell,double x_in[3],double x_out[3]);

        /// \return Xi^2 of the fit.
        void               CalculateChi2                 (void);
        /// \return Xi^2 of the fit.
        double             GetChi2 ( void ) const { return  fXiSqr;}

        ///  \returnNumber degrees of freedom.
        size_t              GetNDF ( void ) const { return fNDF;}

        std::pair< double, size_t> CheckSinglePeakHypothesis (void) const;
        ///  debug printout
        void   CheckPrint( void ) const;

    protected:

        /*! Calculate expected energy and fluctuations in the cell from the Particle.
          Fluctuations here means intrinsic shower dispersion of cell's energy.
          cell_indx - index in the vector of cells
          */
        void                CalculateParticleResponseInCell  (int cell_indx,
                const Cell* &pcell, const Cell* &p_main_cell,
                const std::vector<std::pair<double, double> >& realAmplitudes,
                bool with_derivatives=false);

        /// Used for fine tuning Z position of the particle in the calorimeter
        void                PropagateParticleAtCorrectZ (void);
        /// Development method to work with particle Trace
        void                SortTrace                   (void);
        void                PrintTrace                  (void);
        void                ClearList                   (void);
        void                InitList                    (const Cell &cell);

    // =========================================================================
    // Data Members
    // =========================================================================

    private:

        /// This is our Calorimeter
        const Calorimeter*                          calorimeter;

        /// The Particle which parameters we try to fit or just use in Monte-Carlo case.
        CalorimeterParticle                         fParticle;

        /// Total energy loss.
        double                                      fEnergy_deposit;

        /// Xi^2 of the fit.
        double                                      fXiSqr;

        /// Number degrees of freedom.
        size_t                                      fNDF;

        double                                      e_overlap_;

        /*! List of cells on which the hypothesis (about sort, energy and position of our Particle) can influence.
          This array is determined by FindList() method.
          */
        std::vector<size_t>                         hited_cells;

        /*! Energy deposit in the hited_cells which we assign to our Particle after some
         *         concurrence between several Particles.
         **/
        std::vector<double>                         e_real_;
        std::vector<double>                         de_real_;

        /// Energy deposit in the hited_cells expected from the Particle.
        std::vector<double>                         e_expected_;
        std::vector<double>                         de_expected_;

        /// Not used for this moment. Forseen to use in Xi^2 calculations?
        std::vector<std::pair<double,double> >      shower_fluctuations_one_particle_response;

        /// This mask is used to speed left,right,up,down neighbours search. Is used in SimpleFit() method.
        std::vector<std::pair<int,int> >            fMaskCells;

        /*! Derivatives of expected_one_particle_response_in_cells[]
          Calculated by CalculateDerivativesExpectedResponse() method. Used in Fit() method.
          Derivative d(one_particle_response_in_cells[])/de
          */
        std::vector<double>                         de_one_particle_response_in_cells;

        /// Derivative d(one_particle_response_in_cells[])/dx
        std::vector<double>                         dx_one_particle_response_in_cells;

        /// Derivative d(one_particle_response_in_cells[])/dy
        std::vector<double>                         dy_one_particle_response_in_cells;

        /*! Some information about particle's trace. Will be combined later in Trace class.
          Speaking about trace in cells: the info( cell#, input point, output point) is considered.
          We have to add this information to simplify calculation of Shower response in cells.
          */
        std::vector<size_t>                         cells_trace;
        std::vector<double>                         v_input_trace_x;
        std::vector<double>                         v_input_trace_y;
        std::vector<double>                         v_input_trace_z;
        std::vector<double>                         v_output_trace_x;
        std::vector<double>                         v_output_trace_y;
        std::vector<double>                         v_output_trace_z;
};

////////////////////////////////////////////////////////////////////////////////

} // using namespace Reco

#endif // OneParticleResponse___include
