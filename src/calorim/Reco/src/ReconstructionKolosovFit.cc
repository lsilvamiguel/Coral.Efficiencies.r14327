/*
   $Source: /afs/cern.ch/user/g/gobbo/w0/coral/src/calorim/Reco/src/ReconstructionKolosovFit.cc,v $
   $Date: 2010/04/15 14:20:02 $
   $Revision: 1.4 $
   -------------------------------------------------------------------------

   This file is part of cellular calorimeter reconstruction program.

   Authors:
     Vladimir  Kolosov   ( Kolosov@mx.ihep.su )
     Alexander Zvyagin   ( Alexander.Zviagine@cern.ch, Zvyagin@mx.ihep.su )

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

#include "ReconstructionKolosov.h"

#include <cmath>

#include "mycernlib.h"
#include "CellType.h"
#include "Shower.h"

namespace Reco {

////////////////////////////////////////////////////////////////////////////////

/*!
  \brief Interface to shower fitting
  Call fitting method dependent from fitting_type request
  \param fitting_type {Select fitting option to use}
  \callgraph
  \callergraph
*/
void ReconstructionKolosov::Fit(const Calorimeter::FitMethod& fitting_type, OneParticleResponse& partResponse) {
    if      ( fitting_type == Reco::Calorimeter::Simple ) SimpleFit(partResponse);
    else if ( fitting_type == Reco::Calorimeter::Normal ) Fit      (partResponse);
    else if ( fitting_type == Reco::Calorimeter::NoFit  ) NoFit    (partResponse);
    else {
        std::cerr << " Fitting method is not initialized " << std::endl;
        exit(1);
    }
}

////////////////////////////////////////////////////////////////////////////////

/*!
  \brief Determine particle position from cluster main cell
  \callgraph
  \callergraph
  \remarks {
  Used options:
  energy_leak_for_SimpleFit
  }
*/
void ReconstructionKolosov::NoFit(OneParticleResponse& partResponse) {

    if( partResponse.GetCells().empty() )
        return;

    // Calculate total energy total energy in 9 cell block
    double Ecell, Etotal=0., Etotal9=0.;
    for (size_t i=0; i<partResponse.GetCells().size(); i++ ) {
        Ecell = partResponse.GetCellsEnergy()[i];
        Etotal += Ecell;

        if( abs(partResponse.GetMask()[i].first) <= 1 && abs(partResponse.GetMask()[i].second) <= 1 ) {
            Etotal9 += Ecell;
        }
    }

    // Check that there is energy in 9 cell block
    if( Etotal9 == 0 )
    {
        std::cerr << "ReconstructionKolosov::NoFit():\n"
                                                            << "  Wrong Gamma in NoFit Etotal9=0 " << std::endl;
        return;
    }

    // Get main cell
    assert(partResponse.GetCells().size()>0 && partResponse.GetCells().front()<GetCalorimeter()->NCells());
    size_t main_cell = partResponse.GetCells().front();
    const Cell*  p_main_cell = &GetCalorimeter()->GetCells()[main_cell];

    // Get Size of main cell
    double SizeX = p_main_cell->GetCellType().GetSizeX()/2.;
    double SizeY = p_main_cell->GetCellType().GetSizeY()/2.;

    // Get position of main cell
    double x0 = p_main_cell->GetX();
    double y0 = p_main_cell->GetY();

    // Get position error estimation (half size)
    double sigmax = SizeX/2.;
    double sigmay = SizeY/2.;

    // Set particle position to main cell position
    partResponse.GetParticleToModify().SetX( x0 , sigmax );
    partResponse.GetParticleToModify().SetY( y0 , sigmay );

    double SizeZ = p_main_cell->GetCellType().GetSizeZ()/2.;
    double z0 = p_main_cell->GetZ()-SizeZ;
    double zw = SizeZ;
    double sigmaz = 0;
    partResponse.GetParticleToModify().SetZ( zw+z0 , sigmaz );

    // Set energy for cluster from total energy
    double fraction = 1.+GetCalorimeter()->GetOptions().energy_leak_for_SimpleFit;
    double eg = fraction*Etotal;
    double se = GetCalorimeter()->EnergySigmaInCell(eg, main_cell);
    partResponse.GetParticleToModify().SetE( eg , se );

}

////////////////////////////////////////////////////////////////////////////////

/*!
  \brief Calculate particle position using deterministik hit combination
  \callgraph
  \callergraph
*/
void ReconstructionKolosov::SimpleFit(OneParticleResponse& partResponse) {
    bool debug = false;

    if( debug ) std::cout << " Start SimpleFit debug E = " << partResponse.GetParticle().GetE() << " PID " << partResponse.GetParticle().GetID() <<
        "Emain=" << partResponse.GetCellsEnergy()[0] << std::endl;

    bool check_function = false;
    bool invert_function = false;

    if( partResponse.GetCells().empty() )
        return;

    // Calculate total energy total energy in 9 cell block and direct neighbors of main cell
    double Ecell, Etotal=0., Etotal9=0., Eleft=0., Eright=0., Edown=0., Eup=0.;
    for (size_t i=0; i<partResponse.GetCells().size(); i++) {

        Ecell = partResponse.GetCellsEnergy()[i];
        Etotal += Ecell;

        if( abs(partResponse.GetMask()[i].first) <= 1 && abs(partResponse.GetMask()[i].second) <= 1 ) {
            Etotal9 += Ecell;
            if( partResponse.GetMask()[i].first  == -1 ) Eleft  += Ecell;
            if( partResponse.GetMask()[i].first  ==  1 ) Eright += Ecell;
            if( partResponse.GetMask()[i].second == -1 ) Edown  += Ecell;
            if( partResponse.GetMask()[i].second ==  1 ) Eup    += Ecell;
        }

    }

    if( debug ) std::cout << "  Etotal9 = " << Etotal9 << std::endl;

    // Check that there is energy in 9 cell block
    if( Etotal9 == 0 )
    {
        std::cerr << "ReconstructionKolosov::SimpleFit():\n"
                                                                << "  Wrong Gamma in SimpleFit Etotal9=0 " << std::endl;
        return;
    }

    // Get main cell
    assert(partResponse.GetCells().size()>0 && partResponse.GetCells().front()<GetCalorimeter()->NCells());
    size_t main_cell = partResponse.GetCells().front();
    const Cell*  p_main_cell = &GetCalorimeter()->GetCells()[main_cell];

    // Get Cell size of main cell
    double SizeX = p_main_cell->GetCellType().GetSizeX()/2.;
    double SizeY = p_main_cell->GetCellType().GetSizeY()/2.;

    // Radiation length of main cell
    double RadLeng = p_main_cell->GetCellType().GetRadiationLength();

    // Nuclear interaction length of main cell
    double NuclLeng = p_main_cell->GetCellType().GetNuclearLength();

    if( debug ) std::cout << "  RadLeng = " << RadLeng << " NuclLeng " << NuclLeng << std::endl;




    // Get particle ID
    CalorimeterParticle::ParticleID pid = partResponse.GetParticle().GetID();

    // Get track angels
    double ax = partResponse.GetParticle().GetAngleX();
    double ay = partResponse.GetParticle().GetAngleY();

    // Get position
    double x0 = p_main_cell->GetX();
    double y0 = p_main_cell->GetY();

    // Set variables for coordinate function
    double dirx = 1.;
    double El = Etotal9-Eleft-Eright;
    double Er = Eright;
    if( Eleft > Eright ) {
        dirx = -1.;
        Er = El;
        El = Eleft;
    }
    double diry = 1.;
    double Ed = Etotal9-Edown-Eup;
    double Eu = Eup;
    if( Edown > Eup    ) {
        diry = -1.;
        Eu = Ed;
        Ed = Edown;
    }

    // Get shifted coordinate (by half cell dimension)
    double xnew = x0 + dirx*SizeX;
    double ynew = y0 + diry*SizeY;

    //
    double xx=0;
    double yy=0;


    if( pid <= 3 || (pid == 5 || pid == 6) )   // gamma,electron,positron, (muon)
    {
        xx = CoordinateFunction(RadLeng,El,Er,ax);
        yy = CoordinateFunction(RadLeng,Ed,Eu,ay);
    }
    else if( pid >= 7  )   // hadron (pion_0 can not appear)
    {
        xx = HadronicCoordinateFunction(NuclLeng,RadLeng,El,Er,ax);
        yy = HadronicCoordinateFunction(NuclLeng,RadLeng,Ed,Eu,ay);
    }


    // Set cell center if too far
    if( xx < -SizeX ) xx = -SizeX;
    if( xx >  SizeX ) xx =  SizeX;
    if( yy < -SizeY ) yy = -SizeY;
    if( yy >  SizeY ) yy =  SizeY;

    double xxb(xx),yyb(yy);

    double x1 = xnew + xx - x0;
    double wx1 = El;
    if( Er < El ) wx1= Er;

    double y1 = ynew + yy - y0;
    double wy1= Ed;
    if(Eu < Ed ) wy1= Eu;


    // Repeate with different conditions
    dirx = 1.;
    diry = 1.;

    El = Etotal9-Eleft-Eright;
    Er = Eright;
    Ed = Etotal9-Edown-Eup;
    Eu = Eup;

    if( Eleft < Eright ) {
        dirx = -1.;
        Er = El;
        El = Eleft;
    }

    if( Edown < Eup    ) {
        diry = -1.;
        Eu = Ed;
        Ed = Edown;
    }

    xnew = x0 + dirx*SizeX;
    ynew = y0 + diry*SizeY;

    if( pid <= 3 ||  (pid == 5 || pid == 6) )   // gamma,electron,positron,MUON
    {
        xx = CoordinateFunction(RadLeng,El,Er,ax);
        yy = CoordinateFunction(RadLeng,Ed,Eu,ay);
    }
    else if( pid >= 7  )   // hadron (pion_0 can not appear)
    {
        xx = HadronicCoordinateFunction(NuclLeng,RadLeng,El,Er,ax);
        yy = HadronicCoordinateFunction(NuclLeng,RadLeng,Ed,Eu,ay);
    }

    // Set cell center if too far
    if( xx < -2.*SizeX ) xx = -2.*SizeX;
    if( xx >  2.*SizeX ) xx =  2.*SizeX;
    if( yy < -2.*SizeY ) yy = -2.*SizeY;
    if( yy >  2.*SizeY ) yy =  2.*SizeY;

    double x2=xnew+xx-x0;
    double wx2= El;
    if(Er < El ) wx2= Er;

    double y2=ynew+yy-y0;
    double wy2= Ed;
    if(Eu < Ed ) wy2= Eu;


    // Setting xw yw
    double xw = 0;
    double yw = 0;
    if( check_function )
    {
        if( wx1 > wx2 )
            if(!invert_function )
                xw = x1;
            else
                xw = x2;
        else
            if(!invert_function )
                xw = x2;
            else
                xw = x1;

        if( wy1 > wy2 )
            if(!invert_function )
                yw = y1;
            else
                yw = y2;
        else
            if(!invert_function )
                yw = y2;
            else
                yw = y1;
    }
    else
    {
        if( wx1+wx2 >0 ) xw = (x1*wx1+x2*wx2)/(wx1+wx2);
        if( wy1+wy2 >0 ) yw = (y1*wy1+y2*wy2)/(wy1+wy2);
    }


    // # warning SimpleFit:: Non realistic very BAD SigmaX,Y  setting - to be valid for muons.

    // Calculate sigmas
    double sqe = sqrt(Etotal);
    double sigmax = SizeX/2.;
    double sigmay = SizeY/2.;
    if( pid <= 3 )   // gamma,electron,positron
    {
        sigmax = 6./sqe + fabs(xxb)/10.;
        sigmay = 6./sqe + fabs(yyb)/10.;
    }

    // Set X,Y position of cluster
    double xpart = xw+x0;
    double ypart = yw+y0;
    partResponse.GetParticleToModify().SetX( xpart , sigmax );
    partResponse.GetParticleToModify().SetY( ypart , sigmay );

    // correction of Z can be ommited after tests
    double SizeZ = p_main_cell->GetCellType().GetSizeZ()/2.;
    double z0 = p_main_cell->GetZ()-SizeZ;
    double zw = SizeZ;
    double sigmaz = 0;
    if( pid <= 3 )   // gamma,electron,positron
    {
        zw = ZmidShowerEM(Etotal9)*RadLeng;
        sigmaz = SigmaZmidShowerEM(Etotal9)*RadLeng;
    }
    else if( pid == 5 || pid == 6  )   // muon
    {
        zw = SizeZ;
        sigmaz = SizeZ/1.732;
    }
    else if( pid >= 7  )   // hadron (pion_0 can not appear)
    {
        zw = ZmidShowerHadronic(Etotal9)*NuclLeng;
        sigmaz= SigmaZmidShowerHadronic(Etotal9)*NuclLeng;
    }
    partResponse.GetParticleToModify().SetZ( zw+z0 , sigmaz );


    //  Try some development with "electronic delta recovery"

    bool debug_electronic_delta_recovery = false;

    // Calculate correction for total energy
    double elost_integer = 0.;
    double elost_delta = 0.;
    double main_ecut_delta = 0.;
    double cf_main = 0.;
    for (size_t i=0; i<partResponse.GetCells().size(); i++) {
        double  ercell = partResponse.GetCellsEnergy()[i];
        double  etcell = partResponse.GetExpectedCellsEnergy()[i];
        int jcell = partResponse.GetCells()[i];
        double ecutdelta = GetCalorimeter()->GetEnergyCutSparseMode()[jcell];
        double cf = GetCalorimeter()->GetCellInfo(Calorimeter::CALIB, Calorimeter::OLD, jcell).GetMean();
        if(i == 0 )
        {
            main_ecut_delta = ecutdelta;
            cf_main = cf;
        }
        if( debug_electronic_delta_recovery ) std::cout <<"  " << i <<" Er " << ercell <<" Et " << etcell <<
            " Edelta " << ecutdelta << std::endl;
        if( ercell > 0. && ercell > ecutdelta ) elost_integer += cf;
        if( ercell == 0. && etcell < ecutdelta+cf ) elost_delta += etcell;
    }
    elost_integer /= 2.;
    elost_delta += main_ecut_delta/3.;


    // # warning SimpleFit::   Shower leak corrections valid only for gamma!! Need to be improved!
    if( debug_electronic_delta_recovery )
    {
        std::cout << " Etotal " << Etotal << " elost_integer " << elost_integer << " elost_delta " << elost_delta << std::endl;
    }

    // Correct total energy
    if(GetCalorimeter()->GetOptions().correct_for_digitization)
        Etotal += elost_integer + elost_delta + cf_main/3.;

    // Set particle energy
    double fraction = 1.+GetCalorimeter()->GetOptions().energy_leak_for_SimpleFit;
    double eg = fraction*Etotal;
    double se = GetCalorimeter()->EnergySigmaInCell(eg, main_cell);
    partResponse.GetParticleToModify().SetE( eg , se );

    if( debug ) std::cout << " After SimpleFit debug E = " << partResponse.GetParticle().GetE() << std::endl;
}

////////////////////////////////////////////////////////////////////////////////

/*!
  \brief Determine particle position using iterative algorithm
  \callgraph
  \callergraph
  \remarks {
  Hardcoded:
  bool debug : Debug output
  bool select_to_print : More debug output
  }
*/
void ReconstructionKolosov::Fit(OneParticleResponse& partResponse) {

    bool debug = false;

    bool select_to_print = false;

    if( debug )
    {
        if( partResponse.GetParticle().GetE() > 5. ) select_to_print = true;
    }

    if( partResponse.GetCells().empty() )
        return;


    if( debug && select_to_print )
    {
        std::cout << " START FIT " << std::endl;
        std::cout << "  START E=" << partResponse.GetParticle().GetE() << std::endl;
        std::cout << "  Fit START X=" << partResponse.GetParticle().GetX() << " Y=" << partResponse.GetParticle().GetY() << std::endl;
    }

    int maxiter = 25;
    double sss_old = 0.;
    for( int iter=0; iter != maxiter; iter++ )
    {
        double sss_epsilon = 0.01;

        // Init temporary variables
        double sfe(0.), sfx(0.), sfy(0.);
        double see(0.), sxe(0.), sye(0.);
        double sxx(0.), syx(0.), syy(0.);
        double sss(0.);
        double de(0.), we(0.);
        double dei(0.), dxi(0.), dyi(0.);



        // Calculate expected respons and derivation from expexted response
        partResponse.CalculateDerivativesExpectedResponse();

        // Loop over cells in cluster
        for (size_t i=0; i<partResponse.GetCells().size(); i++) {

            // Get cell index
            int jcell = partResponse.GetCells()[i];

            // Calculate derivation from expected energy to real energy
            de = partResponse.GetCellsEnergy()[i] - partResponse.GetExpectedCellsEnergy()[i];

            // Get error of real energy measurment
            double de_ro =  real_amplitudes[jcell].second;
            // Increase error of real energy measurment for bad cells
            if( GetCalorimeter()->CellIsBad(jcell, Calorimeter::OLD) > 0 )
                de_ro = 4.*de_ro;

            // Calculate weight
            we = 1./(partResponse.GetCellsDispEnergy()[i] + de_ro);

            // Correct for digitalization effects (idf desired)
            if(GetCalorimeter()->GetOptions().correct_for_digitization)
            {
                double ecutdelta = GetCalorimeter()->GetEnergyCutSparseMode()[jcell];
                if( partResponse.GetCellsEnergy()[i] < ecutdelta )
                    we = 1./(partResponse.GetCellsDispEnergy()[i] + de_ro + ecutdelta*ecutdelta );
            }

            // Get derivated energy errors (CalculateDerivativesExpectedResponse())
            dei = partResponse.GetDDE(partResponse.GetCells()[i]);
            dxi = partResponse.GetDDX(partResponse.GetCells()[i]);
            dyi = partResponse.GetDDY(partResponse.GetCells()[i]);

            // Add values to sigmas
            sfe += de *dei*we;
            sfx += de *dxi*we;
            sfy += de *dyi*we;

            see += dei*dei*we;
            sxe += dxi*dei*we;
            sye += dyi*dei*we;

            sxx += dxi*dxi*we;
            syx += dyi*dxi*we;
            syy += dyi*dyi*we;

            sss += de*de*we;
        }

        {
            // Fit parameters 0 -> e, 1 -> x, 2 -> y
            double a[3][3] = {
                { see, sxe, sye },
                { sxe, sxx, syx },
                { sye, syx, syy }
            };
            double b[3] = {sfe, sfx, sfy} ;
            {
                // Variables for external call
                int nf(3), nnf(3), neq(3);
                int ifail(0);
                int ir[1000];

                // Call external function (Matrix equation solver)
                deqinv_(nf,a,nnf,ir,ifail,neq,b);  // Fortran CERNLIB

                if( ifail != 0 )
                    std::cout << " deqinv_ ifail= " << ifail << std::endl;
            }

            // Get changed parameter
            sfe = b[0];
            sfx = b[1];
            sfy = b[2];

        }

        // Apply cut off set to maxima
        if( sfx >  5. ) sfx = 5.;
        if( sfx < -5. ) sfx =-5.;
        if( sfy >  5. ) sfy = 5.;
        if( sfy < -5. ) sfy =-5.;
        if( sfe >  0.1*partResponse.GetParticle().GetE() ) sfe = 0.1*partResponse.GetParticle().GetE();
        if( sfe < -0.1*partResponse.GetParticle().GetE() ) sfe = -0.1*partResponse.GetParticle().GetE();
        //#warning OneParticleResponse::Fit Parmeters errors not implemented yet
        if( partResponse.GetParticle().GetE()+sfe < 0.05 )
            partResponse.GetParticleToModify().SetE( 0.05 , 0.5 );
        else
            partResponse.GetParticleToModify().SetE( partResponse.GetParticle().GetE()+sfe , 0.5 );

        partResponse.GetParticleToModify().SetX( partResponse.GetParticle().GetX()+sfx , 3. );
        partResponse.GetParticleToModify().SetY( partResponse.GetParticle().GetY()+sfy , 3. );

        if( iter > 0 && 1.- sss/sss_old < sss_epsilon )
            break;

        sss_old = sss;

    }

    // Shower leak corrections caused by wrong Shower Profile parametrisation
    //#warning OneParticleResponse::Fit Tune energy by hands NOW +6%
    //      partResponse.GetParticleToModify().SetE( 1.06*partResponse.GetParticle().GetE() , 0.5 );

    double eg = partResponse.GetParticle().GetE();


    // Calculate energy correction
    double elost_integer = 0.;
    double main_ecut_delta = 0.;
    double cf_main = 0.;
    for (size_t i=0; i< partResponse.GetCells().size(); i++) {

        double  ercell = partResponse.GetCellsEnergy()[i];
        int jcell = partResponse.GetCells()[i];
        double ecutdelta = GetCalorimeter()->GetEnergyCutSparseMode()[jcell];
        double cf = GetCalorimeter()->GetCellInfo(Calorimeter::CALIB, Calorimeter::OLD, jcell).GetMean();
        if(i == 0 )
        {
            main_ecut_delta = ecutdelta;
            cf_main = cf;
        }
        if( ercell > 0. && ercell > ecutdelta ) elost_integer += cf;
    }
    elost_integer /= 2.;

    // correct energy if desired
    if(GetCalorimeter()->GetOptions().correct_for_digitization)
        eg += main_ecut_delta/2.+elost_integer+cf_main/3.;

    // Set particle energy
    partResponse.GetParticleToModify().SetE( eg , 0.5 );

    if( debug && select_to_print )
    {
        std::cout << " FINISHED FIT " << std::endl;
        std::cout << "  START E=" << partResponse.GetParticle().GetE() << std::endl;
        std::cout << "  Fit START X=" << partResponse.GetParticle().GetX() << " Y=" << partResponse.GetParticle().GetY() << std::endl;
    }
}

////////////////////////////////////////////////////////////////////////////////

} // namespace Reco

