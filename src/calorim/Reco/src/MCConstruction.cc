/*
   $Source: /afs/cern.ch/user/g/gobbo/w0/coral/src/calorim/Reco/src/MCConstruction.cc,v $
   $Date: 2011/02/01 22:05:52 $
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

// --- Standard C/C++ library ---
#include <cassert>
#include <iostream>
#include <cstdlib>
#include <cmath>

// --- Internal files ----
#include "MCConstruction.h"

#include "Calorimeter.h"
#include "CellDataRaw.h"
#include "CellType.h"
#include "OneParticleResponse.h"

using namespace std;

namespace Reco {

////////////////////////////////////////////////////////////////////////////////

MCConstruction::MCConstruction(const Calorimeter* c) :
    calorimeter_(c) {
    expected_amplitudes.resize(c->NCells(), pair<double, double>(0., 0.));
    expected_times     .resize(c->NCells(), 0.);
}

////////////////////////////////////////////////////////////////////////////////

const std::vector<CellDataRaw>& MCConstruction::DoDigitization(const std::vector<CalorimeterParticle>& particles) {
//    bool debug = true;
    bool debug = false;

    if ( debug )
        cout << "MCConstruction::DoDigitization " << GetCalorimeter()->GetName() << " debug " << endl;

    ClearExpected();

    for (vector<CalorimeterParticle>::const_iterator it=particles.begin(); it!=particles.end(); it++) {
        OneParticleResponse response((*it), GetCalorimeter());

        response.MonteRandomise();

        AddToExpected(response, it->GetTime());
    }

    ExtractExpected(0.);

    return signals_;
}

////////////////////////////////////////////////////////////////////////////////

const std::vector<CellDataRaw>& MCConstruction::DoConstruction(const std::vector<CellDataRaw>& mc_input) {
//    bool debug = true;
    bool debug = false;
    if ( debug )
        cout << "MCConstruction::DoConstruction " << GetCalorimeter()->GetName()
             << " debug real= " << GetCalorimeter()->GetOptions().mc_make_real_digitization << endl;

    ImportExpected(mc_input);

    if ( GetCalorimeter()->GetOptions().mc_smear_response ) {
        if ( GetCalorimeter()->GetOptions().mc_smear_constant   ) AddConstantFluctuations();
        if ( GetCalorimeter()->GetOptions().mc_smear_stochastic ) AddStochasticFluctuations();
        if ( GetCalorimeter()->GetOptions().mc_smear_readout    ) AddReadoutFluctuations();
    }

    bool real_digitization = GetCalorimeter()->GetOptions().mc_make_real_digitization;

    if ( real_digitization ) {
        if ( !GetCalorimeter()->GetOptions().mc_make_fiadc_digitization )
            cerr << "Warning!! in MCConstruction::DoConstruction " << GetCalorimeter()->GetName()
                 << " options.mc_make_fiadc_digitization is disabled ! But sorry there in no alternative now !! " << endl;

        ExtractExpected(GetCalorimeter()->GetOptions().mc_data_extraction_threshold);

        // FIADC elctronics: asuming sparse mode
        int fiadc_delta = (int)GetCalorimeter()->GetOptions().mc_fiadc_sparce_delta;
        std::vector<CellDataRaw> signals_new;
        for (vector<CellDataRaw>::const_iterator it=signals_.begin(); it!=signals_.end(); it++) {
            size_t icell = it->GetCellIdx();
            double e     = it->GetEnergy();
            double e_dig = e/GetCalorimeter()->GetCellInfo(Calorimeter::CALIB, Calorimeter::MC, icell).GetMean();

            // Here we should add electronics noise
            int idig = (int)e_dig;
            if ( idig > fiadc_delta ) {
                if( icell < GetCalorimeter()->NCells() ) {
	                double e_digit = idig;
                    // Amplitudes from calorimeter cells
                    signals_new.push_back(*it);

                    signals_new.back().SetAmplitude( e_digit );
                    signals_new.back().SetEnergy   ( GetCalorimeter()->SignalToEnergy(icell, e_digit, GetCalorimeter()->GetTimeInSpill()) );
                } else{
                    cerr << "MCConstruction::DoConstruction " << GetCalorimeter()->GetName() << " wrong cell " << icell << endl;
	            }
            }
        }

        signals_ = signals_new;
        if ( debug )
            cout << "MCConstruction::DoConstruction data before FilterOutBadCells " << signals_.size() << endl;

        FilterOutBadCells();

        if ( debug )
            cout << "MCConstruction::DoConstruction data after FilterOutBadCells " << signals_.size() << endl;

    } else {
        if ( debug )
            cout << "  Dont Make real Digitization " << endl;

        ExtractExpected(GetCalorimeter()->GetOptions().mc_data_extraction_threshold);

        FilterOutBadCells();

        if ( debug )
            cout << "  Filtered Data  " << signals_.size() << endl;

        for (vector<CellDataRaw>::iterator it=signals_.begin(); it!=signals_.end(); it++ ) {
            size_t icell = it->GetCellIdx();
            if ( icell>=GetCalorimeter()->NCells() )
                throw Exception("MCConstruction::DoConstruction for %s: cell Address %d is out of calorimeter size %d",
                                GetCalorimeter()->GetName().c_str(), icell, GetCalorimeter()->NCells());

            it->SetEnergy( it->GetEnergy()
                           * GetCalorimeter()->GetCellInfo(Calorimeter::CALIB, Calorimeter::OLD, icell).GetMean()
                           / GetCalorimeter()->GetCellInfo(Calorimeter::CALIB, Calorimeter::MC, icell).GetMean() );

            if ( debug )
                cout << "  Amplitude " << it->GetEnergy()
                     << " Cf OlD " << GetCalorimeter()->GetCellInfo(Calorimeter::CALIB, Calorimeter::OLD, icell).GetMean()
                     << " Cf MC "  << GetCalorimeter()->GetCellInfo(Calorimeter::CALIB, Calorimeter::MC,  icell).GetMean()
                     << " Cf default " << GetCalorimeter()->GetOptions().default_calibration << endl;
        }
    }

    return signals_;
}

////////////////////////////////////////////////////////////////////////////////

double MCConstruction::AddParticle(const CalorimeterParticle& particle) {
    bool development = true;
    OneParticleResponse response(particle, GetCalorimeter(), development);

    response.MonteRandomise();

    AddToExpected(response, particle.GetTime());

    return response.GetEnergyDeposit();
}

////////////////////////////////////////////////////////////////////////////////

void MCConstruction::AddConstantFluctuations(void) {
    for (size_t it=0; it<GetCalorimeter()->NCells(); it++) {
// TODO: Monte-Carlo options - Hard coded 100 MeV threshold to randomise Constant term.
        if ( expected_amplitudes[it].first > 0.1 ) {
// TODO: Usage of non CORAL random methods does not seem to be right for reproduceable results
            const double r1 = GetCalorimeter()->GetRandomGaus();

            double Ecell      = expected_amplitudes[it].first;
            double SigmaEcell = Ecell * GetCalorimeter()->GetCells()[it].GetCellType().GetConstantTerm();

            expected_amplitudes[it].first  += r1*SigmaEcell;
            expected_amplitudes[it].second += SigmaEcell*SigmaEcell;

            if (expected_amplitudes[it].first < 0.)
                expected_amplitudes[it].first = 0.;
        }
    }
}

////////////////////////////////////////////////////////////////////////////////

void MCConstruction::AddStochasticFluctuations(void) {
    for (size_t it=0; it<GetCalorimeter()->NCells(); it++) {
// TODO: Monte-Carlo options - Hard coded 5 MeV threshold to randomise photo-electrons light.
        if(expected_amplitudes[it].first > 0.005 ) {
            const double r1 = GetCalorimeter()->GetRandomGaus();

            double Ecell      = expected_amplitudes[it].first;
            double SigmaEcell = sqrt(Ecell) * GetCalorimeter()->GetCells()[it].GetCellType().GetStochasticTerm();

            expected_amplitudes[it].first  += r1*SigmaEcell;
            expected_amplitudes[it].second += SigmaEcell*SigmaEcell;

            if (expected_amplitudes[it].first < 0.)
                expected_amplitudes[it].first = 0.;
        }
    }
}

////////////////////////////////////////////////////////////////////////////////

void MCConstruction::AddReadoutFluctuations(void) {
    for (size_t it=0; it<GetCalorimeter()->NCells(); it++) {
// TODO: Monte-Carlo options - Check ReadOut term.
        const double r1 = GetCalorimeter()->GetRandomGaus();

        double SigmaEcell = GetCalorimeter()->GetCells()[it].GetCellType().GetReadoutTerm();

        if (SigmaEcell > 1.0)
            cout << "Check in:: " << GetCalorimeter()->GetName() << " Very Big ReadOut Term !! " << SigmaEcell << endl;

        expected_amplitudes[it].first  += r1*SigmaEcell;
        expected_amplitudes[it].second += SigmaEcell*SigmaEcell;

        if (expected_amplitudes[it].first < 0.)
            expected_amplitudes[it].first = 0.;
    }
}

////////////////////////////////////////////////////////////////////////////////

void Calorimeter::AddCalibrationSmearing( double s )
{
  double u1,u2,r1,r2,ss;
  ss=s;
  ss=fabs(ss);
  for( size_t i=0; i<NCells(); i++ )
  {
    const double r2 = GetRandomGaus();
    double stat = cells_info[CALIB][OLD][i].GetEntries();
    double mean = (1.+ss*r2)*cells_info[CALIB][OLD][i].GetMean();
    mean = fabs(mean);
    double sigma = cells_info[CALIB][OLD][i].GetSigma();
    cells_info[CALIB][OLD][i].Set(stat, mean, sigma);
  }
}

////////////////////////////////////////////////////////////////////////////////

void MCConstruction::ClearExpected() {
    for (size_t it=0; it<GetCalorimeter()->NCells(); it++) {
        expected_amplitudes[it] = pair<double,double>(0,0);
        expected_times[it]      = 0.;
    }
}

////////////////////////////////////////////////////////////////////////////////

void MCConstruction::AddToExpected(const OneParticleResponse& response, const double& time) {
    size_t i =0;
    for (vector<size_t>::const_iterator it=response.GetCells().begin();
         it!=response.GetCells().end(); it++ ) {
        expected_times[*it] = (expected_times[*it]*expected_amplitudes[*it].first + time*response.GetExpectedCellsEnergy()[i])
                                       / (expected_amplitudes[*it].first + response.GetExpectedCellsEnergy()[i]);

        expected_amplitudes[*it].first  += response.GetExpectedCellsEnergy()[i];
        expected_amplitudes[*it].second += response.GetExpectedCellsDispEnergy()[i];

        i++;
    }
}

////////////////////////////////////////////////////////////////////////////////

void MCConstruction::ExtractExpected(const double& extraction_energy_threshold) {
    signals_.clear();

    for (size_t it=0; it<GetCalorimeter()->NCells(); it++)
        if (expected_amplitudes[it].first > extraction_energy_threshold) {
// TODO: hardcoded 5ns time resolution
            signals_.push_back( CellDataRaw(it,
                                       expected_amplitudes[it].first,
                                       expected_times[it]) );
            signals_.back().SetEnergyErr( expected_amplitudes[it].second );
        }

}

////////////////////////////////////////////////////////////////////////////////

void MCConstruction::ImportExpected(const vector<CellDataRaw>& mc_input) {
    ClearExpected();

    for (vector<CellDataRaw>::const_iterator it=mc_input.begin(); it!=mc_input.end(); it++) {
        // protection against unlogical energies
        if (it->GetEnergy()<=0.)
            continue;

        expected_times     [it->GetCellIdx()]         = (expected_times[it->GetCellIdx()]*expected_amplitudes[it->GetCellIdx()].first
                                                         + it->GetTime()*it->GetEnergy() )
                                                        / (expected_amplitudes[it->GetCellIdx()].first + it->GetEnergy());
        expected_amplitudes[it->GetCellIdx()].first  += it->GetEnergy();

        // this is allowed not to happen:
        // for MC imported from CORAL there is no initial error on the energy,
        // it is only added during further processing
        if (it->HasEnergyErr())
            expected_amplitudes[it->GetCellIdx()].second = sqrt( expected_amplitudes[it->GetCellIdx()].second*expected_amplitudes[it->GetCellIdx()].second
                                                                 + it->GetEnergyErr()*it->GetEnergyErr() );
    }
}

////////////////////////////////////////////////////////////////////////////////

void MCConstruction::FilterOutBadCells() {
    vector<CellDataRaw> signals_new;
    for (vector<CellDataRaw>::iterator it=signals_.begin(); it!=signals_.end(); it++)
        if ( !GetCalorimeter()->CellIsBad(it->GetCellIdx(), Calorimeter::MC) )
            signals_new.push_back(*it);

    signals_ = signals_new;
}

////////////////////////////////////////////////////////////////////////////////

} /// using namespace Reco
