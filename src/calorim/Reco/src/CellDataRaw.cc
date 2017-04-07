/*
   $Source: /afs/cern.ch/user/g/gobbo/w0/coral/src/calorim/Reco/src/CellDataRaw.cc,v $
   $Date: 2010/04/20 15:12:20 $
   $Revision: 1.2 $
*/

// --- Standard C/C++ library ---
#include <cassert>
#include <cmath>
#include <iostream>

// --- Internal files ---
#include "CellDataRaw.h"

namespace Reco {

////////////////////////////////////////////////////////////////////////////////

CellDataRaw::CellDataRaw(const size_t& cellIdx) :
    fCellIdx(cellIdx),
    fHasAmplitude(false),
    fHasEnergy(false), fHasEnergyErr(false),
    fHasTime(false), fHasTimeErr(false) {
}

////////////////////////////////////////////////////////////////////////////////

CellDataRaw::CellDataRaw(const size_t& cellIdx, const double& energy, const double& time) :
    fCellIdx(cellIdx),
    fHasAmplitude(false),
    fEnergy(energy),
    fHasEnergy(true), fHasEnergyErr(false),
    fTime(time),
    fHasTime(true), fHasTimeErr(false) {
}

////////////////////////////////////////////////////////////////////////////////

CellDataRaw::CellDataRaw(const size_t& cellIdx, const CellDataRaw& signal) {
    *this = signal;
    fCellIdx = cellIdx;
}

////////////////////////////////////////////////////////////////////////////////

const size_t& CellDataRaw::GetCellIdx() const {
    return fCellIdx;
}

////////////////////////////////////////////////////////////////////////////////

const double& CellDataRaw::GetAmplitude() const {
    if (!fHasAmplitude)
        std::cerr << "CellDataRaw::GetAmplitude: amplitude not initialized" << std::endl;

    return fAmplitude;
}

////////////////////////////////////////////////////////////////////////////////

const double& CellDataRaw::GetEnergy() const {
    if (!fHasEnergy)
        std::cerr << "CellDataRaw::GetEnergy: energy not initialized" << std::endl;

    return fEnergy;
}

////////////////////////////////////////////////////////////////////////////////

const double& CellDataRaw::GetEnergyErr() const {
    if (!fHasEnergyErr)
        std::cerr << "CellDataRaw::GetEnergyErr: energy error not initialized" << std::endl;

    return fEnergyErr;
}

////////////////////////////////////////////////////////////////////////////////

const double& CellDataRaw::GetTime() const {
    if (!fHasTime)
        std::cerr << "CellDataRaw::GetTime: time not initialized" << std::endl;

    return fTime;
}

////////////////////////////////////////////////////////////////////////////////

const double& CellDataRaw::GetTimeErr() const {
    if (!fHasTimeErr)
        std::cerr << "CellDataRaw::GetTimeErr: time error not initialized" << std::endl;

    return fTimeErr;
}

////////////////////////////////////////////////////////////////////////////////

const bool& CellDataRaw::HasAmplitude() const {
    return fHasAmplitude;
}

////////////////////////////////////////////////////////////////////////////////

const bool& CellDataRaw::HasEnergy() const {
    return fHasEnergy;
}

////////////////////////////////////////////////////////////////////////////////

const bool& CellDataRaw::HasEnergyErr() const {
    return fHasEnergyErr;
}

////////////////////////////////////////////////////////////////////////////////

const bool& CellDataRaw::HasTime() const {
    return fHasTime;
}

////////////////////////////////////////////////////////////////////////////////

const bool& CellDataRaw::HasTimeErr() const {
    return fHasTimeErr;
}

////////////////////////////////////////////////////////////////////////////////

void CellDataRaw::SetAmplitude(const double& amplitude) {
    fAmplitude = amplitude;
    fHasAmplitude = true;
}

////////////////////////////////////////////////////////////////////////////////

void CellDataRaw::SetEnergy(const double& energy) {
    fEnergy = energy;
    fHasEnergy = true;
}

////////////////////////////////////////////////////////////////////////////////

void CellDataRaw::SetEnergyErr(const double& energyErr) {
    fEnergyErr = energyErr;
    fHasEnergyErr = true;
}

////////////////////////////////////////////////////////////////////////////////

void CellDataRaw::SetTime(const double& time) {
    fTime = time;
    fHasTime = true;
}

////////////////////////////////////////////////////////////////////////////////

void CellDataRaw::SetTimeErr(const double& timeErr) {
    fTimeErr = timeErr;
    fHasTimeErr = true;
}

////////////////////////////////////////////////////////////////////////////////

} // namespace Reco
