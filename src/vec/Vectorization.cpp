/*
  Copyright (C) 2021
      Max Planck Institute for Polymer Research & JGU Mainz

  This file is part of ESPResSo++.

  ESPResSo++ is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.

  ESPResSo++ is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/

#include "Vectorization.hpp"
#include "vec/integrator/MDIntegratorVec.hpp"
#include "vec/storage/StorageVec.hpp"

#include "storage/Storage.hpp"

namespace espressopp
{
namespace vec
{
///////////////////////////////////////////////////////////////////////////////////////////////
LOG4ESPP_LOGGER(Vectorization::logger, "Vectorization");

///////////////////////////////////////////////////////////////////////////////////////////////
/// constructor
Vectorization::Vectorization(std::shared_ptr<System> _system,
                             std::shared_ptr<MDIntegrator> _mdintegrator)
    : SystemAccess(_system), mdintegrator(_mdintegrator), vecLevel(1)
{
    connect();
    resetCells();  // immediately retrieve cell information
}

Vectorization::Vectorization(std::shared_ptr<System> _system) : SystemAccess(_system), vecLevel(2)
{
}

///////////////////////////////////////////////////////////////////////////////////////////////
/// destructor
Vectorization::~Vectorization() { disconnect(); }

///////////////////////////////////////////////////////////////////////////////////////////////
/// connect to boost signals in integrator and storage
void Vectorization::connect()
{
    sigResetParticles = getSystem()->storage->onParticlesChanged.connect(
        boost::signals2::at_front,  // call first due to reordering
        std::bind(&Vectorization::resetParticles, this));
    sigResetCells = getSystem()->storage->onCellAdjust.connect(
        boost::signals2::at_back, std::bind(&Vectorization::resetCells, this));
    sigBefCalcForces = mdintegrator->aftInitF.connect(
        boost::signals2::at_back, std::bind(&Vectorization::befCalcForces, this));
    sigUpdateForces = mdintegrator->aftCalcFLocal.connect(
        boost::signals2::at_front, std::bind(&Vectorization::updateForces, this));
}

///////////////////////////////////////////////////////////////////////////////////////////////
/// disconnect boost signals made with connect()
void Vectorization::disconnect()
{
    sigResetParticles.disconnect();
    sigResetCells.disconnect();
    sigBefCalcForces.disconnect();
    sigUpdateForces.disconnect();
}

///////////////////////////////////////////////////////////////////////////////////////////////
/// reset and update cell mapping and neighbor lists from current storage
void Vectorization::resetCells()
{
    if (!getSystem()->storage)
    {
        throw std::runtime_error("System has no storage");
    }
    resetCellsStorage(getSystem()->storage.get());
}

void Vectorization::resetCellsStorage(Storage* storage)
{
    CellList const& localCells = storage->getLocalCells();
    CellList const& realCells = storage->getRealCells();
    Cell* const cell0 = localCells[0];
    const auto realCellIdx = CellListToIdx(realCells, cell0);
    particles.markRealCells(realCellIdx, localCells.size());
    neighborList = CellNeighborList(cell0, localCells, realCellIdx);
    LOG4ESPP_TRACE(logger, "neighborList, ncells: " << neighborList.numCells() << " nnbrs: "
                                                    << neighborList.maxNumNeighbors());
}

///////////////////////////////////////////////////////////////////////////////////////////////
/// set force array/s to zero
void Vectorization::zeroForces()
{
    auto f_zero = [](real* __restrict f, size_t size) {
#ifdef __INTEL_COMPILER
#pragma vector always
#pragma vector aligned
#pragma ivdep
#endif
        for (size_t i = 0; i < size; i++) f[i] = 0.0;
    };
    {
        f_zero(particles.f_x.data(), particles.f_x.size());
        f_zero(particles.f_y.data(), particles.f_y.size());
        f_zero(particles.f_z.data(), particles.f_z.size());
    }
}

///////////////////////////////////////////////////////////////////////////////////////////////
/// reset and update particles from current storage
void Vectorization::resetParticles() { particles.copyFrom(getSystem()->storage->getLocalCells()); }

///////////////////////////////////////////////////////////////////////////////////////////////
/// set force array/s to zero and overwrite particles position data
void Vectorization::befCalcForces()
{
    zeroForces();
    particles.updateFromPositionOnly(getSystem()->storage->getLocalCells());
}

///////////////////////////////////////////////////////////////////////////////////////////////
/// add forces back to storage
void Vectorization::updateForces()
{
    particles.addToForceOnly(getSystem()->storage->getLocalCells());
}

///////////////////////////////////////////////////////////////////////////////////////////////
/// registration with python
void Vectorization::registerPython()
{
    using namespace espressopp::python;

    class_<Vectorization, std::shared_ptr<Vectorization>>(
        "vec_Vectorization", init<std::shared_ptr<System>, std::shared_ptr<MDIntegrator>>())
        .def(init<std::shared_ptr<System>>())
        .add_property("level", &Vectorization::getVecLevel)
        .def_readwrite("storageVec", &Vectorization::storageVec);
}
}  // namespace vec
}  // namespace espressopp