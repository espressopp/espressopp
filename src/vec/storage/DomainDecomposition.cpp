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

#include "vec/storage/DomainDecomposition.hpp"
#include "vec/Vectorization.hpp"
#include "bc/BC.hpp"

const int DD_COMM_TAG = 0xab;

namespace espressopp
{
namespace vec
{
namespace storage
{
LOG4ESPP_LOGGER(DomainDecomposition::logger, "DomainDecomposition");

const Real3D DomainDecomposition::SHIFT_ZERO{0., 0., 0.};

DomainDecomposition::DomainDecomposition(std::shared_ptr<System> system,
                                         const Int3D& _nodeGrid,
                                         const Int3D& _cellGrid,
                                         int _halfCellInt)
    : baseClass(system, _nodeGrid, _cellGrid, _halfCellInt), StorageVec(system)
{
    if (halfCellInt != 1) throw std::runtime_error("vec: Not implemented for halfCellInt!=1.");
    // resetTimers();
    connect();
    resetCells();
    LOG4ESPP_INFO(logger, "DomainDecomposition()");
}

DomainDecomposition::~DomainDecomposition() { disconnect(); }

void DomainDecomposition::connect()
{
    if (!sigResetCells.connected())
        sigResetCells = this->onCellAdjust.connect(
            boost::signals2::at_back, std::bind(&DomainDecomposition::resetCells, this));

    if (!sigLoadCells.connected())
        sigLoadCells = this->onParticlesChanged.connect(
            boost::signals2::at_front, std::bind(&DomainDecomposition::loadCells, this));

    LOG4ESPP_INFO(logger, "DomainDecomposition::connect()");
}

void DomainDecomposition::disconnect()
{
    if (sigResetCells.connected()) sigResetCells.disconnect();
    if (sigLoadCells.connected()) sigLoadCells.disconnect();
    LOG4ESPP_INFO(logger, "DomainDecomposition::disconnect()");
}

/// Copy particles to packed form. To be called at the start of integrator.run
void DomainDecomposition::loadCells()
{
    vectorization->resetParticles();
    if (localParticlesEnabled) localParticlesVec.rebuild(vectorization->particles, uniqueCells);

    prepareGhostBuffers();
}

/// Copy particles back from packed form. To be called at the end of integrator.run
void DomainDecomposition::unloadCells()
{
    vectorization->particles.updateToPositionVelocity(localCells, true);
}

void DomainDecomposition::resetCells()
{
    /// reset realCells and cellNeighborList
    vectorization->resetCellsStorage(this);

    /// Setup ghost communication for this cell grid
    {
        /// convert commCells to commCellIdx
        for (int dir = 0; dir < 6; dir++)
        {
            commCellIdx[dir].reals = CellListToIdx(commCells[dir].reals, localCells[0]);
            commCellIdx[dir].ghosts = CellListToIdx(commCells[dir].ghosts, localCells[0]);
        }
        /// TODO: Check whether sorting improves performance
    }

    /// Determine list of unique cells for this subdomain. Unique cells comprise of real cells
    /// and ghost cells which are not derived from the same subdomain (intra).
    {
        uniqueCells.clear();
        uniqueCells.reserve(localCells.size());
        std::vector<size_t> cells(localCells.size(), 0);

        for (int coord = 0; coord < 3; ++coord)
        {
            const bool doPeriodic = (nodeGrid.getGridSize(coord) == 1);
            for (int lr = 0; lr < 2; ++lr)
            {
                int const dir = 2 * coord + lr;
                if (doPeriodic)
                {
                    auto const& ghosts = commCellIdx[dir].ghosts;
                    for (auto const& gc : ghosts)
                    {
                        cells[gc]++;
                    }
                }
            }
        }

        for (size_t ic = 0; ic < cells.size(); ic++)
        {
            if (cells[ic] == 0) uniqueCells.push_back(ic);
        }
    }
}

///////////////////////////////////////////////////////////////////////////////////////////////
/// preallocate buffers for ghost communication
void DomainDecomposition::prepareGhostBuffers()
{
    const auto& cr = vectorization->particles.cellRange();
    size_t maxReals = 0, maxGhosts = 0;
    for (size_t coord = 0; coord < 3; ++coord)
    {
        for (size_t lr = 0; lr < 2; ++lr)
        {
            size_t const dir = 2 * coord + lr;

            auto f_countParticles = [cr](auto const& cellIdx) {
                size_t total = 0;
                for (auto const& ic : cellIdx) total += (cr[ic + 1] - cr[ic]);
                return total;
            };

            auto& cci = commCellIdx[dir];
            cci.numReals = f_countParticles(cci.reals);
            cci.numGhosts = f_countParticles(cci.ghosts);
            maxReals = std::max(cci.numReals, maxReals);
            maxGhosts = std::max(cci.numGhosts, maxGhosts);
        }
    }

    const size_t preallocReal = maxReals * 3;
    const size_t preallocGhost = maxGhosts * 3;
    buffReal.resize(preallocReal);
    buffGhost.resize(preallocGhost);
}

///////////////////////////////////////////////////////////////////////////////////////////////
void DomainDecomposition::updateGhostsVec() { ghostCommunication_impl<false, true, 0>(); }

///////////////////////////////////////////////////////////////////////////////////////////////
void DomainDecomposition::collectGhostForcesVec() { ghostCommunication_impl<false, false, 0>(); }

///////////////////////////////////////////////////////////////////////////////////////////////
template <bool SIZES_FIRST, bool REAL_TO_GHOSTS, int EXTRA_DATA>
void DomainDecomposition::ghostCommunication_impl()
{
    auto const& comm = *(baseClass::getSystem()->comm);
    auto const boxL = baseClass::getSystem()->bc->getBoxL();

    for (size_t _coord = 0; _coord < 3; ++_coord)
    {
        int coord = REAL_TO_GHOSTS ? _coord : (2 - _coord);
        const real curCoordBoxL = boxL[coord];
        const bool doPeriodic = (nodeGrid.getGridSize(coord) == 1);

        for (size_t lr = 0; lr < 2; ++lr)
        {
            size_t const dir = 2 * coord + lr;
            size_t const oppDir = 2 * coord + (1 - lr);

            Real3D shift(0, 0, 0);
            if (REAL_TO_GHOSTS)
            {
                shift[coord] = nodeGrid.getBoundary(dir) * curCoordBoxL;
            }

            if (doPeriodic)
            {
                if (REAL_TO_GHOSTS)
                {
                    copyRealsToGhostsIntra<ADD_SHIFT>(dir, shift);
                }
                else
                {
                    addGhostForcesToRealsIntra(dir);
                }
            }
            else
            {
                if (REAL_TO_GHOSTS)
                {
                    packCells<PACKED_POSITIONS, ADD_SHIFT>(buffReal, true, dir, shift);
                }
                else
                {
                    packCells<PACKED_FORCES, NO_SHIFT>(buffGhost, false, dir, SHIFT_ZERO);
                }

                {
                    longint recver, sender, countRecv, countSend;
                    real *buffSend, *buffRecv;
                    if (REAL_TO_GHOSTS)
                    {
                        recver = nodeGrid.getNodeNeighborIndex(dir);
                        sender = nodeGrid.getNodeNeighborIndex(oppDir);
                        buffRecv = buffGhost.data();
                        buffSend = buffReal.data();
                        countRecv = commCellIdx[dir].numGhosts * 3;
                        countSend = commCellIdx[dir].numReals * 3;
                    }
                    else
                    {
                        recver = nodeGrid.getNodeNeighborIndex(oppDir);
                        sender = nodeGrid.getNodeNeighborIndex(dir);
                        buffRecv = buffReal.data();
                        buffSend = buffGhost.data();
                        countRecv = commCellIdx[dir].numReals * 3;
                        countSend = commCellIdx[dir].numGhosts * 3;
                    }
                    if (nodeGrid.getNodePosition(coord) % 2 == 0)
                    {
                        comm.send(recver, DD_COMM_TAG, buffSend, countSend);
                        comm.recv(sender, DD_COMM_TAG, buffRecv, countRecv);
                    }
                    else
                    {
                        comm.recv(sender, DD_COMM_TAG, buffRecv, countRecv);
                        comm.send(recver, DD_COMM_TAG, buffSend, countSend);
                    }
                }

                if (REAL_TO_GHOSTS)
                {
                    unpackCells<PACKED_POSITIONS, DATA_INSERT>(buffGhost, false, dir);
                }
                else
                {
                    unpackCells<PACKED_FORCES, DATA_ADD>(buffReal, true, dir);
                }
            }
        }
    }
}

template void DomainDecomposition::ghostCommunication_impl<false, true, 0>();

template void DomainDecomposition::ghostCommunication_impl<false, false, 0>();

///////////////////////////////////////////////////////////////////////////////////////////////
template <DomainDecomposition::AddShift DO_SHIFT>
void DomainDecomposition::copyRealsToGhostsIntra(size_t dir, Real3D const& shift)
{
    const auto& ccr = commCellIdx[dir].reals;
    const auto& ccg = commCellIdx[dir].ghosts;
    const size_t numCells = ccr.size();

    auto& particles = vectorization->particles;
    const auto& cellRange = particles.cellRange();

    {
        auto f_dim = [&](real* __restrict pos, real shift_v) {
            for (size_t ic = 0; ic < numCells; ic++)
            {
                const size_t icr = ccr[ic];
                const size_t icg = ccg[ic];
                const size_t numPart = cellRange[icr + 1] - cellRange[icr];

                const real* __restrict pos_r = pos + cellRange[icr];
                real* __restrict pos_g = pos + cellRange[icg];

                if (DO_SHIFT)
                {
                    ESPP_VEC_PRAGMAS
                    for (size_t ip = 0; ip < numPart; ip++)
                    {
                        pos_g[ip] = pos_r[ip] + shift_v;
                    }
                }
                else
                {
                    ESPP_VEC_PRAGMAS
                    for (size_t ip = 0; ip < numPart; ip++)
                    {
                        pos_g[ip] = pos_r[ip];
                    }
                }
            }
        };

        f_dim(particles.p_x.data(), shift[0]);
        f_dim(particles.p_y.data(), shift[1]);
        f_dim(particles.p_z.data(), shift[2]);
    }
}

template void DomainDecomposition::copyRealsToGhostsIntra<DomainDecomposition::ADD_SHIFT>(
    size_t dir, Real3D const& shift);

template void DomainDecomposition::copyRealsToGhostsIntra<DomainDecomposition::NO_SHIFT>(
    size_t dir, Real3D const& shift);

///////////////////////////////////////////////////////////////////////////////////////////////
void DomainDecomposition::addGhostForcesToRealsIntra(size_t dir)
{
    const auto& ccr = commCellIdx[dir].reals;
    const auto& ccg = commCellIdx[dir].ghosts;
    const size_t numCells = ccr.size();

    auto& particles = vectorization->particles;
    const auto& cellRange = particles.cellRange();

    {
        auto f_dim = [&](real* __restrict f) {
            for (size_t ic = 0; ic < numCells; ic++)
            {
                const size_t icr = ccr[ic];
                const size_t icg = ccg[ic];
                const size_t numPart = cellRange[icr + 1] - cellRange[icr];

                real* __restrict f_r = f + cellRange[icr];
                const real* __restrict f_g = f + cellRange[icg];

                {
                    ESPP_VEC_PRAGMAS
                    for (size_t ip = 0; ip < numPart; ip++)
                    {
                        f_r[ip] += f_g[ip];
                    }
                }
            }
        };
        f_dim(particles.f_x.data());
        f_dim(particles.f_y.data());
        f_dim(particles.f_z.data());
    }
}

///////////////////////////////////////////////////////////////////////////////////////////////
template <DomainDecomposition::PackedData PACKED_DATA, DomainDecomposition::AddShift DO_SHIFT>
void DomainDecomposition::packCells(AlignedVector<real>& sendBuf,
                                    bool commReal,
                                    size_t dir,
                                    Real3D const& shift)
{
    const auto& particles = vectorization->particles;
    const auto& cr = particles.cellRange();
    const auto& cc = commReal ? commCellIdx[dir].reals : commCellIdx[dir].ghosts;
    const auto& numPart = commReal ? commCellIdx[dir].numReals : commCellIdx[dir].numGhosts;

    {
        auto f_pack_dim = [&](size_t dim, const real* __restrict p_ptr, real shift_v) {
            real* __restrict b_ptr = sendBuf.data() + numPart * dim;

            size_t b_off = 0;
            for (const auto& ic : cc)
            {
                real* __restrict b_ptr_c = b_ptr + b_off;
                const real* __restrict p_ptr_c = p_ptr + cr[ic];
                const size_t npart = cr[ic + 1] - cr[ic];

                ESPP_VEC_PRAGMAS
                for (size_t ip = 0; ip < npart; ip++)
                {
                    if (DO_SHIFT == ADD_SHIFT)
                    {
                        b_ptr_c[ip] = p_ptr_c[ip] + shift_v;
                    }
                    else
                    {
                        b_ptr_c[ip] = p_ptr_c[ip];
                    }
                }
                b_off += npart;
            }
        };

        if (PACKED_DATA == PACKED_POSITIONS)
        {
            f_pack_dim(0, particles.p_x.data(), shift[0]);
            f_pack_dim(1, particles.p_y.data(), shift[1]);
            f_pack_dim(2, particles.p_z.data(), shift[2]);
        }
        else
        {
            f_pack_dim(0, particles.f_x.data(), shift[0]);
            f_pack_dim(1, particles.f_y.data(), shift[1]);
            f_pack_dim(2, particles.f_z.data(), shift[2]);
        }
    }
}

template void DomainDecomposition::packCells<DomainDecomposition::PACKED_POSITIONS,
                                             DomainDecomposition::ADD_SHIFT>(
    AlignedVector<real>& sendBuf, bool commReal, size_t dir, Real3D const& shift);

template void
DomainDecomposition::packCells<DomainDecomposition::PACKED_FORCES, DomainDecomposition::NO_SHIFT>(
    AlignedVector<real>& sendBuf, bool commReal, size_t dir, Real3D const& shift);

///////////////////////////////////////////////////////////////////////////////////////////////
template <DomainDecomposition::PackedData PACKED_DATA, DomainDecomposition::DataMode DATA_MODE>
void DomainDecomposition::unpackCells(AlignedVector<real> const& recvBuf, bool commReal, size_t dir)
{
    auto& particles = vectorization->particles;
    const auto& cr = particles.cellRange();
    const auto& cc = commReal ? commCellIdx[dir].reals : commCellIdx[dir].ghosts;
    const auto& numPart = commReal ? commCellIdx[dir].numReals : commCellIdx[dir].numGhosts;

    {
        auto f_unpack_dim = [&](size_t dim, real* __restrict p_ptr) {
            const size_t b_start = (numPart * dim);
            const real* __restrict b_ptr = recvBuf.data() + b_start;

            /// loop over cells
            size_t b_off = 0;
            for (const auto& ic : cc)
            {
                const real* __restrict b_ptr_c = b_ptr + b_off;
                real* __restrict p_ptr_c = p_ptr + cr[ic];
                const size_t npart = cr[ic + 1] - cr[ic];

                ESPP_VEC_PRAGMAS
                for (size_t ip = 0; ip < npart; ip++)
                {
                    if (DATA_MODE == DATA_ADD)
                    {
                        p_ptr_c[ip] += b_ptr_c[ip];
                    }
                    else
                    {
                        p_ptr_c[ip] = b_ptr_c[ip];
                    }
                }
                b_off += npart;
            }
        };

        if (PACKED_DATA == PACKED_POSITIONS)
        {
            f_unpack_dim(0, particles.p_x.data());
            f_unpack_dim(1, particles.p_y.data());
            f_unpack_dim(2, particles.p_z.data());
        }
        else
        {
            f_unpack_dim(0, particles.f_x.data());
            f_unpack_dim(1, particles.f_y.data());
            f_unpack_dim(2, particles.f_z.data());
        }
    }
}

template void DomainDecomposition::unpackCells<DomainDecomposition::PACKED_POSITIONS,
                                               DomainDecomposition::DATA_INSERT>(
    AlignedVector<real> const& recvBuf, bool commReal, size_t dir);

template void
DomainDecomposition::unpackCells<DomainDecomposition::PACKED_FORCES, DomainDecomposition::DATA_ADD>(
    AlignedVector<real> const& recvBuf, bool commReal, size_t dir);

///////////////////////////////////////////////////////////////////////////////////////////////
void DomainDecomposition::registerPython()
{
    using namespace espressopp::python;

    class_<DomainDecomposition, bases<espressopp::storage::DomainDecomposition, StorageVec>,
           boost::noncopyable>("vec_storage_DomainDecomposition",
                               init<std::shared_ptr<System>, const Int3D&, const Int3D&, int>())
        // .def("initChannels", &DomainDecomposition::initChannels)
        // .def("getChannelIndices", &DomainDecomposition::getChannelIndices)
        // .def("connectOffload", &DomainDecomposition::connectOffload)
        // .def("connectedOffload", &DomainDecomposition::connectedOffload)
        // .def("disconnectOffload", &DomainDecomposition::disconnectOffload)
        .def("resetCells", &DomainDecomposition::resetCells)
        // .def("resetTimers", &DomainDecomposition::resetTimers)
        // .def("getTimers", &wrapGetTimers)
        // .def("getTimers2", &wrapGetTimers2)
        ;
}

}  // namespace storage
}  // namespace vec
}  // namespace espressopp

/// WARNING: Storage::localParticles is unused in this implementation
