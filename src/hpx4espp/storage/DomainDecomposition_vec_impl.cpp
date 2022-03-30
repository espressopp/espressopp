/*
  Copyright (C) 2020-2022
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

#include <hpx/config.hpp>
#include <hpx/future.hpp>

#include "hpx4espp/storage/DomainDecomposition.hpp"
#include "hpx4espp/utils/algorithms/for_loop.hpp"
#include "hpx4espp/utils/assert.hpp"
#include "hpx4espp/utils/assert_msg.hpp"
#include "bc/BC.hpp"

const int DD_COMM_TAG = 0xab;

/// manually enable by adding { STATEMENT }
#define HPX4ESPP_DDVI_DEBUG(STATEMENT)
#define HPX4ESPP_DDVI_CHANNELS_DEBUG(STATEMENT)

namespace espressopp
{
namespace hpx4espp
{
namespace storage
{
template <bool SIZES_FIRST, bool REAL_TO_GHOSTS, int EXTRA_DATA>
void DomainDecomposition::ghostCommunication_impl()
{
    HPX4ESPP_DDVI_DEBUG(HPX4ESPP_ASSERT_EQUAL(SIZES_FIRST, false);
                        HPX4ESPP_ASSERT_EQUAL(EXTRA_DATA, 0);)

    auto const& comm = *(getSystem()->comm);

    for (size_t _coord = 0; _coord < 3; ++_coord)
    {
        int coord = REAL_TO_GHOSTS ? _coord : (2 - _coord);
        const real curCoordBoxL = getSystem()->bc->getBoxL()[coord];
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

            auto f_interNode_pack = [this, &shift, dir]() {
                const real time = wallTimer.getElapsedTime();
                auto packNode = [this, dir, &shift](size_t idxCommNode) {
                    if (REAL_TO_GHOSTS)
                    {
                        packCells<PACKED_POSITIONS, ADD_SHIFT>(buffReal, true, dir, idxCommNode,
                                                               shift);
                    }
                    else
                    {
                        packCells<PACKED_FORCES, NO_SHIFT>(buffGhost, false, dir, idxCommNode,
                                                           SHIFT_ZERO);
                    }
                };
                size_t const numCommNodes =
                    REAL_TO_GHOSTS ? commNodesReal[dir].size() : commNodesGhost[dir].size();
                utils::parallelForLoop(0, numCommNodes, packNode);

                timeUpdateGhosts_InterNode_pack += wallTimer.getElapsedTime() - time;
            };

            auto f_interNode_comm = [this, &comm, coord, dir, oppDir]() {
                const real time = wallTimer.getElapsedTime();
                longint recver, sender, countRecv, countSend;
                real *buffSend, *buffRecv;
                if (REAL_TO_GHOSTS)
                {
                    recver = nodeGrid.getNodeNeighborIndex(dir);
                    sender = nodeGrid.getNodeNeighborIndex(oppDir);
                    buffRecv = buffGhost.data();
                    buffSend = buffReal.data();
                    countRecv = nodeRangeGhost[dir].back() * vecModeFactor;
                    countSend = nodeRangeReal[dir].back() * vecModeFactor;
                }
                else
                {
                    recver = nodeGrid.getNodeNeighborIndex(oppDir);
                    sender = nodeGrid.getNodeNeighborIndex(dir);
                    buffRecv = buffReal.data();
                    buffSend = buffGhost.data();
                    countRecv = nodeRangeReal[dir].back() * vecModeFactor;
                    countSend = nodeRangeGhost[dir].back() * vecModeFactor;
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
                timeUpdateGhosts_InterNode_comm += wallTimer.getElapsedTime() - time;
            };

            auto f_interNode_unpack = [this, dir]() {
                const real time = wallTimer.getElapsedTime();

                auto unpackNode = [this, dir](size_t idxCommNode) {
                    if (REAL_TO_GHOSTS)
                    {
                        unpackCells<PACKED_POSITIONS, DATA_INSERT>(buffGhost, false, dir,
                                                                   idxCommNode);
                    }
                    else
                    {
                        unpackCells<PACKED_FORCES, DATA_ADD>(buffReal, true, dir, idxCommNode);
                    }
                };
                size_t const numCommNodes =
                    REAL_TO_GHOSTS ? commNodesGhost[dir].size() : commNodesReal[dir].size();
                utils::parallelForLoop(0, numCommNodes, unpackNode);

                timeUpdateGhosts_InterNode_unpack += wallTimer.getElapsedTime() - time;
            };

            /// Inter-node block
            auto f_interNode = [this, &shift, &f_interNode_pack, &f_interNode_comm,
                                &f_interNode_unpack, coord, dir, oppDir]() {
                f_interNode_pack();
                f_interNode_comm();
                f_interNode_unpack();
            };

            /// Intra-node block (including shifted/periodic)
            auto f_intraNode = [this, &shift, doPeriodic, dir]() {
                const real time = wallTimer.getElapsedTime();

                if (HPX4ESPP_OPTIMIZE_COMM)
                {
                    auto f_pair = [this, dir, &shift](size_t i) {
                        const size_t ip = i / numCommSubs;
                        const size_t is = i % numCommSubs;

                        const auto& pair = subNodePairsIntraPeriodic[dir][ip];

                        {
                            if (REAL_TO_GHOSTS)
                            {
                                copyRealsToGhostsIntra<ADD_SHIFT>(dir, std::get<0>(pair),
                                                                  std::get<1>(pair), is, shift);
                            }
                            else
                            {
                                addGhostForcesToRealsIntra(dir, std::get<0>(pair),
                                                           std::get<1>(pair), is);
                            }
                        }
                        return;
                    };
                    size_t const numPairs = subNodePairsIntraPeriodic[dir].size();
                    utils::parallelForLoop(0, numPairs * numCommSubs, f_pair);
                }
                else
                {
                    auto f_pair = [this, dir, &shift](size_t i) {
                        const size_t ip = i / numCommSubs;
                        const size_t is = i % numCommSubs;

                        const auto& pair = subNodePairsIntra[dir][ip];

                        if (REAL_TO_GHOSTS)
                        {
                            if (std::get<2>(pair))
                                copyRealsToGhostsIntra<ADD_SHIFT>(dir, std::get<0>(pair),
                                                                  std::get<1>(pair), is, shift);
                            else
                            {
                                copyRealsToGhostsIntra<NO_SHIFT>(dir, std::get<0>(pair),
                                                                 std::get<1>(pair), is, SHIFT_ZERO);
                            }
                        }
                        else
                        {
                            addGhostForcesToRealsIntra(dir, std::get<0>(pair), std::get<1>(pair),
                                                       is);
                        }
                    };
                    size_t const numPairs = subNodePairsIntra[dir].size();
                    utils::parallelForLoop(0, numPairs * numCommSubs, f_pair);
                }

                /// NOTE: Time only if comm is involved
                // if (!doPeriodic)
                timeUpdateGhosts_IntraNode += wallTimer.getElapsedTime() - time;
            };

            if ((hpx::threads::get_self_ptr() != nullptr) && !doPeriodic && commAsync)
            {
                {
                    f_interNode_pack();
                    std::array<hpx::future<void>, 2> ret = {hpx::async(f_interNode_comm),
                                                            hpx::async(f_intraNode)};
                    hpx::wait_all(ret);
                    f_interNode_unpack();
                }
            }
            else
            {
                if (!doPeriodic) f_interNode();
                f_intraNode();
            }
        }
    }
}

template void DomainDecomposition::ghostCommunication_impl<false, true, 0>();
template void DomainDecomposition::ghostCommunication_impl<false, false, 0>();

template <DomainDecomposition::AddShift DO_SHIFT>
void DomainDecomposition::copyRealsToGhostsIntra(
    size_t dir, size_t ir, size_t ig, size_t is, Real3D const& shift)
{
    if (HPX4ESPP_OPTIMIZE_COMM)
    {
        const auto& ccr = vdd[ir].commCellsOwned[dir].reals;
        const auto& ccg = vdd[ig].commCellsOwned[dir].ghosts;
        HPX4ESPP_DDVI_DEBUG(HPX4ESPP_ASSERT_EQUAL(ccr.size(), ccg.size());)
        const size_t numCells = ccr.size();

        const auto& pr = virtualStorage[ir].particles;
        auto& pg = virtualStorage[ig].particles;
        const auto& crr = pr.cellRange();
        const auto& crg = pg.cellRange();

        {
            auto f_dim = [&](const real* __restrict ptr_r, real* __restrict ptr_g, real shift_v) {
                for (size_t ic = 0; ic < numCells; ic++)
                {
                    const size_t icr = ccr[ic];
                    const size_t icg = ccg[ic];

                    const size_t numPart = crr[icr + 1] - crr[icr];
                    HPX4ESPP_DDVI_DEBUG(HPX4ESPP_ASSERT_EQUAL(numPart, (crg[icg + 1] - crg[icg]));)

                    const real* __restrict ptr_s_r = ptr_r + crr[icr];
                    real* __restrict ptr_s_g = ptr_g + crg[icg];

                    if (DO_SHIFT)
                    {
                        ESPP_VEC_PRAGMAS
                        for (size_t ip = 0; ip < numPart; ip++)
                        {
                            ptr_s_g[ip] = ptr_s_r[ip] + shift_v;
                        }
                    }
                    else
                    {
                        ESPP_VEC_PRAGMAS
                        for (size_t ip = 0; ip < numPart; ip++)
                        {
                            ptr_s_g[ip] = ptr_s_r[ip];
                        }
                    }
                }
            };
            f_dim(pr.p_x.data(), pg.p_x.data(), shift[0]);
            f_dim(pr.p_y.data(), pg.p_y.data(), shift[1]);
            f_dim(pr.p_z.data(), pg.p_z.data(), shift[2]);
        }
    }
    else
    {
        // if(is!=1) return;

        const auto& ccr = vdd[ir].commCells[dir].reals;
        const auto& ccg = vdd[ig].commCells[dir].ghosts;
        HPX4ESPP_DDVI_DEBUG(HPX4ESPP_ASSERT_EQUAL(ccr.size(), ccg.size());)
        // const size_t numCells = ccr.size();
        const auto& c_range = vdd[ir].commSubCellRanges[dir][is];

        const auto& pr = virtualStorage[ir].particles;
        auto& pg = virtualStorage[ig].particles;
        const auto& crr = pr.cellRange();
        const auto& crg = pg.cellRange();

        {
            auto f_dim = [&](const real* __restrict ptr_r, real* __restrict ptr_g, real shift_v) {
                // for (size_t ic = 0; ic < numCells; ic++)
                for (size_t ic = c_range.first; ic < c_range.second; ic++)
                {
                    const size_t icr = ccr[ic];
                    const size_t icg = ccg[ic];

                    const size_t numPart = crr[icr + 1] - crr[icr];
                    HPX4ESPP_DDVI_DEBUG(HPX4ESPP_ASSERT_EQUAL(numPart, (crg[icg + 1] - crg[icg]));)

                    const real* __restrict ptr_s_r = ptr_r + crr[icr];
                    real* __restrict ptr_s_g = ptr_g + crg[icg];

                    if (DO_SHIFT)
                    {
                        ESPP_VEC_PRAGMAS
                        for (size_t ip = 0; ip < numPart; ip++)
                        {
                            ptr_s_g[ip] = ptr_s_r[ip] + shift_v;
                        }
                    }
                    else
                    {
                        ESPP_VEC_PRAGMAS
                        for (size_t ip = 0; ip < numPart; ip++)
                        {
                            ptr_s_g[ip] = ptr_s_r[ip];
                        }
                    }
                }
            };
            f_dim(pr.p_x.data(), pg.p_x.data(), shift[0]);
            f_dim(pr.p_y.data(), pg.p_y.data(), shift[1]);
            f_dim(pr.p_z.data(), pg.p_z.data(), shift[2]);
        }
    }
}

template void DomainDecomposition::copyRealsToGhostsIntra<DomainDecomposition::ADD_SHIFT>(
    size_t dir, size_t ir, size_t ig, size_t is, Real3D const& shift);

template void DomainDecomposition::copyRealsToGhostsIntra<DomainDecomposition::NO_SHIFT>(
    size_t dir, size_t ir, size_t ig, size_t is, Real3D const& shift);

void DomainDecomposition::addGhostForcesToRealsIntra(size_t dir, size_t ir, size_t ig, size_t is)
{
    if (HPX4ESPP_OPTIMIZE_COMM)
    {
        const auto& ccr = vdd[ir].commCellsOwned[dir].reals;
        const auto& ccg = vdd[ig].commCellsOwned[dir].ghosts;

        HPX4ESPP_DDVI_DEBUG(HPX4ESPP_ASSERT_EQUAL(ccr.size(), ccg.size());)

        const size_t numCells = ccr.size();

        auto& pr = virtualStorage[ir].particles;
        const auto& pg = virtualStorage[ig].particles;
        const auto& crr = pr.cellRange();
        const auto& crg = pg.cellRange();

        {
            auto f_dim = [&](real* __restrict ptr_r, const real* __restrict ptr_g) {
                for (size_t ic = 0; ic < numCells; ic++)
                {
                    const size_t icr = ccr[ic];
                    const size_t icg = ccg[ic];

                    const size_t numPart = crr[icr + 1] - crr[icr];

                    HPX4ESPP_DDVI_DEBUG(HPX4ESPP_ASSERT_EQUAL(numPart, (crg[icg + 1] - crg[icg]));)

                    real* __restrict ptr_s_r = ptr_r + crr[icr];
                    const real* __restrict ptr_s_g = ptr_g + crg[icg];

                    {
                        ESPP_VEC_PRAGMAS
                        for (size_t ip = 0; ip < numPart; ip++)
                        {
                            ptr_s_r[ip] += ptr_s_g[ip];
                        }
                    }
                }
            };
            f_dim(pr.f_x.data(), pg.f_x.data());
            f_dim(pr.f_y.data(), pg.f_y.data());
            f_dim(pr.f_z.data(), pg.f_z.data());
        }
    }
    else
    {
        // if(is!=1) return;

        const auto& ccr = vdd[ir].commCells[dir].reals;
        const auto& ccg = vdd[ig].commCells[dir].ghosts;

        HPX4ESPP_DDVI_DEBUG(HPX4ESPP_ASSERT_EQUAL(ccr.size(), ccg.size());)

        // const size_t numCells = ccr.size();
        const auto& c_range = vdd[ir].commSubCellRanges[dir][is];

        auto& pr = virtualStorage[ir].particles;
        const auto& pg = virtualStorage[ig].particles;
        const auto& crr = pr.cellRange();
        const auto& crg = pg.cellRange();

        {
            auto f_dim = [&](real* __restrict ptr_r, const real* __restrict ptr_g) {
                // for (size_t ic = 0; ic < numCells; ic++)
                for (size_t ic = c_range.first; ic < c_range.second; ic++)
                {
                    const size_t icr = ccr[ic];
                    const size_t icg = ccg[ic];

                    const size_t numPart = crr[icr + 1] - crr[icr];

                    HPX4ESPP_DDVI_DEBUG(HPX4ESPP_ASSERT_EQUAL(numPart, (crg[icg + 1] - crg[icg]));)

                    real* __restrict ptr_s_r = ptr_r + crr[icr];
                    const real* __restrict ptr_s_g = ptr_g + crg[icg];

                    {
                        ESPP_VEC_PRAGMAS
                        for (size_t ip = 0; ip < numPart; ip++)
                        {
                            ptr_s_r[ip] += ptr_s_g[ip];
                        }
                    }
                }
            };
            f_dim(pr.f_x.data(), pg.f_x.data());
            f_dim(pr.f_y.data(), pg.f_y.data());
            f_dim(pr.f_z.data(), pg.f_z.data());
        }
    }
}

template <DomainDecomposition::PackedData PACKED_DATA, DomainDecomposition::AddShift DO_SHIFT>
void DomainDecomposition::packCells(vec::AlignedVector<real>& sendBuf,
                                    bool commReal,
                                    size_t dir,
                                    size_t idxCommNode,
                                    Real3D const& shift)
{
    const size_t nodeStart =
        commReal ? nodeRangeReal[dir][idxCommNode] : nodeRangeGhost[dir][idxCommNode];
    const size_t nodeEnd =
        commReal ? nodeRangeReal[dir][idxCommNode + 1] : nodeRangeGhost[dir][idxCommNode + 1];
    const size_t numPart = nodeEnd - nodeStart;
    const size_t inode =
        commReal ? commNodesReal[dir][idxCommNode] : commNodesGhost[dir][idxCommNode];
    const auto& vs = virtualStorage[inode];
    const auto& cr = vs.particles.cellRange();
    const auto& vd = vdd[inode];
    const auto& cc = commReal ? vd.commCells[dir].reals : vd.commCells[dir].ghosts;

    {
        /// loop over dimensions (x,y,z)
        auto f_pack_dim = [dir, nodeStart, nodeEnd, numPart, inode, &cc, &cr, &sendBuf, &vd, this](
                              size_t dim, const real* __restrict p_ptr, real shift_v) {
            const size_t b_start = (nodeStart * 3) + (numPart * dim);
            real* __restrict b_ptr = sendBuf.data() + b_start;

            /// loop over cells
            size_t b_off = 0;
            for (const auto& ic : cc)
            {
                real* __restrict b_ptr_c = b_ptr + b_off;
                const real* __restrict p_ptr_c = p_ptr + cr[ic];
                const size_t npart = cr[ic + 1] - cr[ic];
                b_off += npart;

                if (HPX4ESPP_OPTIMIZE_COMM)
                {
                    const bool icOwn = (vd.cellGridInfo[ic].subNode == inode);
                    if (!icOwn) continue;
                }

                /// loop over particles
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
            }
        };

        if (PACKED_DATA == PACKED_POSITIONS)
        {
            f_pack_dim(0, vs.particles.p_x.data(), shift[0]);
            f_pack_dim(1, vs.particles.p_y.data(), shift[1]);
            f_pack_dim(2, vs.particles.p_z.data(), shift[2]);
        }
        else
        {
            f_pack_dim(0, vs.particles.f_x.data(), shift[0]);
            f_pack_dim(1, vs.particles.f_y.data(), shift[1]);
            f_pack_dim(2, vs.particles.f_z.data(), shift[2]);
        }
    }
}

template void DomainDecomposition::packCells<DomainDecomposition::PACKED_POSITIONS,
                                             DomainDecomposition::ADD_SHIFT>(
    vec::AlignedVector<real>& sendBuf,
    bool commReal,
    size_t dir,
    size_t idxCommNode,
    Real3D const& shift);

template void
DomainDecomposition::packCells<DomainDecomposition::PACKED_FORCES, DomainDecomposition::NO_SHIFT>(
    vec::AlignedVector<real>& sendBuf,
    bool commReal,
    size_t dir,
    size_t idxCommNode,
    Real3D const& shift);

template <DomainDecomposition::PackedData PACKED_DATA, DomainDecomposition::DataMode DATA_MODE>
void DomainDecomposition::unpackCells(vec::AlignedVector<real> const& recvBuf,
                                      bool commReal,
                                      size_t dir,
                                      size_t idxCommNode)
{
    const size_t nodeStart =
        commReal ? nodeRangeReal[dir][idxCommNode] : nodeRangeGhost[dir][idxCommNode];
    const size_t nodeEnd =
        commReal ? nodeRangeReal[dir][idxCommNode + 1] : nodeRangeGhost[dir][idxCommNode + 1];
    const size_t numPart = nodeEnd - nodeStart;
    const size_t inode =
        commReal ? commNodesReal[dir][idxCommNode] : commNodesGhost[dir][idxCommNode];
    auto& vs = virtualStorage[inode];
    const auto& cr = vs.particles.cellRange();
    const auto& vd = vdd[inode];
    const auto& cc = commReal ? vd.commCells[dir].reals : vd.commCells[dir].ghosts;

    {
        /// loop over dimensions (x,y,z)
        auto f_pack_dim = [dir, nodeStart, nodeEnd, numPart, inode, &cc, &cr, &recvBuf, &vd, this](
                              size_t dim, real* __restrict p_ptr) {
            const size_t b_start = (nodeStart * 3) + (numPart * dim);
            const real* __restrict b_ptr = recvBuf.data() + b_start;

            /// loop over cells
            size_t b_off = 0;
            for (const auto& ic : cc)
            {
                const real* __restrict b_ptr_c = b_ptr + b_off;
                real* __restrict p_ptr_c = p_ptr + cr[ic];
                const size_t npart = cr[ic + 1] - cr[ic];
                b_off += npart;

                if (HPX4ESPP_OPTIMIZE_COMM)
                {
                    const bool icOwn = (vd.cellGridInfo[ic].subNode == inode);
                    if (!icOwn) continue;
                }

                /// loop over particles
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
            }
        };

        if (PACKED_DATA == PACKED_POSITIONS)
        {
            f_pack_dim(0, vs.particles.p_x.data());
            f_pack_dim(1, vs.particles.p_y.data());
            f_pack_dim(2, vs.particles.p_z.data());
        }
        else
        {
            f_pack_dim(0, vs.particles.f_x.data());
            f_pack_dim(1, vs.particles.f_y.data());
            f_pack_dim(2, vs.particles.f_z.data());
        }
    }
}

template void DomainDecomposition::unpackCells<DomainDecomposition::PACKED_POSITIONS,
                                               DomainDecomposition::DATA_INSERT>(
    vec::AlignedVector<real> const& recvBuf, bool commReal, size_t dir, size_t idxCommNode);

template void
DomainDecomposition::unpackCells<DomainDecomposition::PACKED_FORCES, DomainDecomposition::DATA_ADD>(
    vec::AlignedVector<real> const& recvBuf, bool commReal, size_t dir, size_t idxCommNode);

template <bool ALIGNED>
void DomainDecomposition::exchangeGhosts_impl()
{
    auto const& comm = *(getSystem()->comm);
    const int rank = comm.rank();

    constexpr bool sizesFirst = true;
    constexpr bool realToGhosts = true;
    constexpr int extradata = DATA_PROPERTIES;

    constexpr size_t sizePos = sizeof(ParticlePosition);
    constexpr size_t sizePrp = sizeof(ParticleProperties);
    constexpr size_t sizeTotal = sizePos + sizePrp;

    // LOG4ESPP_DEBUG(logger, "do ghost communication " << (sizesFirst ? "with sizes " : "")
    //       << (realToGhosts ? "reals to ghosts " : "ghosts to reals ") << extradata);

    /* direction loop: x, y, z.
      Here we could in principle build in a one sided ghost
      communication, simply by taking the lr loop only over one
      value. */
    for (int _coord = 0; _coord < 3; ++_coord)
    {
        /* inverted processing order for ghost force communication,
          since the corner ghosts have to be collected via several
          nodes. We now add back the corner ghost forces first again
          to ghost forces, which only eventually go back to the real
          particle.
        */
        int coord = realToGhosts ? _coord : (2 - _coord);
        real curCoordBoxL = getSystem()->bc->getBoxL()[coord];

        // lr loop: left right
        for (int lr = 0; lr < 2; ++lr)
        {
            int dir = 2 * coord + lr;
            int oppositeDir = 2 * coord + (1 - lr);

            Real3D shift(0, 0, 0);

            shift[coord] = nodeGrid.getBoundary(dir) * curCoordBoxL;

            // LOG4ESPP_DEBUG(logger, "direction " << dir);

            if (nodeGrid.getGridSize(coord) == 1)
            {
                const real time = wallTimer.getElapsedTime();

                // LOG4ESPP_DEBUG(logger, "local communication");

                // copy operation, we have to receive as many cells as we send
                if (commCells[dir].ghosts.size() != commCells[dir].reals.size())
                {
                    throw std::runtime_error(
                        "DomainDecomposition::doGhostCommunication: send/recv cell structure "
                        "mismatch during local copy");
                }

                auto f_copyRealsToGhosts = [this, &shift, &extradata, &dir](size_t i) {
                    copyRealsToGhosts(*commCells[dir].reals[i], *commCells[dir].ghosts[i],
                                      extradata, shift);
                };
                size_t const numCommCells = commCells[dir].ghosts.size();
                utils::parallelForLoop(0, numCommCells, f_copyRealsToGhosts);

                timeExcGhosts_IntraNode += wallTimer.getElapsedTime() - time;
            }
            else
            {
                // exchange size information, if necessary
                std::vector<longint> sendSizes, recvSizes, sendRanges, recvRanges;
                longint sendTotalBytes, recvTotalBytes;
                {
                    const real time = wallTimer.getElapsedTime();

                    // prepare buffers
                    sendSizes.reserve(commCells[dir].reals.size());
                    for (int i = 0, end = commCells[dir].reals.size(); i < end; ++i)
                    {
                        sendSizes.push_back(commCells[dir].reals[i]->particles.size());
                    }
                    recvSizes.resize(commCells[dir].ghosts.size());

                    // exchange sizes, odd-even rule
                    if (nodeGrid.getNodePosition(coord) % 2 == 0)
                    {
                        comm.send(nodeGrid.getNodeNeighborIndex(dir), DD_COMM_TAG, &(sendSizes[0]),
                                  sendSizes.size());
                        comm.recv(nodeGrid.getNodeNeighborIndex(oppositeDir), DD_COMM_TAG,
                                  &(recvSizes[0]), recvSizes.size());
                    }
                    else
                    {
                        comm.recv(nodeGrid.getNodeNeighborIndex(oppositeDir), DD_COMM_TAG,
                                  &(recvSizes[0]), recvSizes.size());
                        comm.send(nodeGrid.getNodeNeighborIndex(dir), DD_COMM_TAG, &(sendSizes[0]),
                                  sendSizes.size());
                    }

                    // resize according to received information
                    auto f_resize = [this, &recvSizes, &dir](size_t i) {
                        commCells[dir].ghosts[i]->particles.resize(recvSizes[i]);
                    };
                    size_t const numCommCells = commCells[dir].ghosts.size();
                    utils::parallelForLoop(0, numCommCells, f_resize);

                    auto f_sizesToByteRanges = [&sizeTotal](std::vector<longint> const& sizes) {
                        const int nsize = sizes.size();
                        std::vector<longint> ranges;
                        ranges.resize(nsize + 1);
                        longint total = 0;
                        for (int i = 0; i < nsize; i++)
                        {
                            ranges[i] = total;
                            if (ALIGNED)
                            {
                                total += ROUND_TO_CACHE_LINE(sizeTotal * sizes[i]);
                            }
                            else
                            {
                                total += sizeTotal * sizes[i];
                            }
                        }
                        ranges[nsize] = total;
                        return ranges;
                    };
                    sendRanges = f_sizesToByteRanges(sendSizes);
                    recvRanges = f_sizesToByteRanges(recvSizes);
                    sendTotalBytes = sendRanges.back();
                    recvTotalBytes = recvRanges.back();

                    timeExcGhosts_InterNode_sizes += wallTimer.getElapsedTime() - time;
                }

                // prepare send and receive buffers
                longint receiver, sender;
                {
                    const real time = wallTimer.getElapsedTime();

                    auto f_fillBufView = [](auto& bufview, const auto& buffixed,
                                            const auto& ranges) {
                        const size_t size = ranges.size() - 1;
                        if (bufview.size() < size) bufview.resize(size);
                        for (size_t i = 0; i < size; i++)
                        {
                            bufview[i].view(buffixed, ranges[i], ranges[i + 1]);
                        }
                    };

                    outBuf.allocate(sendTotalBytes);
                    f_fillBufView(outBufView, outBuf, sendRanges);
                    inBuf.allocate(recvTotalBytes);
                    f_fillBufView(inBufView, inBuf, recvRanges);

                    {
                        receiver = nodeGrid.getNodeNeighborIndex(dir);
                        sender = nodeGrid.getNodeNeighborIndex(oppositeDir);
                        auto f_pack = [](OutBufferView& buf, Cell& _reals, int extradata,
                                         const Real3D& shift) {
                            ParticleList& reals = _reals.particles;
                            for (ParticleList::iterator src = reals.begin(), end = reals.end();
                                 src != end; ++src)
                            {
                                buf.write(*src, extradata, shift);
                            }
                        };
                        const size_t size = commCells[dir].reals.size();
                        utils::parallelForLoop(
                            0, size, [this, &f_pack, &extradata, &shift, dir](size_t i) {
                                f_pack(outBufView[i], *commCells[dir].reals[i], extradata, shift);
                            });
                    }

                    timeExcGhosts_InterNode_pack += wallTimer.getElapsedTime() - time;
                }

                {
                    const real time = wallTimer.getElapsedTime();

                    // exchange particles, odd-even rule
                    if (nodeGrid.getNodePosition(coord) % 2 == 0)
                    {
                        outBuf.send(receiver, DD_COMM_TAG);
                        inBuf.recv(sender, DD_COMM_TAG);
                    }
                    else
                    {
                        inBuf.recv(sender, DD_COMM_TAG);
                        outBuf.send(receiver, DD_COMM_TAG);
                    }

                    timeExcGhosts_InterNode_comm += wallTimer.getElapsedTime() - time;
                }

                {
                    const real time = wallTimer.getElapsedTime();

                    {
                        auto f_unpack = [](Cell& _ghosts, InBufferView& buf, int extradata) {
                            ParticleList& ghosts = _ghosts.particles;
                            for (ParticleList::iterator dst = ghosts.begin(), end = ghosts.end();
                                 dst != end; ++dst)
                            {
                                buf.read(*dst, extradata);
                                /// NOTE: multithreaded particle map not implemented
                                // if (extradata & DATA_PROPERTIES) {
                                //   updateInLocalParticles(&(*dst), true);
                                // }
                                dst->ghost() = 1;
                            }
                        };
                        const size_t size = commCells[dir].ghosts.size();
                        utils::parallelForLoop(
                            0, size, [this, &f_unpack, &extradata, dir](size_t i) {
                                f_unpack(*commCells[dir].ghosts[i], inBufView[i], extradata);
                            });
                    }

                    timeExcGhosts_InterNode_unpack += wallTimer.getElapsedTime() - time;
                }
            }
        }
    }
    LOG4ESPP_DEBUG(logger, "ghost communication finished");
}

template void DomainDecomposition::exchangeGhosts_impl<true>();
template void DomainDecomposition::exchangeGhosts_impl<false>();

///////////////////////////////////////////////////////////////////////////////////////////////
/// Channels

void DomainDecomposition::prepareChannelAttributes()
{
    /// NOTE: call function during channel initialization and (TODO:) cellAdjust
    const int rank = mpiWorld->rank();

    // std::vector<std::tuple<int,int>> nsubComm(6*numSubNodes,{0,0});
    nsubComm = std::vector<std::tuple<int, int>>(6 * numSubNodes, {0, 0});

    /// use data from commNodesReal/Ghost
    if (!rank) std::cout << "prepareChannelAttributes" << std::endl;
    for (int coord = 0; coord < 3; ++coord)
    {
        const real curCoordBoxL = getSystem()->bc->getBoxL()[coord];
        const bool doPeriodic = (nodeGrid.getGridSize(coord) == 1);
        for (int lr = 0; lr < 2; ++lr)
        {
            int const dir = 2 * coord + lr;
            const auto& thisCommNodes = commNodesReal[dir];

            /// store senders for inter-node (except periodic part)
            if (!doPeriodic)
            {
                for (int ic = 0; ic < thisCommNodes.size(); ic++)
                {
                    const auto ii = thisCommNodes[ic];
                    nsubComm[ii + numSubNodes * dir] = {1, ic};
                }
            }

            HPX4ESPP_DDVI_CHANNELS_DEBUG({
                if (!rank)
                    std::cout << "   coord: " << coord << " lr: " << lr << " commNodesReal[" << dir
                              << "]:";
                if (!rank)
                    for (auto& c : commNodesReal[dir]) std::cout << " " << c;
                if (!rank) std::cout << " size: " << commNodesReal.size() << std::endl;
                if (!rank)
                    std::cout << "   coord: " << coord << " lr: " << lr << " commNodesGhost[" << dir
                              << "]:";
                if (!rank)
                    for (auto& c : commNodesGhost[dir]) std::cout << " " << c;
                if (!rank) std::cout << " size: " << commNodesGhost.size() << std::endl;
            })

            /// store recver (self->dst) nsubs for local comm in remaining nsubComm
            /// periodic case already handled by NodeGrid::getNodeNeighborIndex
            for (size_t i2 = 0; i2 < subNodeGrid[2]; i2++)
            {
                for (size_t i1 = 0; i1 < subNodeGrid[1]; i1++)
                {
                    for (size_t i0 = 0; i0 < subNodeGrid[0]; i0++)
                    {
                        const size_t ii = i0 + subNodeGrid[0] * (i1 + subNodeGrid[1] * i2);
                        auto& ns = nsubComm[ii + numSubNodes * dir];
                        if (std::get<0>(ns)) continue;
                        std::get<1>(ns) = vdd[ii].nodeGridLocal.getNodeNeighborIndex(dir);
                    }
                }
            }

            HPX4ESPP_DDVI_CHANNELS_DEBUG({
                /// sample indexing
                if (!rank)
                    std::cout << "   coord: " << coord << " lr: " << lr << " nsubComm[" << dir
                              << "]:";
                if (!rank)
                    for (int ii = 0; ii < numSubNodes; ii++)
                    {
                        auto const& ns = nsubComm[ii + numSubNodes * dir];
                        std::cout << " " << ii << ":(" << std::get<0>(ns) << "," << std::get<1>(ns)
                                  << ")";

                        if (std::get<0>(ns))
                        {
                            auto const& idx = channels.sendIdx(coord, lr, std::get<1>(ns));
                            std::cout << "[" << std::get<0>(idx) << "," << std::get<1>(idx) << ","
                                      << std::get<2>(idx) << "," << std::get<3>(idx) << "]";
                        }
                    }
                if (!rank) std::cout << " size: " << numSubNodes << std::endl;
                if (!rank) std::cout << std::endl;
            })
        }
    }

    HPX4ESPP_DDVI_CHANNELS_DEBUG(if (!rank) {
        std::cout << "vdd[inode]->nodeGridLocal->getBoundary[dir]" << std::endl;
        for (size_t coord = 0; coord < 3; ++coord)
        {
            for (size_t lr = 0; lr < 2; ++lr)
            {
                size_t const dir = 2 * coord + lr;

                // std::cout << "   coord: " << coord << " lr: " << lr << " dir: "<<dir<<"
                // (inode,val,ns[0]):";
                std::cout << "   coord: " << coord << " lr: " << lr << " dir: " << dir
                          << " (inode,val):";

                for (size_t inode = 0; inode < vdd.size(); inode++)
                {
                    const auto& vd = vdd[inode];
                    const auto& ns =
                        nsubComm[inode + numSubNodes * dir /* (REAL_TO_GHOSTS ? dir : oppDir) */];
                    // std::cout << " (" << inode << "," << vd.nodeGridLocal.getBoundary(dir) << ","
                    // << std::get<0>(ns) << ")"; std::cout << " (" << inode << "," <<
                    // vd.nodeGridLocal.getBoundary(dir) << ")";
                    if (vd.nodeGridLocal.getBoundary(dir) != 0)
                    {
                        std::cout << " (" << inode << "," << vd.nodeGridLocal.getBoundary(dir)
                                  << ")";
                    }
                }
                std::cout << std::endl;
            }
        }
    })

    /// prepare right number of buffers (one for every send operation)
    for (size_t coord = 0; coord < 3; ++coord)
    {
        for (size_t lr = 0; lr < 2; ++lr)
        {
            size_t const dir = 2 * coord + lr;
            sendBufReal[dir].clear();
            sendBufReal[dir].resize(commNodesReal[dir].size());
            sendBufGhost[dir].clear();
            sendBufGhost[dir].resize(commNodesGhost[dir].size());
            // sendBufSizeReal[dir].clear();
            // sendBufSizeReal[dir].resize(commNodesReal[dir].size());
            // sendBufSizeGhost[dir].clear();
            // sendBufSizeGhost[dir].resize(commNodesGhost[dir].size());
        }
    }
}

void DomainDecomposition::prepareGhostBuffers_channel()
{
    HPX4ESPP_ASSERT(channelsInit);
    /// TODO: Check if worth parallelizing across bufSizeList/resize
    /// prepare buffers for every [[send channel]] set of commNodes
    /// same shape as commNodesReal[dir] and commNodesGhost[dir]
    /// indexable by: [dir][isub in commNodes]
    for (size_t coord = 0; coord < 3; ++coord)
    {
        for (size_t lr = 0; lr < 2; ++lr)
        {
            size_t const dir = 2 * coord + lr;
            size_t const oppDir = 2 * coord + (1 - lr);

            auto f_resize_buf = [this](bool commReals, size_t dir) {
                auto& buf = commReals ? sendBufReal[dir] : sendBufGhost[dir];
                // auto& bufSizeList = commReals ? sendBufSizeReal[dir] : sendBufSizeGhost[dir];
                std::vector<size_t> const& cn =
                    commReals ? commNodesReal[dir] : commNodesGhost[dir];
                for (size_t in = 0; in < cn.size(); in++)
                {
                    const size_t& inode = cn[in];
                    const auto& vs = virtualStorage[inode];
                    const auto& cr = vs.particles.cellRange();
                    const auto& vd = vdd[inode];
                    const auto& cc = commReals ? vd.commCells[dir].reals : vd.commCells[dir].ghosts;
                    size_t nodeTotal = 0;
                    for (const auto& ic : cc)
                    {
                        nodeTotal += (cr[ic + 1] - cr[ic]);
                    }
                    const size_t bufSize = vecModeFactor * nodeTotal * sizeof(real);
                    // bufSizeList[in] = bufSize;
                    buf[in].resize(bufSize);
                }
            };
            f_resize_buf(true, dir);
            f_resize_buf(false, dir);
        }
    }
}

template <bool SIZES_FIRST, bool REAL_TO_GHOSTS, int EXTRA_DATA>
void DomainDecomposition::ghostCommunication_channel_impl()
{
    const int rank = mpiWorld->rank();

    for (size_t _coord = 0; _coord < 3; ++_coord)
    {
        int coord = REAL_TO_GHOSTS ? _coord : (2 - _coord);
        const real curCoordBoxL = getSystem()->bc->getBoxL()[coord];
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

            auto f_subNode_send = [this, &shift, coord, lr, dir, oppDir, rank,
                                   doPeriodic](size_t send_node) {
                /// get send mode info from nsubComm
                auto const& ns =
                    nsubComm[send_node + numSubNodes * (REAL_TO_GHOSTS ? dir : oppDir)];

                std::ostringstream errorMsg;
                HPX4ESPP_DDVI_CHANNELS_DEBUG({
                    errorMsg << "REAL_TO_GHOSTS: " << (REAL_TO_GHOSTS ? "true" : "false")
                             << " dir: " << dir << " (coord: " << coord << " lr: " << lr << ")";
                })

                /// internode
                if (std::get<0>(ns))
                {
                    const auto send_node_inter = std::get<1>(ns);
                    const int lr_ch = REAL_TO_GHOSTS ? lr : !lr;

                    /// Check that you are receiving from the correct node
                    HPX4ESPP_DDVI_CHANNELS_DEBUG({
                        const auto& chInfo = channels.sendIdx(coord, lr_ch, send_node_inter);
                        /// ^tuple{lr,sender0,selfid,index(0..numSubNodes)}

                        HPX4ESPP_ASSERT_EQUAL(std::get<0>(chInfo), lr_ch);
                        HPX4ESPP_ASSERT_EQUAL(std::get<1>(chInfo), rank);
                        HPX4ESPP_ASSERT_EQUAL(
                            std::get<2>(chInfo),
                            nodeGrid.getNodeNeighborIndex(REAL_TO_GHOSTS ? dir : oppDir));
                        HPX4ESPP_ASSERT_EQUAL(std::get<3>(chInfo), send_node_inter);
                    })

                    /// which buffer to use
                    auto& sendBuf = REAL_TO_GHOSTS ? sendBufReal[dir][send_node_inter]
                                                   : sendBufGhost[dir][send_node_inter];

                    /// check buffer size (debug)
                    HPX4ESPP_DDVI_CHANNELS_DEBUG({
                        const bool commReal = REAL_TO_GHOSTS;
                        const size_t idxCommNode = send_node_inter;  /// index in commnodes
                        const size_t inode = commReal ? commNodesReal[dir][idxCommNode]
                                                      : commNodesGhost[dir][idxCommNode];
                        HPX4ESPP_ASSERT_EQUAL_MSG(send_node, inode,
                                                  errorMsg.str());  /// index in virtual storage
                        const auto& vs = virtualStorage[inode];
                        const auto& cr = vs.particles.cellRange();
                        const auto& vd = vdd[inode];
                        const auto& cc =
                            commReal ? vd.commCells[dir].reals : vd.commCells[dir].ghosts;
                        size_t nodeTotal = 0;
                        for (const auto& ic : cc)
                        {
                            nodeTotal += cr[ic + 1] - cr[ic];
                        }
                        HPX4ESPP_ASSERT_EQUAL_MSG(vecModeFactor * nodeTotal * sizeof(real),
                                                  sendBuf.size(), errorMsg.str());
                    })

                    /// pack data to corresponding buffer
                    if (REAL_TO_GHOSTS)
                    {
                        packCells_channel<PACKED_POSITIONS, ADD_SHIFT>(sendBuf, true, dir,
                                                                       send_node, shift);
                    }
                    else
                    {
                        packCells_channel<PACKED_FORCES, NO_SHIFT>(sendBuf, false, dir, send_node,
                                                                   SHIFT_ZERO);
                    }

                    /// set channel value
                    /// TODO: Determine most efficient launch policy
                    // channels.sendChannel(coord,lr_ch,send_node_inter).set(sendBuf);
                    // channels.sendChannel(coord,lr_ch,send_node_inter).set(hpx::launch::apply,
                    // sendBuf);
                    return std::move(channels.sendChannel(coord, lr_ch, send_node_inter)
                                         .set(hpx::launch::async, sendBuf));
                }
                /// intranode
                else
                {
                    const auto src_node = send_node;
                    const auto dst_node = std::get<1>(ns);
                    if (REAL_TO_GHOSTS)
                    {
                        if (doPeriodic && (vdd[send_node].nodeGridLocal.getBoundary(dir) != 0))
                        {
                            /// periodic, i.e. in this direction, source node is at the edge and
                            /// there is only one subdomain
                            HPX4ESPP_NOT_IMPLEMENTED("Wrong signature of copyRealsToGhostsIntra");
                            // copyRealsToGhostsIntra<ADD_SHIFT>(dir, src_node, dst_node, shift);
                        }
                        else
                        {
                            HPX4ESPP_NOT_IMPLEMENTED("Wrong signature of copyRealsToGhostsIntra");
                            // copyRealsToGhostsIntra<NO_SHIFT>(dir, src_node, dst_node,
                            // SHIFT_ZERO);
                        }
                    }
                    else
                    {
                        HPX4ESPP_NOT_IMPLEMENTED("Wrong signature of addGhostForcesToRealsIntra");
                        // addGhostForcesToRealsIntra(dir, dst_node, src_node);
                    }
                }
                return std::move(hpx::make_ready_future());
            };

            auto f_send = [this, &f_subNode_send]() {
                std::vector<hpx::future<void>> sendFutures;
                for (size_t send_node = 0; send_node < numSubNodes; send_node++)
                {
                    sendFutures.push_back(hpx::async(f_subNode_send, send_node));
                }
                return std::move(sendFutures);
            };

            auto f_recv = [this, doPeriodic, coord, lr, dir, oppDir, rank]() {
                std::vector<hpx::future<void>> recvFutures;
                auto const& recvNodes = REAL_TO_GHOSTS ? commNodesGhost[dir] : commNodesReal[dir];
                const size_t numRecvNodes = doPeriodic ? 0 : recvNodes.size();
                recvFutures.reserve(numRecvNodes);

                /// push_back recv_futures first only if there are channels to receive from
                for (size_t recv_node = 0; recv_node < numRecvNodes; recv_node++)
                {
                    const int lr_ch = REAL_TO_GHOSTS ? lr : !lr;

                    /// Check that you are receiving from the correct node
                    HPX4ESPP_DDVI_CHANNELS_DEBUG({
                        const auto& chInfo = channels.recvIdx(coord, lr_ch, recv_node);
                        ///^ tuple{lr,sender0,selfid,index(0..numSubNodes)}
                        HPX4ESPP_ASSERT_EQUAL(std::get<0>(chInfo), lr_ch);
                        HPX4ESPP_ASSERT_EQUAL(
                            std::get<1>(chInfo),
                            nodeGrid.getNodeNeighborIndex(REAL_TO_GHOSTS ? oppDir : dir));
                        HPX4ESPP_ASSERT_EQUAL(std::get<2>(chInfo), rank);
                        HPX4ESPP_ASSERT_EQUAL(std::get<3>(chInfo), recv_node);
                    })

                    recvFutures.push_back(
                        channels.recvChannel(coord, lr_ch, recv_node)
                            .get()
                            .then(hpx::launch::async,
                                  // hpx::launch::sync,
                                  [this, recv_node, lr_ch, dir, oppDir, coord,
                                   lr](hpx::future<AlignedVectorChar>&& future) {
                                      AlignedVectorChar recvBuf = future.get();

                                      /// check buffer size
                                      HPX4ESPP_DDVI_CHANNELS_DEBUG({
                                          const bool commReal = !REAL_TO_GHOSTS;
                                          const size_t idxCommNode =
                                              recv_node;  /// index in commnodes
                                          const size_t inode =
                                              commReal ? commNodesReal[dir][idxCommNode]
                                                       : commNodesGhost[dir][idxCommNode];
                                          // HPX4ESPP_ASSERT_EQUAL(send_node,inode); /// index in
                                          // virtual storage
                                          const auto& vs = virtualStorage[inode];
                                          const auto& cr = vs.particles.cellRange();
                                          const auto& vd = vdd[inode];
                                          const auto& cc = commReal ? vd.commCells[dir].reals
                                                                    : vd.commCells[dir].ghosts;
                                          size_t nodeTotal = 0;
                                          for (const auto& ic : cc)
                                          {
                                              nodeTotal += cr[ic + 1] - cr[ic];
                                          }
                                          HPX4ESPP_ASSERT_EQUAL_MSG(
                                              vecModeFactor * nodeTotal * sizeof(real),
                                              recvBuf.size(),
                                              "Failed at"
                                                  << " REAL_TO_GHOSTS: " << REAL_TO_GHOSTS
                                                  << " coord: " << coord << " lr: " << lr);
                                      })

                                      const size_t inode = REAL_TO_GHOSTS
                                                               ? commNodesGhost[dir][recv_node]
                                                               : commNodesReal[dir][recv_node];

                                      /// unpack data from corresponding buffer
                                      if (REAL_TO_GHOSTS)
                                      {
                                          unpackCells_channel<PACKED_POSITIONS, DATA_INSERT>(
                                              recvBuf, false, dir, inode);
                                      }
                                      else
                                      {
                                          unpackCells_channel<PACKED_FORCES, DATA_ADD>(
                                              recvBuf, true, dir, inode);
                                      }
                                  }));
                }

                return std::move(recvFutures);
            };

            if (nodeGrid.getNodePosition(coord) % 2 == 0)
            {
                std::vector<hpx::future<void>> sendFutures = f_send();
                hpx::wait_all(sendFutures);
                std::vector<hpx::future<void>> recvFutures = f_recv();
                hpx::wait_all(recvFutures);
            }
            else
            {
                std::vector<hpx::future<void>> recvFutures = f_recv();
                hpx::wait_all(recvFutures);
                std::vector<hpx::future<void>> sendFutures = f_send();
                hpx::wait_all(sendFutures);
            }
        }
    }
}

template void DomainDecomposition::ghostCommunication_channel_impl<false, true, 0>();
template void DomainDecomposition::ghostCommunication_channel_impl<false, false, 0>();

template <DomainDecomposition::PackedData PACKED_DATA, DomainDecomposition::AddShift DO_SHIFT>
void DomainDecomposition::packCells_channel(
    vec::AlignedVector<char>& sendBuf, bool commReal, size_t dir, size_t inode, Real3D const& shift)
{
    const auto& vs = virtualStorage[inode];
    const auto& cr = vs.particles.cellRange();
    const auto& vd = vdd[inode];
    const auto& cc = commReal ? vd.commCells[dir].reals : vd.commCells[dir].ghosts;

    size_t expNumParts;
    HPX4ESPP_DDVI_CHANNELS_DEBUG({
        const size_t maxBufReals = sendBuf.size() / sizeof(real);
        expNumParts = maxBufReals / vecModeFactor;
        HPX4ESPP_ASSERT_EQUAL(sendBuf.size() % sizeof(real), 0);
        HPX4ESPP_ASSERT_EQUAL(maxBufReals % vecModeFactor, 0);
    })

    {
        /// loop over dimensions (x,y,z)
        size_t b_start = 0;
        auto f_pack_dim = [dir, &cc, &cr, &sendBuf, &b_start, expNumParts](
                              size_t dim, const real* __restrict p_ptr, real shift_v) {
            // const size_t b_start = (nodeStart * 3) + (numPart * dim);
            real* __restrict b_ptr = (real*)(sendBuf.data()) + b_start;

            /// loop over cells
            size_t b_off = 0;
            for (const auto& ic : cc)
            {
                real* __restrict b_ptr_c = b_ptr + b_off;
                const real* __restrict p_ptr_c = p_ptr + cr[ic];
                const size_t npart = cr[ic + 1] - cr[ic];

                /// loop over particles
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
            b_start += b_off;
            HPX4ESPP_DDVI_CHANNELS_DEBUG(HPX4ESPP_ASSERT_EQUAL(expNumParts, b_off))
        };

        if (PACKED_DATA == PACKED_POSITIONS)
        {
            f_pack_dim(0, vs.particles.p_x.data(), shift[0]);
            f_pack_dim(1, vs.particles.p_y.data(), shift[1]);
            f_pack_dim(2, vs.particles.p_z.data(), shift[2]);
        }
        else
        {
            f_pack_dim(0, vs.particles.f_x.data(), shift[0]);
            f_pack_dim(1, vs.particles.f_y.data(), shift[1]);
            f_pack_dim(2, vs.particles.f_z.data(), shift[2]);
        }
    }
}

template void DomainDecomposition::packCells_channel<DomainDecomposition::PACKED_POSITIONS,
                                                     DomainDecomposition::ADD_SHIFT>(
    vec::AlignedVector<char>& sendBuf,
    bool commReal,
    size_t dir,
    size_t inode,
    Real3D const& shift);

template void DomainDecomposition::packCells_channel<DomainDecomposition::PACKED_FORCES,
                                                     DomainDecomposition::NO_SHIFT>(
    vec::AlignedVector<char>& sendBuf,
    bool commReal,
    size_t dir,
    size_t inode,
    Real3D const& shift);

template <DomainDecomposition::PackedData PACKED_DATA, DomainDecomposition::DataMode DATA_MODE>
void DomainDecomposition::unpackCells_channel(vec::AlignedVector<char> const& recvBuf,
                                              bool commReal,
                                              size_t dir,
                                              size_t inode)
{
    auto& vs = virtualStorage[inode];
    const auto& cr = vs.particles.cellRange();
    const auto& vd = vdd[inode];
    const auto& cc = commReal ? vd.commCells[dir].reals : vd.commCells[dir].ghosts;

    size_t expNumParts;
    HPX4ESPP_DDVI_CHANNELS_DEBUG({
        const size_t maxBufReals = recvBuf.size() / sizeof(real);
        expNumParts = maxBufReals / vecModeFactor;
        HPX4ESPP_ASSERT_EQUAL(recvBuf.size() % sizeof(real), 0);
        HPX4ESPP_ASSERT_EQUAL(maxBufReals % vecModeFactor, 0);
    })

    {
        /// loop over dimensions (x,y,z)
        size_t b_start = 0;
        auto f_pack_dim = [dir, inode, &cc, &cr, &recvBuf, &b_start, expNumParts](
                              size_t dim, real* __restrict p_ptr) {
            const real* __restrict b_ptr = (real*)(recvBuf.data()) + b_start;

            /// loop over cells
            size_t b_off = 0;
            for (const auto& ic : cc)
            {
                const real* __restrict b_ptr_c = b_ptr + b_off;
                real* __restrict p_ptr_c = p_ptr + cr[ic];
                const size_t npart = cr[ic + 1] - cr[ic];

                /// loop over particles
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
            b_start += b_off;
            HPX4ESPP_DDVI_CHANNELS_DEBUG(HPX4ESPP_ASSERT_EQUAL(expNumParts, b_off))
        };

        if (PACKED_DATA == PACKED_POSITIONS)
        {
            f_pack_dim(0, vs.particles.p_x.data());
            f_pack_dim(1, vs.particles.p_y.data());
            f_pack_dim(2, vs.particles.p_z.data());
        }
        else
        {
            f_pack_dim(0, vs.particles.f_x.data());
            f_pack_dim(1, vs.particles.f_y.data());
            f_pack_dim(2, vs.particles.f_z.data());
        }
    }
}

template void DomainDecomposition::unpackCells_channel<DomainDecomposition::PACKED_POSITIONS,
                                                       DomainDecomposition::DATA_INSERT>(
    vec::AlignedVector<char> const& recvBuf, bool commReal, size_t dir, size_t inode);

template void DomainDecomposition::unpackCells_channel<DomainDecomposition::PACKED_FORCES,
                                                       DomainDecomposition::DATA_ADD>(
    vec::AlignedVector<char> const& recvBuf, bool commReal, size_t dir, size_t inode);

}  // namespace storage
}  // namespace hpx4espp
}  // namespace espressopp
