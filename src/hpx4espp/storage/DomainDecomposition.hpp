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

#ifndef HPX4ESPP_STORAGE_DOMAINDECOMPOSITION_HPP
#define HPX4ESPP_STORAGE_DOMAINDECOMPOSITION_HPP

#include <iostream>
#include "log4espp.hpp"
#include "storage/DomainDecomposition.hpp"
#include "StorageHPX.hpp"
#include "Channels.hpp"
#include "esutil/Timer.hpp"
#include "hpx4espp/BufferView.hpp"

// If false, revert to original domain decomp
// #define HPX4ESPP_OPTIMIZE_COMM (1)

/// Forward declaration
namespace espressopp
{
namespace integrator
{
class MDIntegrator;
}
}  // namespace espressopp

namespace espressopp
{
namespace hpx4espp
{
namespace storage
{
struct VirtualDomainDecomposition
{
    espressopp::storage::NodeGrid nodeGrid;       /// index wrt entire system domain
    espressopp::storage::NodeGrid nodeGridLocal;  /// index wrt locality subdomain

    espressopp::CellGrid cellGrid;

    struct CellInfo
    {
        enum CellType
        {
            CELL_REAL,
            CELL_GHOST_EXTERNAL,
            CELL_GHOST_INTERNAL
        };
        CellType cellType;
        size_t subNode;
        size_t subNodeCell;
    };

    std::vector<CellInfo> cellGridInfo;

    struct CommCellIdx
    {
        std::vector<size_t> reals;
        std::vector<size_t> ghosts;
    };

    std::array<CommCellIdx, 6> commCells;
    std::array<CommCellIdx, 6> commCellsOwned;

    std::array<std::vector<size_t>, 6> numCommSubCells;
    std::array<std::vector<std::pair<size_t, size_t>>, 6> commSubCellRanges;
    std::array<std::vector<size_t>, 6> subCommRangeReal;
    std::array<std::vector<size_t>, 6> subCommRangeGhost;

    LocalCellList vGhostCells;
    CellList vLocalCells;
    CellList vRealCells;
};

typedef VirtualDomainDecomposition VDD;

class DomainDecomposition : public espressopp::storage::DomainDecomposition, public StorageHPX
{
public:
    typedef espressopp::storage::DomainDecomposition baseClass;

    bool HPX4ESPP_OPTIMIZE_COMM;

    DomainDecomposition(shared_ptr<System> system,
                        const Int3D& nodeGrid,
                        const Int3D& cellGrid,
                        int halfCellInt,
                        const Int3D& subCellGrid,
                        int numCommSubs,
                        bool commAsync,
                        bool excgAligned,
                        bool commUseChannels,
                        bool decompUseParFor,
                        bool _HPX4ESPP_OPTIMIZE_COMM);

    ~DomainDecomposition();

    void initChannels();

    /** Copy particles to packed form. To be called at the start of integrator.run */
    void loadCells();

    /** Copy particles back from packed form. To be called at the end of integrator.run */
    void unloadCells();

    void prepareGhostBuffers();

    void updateGhostsBlocking();

    void collectGhostForcesBlocking();

    void connect();

    void disconnect();

    void resetVirtualStorage();

    void decomposeHPX();

    static void registerPython();

protected:
    std::array<size_t, 3> subCellGrid;

    std::array<size_t, 3> subNodeGrid;

    size_t numSubNodes;

    size_t numSubCells;

    std::vector<std::vector<size_t>> resortWaves;

    std::vector<VDD> vdd;

    std::array<std::vector<size_t>, 6> commNodesReal;
    std::array<std::vector<size_t>, 6> commNodesGhost;
    std::array<std::vector<size_t>, 6> nodeRangeReal;
    std::array<std::vector<size_t>, 6> nodeRangeGhost;
    std::array<std::vector<std::tuple<size_t, size_t, bool>>, 6> subNodePairsIntra;
    std::array<std::vector<std::tuple<size_t, size_t, bool>>, 6> subNodePairsIntraPeriodic;
    size_t maxReal, maxGhost;
    vec::AlignedVector<real> buffReal, buffGhost;

    /// range of cells for each node
    std::array<std::vector<size_t>, 6> nodeCommCellRange;
    ParticleList decompSendBuf, decompRecvBuf;

    template <bool SIZES_FIRST, bool REAL_TO_GHOSTS, int EXTRA_DATA>
    void ghostCommunication_impl();
    bool commAsync;

    void prepareChannelAttributes();

    void prepareGhostBuffers_channel();

    template <bool SIZES_FIRST, bool REAL_TO_GHOSTS, int EXTRA_DATA>
    void ghostCommunication_channel_impl();

    std::vector<std::tuple<int, int>> nsubComm;
    std::array<std::vector<AlignedVectorChar>, 6> sendBufReal;
    std::array<std::vector<AlignedVectorChar>, 6> sendBufGhost;
    std::array<std::vector<size_t>, 6> sendBufSizeReal;
    std::array<std::vector<size_t>, 6> sendBufSizeGhost;

    enum PackedData
    {
        PACKED_POSITIONS = 0,
        PACKED_FORCES = 1
    };
    enum DataMode
    {
        DATA_INSERT = 0,
        DATA_ADD = 1
    };
    enum AddShift
    {
        NO_SHIFT = 0,
        ADD_SHIFT = 1
    };
    static const espressopp::Real3D SHIFT_ZERO;

    template <PackedData PACKED_DATA, AddShift DO_SHIFT>
    void packCells_channel(vec::AlignedVector<char>& sendBuf,
                           bool commReal,
                           size_t dir,
                           size_t inode,
                           Real3D const& shift);

    template <PackedData PACKED_DATA, DataMode DATA_MODE>
    void unpackCells_channel(vec::AlignedVector<char> const& recvBuf,
                             bool commReal,
                             size_t dir,
                             size_t inode);

    const size_t numCommSubs;

public:
    python::object getChannelIndices() const;

protected:
    bool commUseChannels;
    bool channelsInit = false;
    Channels channels;

    template <AddShift DO_SHIFT>
    void copyRealsToGhostsIntra(size_t dir, size_t ir, size_t ig, size_t is, Real3D const& shift);

    void addGhostForcesToRealsIntra(size_t dir, size_t ir, size_t ig, size_t is);

    template <PackedData PACKED_DATA, AddShift DO_SHIFT>
    void packCells(vec::AlignedVector<real>& sendBuf,
                   bool commReal,
                   size_t dir,
                   size_t idxCommNode,
                   Real3D const& shift);

    template <PackedData PACKED_DATA, DataMode DATA_MODE>
    void unpackCells(vec::AlignedVector<real> const& recvBuf,
                     bool commReal,
                     size_t dir,
                     size_t idxCommNode);

    static constexpr size_t vecModeFactor = 3;

    /** Members extending base DomainDecomposition class */
protected:
    BufferFixed inBuf, outBuf;
    std::vector<InBufferView> inBufView;
    std::vector<OutBufferView> outBufView;
    // std::vector<BufferView> inBufView, outBufView;

    template <bool ALIGNED>
    void exchangeGhosts_impl();
    bool excgAligned = false;

    void exchangeGhosts_SingleNode();

    bool decompUseParFor;

    void decomposeRealParticlesHPXCellTask_MultiNode_NonPeriodic();

    void decomposeRealParticlesHPXParFor_MultiNode_NonPeriodic();

    /** Utility functions for connecting to integrator to do stepwise offload */
public:
    void connectOffload(boost::shared_ptr<espressopp::integrator::MDIntegrator> mdintegrator);

    void disconnectOffload();

    bool connectedOffload() const { return offload; }

protected:
    void resetParticles();

    void resetCells();

    void befCalcForces();

    void updateForces();

    // signals that connect to integrator
    boost::signals2::connection sigBefCalcForces;
    boost::signals2::connection sigUpdateForces;

    // signals that connect to system
    boost::signals2::connection sigResetParticles;
    boost::signals2::connection sigResetCells;

    bool offload = false;

    /** Timer-related  */
public:
    void loadTimers(real t[2]);

    void resetTimers();

protected:
    real timeUpdateGhostsBlocking_comm;
    real timeCollectGhostForcesBlocking_comm;
    real timeParticlesCopyFrom;
    real timeParticlesUpdateTo;
    real timeDecomposeInvGhosts;
    real timeDecomposeReal;
    real timeDecomposeExcGhosts;
    real timeDecomposeSignal;
    real timeDecomposeRealResort;
    real timeDecomposeRealComm;
    real timeLoadLoop;
    real timeLoadPrepGhost;
    real timeLoadSignal;
    real timeUnloadLoop;
    real timeUnloadSignal;
    real timeUpdateGhosts_InterNode_pack;
    real timeUpdateGhosts_InterNode_comm;
    real timeUpdateGhosts_InterNode_unpack;
    real timeUpdateGhosts_IntraNode;
    real timeExcGhosts_InterNode_sizes;
    real timeExcGhosts_InterNode_pack;
    real timeExcGhosts_InterNode_comm;
    real timeExcGhosts_InterNode_unpack;
    real timeExcGhosts_IntraNode;

    espressopp::esutil::WallTimer wallTimer;

    /** Provides interace to testing functions */
public:
    std::vector<VDD> getVDD() const { return vdd; }
    std::vector<VDD::CommCellIdx> getCommCellsAsIdx() const;

private:
    static LOG4ESPP_DECL_LOGGER(logger);
};

}  // namespace storage
}  // namespace hpx4espp
}  // namespace espressopp

#endif  // HPX4ESPP_STORAGE_DOMAINDECOMPOSITION_HPP
