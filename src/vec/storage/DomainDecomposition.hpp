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

#ifndef VEC_STORAGE_DOMAINDECOMPOSITION_HPP
#define VEC_STORAGE_DOMAINDECOMPOSITION_HPP

#include "StorageVec.hpp"
#include "vec/ParticleArray.hpp"
#include "vec/include/types.hpp"

#include <iostream>
#include "log4espp.hpp"
#include "storage/DomainDecomposition.hpp"
#include "esutil/Timer.hpp"
// #include "vec/BufferView.hpp"

#include "types.hpp"

/// Forward declaration
namespace espressopp { namespace integrator {
  class MDIntegrator;
}}

namespace espressopp { namespace vec {
  namespace storage {

#if 0
    struct VirtualDomainDecomposition
    {
      espressopp::storage::NodeGrid nodeGrid; /// index wrt entire system domain
      espressopp::storage::NodeGrid nodeGridLocal; /// index wrt locality subdomain

      espressopp::CellGrid cellGrid;

      struct CommCellIdx {
        std::vector<size_t> reals;
        std::vector<size_t> ghosts;
      };

      std::array< CommCellIdx, 6 > commCells;

      LocalCellList vGhostCells;
      CellList vLocalCells;
      CellList vRealCells;
    };

    typedef VirtualDomainDecomposition VDD;
#endif

    class DomainDecomposition:
      public espressopp::storage::DomainDecomposition,
      public StorageVec
    {
    public:
      typedef espressopp::storage::DomainDecomposition baseClass;

      DomainDecomposition(
        shared_ptr< Vectorization > vectorization,
        const Int3D& nodeGrid,
        const Int3D& cellGrid,
        int halfCellInt
        );

      ~DomainDecomposition();

      void connect();

      void disconnect();

      void resetCells();

      void loadCells();

      void unloadCells();

      void updateGhostsVec();

      void collectGhostForcesVec();

      static void registerPython();

    protected:

      boost::signals2::connection sigResetCells;

      const Mode vecMode;

      /////////////////////////////////////////////////////////////////////////////////////////////
      //// members involved in ghost communication
      struct CommCellIdx {
        std::vector<size_t> reals;
        std::vector<size_t> ghosts;
        size_t numReals;  /// buffer size of reals
        size_t numGhosts; /// buffer size of ghosts
      };

      std::array< CommCellIdx, 6 > commCellIdx;

      const size_t vecModeFactor;

      AlignedVector< real > buffReal, buffGhost;

      void prepareGhostBuffers();

      template < bool SIZES_FIRST, bool REAL_TO_GHOSTS, int EXTRA_DATA >
      void ghostCommunication_impl();

      enum PackedData { PACKED_POSITIONS=0, PACKED_FORCES=1 };
      enum DataMode { DATA_INSERT=0, DATA_ADD=1 };
      enum AddShift { NO_SHIFT=0, ADD_SHIFT=1 };
      static const espressopp::Real3D SHIFT_ZERO;

      template< AddShift DO_SHIFT >
      void copyRealsToGhostsIntra(size_t dir, size_t ir, size_t ig, Real3D const& shift);

      void addGhostForcesToRealsIntra(size_t dir, size_t ir, size_t ig);

      template<
        PackedData PACKED_DATA,
        AddShift DO_SHIFT >
      void packCells(
        AlignedVector<real> & sendBuf,
        bool commReal,
        size_t dir,
        size_t idxCommNode,
        Real3D const& shift
        );

      template<
        PackedData PACKED_DATA,
        DataMode DATA_MODE >
      void unpackCells(
        AlignedVector<real> const& recvBuf,
        bool commReal,
        size_t dir,
        size_t idxCommNode
        );

      /////////////////////////////////////////////////////////////////////////////////////////////

#if 0

      void updateGhostsBlocking();

      void collectGhostForcesBlocking();

      void resetVirtualStorage();

      void decomposeVec();

    protected:

      std::array<size_t,3> subCellGrid;

      std::array<size_t,3> subNodeGrid;

      size_t numSubNodes;

      size_t numSubCells;

      std::vector<std::vector<size_t>> resortWaves;

      std::vector<VDD> vdd;

      std::array< std::vector< size_t >, 6 > commNodesReal;
      std::array< std::vector< size_t >, 6 > commNodesGhost;
      std::array< std::vector< size_t >, 6 > nodeRangeReal;
      std::array< std::vector< size_t >, 6 > nodeRangeGhost;
      std::array< std::vector< std::tuple<size_t,size_t,bool> >, 6 > subNodePairsIntra;
      size_t maxReal, maxGhost;
      AlignedVector< real > buffReal, buffGhost;

      /// range of cells for each node
      std::array< std::vector< size_t >, 6 > nodeCommCellRange;
      ParticleList decompSendBuf, decompRecvBuf;

      bool commAsync;

      void prepareChannelAttributes();

      void prepareGhostBuffers_channel();

      template < bool SIZES_FIRST, bool REAL_TO_GHOSTS, int EXTRA_DATA >
      void ghostCommunication_channel_impl();

      std::vector<std::tuple<int,int>> nsubComm;
      std::array<std::vector<AlignedVectorChar>, 6> sendBufReal;
      std::array<std::vector<AlignedVectorChar>, 6> sendBufGhost;
      std::array<std::vector<size_t>, 6> sendBufSizeReal;
      std::array<std::vector<size_t>, 6> sendBufSizeGhost;

      template<
        PackedData PACKED_DATA,
        AddShift DO_SHIFT >
      void packCells_channel(
        AlignedVector<char> & sendBuf,
        bool commReal,
        size_t dir,
        size_t inode,
        Real3D const& shift
        );

      template<
        PackedData PACKED_DATA,
        DataMode DATA_MODE >
      void unpackCells_channel(
        AlignedVector<char> const& recvBuf,
        bool commReal,
        size_t dir,
        size_t inode
        );

    public:
      python::object getChannelIndices() const;

    protected:
      bool commUseChannels;
      bool channelsInit=false;
      Channels channels;


    /** Members extending base DomainDecomposition class */
    protected:

      BufferFixed inBuf, outBuf;
      std::vector<InBufferView> inBufView;
      std::vector<OutBufferView> outBufView;
      // std::vector<BufferView> inBufView, outBufView;

      template< bool ALIGNED >
      void exchangeGhosts_impl();
      bool excgAligned = false;

      void exchangeGhosts_SingleNode();

      bool decompUseParFor;

      void decomposeRealParticlesVecCellTask_MultiNode_NonPeriodic();

      void decomposeRealParticlesVecParFor_MultiNode_NonPeriodic();

    /** Utility functions for connecting to integrator to do stepwise offload */
    public:
      void connectOffload(
        boost::shared_ptr<espressopp::integrator::MDIntegrator> mdintegrator
        );

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

#endif

    private:
      static LOG4ESPP_DECL_LOGGER(logger);
    };

  }
}}

#endif//VEC_STORAGE_DOMAINDECOMPOSITION_HPP
