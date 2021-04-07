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

namespace espressopp { namespace vec {
  namespace storage {

    LOG4ESPP_LOGGER(DomainDecomposition::logger, "DomainDecomposition");

    // const Real3D DomainDecomposition::SHIFT_ZERO {0.,0.,0.};

    DomainDecomposition::DomainDecomposition(
      shared_ptr< Vectorization > vectorization,
      const Int3D& _nodeGrid,
      const Int3D& _cellGrid,
      int _halfCellInt
      )
      : baseClass(vectorization->getSystem(), _nodeGrid, _cellGrid, _halfCellInt),
        StorageVec(vectorization), vecMode(vectorization->getVecMode()),
        vecModeFactor((vecMode==ESPP_VEC_SOA) ? 3 : 4)
    {
      if(halfCellInt!=1)
        throw std::runtime_error("vec: Not implemented for halfCellInt!=1.");
      // resetTimers();
      connect();
      resetCells();
      LOG4ESPP_INFO(logger, "DomainDecomposition()");
    }

    DomainDecomposition::~DomainDecomposition()
    {
      disconnect();
    }

    void DomainDecomposition::connect()
    {
      if(!sigResetCells.connected())
        sigResetCells = this->onCellAdjust.connect(
                      boost::signals2::at_back,
                      boost::bind(&DomainDecomposition::resetCells, this));
      LOG4ESPP_INFO(logger, "DomainDecomposition::connect()");
    }

    void DomainDecomposition::disconnect()
    {
      if(sigResetCells.connected())
        sigResetCells.disconnect();
      LOG4ESPP_INFO(logger, "DomainDecomposition::disconnect()");
    }

    /// Copy particles to packed form. To be called at the start of integrator.run
    void DomainDecomposition::loadCells()
    {
      vectorization->resetParticles();
      prepareGhostBuffers();
      onLoadCells();
    }

    /// Copy particles back from packed form. To be called at the end of integrator.run
    void DomainDecomposition::unloadCells()
    {
      vectorization->particles.updateToPositionVelocity(localCells, true);
      onUnloadCells();
    }

    void DomainDecomposition::resetCells()
    {
      /// reset realCells and cellNeighborList
      vectorization->resetCells(this);

      /// Setup ghost communication for this cell grid
      {
        /// convert commCells to commCellsIdx
        for(int dir=0; dir<6; dir++){
          commCellIdx[dir].reals  = CellListToIdx(commCells[dir].reals,  localCells[0]);
          commCellIdx[dir].ghosts = CellListToIdx(commCells[dir].ghosts, localCells[0]);
        }
        /// TODO: Check whether sorting improves performance
      }
    }

    ///////////////////////////////////////////////////////////////////////////////////////////////
    /// preallocate buffers for ghost communication
    void DomainDecomposition::prepareGhostBuffers()
    {
      const auto& cr = vectorization->particles.cellRange();
      size_t maxReals = 0, maxGhosts = 0;
      for (size_t coord = 0; coord < 3; ++coord) {
        for (size_t lr = 0; lr < 2; ++lr) {
          size_t const dir    = 2 * coord + lr;
          size_t const oppDir = 2 * coord + (1-lr);

          auto f_countParticles = [cr](auto const& cellIdx) {
            size_t total = 0;
            for(auto const& ic: cellIdx)
              total += (cr[ic+1]-cr[ic]);
            return total;
          };

          auto& cci = commCellIdx[dir];
          cci.numReals  = f_countParticles(cci.reals);
          cci.numGhosts = f_countParticles(cci.ghosts);
          maxReals  = std::max(cci.numReals,  maxReals);
          maxGhosts = std::max(cci.numGhosts, maxGhosts);
        }
      }

      const size_t preallocReal  = maxReals  * vecModeFactor;
      const size_t preallocGhost = maxGhosts * vecModeFactor;
      buffReal.resize(preallocReal);
      buffGhost.resize(preallocGhost);
    }

    ///////////////////////////////////////////////////////////////////////////////////////////////
    void DomainDecomposition::updateGhostsVec()
    {

    }

    ///////////////////////////////////////////////////////////////////////////////////////////////
    void DomainDecomposition::collectGhostForcesVec()
    {

    }

    ///////////////////////////////////////////////////////////////////////////////////////////////
    template < bool SIZES_FIRST, bool REAL_TO_GHOSTS, int EXTRA_DATA >
    void DomainDecomposition::ghostCommunication_impl()
    {

    }

    template void DomainDecomposition::ghostCommunication_impl<false, true, 0>();

    template void DomainDecomposition::ghostCommunication_impl<false, false, 0>();

    ///////////////////////////////////////////////////////////////////////////////////////////////
    template< DomainDecomposition::AddShift DO_SHIFT >
    void DomainDecomposition::copyRealsToGhostsIntra(
      size_t dir, size_t ir, size_t ig, Real3D const& shift
      )
    {

    }

    template void DomainDecomposition::copyRealsToGhostsIntra< DomainDecomposition::ADD_SHIFT >(
      size_t dir, size_t ir, size_t ig, Real3D const& shift);

    template void DomainDecomposition::copyRealsToGhostsIntra< DomainDecomposition::NO_SHIFT >(
      size_t dir, size_t ir, size_t ig, Real3D const& shift);

    ///////////////////////////////////////////////////////////////////////////////////////////////
    void DomainDecomposition::addGhostForcesToRealsIntra(size_t dir, size_t ir, size_t ig)
    {

    }

    ///////////////////////////////////////////////////////////////////////////////////////////////
    template<
      DomainDecomposition::PackedData PACKED_DATA,
      DomainDecomposition::AddShift DO_SHIFT >
    void DomainDecomposition::packCells(
      AlignedVector<real> & sendBuf,
      bool commReal,
      size_t dir,
      size_t idxCommNode,
      Real3D const& shift
      )
    {

    }

    template void DomainDecomposition::packCells<
      DomainDecomposition::PACKED_POSITIONS,
      DomainDecomposition::ADD_SHIFT>
    (
      AlignedVector<real> & sendBuf,
      bool commReal,
      size_t dir,
      size_t idxCommNode,
      Real3D const& shift
    );

    template void DomainDecomposition::packCells<
      DomainDecomposition::PACKED_FORCES,
      DomainDecomposition::NO_SHIFT>
    (
      AlignedVector<real> & sendBuf,
      bool commReal,
      size_t dir,
      size_t idxCommNode,
      Real3D const& shift
    );

    ///////////////////////////////////////////////////////////////////////////////////////////////
    template<
      DomainDecomposition::PackedData PACKED_DATA,
      DomainDecomposition::DataMode DATA_MODE >
    void DomainDecomposition::unpackCells(
      AlignedVector<real> const& recvBuf,
      bool commReal,
      size_t dir,
      size_t idxCommNode
      )
    {

    }

    template void DomainDecomposition::unpackCells<
      DomainDecomposition::PACKED_POSITIONS,
      DomainDecomposition::DATA_INSERT>
    (
      AlignedVector<real> const& recvBuf,
      bool commReal,
      size_t dir,
      size_t idxCommNode
    );

    template void DomainDecomposition::unpackCells<
      DomainDecomposition::PACKED_FORCES,
      DomainDecomposition::DATA_ADD>
    (
      AlignedVector<real> const& recvBuf,
      bool commReal,
      size_t dir,
      size_t idxCommNode
    );

    ///////////////////////////////////////////////////////////////////////////////////////////////
    void DomainDecomposition::registerPython()
    {
      using namespace espressopp::python;

      class_< DomainDecomposition, bases<espressopp::storage::DomainDecomposition, StorageVec >, boost::noncopyable >
        ("vec_storage_DomainDecomposition", init< shared_ptr< Vectorization >, const Int3D&,
            const Int3D&, int >())
        // .def("initChannels", &DomainDecomposition::initChannels)
        // .def("getChannelIndices", &DomainDecomposition::getChannelIndices)
        // .def("connectOffload", &DomainDecomposition::connectOffload)
        // .def("connectedOffload", &DomainDecomposition::connectedOffload)
        // .def("disconnectOffload", &DomainDecomposition::disconnectOffload)
        // .def("resetCells", &DomainDecomposition::resetCells)
        // .def("resetTimers", &DomainDecomposition::resetTimers)
        // .def("getTimers", &wrapGetTimers)
        // .def("getTimers2", &wrapGetTimers2)
        ;
    }

  }
}}

/// WARNING: Storage::localParticles is unused in this implementation
