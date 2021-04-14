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

#include "vec/include/types.hpp"
#include "vec/ParticleArray.hpp"
#include "StorageVec.hpp"

#include <iostream>
#include "log4espp.hpp"
#include "storage/DomainDecomposition.hpp"
#include "esutil/Timer.hpp"

#include "types.hpp"

namespace espressopp { namespace vec {
  namespace storage {

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
        int halfCellInt,
        bool rebuildLocalParticles
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

      bool rebuildLocalParticles;

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
      void copyRealsToGhostsIntra(size_t dir, Real3D const& shift);

      void addGhostForcesToRealsIntra(size_t dir);

      template<
        PackedData PACKED_DATA,
        AddShift DO_SHIFT >
      void packCells(
        AlignedVector<real> & sendBuf,
        bool commReal,
        size_t dir,
        Real3D const& shift
        );

      template<
        PackedData PACKED_DATA,
        DataMode DATA_MODE >
      void unpackCells(
        AlignedVector<real> const& recvBuf,
        bool commReal,
        size_t dir
        );

      /////////////////////////////////////////////////////////////////////////////////////////////

    private:
      static LOG4ESPP_DECL_LOGGER(logger);
    };

  }
}}

#endif//VEC_STORAGE_DOMAINDECOMPOSITION_HPP
