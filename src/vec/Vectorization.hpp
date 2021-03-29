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

#ifndef VEC_VECTORIZATION_HPP
#define VEC_VECTORIZATION_HPP

#include "SystemAccess.hpp"
#include "integrator/MDIntegrator.hpp"
#include "storage/DomainDecomposition.hpp"
#include "ParticleArray.hpp"
#include "CellNeighborList.hpp"

#include "types.hpp"
#include "log4espp.hpp"
#include "boost/signals2.hpp"

#include "vec/include/types.hpp"

namespace espressopp {
  namespace vec {

    ///////////////////////////////////////////////////////////////////////////////////////////////
    /// facilitates offloading of particle data to vectorization-friendly form
    class Vectorization
      : public SystemAccess
      , public boost::enable_shared_from_this<Vectorization>
    {
      typedef espressopp::storage::Storage Storage;
      typedef espressopp::vec::storage::StorageVec StorageVec;

      typedef espressopp::integrator::MDIntegrator MDIntegrator;
      typedef espressopp::vec::integrator::MDIntegratorVec MDIntegratorVec;

    public:
      Vectorization(
        shared_ptr<System> system,
        shared_ptr<Storage> storage,
        shared_ptr<MDIntegrator> mdintegrator,
        Mode vecMode = ESPP_VEC_MODE_DEFAULT);

      Vectorization(
        shared_ptr<System> system,
        shared_ptr<StorageVec> storageVec,
        shared_ptr<MDIntegratorVec> mdintegratorVec,
        Mode vecMode = ESPP_VEC_MODE_DEFAULT);

      ~Vectorization();

      void connect_level_1();
      void connect_level_2();
      void disconnect();

      ParticleArray particles;
      CellNeighborList neighborList;

      int getVecLevel() {
        return vecLevel;
      }

      static void registerPython();

    protected:
      const Mode vecMode;
      const int vecLevel;

      void resetParticles();
      void befCalcForces();
      void updatePositions();
      void updateForces();

      void resetCells();

      shared_ptr<Storage> storage;
      shared_ptr<StorageVec> storageVec;

      shared_ptr<MDIntegrator> mdintegrator;
      shared_ptr<MDIntegratorVec> mdintegratorVec;

      // signals that connect to integrator
      boost::signals2::connection sigBefCalcForces;
      boost::signals2::connection sigUpdateForces;

      // signals that connect to system
      boost::signals2::connection sigResetParticles;
      boost::signals2::connection sigResetCells;

      static LOG4ESPP_DECL_LOGGER(logger);
    };

  }
}

#endif // VEC_VECTORIZATION_HPP
