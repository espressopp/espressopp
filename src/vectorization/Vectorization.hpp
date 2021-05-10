/*
  Copyright (C) 2019-2020
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

#ifndef _VECTORIZATION_VECTORIZATION_HPP
#define _VECTORIZATION_VECTORIZATION_HPP

#include "SystemAccess.hpp"
#include "integrator/MDIntegrator.hpp"
#include "storage/DomainDecomposition.hpp"
#include "ParticleArray.hpp"
#include "CellNeighborList.hpp"

#include "types.hpp"
#include "log4espp.hpp"
#include "boost/signals2.hpp"

namespace espressopp {
  namespace vectorization {

    ///////////////////////////////////////////////////////////////////////////////////////////////
    /// facilitates offloading of particle data to vectorization-friendly form
    class Vectorization : public SystemAccess
    {
    public:
      Vectorization(std::shared_ptr<System> system, std::shared_ptr<integrator::MDIntegrator> mdintegrator,
        Mode mode = ESPP_VEC_MODE_DEFAULT);
      ~Vectorization();

      void connect();
      void disconnect();

      ParticleArray& getParticleArray() { return particleArray; }
      CellNeighborList const& getNeighborList() const { return neighborList; }

      static void registerPython();

    private:
      ParticleArray particleArray;
      Mode mode;
      void resetParticles();
      void befCalcForces();
      void updatePositions();
      void updateForces();

      CellNeighborList neighborList;
      void resetCells();

      std::shared_ptr<integrator::MDIntegrator> mdintegrator;

      // signals that connect to integrator
      boost::signals2::connection sigBefCalcForces;
      boost::signals2::connection sigUpdateForces;

      // signals that connect to system
      boost::signals2::connection sigResetParticles;
      boost::signals2::connection sigResetCells;

      std::shared_ptr<storage::DomainDecomposition> decomp;

      static LOG4ESPP_DECL_LOGGER(logger);
    };

  }
}

#endif
