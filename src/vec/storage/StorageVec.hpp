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

#ifndef VEC_STORAGEVEC_HPP
#define VEC_STORAGEVEC_HPP

#include "vec/include/types.hpp"
#include "vec/ParticleArray.hpp"
#include "vec/Vectorization.hpp"
#include "LocalParticles.hpp"

#include "SystemAccess.hpp"
#include "python.hpp"
#include "log4espp.hpp"

#include <boost/signals2.hpp>

#define VEC_PARTICLE_NOT_FOUND (std::numeric_limits<size_t>::max())

namespace espressopp { namespace vec {
  namespace storage {

    class StorageVec
      : protected SystemAccess
    {
    public:
      StorageVec(shared_ptr<System> system);

      virtual void loadCells() = 0;

      virtual void unloadCells() = 0;

      virtual void updateGhostsVec() = 0;

      virtual void collectGhostForcesVec() = 0;

      boost::signals2::signal<void ()> onLoadCells;

      boost::signals2::signal<void ()> onUnloadCells;

      static void registerPython();

    protected:

      LocalParticles localParticlesVec;
      std::vector<size_t> uniqueCells;

    public:
      inline size_t lookupLocalParticleVec(size_t id)
      {
        const auto it = localParticlesVec.find(id);
        return (it != localParticlesVec.end()) ?
          it->second : VEC_PARTICLE_NOT_FOUND;
      }

      inline size_t lookupRealParticleVec(size_t id)
      {
        const auto it = localParticlesVec.find(id);

        size_t ret = VEC_PARTICLE_NOT_FOUND;
        if(it!=localParticlesVec.end()) {
          if(!(getSystem()->vectorization->particles.ghost[it->second])) {
            ret = it->second;
          }
        }
        return ret;
      }

    private:
      static LOG4ESPP_DECL_LOGGER(logger);
    };

  }
}}

#endif//VEC_STORAGEVEC_HPP
