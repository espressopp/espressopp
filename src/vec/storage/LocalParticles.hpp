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

#ifndef VEC_STORAGE_LOCALPARTICLES_HPP
#define VEC_STORAGE_LOCALPARTICLES_HPP

#include "vec/include/types.hpp"
#include "vec/include/simdconfig.hpp"

#include <unordered_map>

namespace espressopp { namespace vec {
  namespace storage {

    typedef std::unordered_map<lint, size_t> LocalParticlesBase;

    ///////////////////////////////////////////////////////////////////////////
    /// Provides a lookup from particle id to (cellIdx, particleIdx)
    /// where particleIdx is the relative position in the cell list.
    /// This version resets everytime new particles are loaded so it will only
    /// store the index in ParticleArray.
    struct LocalParticles :
      public LocalParticlesBase
    {
    public:
      typedef LocalParticlesBase base;

      /// Default constructor
      LocalParticles(){}

      /// Rebuild using all particles in uniqueCells, which should prevent any particles to be
      /// inserted multiple times
      void rebuild(ParticleArray const& pa, std::vector<size_t> const& uniqueCells);
    };

  }
}}

#endif//VEC_STORAGE_LOCALPARTICLES_HPP
