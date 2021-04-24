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

// ESPP_CLASS
#ifndef VEC_PARTICLEARRAY_HPP
#define VEC_PARTICLEARRAY_HPP

#include "vec/include/simdconfig.hpp"
#include "types.hpp"
#include "Particle.hpp"

namespace espressopp { namespace vec {

  enum ParticleElements {
    PARTICLE_PROPERTIES=1,
    PARTICLE_POSITION=2,
    PARTICLE_MOMENTUM=4,
    PARTICLE_LOCAL=8,
    PARTICLE_FORCE=16,
    PARTICLE_ALL=31,
    PARTICLE_POSITION_ONLY=32,
    PARTICLE_VELOCITY_ONLY=64,
    PARTICLE_FORCE_ONLY=128
  };

  /*
      Optimizations:
        - store particle data separately as arrays of individual attributes (SOA)
        - use aligned vector everywhere
        - pad ends of vector w/ pseudo particles to keep cell length a multiple of the cpu vector width
   */
  class ParticleArray
  {
  public:
    ParticleArray();

    void markRealCells(CellList const& realcells, const Cell* cell0, size_t numLocalCells);
    void markRealCells(std::vector<size_t> const& realcells, size_t numLocalCells);

    void copyFrom(CellList const& srcCells);

    void updateFromPositionVelocity(CellList const& srcCells, bool realOnly);
    void updateToPositionVelocity(CellList & srcCells, bool realOnly) const;

    void updateFromPositionOnly(CellList const& srcCells);
    void updateToPositionOnly(CellList & srcCells) const;
    void updateFromForceOnly(CellList const& srcCells, bool realOnly);
    void updateToForceOnly(CellList & srcCells, bool realOnly) const;

    void addToForceOnly(CellList & srcCells) const;

    inline size_t size() const { return size_; }
    inline size_t numCells() const { return sizes_.size(); }
    inline std::vector<size_t> const& cellRange() const { return cellRange_; }
    inline std::vector<size_t> const& sizes() const { return sizes_; }
    inline std::vector<size_t> const& realCells() const { return realCells_; }
    inline std::vector<size_t> const& ghostCells() const { return ghostCells_; }
    bool checkSizes() const;
    void verify(CellList const& srcCells) const;

    AlignedVector< size_t > id;
    AlignedVector< real > mass;
    AlignedVector< real > q;
    AlignedVector< bool > ghost;

    /* SOA */

    AlignedVector< real > p_x;
    AlignedVector< real > p_y;
    AlignedVector< real > p_z;
    AlignedVector< lint > type;

    AlignedVector< real > v_x;
    AlignedVector< real > v_y;
    AlignedVector< real > v_z;

    AlignedVector< real > f_x;
    AlignedVector< real > f_y;
    AlignedVector< real > f_z;

  protected:
    /// start=cellRange_[i] to end=cellRange_[i+1] for cell[i] including padding
    std::vector<size_t> cellRange_;
    /// actual size of particle array
    std::vector<size_t> sizes_;
    /// indices of real cells
    std::vector<size_t> realCells_;
    /// indices of ghost cells
    std::vector<size_t> ghostCells_;
    /// number of local cells in source storage
    std::size_t numLocalCells_ = 0;
    std::size_t size_ = 0;
    std::size_t data_size_ = 0; // including padding
    std::size_t reserve_size_ = 0; // actual size of vectors (allows for re-use of arrays)
    static const std::size_t chunk_size_ = ESPP_VECTOR_WIDTH;
    inline std::size_t calc_data_size(size_t const& size) { return (1 + ((size - 1) / chunk_size_)) * chunk_size_; }

    void updateFrom(std::vector<Particle> const& particlelist, size_t start);
    void markGhostCells();

  public:

    inline Real3D getPosition(size_t i) const
    {
      return Real3D(p_x[i],p_y[i],p_z[i]);
    }

    inline lint getType(size_t i) const
    {
      return type[i];
    }

    inline void addForce(size_t i, Real3D const& ff)
    {
      f_x[i] += ff[0];
      f_y[i] += ff[1];
      f_z[i] += ff[2];
    }

    inline void subForce(size_t i, Real3D const& ff)
    {
      f_x[i] -= ff[0];
      f_y[i] -= ff[1];
      f_z[i] -= ff[2];
    }

  };
}}

#endif//VEC_PARTICLEARRAY_HPP
