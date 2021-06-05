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

// ESPP_CLASS
#ifndef _PARTICLEARRAY_HPP
#define _PARTICLEARRAY_HPP

#include "types.hpp"
#include "simdconfig.hpp"
#include "Particle.hpp"

#define ESPP_AOS
#undef ESPP_AOS  // USE SOA

#ifndef ESPP_SOA
#define ESPP_SOA
#endif

#ifdef ESPP_AOS
#undef ESPP_SOA
#endif

namespace espressopp
{
namespace vectorization
{
typedef std::uint64_t ulongint;

enum ParticleElements
{
    PARTICLE_PROPERTIES = 1,
    PARTICLE_POSITION = 2,
    PARTICLE_MOMENTUM = 4,
    PARTICLE_LOCAL = 8,
    PARTICLE_FORCE = 16,
    PARTICLE_ALL = 31,
    PARTICLE_POSITION_ONLY = 32,
    PARTICLE_VELOCITY_ONLY = 64,
    PARTICLE_FORCE_ONLY = 128
};

enum Mode
{
    ESPP_VEC_SOA = 0,
    ESPP_VEC_AOS = 1
};

#if defined(ESPP_AOS)
#define ESPP_VEC_MODE_DEFAULT ESPP_VEC_AOS
#else
#define ESPP_VEC_MODE_DEFAULT ESPP_VEC_SOA
#endif

struct Real3DInt
{
    real x;
    real y;
    real z;
    ulongint t;

    Real3DInt()
    {
        // for padding to work, real must be 64 bit floating point type
        static_assert(sizeof(real) == sizeof(ulongint),
                      "Size of real not equal to ulongint. Check esconfig.hpp for typedefs.");
    }

    Real3DInt(real x, real y, real z, ulongint t) : x(x), y(y), z(z), t(t) {}

    Real3DInt(Real3D const& r, ulongint t) : x(r[0]), y(r[1]), z(r[2]), t(t) {}

    Real3D const& to_Real3D() const { return *((Real3D*)(this)); }
    void operator=(Real3D const& src)
    {
        x = src[0];
        y = src[1];
        z = src[2];
    }
};

struct Real4D
{
    real x;
    real y;
    real z;
    real w;

    Real4D() {}

    Real4D(real x, real y, real z, real w) : x(x), y(y), z(z), w(w) {}

    Real4D(Real3D const& r, real w) : x(r[0]), y(r[1]), z(r[2]), w(w) {}

    Real3D const& to_Real3D() const { return *((Real3D*)(this)); }
    void operator=(Real3D const& src)
    {
        x = src[0];
        y = src[1];
        z = src[2];
    }
};

/*
    Optimizations:
      - store particle data separately as arrays of individual attributes (SOA)
      - store position and type in packed format (x,y,z,type)
      - use aligned vector everywhere
      - pad ends of vector w/ fake particles to keep cell length a multiple of the cpu vector width
 */
class ParticleArray
{
public:
    void copyFrom(CellList const& srcCells, Mode mode = ESPP_VEC_MODE_DEFAULT);

    void updateFromPositionOnly(CellList const& srcCells);
    void addToForceOnly(CellList& srcCells) const;

    std::vector<size_t> const& cellRange() const { return cellRange_; }
    std::vector<size_t> const& sizes() const { return sizes_; }
    bool checkSizes();
    inline bool mode_aos() const { return mode == ESPP_VEC_AOS; }

    AlignedVector<Real3DInt> position; /* p_x,p_y,p_z,p_type */
    AlignedVector<Real4D> force;       /* f_x,f_y,f_z,padding */

    AlignedVector<real> p_x;
    AlignedVector<real> p_y;
    AlignedVector<real> p_z;
    AlignedVector<real> f_x;
    AlignedVector<real> f_y;
    AlignedVector<real> f_z;
    AlignedVector<ulongint> type;

protected:
    /// start=cellRange_[i] to end=cellRange_[i+1] for cell[i] including padding
    std::vector<size_t> cellRange_;
    /// actual size of particle array
    std::vector<size_t> sizes_;

    Mode mode;

    std::size_t size_ = 0;
    std::size_t data_size_ = 0;     // including padding
    std::size_t reserve_size_ = 0;  // actual size of vectors (allows for re-use of arrays)
    static const std::size_t chunk_size_ = ESPP_VECTOR_WIDTH;
    inline std::size_t calc_data_size(size_t const& size)
    {
        return (1 + ((size - 1) / chunk_size_)) * chunk_size_;
    }

    void updateFrom(std::vector<Particle> const& particlelist, size_t start);
};
}  // namespace vectorization
}  // namespace espressopp

#define ESPP_PARTICLEARRAY_AOS_APPLY(COMMAND) \
    position.COMMAND;                         \
    force.COMMAND;                            \
    /* */

#define ESPP_PARTICLEARRAY_SOA_APPLY(COMMAND) \
    p_x.COMMAND;                              \
    p_y.COMMAND;                              \
    p_z.COMMAND;                              \
    f_x.COMMAND;                              \
    f_y.COMMAND;                              \
    f_z.COMMAND;                              \
    type.COMMAND;                             \
    /* */

#endif  // _PARTICLEARRAY_HPP
