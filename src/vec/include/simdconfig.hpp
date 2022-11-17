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

#ifndef VEC_INCLUDE_SIMDCONFIG_HPP
#define VEC_INCLUDE_SIMDCONFIG_HPP

#include <vector>
#include <cstdint>
#include <boost/align/aligned_allocator.hpp>
#include "include/esconfig.hpp"

#ifdef __INTEL_COMPILER
#define ESPP_VEC_PRAGMAS _Pragma("ivdep")
#else
#define ESPP_VEC_PRAGMAS _Pragma("GCC ivdep")
#endif
// TODO: Test and define for more compilers

namespace espressopp
{
namespace vec
{
/// Padding and alignment to 64-byte boundaries
constexpr size_t ESPP_VECTOR_ALIGNMENT = 64;

/// Number of espressopp::real elements in the SIMD width
constexpr size_t ESPP_VECTOR_WIDTH = ESPP_VECTOR_ALIGNMENT / sizeof(espressopp::real);

constexpr size_t ESPP_FIT_TO_VECTOR_WIDTH(size_t SIZE)
{
    return ((((SIZE) + ESPP_VECTOR_WIDTH - 1) / ESPP_VECTOR_WIDTH) * ESPP_VECTOR_WIDTH);
}

/// Vector allocated aligned to given Alignment (64 bytes by default)
template <typename T, std::size_t Alignment = ESPP_VECTOR_ALIGNMENT>
using AlignedVector = std::vector<T, boost::alignment::aligned_allocator<T, Alignment>>;

static_assert(std::is_same<real, double>::value || std::is_same<real, float>::value,
              "Only float and double are allowed for espressopp::real. Otherwise, manually choose "
              "a value for large_pos.");

/// Represents a very large number for padding positions of "fake" particles = sqrt(max/3).
constexpr real large_pos = std::is_same<real, double>::value ? 7.74099e150 : 7.74099e15;

}  // namespace vec
}  // namespace espressopp

#endif  // VEC_INCLUDE_SIMDCONFIG_HPP
