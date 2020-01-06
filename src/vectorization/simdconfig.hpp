/*
  Copyright (C) 2019
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

#ifndef _VECTORIZATION_SIMDCONFIG_HPP
#define _VECTORIZATION_SIMDCONFIG_HPP

// ESPP_VECTOR_WIDTH assumes real = double and ulongint = uint64_t
#if defined(__AVX512F__)
 #if __AVX512F__
  #define ESPP_VECTOR_ALIGNMENT 64
  #define ESPP_VECTOR_WIDTH      8
  #define ESPP_VECTOR_MASK
 #endif
#elif defined(__AVX2__)
 #if __AVX2__
  #define ESPP_VECTOR_ALIGNMENT 32
  #define ESPP_VECTOR_WIDTH      4
 #endif
#else
 #define ESPP_VECTOR_ALIGNMENT 16
 #define ESPP_VECTOR_WIDTH      2
#endif

#define ESPP_FIT_TO_VECTOR_WIDTH(SIZE) ((((SIZE)+ESPP_VECTOR_WIDTH-1)/ESPP_VECTOR_WIDTH)*ESPP_VECTOR_WIDTH)

#include <vector>
#include <boost/align/aligned_allocator.hpp>

namespace espressopp {
  namespace vectorization {

    // aligned
    template <typename T, std::size_t Alignment=ESPP_VECTOR_ALIGNMENT>
    using AlignedVector = std::vector<T, boost::alignment::aligned_allocator<T,Alignment>>;

    // represents a very large number for padding positions of "fake" particles = sqrt(max/3)
    // compatible only with real = double
    static const real large_pos = 7.74099e150;

  }
}

#endif
