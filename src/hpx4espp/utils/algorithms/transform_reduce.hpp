/*
  Copyright (C) 2020-2022
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

#ifndef HPX4ESPP_UTILS_ALGORITHMS_TRANSFORM_REDUCE_HPP
#define HPX4ESPP_UTILS_ALGORITHMS_TRANSFORM_REDUCE_HPP

#include "hpx4espp/include/hpx_version.hpp"

#if (HPX4ESPP_HPX_VERSION_FULL >= 10500)
#include <hpx/threading_base/thread_data.hpp>
#else
#include <hpx/runtime/threads/thread_data_fwd.hpp>
#endif
#include <hpx/parallel/algorithms/transform_reduce.hpp>

namespace espressopp
{
namespace hpx4espp
{
namespace utils
{
/// Uses HPX parallel transform reduce when running inside an HPX thread.
/// Otherwise, a regular serial implementation is used instead.
template <typename FwdIter, typename T, typename Reduce, typename Convert>
inline T parallelTransformReduce(
    FwdIter first, FwdIter last, T init, Reduce&& red_op, Convert&& conv_op)
{
    if (hpx::threads::get_self_ptr() != nullptr)
    {
        return hpx::parallel::transform_reduce(hpx::parallel::execution::par, first, last,
                                               std::move(init), std::forward<Reduce>(red_op),
                                               std::forward<Convert>(conv_op));
    }
    else
    {
        return hpx::parallel::transform_reduce(hpx::parallel::execution::seq, first, last,
                                               std::move(init), std::forward<Reduce>(red_op),
                                               std::forward<Convert>(conv_op));
    }
}

}  // namespace utils
}  // namespace hpx4espp
}  // namespace espressopp

#endif  // HPX4ESPP_UTILS_ALGORITHMS_TRANSFORM_REDUCE_HPP
