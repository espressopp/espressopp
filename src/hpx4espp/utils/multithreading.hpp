/*
  Copyright (C) 2020
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

#ifndef HPX4ESPP_UTILS_MULTITHREADING_HPP
#define HPX4ESPP_UTILS_MULTITHREADING_HPP

#include "hpx4espp/include/hpx_version.hpp"

#if (HPX4ESPP_HPX_VERSION_FULL >= 10500)
#include <hpx/threading_base/thread_data.hpp>
#include <hpx/include/run_as.hpp>
#else
#include <hpx/runtime/threads/thread_data_fwd.hpp>
#include <hpx/runtime/threads/run_as_hpx_thread.hpp>
#endif

#include <hpx/include/parallel_for_loop.hpp>

namespace espressopp
{
namespace hpx4espp
{
namespace utils
{
bool rtsIsRunning();

size_t getNumThreads();

inline bool isOutsideHPXThread() { return hpx::threads::get_self_ptr() == nullptr; }

inline bool isInsideHPXThread() { return hpx::threads::get_self_ptr() != nullptr; }

/// Runs a function inside an HPX thread when the runtime system is active
/// Otherwise, the function is called directly on the main thread
template <typename F, typename... Ts>
inline typename hpx::util::invoke_result<F, Ts...>::type runAsHPXThread(F const& f, Ts&&... vs)
{
    if (rtsIsRunning())
    {
        // std::cout << "NOTE: runAsHPXThread called and RTS is running" << std::endl;
        return hpx::threads::run_as_hpx_thread(f, std::forward<Ts>(vs)...);
    }
    else
    {
        // std::cout << "WARNING: runAsHPXThread called but RTS is not running" << std::endl;
        return f(std::forward<Ts>(vs)...);
    }
}

}  // namespace utils
}  // namespace hpx4espp
}  // namespace espressopp

/// Requires #include <hpx4espp/include/error.hpp>
#define HPX4ESPP_ASSERT_IS_INSIDE_HPX_THREAD(FUNCTION)                                 \
    {                                                                                  \
        if (!espressopp::hpx4espp::utils::isInsideHPXThread())                         \
            HPX4ESPP_THROW_EXCEPTION(                                                  \
                hpx::assertion_failure, __FUNCTION__,                                  \
                "HPX4ESPP_CHECK_IS_INSIDE_HPX_THREAD: Not running inside HPX thread"); \
    }                                                                                  \
    /* */

#endif  // HPX4ESPP_UTILS_MULTITHREADING_HPP
