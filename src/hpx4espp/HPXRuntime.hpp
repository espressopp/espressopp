/*
  Copyright (C) 2020-2022
      Max Planck Institute for Polymer Research & JGU Mainz
  Copyright (C) 2018
      Hartmut Kaiser (for manage_hpx_runtime)

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

#ifndef HPX4ESPP_HPXRUNTIME_HPP
#define HPX4ESPP_HPXRUNTIME_HPP

#include "log4espp.hpp"

#include <hpx/runtime.hpp>
#include <hpx/synchronization/spinlock.hpp>
#include <hpx/synchronization/condition_variable.hpp>

#include <condition_variable>

namespace espressopp
{
namespace hpx4espp
{
namespace detail
{
/////////////////////////////////////////////////////////////////////////////
/**
    Initializes a console instance of HPX.

    Taken from:
      - phylanx/python/src/init_hpx.cpp
      - hpx/examples/quickstart/init_globally.cpp
*/
struct manage_hpx_runtime
{
    /** Initializes the HPX runtime by calling \a hpx::start */
    manage_hpx_runtime(bool disable_tcp, size_t threads);

    /** Terminates the HPX runtime system by calling \a hpx::stop */
    ~manage_hpx_runtime();

    /** Register an OS thread as an HPX thread */
    void register_thread(char const* name);

    /** Unregisters an HPX thread */
    void unregister_thread();

protected:
    /** Builds the list of all localities and upon termination, calls
        \a hpx::finalize from the console locality */
    int hpx_main(int argc, char** argv);

private:
    hpx::spinlock mtx_;
    hpx::condition_variable_any cond_;
    std::mutex startup_mtx_;
    std::condition_variable startup_cond_;
    bool running_;
    hpx::runtime* rts_;
};
}  // namespace detail

/** Wrapper that exposes the HPX runtime to python */
class HPXRuntime
{
public:
    HPXRuntime();
    ~HPXRuntime();
    void start(bool disable_tcp, size_t threads);
    void stop();
    bool networkingEnabled() const;
    size_t getNumOverallThreads() const;
    size_t getNumThreads() const;
    static bool isRunning();
    std::string parcelport() const;
    static void registerPython();

private:
    static LOG4ESPP_DECL_LOGGER(logger);
};
}  // namespace hpx4espp
}  // namespace espressopp

#endif  // HPX4ESPP_HPXRUNTIME_HPP
