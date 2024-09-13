/*
  Copyright (C) 2020-2022-2022
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

#include "hpx4espp/HPXRuntime.hpp"
#include "hpx4espp/include/logging.hpp"
#include "hpx4espp/include/hpx_version.hpp"

#include <hpx/hpx_start.hpp>
#include <hpx/runtime_distributed.hpp>
#include <hpx/modules/runtime_local.hpp>

#if (HPX4ESPP_HPX_VERSION_FULL >= 10701)
#include <hpx/modules/agas.hpp>
#else
#include <hpx/runtime/agas/addressing_service.hpp>
#endif

#include "python.hpp"

int __argc = 0;
char** __argv = nullptr;

void set_argc_argv(int argc, char* argv[], char* env[])
{
    __argc = argc;
    __argv = argv;
}

/** parse argc and argv at program startup */
__attribute__((section(".init_array"))) void (*set_global_argc_argv)(int,
                                                                     char*[],
                                                                     char*[]) = &set_argc_argv;

espressopp::hpx4espp::detail::manage_hpx_runtime* rts = nullptr;

namespace espressopp
{
namespace hpx4espp
{
namespace detail
{
/////////////////////////////////////////////////////////////////////////////////////////////

manage_hpx_runtime::manage_hpx_runtime(bool disable_tcp, size_t threads)
    : running_(false), rts_(nullptr)
{
    if (__argv == nullptr && __argc == 0)
    {
        std::cerr << "argc and argv were not parsed correctly." << std::endl;
        std::abort();
    }

    std::vector<std::string> cfg = {
        // make sure hpx_main is always executed regardless of locality mode
        "hpx.run_hpx_main!=1",
        // allow for unknown command line options
        "hpx.commandline.allow_unknown!=1",
        // disable HPX' short options
        "hpx.commandline.aliasing!=0",
    };

    if (disable_tcp)
    {
        // disable the TCP parcelport and exclusively use MPI
        cfg.push_back("hpx.parcel.tcp.enable!=0");
    }

    if (threads > 0)
    {
        std::string cfg_threads = "hpx.os_threads!=" + std::to_string(threads);
        cfg.push_back(cfg_threads);
    }

    using hpx::placeholders::_1;
    using hpx::placeholders::_2;
    hpx::util::function_nonser<int(int, char**)> start_function =
        hpx::bind(&manage_hpx_runtime::hpx_main, this, hpx::placeholders::_1, hpx::placeholders::_2);

    hpx::init_params init_args;
    init_args.cfg = cfg;
    init_args.mode = hpx::runtime_mode::console;

    HPX4ESPP_DEBUG_MSG("Starting HPX runtime")
    if (!hpx::start(start_function, __argc, __argv, init_args))
    {
        HPX4ESPP_DEBUG_MSG("Something went wrong");
        // Something went wrong while initializing the runtime.
        // This early we can't generate any output, just bail out.
        std::abort();
    }
    HPX4ESPP_DEBUG_MSG("HPX runtime started")
    // Wait for the main HPX thread (hpx_main below) to have started running
    std::unique_lock<std::mutex> lk(startup_mtx_);
    while (!running_) startup_cond_.wait(lk);
}

manage_hpx_runtime::~manage_hpx_runtime()
{
    // notify hpx_main to tear down the runtime
    {
        std::lock_guard<hpx::lcos::local::spinlock> lk(mtx_);
        rts_ = nullptr;  // reset pointer
    }
    cond_.notify_one();  // signal exit

    // wait for the runtime to exit
    hpx::stop();
    HPX4ESPP_DEBUG_MSG("HPX runtime terminated")
}

// registration of external (to HPX) threads
void manage_hpx_runtime::register_thread(char const* name) { hpx::register_thread(rts_, name); }

void manage_hpx_runtime::unregister_thread() { hpx::unregister_thread(rts_); }

int manage_hpx_runtime::hpx_main(int argc, char** argv)
{
    // Store a pointer to the runtime here.
    rts_ = hpx::get_runtime_ptr();

    // Signal to constructor that thread has started running.
    {
        std::lock_guard<std::mutex> lk(startup_mtx_);
        running_ = true;
    }
    startup_cond_.notify_one();

    // Now, wait for destructor to be called.
    {
        std::unique_lock<hpx::lcos::local::spinlock> lk(mtx_);
        if (rts_ != nullptr) cond_.wait(lk);
    }

#if defined(HPX_HAVE_NETWORKING)
    if (hpx::get_locality_id() == 0)
    {
        return hpx::finalize();
    }
    else
    {
        return 0;
    }
#else
    return hpx::finalize();
#endif
}
}  // namespace detail

LOG4ESPP_LOGGER(HPXRuntime::logger, "HPXRuntime");

HPXRuntime::HPXRuntime() {}

HPXRuntime::~HPXRuntime() {}

void HPXRuntime::start(bool disable_tcp, size_t threads)
{
    if (rts == nullptr) rts = new detail::manage_hpx_runtime(disable_tcp, threads);
}

void HPXRuntime::stop()
{
    detail::manage_hpx_runtime* r = rts;
    rts = nullptr;
    if (r != nullptr) delete r;
}

bool HPXRuntime::networkingEnabled() const { return hpx::is_networking_enabled(); }

size_t HPXRuntime::getNumOverallThreads() const
{
    if (isRunning())
    {
        return hpx::agas::get_num_overall_threads(hpx::launch::sync);
    }
    else
    {
        int size = 0;
        MPI_Comm_size(MPI_COMM_WORLD, &size);
        return size_t(size);
    }
}

size_t HPXRuntime::getNumThreads() const
{
    if (isRunning())
    {
        return hpx::get_os_thread_count();
    }
    else
    {
        return size_t(1);
    }
}

bool HPXRuntime::isRunning() { return rts != nullptr; }

std::string HPXRuntime::parcelport() const
{
#ifdef HPX_HAVE_NETWORKING
    if (rts != nullptr && hpx::is_networking_enabled())
        return hpx::get_runtime_distributed()
            .get_parcel_handler()
            .get_bootstrap_parcelport()
            ->type();
#endif
    return "none";
}

void HPXRuntime::registerPython()
{
    using namespace espressopp::python;

    class_<HPXRuntime, std::shared_ptr<HPXRuntime> >("HPXRuntime", init<>())
        .def("start", &HPXRuntime::start)
        .def("stop", &HPXRuntime::stop)
        .def("getNumOverallThreads", &HPXRuntime::getNumOverallThreads)
        .def("getNumThreads", &HPXRuntime::getNumThreads)
        .def("networkingEnabled", &HPXRuntime::networkingEnabled)
        .def("parcelport", &HPXRuntime::parcelport);

    def("hpx4espp_isRunning", &HPXRuntime::isRunning);
}
}  // namespace hpx4espp
}  // namespace espressopp
