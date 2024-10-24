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

#ifndef HPX4ESPP_INCLUDE_LOGGING_HPP
#define HPX4ESPP_INCLUDE_LOGGING_HPP

/// Some debugging messages for code dev
#ifdef HPX4ESPP_DEBUG

#include <string.h>
#include <iostream>
#include "mpi.h"

#define HPX4ESPP_DEBUG_MSG(MSG)                                                              \
    {                                                                                        \
        const char* _fn_ = (strrchr(__FILE__, '/') ? strrchr(__FILE__, '/') + 1 : __FILE__); \
        int _rank_, _ec_;                                                                    \
        _ec_ = MPI_Comm_rank(MPI_COMM_WORLD, &_rank_);                                       \
        if (_ec_) MPI_Abort(MPI_COMM_WORLD, _ec_);                                           \
        std::cout << "[" << _rank_ << "] [" << _fn_ << ":" << __LINE__ << "] " << MSG        \
                  << std::endl;                                                              \
    }                                                                                        \
    /*  */

// #include <hpx/include/iostreams.hpp>
// #include <sstream>

// #define HPX4ESPP_DEBUG_MSG_THREAD(MSG) {                                                        \
  //   const char* filename_ = (strrchr(__FILE__, '/') ? strrchr(__FILE__, '/') + 1 : __FILE__);     \
  //   std::stringstream ss_; ss_ << MSG;                                                            \
  //   hpx::util::format_to(hpx::cout, "[{1}:{2}] {3}",filename_,__LINE__,ss_.str()) << hpx::flush;  \
  // }                                                                                               \
  // /*  */

#else
#define HPX4ESPP_DEBUG_MSG(MSG) \
    {                           \
    }
#endif

#endif  // HPX4ESPP_INCLUDE_LOGGING_HPP
