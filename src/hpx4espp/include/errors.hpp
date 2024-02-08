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

#ifndef HPX4ESPP_INCLUDE_ERRORS_HPP
#define HPX4ESPP_INCLUDE_ERRORS_HPP

#include <hpx/errors/error.hpp>  // list of error codes
#include <hpx/errors/throw_exception.hpp>
#include <string.h>
#include <sstream>
#include "mpi.h"

#define HPX4ESPP_THROW_EXCEPTION(errcode, f, msg)                                                 \
    {                                                                                             \
        const char* filename_ = (strrchr(__FILE__, '/') ? strrchr(__FILE__, '/') + 1 : __FILE__); \
        int _rank_, _ec_;                                                                         \
        _ec_ = MPI_Comm_rank(MPI_COMM_WORLD, &_rank_);                                            \
        if (_ec_) MPI_Abort(MPI_COMM_WORLD, _ec_);                                                \
        std::stringstream ss_;                                                                    \
        ss_ << "MPI rank " << _rank_ << ", file \"" << __FILE__ << "\", line " << __LINE__        \
            << ", in " << f << "\n  " << msg;                                                     \
        std::cerr << "\nHPX4ESPP_THROW_EXCEPTION: " << ss_.str() << ": HPX("                      \
                  << hpx::get_error_name(errcode) << ")\n"                                        \
                  << std::endl;                                                                   \
        HPX_THROW_EXCEPTION(errcode, f, ss_.str());                                               \
    }                                                                                             \
    /*  */

#define HPX4ESPP_NOT_IMPLEMENTED(msg) \
    HPX4ESPP_THROW_EXCEPTION(hpx::not_implemented, __FUNCTION__, msg)

#endif  // HPX4ESPP_INCLUDE_ERRORS_HPP
