/*
  Copyright (C) 2021
      Sebastian Eibl, Max Planck Computing & Data Facility

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

#pragma once

#include "checks.hpp"
#include <hdf5.h>
#include <hdf5_hl.h>

namespace espressopp
{
namespace io
{
template <typename T>
hid_t typeToHDF5()
{
    CHECK_TRUE(false, "This type is not supported!");
}
template <>
inline hid_t typeToHDF5<int8_t>()
{
    return H5T_NATIVE_INT8;
}
template <>
inline hid_t typeToHDF5<int16_t>()
{
    return H5T_NATIVE_INT16;
}
template <>
inline hid_t typeToHDF5<int32_t>()
{
    return H5T_NATIVE_INT32;
}
template <>
inline hid_t typeToHDF5<int64_t>()
{
    return H5T_NATIVE_INT64;
}
template <>
inline hid_t typeToHDF5<double>()
{
    return H5T_NATIVE_DOUBLE;
}

template <typename T>
auto CHECK_HDF5(const T& status)
{
    if (status < 0)
    {
        H5Eprint(H5E_DEFAULT, stderr);
        CHECK_GREATER_EQUAL(status, 0);
    }
    return status;
}

}  // namespace io
}  // namespace espressopp