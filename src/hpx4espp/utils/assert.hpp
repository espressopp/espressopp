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

#ifndef HPX4ESPP_UTILS_ASSERT_HPP
#define HPX4ESPP_UTILS_ASSERT_HPP

#include <hpx4espp/include/errors.hpp>

#define HPX4ESPP_ASSERT(VALUE)                                                        \
    {                                                                                 \
        if (!(VALUE))                                                                 \
            HPX4ESPP_THROW_EXCEPTION(                                                 \
                hpx::assertion_failure, __FUNCTION__,                                 \
                "HPX4ESPP_ASSERT FAILED. " #VALUE "=" << (VALUE ? "true" : "false")); \
    }                                                                                 \
    /* */

#define HPX4ESPP_ASSERT_EQUAL(VALUE1, VALUE2)                                                  \
    {                                                                                          \
        if (!((VALUE1) == (VALUE2)))                                                           \
            HPX4ESPP_THROW_EXCEPTION(hpx::assertion_failure, __FUNCTION__,                     \
                                     "HPX4ESPP_ASSERT_EQUAL FAILED."                           \
                                         << " " #VALUE1 "=" << (VALUE1) << " and " #VALUE2 "=" \
                                         << (VALUE2) << " are not equal.");                    \
    }                                                                                          \
    /* */

#define HPX4ESPP_ASSERT_NEQ(VALUE1, VALUE2)                                                    \
    {                                                                                          \
        if (!((VALUE1) != (VALUE2)))                                                           \
            HPX4ESPP_THROW_EXCEPTION(hpx::assertion_failure, __FUNCTION__,                     \
                                     "HPX4ESPP_ASSERT_NEQ FAILED."                             \
                                         << " " #VALUE1 "=" << (VALUE1) << " and " #VALUE2 "=" \
                                         << (VALUE2) << " are not unequal.");                  \
    }                                                                                          \
    /* */

#define HPX4ESPP_ASSERT_LT(VALUE1, VALUE2)                                                     \
    {                                                                                          \
        if (!((VALUE1) < (VALUE2)))                                                            \
            HPX4ESPP_THROW_EXCEPTION(hpx::assertion_failure, __FUNCTION__,                     \
                                     "HPX4ESPP_ASSERT_LT FAILED."                              \
                                         << " " #VALUE1 "=" << (VALUE1)                        \
                                         << " not less than " #VALUE2 "=" << (VALUE2) << "."); \
    }                                                                                          \
    /* */

#define HPX4ESPP_ASSERT_GT(VALUE1, VALUE2)                                                        \
    {                                                                                             \
        if (!((VALUE1) > (VALUE2)))                                                               \
            HPX4ESPP_THROW_EXCEPTION(hpx::assertion_failure, __FUNCTION__,                        \
                                     "HPX4ESPP_ASSERT_GT FAILED."                                 \
                                         << " " #VALUE1 "=" << (VALUE1)                           \
                                         << " not greater than " #VALUE2 "=" << (VALUE2) << "."); \
    }                                                                                             \
    /* */

#define HPX4ESPP_ASSERT_LEQ(VALUE1, VALUE2)                                                        \
    {                                                                                              \
        if (!((VALUE1) <= (VALUE2)))                                                               \
            HPX4ESPP_THROW_EXCEPTION(hpx::assertion_failure, __FUNCTION__,                         \
                                     "HPX4ESPP_ASSERT_LEQ FAILED."                                 \
                                         << " " #VALUE1 "=" << (VALUE1)                            \
                                         << " not less than nor equal to " #VALUE2 "=" << (VALUE2) \
                                         << ".");                                                  \
    }                                                                                              \
    /* */

#define HPX4ESPP_ASSERT_GEQ(VALUE1, VALUE2)                                               \
    {                                                                                     \
        if (!((VALUE1) >= (VALUE2)))                                                      \
            HPX4ESPP_THROW_EXCEPTION(hpx::assertion_failure, __FUNCTION__,                \
                                     "HPX4ESPP_ASSERT_GEQ FAILED."                        \
                                         << " " #VALUE1 "=" << (VALUE1)                   \
                                         << " not greater than nor equal to " #VALUE2 "=" \
                                         << (VALUE2) << ".");                             \
    }                                                                                     \
    /* */

#endif  // HPX4ESPP_UTILS_ASSERT_HPP
