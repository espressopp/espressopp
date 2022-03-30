/*
  Copyright (C) 2021-2022
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

#ifndef HPX4ESPP_UTILS_ASSERTMSG_HPP
#define HPX4ESPP_UTILS_ASSERTMSG_HPP

#include <hpx4espp/include/errors.hpp>

#define HPX4ESPP_ASSERT_MSG(VALUE, MSG)                                                 \
    {                                                                                   \
        if (!(VALUE))                                                                   \
            HPX4ESPP_THROW_EXCEPTION(hpx::assertion_failure, __FUNCTION__,              \
                                     "HPX4ESPP_ASSERT FAILED. " #VALUE "="              \
                                         << (VALUE ? "true" : "false") << ". " << MSG); \
    }                                                                                   \
    /* */

#define HPX4ESPP_ASSERT_EQUAL_MSG(VALUE1, VALUE2, MSG)                                         \
    {                                                                                          \
        if (!((VALUE1) == (VALUE2)))                                                           \
            HPX4ESPP_THROW_EXCEPTION(hpx::assertion_failure, __FUNCTION__,                     \
                                     "HPX4ESPP_ASSERT_EQUAL FAILED."                           \
                                         << " " #VALUE1 "=" << (VALUE1) << " and " #VALUE2 "=" \
                                         << (VALUE2) << " are not equal. " << MSG);            \
    }                                                                                          \
    /* */

#define HPX4ESPP_ASSERT_NEQ_MSG(VALUE1, VALUE2, MSG)                                           \
    {                                                                                          \
        if (!((VALUE1) != (VALUE2)))                                                           \
            HPX4ESPP_THROW_EXCEPTION(hpx::assertion_failure, __FUNCTION__,                     \
                                     "HPX4ESPP_ASSERT_NEQ FAILED."                             \
                                         << " " #VALUE1 "=" << (VALUE1) << " and " #VALUE2 "=" \
                                         << (VALUE2) << " are not unequal. " << MSG);          \
    }                                                                                          \
    /* */

#define HPX4ESPP_ASSERT_LT_MSG(VALUE1, VALUE2, MSG)                                           \
    {                                                                                         \
        if (!((VALUE1) < (VALUE2)))                                                           \
            HPX4ESPP_THROW_EXCEPTION(hpx::assertion_failure, __FUNCTION__,                    \
                                     "HPX4ESPP_ASSERT_LT FAILED."                             \
                                         << " " #VALUE1 "=" << (VALUE1)                       \
                                         << " not less than " #VALUE2 "=" << (VALUE2) << ". " \
                                         << MSG);                                             \
    }                                                                                         \
    /* */

#define HPX4ESPP_ASSERT_GT_MSG(VALUE1, VALUE2, MSG)                                              \
    {                                                                                            \
        if (!((VALUE1) > (VALUE2)))                                                              \
            HPX4ESPP_THROW_EXCEPTION(hpx::assertion_failure, __FUNCTION__,                       \
                                     "HPX4ESPP_ASSERT_GT FAILED."                                \
                                         << " " #VALUE1 "=" << (VALUE1)                          \
                                         << " not greater than " #VALUE2 "=" << (VALUE2) << ". " \
                                         << MSG);                                                \
    }                                                                                            \
    /* */

#define HPX4ESPP_ASSERT_LEQ_MSG(VALUE1, VALUE2, MSG)                                               \
    {                                                                                              \
        if (!((VALUE1) <= (VALUE2)))                                                               \
            HPX4ESPP_THROW_EXCEPTION(hpx::assertion_failure, __FUNCTION__,                         \
                                     "HPX4ESPP_ASSERT_LEQ FAILED."                                 \
                                         << " " #VALUE1 "=" << (VALUE1)                            \
                                         << " not less than nor equal to " #VALUE2 "=" << (VALUE2) \
                                         << ". " << MSG);                                          \
    }                                                                                              \
    /* */

#define HPX4ESPP_ASSERT_GEQ_MSG(VALUE1, VALUE2, MSG)                                      \
    {                                                                                     \
        if (!((VALUE1) >= (VALUE2)))                                                      \
            HPX4ESPP_THROW_EXCEPTION(hpx::assertion_failure, __FUNCTION__,                \
                                     "HPX4ESPP_ASSERT_GEQ FAILED."                        \
                                         << " " #VALUE1 "=" << (VALUE1)                   \
                                         << " not greater than nor equal to " #VALUE2 "=" \
                                         << (VALUE2) << ". " << MSG);                     \
    }                                                                                     \
    /* */

#endif  // HPX4ESPP_UTILS_ASSERTMSG_HPP
