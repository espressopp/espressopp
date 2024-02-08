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

#include "python.hpp"
#include "HPXUtils.hpp"

#include <stdlib.h>
#include <stdio.h>
#include <string.h>

namespace espressopp
{
namespace hpx4espp
{
namespace utils
{
namespace HPXUtils
{
bool enabled()
{
#if HPX4ESPP_ENABLED
    return true;
#else
    return false;
#endif
}

bool have_networking()
{
#ifdef HPX_HAVE_NETWORKING
    return true;
#else
    return false;
#endif
}

int get_cpu_id()
{
#ifdef __linux__
    {
        /// Adapted from
        /// https://hpcf.umbc.edu/general-productivity/checking-which-cpus-are-used-by-your-program/
        FILE* procfile = fopen("/proc/self/stat", "r");
        long to_read = 8192;
        char buffer[to_read];
        int read = fread(buffer, sizeof(char), to_read, procfile);
        fclose(procfile);

        // Field with index 38 (zero-based counting) is the one we want
        char* line = strtok(buffer, " ");
        for (int i = 1; i < 38; i++)
        {
            line = strtok(NULL, " ");
        }
        line = strtok(NULL, " ");
        int cpu_id = atoi(line);
        return cpu_id;
    }
#else
    {
        std::cerr << "WARNING: get_cpu_id implemented only for linux" << std::endl;
        return -1;
    }
#endif
}

void registerPython()
{
    using namespace espressopp::python;
    def("hpx4espp_enabled", enabled);
    def("hpx4espp_have_networking", have_networking);
    def("hpx4espp_get_cpu_id", get_cpu_id);
}
}  // namespace HPXUtils
}  // namespace utils
}  // namespace hpx4espp
}  // namespace espressopp
