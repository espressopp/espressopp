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

#include <functional>

#include "System.hpp"
#include "checks.hpp"
#include "hdf5.hpp"
#include "types.hpp"

namespace espressopp
{
namespace io
{
class RestoreH5MDParallel
{
public:
    RestoreH5MDParallel(const shared_ptr<System>& system, const std::string& filename)
        : system_(system), filename_(filename)
    {
    }

    void restore();

    std::string author = "xxx";
    std::string particleGroupName = "atoms";

    bool restoreId = true;
    bool restoreType = true;
    bool restoreMass = true;
    bool restoreQ = true;
    bool restoreGhost = true;
    bool restorePosition = true;
    bool restoreVelocity = true;
    bool restoreForce = true;

    std::string idDataset = "id";
    std::string typeDataset = "type";
    std::string massDataset = "mass";
    std::string qDataset = "charge";
    std::string ghostDataset = "ghost";
    std::string positionDataset = "position";
    std::string velocityDataset = "velocity";
    std::string forceDataset = "force";

    static void registerPython();

private:
    void updateCache();

    template <typename T>
    void readParallel(hid_t fileId, const std::string& name, std::vector<T>& data);

    shared_ptr<System> system_ = nullptr;
    std::string filename_ = "";  ///<  output filename

    MPI_Comm comm = MPI_COMM_NULL;
    int rank = -1;
    int numProcesses = -1;
};

}  // namespace io
}  // namespace espressopp