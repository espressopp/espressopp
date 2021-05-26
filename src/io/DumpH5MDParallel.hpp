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
class DumpH5MDParallel
{
public:
    DumpH5MDParallel(const shared_ptr<System>& system, const std::string& filename)
        : system_(system), filename_(filename)
    {
    }

    void dump();

    std::string author = "xxx";
    std::string particleGroupName = "atoms";

    bool dumpId = true;
    bool dumpType = true;
    bool dumpMass = true;
    bool dumpQ = true;
    bool dumpGhost = true;
    bool dumpPosition = true;
    bool dumpVelocity = true;
    bool dumpForce = true;

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

    void writeHeader(hid_t fileId) const;
    void writeBox(hid_t fileId);
    void writeId(hid_t fileId);
    void writeType(hid_t fileId);
    void writeMass(hid_t fileId);
    void writeQ(hid_t fileId);
    void writeGhost(hid_t fileId);
    void writePosition(hid_t fileId);
    void writeVelocity(hid_t fileId);
    void writeForce(hid_t fileId);

    template <typename T>
    void writeParallel(hid_t fileId,
                       const std::string& name,
                       const std::vector<hsize_t>& globalDims,
                       const std::vector<hsize_t>& localDims,
                       const std::vector<T>& data);

    shared_ptr<System> system_ = nullptr;
    std::string filename_ = "";  ///< output filename

    MPI_Comm comm = MPI_COMM_NULL;
    int rank = -1;
    int numProcesses = -1;

    int64_t numLocalParticles = -1;
    int64_t numTotalParticles = -1;
    /// Offset of the local particle chunk in the global particle array.
    int64_t particleOffset = -1;
};

}  // namespace io
}  // namespace espressopp