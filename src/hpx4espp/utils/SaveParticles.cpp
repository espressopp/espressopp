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

#include <hpx/config.hpp>

#define HPX4ESPP_DEBUG
#include "hpx4espp/include/logging.hpp"

#include "hpx4espp/utils/assert.hpp"
#include "hpx4espp/utils/assert_msg.hpp"
#include "SaveParticles.hpp"

#include "python.hpp"
#include <fstream>

namespace espressopp
{
namespace hpx4espp
{
namespace utils
{
namespace SaveParticles
{
/// output - subnode, icell, particle id, position, velocity, force
/// filterLevel: 0 - all, 1 - owned, 2 - real
void SAVE_PARTICLE_ARRAY(boost::mpi::communicator const& comm,
                         std::vector<storage::VirtualStorage> const& vss,
                         std::string label,
                         int filterLevel)
{
    auto getEnvStr = [](std::string const& key)
    {
        char* val = std::getenv(key.c_str());
        return val == NULL ? std::string("") : std::string(val);
    };

    /// parse output filename
    const auto OUT_BASE = getEnvStr("OUT_BASE");
    HPX4ESPP_ASSERT_GT(OUT_BASE.size(), 0);
    std::ostringstream outFileOss;
    outFileOss << OUT_BASE << "---SAVE-" << label << ".csv";
    const auto outFile = outFileOss.str();

    HPX4ESPP_ASSERT_EQUAL_MSG(comm.size(), 1, "implemented only for 1 MPI rank");
    const auto& rank = comm.rank();

    struct SaveParticle
    {
        size_t rank, inode, icell, id;
        real p_x, p_y, p_z, v_x, v_y, v_z, f_x, f_y, f_z;
    };
    std::vector<SaveParticle> saveParticles;

    auto appendCell = [&](const auto& inode, size_t ic)
    {
        const auto& vs = vss[inode];
        const auto& pa = vs.particles;

        const auto cellStart = pa.cellRange()[ic];
        const auto cellEnd = cellStart + pa.sizes()[ic];

        for (size_t ip = cellStart; ip < cellEnd; ip++)
        {
            saveParticles.push_back({static_cast<size_t>(rank), inode, ic, pa.id[ip], pa.p_x[ip],
                                     pa.p_y[ip], pa.p_z[ip], pa.v_x[ip], pa.v_y[ip], pa.v_z[ip],
                                     pa.f_x[ip], pa.f_y[ip], pa.f_z[ip]});
        }
    };

    for (size_t inode = 0; inode < vss.size(); inode++)
    {
        const auto& vs = vss[inode];
        if (filterLevel == 0)
        {
            /// loop over all cells
            for (size_t ic = 0; ic < vs.localCells.size(); ic++)
            {
                appendCell(inode, ic);
            }
        }
        else if (filterLevel == 1)
        {
            /// loop over owned cells
            for (const auto ic : vs.ownCellsIdx)
            {
                appendCell(inode, ic);
            }
        }
        else if (filterLevel == 2)
        {
            /// loop over real cells
            for (const auto ic : vs.particles.realCells())
            {
                appendCell(inode, ic);
            }
        }
    }

    /// TODO: Gather to rank 0

    auto saveToFile = [](const auto& saveParticles, const auto& outFile)
    {
        HPX4ESPP_DEBUG_MSG("Saving " << saveParticles.size() << " particles to " << outFile);

        // verify that file does not yet exist
        {
            auto fileExists = [](const std::string& filename) -> bool
            {
                struct stat buf;
                if (stat(filename.c_str(), &buf) != -1)
                {
                    return true;
                }
                return false;
            };

            HPX4ESPP_ASSERT_EQUAL_MSG(fileExists(outFile), false,
                                      "File " << outFile << " already exists");
        }

        std::ofstream out(outFile);

        // save the items line by line with space separation
        size_t iline = 0;
        out << "line"
               ",rank"
               ",inode"
               ",icell"
               ",id"
               ",p_x"
               ",p_y"
               ",p_z"
               ",v_x"
               ",v_y"
               ",v_z"
               ",f_x"
               ",f_y"
               ",f_z"
            << std::endl;
        for (const auto& s : saveParticles)
        {
            out << (iline++) << "," << s.rank << "," << s.inode << "," << s.icell << "," << s.id
                << "," << s.p_x << "," << s.p_y << "," << s.p_z << "," << s.v_x << "," << s.v_y
                << "," << s.v_z << "," << s.f_x << "," << s.f_y << "," << s.f_z << std::endl;
        }
    };

    if (!rank) saveToFile(saveParticles, outFile);
}

}  // namespace SaveParticles
}  // namespace utils
}  // namespace hpx4espp
}  // namespace espressopp
