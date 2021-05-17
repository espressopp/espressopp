/*
  Copyright (C) 2021
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

#include "LocalParticles.hpp"
#include "vec/ParticleArray.hpp"

namespace espressopp
{
namespace vec
{
namespace storage
{
void LocalParticles::rebuild(ParticleArray const& pa, std::vector<size_t> const& uniqueCells)
{
    clear();

    auto const numCells = pa.numCells();
    auto const& cellRange = pa.cellRange();
    auto const& cellSizes = pa.sizes();
    auto const& pids = pa.id;

    for (auto const& ic : uniqueCells)
    {
        auto const start = cellRange[ic];
        auto const size = cellSizes[ic];
        auto const end = start + size;
        for (auto ip = start; ip < end; ip++)
        {
            auto const pid = pids[ip];
            insert({pid, ip});
        }
    }
}

}  // namespace storage
}  // namespace vec
}  // namespace espressopp
