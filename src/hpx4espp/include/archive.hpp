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

#ifndef HPX4ESPP_INCLUDE_ARCHIVE_HPP
#define HPX4ESPP_INCLUDE_ARCHIVE_HPP

namespace espressopp
{
namespace hpx4espp
{
/// Store T as as char array
template <typename T>
struct Archive
{
    Archive() {}

    Archive(T const& p) { *(T*)(archive.data()) = p; }

    T const& get() { return *((T*)(archive.data())); }

protected:
    std::array<char, sizeof(T)> archive;
};

}  // namespace hpx4espp
}  // namespace espressopp

#endif  // HPX4ESPP_INCLUDE_ARCHIVE_HPP