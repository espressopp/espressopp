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

#ifndef HPX4ESPP_STORAGE_STORAGEHPX_HPP
#define HPX4ESPP_STORAGE_STORAGEHPX_HPP

#include "log4espp.hpp"
#include "VirtualStorage.hpp"

#include <boost/signals2.hpp>

#include <iostream>

namespace espressopp
{
namespace hpx4espp
{
namespace storage
{
class StorageHPX
{
public:
    StorageHPX() {}

    /** Initialize connections between subdomains */
    virtual void initChannels() = 0;

    /** Copy particles to packed form. To be called at the start of integrator.run */
    virtual void loadCells() = 0;

    /** Copy particles back from packed form. To be called at the end of integrator.run */
    virtual void unloadCells() = 0;

    virtual void updateGhostsBlocking() = 0;

    virtual void collectGhostForcesBlocking() = 0;

    virtual void connect() = 0;

    virtual void disconnect() = 0;

    virtual void resetVirtualStorage() = 0;

    virtual void decomposeHPX() = 0;

    boost::signals2::signal<void()> onLoadCells;

    boost::signals2::signal<void()> onUnloadCells;

    static void registerPython();

public:
    std::vector<VirtualStorage> virtualStorage;

    boost::signals2::connection sigResetVirtualStorage;

private:
    static LOG4ESPP_DECL_LOGGER(logger);
};

}  // namespace storage
}  // namespace hpx4espp
}  // namespace espressopp

#endif  // HPX4ESPP_STORAGE_STORAGEHPX_HPP
