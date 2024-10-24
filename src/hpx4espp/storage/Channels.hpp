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

#ifndef HPX4ESPP_STORAGE_CHANNELS_HPP
#define HPX4ESPP_STORAGE_CHANNELS_HPP

#include "hpx4espp/include/hpx_version.hpp"
#if (HPX4ESPP_HPX_VERSION_FULL >= 10500)
#include <hpx/channel.hpp>
#include <hpx/lcos_local/channel.hpp>
#else
#include <hpx/lcos/channel.hpp>
#include <hpx/local_lcos/channel.hpp>
#endif

#include <hpx/include/lcos.hpp>

#include "vec/include/simdconfig.hpp"
#include "storage/NodeGrid.hpp"
#include "python.hpp"

namespace espressopp
{
namespace hpx4espp
{
namespace storage
{
typedef espressopp::vec::AlignedVector<char> AlignedVectorChar;
typedef hpx::lcos::channel<AlignedVectorChar> ChannelType;
typedef std::tuple<int, int, int, int> IndexType;

class Channels
{
public:
    Channels(){};
    Channels(espressopp::storage::NodeGrid const& nodeGrid,
             std::array<size_t, 3> const& subNodeGrid);
    python::object getChannelIndices() const;
    void verifyChannels();

    ChannelType& sendChannel(size_t coord, size_t lr, size_t isub);
    ChannelType& recvChannel(size_t coord, size_t lr, size_t isub);

    IndexType const& sendIdx(size_t coord, size_t lr, size_t isub) const;
    IndexType const& recvIdx(size_t coord, size_t lr, size_t isub) const;

protected:
    std::vector<ChannelType> sendChannels;
    std::vector<ChannelType> recvChannels;

    std::vector<AlignedVectorChar> sendBufs;

    std::vector<IndexType> sendIdxs;
    std::vector<IndexType> recvIdxs;

    std::string getChannelName(int lr, int src, int dst, int isub) const;
    std::string getChannelName(IndexType const& params) const;

    std::array<size_t, 3> numSubNodes;
    std::array<size_t, 3> shiftSubNodes;
};

}  // namespace storage
}  // namespace hpx4espp
}  // namespace espressopp

typedef espressopp::hpx4espp::storage::AlignedVectorChar hpx4espp_AlignedVectorChar;

HPX_REGISTER_CHANNEL_DECLARATION(hpx4espp_AlignedVectorChar);

#endif  // HPX4ESPP_STORAGE_CHANNELS_HPP
