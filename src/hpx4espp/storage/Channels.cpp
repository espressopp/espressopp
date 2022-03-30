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

#include <hpx/include/run_as.hpp>
#include <hpx/include/components.hpp>
#include "mpi.hpp"
#include "hpx4espp/utils/assert.hpp"

#include "Channels.hpp"

HPX_REGISTER_CHANNEL(hpx4espp_AlignedVectorChar, hpx4espp_AlignedVectorChar_comm);

#define CHIDX(COORD, LR, ISUB) ((ISUB) + numSubNodes[(COORD)] * (LR) + 2 * shiftSubNodes[(COORD)])

namespace espressopp
{
namespace hpx4espp
{
namespace storage
{
Channels::Channels(espressopp::storage::NodeGrid const& nodeGrid,
                   std::array<size_t, 3> const& subNodeGrid)
{
    const int rank = mpiWorld->rank();
    size_t totSubNodes = 0;
    for (int coord = 0; coord < 3; ++coord)
    {
        /// correctly skip self-interactions
        if (nodeGrid.getGridSize(coord) == 1)
        {
            numSubNodes[coord] = 0;
            shiftSubNodes[coord] = totSubNodes;
            totSubNodes += numSubNodes[coord];
        }
        else
        {
            auto subCommGrid = subNodeGrid;
            subCommGrid[coord] = 1;
            numSubNodes[coord] = subNodeGrid[0] * subNodeGrid[1] * subNodeGrid[2];
            shiftSubNodes[coord] = totSubNodes;
            totSubNodes += numSubNodes[coord];

            const int dir = 2 * coord + 0;
            const int oppDir = 2 * coord + 1;
            const int selfid = rank;
            int lr;

            /// senders
            {
                lr = 0;
                const int recver0 = nodeGrid.getNodeNeighborIndex(dir);
                for (int in = 0; in < numSubNodes[coord]; in++)
                    sendIdxs.push_back({lr, selfid, recver0, in});

                lr = 1;
                const int recver1 = nodeGrid.getNodeNeighborIndex(oppDir);
                for (int in = 0; in < numSubNodes[coord]; in++)
                    sendIdxs.push_back({lr, selfid, recver1, in});
            }

            /// recvers
            {
                lr = 0;
                const int sender0 = nodeGrid.getNodeNeighborIndex(oppDir);
                for (int in = 0; in < numSubNodes[coord]; in++)
                    recvIdxs.push_back({lr, sender0, selfid, in});

                lr = 1;
                const int sender1 = nodeGrid.getNodeNeighborIndex(dir);
                for (int in = 0; in < numSubNodes[coord]; in++)
                    recvIdxs.push_back({lr, sender1, selfid, in});
            }
        }
    }

    /// setup receiver channel
    {
        for (auto const& c : recvIdxs)
        {
            std::string name = getChannelName(c);
            recvChannels.push_back(ChannelType(hpx::find_here()));
            hpx::register_with_basename(name, recvChannels.back(), 0);
        }
    }

    /// setup sender channel
    {
        for (auto const& c : sendIdxs)
        {
            std::string name = getChannelName(c);
            sendChannels.push_back(hpx::find_from_basename<ChannelType>(name, 0));
        }
    }

    /// verify channels
    {
        verifyChannels();
    }
}

void Channels::verifyChannels()
{
    hpx::threads::run_as_hpx_thread([this] {
        const int rank = mpiWorld->rank();

        /// send the channel's data
        for (int coord = 0; coord < 3; ++coord)
        {
            const int nsub = numSubNodes[coord];
            for (int lr = 0; lr < 2; ++lr)
            {
                for (int isub = 0; isub < nsub; isub++)
                {
                    auto const& idx = sendIdxs[CHIDX(coord, lr, isub)];
                    AlignedVectorChar avc(sizeof(int) * 4);
                    int* buf = (int*)(avc.data());
                    buf[0] = std::get<0>(idx);
                    buf[1] = std::get<1>(idx);
                    buf[2] = std::get<2>(idx);
                    buf[3] = std::get<3>(idx);
                    // sendChannels[CHIDX(coord,lr,isub)].set(hpx::launch::sync, avc);
                    sendChannels[CHIDX(coord, lr, isub)].set(hpx::launch::apply, std::move(avc));
                    // sendChannels[CHIDX(coord,lr,isub)].set(hpx::launch::apply, avc);
                }
            }
        }

        std::vector<hpx::future<AlignedVectorChar>> futures;
        futures.reserve(2 * (numSubNodes[0] + numSubNodes[1] + numSubNodes[2]));
        for (int coord = 0; coord < 3; ++coord)
        {
            const int nsub = numSubNodes[coord];
            for (int lr = 0; lr < 2; ++lr)
            {
                for (int isub = 0; isub < nsub; isub++)
                {
                    // std::future<AlignedVectorChar> buf =
                    // recvChannels[CHIDX(coord,lr,isub)].get();
                    futures.push_back(recvChannels[CHIDX(coord, lr, isub)].get());
                }
            }
        }

        std::vector<AlignedVectorChar> bufs = hpx::util::unwrap(futures);
        for (int coord = 0; coord < 3; ++coord)
        {
            const int nsub = numSubNodes[coord];
            for (int lr = 0; lr < 2; ++lr)
            {
                for (int isub = 0; isub < nsub; isub++)
                {
                    const int chidx = CHIDX(coord, lr, isub);
                    auto const& avc = bufs[chidx];
                    HPX4ESPP_ASSERT_EQUAL(avc.size(), sizeof(int) * 4);
                    auto const& idx = recvIdxs[chidx];
                    const int* buf = (int*)(avc.data());
                    HPX4ESPP_ASSERT_EQUAL(buf[0], std::get<0>(idx));
                    HPX4ESPP_ASSERT_EQUAL(buf[1], std::get<1>(idx));
                    HPX4ESPP_ASSERT_EQUAL(buf[2], std::get<2>(idx));
                    HPX4ESPP_ASSERT_EQUAL(buf[3], std::get<3>(idx));
                }
            }
        }
    });
}

ChannelType& Channels::sendChannel(size_t coord, size_t lr, size_t isub)
{
    return sendChannels[CHIDX(coord, lr, isub)];
}

ChannelType& Channels::recvChannel(size_t coord, size_t lr, size_t isub)
{
    return recvChannels[CHIDX(coord, lr, isub)];
}

IndexType const& Channels::sendIdx(size_t coord, size_t lr, size_t isub) const
{
    return sendIdxs[CHIDX(coord, lr, isub)];
}

IndexType const& Channels::recvIdx(size_t coord, size_t lr, size_t isub) const
{
    return recvIdxs[CHIDX(coord, lr, isub)];
}

std::string Channels::getChannelName(int lr, int src, int dst, int isub) const
{
    std::stringstream ss;
    ss << "/hpx4espp_channel/" << lr << "/" << src << "/" << dst << "/" << isub << "/";
    return std::move(ss.str());
}

std::string Channels::getChannelName(IndexType const& p) const
{
    return std::move(
        getChannelName(std::get<0>(p), std::get<1>(p), std::get<2>(p), std::get<3>(p)));
}

/// store channel parameters into list of tuples
python::object Channels::getChannelIndices() const
{
    auto f_append_idx = [this](auto& l, auto const& idxList) {
        HPX4ESPP_ASSERT_EQUAL(idxList.size(),
                              2 * (numSubNodes[0] + numSubNodes[1] + numSubNodes[2]));

        auto s_coord = python::list();
        for (int coord = 0; coord < 3; ++coord)
        {
            const int nsub = numSubNodes[coord];
            auto s_lr = python::list();
            for (int lr = 0; lr < 2; ++lr)
            {
                auto s_isub = python::list();
                for (int isub = 0; isub < nsub; isub++)
                {
                    auto const& c = idxList.at(CHIDX(coord, lr, isub));
                    auto ll = python::list();
                    ll.append(std::get<0>(c));
                    ll.append(std::get<1>(c));
                    ll.append(std::get<2>(c));
                    ll.append(std::get<3>(c));
                    s_isub.append(ll);
                }
                s_lr.append(s_isub);
            }
            s_coord.append(s_lr);
        }
        l.append(s_coord);
    };

    auto l = python::list();
    f_append_idx(l, sendIdxs);
    f_append_idx(l, recvIdxs);
    return l;
}
}  // namespace storage
}  // namespace hpx4espp
}  // namespace espressopp
