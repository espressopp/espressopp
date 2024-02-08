/*
  Copyright (C) 2019-2022
      Max Planck Institute for Polymer Research & JGU Mainz
  Copyright (C) 2012,2013
      Max Planck Institute for Polymer Research
  Copyright (C) 2008,2009,2010,2011
      Max-Planck-Institute for Polymer Research & Fraunhofer SCAI

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
#include <hpx/include/parallel_for_loop.hpp>

#include "VerletList.hpp"
#include "hpx4espp/include/logging.hpp"
#include "hpx4espp/include/errors.hpp"
#include "hpx4espp/utils/algorithms/for_loop.hpp"
#include "hpx4espp/utils/multithreading.hpp"
#include "hpx4espp/utils/assert.hpp"
#include "hpx4espp/utils/assert_msg.hpp"

#include "python.hpp"
#include "Real3D.hpp"
#include "Particle.hpp"
#include "Cell.hpp"
#include "System.hpp"
#include "storage/Storage.hpp"
#include "bc/BC.hpp"
#include "iterator/CellListAllPairsIterator.hpp"

#include <atomic>

namespace espressopp
{
namespace hpx4espp
{
using namespace espressopp::iterator;
using storage::VirtualStorage;

LOG4ESPP_LOGGER(VerletList::theLogger, "VerletList");

/*-------------------------------------------------------------*/

// cut is a cutoff (without skin)
VerletList::VerletList(shared_ptr<System> system,
                       shared_ptr<StorageHPX> storageHPX,
                       real _cut,
                       bool rebuildVL,
                       bool resortOnLoad)
    : SystemAccess(system), storageHPX(storageHPX)
{
    LOG4ESPP_INFO(theLogger, "construct VerletList, cut = " << _cut);

    if (!system->storage)
    {
        throw std::runtime_error("system has no storage");
    }

    cut = _cut;
    cutVerlet = cut + system->getSkin();
    cutsq = cutVerlet * cutVerlet;
    builds = 0;
    max_type = 0;

    resetTimers();
    if (rebuildVL) rebuild();  // not called if exclutions are provided

    if (resortOnLoad)
    {
        // make a connection to StorageHPX to invoke rebuild on loadCells
        connectionResort = storageHPX->onLoadCells.connect(boost::bind(&VerletList::rebuild, this));
    }
    else
    {
        // make a connection to System to invoke rebuild on resort
        connectionResort =
            system->storage->onParticlesChanged.connect(boost::bind(&VerletList::rebuild, this));
    }
}

real VerletList::getVerletCutoff() { return cutVerlet; }

void VerletList::connect()
{
    // make a connection to System to invoke rebuild on resort
    connectionResort =
        getSystem()->storage->onParticlesChanged.connect(boost::bind(&VerletList::rebuild, this));
}

void VerletList::disconnect()
{
    // disconnect from System to avoid rebuild on resort
    connectionResort.disconnect();
}

/*-------------------------------------------------------------*/

void VerletList::rebuild()
{
    HPX4ESPP_DEBUG_MSG("VerletList::rebuild()");
    timer.reset();
    real currTime = timer.getElapsedTime();

    cutVerlet = cut + getSystem()->getSkin();
    cutsq = cutVerlet * cutVerlet;

    if (neighborLists.size() != getVirtualStorage().size())
    {
        neighborLists.clear();
        neighborLists.resize(getVirtualStorage().size());
    }

    num_pairs = 0;
    const size_t nvs = getVirtualStorage().size();

    auto f = [this](size_t i) {
        const auto& vs = getVirtualStorage().at(i);
        auto& nls = neighborLists.at(i);

        int num_pairs = 0;
        size_t max_type = 0;

        auto f_rebuild = [&, this](vec::VerletList::NeighborList& nl, const auto& cnl) {
            nl.reset();
            if (vs.particles.size()) nl.rebuild<1>(cutsq, cnl, vs.particles);
            num_pairs += 2 * nl.num_pairs;
            max_type = std::max(max_type, nl.max_type);
        };

        f_rebuild(nls.realNbrs, vs.realNbrs);
        f_rebuild(nls.externalNbrs, vs.externalNbrs);

        // f_rebuild(nls.realNbrs, vs.cellNeighborList);

        const auto nin = vs.internalNbrs.size();
        if (nls.internalNbrs.size() != nin)
        {
            nls.internalNbrs.clear();
            nls.internalNbrs.resize(nin);
        }

        for (size_t in = 0; in < nin; in++)
        {
            const auto& nnode_cnl = vs.internalNbrs.at(in);
            const auto& cnl = nnode_cnl.second;
            const auto nnode = nnode_cnl.first;
            auto& nnode_nl = nls.internalNbrs.at(in);
            nnode_nl.first = nnode;
            auto& nl = nnode_nl.second;
            const auto& vsNbr = getVirtualStorage().at(nnode);
            nl.reset();
            if (vs.particles.size())
            {
                nl.rebuildMulti<1, 0>(cutsq, cnl, vs.particles, vsNbr.particles);
                num_pairs += nl.num_pairs;
                max_type = std::max(max_type, nl.max_type);
            }
        }
        nls.num_pairs = num_pairs;
        nls.max_type = max_type;
    };

    utils::parallelForLoop(0, nvs, f);

    /// TODO: verify that this operation is relatively fast ~O(nvs)
    /// or perform reduction using atomics inside prev loop
    for (size_t i = 0; i < nvs; i++)
    {
        const auto& nl = neighborLists.at(i);
        num_pairs += nl.num_pairs;
        max_type = std::max(max_type, nl.max_type);
    }

    timeRebuild += timer.getElapsedTime() - currTime;
    builds++;

    HPX4ESPP_DEBUG_MSG("VerletList::rebuild()...DONE");
}

/*-------------------------------------------------------------*/

int VerletList::totalSize() const
{
    System& system = getSystemRef();
    int size = localSize();
    int allsize;

    mpi::all_reduce(*system.comm, size, allsize, std::plus<int>());
    return allsize;
}

int VerletList::localSize() const
{
    System& system = getSystemRef();
    // return vlPairs.size();
    return num_pairs;
}

python::list VerletList::getSubSizes() const
{
    python::list sizes;
    for (const auto& mnl : neighborLists)
    {
        sizes.append(mnl.realNbrs.num_pairs);
    }
    return sizes;
}

#if 0
  python::tuple VerletList::getPair(int i) {
    if (i <= 0 || i > vlPairs.size()) {
      std::cout << "ERROR VerletList pair " << i << " does not exists" << std::endl;
      return python::make_tuple();
    } else {
      return python::make_tuple(vlPairs[i-1].first->id(), vlPairs[i-1].second->id());
    }
  }


  bool VerletList::exclude(longint pid1, longint pid2) {

      throw std::runtime_error("Exclusions in VerletList not implemented.");

      exList.insert(std::make_pair(pid1, pid2));

      return true;
  }
#endif

/*-------------------------------------------------------------*/

VerletList::~VerletList()
{
    LOG4ESPP_INFO(theLogger, "~VerletList");

    if (!connectionResort.connected())
    {
        connectionResort.disconnect();
    }
}

/*-------------------------------------------------------------*/

void VerletList::resetTimers() { timeRebuild = 0.0; }

void VerletList::loadTimers(real* t) { t[0] = timeRebuild; }

static boost::python::object wrapGetTimers(class VerletList* obj)
{
    real tms[1];
    obj->loadTimers(tms);
    return boost::python::make_tuple(tms[0]);
}

/****************************************************
** REGISTRATION WITH PYTHON
****************************************************/

void VerletList::registerPython()
{
    using namespace espressopp::python;

    bool (VerletList::*pyExclude)(longint pid1, longint pid2) = &VerletList::exclude;

    class_<VerletList, shared_ptr<VerletList>>(
        "hpx4espp_VerletList", init<shared_ptr<System>, shared_ptr<StorageHPX>, real, bool, bool>())
        .add_property("system", &SystemAccess::getSystem)
        .add_property("builds", &VerletList::getBuilds, &VerletList::setBuilds)
        .def("totalSize", &VerletList::totalSize)
        .def("localSize", &VerletList::localSize)
        .def("getSubSizes", &VerletList::getSubSizes)
        // .def("getPair", &VerletList::getPair)
        // .def("exclude", pyExclude)
        .def("rebuild", &VerletList::rebuild)
        // .def("rebuildPairs", &VerletList::rebuildPairs)
        .def("connect", &VerletList::connect)
        .def("disconnect", &VerletList::disconnect)
        .def("getVerletCutoff", &VerletList::getVerletCutoff)
        .def("resetTimers", &VerletList::resetTimers)
        .def("getTimers", wrapGetTimers)
        .def("preallocFactor", &VerletList::preallocFactor);
}

}  // namespace hpx4espp
}  // namespace espressopp
