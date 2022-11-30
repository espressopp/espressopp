/*
  Copyright (C) 2021
      Max Planck Institute for Polymer Research & JGU Mainz
  Copyright (C) 2016
      Jakub Krajniak (jkrajniak at gmail.com)
  Copyright (C) 2012,2013,2016
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

#include "FixedPairList.hpp"
#include "vec/storage/StorageVec.hpp"

#include "python.hpp"
#include "storage/Storage.hpp"
#include "boost/bind.hpp"
#include "boost/serialization/vector.hpp"
#include "esutil/Error.hpp"
#include "Buffer.hpp"
#include "checks.hpp"

#include <sstream>

using namespace std;

namespace espressopp
{
namespace vec
{
LOG4ESPP_LOGGER(FixedPairList::theLogger, "FixedPairList");

FixedPairList::FixedPairList(std::shared_ptr<espressopp::storage::Storage> storage) : globalPairs()
{
    LOG4ESPP_INFO(theLogger, "construct FixedPairList");

    if (!storage->getSystem()->vectorization)
    {
        throw std::runtime_error("system has no vectorization");
    }
    vectorization = storage->getSystem()->vectorization;

    if (!(vectorization->storageVec))
        throw std::runtime_error("vectorization->storageVec cannot be null");
    auto& storageVec = vectorization->storageVec;
    storageVec->enableLocalParticles();

    sigBeforeSend = storage->beforeSendParticles.connect(std::bind(
        &FixedPairList::beforeSendParticles, this, std::placeholders::_1, std::placeholders::_2));
    sigAfterRecv = storage->afterRecvParticles.connect(std::bind(
        &FixedPairList::afterRecvParticles, this, std::placeholders::_1, std::placeholders::_2));
    sigOnParticlesChanged =
        storage->onParticlesChanged.connect(std::bind(&FixedPairList::onParticlesChanged, this));

    if (vectorization->getVecLevel() != 2)
    {
        throw std::runtime_error("espressopp::vec::FixedPairList can only be used for vecLevel=2");
    }
}

FixedPairList::~FixedPairList()
{
    LOG4ESPP_INFO(theLogger, "~FixedPairList");
    sigBeforeSend.disconnect();
    sigAfterRecv.disconnect();
    sigOnParticlesChanged.disconnect();
}

real FixedPairList::getLongtimeMaxBondSqr() { return longtimeMaxBondSqr; }

void FixedPairList::setLongtimeMaxBondSqr(real d) { longtimeMaxBondSqr = d; }

void FixedPairList::resetLongtimeMaxBondSqr() { longtimeMaxBondSqr = 0.0; }

bool FixedPairList::add(size_t pid1, size_t pid2)
{
    if (pid1 > pid2) std::swap(pid1, pid2);

    auto const& storageVec = vectorization->storageVec;
    size_t const p1 = storageVec->lookupRealParticleVec(pid1);
    size_t const p2 = storageVec->lookupLocalParticleVec(pid2);

    if (p1 == VEC_PARTICLE_NOT_FOUND)
    {
        // Particle does not exist here, return false
        LOG4ESPP_DEBUG(theLogger, "Leaving add with returnVal " << false);
        return false;
    }
    else
    {
        if (p2 == VEC_PARTICLE_NOT_FOUND)
        {
            LOG4ESPP_DEBUG(theLogger, "Particle p2 " << pid2 << " not found");
        }
    }

    {
        // add the pair locally
        this->push_back({p1, p2});

        // ADD THE GLOBAL PAIR
        // see whether the particle already has pairs
        const auto equalRange = globalPairs.equal_range(pid1);
        if (equalRange.first == globalPairs.end())
        {
            // if it hasn't, insert the new pair
            globalPairs.insert(std::make_pair(pid1, pid2));
        }
        else
        {
            // otherwise test whether the pair already exists
            for (auto it = equalRange.first; it != equalRange.second; ++it)
            {
                if (it->second == pid2)
                {
                    throw std::runtime_error("Duplicate pair inserted");
                }
            }
            // if not, insert the new pair
            globalPairs.insert(equalRange.first, std::make_pair(pid1, pid2));
        }
        LOG4ESPP_INFO(theLogger, "added fixed pair to global pair list");
    }
    LOG4ESPP_DEBUG(theLogger, "Leaving add with returnVal " << true);
    return true;
}

python::list FixedPairList::getBonds()
{
    python::tuple bond;
    python::list bonds;
    for (auto it = globalPairs.cbegin(); it != globalPairs.cend(); it++)
    {
        bond = python::make_tuple(it->first, it->second);
        bonds.append(bond);
    }

    return bonds;
}

std::vector<size_t> FixedPairList::getPairList()
{
    std::vector<size_t> ret;
    for (auto it = globalPairs.cbegin(); it != globalPairs.cend(); it++)
    {
        ret.push_back(it->first);
        ret.push_back(it->second);
    }
    return ret;
}

python::list FixedPairList::getAllBonds()
{
    std::vector<size_t> local_bonds = getPairList();
    std::vector<std::vector<size_t> > global_bonds;
    python::list bonds;

    for (auto it = globalPairs.cbegin(); it != globalPairs.cend(); it++)
    {
        local_bonds.push_back(it->first);
        local_bonds.push_back(it->second);
    }
    System& system = vectorization->getSystemRef();
    if (system.comm->rank() == 0)
    {
        mpi::gather(*system.comm, local_bonds, global_bonds, 0);
        for (auto it = global_bonds.begin(); it != global_bonds.end(); ++it)
        {
            CHECK_EQUAL(it->size() % 2, 0, "Size needs to be a multiple of two!");
            for (auto iit = it->begin(); iit != it->end(); iit += 2)
            {
                size_t pid1 = *(iit);
                size_t pid2 = *(iit + 1);
                bonds.append(python::make_tuple(pid1, pid2));
            }
        }
    }
    else
    {
        mpi::gather(*system.comm, local_bonds, global_bonds, 0);
    }
    return bonds;
}

void FixedPairList::beforeSendParticles(ParticleList& pl, OutBuffer& buf)
{
    std::vector<size_t> toSend;
    // loop over the particle list
    for (ParticleList::Iterator pit(pl); pit.isValid(); ++pit)
    {
        size_t pid = pit->id();

        LOG4ESPP_DEBUG(theLogger, "send particle with pid " << pid << ", find pairs");

        // find all pairs that involve this particle

        int n = globalPairs.count(pid);

        if (n > 0)
        {
            const auto equalRange = globalPairs.equal_range(pid);

            // first write the pid of the first particle
            // then the number of partners
            // and then the pids of the partners
            toSend.reserve(toSend.size() + n + 1);
            toSend.push_back(pid);
            toSend.push_back(n);
            for (auto it = equalRange.first; it != equalRange.second; ++it)
            {
                toSend.push_back(it->second);
                LOG4ESPP_DEBUG(theLogger,
                               "send global bond: pid " << pid << " and partner " << it->second);
            }

            // delete all of these pairs from the global list
            globalPairs.erase(equalRange.first, equalRange.second);
            // std::cout << "erasing particle " << pid << " from here" << std::endl;
        }
    }
    // send the list
    buf.write(toSend);
    LOG4ESPP_INFO(theLogger, "prepared fixed pair list before send particles");
}

void FixedPairList::afterRecvParticles(ParticleList& pl, InBuffer& buf)
{
    std::vector<size_t> received;
    auto it = globalPairs.begin();
    // receive the bond list
    buf.read(received);
    int size = received.size();
    int i = 0;
    while (i < size)
    {
        // unpack the list
        size_t pid1 = received[i++];
        int n = received[i++];
        LOG4ESPP_DEBUG(theLogger, "recv particle " << pid1 << ", has " << n << " global pairs");
        for (; n > 0; --n)
        {
            size_t pid2 = received[i++];
            // add the bond to the global list
            LOG4ESPP_DEBUG(theLogger, "received pair " << pid1 << " , " << pid2);
            it = globalPairs.insert(it, std::make_pair(pid1, pid2));
        }
    }
    if (i != size)
    {
        LOG4ESPP_ERROR(theLogger, "ATTENTION:  recv particles might have read garbage\n");
    }
    LOG4ESPP_INFO(theLogger, "received fixed pair list after receive particles");
}

void FixedPairList::onParticlesChanged()
{
    LOG4ESPP_INFO(theLogger, "rebuild local bond list from global\n");

    auto const& storageVec = vectorization->storageVec;

    System& system = vectorization->getSystemRef();
    esutil::Error err(system.comm);

    this->clear();
    size_t lastpid1 = VEC_PARTICLE_NOT_FOUND;
    size_t p1 = VEC_PARTICLE_NOT_FOUND;
    size_t p2 = VEC_PARTICLE_NOT_FOUND;
    for (auto it = globalPairs.cbegin(); it != globalPairs.cend(); ++it)
    {
        if (it->first != lastpid1)
        {
            p1 = storageVec->lookupRealParticleVec(it->first);
            if (p1 == VEC_PARTICLE_NOT_FOUND)
            {
                std::stringstream msg;
                msg << "onParticlesChanged error. Fixed Pair List particle p1 " << it->first
                    << " does not exist here. "
                    << "Checking pair (" << it->first << "," << it->second << ")";
                err.setException(msg.str());
            }
            lastpid1 = it->first;
        }
        p2 = storageVec->lookupLocalParticleVec(it->second);
        if (p2 == VEC_PARTICLE_NOT_FOUND)
        {
            std::stringstream msg;
            msg << "onParticlesChanged error. Fixed Pair List particle p2 " << it->second
                << " does not exist here. "
                << "Checking pair (" << it->first << "," << it->second << ")";
            err.setException(msg.str());
        }
        this->push_back({p1, p2});
    }
    err.checkException();

    LOG4ESPP_INFO(theLogger, "regenerated local fixed pair list from global list");
}

void FixedPairList::remove()
{
    this->clear();
    globalPairs.clear();
    sigBeforeSend.disconnect();
    sigAfterRecv.disconnect();
    sigOnParticlesChanged.disconnect();
}

int FixedPairList::totalSize()
{
    int local_size = globalPairs.size();
    int global_size;
    System& system = vectorization->getSystemRef();
    mpi::all_reduce(*system.comm, local_size, global_size, std::plus<int>());
    return global_size;
}

/****************************************************
** REGISTRATION WITH PYTHON
****************************************************/

void FixedPairList::registerPython()
{
    using namespace espressopp::python;

    bool (FixedPairList::*pyAdd)(size_t pid1, size_t pid2) = &FixedPairList::add;
    // bool (FixedPairList::*pyAdd)(pvec pids) = &FixedPairList::add;

    class_<FixedPairList, std::shared_ptr<FixedPairList> >(
        "vec_FixedPairList", init<std::shared_ptr<espressopp::storage::Storage> >())
        .def("add", pyAdd)
        .def("size", &FixedPairList::size)
        .def("totalSize", &FixedPairList::totalSize)
        .def("getBonds", &FixedPairList::getBonds)
        .def("remove", &FixedPairList::remove)
        .def("getAllBonds", &FixedPairList::getAllBonds)
        .def("resetLongtimeMaxBondSqr", &FixedPairList::resetLongtimeMaxBondSqr)
        .def("getLongtimeMaxBondSqr", &FixedPairList::getLongtimeMaxBondSqr);
}

}  // namespace vec
}  // namespace espressopp
