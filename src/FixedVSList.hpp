/*
  Copyright (C) 2014-2017
      Jakub Krajniak (jkrajniak at gmail.com)

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

// ESPP_CLASS
#ifndef _FixedVSList_HPP
#define _FixedVSList_HPP

#include "log4espp.hpp"

#include "Particle.hpp"
#include "esutil/ESPPIterator.hpp"
#include <boost/unordered_map.hpp>
#include <boost/signals2.hpp>

namespace espressopp
{
class FixedVSList
{
protected:
    shared_ptr<storage::Storage> storage;

public:
    typedef std::vector<size_t> tuple;
    typedef boost::unordered_map<size_t, tuple> GlobalTuples;
    GlobalTuples globalTuples;

    FixedVSList(shared_ptr<storage::Storage> _storage);
    ~FixedVSList();

    void add(size_t pid) { tmppids.push_back(pid); }
    void addTs()
    {
        addT(tmppids);
        tmppids.clear();
    }

    tuple& getATParticleIds(size_t pidK) { return globalTuples[pidK]; }

    static void registerPython();

private:
    tuple tmppids;
    bool addT(tuple pids);  // add tuple
    static LOG4ESPP_DECL_LOGGER(theLogger);
};
}  // namespace espressopp

#endif
