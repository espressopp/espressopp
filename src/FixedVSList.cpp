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

#include "python.hpp"
#include "FixedVSList.hpp"

#include <boost/bind.hpp>
#include "storage/Storage.hpp"
#include "Buffer.hpp"

#include "System.hpp"
#include "bc/BC.hpp"

#include "esutil/Error.hpp"

namespace espressopp
{
LOG4ESPP_LOGGER(FixedVSList::theLogger, "FixedVSList");

FixedVSList::FixedVSList(shared_ptr<storage::Storage> _storage) : storage(_storage), globalTuples()
{
    LOG4ESPP_INFO(theLogger, "construct FixedVSList");
}

FixedVSList::~FixedVSList() { LOG4ESPP_INFO(theLogger, "~FixedVSList"); }

bool FixedVSList::addT(tuple pids)
{
    bool returnVal = true;

    tuple::iterator it = pids.begin();
    longint pid_cg = *it;
    std::vector<longint> pidstmp;

    for (++it; it != pids.end(); ++it)
    {
        pidstmp.push_back(*it);
    }
    globalTuples.insert(make_pair(pid_cg, pidstmp));

    return returnVal;
}

/****************************************************
** REGISTRATION WITH PYTHON
****************************************************/
void FixedVSList::registerPython()
{
    using namespace espressopp::python;

    void (FixedVSList::*pyAdd)(longint pid) = &FixedVSList::add;

    class_<FixedVSList, shared_ptr<FixedVSList>, boost::noncopyable>(
        "FixedVSList", init<shared_ptr<storage::Storage> >())
        .def("add", pyAdd)
        .def("addTs", &FixedVSList::addTs);
}

}  // namespace espressopp
