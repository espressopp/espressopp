/*
  Copyright (C) 2016
      Jakub Krajniak <jkrajniak at gmail.com>

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
#include "VerletListHybrid.hpp"
#include "Real3D.hpp"
#include "Particle.hpp"
#include "Cell.hpp"
#include "System.hpp"
#include "storage/Storage.hpp"
#include "bc/BC.hpp"
#include "iterator/CellListAllPairsIterator.hpp"

namespace espressopp
{
using namespace espressopp::iterator;

/** Implementation of VerletListHybrid **/
LOG4ESPP_LOGGER(VerletListHybridAT::theLogger, "VerletListHybridAT");
LOG4ESPP_LOGGER(VerletListHybridCG::theLogger, "VerletListHybridCG");

void VerletListHybridAT::checkPair(Particle &pt1, Particle &pt2)
{
    if (pt1.vp() || pt2.vp())
    {  // only AT pairs
        LOG4ESPP_DEBUG(theLogger, " at skip pair " << pt1.id() << "-" << pt2.id()
                                                   << " vp1=" << pt1.vp() << " vp2=" << pt2.vp());
        return;
    }

    Real3D d = pt1.position() - pt2.position();
    real distsq = d.sqr();

    LOG4ESPP_TRACE(theLogger, "p1: " << pt1.id() << " @ " << pt1.position() << " - p2: " << pt2.id()
                                     << " @ " << pt2.position() << " -> distsq = " << distsq
                                     << " cutsq=" << cutsq);

    if (distsq > cutsq) return;

    // see if it's in the exclusion list (both directions)
    if (exList.count(std::make_pair(pt1.id(), pt2.id())) == 1) return;
    if (exList.count(std::make_pair(pt2.id(), pt1.id())) == 1) return;

    vlPairs.add(pt1, pt2);  // add pair to Verlet List
}

void VerletListHybridCG::checkPair(Particle &pt1, Particle &pt2)
{
    if ((!pt1.vp()) || (!pt2.vp()))
    {  // only CG pairs
        LOG4ESPP_DEBUG(theLogger, "cg skip pair " << pt1.id() << "-" << pt2.id()
                                                  << " vp1=" << pt1.vp() << " vp2=" << pt2.vp());
        return;
    }

    Real3D d = pt1.position() - pt2.position();
    real distsq = d.sqr();

    LOG4ESPP_TRACE(theLogger, "p1: " << pt1.id() << " @ " << pt1.position() << " - p2: " << pt2.id()
                                     << " @ " << pt2.position() << " -> distsq = " << distsq);

    if (distsq > cutsq) return;

    // see if it's in the exclusion list (both directions)
    if (exList.count(std::make_pair(pt1.id(), pt2.id())) == 1) return;
    if (exList.count(std::make_pair(pt2.id(), pt1.id())) == 1) return;

    vlPairs.add(pt1, pt2);  // add pair to Verlet List
}

/****************************************************
** REGISTRATION WITH PYTHON
****************************************************/
void VerletListHybridAT::registerPython()
{
    using namespace espressopp::python;

    bool (VerletListHybridAT::*pyExclude)(longint pid1, longint pid2) =
        &VerletListHybridAT::exclude;

    class_<VerletListHybridAT, shared_ptr<VerletListHybridAT>, bases<VerletList> >(
        "VerletListHybridAT", init<shared_ptr<System>, real, bool>())
        .add_property("system", &SystemAccess::getSystem)
        .add_property("builds", &VerletListHybridAT::getBuilds, &VerletListHybridAT::setBuilds)
        .def("totalSize", &VerletListHybridAT::totalSize)
        .def("localSize", &VerletListHybridAT::localSize)
        .def("getPair", &VerletListHybridAT::getPair)
        .def("exclude", pyExclude)
        .def("rebuild", &VerletListHybridAT::rebuild)
        .def("connect", &VerletListHybridAT::connect)
        .def("disconnect", &VerletListHybridAT::disconnect)
        .def("getVerletCutoff", &VerletListHybridAT::getVerletCutoff);
}

void VerletListHybridCG::registerPython()
{
    using namespace espressopp::python;

    bool (VerletListHybridCG::*pyExclude)(longint pid1, longint pid2) =
        &VerletListHybridCG::exclude;

    class_<VerletListHybridCG, shared_ptr<VerletListHybridCG>, bases<VerletList> >(
        "VerletListHybridCG", init<shared_ptr<System>, real, bool>())
        .add_property("system", &SystemAccess::getSystem)
        .add_property("builds", &VerletListHybridCG::getBuilds, &VerletListHybridCG::setBuilds)
        .def("totalSize", &VerletListHybridCG::totalSize)
        .def("localSize", &VerletListHybridCG::localSize)
        .def("getPair", &VerletListHybridCG::getPair)
        .def("exclude", pyExclude)
        .def("rebuild", &VerletListHybridCG::rebuild)
        .def("connect", &VerletListHybridCG::connect)
        .def("disconnect", &VerletListHybridCG::disconnect)
        .def("getVerletCutoff", &VerletListHybridCG::getVerletCutoff);
}

}  // namespace espressopp
