/*
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

#include "python.hpp"
#include "VerletList.hpp"
#include "Real3D.hpp"
#include "Particle.hpp"
#include "Cell.hpp"
#include "System.hpp"
#include "storage/Storage.hpp"
#include "bc/BC.hpp"
#include "iterator/CellListAllPairsIterator.hpp"

namespace espressopp {

  using namespace espressopp::iterator;

  LOG4ESPP_LOGGER(VerletList::theLogger, "VerletList");

/*-------------------------------------------------------------*/

  // cut is a cutoff (without skin)
  VerletList::VerletList(shared_ptr<System> system, real _cut, bool rebuildVL) : SystemAccess(system)
  {
    LOG4ESPP_INFO(theLogger, "construct VerletList, cut = " << _cut);
  
    if (!system->storage) {
       throw std::runtime_error("system has no storage");
    }

    cut = _cut;
    cutVerlet = cut + system -> getSkin();
    cutsq = cutVerlet * cutVerlet;
    builds = 0;
    max_type = 0;

    if (rebuildVL) rebuild(); // not called if exclutions are provided

  
    // make a connection to System to invoke rebuild on resort
    connectionResort = system->storage->onParticlesChanged.connect(
        boost::bind(&VerletList::rebuild, this));
  }
  
  real VerletList::getVerletCutoff(){
    return cutVerlet;
  }
  
  void VerletList::connect()
  {

  // make a connection to System to invoke rebuild on resort
  connectionResort = getSystem()->storage->onParticlesChanged.connect(
      boost::bind(&VerletList::rebuild, this));
  }

  void VerletList::disconnect()
  {

  // disconnect from System to avoid rebuild on resort
  connectionResort.disconnect();
  }

  /*-------------------------------------------------------------*/
  
  void VerletList::rebuild()
  {
    //real cutVerlet = cut + getSystem() -> getSkin();
    cutVerlet = cut + getSystem() -> getSkin();
    cutsq = cutVerlet * cutVerlet;
    
    vlPairs.clear();

    // add particles to adress zone
    CellList cl = getSystem()->storage->getRealCells();
    LOG4ESPP_DEBUG(theLogger, "local cell list size = " << cl.size());
    for (CellListAllPairsIterator it(cl); it.isValid(); ++it) {
      checkPair(*it->first, *it->second);
      LOG4ESPP_DEBUG(theLogger, "checking particles " << it->first->id() << " and " << it->second->id());
    }
    
    builds++;
    LOG4ESPP_DEBUG(theLogger, "rebuilt VerletList (count=" << builds << "), cutsq = " << cutsq
                 << " local size = " << vlPairs.size());
  }
  

  /*-------------------------------------------------------------*/
  
  void VerletList::checkPair(Particle& pt1, Particle& pt2)
  {

    Real3D d = pt1.position() - pt2.position();
    real distsq = d.sqr();

    LOG4ESPP_TRACE(theLogger, "p1: " << pt1.id()
                   << " @ " << pt1.position() 
		   << " - p2: " << pt2.id() << " @ " << pt2.position()
		   << " -> distsq = " << distsq);

    if (distsq > cutsq) return;

    // see if it's in the exclusion list (both directions)
    if (exList.count(std::make_pair(pt1.id(), pt2.id())) == 1) return;
    if (exList.count(std::make_pair(pt2.id(), pt1.id())) == 1) return;

    max_type = std::max(max_type, std::max(pt1.type(), pt2.type()));
    vlPairs.add(pt1, pt2); // add pair to Verlet List
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
    return vlPairs.size();
  }

  python::tuple VerletList::getPair(int i) {
	  if (i <= 0 || i > vlPairs.size()) {
	    std::cout << "ERROR VerletList pair " << i << " does not exists" << std::endl;
	    return python::make_tuple();
	  } else {
	    return python::make_tuple(vlPairs[i-1].first->id(), vlPairs[i-1].second->id());
	  }
  }


  bool VerletList::exclude(longint pid1, longint pid2) {

      exList.insert(std::make_pair(pid1, pid2));

      return true;
  }
  

  /*-------------------------------------------------------------*/
  
  VerletList::~VerletList()
  {
    LOG4ESPP_INFO(theLogger, "~VerletList");
  
    if (!connectionResort.connected()) {
      connectionResort.disconnect();
    }
  }
  
  /****************************************************
  ** REGISTRATION WITH PYTHON
  ****************************************************/
  
  void VerletList::registerPython() {
    using namespace espressopp::python;

    bool (VerletList::*pyExclude)(longint pid1, longint pid2)
          = &VerletList::exclude;


    class_<VerletList, shared_ptr<VerletList> >
      ("VerletList", init< shared_ptr<System>, real, bool >())
      .add_property("system", &SystemAccess::getSystem)
      .add_property("builds", &VerletList::getBuilds, &VerletList::setBuilds)
      .def("totalSize", &VerletList::totalSize)
      .def("localSize", &VerletList::localSize)
      .def("getPair", &VerletList::getPair)
      .def("exclude", pyExclude)
      .def("rebuild", &VerletList::rebuild)
      .def("connect", &VerletList::connect)
      .def("disconnect", &VerletList::disconnect)
    
      .def("getVerletCutoff", &VerletList::getVerletCutoff)
      ;
  }

}
