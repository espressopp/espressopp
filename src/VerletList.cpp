/*
  Copyright (C) 2020
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
  VerletList::VerletList(std::shared_ptr<System> system, real _cut, bool rebuildVL, bool useBuffers, bool useSOA)
    : SystemAccess(system), useBuffers(useBuffers), useSOA(useSOA)
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

    resetTimers();
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
    timer.reset();
    real currTime = timer.getElapsedTime();

    //real cutVerlet = cut + getSystem() -> getSkin();
    cutVerlet = cut + getSystem() -> getSkin();
    cutsq = cutVerlet * cutVerlet;
    
    vlPairs.clear();

    if(useBuffers) {
      rebuildUsingBuffers(exList.size(), useSOA);
    } else {
      // add particles to adress zone
      CellList cl = getSystem()->storage->getRealCells();
      LOG4ESPP_DEBUG(theLogger, "local cell list size = " << cl.size());
      for (CellListAllPairsIterator it(cl); it.isValid(); ++it) {
        checkPair(*it->first, *it->second);
        LOG4ESPP_DEBUG(theLogger, "checking particles " << it->first->id() << " and " << it->second->id());
      }
    }

    builds++;
    timeRebuild += timer.getElapsedTime() - currTime;
    LOG4ESPP_DEBUG(theLogger, "rebuilt VerletList (count=" << builds << "), cutsq = " << cutsq
                 << " local size = " << vlPairs.size());
  }

  /*-------------------------------------------------------------*/

  template< bool USE_EXCLUSION_LIST, bool USE_SOA >
  void VerletList::_rebuildUsingBuffers()
  {
    const CellList& realCells = getSystem()->storage->getRealCells();
    const size_t numRealCells = realCells.size();

    // stores the range of neighbor particles belonging to cell i: with end=c_range[i]
    std::vector<int> c_range;
    c_range.reserve(numRealCells);

    // get the number of particles in all neighbor cells
    size_t c_reserve = 0;
    for(size_t icell=0; icell<numRealCells; icell++) {
      size_t row_reserve = 0;
      for(NeighborCellInfo& nc: realCells[icell]->neighborCells) {
        if(!nc.useForAllPairs) {
          row_reserve += nc.cell->particles.size();
        }
      }
      c_reserve += row_reserve;
      c_range.push_back(c_reserve);
    }

    // resize buffer
    if(c_reserve > c_p.size()) {
      size_t c_resize = 2*c_reserve;
      c_p.resize(c_resize);
      if(USE_SOA){
        c_x.resize(c_resize);
        c_y.resize(c_resize);
        c_z.resize(c_resize);
      } else {
        c_pos.resize(c_resize);
      }
      c_type.resize(c_resize);
    }

    // exclusion list may have been added later so check size of c_id separately
    if(USE_EXCLUSION_LIST){
      if(c_reserve > c_id.size()) {
        size_t c_resize = 2*c_reserve;
        c_id.resize(c_resize);
      }
    }

    // fill buffer
    size_t ip = 0;
    for(size_t icell=0; icell<numRealCells; icell++) {
      size_t end = c_range[icell];
      for(NeighborCellInfo& nc: realCells[icell]->neighborCells) {
        if(!nc.useForAllPairs) {
          for(Particle& p: nc.cell->particles) {
            c_p[ip] = &p;
            if(USE_EXCLUSION_LIST)
              c_id[ip] = p.id();
            c_type[ip] = p.type();
            if(USE_SOA){
              const Real3D& pos = p.position();
              c_x[ip] = pos[0];
              c_y[ip] = pos[1];
              c_z[ip] = pos[2];
            } else {
              c_pos[ip] = p.position();
            }
            ip++;
          }
        }
      }
      if(ip!=end) std::runtime_error("Range mismatch.");
    }

    // rebuild neighbor list
    size_t start=0;
    for(size_t icell=0; icell<numRealCells; icell++) {
      size_t end=c_range[icell];
      ParticleList& particles = realCells[icell]->particles;
      size_t numParticles = particles.size();
      for(size_t p1=0; p1<numParticles; p1++) {
        Particle& part1 = particles[p1];

        // self-loop
        for(size_t p2=p1+1; p2<numParticles; p2++) {
          Particle& part2 = particles[p2];
          checkPair(part1, part2);
        }

        Real3D p1_pos;
        real x1, y1, z1;
        if(USE_SOA) {
          const Real3D& pos = part1.position();
          x1 = pos[0];
          y1 = pos[1];
          z1 = pos[2];
        } else {
          p1_pos = part1.position();
        }
        size_t id1;
        if(USE_EXCLUSION_LIST) {
          id1 = part1.id();
        }
        const size_t type1 = part1.type();

        // neighbor-loop
        for(size_t p2=start; p2<end; p2++) {
          real distsq;
          if(USE_SOA) {
            real dist_x = x1 - c_x[p2];
            real dist_y = y1 - c_y[p2];
            real dist_z = z1 - c_z[p2];
            distsq = dist_x*dist_x + dist_y*dist_y + dist_z*dist_z;
          } else {
            Real3D d = p1_pos - c_pos[p2];
            distsq = d.sqr();
          }

          if (distsq > cutsq) continue;

          if(USE_EXCLUSION_LIST) {
            size_t const& id2 = c_id[p2];
            if (exList.count(std::make_pair(id1, id2)) == 1) continue;
            if (exList.count(std::make_pair(id2, id1)) == 1) continue;
          }

          max_type = std::max(max_type, std::max(type1, c_type[p2]));
          vlPairs.add(&part1,c_p[p2]);
        }
      }
      start = end;
    }
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
  
  /*-------------------------------------------------------------*/

  void VerletList::resetTimers()
  {
    timeRebuild = 0.0;
  }

  void VerletList::loadTimers(real* t)
  {
    t[0] = timeRebuild;
  }

  static boost::python::object wrapGetTimers(class VerletList* obj)
  {
    real tms[1];
    obj->loadTimers(tms);
    return boost::python::make_tuple(tms[0]);
  }

  /****************************************************
  ** REGISTRATION WITH PYTHON
  ****************************************************/
  
  void VerletList::registerPython() {
    using namespace espressopp::python;

    bool (VerletList::*pyExclude)(longint pid1, longint pid2)
          = &VerletList::exclude;


    class_<VerletList, std::shared_ptr<VerletList> >
      ("VerletList", init< std::shared_ptr<System>, real, bool, bool, bool >())
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
      .def("resetTimers", &VerletList::resetTimers)
      .def("getTimers", wrapGetTimers)
      ;
  }

}
