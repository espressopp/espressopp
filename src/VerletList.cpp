/*
  Copyright (C) 2017
      Jakub Krajniak (jkrajniak at gmail.com)
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

LOG4ESPP_LOGGER(DynamicExcludeList::theLogger, "DynamicExcludeList");

DynamicExcludeList::DynamicExcludeList(shared_ptr<integrator::MDIntegrator> integrator):
    integrator_(integrator) {
  LOG4ESPP_INFO(theLogger, "construct of DynamicExcludeList");
  exList = boost::make_shared<ExcludeList>();
  connect();
  exList_remove.clear();
  exList_add.clear();
  is_dirty_ = false;

  system_ = integrator->getSystem();
}

DynamicExcludeList::~DynamicExcludeList() {
  disconnect();
}

void DynamicExcludeList::connect() {
  LOG4ESPP_INFO(theLogger, "Connected to integrator");
  befIntP = integrator_->befIntP.connect(boost::bind(&DynamicExcludeList::updateList, this));
  runInit = integrator_->runInit.connect(boost::bind(&DynamicExcludeList::updateList, this));
}

void DynamicExcludeList::disconnect() {
  LOG4ESPP_INFO(theLogger, "Disconnected from integrator");
  befIntP.disconnect();
  runInit.disconnect();
}

void DynamicExcludeList::observe_tuple(shared_ptr<FixedPairList> fpl) {
  fpl->onTupleAdded.connect(
      boost::bind(&DynamicExcludeList::exclude, this, _1, _2));

  fpl->onTupleRemoved.connect(
      boost::bind(&DynamicExcludeList::unexclude, this, _1, _2));
}

void DynamicExcludeList::observe_triple(shared_ptr<FixedTripleList> ftl) {
  ftl->onTupleAdded.connect(
      boost::bind(&DynamicExcludeList::exclude, this, _1, _3));
  ftl->onTupleRemoved.connect(
      boost::bind(&DynamicExcludeList::unexclude, this, _1, _3));
}

void DynamicExcludeList::observe_quadruple(shared_ptr<FixedQuadrupleList> fql) {
  fql->onTupleAdded.connect(
      boost::bind(&DynamicExcludeList::exclude, this, _1, _3));
  fql->onTupleAdded.connect(
      boost::bind(&DynamicExcludeList::exclude, this, _1, _4));
  fql->onTupleAdded.connect(
      boost::bind(&DynamicExcludeList::exclude, this, _2, _4));
  fql->onTupleRemoved.connect(
      boost::bind(&DynamicExcludeList::unexclude, this, _1, _3));
  fql->onTupleRemoved.connect(
      boost::bind(&DynamicExcludeList::unexclude, this, _1, _4));
  fql->onTupleRemoved.connect(
      boost::bind(&DynamicExcludeList::unexclude, this, _2, _4));
}

void DynamicExcludeList::updateList() {
  LOG4ESPP_INFO(theLogger, "Update dynamic list.");
  // check is dirty flag on all processes.
  bool global_is_dirty;
  mpi::all_reduce(*(system_->comm), is_dirty_, global_is_dirty, std::logical_or<bool>());

  if (!global_is_dirty)  // skip update
    return;

  // Collect state from all CPUs.
  std::vector<longint> out_buffer;

  // Prepare output.
  out_buffer.reserve(2 + exList_remove.size() + exList_add.size());
  out_buffer.push_back(exList_remove.size() / 2);
  out_buffer.push_back(exList_add.size() / 2);
  out_buffer.insert(out_buffer.end(), exList_remove.begin(), exList_remove.end());
  out_buffer.insert(out_buffer.end(), exList_add.begin(), exList_add.end());

  // Gather everywhere updates and apply.
  std::vector<std::vector<longint> > in_buffer;
  mpi::all_gather(*(integrator_->getSystem()->comm), out_buffer, in_buffer);
  //Update list.
  LOG4ESPP_DEBUG(theLogger, "update data from " << in_buffer.size());

  for (std::vector<std::vector<longint> >::iterator it = in_buffer.begin(); it != in_buffer.end(); it++) {
    for (std::vector<longint>::iterator itm = it->begin(); itm != it->end();) {
      longint remove_size = *(itm++);
      longint add_size = *(itm++);
      LOG4ESPP_DEBUG(theLogger, "remove_size=" << remove_size << " add_size=" << add_size);
      for (int i = 0; i < remove_size; i++) {
        longint f1 = *(itm++);
        longint f2 = *(itm++);
        exList->erase(std::make_pair(f1, f2));
        exList->erase(std::make_pair(f2, f1));
        onPairUnexclude(f1, f2);
      }
      for (int i = 0; i < add_size; i++) {
        longint f1 = *(itm++);
        longint f2 = *(itm++);
        exList->insert(std::make_pair(f1, f2));
        exList->insert(std::make_pair(f2, f1));
        onPairExclude(f1, f2);
      }
    }
  }
  exList_remove.clear();
  exList_add.clear();
  is_dirty_ = false;

  // Rebuild list.
  onListUpdated();

  LOG4ESPP_DEBUG(theLogger, "leave DynamicExcludeList::updateList");
}

python::list DynamicExcludeList::getList() {
  python::list return_list;
  for (ExcludeList::iterator it = exList->begin(); it != exList->end(); ++it) {
    return_list.append(python::make_tuple(it->first, it->second));
  }
  return return_list;
}

void DynamicExcludeList::exclude(longint pid1, longint pid2) {
  LOG4ESPP_INFO(theLogger, "new exclude pair " << pid1 << "-" << pid2);
  exList_add.push_back(pid1);
  exList_add.push_back(pid2);
  is_dirty_ = true;
  onPairExclude(pid1, pid2);
}

void DynamicExcludeList::unexclude(longint pid1, longint pid2) {
  LOG4ESPP_INFO(theLogger, "removed exclude pair " << pid1 << "-" << pid2);
  exList_remove.push_back(pid1);
  exList_remove.push_back(pid2);
  is_dirty_ = true;
  onPairUnexclude(pid1, pid2);
}

void DynamicExcludeList::registerPython() {
  using namespace espressopp::python;

  class_<DynamicExcludeList, shared_ptr<DynamicExcludeList>, boost::noncopyable >
      ("DynamicExcludeList", init< shared_ptr<integrator::MDIntegrator> >())
       .add_property("size", &DynamicExcludeList::getSize)
       .def("exclude", &DynamicExcludeList::exclude)
       .def("unexclude", &DynamicExcludeList::unexclude)
       .def("observe_tuple", &DynamicExcludeList::observe_tuple)
       .def("observe_triple", &DynamicExcludeList::observe_triple)
       .def("observe_quadruple", &DynamicExcludeList::observe_quadruple)
       .def("get_list", &DynamicExcludeList::getList)
       .def("update", &DynamicExcludeList::updateList)
       .def("connect", &DynamicExcludeList::connect)
       .def("disconnect", &DynamicExcludeList::disconnect);
}

/** Implementation of VerletList **/

  LOG4ESPP_LOGGER(VerletList::theLogger, "VerletList");
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

    exList = boost::make_shared<ExcludeList>();
    isDynamicExList = false;

    if (rebuildVL) rebuild(); // not called if exclutions are provided

    // make a connection to System to invoke rebuild on resort
    connect();
  }
  
  VerletList::VerletList(shared_ptr<System> system, real _cut,
                         shared_ptr<DynamicExcludeList> dynamicExList_, bool rebuildVL):
      SystemAccess(system), dynamicExcludeList(dynamicExList_) {
    LOG4ESPP_INFO(theLogger, "construct VerletList with dynamic exclusion list, cut = " << _cut);

    if (!system->storage) {
      throw std::runtime_error("system has no storage");
    }

    cut = _cut;
    cutVerlet = cut + system -> getSkin();
    cutsq = cutVerlet * cutVerlet;
    builds = 0;

    exList = dynamicExList_->getExList();

    isDynamicExList = true;

    dynamicExList_->onListUpdated.connect(boost::bind(&VerletList::rebuild, this));

    // proxy signals from DynamicExclude list to signals of VerletList.
    dynamicExcludeList->onPairExclude.connect(onPairExclude);
    dynamicExcludeList->onPairUnexclude.connect(onPairUnexclude);

    if (rebuildVL) rebuild(); // not called if exclutions are provided

    // make a connection to System to invoke rebuild on resort
    connect();

    resetTimers();
  }
  real VerletList::getVerletCutoff(){
    return cutVerlet;
  }

  void VerletList::setVerletCutoff(real _cut) {
    cut = _cut;
    cutVerlet = cut + getSystem()->getSkin();
    cutsq = cutVerlet * cutVerlet;
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
    real time0 = wallTimer.getElapsedTime();
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
    timeRebuild_ += wallTimer.getElapsedTime() - time0;
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
    if (exList->count(std::make_pair(pt1.id(), pt2.id())) == 1) return;
    //if (exList->count(std::make_pair(pt2.id(), pt1.id())) == 1) return;

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
      if (isDynamicExList) {
        dynamicExcludeList->exclude(pid1, pid2);
      } else {
        exList->insert(std::make_pair(pid1, pid2));
        exList->insert(std::make_pair(pid2, pid1));
        onPairExclude(pid1, pid2);
      }
      return true;
  }

  bool VerletList::unexclude(longint pid1, longint pid2) {
    if (isDynamicExList) {
      dynamicExcludeList->unexclude(pid1, pid2);
    } else {
      exList->erase(std::make_pair(pid1, pid2));
      exList->erase(std::make_pair(pid2, pid1));
      onPairUnexclude(pid1, pid2);
    }
  }

  longint VerletList::excludeListSize() const {
    return exList->size();
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


    class_<VerletList, shared_ptr<VerletList>, boost::noncopyable >
      ("VerletList", init< shared_ptr<System>, real, bool >())
      .def(init<shared_ptr<System>, real, shared_ptr<DynamicExcludeList>, bool>())
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
      .def("setVerletCutoff", &VerletList::setVerletCutoff)
      .def("get_timers", &VerletList::getTimers)
      .def("excludeListSize", &VerletList::excludeListSize)
      ;
  }

}
