#include "python.hpp"
#include "VerletList.hpp"
#include "Real3D.hpp"
#include "Particle.hpp"
#include "Cell.hpp"
#include "System.hpp"
#include "storage/Storage.hpp"
#include "bc/BC.hpp"
#include "iterator/CellListAllPairsIterator.hpp"

namespace espresso {

  using namespace espresso::iterator;

  LOG4ESPP_LOGGER(VerletList::theLogger, "VerletList");

/*-------------------------------------------------------------*/

  VerletList::VerletList(shared_ptr<System> system, real cut, bool rebuildVL) : SystemAccess(system)
  {
    LOG4ESPP_INFO(theLogger, "construct VerletList, cut = " << cut);
  
    if (!system->storage) {
       throw std::runtime_error("system has no storage");
    }

    real skin = system->skin;
    real cutVerlet = cut + skin;
    cutsq = cutVerlet * cutVerlet;
    builds = 0;

    if (rebuildVL) rebuild(); // not called if exclutions are provided

  
    // make a connection to System to invoke rebuild on resort
    connectionResort = system->storage->onParticlesChanged.connect(
        boost::bind(&VerletList::rebuild, this));
  }
  
  /*-------------------------------------------------------------*/
  
  void VerletList::rebuild()
  {
    vlPairs.clear();

    // add particles to adress zone
    CellList cl = getSystem()->storage->getRealCells();
    for (CellListAllPairsIterator it(cl); it.isValid(); ++it) {
      checkPair(*it->first, *it->second);
    }

    LOG4ESPP_INFO(theLogger, "rebuilt VerletList, cutsq = " << cutsq 
                 << " local size = " << vlPairs.size());
    builds++;
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

    vlPairs.add(pt1, pt2); // add pair to Verlet List
  }
  
  /*-------------------------------------------------------------*/
  
  int VerletList::totalSize() const
  {
    System& system = getSystemRef();
    int size = vlPairs.size();
    int allsize;
  
    mpi::all_reduce(*system.comm, size, allsize, std::plus<int>());
    return allsize;
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
    using namespace espresso::python;

    bool (VerletList::*pyExclude)(longint pid1, longint pid2)
          = &VerletList::exclude;


    class_<VerletList, shared_ptr<VerletList> >
      ("VerletList", init< shared_ptr<System>, real, bool >())
      .add_property("system", &SystemAccess::getSystem)
      .add_property("builds", &VerletList::getBuilds, &VerletList::setBuilds)
      .def("totalSize", &VerletList::totalSize)
      .def("exclude", pyExclude)
      .def("rebuild", &VerletList::rebuild)
      ;
  }

}
