#include "python.hpp"
#include "VerletListTriple.hpp"
#include "Real3D.hpp"
#include "Particle.hpp"
#include "Cell.hpp"
#include "System.hpp"
#include "storage/Storage.hpp"
#include "bc/BC.hpp"
#include "iterator/CellListAllPairsIterator.hpp"

namespace espresso {

  using namespace espresso::iterator;

  LOG4ESPP_LOGGER(VerletListTriple::theLogger, "VerletListTriple");

/*-------------------------------------------------------------*/

  VerletListTriple::VerletListTriple(shared_ptr<System> system, real cut, bool rebuildVL) : SystemAccess(system)
  {
    LOG4ESPP_INFO(theLogger, "construct VerletListTriple, cut = " << cut);
  
    if (!system->storage) {
       throw std::runtime_error("system has no storage");
    }

    real cutVerlet = cut;
    cutsq = cutVerlet * cutVerlet;
    builds = 0;

    if (rebuildVL) rebuild(); // not called if exclutions are provided

  
    // make a connection to System to invoke rebuild on resort
    connectionResort = system->storage->onParticlesChanged.connect(
        boost::bind(&VerletListTriple::rebuild, this));
  }
  
  void VerletListTriple::connect()
  {

  // make a connection to System to invoke rebuild on resort
  connectionResort = getSystem()->storage->onParticlesChanged.connect(
      boost::bind(&VerletListTriple::rebuild, this));
  }

  void VerletListTriple::disconnect()
  {

  // disconnect from System to avoid rebuild on resort
  connectionResort.disconnect();
  }

  /*-------------------------------------------------------------*/
  
  void VerletListTriple::rebuild()
  {
    vlTriples.clear();

    // add particles to adress zone
    CellList cl = getSystem()->storage->getRealCells();
    LOG4ESPP_DEBUG(theLogger, "local cell list size = " << cl.size());
    for (CellListAllTriplesIterator it(cl); it.isValid(); ++it) {
      checkTriple(*it->first, *it->second, *it->third);
      LOG4ESPP_DEBUG(theLogger, "checking particles " << it->first->id() << " , " << it->second->id() << " and " << it->third());
    }

    builds++;
    LOG4ESPP_DEBUG(theLogger, "rebuilt VerletListTriple (count=" << builds << "), cutsq = " << cutsq
                 << " local size = " << vlTriples.size());
  }
  

  /*-------------------------------------------------------------*/
  
  void VerletListTriple::checkTriple(Particle& pt1, Particle& pt2, Particle &pt3)
  {

	Real3D d12 = pt1.position() - pt2.position();
	Real3D d13 = pt1.position() - pt3.position();
	Real3D d23 = pt2.position() - pt3.position();
    real distsq = max(d12, d13, d23).sqr();

    LOG4ESPP_TRACE(theLogger, "max(d12, d13, d23).sqr()" << distsq);

    if (distsq > cutsq) return;

    // see if it's in the exclusion list (all three possibilities)
    if (exList.count(std::make_tuple(pt1.id(), pt2.id(), pt3.id())) == 1) return;
    if (exList.count(std::make_tuple(pt2.id(), pt3.id(), pt1.id())) == 1) return;
    if (exList.count(std::make_tuple(pt3.id(), pt1.id(), pt2.id())) == 1) return;

    vlTriples.add(pt1, pt2, pt3); // add pair to Verlet List Triple
  }
  
  /*-------------------------------------------------------------*/
  
  int VerletListTriple::totalSize() const
  {
    System& system = getSystemRef();
    int size = vlTriples.size();
    int allsize;
  
    mpi::all_reduce(*system.comm, size, allsize, std::plus<int>());
    return allsize;
  }

  python::tuple VerletListTriple::getTriple(int i)
  {
	if (i > 0 && i <= vlTriples.size()) {
	   return python::make_tuple(vlTriples[i-1].first->id(), vlTriples[i-1].second->id(), vlTriples[i-1].third->id());
	}
  }


  bool VerletListTriple::exclude(longint pid1, longint pid2, longint pid3) {

      exList.insert(std::make_tuple(pid1, pid2, pid3));

      return true;
  }
  

  /*-------------------------------------------------------------*/
  
  VerletListTriple::~VerletListTriple()
  {
    LOG4ESPP_INFO(theLogger, "~VerletListTriple");
  
    if (!connectionResort.connected()) {
      connectionResort.disconnect();
    }
  }
  
  /****************************************************
  ** REGISTRATION WITH PYTHON
  ****************************************************/
  
  void VerletListTriple::registerPython() {
    using namespace espresso::python;

    bool (VerletListTriple::*pyExclude)(longint pid1, longint pid2, longint pid3)
          = &VerletListTriple::exclude;


    class_<VerletListTriple, shared_ptr<VerletListTriple> >
      ("VerletListTriple", init< shared_ptr<System>, real, bool >())
      .add_property("system", &SystemAccess::getSystem)
      .add_property("builds", &VerletListTriple::getBuilds, &VerletListTriple::setBuilds)
      .def("totalSize", &VerletListTriple::totalSize)
      .def("getPair", &VerletListTriple::getTriple)
      .def("exclude", pyExclude)
      .def("rebuild", &VerletListTriple::rebuild)
      .def("connect", &VerletListTriple::connect)
      .def("disconnect", &VerletListTriple::disconnect)
      ;
  }

}
