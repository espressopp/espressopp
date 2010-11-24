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

  VerletList::VerletList(shared_ptr< System > system, real cut) : SystemAccess(system)
  {
    LOG4ESPP_INFO(theLogger, "construct VerletList, cut = " << cut);
  
    if (!system->storage) {
       throw std::runtime_error("system has no storage");
    }

    real cutVerlet = cut + system->skin;
    cutsq = cutVerlet * cutVerlet;
    builds = 0;
    rebuild();
  
    // make a connection to System to invoke rebuild on resort
    connectionResort = system->storage->onParticlesChanged.connect(
        boost::bind(&VerletList::rebuild, this));
  }
  
  /*-------------------------------------------------------------*/
  
  void VerletList::rebuild()
  {
    myList.clear();
  
    CellList cl = getSystem()->storage->getRealCells();
    for (CellListAllPairsIterator it(cl); it.isValid(); ++it) {
      checkPair(*it->first, *it->second);
    }
  
    LOG4ESPP_INFO(theLogger, "rebuilt VerletList, cutsq = " << cutsq 
                 << " local size = " << myList.size());
    builds++;
  }
  
  /*-------------------------------------------------------------*/
  
  void VerletList::checkPair(Particle& pt1, Particle& pt2)
    {
    //Real3DRef pos1 = pt1.r.p;
    //Real3DRef pos2 = pt2.r.p;

    real pos1[3] = {pt1.r.p[0], pt1.r.p[1], pt1.r.p[2]};
    real pos2[3] = {pt2.r.p[0], pt2.r.p[1], pt2.r.p[2]};
    
    //Real3D d;
    real d[3];

    real distsq;

    d[0] = pos1[0];  
    d[1] = pos1[1];  
    d[2] = pos1[2];  

    d[0] -= pos2[0];
    d[1] -= pos2[1];
    d[2] -= pos2[2];

    distsq = d[0]*d[0] + d[1]*d[1] + d[2]*d[2];
  
    LOG4ESPP_TRACE(theLogger, "p1: " << pt1.p.id << " @ " << pos1
		   << " - p2: " << pt2.p.id << " @ " << pos2
		   << " -> distsq = " << distsq);
  
    if (distsq > cutsq) return;
    
    myList.add(pt1, pt2);
  }
  
  /*-------------------------------------------------------------*/
  
  int VerletList::totalSize() const
  {
    System& system = getSystemRef();
    int size = myList.size();
    int allsize;
  
    mpi::all_reduce(*system.comm, size, allsize, std::plus<int>());
    return allsize;
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

    class_<VerletList, shared_ptr<VerletList> >
      ("VerletList", init< shared_ptr<System>, real >())
      .add_property("system", &SystemAccess::getSystem)
      .add_property("builds", &VerletList::getBuilds, &VerletList::setBuilds)
      .def("totalSize", &VerletList::totalSize)
      ;
  }

}
