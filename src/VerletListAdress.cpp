#include "python.hpp"
#include "VerletListAdress.hpp"
#include "Real3D.hpp"
#include "Particle.hpp"
#include "Cell.hpp"
#include "System.hpp"
#include "storage/Storage.hpp"
#include "bc/BC.hpp"
#include "iterator/CellListAllPairsIterator.hpp"

namespace espresso {

  using namespace espresso::iterator;

  LOG4ESPP_LOGGER(VerletListAdress::theLogger, "VerletList");

/*-------------------------------------------------------------*/

  VerletListAdress::VerletListAdress(shared_ptr<System> system, real cut, bool rebuildVL) : SystemAccess(system)
  {
    LOG4ESPP_INFO(theLogger, "construct VerletList, cut = " << cut);
  
    if (!system->storage) {
       throw std::runtime_error("system has no storage");
    }

    skin = system->skin;
    real cutVerlet = cut + skin;
    cutsq = cutVerlet * cutVerlet;
    builds = 0;

    adresscut = cutVerlet; //TODO now it's fixed
    adrsq = adresscut * adresscut;

    std::cout << "\n------constructor----- ";
    if (rebuildVL) rebuild(); // not called if exclutions are provided

  
    // make a connection to System to invoke rebuild on resort
    connectionResort = system->storage->onParticlesChanged.connect(
        boost::bind(&VerletListAdress::rebuild, this));
  }
  
  /*-------------------------------------------------------------*/
  
  void VerletListAdress::rebuild()
  {
    vlPairs.clear();
    adrZone.clear(); // particles in adress zone
    adrPairs.clear(); // pairs in adress zone
    atmPositions.clear(); // clear position pointers


    std::cout << "\n-- VL Rebuild --\n";

    // add particles to adress zone
    CellList cl = getSystem()->storage->getRealCells();
    for (CellListAllPairsIterator it(cl); it.isValid(); ++it) {
      if ((*it->first).type() > 100) continue;
      if ((*it->second).type() > 100) continue;
      isPairInAdrZone(*it->first, *it->second);
    }

    // add particles to VL
    for (CellListAllPairsIterator it(cl); it.isValid(); ++it) {
      if ((*it->first).type() > 100) continue;
      if ((*it->second).type() > 100) continue;
      checkPair(*it->first, *it->second);
    }

    std::cout << "verlet list pairs (vlPairs size " << vlPairs.size() << "):\n";
    //std::cout << "\n\n";



    // AdResS testing
    /*
    // print VL pairs of atomistic particles
    for (std::set<longint>::iterator it = atmList.begin(); it != atmList.end(); ++it) {
       std::cout << *it << " interacts with:\n";
       for (PairList::Iterator it2(vlPairs); it2.isValid(); ++it2) {
             if ((*it2->first).id() == *it) {
                 Real3D d = (*it2->first).position() - (*it2->second).position();
                 real distsq = d.sqr();
                 std::cout << " " << (*it2->second).id() << " (d = " << sqrt(distsq) << ")\n";
             }

             else if ((*it2->second).id() == *it) {
               Real3D d = (*it2->first).position() - (*it2->second).position();
               real distsq = d.sqr();
               std::cout << " " << (*it2->first).id() << " (d = " << sqrt(distsq) << ")*\n";
             }
       }
       std::cout << "\n";
    }*/


    // print particles in adress zone
    std::cout << "in adress zone (adrZone size " << adrZone.size() <<  "):\n";
    for (std::set<Particle*>::iterator it = adrZone.begin(); it != adrZone.end(); ++it) {
        std::cout << (*it)->id() << "-";
        std::cout << (*it)->ghost() << " ";
    }
    std::cout << "\n\n";


    // print adrPairs
    std::cout << "adress pairs (adrPairs size " << adrPairs.size() << "):\n";
    for (PairList::Iterator it(adrPairs); it.isValid(); ++it) {
        std::cout << "(" << (*it->first).id() << "-" << (*it->first).ghost() <<
                ", " << (*it->second).id() << "-" << (*it->second).ghost() << ") ";
    }
    std::cout << "\n\n";



    LOG4ESPP_INFO(theLogger, "rebuilt VerletList, cutsq = " << cutsq 
                 << " local size = " << vlPairs.size());
    builds++;
  }
  

  /*-------------------------------------------------------------*/

  // add particles to adress zone
  void VerletListAdress::isPairInAdrZone(Particle& pt1, Particle& pt2) {

      Real3D d = pt1.position() - pt2.position();
      real distsq = d.sqr();

      // AdResS stuff
      //if (distsq > adrsq) return;
      // check if one of the particles is an atomistic particle
      if (atmList.count(pt1.id()) == 1) {
          adrZone.insert(&pt1);
          atmPositions.push_back(&pt1.position());
          if (distsq > adrsq) return;
          else adrZone.insert(&pt2);
      }
      if (atmList.count(pt2.id()) == 1) {
          adrZone.insert(&pt2);
          atmPositions.push_back(&pt2.position());
          if (distsq > adrsq) return;
          else adrZone.insert(&pt1);
      }

          /*
          // insert and overwrite position
          std::pair<std::map<longint, Real3D>::iterator,bool> res;
          res = adrZone.insert(std::make_pair(pt1.id(), pt1.position()));
          if(!res.second) // key was already in map change it.
              res.first->second = pt1.position();
          res = adrZone.insert(std::make_pair(pt2.id(), pt2.position()));
          if(!res.second) // key was already in map change it.
              res.first->second = pt2.position();
          */
  }

  /*-------------------------------------------------------------*/
  
  void VerletListAdress::checkPair(Particle& pt1, Particle& pt2)
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
    // see if it's in the adress zone
    if (adrZone.count(&pt1) == 1 || adrZone.count(&pt2) == 1) {
        adrPairs.add(pt1, pt2); // add to adress pairs
    }
    else {
        //std::cout << "(" << pt1.id() << "-" << pt1.ghost() << ", " << pt2.id() << "-" << pt2.ghost() << ") ";
        vlPairs.add(pt1, pt2); // add pair to Verlet List
    }
  }
  
  /*-------------------------------------------------------------*/
  
  int VerletListAdress::totalSize() const
  {
    System& system = getSystemRef();
    int size = vlPairs.size();
    int allsize;
  
    mpi::all_reduce(*system.comm, size, allsize, std::plus<int>());
    return allsize;
  }


  bool VerletListAdress::exclude(longint pid1, longint pid2) {

      exList.insert(std::make_pair(pid1, pid2));

      return true;
  }
  
  void VerletListAdress::addAtParticle(longint pid) {
        atmList.insert(pid);
  }

  /*-------------------------------------------------------------*/
  
  VerletListAdress::~VerletListAdress()
  {
    LOG4ESPP_INFO(theLogger, "~VerletList");
  
    if (!connectionResort.connected()) {
      connectionResort.disconnect();
    }
  }
  
  /****************************************************
  ** REGISTRATION WITH PYTHON
  ****************************************************/
  
  void VerletListAdress::registerPython() {
    using namespace espresso::python;

    bool (VerletListAdress::*pyExclude)(longint pid1, longint pid2)
          = &VerletListAdress::exclude;

    void (VerletListAdress::*pyAddAtParticle)(longint pid)
          = &VerletListAdress::addAtParticle;

    class_<VerletListAdress, shared_ptr<VerletList> >
      ("VerletListAdress", init< shared_ptr<System>, real, bool >())
      .add_property("system", &SystemAccess::getSystem)
      .add_property("builds", &VerletListAdress::getBuilds, &VerletListAdress::setBuilds)
      .def("totalSize", &VerletListAdress::totalSize)
      .def("exclude", pyExclude)
      .def("addAtParticle", pyAddAtParticle)
      .def("rebuild", &VerletListAdress::rebuild)
      ;
  }

}
