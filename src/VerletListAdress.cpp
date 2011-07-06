#include "python.hpp"
#include "VerletListAdress.hpp"
#include "Real3D.hpp"
#include "Particle.hpp"
#include "Cell.hpp"
#include "System.hpp"
#include "storage/Storage.hpp"
#include "bc/BC.hpp"
#include "iterator/CellListAllPairsIterator.hpp"
#include "iterator/CellListIterator.hpp" // TODO remove later, used just for testing

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
      adrPositions.clear(); // clear position pointers


      std::cout << "\n-- VL Rebuild --\n";

      // add particles to adress zone
      /*
      CellList cl = getSystem()->storage->getRealCells();
      int count = 0;
      for (CellListAllPairsIterator it(cl); it.isValid(); ++it) {
        std::cout << "p1: " << (*it->first).id() << "-" << (*it->first).ghost()
                << " p2: " << (*it->second).id() << "-" << (*it->second).ghost() << "\n";
        ++count;
        if ((*it->first).type() >= atType) continue;  // only check if VP/CG particle!
        if ((*it->second).type() >= atType) continue; // only check if VP/CG particle!
        isPairInAdrZone(*it->first, *it->second);
      }
      std::cout << "pairs: " << count << std::endl;
      */

      /*
      std::cout << "particles of all local cells:\n";
      int count = 0;
      CellList localcells = getSystem()->storage->getLocalCells();
      Cell* cellp;
      for (CellListIterator it(localcells); it.isValid(); ++it) {
          cellp = getSystem()->storage->mapPositionToCell(it->position());
          std::cout << it->id() << "-" << it->ghost() << " " << it->position()
                  << " in cell " << cellp - (getSystem()->storage->getFirstCell()) << "\n";
          ++count;
      }
      std::cout << "(" << count <<" particles)\n";
      */


      // add particles to adress pairs and VL
      CellList cl = getSystem()->storage->getRealCells();
      for (CellListAllPairsIterator it(cl); it.isValid(); ++it) {
        if ((*it->first).type() >= atType) continue;  // only check if VP/CG particle!
        if ((*it->second).type() >= atType) continue; // only check if VP/CG particle!
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
      std::cout << "\nin adress zone (adrZone size " << adrZone.size() <<  "):\n";
      for (std::set<Particle*>::iterator it = adrZone.begin(); it != adrZone.end(); ++it) {
          std::cout << (*it)->id() << "-";
          std::cout << (*it)->ghost() << " (";
          std::cout << (*it)->position() << ")\n";

      }
      std::cout << "\n";


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

        //if (distsq > adrsq) return;
        // check if one of the particles is an atomistic particle
        //std::cout << "isPairInAdrZone: " << pt1.id() << ", " << pt2.id() << " " << sqrt(distsq) << "\n";
        if (adrList.count(pt1.id()) == 1) {
            adrZone.insert(&pt1);
            //std::cout << "Add " << pt1.id() << " to adr zone, pos: " << pt1.position() <<"\n";
            adrPositions.insert(&pt1.position());
            if (distsq < adrsq) adrZone.insert(&pt2);
            //std::cout << "Add " << pt2.id() << " to adr zone\n";
        }
        if (adrList.count(pt2.id()) == 1) {
            adrZone.insert(&pt2); // it's a set, so no duplicates possible
            //std::cout << "Add " << pt2.id() << " to adr zone, pos: " << pt2.position() <<"\n";
            adrPositions.insert(&pt2.position());
            if (distsq < adrsq) adrZone.insert(&pt1); // it's a set, so no duplicates possible
            //std::cout << "Add " << pt1.id() << " to adr zone\n";
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

    void VerletListAdress::addAdrParticle(longint pid) {
          adrList.insert(pid);
    }

    // types above this number are considered atomistic
    void VerletListAdress::setAtType(size_t type) {
        atType = type;
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

      void (VerletListAdress::*pyAddAdrParticle)(longint pid)
            = &VerletListAdress::addAdrParticle;

      void (VerletListAdress::*pySetAtType)(size_t type)
            = &VerletListAdress::setAtType;

      class_<VerletListAdress, shared_ptr<VerletList> >
        ("VerletListAdress", init< shared_ptr<System>, real, bool >())
        .add_property("system", &SystemAccess::getSystem)
        .add_property("builds", &VerletListAdress::getBuilds, &VerletListAdress::setBuilds)
        .def("totalSize", &VerletListAdress::totalSize)
        .def("exclude", pyExclude)
        .def("addAdrParticle", pyAddAdrParticle)
        .def("rebuild", &VerletListAdress::rebuild)
        .def("setAtType", pySetAtType)
        ;
    }

}
