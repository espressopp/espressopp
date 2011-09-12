#include "python.hpp"
#include "VerletListAdress.hpp"
#include "Real3D.hpp"
#include "Particle.hpp"
#include "Cell.hpp"
#include "System.hpp"
#include "storage/Storage.hpp"
#include "bc/BC.hpp"
#include "iterator/CellListAllPairsIterator.hpp"
#include "iterator/CellListIterator.hpp"

namespace espresso {

  using namespace espresso::iterator;

  LOG4ESPP_LOGGER(VerletListAdress::theLogger, "VerletList");

  /*-------------------------------------------------------------*/

    VerletListAdress::VerletListAdress(shared_ptr<System> system, real cut, bool rebuildVL, real _dEx, real _dHy) : SystemAccess(system)
    {
      LOG4ESPP_INFO(theLogger, "construct VerletList, cut = " << cut);

      if (!system->storage) {
         throw std::runtime_error("system has no storage");
      }

      skin = system->skin;
      real cutVerlet = cut + skin;
      cutsq = cutVerlet * cutVerlet;
      builds = 0;

      // AdResS stuff
      dEx = _dEx;
      dHy = _dHy;
      adrCenterSet = false;
      adresscut = dEx + dHy;
      adrsq = adresscut * adresscut;

      //std::cout << getSystem()->comm->rank() << ": " << "------constructor----- \n";
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
      cgZone.clear(); // particles in CG zone
      adrPairs.clear(); // pairs in adress zone


      //std::cout << getSystem()->comm->rank() << ": " << "-- VL Rebuild --\n";

      /*
      // add particles to adress zone -- not used anymore
      CellList cl = getSystem()->storage->getRealCells();
      //int count = 0;
      for (CellListAllPairsIterator it(cl); it.isValid(); ++it) {
        std::cout << "p1: " << (*it->first).id() << "-" << (*it->first).ghost()
                << " p2: " << (*it->second).id() << "-" << (*it->second).ghost() << "\n";
        //++count;
        //if ((*it->first).type() >= atType) continue;  // only check if VP/CG particle!
        //if ((*it->second).type() >= atType) continue; // only check if VP/CG particle!
        //isPairInAdrZone(*it->first, *it->second);
      }*/


      // get local cells
      CellList localcells = getSystem()->storage->getLocalCells();

      // if adrCenter is not set, the center of adress zone moves along with some particles
      if (!adrCenterSet) {
          // loop over all VP particles (reals and ghosts) on node
          //std::cout << "particles of all local cells:\n";
          //int count = 0;
          //Cell* cellp;
          adrPositions.clear(); // clear position pointers
          for (CellListIterator it(localcells); it.isValid(); ++it) {

              /*cellp = getSystem()->storage->mapPositionToCell(it->position());
              ++count;
              std::cout << it->id() << "-" << it->ghost() << " " << it->position()
                      << " in cell " << cellp - (getSystem()->storage->getFirstCell()) << "\n";*/

              if (adrList.count(it->id()) == 1) {
                  //std::cout << getSystem()->comm->rank() << ": " << " adding particle position (" << it->position() << ") to adrPositions and adrZone\n";
                  adrPositions.push_back(&(it->position()));
                  adrZone.insert(&(*it));
              }
          }
          //std::cout << "(" << count <<" particles)\n";

          // again, loop over all VP particles and check if they are close enough to adrPositions and add to adrZone
          //std::cout << "\nAdding particles to adrZone ...\n";
          Real3D dist;
          real distsq;
          for (CellListIterator it(localcells); it.isValid(); ++it) {

                /*cellp = getSystem()->storage->mapPositionToCell(it->position());
                std::cout << it->id() << "-" << it->ghost() << " " << it->position()
                        << " in cell " << cellp - (getSystem()->storage->getFirstCell()) << "\n";*/

                // loop over positions
                for (std::vector<Real3D*>::iterator it2 = adrPositions.begin(); it2 != adrPositions.end(); ++it2){
                    dist = it->getPos() - **it2;
                    distsq = dist.sqr();

                    //std::cout << "distance " << sqrt(distsq) << "\n";
                    if (distsq < adrsq) {
                        adrZone.insert(&(*it));
                        //std::cout << " added " << it->getId() << "-" << it->ghost() <<  "\n";
                        //std::cout << " adding particle " << it->getId() << "-" << it->ghost() << " to adrZone\n";
                        break; // do not need to loop further
                    }
                }
          }
      }
      // center of adress zone is fixed
      else {
          Real3D dist;
          real distsq;
          for (CellListIterator it(localcells); it.isValid(); ++it) {
              dist = it->getPos() - adrCenter;
              distsq = dist.sqr();
              //std::cout << "distance " << sqrt(distsq) << "\n";
              if (distsq < adrsq) {
                  adrZone.insert(&(*it));
                  //std::cout << " added " << it->getId() << "-" << it->ghost() <<  "\n";
                  //std::cout << " adding particle " << it->getId() << "-" << it->ghost() << " to adrZone\n";
              }
          }
      }


      // add particles to adress pairs and VL
      CellList cl = getSystem()->storage->getRealCells();
      //std::cout << "local cell list size = " << cl.size() << "\n";
      for (CellListAllPairsIterator it(cl); it.isValid(); ++it) {
        //if ((*it->first).type() >= atType) continue;  // only check if VP/CG particle!
        //if ((*it->second).type() >= atType) continue; // only check if VP/CG particle!
        //std::cout << "iterating over " << (*it->first).id() << ", " << (*it->second).id() << "\n";
        checkPair(*it->first, *it->second);
      }

      //std::cout << getSystem()->comm->rank() << ": " << "verlet list pairs (vlPairs size " << vlPairs.size() << "):\n";
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
      //std::cout << getSystem()->comm->rank() << ": " << "in adress zone (adrZone size " << adrZone.size() <<  "):\n";
      /*for (std::set<Particle*>::iterator it = adrZone.begin(); it != adrZone.end(); ++it) {
          std::cout << (*it)->id() << "-";
          std::cout << (*it)->ghost() << " (";
          std::cout << (*it)->position() << ")\n";
      }
      std::cout << "\n";*/


      // print adrPairs
      //std::cout << getSystem()->comm->rank() << ": " << "adress pairs (adrPairs size " << adrPairs.size() << "):\n";
      /*for (PairList::Iterator it(adrPairs); it.isValid(); ++it) {
          std::cout << "(" << (*it->first).id() << "-" << (*it->first).ghost() <<
                  ", " << (*it->second).id() << "-" << (*it->second).ghost() << ") ";
      }
      std::cout << "\n\n";*/



      LOG4ESPP_INFO(theLogger, "rebuilt VerletList, cutsq = " << cutsq
                   << " local size = " << vlPairs.size());
      builds++;
    }


    /*-------------------------------------------------------------*/

    // add particles to adress zone -- not used anymore
    /*
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

        // insert and overwrite position
        //std::pair<std::map<longint, Real3D>::iterator,bool> res;
        //res = adrZone.insert(std::make_pair(pt1.id(), pt1.position()));
        //if(!res.second) // key was already in map change it.
        //   res.first->second = pt1.position();
        //res = adrZone.insert(std::make_pair(pt2.id(), pt2.position()));
        //if(!res.second) // key was already in map change it.
        //    res.first->second = pt2.position();

    }*/


    /*-------------------------------------------------------------*/

    void VerletListAdress::checkPair(Particle& pt1, Particle& pt2)
    {

      Real3D d = pt1.position() - pt2.position();
      real distsq = d.sqr();

      LOG4ESPP_TRACE(theLogger, "p1: " << pt1.id()
                     << " @ " << pt1.position()
             << " - p2: " << pt2.id() << " @ " << pt2.position()
             << " -> distsq = " << distsq);

      //std::cout << "checkPair(" << pt1.id() << ", " << pt2.id() << ")\n";

      // see if it's in the exclusion list (both directions)
      if (exList.count(std::make_pair(pt1.id(), pt2.id())) == 1) return;
      if (exList.count(std::make_pair(pt2.id(), pt1.id())) == 1) return;
      // see if it's in the adress zone
      if (adrZone.count(&pt1) == 1 || adrZone.count(&pt2) == 1) {
          adrPairs.add(pt1, pt2); // add to adress pairs
          //std::cout << "adding pair (" << pt1.id() << ", " << pt2.id() << ")\n";
      }
      else {
          if (distsq > cutsq) return;
          //std::cout << "not adding, adding to VL (" << pt1.id() << "-" << pt1.ghost() << ", " << pt2.id() << "-" << pt2.ghost() << ")\n";
          vlPairs.add(pt1, pt2); // add pair to Verlet List
          cgZone.insert(&pt1);
          cgZone.insert(&pt2);
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

    void VerletListAdress::setAdrCenter(real x, real y, real z){
        adrCenter = Real3D(x, y, z);
        adrCenterSet = true;
        adrPositions.push_back(&adrCenter);
    }

    /* not used anymore
    // types above this number are considered atomistic
    void VerletListAdress::setAtType(size_t type) {
        atType = type;
    }*/

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

      void (VerletListAdress::*pySetAdrCenter)(real x, real y, real z)
                  = &VerletListAdress::setAdrCenter;

      /*void (VerletListAdress::*pySetAtType)(size_t type)
            = &VerletListAdress::setAtType;*/

      class_<VerletListAdress, shared_ptr<VerletList> >
        ("VerletListAdress", init< shared_ptr<System>, real, bool, real, real>())
        .add_property("system", &SystemAccess::getSystem)
        .add_property("builds", &VerletListAdress::getBuilds, &VerletListAdress::setBuilds)
        .def("totalSize", &VerletListAdress::totalSize)
        .def("exclude", pyExclude)
        .def("addAdrParticle", pyAddAdrParticle)
        .def("setAdrCenter", pySetAdrCenter)
        .def("rebuild", &VerletListAdress::rebuild)
        //.def("setAtType", pySetAtType)
        ;
    }

}
