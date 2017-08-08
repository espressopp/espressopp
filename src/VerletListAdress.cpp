/*
  Copyright (C) 2012,2013,2014,2015,2016,2017(H)
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
#include "VerletListAdress.hpp"
#include "Real3D.hpp"
#include "Particle.hpp"
#include "Cell.hpp"
#include "System.hpp"
#include "storage/Storage.hpp"
#include "bc/BC.hpp"
#include "iterator/CellListAllPairsIterator.hpp"
#include "iterator/CellListIterator.hpp"

namespace espressopp {

  using namespace espressopp::iterator;

  LOG4ESPP_LOGGER(VerletListAdress::theLogger, "VerletList");

  /*-------------------------------------------------------------*/

    VerletListAdress::VerletListAdress(shared_ptr<System> system, real cut, real adrCut,
                                       bool rebuildVL, real _dEx, real _dHy)
    :SystemAccess(system) {
      LOG4ESPP_INFO(theLogger, "construct VerletList, cut = " << cut);

      if (!system->storage) {
         throw std::runtime_error("system has no storage");
      }
      skin = system->getSkin();
      cutverlet = cut + skin;
      cutsq = cutverlet * cutverlet;
      builds = 0;
      dEx = _dEx;
      dHy = _dHy;
      adrCenterSet = false;
      real adressSize = dEx + dHy + skin; // adress region size
      if (dEx + dHy == 0) adressSize = 0; // 0 should be 0
      adrsq = adressSize * adressSize;
      adrCutverlet = adrCut + skin;
      adrcutsq = adrCutverlet*adrCutverlet;

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
      const bc::BC& bc = *getSystemRef().bc;

      // get local cells
      CellList localcells = getSystem()->storage->getLocalCells();

      // if adrCenter is not set, the center of adress zone moves along with some particles
      if (!adrCenterSet) {
          // loop over all VP particles (reals and ghosts) on node
          for (CellListIterator it(localcells); it.isValid(); ++it) {

                // loop over positions
                for (std::vector<Real3D*>::iterator it2 = adrPositions.begin(); it2 != adrPositions.end(); ++it2){
                    Real3D dist;
                    real distsq;
                    bc.getMinimumImageVectorBox(dist, it->getPos(), **it2);

                    if (getAdrRegionType()){ // spherical adress region
                       distsq=dist.sqr();
                    }
                    else {  // slab-type adress region
                       distsq=dist[0]*dist[0];
                    }
                    if (distsq <= adrsq) {
                        adrZone.insert(&(*it));
                        break; // do not need to loop further
                    }
                }
                // if not near enough to any adrPositions, put in cgZone
                if (adrZone.count(&(*it)) == 0) {
                    cgZone.insert(&(*it));
                }
          }
      }
      // center of adress zone is fixed
      else {
          for (CellListIterator it(localcells); it.isValid(); ++it) {
              Real3D dist;
              real distsq;
              bc.getMinimumImageVectorBox(dist, it->getPos(), adrCenter);
              //Real3D dist = it->getPos() - adrCenter;                                   					        // CHANGE FOR X SPLIT VS SPHERE
              //real dist = it->getPos()[0] - adrCenter[0];                                  					        // CHANGE FOR X SPLIT VS SPHERE
              //real distsq = dist.sqr();                                                                                                 // CHANGE FOR X SPLIT VS SPHERE
              //real distsq = dist*dist;                                                                                                // CHANGE FOR X SPLIT VS SPHERE
              if (getAdrRegionType()){ // spherical adress region
                distsq=dist.sqr();
              }
              else {  // slab-type adress region
                distsq=dist[0]*dist[0];
              }
              if (distsq <= adrsq) {
                  adrZone.insert(&(*it));
              }
              else {
                  cgZone.insert(&(*it));
              }
          }

      }

      // add particles to adress pairs and VL
      CellList cl = getSystem()->storage->getRealCells();
      int count=0;

      for (CellListAllPairsIterator it(cl); it.isValid(); ++it) {
        checkPair(*it->first, *it->second);
        count+=1;
      }

      LOG4ESPP_INFO(theLogger, "rebuilt VerletList, cutsq = " << cutsq
                   << " local size = " << vlPairs.size());
      builds++;
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

      // see if it's in the exclusion list (both directions, CG particles only)
      if (exList.count(std::make_pair(pt1.id(), pt2.id())) == 1) return;
      if (exList.count(std::make_pair(pt2.id(), pt1.id())) == 1) return;
      // see if it's in the adress zone
      if (adrZone.count(&pt1) == 1 || adrZone.count(&pt2) == 1) {
          if (distsq > adrcutsq) return;
          adrPairs.add(pt1, pt2); // add to adress pairs
      }
      else {
          if (distsq > cutsq) return;
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

    void VerletListAdress::setAdrCenter(real x, real y, real z){
        adrCenter = Real3D(x, y, z);
        adrCenterSet = true;
        adrPositions.push_back(&adrCenter);
    }

    Real3D VerletListAdress::getAdrCenter(){return adrCenter;}

    void VerletListAdress::setAdrRegionType(bool _sphereAdr){
        sphereAdr = _sphereAdr;
    }

    bool VerletListAdress::getAdrRegionType(){
        return sphereAdr;
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
      using namespace espressopp::python;

      bool (VerletListAdress::*pyExclude)(longint pid1, longint pid2)
            = &VerletListAdress::exclude;

      void (VerletListAdress::*pyAddAdrParticle)(longint pid)
            = &VerletListAdress::addAdrParticle;

      void (VerletListAdress::*pySetAdrCenter)(real x, real y, real z)
                  = &VerletListAdress::setAdrCenter;

      void (VerletListAdress::*pySetAdrRegionType)(bool _sphereAdr)
            = &VerletListAdress::setAdrRegionType;

      class_<VerletListAdress, shared_ptr<VerletList> >
        ("VerletListAdress", init< shared_ptr<System>, real, real, bool, real, real>())
        .add_property("system", &SystemAccess::getSystem)
        .add_property("builds", &VerletListAdress::getBuilds, &VerletListAdress::setBuilds)
        .def("totalSize", &VerletListAdress::totalSize)
        .def("exclude", pyExclude)
        .def("addAdrParticle", pyAddAdrParticle)
        .def("setAdrCenter", pySetAdrCenter)
        .def("setAdrRegionType", pySetAdrRegionType)
        .def("rebuild", &VerletListAdress::rebuild)
        ;
    }

}
