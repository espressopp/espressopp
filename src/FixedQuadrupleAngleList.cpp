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
#include <sstream>
#include "FixedQuadrupleAngleList.hpp"
#include <boost/bind.hpp>
#include "storage/Storage.hpp"
#include "bc/BC.hpp"
#include "Buffer.hpp"
#include "esutil/Error.hpp"

namespace espressopp {
  
  /* We assume that the particle with pid2 is a real one, others can be ghosts.
   * It is in order to minimize the cize of the cell, because otherwise user should 
   * be sure that
   * (r4-r1) < cell size
   * distance between first and the last particle should be smaller then the cell size.
   * Now it is (r3-r1)
   * 
   * That means that the structure of our quadruple in our local list will be the next
   * multimap< pid2, pair<Triple < pid1, pid3, pid4 >, angle> >
   */

  LOG4ESPP_LOGGER(FixedQuadrupleAngleList::theLogger, "FixedQuadrupleAngleList");

  FixedQuadrupleAngleList::FixedQuadrupleAngleList(shared_ptr< storage::Storage > _storage) 
    : storage(_storage), quadruplesAngles()
  {
    LOG4ESPP_INFO(theLogger, "construct FixedQuadrupleAngleList");

    con1 = storage->beforeSendParticles.connect
      (boost::bind(&FixedQuadrupleAngleList::beforeSendParticles, this, _1, _2));
    con2 = storage->afterRecvParticles.connect
      (boost::bind(&FixedQuadrupleAngleList::afterRecvParticles, this, _1, _2));
    con3 = storage->onParticlesChanged.connect
      (boost::bind(&FixedQuadrupleAngleList::onParticlesChanged, this));
  }

  FixedQuadrupleAngleList::~FixedQuadrupleAngleList() {
    LOG4ESPP_INFO(theLogger, "~FixedQuadrupleAngleList");

    con1.disconnect();
    con2.disconnect();
    con3.disconnect();
  }

  bool FixedQuadrupleAngleList::
  add(longint pid1, longint pid2, longint pid3, longint pid4) {
    bool returnVal = true;
    System& system = storage->getSystemRef();
    esutil::Error err(system.comm);
    
    // ADD THE LOCAL QUADRUPLET
    Particle *p1 = storage->lookupLocalParticle(pid1);
    Particle *p2 = storage->lookupRealParticle(pid2);
    Particle *p3 = storage->lookupLocalParticle(pid3);
    Particle *p4 = storage->lookupLocalParticle(pid4);
    
    // at first we check the real particle
    if (!p2){
      returnVal = false;
      // Particle does not exist here, return false
    }
    else{
      if (!p1) {
        std::stringstream msg;
        msg << "adding error: quadruple particle p1 " << pid1 << " does not exists here and cannot be added";
        err.setException( msg.str() );
      }
      if (!p3) {
        std::stringstream msg;
        msg << "adding error: quadruple particle p3 " << pid3 << " does not exists here and cannot be added";
        err.setException( msg.str() );
      }
      if (!p4) {
        std::stringstream msg;
        msg << "adding error: quadruple particle p4 " << pid4 << " does not exists here and cannot be added";
        err.setException( msg.str() );
      }
    }
    
    err.checkException();
    
    if(returnVal){
      // add the quadruple locally
      this->add(p1, p2, p3, p4);

      Real3D pos1 = p1->position();
      Real3D pos2 = p2->position();
      Real3D pos3 = p3->position();
      Real3D pos4 = p4->position();

      Real3D r21(0.0,0.0,0.0);
      system.bc->getMinimumImageVectorBox(r21, pos2, pos1);
      Real3D r32(0.0,0.0,0.0);
      system.bc->getMinimumImageVectorBox(r32, pos3, pos2);
      Real3D r43(0.0,0.0,0.0);
      system.bc->getMinimumImageVectorBox(r43, pos4, pos3);

      Real3D rijjk = r21.cross(r32); // [r21 x r32]
      Real3D rjkkn = r32.cross(r43); // [r32 x r43]

      real rijjk_sqr = rijjk.sqr();
      real rjkkn_sqr = rjkkn.sqr();

      real rijjk_abs = sqrt(rijjk_sqr);
      real rjkkn_abs = sqrt(rjkkn_sqr);

      real inv_rijjk = 1.0 / rijjk_abs;
      real inv_rjkkn = 1.0 / rjkkn_abs;

      // cosine between planes
      real cos_phi = (rijjk * rjkkn) * (inv_rijjk * inv_rjkkn);
      if (cos_phi > 1.0) cos_phi = 1.0;
      else if (cos_phi < -1.0) cos_phi = -1.0;

      quadruplesAngles.insert(std::make_pair(pid2, 
              std::make_pair(Triple<longint, longint, longint>(pid1, pid3, pid4),
                            acos(cos_phi))));
    }
    LOG4ESPP_INFO(theLogger, "added fixed quadruple to global quadruple list");
    
    return returnVal;
  }

  // anyway it returns (pid1, pid2, pid3, pid4)
  python::list FixedQuadrupleAngleList::getQuadruples(){
	python::tuple quadruple;
	python::list quadruples;
	for (QuadruplesAngles::const_iterator it=quadruplesAngles.begin(); it != quadruplesAngles.end(); it++) {
      quadruple = python::make_tuple(it->second.first.first,
              it->first,
              it->second.first.second,
              it->second.first.third);
      quadruples.append(quadruple);
    }

	return quadruples;
  }
  
  // anyway it returns (pid1, pid2, pid3, pid4, angle)
  python::list FixedQuadrupleAngleList::getQuadruplesAngles(){
	python::tuple quadruple;
	python::list quadruples;
	for (QuadruplesAngles::const_iterator it=quadruplesAngles.begin(); it != quadruplesAngles.end(); it++) {
      quadruple = python::make_tuple(it->second.first.first,
              it->first,
              it->second.first.second,
              it->second.first.third,
              it->second.second);
      quadruples.append(quadruple);
    }
	return quadruples;
  }
  
  real FixedQuadrupleAngleList::getAngle(int pid1, int pid2, int pid3, int pid4){
    real returnVal = -3;
    
    QuadruplesAngles::iterator itr;
	QuadruplesAngles::iterator lastElement;
	
	// locate an iterator to the first pair object associated with key
	itr = quadruplesAngles.find(pid2);
	if (itr == quadruplesAngles.end())
      return returnVal; // no elements associated with key, so return immediately

	// get an iterator to the element that is one past the last element associated with key
	lastElement = quadruplesAngles.upper_bound(pid2);
    
    Triple<longint, longint, longint> neededTriple = Triple<longint, longint, longint>(pid1, pid3, pid4);

	for ( ; itr != lastElement; ++itr){
      if(neededTriple==itr->second.first){
        returnVal = itr->second.second;
        break;
      }
    }
	return returnVal;
  }
  
  void FixedQuadrupleAngleList::
  beforeSendParticles(ParticleList& pl, OutBuffer& buf) {
    
    std::vector< longint > toSendInt;
    std::vector< real > toSendReal;
    // loop over the particle list
    for (ParticleList::Iterator pit(pl); pit.isValid(); ++pit) {
      longint pid = pit->id();
      
      // find all quadruples that involve this particle
      int n = quadruplesAngles.count(pid);

      if (n > 0) {
        std::pair<QuadruplesAngles::const_iterator, 
                QuadruplesAngles::const_iterator> equalRange 
                = quadruplesAngles.equal_range(pid);

        // first write the pid of this particle
        // then the number of partners (n)
        // and then the pids of the partners
        toSendInt.reserve(toSendInt.size()+3*n+1);
        toSendReal.reserve(toSendReal.size()+n);
        toSendInt.push_back(pid);
        toSendInt.push_back(n);
        for (QuadruplesAngles::const_iterator it = equalRange.first;
            it != equalRange.second; ++it) {
          toSendInt.push_back(it->second.first.first);
          toSendInt.push_back(it->second.first.second);
          toSendInt.push_back(it->second.first.third);
          toSendReal.push_back(it->second.second);
        }

        // delete all of these quadruples from the global list
        quadruplesAngles.erase(pid);
      }
    }
    // send the list
    buf.write(toSendInt);
    buf.write(toSendReal);
    LOG4ESPP_INFO(theLogger, "prepared fixed quadruple list before send particles");
  }

  void FixedQuadrupleAngleList::
  afterRecvParticles(ParticleList &pl, InBuffer& buf) {

    std::vector< longint > receivedInt;
    std::vector< real > receivedReal;
    int n;
    longint pid1, pid2, pid3, pid4;
    real angleVal;
    QuadruplesAngles::iterator it = quadruplesAngles.begin();
    // receive the quadruple list
    buf.read(receivedInt);
    buf.read(receivedReal);
    int size = receivedInt.size(); int i = 0; int j = 0;
    while (i < size) {
      // unpack the list
      pid2 = receivedInt[i++];
      n = receivedInt[i++];
      for (; n > 0; --n) {
        pid1 = receivedInt[i++];
        pid3 = receivedInt[i++];
        pid4 = receivedInt[i++];
        angleVal = receivedReal[j++];
        // add the quadruple to the global list
        it = quadruplesAngles.insert(it, std::make_pair(pid2,
              std::make_pair(Triple<longint, longint, longint>(pid1, pid3, pid4),
                             angleVal)));
      }
    }
    if (i != size) {
      printf("ATTETNTION:  recv particles might have read garbage\n");
    }
    LOG4ESPP_INFO(theLogger, "received fixed quadruple list after receive particles");
  }

  void FixedQuadrupleAngleList::onParticlesChanged() {
    System& system = storage->getSystemRef();
    esutil::Error err(system.comm);
    
    // (re-)generate the local quadruple list from the global list
    this->clear();
    longint lastpid2 = -1;
    Particle *p1;
    Particle *p2;
    Particle *p3;
    Particle *p4;
    for (QuadruplesAngles::const_iterator it = quadruplesAngles.begin();
            it != quadruplesAngles.end(); ++it) {
    
      if (it->first != lastpid2) {
        p2 = storage->lookupRealParticle(it->first);
        if (p2 == NULL) {
          std::stringstream msg;
          msg << "quadruple particle p2 " << it->first << " does not exists here";
          err.setException( msg.str() );
        }
        lastpid2 = it->first;
      }
        
      p1 = storage->lookupLocalParticle(it->second.first.first);
      if(p1 == NULL) {
        std::stringstream msg;
        msg << "quadruple particle p1 " << it->second.first.first << " does not exists here";
        err.setException( msg.str() );
      }
        
      p3 = storage->lookupLocalParticle(it->second.first.second);
      if (p3 == NULL) {
        std::stringstream msg;
        msg << "quadruple particle p3 " << it->second.first.second << " does not exists here";
        err.setException( msg.str() );
      }
      p4 = storage->lookupLocalParticle(it->second.first.third);
      if (p4 == NULL) {
        std::stringstream msg;
        msg << "quadruple particle p4 " << it->second.first.third << " does not exists here";
        err.setException( msg.str() );
      }

      this->add(p1, p2, p3, p4);
    }
    err.checkException();
    
    LOG4ESPP_INFO(theLogger, "regenerated local fixed quadruple list from global list");
  }


  /****************************************************
  ** REGISTRATION WITH PYTHON
  ****************************************************/

  void FixedQuadrupleAngleList::registerPython() {

    using namespace espressopp::python;

    bool (FixedQuadrupleAngleList::*pyAdd)(longint pid1, longint pid2,
           longint pid3, longint pid4) = &FixedQuadrupleAngleList::add;

    class_< FixedQuadrupleAngleList, shared_ptr< FixedQuadrupleAngleList > >
      ("FixedQuadrupleAngleList", init< shared_ptr< storage::Storage > >())
      .def("add", pyAdd)
      .def("size", &FixedQuadrupleAngleList::size)
      .def("getQuadruples",  &FixedQuadrupleAngleList::getQuadruples)
      .def("getQuadruplesAngles",  &FixedQuadrupleAngleList::getQuadruplesAngles)
      .def("getAngle",  &FixedQuadrupleAngleList::getAngle)
     ;
  }
}
