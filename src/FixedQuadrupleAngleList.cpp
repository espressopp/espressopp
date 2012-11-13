#include "python.hpp"
#include <sstream>
#include "FixedQuadrupleAngleList.hpp"
#include <boost/bind.hpp>
#include "storage/Storage.hpp"
#include "bc/BC.hpp"
#include "Buffer.hpp"

namespace espresso {

  LOG4ESPP_LOGGER(FixedQuadrupleAngleList::theLogger, "FixedQuadrupleAngleList");

  FixedQuadrupleAngleList::FixedQuadrupleAngleList(shared_ptr< System > _system) 
    : SystemAccess(_system), quadruplesAngles()
  {
    LOG4ESPP_INFO(theLogger, "construct FixedQuadrupleAngleList");

    System& system = getSystemRef();

    con1 = system.storage->beforeSendParticles.connect
      (boost::bind(&FixedQuadrupleAngleList::beforeSendParticles, this, _1, _2));
    con2 = system.storage->afterRecvParticles.connect
      (boost::bind(&FixedQuadrupleAngleList::afterRecvParticles, this, _1, _2));
    con3 = system.storage->onParticlesChanged.connect
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
    // here we assume pid1 < pid2 < pid3 < pid4
    System& system = getSystemRef();

    // ADD THE LOCAL QUADRUPLET
    Particle *p1 = system.storage->lookupRealParticle(pid1);
    Particle *p2 = system.storage->lookupLocalParticle(pid2);
    Particle *p3 = system.storage->lookupLocalParticle(pid3);
    Particle *p4 = system.storage->lookupLocalParticle(pid4);
    if (!p1)
      // Particle does not exist here, return false
      return false;
    if (!p2) {
      std::stringstream err;
      err << "quadruple particle p2 " << pid2 << " does not exists here and cannot be added";
      std::runtime_error(err.str());
    }
    if (!p3) {
      std::stringstream err;
      err << "quadruple particle p3 " << pid3 << " does not exists here and cannot be added";
      std::runtime_error(err.str());
    }
    if (!p4) {
      std::stringstream err;
      err << "quadruple particle p4 " << pid4 << " does not exists here and cannot be added";
      std::runtime_error(err.str());
    }
    // add the quadruple locally
    this->add(p1, p2, p3, p4);
    //printf("me = %d: pid1 %d, pid2 %d, pid3 %d\n", mpiWorld->rank(), pid1, pid2, pid3);

    Real3D pos1 = p1->position();
    Real3D pos2 = p2->position();
    Real3D pos3 = p3->position();
    Real3D pos4 = p4->position();
    
    Real3D r21 = system.bc->getMinimumImageVector( pos2, pos1);
    Real3D r32 = system.bc->getMinimumImageVector( pos3, pos2);
    Real3D r43 = system.bc->getMinimumImageVector( pos4, pos3);
    
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

    quadruplesAngles.insert(std::make_pair(pid1, 
            std::make_pair(Triple<longint, longint, longint>(pid2, pid3, pid4),
                           acos(cos_phi))));

    LOG4ESPP_INFO(theLogger, "added fixed quadruple to global quadruple list");
    return true;
  }

  python::list FixedQuadrupleAngleList::getQuadruples(){
	python::tuple quadruple;
	python::list quadruples;
	for (QuadruplesAngles::const_iterator it=quadruplesAngles.begin(); it != quadruplesAngles.end(); it++) {
      quadruple = python::make_tuple(it->first, it->second.first.first, it->second.first.second, it->second.first.third);
      quadruples.append(quadruple);
    }

	return quadruples;
  }
  python::list FixedQuadrupleAngleList::getQuadruplesAngles(){
	python::tuple quadruple;
	python::list quadruples;
	for (QuadruplesAngles::const_iterator it=quadruplesAngles.begin(); it != quadruplesAngles.end(); it++) {
      quadruple = python::make_tuple(it->first,
              it->second.first.first,
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
	itr = quadruplesAngles.find(pid1);
	if (itr == quadruplesAngles.end())
		return returnVal; // no elements associated with key, so return immediately

	// get an iterator to the element that is one past the last element associated with key
	lastElement = quadruplesAngles.upper_bound(pid2);
    
    Triple<longint, longint, longint> neededTriple = Triple<longint, longint, longint>(pid2, pid3, pid4);

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
      pid1 = receivedInt[i++];
      n = receivedInt[i++];
      for (; n > 0; --n) {
        pid2 = receivedInt[i++];
        pid3 = receivedInt[i++];
        pid4 = receivedInt[i++];
        angleVal = receivedReal[j++];
        // add the quadruple to the global list
        it = quadruplesAngles.insert(it, std::make_pair(pid1,
              std::make_pair(Triple<longint, longint, longint>(pid2, pid3, pid4),
                             angleVal)));
      }
    }
    if (i != size) {
      printf("ATTETNTION:  recv particles might have read garbage\n");
    }
    LOG4ESPP_INFO(theLogger, "received fixed quadruple list after receive particles");
  }

  void FixedQuadrupleAngleList::onParticlesChanged() {
    System& system = getSystemRef();
    
    // (re-)generate the local quadruple list from the global list
    //printf("FixedQuadrupleAngleList: rebuild local quadruple list from global\n");
    this->clear();
    longint lastpid1 = -1;
    Particle *p1;
    Particle *p2;
    Particle *p3;
    Particle *p4;
    for (QuadruplesAngles::const_iterator it = quadruplesAngles.begin(); it != quadruplesAngles.end(); ++it) {
      //it->first, it->second.first, it->second.second, it->second.third);
      if (it->first != lastpid1) {
        //p1 = system.storage->lookupRealParticle(it->first);
        p1 = system.storage->lookupRealParticle(it->first);
        if(p1 == NULL) {
          std::stringstream err;
          err << "quadruple particle p1 " << it->first << " does not exists here";
          std::runtime_error(err.str());
        }
        lastpid1 = it->first;
      }
      p2 = system.storage->lookupLocalParticle(it->second.first.first);
      if (p2 == NULL) {
        std::stringstream err;
        err << "quadruple particle p2 " << it->second.first.first << " does not exists here";
        std::runtime_error(err.str());
      }
      p3 = system.storage->lookupLocalParticle(it->second.first.second);
      if (p3 == NULL) {
        std::stringstream err;
        err << "quadruple particle p3 " << it->second.first.second << " does not exists here";
        std::runtime_error(err.str());
      }
      p4 = system.storage->lookupLocalParticle(it->second.first.third);
      if (p4 == NULL) {
        std::stringstream err;
        err << "quadruple particle p4 " << it->second.first.third << " does not exists here";
        std::runtime_error(err.str());
      }
      this->add(p1, p2, p3, p4);
    }
    LOG4ESPP_INFO(theLogger, "regenerated local fixed quadruple list from global list");
  }


  /****************************************************
  ** REGISTRATION WITH PYTHON
  ****************************************************/

  void FixedQuadrupleAngleList::registerPython() {

    using namespace espresso::python;

    bool (FixedQuadrupleAngleList::*pyAdd)(longint pid1, longint pid2,
           longint pid3, longint pid4) = &FixedQuadrupleAngleList::add;

    class_< FixedQuadrupleAngleList, shared_ptr< FixedQuadrupleAngleList > >
      ("FixedQuadrupleAngleList", init< shared_ptr< System > >())
      .def("add", pyAdd)
      .def("size", &FixedQuadrupleAngleList::size)
      .def("getQuadruples",  &FixedQuadrupleAngleList::getQuadruples)
      .def("getQuadruplesAngles",  &FixedQuadrupleAngleList::getQuadruplesAngles)
      .def("getAngle",  &FixedQuadrupleAngleList::getAngle)
     ;
  }
}
