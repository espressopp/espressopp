#include "python.hpp"
#define LOG4ESPP_LEVEL_DEBUG

#include "FixedQuadrupleList.hpp"

#include <vector>
#include <utility>
#include <algorithm>
#include <boost/bind.hpp>
#include "storage/Storage.hpp"
#include "Buffer.hpp"

using namespace std;

namespace espresso {

  LOG4ESPP_LOGGER(FixedQuadrupleList::theLogger, "FixedQuadrupleList");

  FixedQuadrupleList::FixedQuadrupleList(shared_ptr< storage::Storage > _storage) 
    : storage(_storage), globalQuadruples()
  {
    LOG4ESPP_INFO(theLogger, "construct FixedQuadrupleList");

    con1 = storage->beforeSendParticles.connect
      (boost::bind(&FixedQuadrupleList::beforeSendParticles, this, _1, _2));
    con2 = storage->afterRecvParticles.connect
      (boost::bind(&FixedQuadrupleList::afterRecvParticles, this, _1, _2));
    con3 = storage->onParticlesChanged.connect
      (boost::bind(&FixedQuadrupleList::onParticlesChanged, this));
  }

  FixedQuadrupleList::~FixedQuadrupleList() {

    LOG4ESPP_INFO(theLogger, "~FixedQuadrupleList");

    con1.disconnect();
    con2.disconnect();
    con3.disconnect();
  }

  bool FixedQuadrupleList::
  add(longint pid1, longint pid2, longint pid3, longint pid4) {
    // here we assume pid1 < pid2 < pid3 < pid4

    // ADD THE LOCAL QUADRUPLET
    Particle *p1 = storage->lookupRealParticle(pid1);
    Particle *p2 = storage->lookupLocalParticle(pid2);
    Particle *p3 = storage->lookupLocalParticle(pid3);
    Particle *p4 = storage->lookupLocalParticle(pid4);
    if (!p1)
      // Particle does not exist here, return false
      return false;
    if (!p2)
      // TODO: Second particle does not exist here, throw exception! (Why? JH)
      return false;
    if (!p3)
      // TODO: Third particle does not exist here, throw exception! (Why? JH)
      return false;
    if (!p4)
      // TODO: Fourth particle does not exist here, throw exception! (Why? JH)
      return false;
    // add the quadruple locally
    this->add(p1, p2, p3, p4);
    //printf("me = %d: pid1 %d, pid2 %d, pid3 %d\n", mpiWorld->rank(), pid1, pid2, pid3);

    // ADD THE GLOBAL PAIR
    // see whether the particle already has quadruples
    std::pair<GlobalQuadruples::const_iterator,
              GlobalQuadruples::const_iterator> equalRange 
      = globalQuadruples.equal_range(pid1);
    if (equalRange.first == globalQuadruples.end()) {
      // if it hasn't, insert the new quadruple
      globalQuadruples.insert(std::make_pair(pid1,
        Triple<longint, longint, longint>(pid2, pid3, pid4)));
    }
    else {
      // otherwise test whether the quadruple already exists
      for (GlobalQuadruples::const_iterator it = equalRange.first;
	   it != equalRange.second; ++it)
	if (it->second == Triple<longint, longint, longint>(pid2, pid3, pid4))
	  // TODO: Quadruple already exists, generate error!
	  ;
      // if not, insert the new quadruple
      globalQuadruples.insert(equalRange.first,
        std::make_pair(pid1, Triple<longint, longint, longint>(pid2, pid3, pid4)));
    }
    LOG4ESPP_INFO(theLogger, "added fixed quadruple to global quadruple list");
    return true;
  }

  void FixedQuadrupleList::
  beforeSendParticles(ParticleList& pl, OutBuffer& buf) {
    
    vector< longint > toSend;
    // loop over the particle list
    for (ParticleList::Iterator pit(pl); pit.isValid(); ++pit) {
      longint pid = pit->p.id;
      
      // LOG4ESPP_DEBUG(theLogger, "send particle with pid " << pid << ", find quadruples");
      //printf ("me = %d: send particle with pid %d find quadruples\n", mpiWorld->rank(), pid);

      // find all quadruples that involve this particle
      int n = globalQuadruples.count(pid);
      //printf ("me = %d: send particle with pid %d, has %d global quadruples\n", 
                //mpiWorld->rank(), pid, n);

      if (n > 0) {
	std::pair<GlobalQuadruples::const_iterator, 
	  GlobalQuadruples::const_iterator> equalRange 
	  = globalQuadruples.equal_range(pid);

	// first write the pid of this particle
	// then the number of partners (n)
	// and then the pids of the partners
	toSend.reserve(toSend.size()+3*n+1);
	toSend.push_back(pid);
	toSend.push_back(n);
	for (GlobalQuadruples::const_iterator it = equalRange.first;
	     it != equalRange.second; ++it) {
	  toSend.push_back(it->second.first);
	  toSend.push_back(it->second.second);
	  toSend.push_back(it->second.third);
          //printf ("send global quadruple: pid %d and partner %d\n", pid, it->second.first);
          //printf ("send global quadruple: pid %d and partner %d\n", pid, it->second.second);
          //printf ("send global quadruple: pid %d and partner %d\n", pid, it->second.third);
        }

	// delete all of these quadruples from the global list
	globalQuadruples.erase(equalRange.first, equalRange.second);
      }
    }
    // send the list
    buf.write(toSend);
    LOG4ESPP_INFO(theLogger, "prepared fixed quadruple list before send particles");
  }

  void FixedQuadrupleList::
  afterRecvParticles(ParticleList &pl, InBuffer& buf) {

    vector< longint > received;
    int n;
    longint pid1, pid2, pid3, pid4;
    GlobalQuadruples::iterator it = globalQuadruples.begin();
    // receive the quadruple list
    buf.read(received);
    int size = received.size(); int i = 0;
    while (i < size) {
      // unpack the list
      pid1 = received[i++];
      n = received[i++];
      //printf ("me = %d: recv particle with pid %d, has %d global quadruples\n",
                //mpiWorld->rank(), pid1, n);
      for (; n > 0; --n) {
	pid2 = received[i++];
	pid3 = received[i++];
	pid4 = received[i++];
	// add the quadruple to the global list
        //printf("received quadruple %d %d %d %d, add quadruple to global list\n", pid1, pid2, pid3, pid4);
	it = globalQuadruples.insert(it, std::make_pair(pid1,
          Triple<longint, longint, longint>(pid2, pid3, pid4)));
      }
    }
    if (i != size) {
      printf("ATTETNTION:  recv particles might have read garbage\n");
    }
    LOG4ESPP_INFO(theLogger, "received fixed quadruple list after receive particles");
  }

  void FixedQuadrupleList::onParticlesChanged() {
    
    // (re-)generate the local quadruple list from the global list
    //printf("FixedQuadrupleList: rebuild local quadruple list from global\n");
    this->clear();
    longint lastpid1 = -1;
    Particle *p1;
    Particle *p2;
    Particle *p3;
    Particle *p4;
    for (GlobalQuadruples::const_iterator it = globalQuadruples.begin();
	 it != globalQuadruples.end(); ++it) {
      //printf("lookup global quadruple %d %d %d %d\n",
        //it->first, it->second.first, it->second.second, it->second.third);
      if (it->first != lastpid1) {
	p1 = storage->lookupRealParticle(it->first);
        if (p1 == NULL) {
          printf("SERIOUS ERROR: particle %d not available\n", it->first);
        }
	lastpid1 = it->first;
      }
      p2 = storage->lookupLocalParticle(it->second.first);
      if (p2 == NULL) {
        printf("SERIOUS ERROR: 2nd particle %d not available\n", it->second.first);
      }
      p3 = storage->lookupLocalParticle(it->second.second);
      if (p3 == NULL) {
        printf("SERIOUS ERROR: 3rd particle %d not available\n", it->second.second);
      }
      p4 = storage->lookupLocalParticle(it->second.third);
      if (p4 == NULL) {
        printf("SERIOUS ERROR: 4th particle %d not available\n", it->second.third);
      }
      this->add(p1, p2, p3, p4);
    }
    LOG4ESPP_INFO(theLogger, "regenerated local fixed quadruple list from global list");
  }

  /****************************************************
  ** REGISTRATION WITH PYTHON
  ****************************************************/

  void FixedQuadrupleList::registerPython() {

    using namespace espresso::python;

    bool (FixedQuadrupleList::*pyAdd)(longint pid1, longint pid2,
           longint pid3, longint pid4) = &FixedQuadrupleList::add;

    class_< FixedQuadrupleList, shared_ptr< FixedQuadrupleList > >
      ("FixedQuadrupleList", init< shared_ptr< storage::Storage > >())
      .def("add", pyAdd)
      ;
  }
}
