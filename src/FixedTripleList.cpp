#include "python.hpp"
#define LOG4ESPP_LEVEL_DEBUG

#include "FixedTripleList.hpp"

#include <vector>
#include <utility>
#include <algorithm>
#include <boost/bind.hpp>
#include "storage/Storage.hpp"

using namespace std;

namespace espresso {

  LOG4ESPP_LOGGER(FixedTripleList::theLogger, "FixedTripleList");

  FixedTripleList::FixedTripleList(shared_ptr< storage::Storage > _storage) 
    : storage(_storage), globalTriples()
  {
    LOG4ESPP_INFO(theLogger, "construct FixedTripleList");

    con1 = storage->beforeSendParticles.connect
      (boost::bind(&FixedTripleList::beforeSendParticles, this, _1, _2));
    con2 = storage->afterRecvParticles.connect
      (boost::bind(&FixedTripleList::afterRecvParticles, this, _1, _2));
    con3 = storage->onParticlesChanged.connect
      (boost::bind(&FixedTripleList::onParticlesChanged, this));
  }

  FixedTripleList::~FixedTripleList() {

    LOG4ESPP_INFO(theLogger, "~FixedTripleList");

    con1.disconnect();
    con2.disconnect();
    con3.disconnect();
  }

  bool FixedTripleList::
  add(longint pid1, longint pid2, longint pid3) {
    if (pid3 < pid2)
      std::swap(pid2, pid3);
    if (pid2 < pid1)
      std::swap(pid1, pid2); 
    if (pid3 < pid2)
      std::swap(pid2, pid3);

    // ADD THE LOCAL TRIPLET
    Particle *p1 = storage->lookupRealParticle(pid1);
    Particle *p2 = storage->lookupLocalParticle(pid2);
    Particle *p3 = storage->lookupLocalParticle(pid3);
    if (!p1)
      // Particle does not exist here, return false
      return false;
    if (!p2)
      // TODO: Second particle does not exist here, throw exception! (Why? JH)
      return false;
    if (!p3)
      // TODO: Third particle does not exist here, throw exception! (Why? JH)
      return false;
    // add the triple locally
    this->add(p1, p2, p3);
    printf("me = %d: pid1 %d, pid2 %d, pid3 %d\n", mpiWorld->rank(), pid1, pid2, pid3);

    // ADD THE GLOBAL PAIR
    // see whether the particle already has triples
    std::pair<GlobalTriples::const_iterator, GlobalTriples::const_iterator> equalRange 
      = globalTriples.equal_range(pid1);
    if (equalRange.first == globalTriples.end()) {
      // if it hasn't, insert the new triple
      globalTriples.insert(std::make_pair(pid1, std::pair<longint, longint>(pid2, pid3)));
    }
    else {
      // otherwise test whether the triple already exists
      for (GlobalTriples::const_iterator it = equalRange.first;
	   it != equalRange.second; ++it)
	if (it->second == std::pair<longint, longint>(pid2, pid3))
	  // TODO: Triple already exists, generate error!
	  ;
      // if not, insert the new triple
      globalTriples.insert(equalRange.first, std::make_pair(pid1, std::pair<longint, longint>(pid2, pid3)));
    }
    LOG4ESPP_INFO(theLogger, "added fixed triple to global triple list");
    return true;
  }

  void FixedTripleList::
  beforeSendParticles(ParticleList& pl, 
		      mpi::packed_oarchive& ar) {
    /*
    vector< longint > toSend;
    // loop over the particle list
    for (ParticleList::Iterator pit(pl);
	 pit.isValid(); ++pit) {
      longint pid = pit->p.id;
      
      // LOG4ESPP_DEBUG(theLogger, "send particle with pid " << pid << ", find triples");
      printf ("me = %d: send particle with pid %d find triples\n", mpiWorld->rank(), pid);

      // find all triples that involve this particle
      
      int n = globalTriples.count(pid);
      printf ("me = %d: send particle with pid %d, has %d global triples\n", 
                mpiWorld->rank(), pid, n);

      if (n > 0) {
	triple<GlobalTriples::const_iterator, 
	  GlobalTriples::const_iterator, GlobalTriples::const_iterator> equalRange 
	  = globalTriples.equal_range(pid);

	// first write the pid of the first particle
	// then the number of partners
	// and then the pids of the partners
	toSend.reserve(toSend.size()+n+1);
	toSend.push_back(pid);
	toSend.push_back(n);
	for (GlobalPairs::const_iterator it = equalRange.first; 
	     it != equalRange.second; ++it) {
	  toSend.push_back(it->second);
          printf ("send global triple: pid %d and partner %d\n", pid, it->second);
        }

	// delete all of these triples from the global list
	globalTriples.erase(equalRange.first, equalRange.second, equalRange.third);
      }
    }
    // send the list
    ar << toSend;
    LOG4ESPP_INFO(theLogger, "prepared fixed triple list before send particles");
    */
  }

  void FixedTripleList::
  afterRecvParticles(ParticleList &pl, 
		     mpi::packed_iarchive& ar) {
    /*
    vector< longint > received;
    longint n;
    longint pid1, pid2;
    GlobalTriples::iterator it = globalTriples.begin();
    // receive the triple list
    ar >> received;
    int size = received.size(); int i = 0;
    while (i < size) {
      // unpack the list
      pid1 = received[i++];
      n = received[i++];
      printf ("me = %d: recv particle with pid %d, has %d global triples\n",
                mpiWorld->rank(), pid1, n);
      for (; n > 0; --n) {
	pid2 = received[i++];
	// add the triple to the global list
        printf("received triple %d %d, add triple to global list\n", pid1, pid2);
	it = globalTriples.insert(it, triple(pid1, pid2, pid3));
      }
    }
    if (i != size) {
      printf("ATTETNTION:  recv particles might have read garbage\n");
    }
    LOG4ESPP_INFO(theLogger, "received fixed triple list after receive particles");
    */
  }

  void FixedTripleList::
  onParticlesChanged() {
    /*
    // (re-)generate the local triple list from the global list
    printf("FixedTripleList: rebuild local triple list from global\n");
    this->clear();
    longint lastpid1 = -1;
    Particle *p1;
    Particle *p2;
    Particle *p3;
    for (GlobalPairs::const_iterator it = globalTriples.begin();
	 it != globalTriples.end(); ++it) {
      printf("lookup global triple %d %d\n", it->first, it->second);
      if (it->first != lastpid1) {
	p1 = storage->lookupRealParticle(it->first);
        if (p1 == NULL) {
          printf("SERIOUS ERROR: particle %d not available\n", it->first);
        }
	lastpid1 = it->first;
      }
      p2 = storage->lookupLocalParticle(it->second);
      if (p2 == NULL) {
        printf("SERIOUS ERROR: 2nd particle %d not available\n", it->second);
      }
      this->add(p1, p2);
    }
    LOG4ESPP_INFO(theLogger, "regenerated local fixed triple list from global list");
    */
  }

  /****************************************************
  ** REGISTRATION WITH PYTHON
  ****************************************************/

  void FixedTripleList::registerPython() {

    using namespace espresso::python;

    bool (FixedTripleList::*pyAdd)(longint pid1, longint pid2, longint pid3) 
      = &FixedTripleList::add;

    class_< FixedTripleList, shared_ptr< FixedTripleList > >
      ("FixedTripleList", init< shared_ptr< storage::Storage > >())
      .def("add", pyAdd)
      ;
  }
}
