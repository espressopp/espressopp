#include "python.hpp"
#define LOG4ESPP_LEVEL_DEBUG

#include "FixedPairList.hpp"

#include <vector>
#include <utility>
#include <algorithm>
#include <boost/bind.hpp>
#include "storage/Storage.hpp"

using namespace std;

namespace espresso {

  LOG4ESPP_LOGGER(FixedPairList::theLogger, "FixedPairList");

  FixedPairList::FixedPairList(shared_ptr< storage::Storage > _storage) 
    : storage(_storage), globalPairs()
  {
    LOG4ESPP_INFO(theLogger, "construct FixedPairList");

    con1 = storage->beforeSendParticles.connect
      (boost::bind(&FixedPairList::beforeSendParticles, this, _1, _2));
    con2 = storage->afterRecvParticles.connect
      (boost::bind(&FixedPairList::afterRecvParticles, this, _1, _2));
    con3 = storage->onResortParticles.connect
      (boost::bind(&FixedPairList::onResortParticles, this));
  }

  FixedPairList::~FixedPairList() {

    LOG4ESPP_INFO(theLogger, "~FixedPairList");

    con1.disconnect();
    con2.disconnect();
    con3.disconnect();
  }

  bool FixedPairList::
  add(longint pid1, longint pid2) {
    if (pid1 > pid2)
      std::swap(pid1, pid2);

    // ADD THE LOCAL PAIR
    Particle *p1 = storage->lookupRealParticle(pid1);
    Particle *p2 = storage->lookupLocalParticle(pid2);
    if (!p1)
      // Particle does not exist here, return false
      return false;
    if (!p2)
      // TODO: Second particle does not exist here, throw exception!
      return false;
    // add the pair locally
    this->add(p1, p2);

    // ADD THE GLOBAL PAIR
    // see whether the particle already has pairs
    std::pair<GlobalPairs::const_iterator, 
      GlobalPairs::const_iterator> equalRange 
      = globalPairs.equal_range(pid1);
    if (equalRange.first == globalPairs.end()) {
      // if it hasn't, insert the new pair
      globalPairs.insert(make_pair(pid1, pid2));
    }
    else {
      // otherwise test whether the pair already exists
      for (GlobalPairs::const_iterator it = equalRange.first;
	   it != equalRange.second; ++it)
	if (it->second == pid2)
	  // TODO: Pair already exists, generate error!
	  ;
      // if not, insert the new pair
      globalPairs.insert(equalRange.first, make_pair(pid1, pid2));
    }
    LOG4ESPP_INFO(theLogger, "added fixed pair to global pair list");
    return true;
  }

  void FixedPairList::
  beforeSendParticles(ParticleList& pl, 
		      mpi::packed_oarchive& ar) {
    vector< longint > toSend;
    // loop over the particle list
    for (ParticleList::Iterator pit(pl);
	 pit.isValid(); ++pit) {
      longint pid = pit->p.id;
      
      // find all pairs that involve this particle
      
      int n = globalPairs.count(pid);
      if (n > 0) {
	std::pair<GlobalPairs::const_iterator, 
	  GlobalPairs::const_iterator> equalRange 
	  = globalPairs.equal_range(pid);

	// first write the pid of the first particle
	// then the number of partners
	// and then the pids of the partners
	toSend.reserve(toSend.size()+n+1);
	toSend.push_back(pid);
	toSend.push_back(n);
	for (GlobalPairs::const_iterator it = equalRange.first; 
	     it != equalRange.second; ++it)
	  toSend.push_back(it->second);

	// delete all of these pairs from the global list
	globalPairs.erase(equalRange.first, equalRange.second);
      }
    }
    // send the list
    ar << toSend;
    LOG4ESPP_INFO(theLogger, "prepared fixed pair list before send particles");
  }

  void FixedPairList::
  afterRecvParticles(ParticleList &pl, 
		     mpi::packed_iarchive& ar) {
    vector< longint > received;
    longint n;
    longint pid1, pid2;
    GlobalPairs::iterator it = globalPairs.begin();
    // receive the bond list
    ar >> received;
    for (unsigned int i = 0; i < received.size(); ++i) {
      // unpack the list
      pid1 = received[i++];
      n = received[i++];
      for (; n > 0; --n) {
	pid2 = received[i++];
	// add the bond to the global list
	it = globalPairs.insert(it, make_pair(pid1, pid2));
      }
    }
    LOG4ESPP_INFO(theLogger, "received fixed pair list after receive particles");
  }

  void FixedPairList::
  onResortParticles() {
    // (re-)generate the local bond list from the global list
    this->clear();
    longint lastpid1 = -1;
    Particle *p1;
    Particle *p2;
    for (GlobalPairs::const_iterator it = globalPairs.begin();
	 it != globalPairs.end(); ++it) {
      if (it->first != lastpid1) {
	p1 = storage->lookupRealParticle(it->first);
	lastpid1 = it->first;
      }
      p2 = storage->lookupLocalParticle(it->second);
      this->add(p1, p2);
    }
    LOG4ESPP_INFO(theLogger, "regenerated local fixed pair list from global list");
  }

  /****************************************************
  ** REGISTRATION WITH PYTHON
  ****************************************************/

  void FixedPairList::registerPython() {

    using namespace espresso::python;

    bool (FixedPairList::*pyAdd)(longint pid1, longint pid2) 
      = &FixedPairList::add;

    class_< FixedPairList, shared_ptr< FixedPairList > >
      ("FixedPairList", init< shared_ptr< storage::Storage > >())
      .def("add", pyAdd)
      ;
  }
}
