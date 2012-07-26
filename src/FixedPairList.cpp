#include "python.hpp"

#include "FixedPairList.hpp"

//#include <vector>
//#include <utility>
//#include <algorithm>
#include <boost/bind.hpp>
#include "storage/Storage.hpp"
#include "Buffer.hpp"

//using namespace std;

namespace espresso {

  /*
  FixedPairList::FixedPairList(shared_ptr< storage::Storage > _storage)
  : FixedListComm (_storage){}

  FixedPairList::~FixedPairList() {
    //std::cout << "~fixedpairlist" << std::endl;
    //FixedListComm::~FixedListComm();
  }*/


  LOG4ESPP_LOGGER(FixedPairList::theLogger, "FixedPairList");


  FixedPairList::FixedPairList(shared_ptr< storage::Storage > _storage) 
    : storage(_storage), globalPairs()
  {
    LOG4ESPP_INFO(theLogger, "construct FixedPairList");

    con1 = storage->beforeSendParticles.connect
      (boost::bind(&FixedPairList::beforeSendParticles, this, _1, _2));
    con2 = storage->afterRecvParticles.connect
      (boost::bind(&FixedPairList::afterRecvParticles, this, _1, _2));
    con3 = storage->onParticlesChanged.connect
      (boost::bind(&FixedPairList::onParticlesChanged, this));
  }

  FixedPairList::~FixedPairList() {

    LOG4ESPP_INFO(theLogger, "~FixedPairList");

    con1.disconnect();
    con2.disconnect();
    con3.disconnect();
  }


  /*
  bool FixedPairList::add(longint pid1, longint pid2) {
    std::vector<longint> tmp;
    tmp.push_back(pid2);
    tmp.push_back(pid1); // this is used as key

    return FixedListComm::add(tmp);
  }*/


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
      globalPairs.insert(std::make_pair(pid1, pid2));
    }
    else {
      // otherwise test whether the pair already exists
      for (GlobalPairs::const_iterator it = equalRange.first; it != equalRange.second; ++it) {
	    if (it->second == pid2) {
	      // TODO: Pair already exists, generate error!
	      ;
	    }
      }
      // if not, insert the new pair
      globalPairs.insert(equalRange.first, std::make_pair(pid1, pid2));
    }
    LOG4ESPP_INFO(theLogger, "added fixed pair to global pair list");
    return true;
  }

  python::list FixedPairList::getBonds()
  {
	python::tuple bond;
	python::list bonds;
	for (GlobalPairs::const_iterator it=globalPairs.begin(); it != globalPairs.end(); it++) {
      bond = python::make_tuple(it->first, it->second);
      bonds.append(bond);
    }

	return bonds;
  }

  void FixedPairList::
  beforeSendParticles(ParticleList& pl, 
		      OutBuffer& buf) {
    std::vector< longint > toSend;
    // loop over the particle list
    for (ParticleList::Iterator pit(pl); pit.isValid(); ++pit) {
      longint pid = pit->id();
      
      LOG4ESPP_DEBUG(theLogger, "send particle with pid " << pid << ", find pairs");

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
             it != equalRange.second; ++it) {
          toSend.push_back(it->second);
          LOG4ESPP_DEBUG(theLogger, "send global bond: pid "
                       << pid << " and partner " << it->second);
        }

        // delete all of these pairs from the global list
        globalPairs.erase(equalRange.first, equalRange.second);
      }
    }
    // send the list
    buf.write(toSend);
    LOG4ESPP_INFO(theLogger, "prepared fixed pair list before send particles");
  }

  void FixedPairList::
  afterRecvParticles(ParticleList &pl, 
		     InBuffer& buf) {
    std::vector< longint > received;
    int n;
    longint pid1, pid2;
    GlobalPairs::iterator it = globalPairs.begin();
    // receive the bond list
    buf.read(received);
    int size = received.size(); int i = 0;
    while (i < size) {
      // unpack the list
      pid1 = received[i++];
      n = received[i++];
      LOG4ESPP_DEBUG(theLogger, "recv particle " << pid1 << 
                                ", has " << n << " global pairs");
      for (; n > 0; --n) {
	pid2 = received[i++];
	// add the bond to the global list
        LOG4ESPP_DEBUG(theLogger, "received pair " << pid1 << " , " << pid2);
	it = globalPairs.insert(it, std::make_pair(pid1, pid2));
      }
    }
    if (i != size) {
      LOG4ESPP_ERROR(theLogger, 
        "ATTETNTION:  recv particles might have read garbage\n");
    }
    LOG4ESPP_INFO(theLogger, "received fixed pair list after receive particles");
  }

  void FixedPairList::
  onParticlesChanged() {
    LOG4ESPP_INFO(theLogger, "rebuild local bond list from global\n");
    this->clear();
    longint lastpid1 = -1;
    Particle *p1;
    Particle *p2;
    for (GlobalPairs::const_iterator it = globalPairs.begin();
	 it != globalPairs.end(); ++it) {
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
    LOG4ESPP_INFO(theLogger, "regenerated local fixed pair list from global list");
  }

  /****************************************************
  ** REGISTRATION WITH PYTHON
  ****************************************************/

  void FixedPairList::registerPython() {

    using namespace espresso::python;

    bool (FixedPairList::*pyAdd)(longint pid1, longint pid2)
      = &FixedPairList::add;
    //bool (FixedPairList::*pyAdd)(pvec pids) = &FixedPairList::add;

    class_<FixedPairList, shared_ptr<FixedPairList> >
      ("FixedPairList", init <shared_ptr<storage::Storage> >())
      .def("add", pyAdd)
      .def("size", &FixedPairList::size)
      .def("getBonds",  &FixedPairList::getBonds)
      ;
  }
}
