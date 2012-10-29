#include "python.hpp"
#include <sstream>
#include "FixedTripleCosList.hpp"

//#include <vector>
//#include <utility>
//#include <algorithm>
#include <boost/bind.hpp>
#include "storage/Storage.hpp"
#include "Buffer.hpp"

namespace espresso {

  /*
  FixedTripleCosList::FixedTripleCosList(shared_ptr< storage::Storage > _storage)
  : FixedListComm (_storage){}
  */


  LOG4ESPP_LOGGER(FixedTripleCosList::theLogger, "FixedTripleCosList");

  FixedTripleCosList::FixedTripleCosList(shared_ptr< storage::Storage > _storage)
    : storage(_storage), triplesCos()
  {
    LOG4ESPP_INFO(theLogger, "construct FixedTripleCosList");

    con1 = storage->beforeSendParticles.connect
      (boost::bind(&FixedTripleCosList::beforeSendParticles, this, _1, _2));
    con2 = storage->afterRecvParticles.connect
      (boost::bind(&FixedTripleCosList::afterRecvParticles, this, _1, _2));
    con3 = storage->onParticlesChanged.connect
      (boost::bind(&FixedTripleCosList::onParticlesChanged, this));
  }

  FixedTripleCosList::~FixedTripleCosList() {

    LOG4ESPP_INFO(theLogger, "~FixedTripleCosList");

    con1.disconnect();
    con2.disconnect();
    con3.disconnect();
  }

  bool FixedTripleCosList::
  add(longint pid1, longint pid2, longint pid3) {
    // ADD THE LOCAL TRIPLET
    Particle *p1 = storage->lookupLocalParticle(pid1);
    Particle *p2 = storage->lookupRealParticle(pid2);
    Particle *p3 = storage->lookupLocalParticle(pid3);

    // middle particle is the reference particle and must exist here
    if (!p2)
      // particle does not exists here (some other CPU must have it)
      return false;


    if (!p1) {
      std::stringstream err;
      err << "triple particle p1 " << pid1 << " does not exists here and cannot be added";
      std::runtime_error(err.str());
    }

    if (!p3) {
      std::stringstream err;
      err << "triple particle p3 " << pid1 << " does not exists here and cannot be added";
      std::runtime_error(err.str());
    }

    // add the triple locally
    this->add(p1, p2, p3);
    //printf("me = %d: pid1 %d, pid2 %d, pid3 %d\n", mpiWorld->rank(), pid1, pid2, pid3);
    
    Real3D pos1 = p1->position();
    Real3D pos2 = p2->position();
    Real3D pos3 = p3->position();
    
    Real3D r12 = pos2 - pos1;
    Real3D r32 = pos2 - pos3;
    
    real cosVal = r12*r32 / r12.abs() / r32.abs();

    // ADD THE GLOBAL TRIPLET
    triplesCos.insert(std::make_pair(pid2, std::make_pair(std::make_pair(pid1, pid3), cosVal)));
    LOG4ESPP_INFO(theLogger, "added fixed triple to global triple list");
    return true;
  }

  python::list FixedTripleCosList::getTriples(){
	python::tuple triple;
	python::list triples;
	for (TriplesCos::const_iterator it=triplesCos.begin(); it!=triplesCos.end(); it++) {
      triple = python::make_tuple(it->first, it->second.first.first, it->second.first.second);
      triples.append(triple);
    }

	return triples;
  }
  
  real FixedTripleCosList::getCos(int pid1, int pid2, int pid3){
    real returnVal = -3;
    
    TriplesCos::iterator itr;
	TriplesCos::iterator lastElement;
	
	// locate an iterator to the first pair object associated with key
	itr = triplesCos.find(pid2);
	if (itr == triplesCos.end())
		return returnVal; // no elements associated with key, so return immediately

	// get an iterator to the element that is one past the last element associated with key
	lastElement = triplesCos.upper_bound(pid2);
    
    std::pair<longint, longint> neededPair = std::make_pair(pid1, pid3);

	for ( ; itr != lastElement; ++itr){
      if(neededPair==itr->second.first){
        returnVal = itr->second.second;
        break;
      }
    }
    
	return returnVal;
  }
  
  void FixedTripleCosList::setCos(int pid1, int pid2, int pid3){
    TriplesCos::iterator itr;
	TriplesCos::iterator lastElement;
	
	// locate an iterator to the first pair object associated with key
	itr = triplesCos.find(pid2);
	if (itr == triplesCos.end())
      std::cout<<"Warning! triplesCos was not modified because the triple "
              "doesn't exists here"<<std::endl;
      return; // no elements associated with key, so return immediately

	// get an iterator to the element that is one past the last element associated with key
	lastElement = triplesCos.upper_bound(pid2);

    std::pair<longint, longint> neededPair = std::make_pair(pid1, pid3);
	for ( ; itr != lastElement; ++itr){
      
      if(neededPair==itr->second.first){
        Particle *p1 = storage->lookupLocalParticle(pid1);
        Particle *p2 = storage->lookupRealParticle(pid2);
        Particle *p3 = storage->lookupLocalParticle(pid3);

        Real3D pos1 = p1->position();
        Real3D pos2 = p2->position();
        Real3D pos3 = p3->position();

        Real3D r12 = pos2 - pos1;
        Real3D r32 = pos2 - pos3;

        real cosVal = r12*r32 / r12.abs() / r32.abs();
        
        itr->second.second = cosVal;

        break;
      }
    }
  }

  void FixedTripleCosList::beforeSendParticles(ParticleList& pl, OutBuffer& buf) {
    std::vector< longint > toSend;
    // loop over the particle list
    for (ParticleList::Iterator pit(pl); pit.isValid(); ++pit) {
      longint pid = pit->id();
      
      // LOG4ESPP_DEBUG(theLogger, "send particle with pid " << pid << ", find triples");
      //printf ("me = %d: send particle with pid %d find triples\n", mpiWorld->rank(), pid);

      // find all triples that involve this particle
      int n = triplesCos.count(pid);
      //printf ("me = %d: send particle with pid %d, has %d global triples\n", 
                //mpiWorld->rank(), pid, n);

      if (n > 0) {
        std::pair<TriplesCos::const_iterator,
        TriplesCos::const_iterator> equalRange = triplesCos.equal_range(pid);

        // first write the pid of this particle
        // then the number of partners (n)
        // and then the pids of the partners
        toSend.reserve(toSend.size()+2*n+1);
        toSend.push_back(pid);
        toSend.push_back(n);
        for(TriplesCos::const_iterator it = equalRange.first; it!=equalRange.second; ++it) {
          toSend.push_back(it->second.first.first);
          toSend.push_back(it->second.first.second);
        }

        // delete all of these triples from the global list
        triplesCos.erase(pid);
      }
    }
    // send the list
    buf.write(toSend);
    LOG4ESPP_INFO(theLogger, "prepared fixed triple list before send particles");
  }

  void FixedTripleCosList::afterRecvParticles(ParticleList &pl, InBuffer& buf) {

    std::vector< longint > received;
    int n;
    longint pid1, pid2, pid3;
    TriplesCos::iterator it = triplesCos.begin();
    // receive the triple list
    buf.read(received);
    int size = received.size(); int i = 0;
    while (i < size) {
      // unpack the list
      pid2 = received[i++];
      n = received[i++];
      //printf ("me = %d: recv particle with pid %d, has %d global triples\n",
                //mpiWorld->rank(), pid1, n);
      for (; n > 0; --n) {
        pid1 = received[i++];
        pid3 = received[i++];
        // add the triple
        it = triplesCos.insert(it, std::make_pair(pid2, std::make_pair(std::make_pair(pid1, pid3), 0.0)));
      }
    }
    if (i != size) {
      printf("ATTETNTION:  recv particles might have read garbage\n");
    }
    LOG4ESPP_INFO(theLogger, "received fixed triple list after receive particles");
  }

  void FixedTripleCosList::onParticlesChanged() {
    // (re-)generate the local triple list from the global list
    //printf("FixedTripleCosList: rebuild local triple list from global\n");
    this->clear();
    longint lastpid2 = -1;
    Particle *p1;
    Particle *p2;
    Particle *p3;
    for (TriplesCos::const_iterator it = triplesCos.begin(); it!=triplesCos.end(); ++it) {
      //printf("lookup global triple %d %d %d\n", it->first, it->second.first, it->second.second);
      if (it->first != lastpid2) {
        p2 = storage->lookupRealParticle(it->first);
        if (p2 == NULL) {
          std::stringstream err;
          err << "triple particle p2 " << it->first << " does not exists here";
          std::runtime_error(err.str());
        }
	    lastpid2 = it->first;
      }
      p1 = storage->lookupLocalParticle(it->second.first.first);
      if (p1 == NULL) {
        std::stringstream err;
        err << "triple particle p1 " << it->second.first.first << " does not exists here";
        std::runtime_error(err.str());
      }
      p3 = storage->lookupLocalParticle(it->second.first.second);
      if (p3 == NULL) {
        std::stringstream err;
        err << "triple particle p3 " << it->second.first.second << " does not exists here";
        std::runtime_error(err.str());
      }
      this->add(p1, p2, p3);
      
      setCos(it->second.first.first,it->first,it->second.first.second);
    }
    LOG4ESPP_INFO(theLogger, "regenerated local fixed triple list from global list");
  }

  /****************************************************
  ** REGISTRATION WITH PYTHON
  ****************************************************/

  void FixedTripleCosList::registerPython() {

    using namespace espresso::python;

    bool (FixedTripleCosList::*pyAdd)(longint pid1, longint pid2, longint pid3)
      = &FixedTripleCosList::add;
    //bool (FixedTripleCosList::*pyAdd)(pvec pids)
    //      = &FixedTripleCosList::add;

    class_< FixedTripleCosList, shared_ptr< FixedTripleCosList > >
      ("FixedTripleCosList", init< shared_ptr< storage::Storage > >())
      .def("add", pyAdd)
      .def("size", &FixedTripleCosList::size)
      .def("getTriples",  &FixedTripleCosList::getTriples)
      .def("getCos",  &FixedTripleCosList::getCos)
     ;
  }
}
