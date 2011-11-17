#include "python.hpp"

#include "FixedTripleListAdress.hpp"

//#include <vector>
//#include <utility>
//#include <algorithm>
#include <boost/bind.hpp>
#include "storage/Storage.hpp"
#include "Buffer.hpp"

//using namespace std;

namespace espresso {

  LOG4ESPP_LOGGER(FixedTripleListAdress::theLogger, "FixedTripleListAdress");


  FixedTripleListAdress::
  FixedTripleListAdress(shared_ptr< storage::Storage > _storage,
          shared_ptr<FixedTupleList> _fixedtupleList)
          : FixedTripleList(_storage), fixedtupleList(_fixedtupleList) {
    LOG4ESPP_INFO(theLogger, "construct FixedTripleListAdress");

    con = fixedtupleList->beforeSendATParticles.connect
          (boost::bind(&FixedTripleListAdress::beforeSendATParticles, this, _1, _2));
  }

  FixedTripleListAdress::~FixedTripleListAdress() {

    LOG4ESPP_INFO(theLogger, "~FixedTripleListAdress");

    con.disconnect();
  }


  // override parent function (use lookupAdrATParticle)
  bool FixedTripleListAdress::add(longint pid1, longint pid2, longint pid3) {

    if (pid1 > pid2)
      std::swap(pid1, pid2);

    // ADD THE LOCAL TRIPLES
    Particle *p1 = storage->lookupAdrATParticle(pid1);
    Particle *p2 = storage->lookupAdrATParticle(pid2);
    Particle *p3 = storage->lookupAdrATParticle(pid3);
    if (!p1)
      // Particle does not exist here, return false
      return false;
    if (!p2)
      // TODO: Second particle does not exist here, throw exception!
      return false;
    if (!p3)
      // TODO: Third particle does not exist here, throw exception!
      return false;
    // add the triple locally
    this->add(p1, p2, p3);

    // ADD THE GLOBAL TRIPLET
    // see whether the particle already has triples
    std::pair<GlobalTriples::const_iterator,
              GlobalTriples::const_iterator> equalRange
      = globalTriples.equal_range(pid2);
    if (equalRange.first == globalTriples.end()) {
      // if it hasn't, insert the new triple
      globalTriples.insert(std::make_pair(pid2,
                           std::pair<longint, longint>(pid1, pid3)));
    }
    else {
      // otherwise test whether the triple already exists
      for (GlobalTriples::const_iterator it = equalRange.first;
           it != equalRange.second; ++it)
        if (it->second == std::pair<longint, longint>(pid1, pid3))
        // TODO: Triple already exists, generate error!
      ;
      // if not, insert the new triple
      globalTriples.insert(equalRange.first, std::make_pair(pid2,
                           std::pair<longint, longint>(pid1, pid3)));
    }
    LOG4ESPP_INFO(theLogger, "added fixed pair to global pair list");
    return true;
  }


  void FixedTripleListAdress::beforeSendATParticles(std::vector<longint>& atpl,
          OutBuffer& buf) {

        //std::cout << "beforeSendATParticles() fixed pl (size " << atpl.size() << ")\n";

        std::vector< longint > toSend;

        // loop over the VP particle list
        for (std::vector<longint>::iterator it = atpl.begin();
                it != atpl.end(); ++it) {
          longint pid = *it;

          LOG4ESPP_DEBUG(theLogger, "send particle with pid " << pid << ", find pairs");

          // find all triples that involve this particle
          int n = globalTriples.count(pid);

          if (n > 0) {
            std::pair<GlobalTriples::const_iterator,
            GlobalTriples::const_iterator> equalRange
                = globalTriples.equal_range(pid);

            // first write the pid of this particle
            // then the number of partners (n)
            // and then the pids of the partners
            toSend.reserve(toSend.size()+2*n+1);
            toSend.push_back(pid);
            toSend.push_back(n);
            for (GlobalTriples::const_iterator it = equalRange.first;
                    it != equalRange.second; ++it) {
                toSend.push_back(it->second.first);
                toSend.push_back(it->second.second);
            }

            // delete all of these triples from the global list
            globalTriples.erase(equalRange.first, equalRange.second);
          }
        }

        // send the list
        buf.write(toSend);
        LOG4ESPP_INFO(theLogger, "prepared fixed triple list before send particles");
  }


  // override parent function, this one should be empty
  void FixedTripleListAdress::beforeSendParticles(ParticleList& pl,
                                                    OutBuffer& buf) {
        //std::cout << storage->getSystem()->comm->rank() << ": beforeSendParticles() fixed pl (size " << pl.size() << ")\n";
  }


  // override parent function (use lookupAdrATParticle())
  void FixedTripleListAdress::onParticlesChanged() {

    LOG4ESPP_INFO(theLogger, "rebuild local bond list from global\n");

    this->clear();
    longint lastpid2 = -1;
    Particle *p1;
    Particle *p2;
    Particle *p3;

    for (GlobalTriples::const_iterator it = globalTriples.begin();
        it != globalTriples.end(); ++it) {

          //printf("lookup global triple %d %d %d\n", it->first, it->second.first, it->second.second);
          if (it->first != lastpid2) {
            p2 = storage->lookupAdrATParticle(it->first);
            if (p2 == NULL) {
              printf("SERIOUS ERROR: particle %d not available\n", it->first);
            }
            lastpid2 = it->first;
          }

          p1 = storage->lookupAdrATParticle(it->second.first);
          if (p1 == NULL) {
            printf("SERIOUS ERROR: 2nd particle %d not available\n", it->second.first);
          }

          p3 = storage->lookupAdrATParticle(it->second.second);
          if (p3 == NULL) {
            printf("SERIOUS ERROR: 3rd particle %d not available\n", it->second.second);
          }

          this->add(p1, p2, p3);
        }
    LOG4ESPP_INFO(theLogger, "regenerated local fixed pair list from global list");
  }

  /****************************************************
  ** REGISTRATION WITH PYTHON
  ****************************************************/

  void FixedTripleListAdress::registerPython() {

    using namespace espresso::python;

    bool (FixedTripleListAdress::*pyAdd)(longint pid1, longint pid2, longint pid3)
      = &FixedTripleListAdress::add;

    class_<FixedTripleListAdress, shared_ptr<FixedTripleListAdress> >
      ("FixedTripleListAdress",
              init <shared_ptr<storage::Storage>,
                     shared_ptr<FixedTupleList> >())
      .def("add", pyAdd)
      ;
  }

}
