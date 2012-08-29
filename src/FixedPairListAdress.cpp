#include <sstream>
#include "python.hpp"

#include "FixedPairListAdress.hpp"

//#include <vector>
//#include <utility>
//#include <algorithm>
#include <boost/bind.hpp>
#include "storage/Storage.hpp"
#include "Buffer.hpp"

//using namespace std;

namespace espresso {

  LOG4ESPP_LOGGER(FixedPairListAdress::theLogger, "FixedPairListAdress");


  FixedPairListAdress::
  FixedPairListAdress(shared_ptr< storage::Storage > _storage,
          shared_ptr<FixedTupleList> _fixedtupleList)
          : FixedPairList(_storage), fixedtupleList(_fixedtupleList) {
    LOG4ESPP_INFO(theLogger, "construct FixedPairListAdress");

    con = fixedtupleList->beforeSendATParticles.connect
          (boost::bind(&FixedPairListAdress::beforeSendATParticles, this, _1, _2));
  }

  FixedPairListAdress::~FixedPairListAdress() {

    LOG4ESPP_INFO(theLogger, "~FixedPairListAdress");

    con.disconnect();
  }


  // override parent function (use lookupAdrATParticle)
  bool FixedPairListAdress::add(longint pid1, longint pid2) {

    if (pid1 > pid2)
      std::swap(pid1, pid2);

    // ADD THE LOCAL PAIR
    Particle *p1 = storage->lookupAdrATParticle(pid1);
    Particle *p2 = storage->lookupAdrATParticle(pid2);
    if (!p1)
      // Particle does not exist here, return false
      return false;
    if (!p2) {
      std::stringstream err;
      err << "atomistic bond particle p2 " << pid2 << " does not exists here and cannot be added";
      std::runtime_error(err.str());
    }
    // add the pair locally
    this->add(p1, p2);

    globalPairs.insert(std::make_pair(pid1, pid2));
    LOG4ESPP_INFO(theLogger, "added fixed pair to global pair list");
    return true;
  }


  void FixedPairListAdress::beforeSendATParticles(std::vector<longint>& atpl,
          OutBuffer& buf) {

        //std::cout << "beforeSendATParticles() fixed pl (size " << atpl.size() << ")\n";

        std::vector< longint > toSend;

        // loop over the VP particle list
        for (std::vector<longint>::iterator it = atpl.begin();
                it != atpl.end(); ++it) {
          longint pid = *it;

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
            globalPairs.erase(pid);
          }
        }

        // send the list
        buf.write(toSend);
        LOG4ESPP_INFO(theLogger, "prepared fixed pair list before send particles");
  }


  // override parent function, this one should be empty
  void FixedPairListAdress::beforeSendParticles(ParticleList& pl,
                                                    OutBuffer& buf) {
        //std::cout << storage->getSystem()->comm->rank() << ": beforeSendParticles() fixed pl (size " << pl.size() << ")\n";
  }


  // override parent function (use lookupAdrATParticle())
  void FixedPairListAdress::onParticlesChanged() {

    LOG4ESPP_INFO(theLogger, "rebuild local bond list from global\n");

    this->clear();
    longint lastpid1 = -1;
    Particle *p1;
    Particle *p2;

    for (GlobalPairs::const_iterator it = globalPairs.begin();
	 it != globalPairs.end(); ++it) {

        if (it->first != lastpid1) {
            p1 = storage->lookupAdrATParticle(it->first);
            if (p1 == NULL) {
                std::stringstream err;
                err << "atomistic bond particle p1 " << it->first << " does not exists here";
                std::runtime_error(err.str());
            }
            lastpid1 = it->first;
        }

        p2 = storage->lookupAdrATParticle(it->second);
        if (p2 == NULL) {
            std::stringstream err;
            err << "atomistic bond particle p2 " << it->second << " does not exists here";
            std::runtime_error(err.str());
        }

        //std::cout << " adding (" << p1->getId() << ", " << p2->getId() << ")\n";
        this->add(p1, p2);
    }
    LOG4ESPP_INFO(theLogger, "regenerated local fixed pair list from global list");
  }

  /****************************************************
  ** REGISTRATION WITH PYTHON
  ****************************************************/

  void FixedPairListAdress::registerPython() {

    using namespace espresso::python;

    bool (FixedPairListAdress::*pyAdd)(longint pid1, longint pid2)
      = &FixedPairListAdress::add;

    class_<FixedPairListAdress, shared_ptr<FixedPairListAdress> >
      ("FixedPairListAdress",
              init <shared_ptr<storage::Storage>,
                     shared_ptr<FixedTupleList> >())
      .def("add", pyAdd)
      ;
  }

}
