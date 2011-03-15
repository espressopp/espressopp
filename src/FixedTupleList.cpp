#include "python.hpp"

#include "FixedTupleList.hpp"

//#include <vector>
//#include <utility>
//#include <algorithm>
#include <boost/bind.hpp>
#include "storage/Storage.hpp"
#include "Buffer.hpp"

//using namespace std;

namespace espresso {

    LOG4ESPP_LOGGER(FixedTupleList::theLogger, "FixedTupleList");

    FixedTupleList::FixedTupleList(shared_ptr<storage::Storage> _storage)
        : storage(_storage), globalTuples(){

        LOG4ESPP_INFO(theLogger, "construct FixedTupleList");

        con1 = storage->beforeSendParticles.connect
          (boost::bind(&FixedTupleList::beforeSendParticles, this, _1, _2));
        con2 = storage->afterRecvParticles.connect
          (boost::bind(&FixedTupleList::afterRecvParticles, this, _1, _2));
        con3 = storage->onParticlesChanged.connect
          (boost::bind(&FixedTupleList::onParticlesChanged, this));
    }

    FixedTupleList::~FixedTupleList() {

        LOG4ESPP_INFO(theLogger, "~FixedTupleList");

        con1.disconnect();
        con2.disconnect();
        con3.disconnect();
    }


    bool FixedTupleList::add(tuple pids) {

        // ADD THE LOCAL PARTICLES
        Particle* p;
        std::vector<Particle*> tmp; // temporary vector
        std::vector<longint> pidstmp; // temporary vector
        longint pidK; // the pid used as key

        tuple::iterator it = pids.begin();
        p = storage->lookupRealParticle(*it);
        if (!p) // Particle does not exist here, return false
            return false;
        pidK = *it; // first pid is key
        tmp.push_back(p);

        for (++it; it != pids.end(); ++it) {
            p = storage->lookupLocalParticle(*it);
            if (!p) // Particle does not exist here, return false
                return false;
            tmp.push_back(p);
            pidstmp.push_back(*it); // pidK is not in this vector
        }
        this->add(tmp);
        tmp.clear();
        pids.clear();

        // ADD THE GLOBAL PARTICLES
        std::pair<GlobalTuples::const_iterator, GlobalTuples::const_iterator>
        equalRange = globalTuples.equal_range(pidK);
        // see whether the particle already has partners
        if (equalRange.first == globalTuples.end()) {
            // it hasn't, insert the new pair
            globalTuples.insert(make_pair(pidK, pidstmp));
        }
        else {// otherwise test whether the pair already exists
            for (GlobalTuples::const_iterator it = equalRange.first;
            it != equalRange.second; ++it) {
                if (it->second == pidstmp) {
                   // TODO: Pair already exists, generate error!
                }
                // if not, insert the new pair
                globalTuples.insert(equalRange.first, make_pair(pidK, pidstmp));
            }
        }
        LOG4ESPP_INFO(theLogger, "added fixed tuple to global tuples");

        return true;
    }

    void FixedTupleList::beforeSendParticles
                                    (ParticleList& pl, OutBuffer& buf) {

        std::vector<longint> toSend;

        // loop over the particle list
        for (ParticleList::Iterator pit(pl); pit.isValid(); ++pit) {
            longint pidK = pit->id();
            LOG4ESPP_DEBUG(theLogger, "send particle with pid " << pidK << ", find pairs");

            // find all particles that involve this particle id
            int n = globalTuples.count(pidK);
            if (n > 0) {
                std::pair<GlobalTuples::const_iterator,
                 GlobalTuples::const_iterator>
                equalRange = globalTuples.equal_range(pidK);

                // first write the pid of the first particle
                toSend.push_back(pidK);
                // then the number of partners
                toSend.push_back(n);
                // and then the size and pids of the partners
                for (GlobalTuples::const_iterator it = equalRange.first;
                it != equalRange.second; ++it) {
                    // write the size of the vector first
                    int s = it->second.size();
                    toSend.push_back(s);

                    // iterate through vector and add pids
                    for (tuple::const_iterator it2 = it->second.begin();
                    it2 != it->second.end(); ++it2) {
                        toSend.push_back(*it2);
                        LOG4ESPP_DEBUG(theLogger, "send global bond: pid "
                               << pidK << " and its vector");
                    }
                }

                // delete all of these pairs from the global list
                globalTuples.erase(equalRange.first, equalRange.second);
            }
        }
        // send the list
        buf.write(toSend);
        LOG4ESPP_INFO(theLogger, "prepared fixed tuple before send particles");
    }

    void FixedTupleList::afterRecvParticles
                                    (ParticleList &pl,InBuffer& buf) {

        std::vector<longint> received, pids;
        int n, s, ts;
        longint pidK;
        GlobalTuples::iterator it = globalTuples.begin();


        // receive the bond list
        buf.read(received);
        int size = received.size();
        int i = 0;
        while (i < size) {
            // unpack the list
            pidK = received[i++];
            n = received[i++];
            LOG4ESPP_DEBUG(theLogger,
                    "recv particle " << pidK << ", has " << n << " global pairs");
            for (; n > 0; --n) {
                s = received[i++];
                for (ts = s; ts > 0; --ts) { // read the pids vector
                    pids.push_back(received[i++]);
                }
                // add pids vector to global tuples
                LOG4ESPP_DEBUG(theLogger, "received vector for pid " << pidK);
                it = globalTuples.insert(it, std::make_pair(pidK, pids));
                pids.clear();
            }
        }
        if (i != size) {
            LOG4ESPP_ERROR(theLogger,
                    "recv particles might have read garbage\n");
        }
        LOG4ESPP_INFO(theLogger,
                "received fixed particle list after receive particles");
    }

    void FixedTupleList::onParticlesChanged() {
        LOG4ESPP_INFO(theLogger, "rebuild local particle list from global tuples\n");
        this->clear();
        longint lastpidK = -1;

        Particle* p;
        std::vector<Particle*> tmp;

        GlobalTuples::const_iterator it = globalTuples.begin();

        // iterate through keys of map
        for (;it != globalTuples.end(); ++it) {
            if (it->first != lastpidK) { // don't check same pidK twice
                lastpidK = it->first;
                p = storage->lookupRealParticle(lastpidK);
                if (p == NULL)
                    printf("SERIOUS ERROR: particle %d not available\n", lastpidK);
                tmp.push_back(p);
            }
            // iterate through vector in map
            for (tuple::const_iterator it2 = it->second.begin(); it2 != it->second.end(); ++it2) {
                p = storage->lookupLocalParticle(*it2);
                if (p == NULL)
                  printf("SERIOUS ERROR: particle %d not available\n", *it2);
                tmp.push_back(p);
            }
            // add the particles
            this->add(tmp);
            tmp.clear();
        }
        LOG4ESPP_INFO(theLogger, "regenerated local fixed list from global tuples");
    }
}
