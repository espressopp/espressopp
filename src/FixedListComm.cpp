#include "FixedListComm.hpp"


#include <boost/bind.hpp>
#include "storage/Storage.hpp"
#include "Buffer.hpp"

//using namespace std;

namespace espresso {

    LOG4ESPP_LOGGER(FixedListComm::theLogger, "FixedListComm");

    FixedListComm::FixedListComm(shared_ptr <storage::Storage> _storage)
        : storage(_storage), globalLists(){

        //std::cout << "fixedlist" << std::endl;

        LOG4ESPP_INFO(theLogger, "construct FixedPairList");
        con1 = storage->beforeSendParticles.connect
          (boost::bind(&FixedListComm::beforeSendParticles, this, _1, _2));
        con2 = storage->afterRecvParticles.connect
          (boost::bind(&FixedListComm::afterRecvParticles, this, _1, _2));
        con3 = storage->onParticlesChanged.connect
          (boost::bind(&FixedListComm::onParticlesChanged, this));
    }

    FixedListComm::~FixedListComm() {

        //std::cout << "~fixedlist" << std::endl;

        LOG4ESPP_INFO(theLogger, "~FixedListComm");

        con1.disconnect();
        con2.disconnect();
        con3.disconnect();
    }

    bool FixedListComm::add(pvec pids) {

        //if (pid1 > pid2) std::swap(pid1, pid2);


        // ADD THE LOCAL PARTICLES
        Particle* p;
        std::vector<Particle*> tmp;
        for (pvec::iterator it = pids.begin(); it!=pids.end(); ++it) {
            p = storage->lookupRealParticle(*it);
            //Particle* p2 = storage->lookupLocalParticle(pid2);
            if (!p)
                // Particle does not exist here, return false
                return false;
            tmp.push_back(p);
        }

        //add(tmp);
        if (pids.size() == 2) this->PairList::add(tmp.at(0), tmp.at(1));
        else if (pids.size() == 3) this->TripleList::add(tmp.at(0), tmp.at(1), tmp.at(2));
        else if (pids.size() == 4) this->QuadrupleList::add(tmp.at(0), tmp.at(1), tmp.at(2), tmp.at(3));
        tmp.clear();

        // ADD THE GLOBAL PARTICLES
        std::pair<GlobalList::const_iterator, GlobalList::const_iterator>
        equalRange = globalLists.equal_range(pids[0]);
        // see whether the particle already has pairs
        if (equalRange.first == globalLists.end()) {
            // it hasn't, insert the new pair
            //globalLists.insert(make_pair(pid1, pid2));
            globalLists.insert(make_pair(pids.at(0), pids));
        }
        else {// otherwise test whether the pair already exists
            for (GlobalList::const_iterator it = equalRange.first;
            it != equalRange.second; ++it) {
                if (it->second == pids) {
                   // TODO: Pair already exists, generate error!
                }
                // if not, insert the new pair
                globalLists.insert(equalRange.first, make_pair(pids.at(0), pids));
            }
        }
        LOG4ESPP_INFO(theLogger, "added fixed pair to global pair list");
        return true;
    }

    void FixedListComm::beforeSendParticles
        (ParticleList& pl, OutBuffer& buf) {
        std::vector<longint> toSend;


        // loop over the particle list
        for (ParticleList::Iterator pit(pl); pit.isValid(); ++pit) {
            longint pid = pit->id();
            LOG4ESPP_DEBUG(theLogger, "send particle with pid " << pid << ", find pairs");

            // find all particles that involve this particle id
            int n = globalLists.count(pid);
            if (n > 0) {
                std::pair<GlobalList::const_iterator,
                 GlobalList::const_iterator>
                equalRange = globalLists.equal_range(pid);
                // get the length of the vector in the map
                int l = equalRange.first->second.size();
                toSend.reserve(toSend.size()+3+n*l);
                // first write the pid of the first particle
                toSend.push_back(pid);
                // then the number of partners
                toSend.push_back(n);
                // then the size of the vector
                toSend.push_back(l);
                // and then the pids of the partners
                for (GlobalList::const_iterator it = equalRange.first;
                it != equalRange.second; ++it) {
                    // iterate through vector to add pids
                    for (pvec::const_iterator it2=it->second.begin();
                    it2!=it->second.end(); ++it2) {
                        toSend.push_back(*it2);
                        LOG4ESPP_DEBUG(theLogger, "send global bond: pid "
                               << pid << " and its vector");
                    }
                }

                // delete all of these pairs from the global list
                globalLists.erase(equalRange.first, equalRange.second);
            }
        }
        // send the list
        buf.write(toSend);
        LOG4ESPP_INFO(theLogger, "prepared fixed pair list before send particles");
    }

    void FixedListComm::afterRecvParticles
        (ParticleList &pl,InBuffer& buf) {
        std::vector<longint> received, pids;
        int n, l, tl;
        longint pid1;
        GlobalList::iterator it = globalLists.begin();


        // receive the bond list
        buf.read(received);
        int size = received.size();
        int i = 0;
        while (i < size) {
            // unpack the list
            pid1 = received[i++];
            n = received[i++];
            l = received[i++];
            LOG4ESPP_DEBUG(theLogger,
                    "recv particle " << pid1 << ", has " << n << " global pairs");
            for (; n > 0; --n) {
                for (tl = l; tl > 0; --tl) {
                    pids.push_back(received[i++]);
                }
                // add the bond to the global list
                LOG4ESPP_DEBUG(theLogger, "received vector for pid " << pid1);
                it = globalLists.insert(it, std::make_pair(pid1, pids));
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

    void FixedListComm::onParticlesChanged() {
        LOG4ESPP_INFO(theLogger, "rebuild local bond list from global\n");
        //this->clear();
        longint lastpid1 = -1;

        Particle* p;
        std::vector<Particle*> tmp;

        GlobalList::const_iterator it = globalLists.begin();
        int size = it->second.size();
        if (size == 2) this->PairList::clear();
        else if (size == 3) this->TripleList::clear();
        else if (size == 4) this->QuadrupleList::clear();


        // iterate through keys of map
        for (;it != globalLists.end(); ++it) {
            if (it->first != lastpid1) { // don't check same pid twice
                p = storage->lookupRealParticle(it->first);
                if (p == NULL) {
                    printf("SERIOUS ERROR: particle %d not available\n", it->first);
                }
                lastpid1 = it->first;
            }
            // iterate through vector in map
            for (pvec::const_iterator it2 = it->second.begin(); it2!=it->second.end(); ++it2) {
                p = storage->lookupLocalParticle(*it2);
                if (p == NULL) {
                  printf("SERIOUS ERROR: particle %d not available\n", *it2);
                }
                tmp.push_back(p);
            }
            // add the particles
            if (size == 2) this->PairList::add(tmp.at(0), tmp.at(1));
            else if (size == 3) this->TripleList::add(tmp.at(0), tmp.at(1), tmp.at(2));
            else if (size == 4) this->QuadrupleList::add(tmp.at(0), tmp.at(1), tmp.at(2), tmp.at(3));
            tmp.clear();
        }
        LOG4ESPP_INFO(theLogger, "regenerated local fixed list from global list");
    }


} // ns espresso
