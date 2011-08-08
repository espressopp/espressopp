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
        con4 = storage->onTuplesChanged.connect
          (boost::bind(&FixedTupleList::onParticlesChanged, this));

        //storage->setFixedTuples(this);

    }

    FixedTupleList::~FixedTupleList() {

        LOG4ESPP_INFO(theLogger, "~FixedTupleList");

        con1.disconnect();
        con2.disconnect();
        con3.disconnect();
        con4.disconnect();
    }


    bool FixedTupleList::addT(tuple pids) {


        // ADD THE LOCAL PARTICLES (pointers)
        Particle* vp, *at;
        std::vector<Particle*> tmp; // temporary vector
        std::vector<longint> pidstmp; // temporary vector
        longint pidK; // the pid used as key

        tuple::iterator it = pids.begin();
        vp = storage->lookupRealParticle(*it);
        if (!vp) { // Particle does not exist here, return false
            //std::cout << "particle " << *it << " not found in localParticles \n";
            return false;
        }
        pidK = *it; // first pid is key
        //std::cout << "Add key: " << *it << "\n";

        for (++it; it != pids.end(); ++it) {

            at = storage->lookupAdrATParticle(*it);
            if (!at) { // Particle does not exist here, return false
                std::cout << "ERROR: AT particle " << *it << " not found in localAdrATParticles \n";
                return false;
            }
            tmp.push_back(at);
            //std::cout << " add: " << *it << "\n";
            pidstmp.push_back(*it); // pidK is not in this vector
        }
        this->add(vp, tmp); // add to TupleList


        // ADD THE GLOBAL PARTICLES (ids)
        globalTuples.insert(make_pair(pidK, pidstmp));
        LOG4ESPP_INFO(theLogger, "added fixed tuple to global tuples");

        tmp.clear();
        pids.clear();
        pidstmp.clear();

        //std::cout << "\n";

        return true;
    }

    /* send global tuple information */
    void FixedTupleList::beforeSendParticles
                                    (ParticleList& pl, OutBuffer& buf) {

        //std::cout << "beforeSendParticles\n";

        //std::vector<longint> toSend;

        // loop over the particle list
        for (ParticleList::Iterator pit(pl); pit.isValid(); ++pit) {
            longint pidK = pit->id();
            LOG4ESPP_DEBUG(theLogger, "send particle with pid " << pidK << ", find tuples");

            // find particle that involves this particle id
            GlobalTuples::const_iterator it = globalTuples.find(pidK);
            if (it != globalTuples.end()) {

                // first write the pid of the first particle
                //toSend.push_back(pidK);
                //std::cout << "write pidK "<< pidK << "\n";
            	buf.write(pidK);

				// write the size of the vector
				int s = it->second.size();
				//toSend.push_back(s);
				//std::cout << "write s "<< s << "\n";
				buf.write(s);

				// iterate through vector and add pids
				for (tuple::const_iterator it2 = it->second.begin();
				 it2 != it->second.end(); ++it2) {
					//toSend.push_back(*it2);
					//std::cout << " write pid "<< *it2 << " (";

					Particle* tp = storage->lookupAdrATParticle(*it2);
					//std::cout << " write " << tp->getId() << " ("  << tp->getPos() << ")\n";
					buf.write(*tp);

					// remove AT particle from storage
					storage->removeAdrATParticle(*it2);
					//std::cout << "removing AT particle " << *it2 <<"\n";
				}

                // delete this pid from the global list
                globalTuples.erase(pidK);

            }
        }

        // send the list
        //buf.write(toSend);
    }

    /* recieve and rebuild global tuple information */
    void FixedTupleList::afterRecvParticles
                                    (ParticleList &pl, InBuffer& buf) {

        //std::cout << "afterRecvParticles\n";
        /*
        std::vector<longint> received, pids;
        int n;
        longint pidK;
        GlobalTuples::iterator it = globalTuples.begin();


        // receive the tuple list
        buf.read(received);
        int size = received.size();

        int i = 0;
        while (i < size) {
            // unpack the list
            pidK = received[i++];
            //std::cout << "receive pidK "<< pidK << "\n";

            n = received[i++];
            //std::cout << "receive n "<< n << "\n";

            for (; n > 0; --n) {
            	LOG4ESPP_DEBUG(theLogger, "received vector for pid " << pidK);
                //std::cout << "receive pid "<< received[i] << "\n";
                storage->addAdrATParticleFTPL(received[i]); // add AT particle to storage
                pids.push_back(received[i++]);
            }

            // add pids vector to global tuples
            it = globalTuples.insert(it, std::make_pair(pidK, pids));
            pids.clear();
        }

        if (i != size) {
            LOG4ESPP_ERROR(theLogger,
                    "recv particles might have read garbage\n");
        }

        LOG4ESPP_INFO(theLogger,
                "received fixed particle list after receive particles");
        */


        std::vector<longint> pids;
		int size, i, n;
		longint pidK;
		GlobalTuples::iterator it = globalTuples.begin();

		size = pl.size();

		for (i = 0; i < size; ++i) {
		    //std::cout << "i: " << i << "\n";

            // receive the tuple list
            //std::cout << "receive pidK: ";
            buf.read(pidK);
            //std::cout << pidK << "\n";

            //std::cout << "receive n: ";
            buf.read(n);
            //std::cout << n << "\n";

            for (; n > 0; --n) {
                LOG4ESPP_DEBUG(theLogger, "received vector for pid " << pidK);
                //std::cout << "receive pid "<< received[i] << "\n";
                //storage->addAdrATParticleFTPL(received[i]); // add AT particle to storage
                //pids.push_back(received[i++]);
                //Particle *p = storage->addAdrATParticleFTPL();
                Particle p;
                //std::cout << " read *p : ";
                buf.read(p);
                storage->addAdrATParticleFTPL(p);
                //std::cout << p.getId() << " (" << p.getPos() << ")\n";
                pids.push_back(p.id());
            }

            // add pids vector to global tuples
            it = globalTuples.insert(it, std::make_pair(pidK, pids));
            pids.clear();
		}

    }

    void FixedTupleList::onParticlesChanged() {

        //std::cout << "onParticlesChanged\n";

        LOG4ESPP_INFO(theLogger, "rebuild local particle list from global tuples\n");

        this->clear();
        //std::cout << " ---- CLEAR TUPLES ----  \n\n";

        Particle* vp, * at;
        std::vector<Particle*> tmp;

        GlobalTuples::const_iterator it = globalTuples.begin();

        // iterate through keys of map
        for (;it != globalTuples.end(); ++it) {
            vp = storage->lookupRealParticle(it->first);
            //std::cout << it->first << "\n";
            if (vp == NULL) {
            	printf("SERIOUS ERROR: VP particle %d not available\n", it->first);
            	exit(1);
            	return;
            }

            // iterate through vector in map
            for (tuple::const_iterator it2 = it->second.begin(); it2 != it->second.end(); ++it2) {
                at = storage->lookupAdrATParticle(*it2);
                //std::cout << "  " << *it2 << "\n";
                if (at == NULL) {
                	printf("SERIOUS ERROR: AT particle %d not available\n", *it2);
                	exit(1);
                	return;
                }
                tmp.push_back(at);
            }
            // add the particles
            this->add(vp, tmp);
            tmp.clear();
        }
        LOG4ESPP_INFO(theLogger, "regenerated local fixed list from global tuples");
        //std::cout << "\n";
    }

    /****************************************************
    ** REGISTRATION WITH PYTHON
    ****************************************************/

    void FixedTupleList::registerPython() {

      using namespace espresso::python;

      void (FixedTupleList::*pyAdd)(longint pid) = &FixedTupleList::add;

      class_<FixedTupleList, shared_ptr<FixedTupleList> >
        ("FixedTupleList", init<shared_ptr<storage::Storage> >())
        .def("add", pyAdd)
        .def("addTs", &FixedTupleList::addTs)
        ;
    }

}
