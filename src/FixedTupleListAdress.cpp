/*
  Copyright (C) 2012,2013
      Max Planck Institute for Polymer Research
  Copyright (C) 2008,2009,2010,2011
      Max-Planck-Institute for Polymer Research & Fraunhofer SCAI
  
  This file is part of ESPResSo++.
  
  ESPResSo++ is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.
  
  ESPResSo++ is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.
  
  You should have received a copy of the GNU General Public License
  along with this program.  If not, see <http://www.gnu.org/licenses/>. 
*/

#include "python.hpp"

#include "FixedTupleListAdress.hpp"

//#include <vector>
//#include <utility>
//#include <algorithm>
#include <boost/bind.hpp>
#include "storage/Storage.hpp"
#include "Buffer.hpp"

#include "System.hpp"
#include "bc/BC.hpp"

#include "esutil/Error.hpp"
//using namespace std;

namespace espressopp {

    LOG4ESPP_LOGGER(FixedTupleListAdress::theLogger, "FixedTupleListAdress");

    FixedTupleListAdress::FixedTupleListAdress(shared_ptr<storage::Storage> _storage)
        : storage(_storage), globalTuples(){

        LOG4ESPP_INFO(theLogger, "construct FixedTupleListAdress");

        sigBeforeSend = storage->beforeSendParticles.connect
           (boost::bind(&FixedTupleListAdress::beforeSendParticles, this, _1, _2));
        sigAfterRecv = storage->afterRecvParticles.connect
           (boost::bind(&FixedTupleListAdress::afterRecvParticles, this, _1, _2));
        sigOnTupleChanged = storage->onTuplesChanged.connect
           (boost::bind(&FixedTupleListAdress::onParticlesChanged, this));
    }

    FixedTupleListAdress::~FixedTupleListAdress() {

        LOG4ESPP_INFO(theLogger, "~FixedTupleListAdress");
        sigBeforeSend.disconnect();
        sigAfterRecv.disconnect();
        sigOnTupleChanged.disconnect();
    }

    bool FixedTupleListAdress::addT(tuple pids) {
        bool returnVal = true;
        System& system = storage->getSystemRef();
        esutil::Error err(system.comm);

        // ADD THE LOCAL PARTICLES (pointers)
        Particle* vp, *at;
        std::vector<Particle*> tmp; // temporary vector
        std::vector<longint> pidstmp; // temporary vector
        longint pidK; // the pid used as key

        tuple::iterator it = pids.begin();
        vp = storage->lookupRealParticle(*it);
        if (!vp) { // Particle does not exist here, return false
            //std::cout << "particle " << *it << " not found in localParticles \n";
            returnVal = false;
        }
        else{
          pidK = *it; // first pid is key
          //std::cout << "Add key: " << *it << "\n";

          for (++it; it != pids.end(); ++it) {

              at = storage->lookupAdrATParticle(*it);
              if (!at) { // Particle does not exist here, return false
                  std::stringstream msg;
                  msg << "ERROR: AT particle " << *it << " not found in localAdrATParticles \n";
                  err.setException( msg.str() );
                  returnVal = false;
                  break;
              }
              tmp.push_back(at);
              //std::cout << " add: " << *it << "\n";
              pidstmp.push_back(*it); // pidK is not in this vector
          }
        }
        err.checkException();
        
        if(returnVal){
            this->add(vp, tmp); // add to TupleList

            // ADD THE GLOBAL PARTICLES (ids)
            globalTuples.insert(make_pair(pidK, pidstmp));
        }
        LOG4ESPP_INFO(theLogger, "added fixed tuple to global tuples");

        tmp.clear();
        pids.clear();
        pidstmp.clear();

        //std::cout << "\n";

        return returnVal;
    }

    /* send global tuple information */
    void FixedTupleListAdress::beforeSendParticles
                                    (ParticleList& pl, OutBuffer& buf) {

        std::vector<longint> atpl;

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
				atpl.reserve(s);

				// iterate through vector and add pids
				//std::cout << storage->getRank() << ": removing AT particles ";
				for (tuple::const_iterator it2 = it->second.begin();
				 it2 != it->second.end(); ++it2) {
					//toSend.push_back(*it2);
					//std::cout << " write pid "<< *it2 << " (";

					Particle* tp = storage->lookupAdrATParticle(*it2);
					//std::cout << " write " << tp->getId() << " ("  << tp->getPos() << ")\n";
					buf.write(*tp);
					atpl.push_back(*it2);

					// remove AT particle from storage
					storage->removeAdrATParticle(*it2);
					//std::cout << " " << *it2;
				}
				//std::cout << "\n";

                // delete this pid from the global list
                globalTuples.erase(pidK);

            }
        }

        beforeSendATParticles(atpl, buf);
    }

    /* recieve and rebuild global tuple information */
    void FixedTupleListAdress::afterRecvParticles
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
            //std::cout << pidK << " at ";


            /*
            // testing
            Particle* vp = storage->lookupRealParticle(pidK);
            Real3D vpp = vp->position();

            // see where VP is folded
            Real3D vpp_old = vpp;
            Real3D vpp_new = vpp;
            Real3D moved;
            int dir;
            Int3D image(0,0,0);

            for (dir = 0; dir < 3; ++dir) {
                storage->getSystem()->bc->foldCoordinate(vpp_new, image, dir);
                if (vpp_new[dir] != vpp_old[dir]) {
                    moved[dir] = vpp_old[dir] - vpp_new[dir];
                    break; // do not continue looping
                }
            }
            */


            //std::cout << "receive n: ";
            buf.read(n);
            //std::cout << n << "\n";

            //std::cout << storage->getRank() << ": add AT particles ";
            for (; n > 0; --n) {
                LOG4ESPP_DEBUG(theLogger, "received vector for pid " << pidK);
                /*storage->addAdrATParticleFTPL(received[i]); // add AT particle to storage
                pids.push_back(received[i++]);
                Particle *p = storage->addAdrATParticleFTPL();*/

                Particle p;
                //std::cout << " read *p : ";
                buf.read(p);

                //std::cout << " AT particle " << p.id() << " at " << p.position() << " ";
                //p.position()[dir] = p.position()[dir] - moved[dir];
                //std::cout << "--> moved to " << p.position() << "\n";

                storage->addAdrATParticleFTPL(p);

                //std::cout << p.getId() << " at " << p.position() << "\n";
                pids.push_back(p.id());
            }
            //std::cout << "\n";

            // add pids vector to global tuples
            it = globalTuples.insert(it, std::make_pair(pidK, pids));
            pids.clear();
		}
        afterRecvATParticles(pl, buf);
    }

    void FixedTupleListAdress::onParticlesChanged() {
      
      // TODO errors should be thrown in a nice way

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
            //std::cout << storage->getRank() << ": loopup for AT particle: ";
            for (tuple::const_iterator it2 = it->second.begin(); it2 != it->second.end(); ++it2) {
                at = storage->lookupAdrATParticle(*it2);
                if (at == NULL) {
                	printf("SERIOUS ERROR: AT particle %d not available\n", *it2);
                	exit(1);
                	return;
                }


                // fold AT coordinates to follow VP if necessary
                real dif;
                Real3D boxL = storage->getSystem()->bc->getBoxL();
                for (int dir = 0; dir < 3; ++dir) {
                    dif = vp->position()[dir] - at->position()[dir];
                    if (dif > boxL[dir]/2) {
                        //std::cout << "VP " << vp->getId() << " at " << vp->position() << "\n";
                        //std::cout << " moving AT " << at->getId() << " " << at->position();
                        at->position()[dir] = at->position()[dir] + boxL[dir];
                        //std::cout << " --> " << at->position() << "\n";
                        at->image()[dir] = vp->image()[dir];
                    }
                    else if (dif < -boxL[dir]/2) {
                        //std::cout << "VP " << vp->getId() << " at " << vp->position() << "\n";
                        //std::cout << " moving AT " << at->getId() << " " << at->position();
                        at->position()[dir] = at->position()[dir] - boxL[dir];
                        //std::cout << " --> " << at->position() << "\n";
                        at->image()[dir] = vp->image()[dir];
                    }
                }


                //std::cout << " " << *it2;
                tmp.push_back(at);
            }
            //std::cout << "\n";
            // add the particles to tuples
            this->add(vp, tmp);
            tmp.clear();
        }
        LOG4ESPP_INFO(theLogger, "regenerated local fixed list from global tuples");
        //std::cout << "\n";
    }

    /****************************************************
    ** REGISTRATION WITH PYTHON
    ****************************************************/

    void FixedTupleListAdress::registerPython() {

      using namespace espressopp::python;

      void (FixedTupleListAdress::*pyAdd)(longint pid) = &FixedTupleListAdress::add;

      class_<FixedTupleListAdress, shared_ptr<FixedTupleListAdress>, boost::noncopyable >
        ("FixedTupleListAdress", init<shared_ptr<storage::Storage> >())
        .def("add", pyAdd)
        .def("addTs", &FixedTupleListAdress::addTs)
        ;
    }

}
