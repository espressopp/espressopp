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
//#include <iostream>
#include "ParticleGroup.hpp"
#include "storage/Storage.hpp"

namespace espressopp {

    ParticleGroup::ParticleGroup(shared_ptr <storage::Storage> _storage)
    : storage(_storage) {
        con_send = storage->beforeSendParticles.connect
                (boost::bind(&ParticleGroup::beforeSendParticles, this, _1, _2));
        con_recv = storage->afterRecvParticles.connect
                (boost::bind(&ParticleGroup::afterRecvParticles, this, _1, _2));
        con_changed = storage->onParticlesChanged.connect
                (boost::bind(&ParticleGroup::onParticlesChanged, this));
    }


    ParticleGroup::~ParticleGroup() {
        con_send.disconnect();
        con_recv.disconnect();
        con_changed.disconnect();
    }


    void ParticleGroup::add(longint pid) {
        // check unique
        //particles[pid] = pid;
    	particles.insert(pid);
        Particle *p1 = storage->lookupRealParticle(pid);
        if (p1) active[pid] = p1;
    }


    bool ParticleGroup::has(longint pid) {
    	return particles.find(pid) != particles.end();
    }

    // for debugging purpose
    void ParticleGroup::print() {
        std::cout << "####### I have " << active.size() << " active particles" << std::endl;
        for(iterator i=begin(); i!=end(); ++i ) {
            std::cout << i->getId() << " ";
        }
        std::cout << std::endl;
        //for (std::map<longint,longint>::iterator iter = particles.begin();
        for (std::set<longint>::iterator iter = particles.begin();
        									iter != particles.end(); ++iter) {
            //std::cout << iter->first << " ";
        	std::cout << *iter << " ";
        }
        std::cout << std::endl;
    }


    void ParticleGroup::beforeSendParticles(ParticleList& pl, class OutBuffer& buf) {
        // remove all particles that move to a different node
        for (ParticleList::Iterator pit(pl); pit.isValid(); ++pit) {
            longint pid = pit->id();

            std::map<longint, Particle*>::iterator p;
            p = active.find(pid);
            if (p != active.end())
                active.erase(p);
        }
    }


    void ParticleGroup::afterRecvParticles(ParticleList& pl, class InBuffer& buf) {
        // add all particles that moved to this node
        for (ParticleList::Iterator pit(pl);
                pit.isValid(); ++pit) {
            longint pid = pit->id();

            //std::map<longint, longint>::iterator p;
            std::set<longint>::iterator p;
            p = particles.find(pid);
            if (p != particles.end())
                active[pid] = NULL; // will be set in onchanged
        }
    }


    void ParticleGroup::onParticlesChanged() {
        std::map<longint, Particle*>::iterator p;
        std::list<longint> remove;

        // update all pointers
        for (p = active.begin(); p != active.end(); ++p) {
            if(!(p->second = storage->lookupRealParticle(p->first))) {
                fprintf(stderr, "ParticleGroup: non local particle\n");
                remove.push_back(p->first); // \todo do we need to check for ghosts here
            }
        }
        // remove ghosts, check how iterators behave on erase,
        // eventually this can be done in the loop
        while(!remove.empty()) {
            active.erase(remove.front());
            remove.pop_front();
        }
    }


    void ParticleGroup::registerPython() {

        using namespace espressopp::python;

        class_< ParticleGroup, shared_ptr <ParticleGroup> >
                ("ParticleGroup", init <shared_ptr <storage::Storage> >())
                .def("add", &ParticleGroup::add)
                .def("show", &ParticleGroup::print)
                .def("has", &ParticleGroup::has)
                .def("size", &ParticleGroup::size)
                ;
    }
}
