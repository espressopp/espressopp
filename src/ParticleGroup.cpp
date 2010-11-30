#include "python.hpp"
#include <iostream>
#include "ParticleGroup.hpp"
#include "storage/Storage.hpp"

namespace espresso {

    ParticleGroup::ParticleGroup(shared_ptr< storage::Storage > _storage)
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
        particles[pid]=pid;
        Particle *p1 = storage->lookupRealParticle(pid);
        if (p1)
            active[pid] = p1;
    }

    // for debugging purpose

    void ParticleGroup::print() {
        std::cout << "####### I have " << active.size() << " active particles" << std::endl;
        for(iterator i=begin(); i!=end(); ++i )
            std::cout << (*i)->getId() << " ";
        std::cout << std::endl;
        for (std::map<longint,longint>::iterator iter = particles.begin();
                iter != particles.end(); ++iter) {
            std::cout << iter->first << " ";
        }
        std::cout << std::endl;
    }

    void ParticleGroup::beforeSendParticles(ParticleList& pl,
            class OutBuffer& buf) {
        // loop over the particle list
        for (ParticleList::Iterator pit(pl);
                pit.isValid(); ++pit) {
            longint pid = pit->id();

            std::map<longint, Particle*>::iterator p;
            p = active.find(pid);
            if (p != active.end())
                active.erase(p);
        }

    }

    void ParticleGroup::afterRecvParticles(ParticleList& pl,
            class InBuffer& buf) {
        for (ParticleList::Iterator pit(pl);
                pit.isValid(); ++pit) {
            longint pid = pit->id();

            std::map<longint, longint>::iterator p;
            p = particles.find(pid);
            if (p != particles.end())
                active[pid] = NULL; // will by set in onchanged
        }
    }

    void ParticleGroup::onParticlesChanged() {
        std::map<longint, Particle*>::iterator p;
        std::list<longint> remove;
        for (p = active.begin(); p != active.end(); ++p) {
            // TODO this should be ok but check
            if(!(p->second = storage->lookupRealParticle(p->first)))
                    remove.push_back(p->first);
        }
        // needed to remove ghosts
        while(!remove.empty()) {
            active.erase(remove.front());
            remove.pop_front();
        }


    }

    void ParticleGroup::registerPython() {

        using namespace espresso::python;

        class_< ParticleGroup, shared_ptr< ParticleGroup > >
                ("ParticleGroup", init< shared_ptr< storage::Storage > >())
                .def("add", &ParticleGroup::add)
                .def("show", &ParticleGroup::print)
                ;
    }
}
