#include "python.hpp"
#include <iostream>
#include "ParticleGroup.hpp"
#include "storage/Storage.hpp"

namespace espresso {

ParticleGroup::ParticleGroup(shared_ptr< storage::Storage > _storage)
  : storage(_storage)
{
    con_send = storage->beforeSendParticles.connect
      (boost::bind(&ParticleGroup::beforeSendParticles, this, _1, _2));
    con_recv = storage->afterRecvParticles.connect
      (boost::bind(&ParticleGroup::afterRecvParticles, this, _1, _2));
    con_changed = storage->onParticlesChanged.connect
      (boost::bind(&ParticleGroup::onParticlesChanged, this));
}

ParticleGroup::~ParticleGroup()
{
    con_send.disconnect();
    con_recv.disconnect();
    con_changed.disconnect();
}

  void ParticleGroup::add(longint pid)
  {
      // check unique
      particles.push_back(pid);
      Particle *p1 = storage->lookupRealParticle(pid);
      if(p1)
          active[pid] = p1;
  }

  // for debugging purpose
  void ParticleGroup::print()
  {
    std::cout << "I have " << active.size() << " active particles" << std::endl;
    for(std::list<longint>::iterator iter = particles.begin();
            iter!=particles.end(); ++iter) {
        std::cout << *iter << " ";
    }
    std::cout << std::endl;
  }
    void ParticleGroup::beforeSendParticles(ParticleList& pl,
			     class OutBuffer& buf)
    {
        std::cout << "before send\n";
    }
    void ParticleGroup::afterRecvParticles(ParticleList& pl,
			    class InBuffer& buf)
    {
        std::cout << "after recv\n";
    }
    void ParticleGroup::onParticlesChanged()
    {
        std::cout << "changed\n";
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
