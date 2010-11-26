#include "python.hpp"
#include <iostream>
#include "ParticleGroup.hpp"

namespace espresso {

ParticleGroup::ParticleGroup(shared_ptr< storage::Storage > _storage)
{
/*    con_send = storage->beforeSendParticles.connect
      (boost::bind(&FixedPairList::beforeSendParticles, this, _1, _2));
    con_recv = storage->afterRecvParticles.connect
      (boost::bind(&FixedPairList::afterRecvParticles, this, _1, _2));
    con_changed = storage->onParticlesChanged.connect
      (boost::bind(&FixedPairList::onParticlesChanged, this));*/
}

ParticleGroup::~ParticleGroup()
{
    /*con_send.disconnect();
    con_recv.disconnect();
    con_changed.disconnect();*/
}

  void ParticleGroup::add(longint pid)
  {
      // check unique
      particles.push_back(pid);
  }

  // for debugging purpose
  void ParticleGroup::print()
  {
    for(std::list<longint>::iterator iter = particles.begin();
            iter!=particles.end(); ++iter) {
        std::cout << *iter << " ";
    }
    std::cout << std::endl;
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
