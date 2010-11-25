#include "ParticleGroup.h"

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

  void ParticleGroup::registerPython() {

    using namespace espresso::python;

    bool (ParticleGroup::*pyAdd)(longint pid)
      = &ParticleGroup::add;


    class_< ParticleGroup, shared_ptr< ParticleGroup > >
      ("ParticleGroup", init< shared_ptr< storage::Storage > >())
      .def("add", pyAdd)
      ;

    bool (ParticleGroup::*pyPrint)()
      = &ParticleGroup::print;


    class_< ParticleGroup, shared_ptr< ParticleGroup > >
      ("ParticleGroup", init< shared_ptr< storage::Storage > >())
      .def("add", pyPrint)
      ;
  }
}
