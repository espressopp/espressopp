#include "Storage.hpp"
#include "python.hpp"

#include <boost/foreach.hpp>
#include <algorithm>
#include <stdexcept>
#include "Property.hpp"

using namespace espresso;
using namespace espresso::storage;
using namespace espresso::python;

Storage::Storage(bc::BC::SelfPtr _bc): bc(_bc) {}

Storage::~Storage() {}

void Storage::positionPropertyModified() {
}

void Storage::_addParticle(ParticleId id) {
  addParticle(id);
  positionPropertyModified();
}

void Storage::deleteProperty(PropertyId id) {
  getTupleVector().eraseProperty(id);
}

Storage::SelfPtr
Storage::getStorage() { return shared_from_this(); }
Storage::SelfPtr 
Storage::getLeftStorage() { return getStorage(); }
Storage::SelfPtr 
Storage::getRightStorage() { return getStorage(); }

// Implementation of the particles::Set interface

bool Storage::contains(ParticleHandle) { return true; }
bool Storage::contains(ParticleId) { return true; }

// Implementation of the pairs::Set interface

namespace {
  // basic pair traversing scheme
  struct InnerPairTraverser : particles::Computer {
    bool cont;
    pairs::Computer &computer;
    bc::BC &bc;

    storage::ParticleHandle p1;
    storage::PropertyHandle< Real3D > pos;

    InnerPairTraverser(pairs::Computer &_computer, 
		       Storage &_store) 
      : cont(true), 
	computer(_computer), bc(*_store.getBoundaryConditions())
    {}

    void setP1(ParticleHandle _p1) { p1 = _p1; }

    void prepare(Storage::SelfPtr storage) {
      pos = storage->getPositionPropertyHandle();
    }

    bool apply(ParticleHandle p2) {
      if (p1 == p2) return false;
      Real3D pos1 = pos[p1];
      Real3D pos2 = pos[p2];
      Real3D dist = bc.getDist(pos1, pos2);

      cont = computer.apply(dist, p1, p2);
      return cont;
    }
  };

  struct OuterPairTraverser : particles::Computer {
    InnerPairTraverser innerTraverser;
    Storage &store;

    OuterPairTraverser(InnerPairTraverser &_innerTraverser, Storage &_store) 
      : innerTraverser(_innerTraverser), store(_store)
    {}

    void prepare(Storage::SelfPtr set) {}

    bool apply(ParticleHandle p1) {
      innerTraverser.setP1(p1);
      store.foreach(innerTraverser);
      return innerTraverser.cont;
    }
  };
}

bool Storage::foreachPairApply(pairs::Computer &computer) {
  computer.prepare(getStorage(), getStorage());
  InnerPairTraverser innerTraverser(computer, *this);
  OuterPairTraverser outerTraverser(innerTraverser, *this);
  bool res = this->foreach(outerTraverser);
  computer.finalize();
  return res;
}

void Storage::setIdProperty(boost::shared_ptr< Property< ParticleId > >prop) {
  if (particleIdProperty != PropertyId()) {
    throw std::runtime_error("Storage::setIdProperty called more than once");
  }
  if (prop.get()) {
    particleIdProperty = prop->id;
  }
  else {
    particleIdProperty = addProperty< Property< ParticleId > >();
  }
}

void Storage::setPositionProperty(boost::shared_ptr< Property< Real3D > >prop) {
  if (prop.get()) {
    particlePosProperty = prop->id;
  }
  else {
    particlePosProperty = addProperty< Property< Real3D > >();
  }
  positionPropertyModified();
}

IdPropertyHandle
Storage::getIdPropertyHandle()
{ 
  if (particleIdProperty == PropertyId()) {
    throw std::runtime_error("Storage::getIdPropertyHandle: no id-property configured");
  }
  return getTupleVector().getProperty< ParticleId >(particleIdProperty);
}

PosPropertyHandle
Storage::getPositionPropertyHandle()
{
  if (particlePosProperty == PropertyId()) {
    throw std::runtime_error("Storage::getPositionPropertyHandle: no pos-property configured");
  }
  return getTupleVector().getProperty< Real3D >(particlePosProperty);
}

//////////////////////////////////////////////////
// REGISTRATION WITH PYTHON
//////////////////////////////////////////////////

void Storage::registerPython() {
  using namespace espresso::python;

  class_< Storage, boost::noncopyable, bases < particles::Set, pairs::Set > >
    ("storage_Storage", no_init)
    .def("addParticle", &Storage::_addParticle)
    .def("deleteParticle", &Storage::deleteParticle)
    .def("setIdProperty", &Storage::setIdProperty)
    .def("setPositionProperty", &Storage::setPositionProperty)
    ;
}

