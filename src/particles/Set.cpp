#include "particles/Set.hpp"
#include "particles/Computer.hpp"
#include "particles/ParticleHandle.hpp"
#include "particles/Storage.hpp"

#include <boost/bind.hpp>
#include <python.hpp>

// Implementation of a generic isMember

using namespace espresso::particles;

namespace {
  class Found {};
  void findMember(const ConstParticleHandle p1, 
		  const ConstParticleHandle p2) {
    if (p1 == p2) throw new Found();
  }
}

bool
Set::isMember(const ConstParticleHandle p) const {
  try {
    foreach(boost::bind(findMember, p, _1));
  } catch (Found) {
    return true;
  }
  return false;
}

bool
Set::isMember(ParticleId pid) {
  return isMember(getStorage()->getParticleHandle(pid));
}

void
Set::foreach(Computer &computer) {
  computer.prepare(getStorage());
  try {
    foreach(boost::bind(&Computer::apply, &computer, _1));
  } catch (ForeachBreak) {
    computer.finalize();
    throw;
  } 
  computer.finalize();
}
  
void
Set::foreach(ConstComputer &computer) const {
  computer.prepare(getStorage());
  try {
    foreach(boost::bind(&ConstComputer::apply, &computer, _1));
  } catch (ForeachBreak &exc) {
    computer.finalize();
    throw;
    }
  computer.finalize();
}

void
Set::foreach(const Computer::SelfPtr computer) 
{ foreach(*computer); }

const Storage::SelfPtr
Set::getStorage() const { 
  return const_pointer_cast< Storage, const Storage >(getStorage()); 
}

//////////////////////////////////////////////////
// REGISTRATION WITH PYTHON
//////////////////////////////////////////////////
void
Set::registerPython() {
  using namespace espresso::python;
  
  void (Set::*pyForeach)(const Computer::SelfPtr computer) 
    = &Set::foreach;

  bool (Set::*pyIsMember)(ParticleId pid)
    = &Set::isMember;

  class_< Set, boost::noncopyable >
    ("particles_Set", no_init)
    .def("foreach", pyForeach)
    .def("isMember", pyIsMember)
    ;
}
