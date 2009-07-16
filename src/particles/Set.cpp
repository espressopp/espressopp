#include "python.hpp"
#include "boost/bind.hpp"

#include "particles/Set.hpp"
#include "particles/Computer.hpp"


using namespace espresso::particles;

// Implementation of a generic isMember

namespace {
  class Found {};
  void findMember(const ConstParticleHandle p1, const ConstParticleHandle p2) {
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

void
Set::foreach(Computer &computer) {
  computer.prepare();
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
  computer.prepare();
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

//////////////////////////////////////////////////
// REGISTRATION WITH PYTHON
//////////////////////////////////////////////////
void
Set::registerPython() {
  using namespace espresso::python;
  
  void (Set::*py_foreach)(const Computer::SelfPtr computer) 
    = &Set::foreach;
  
  class_< Set, boost::noncopyable >
    ("particles_Set", no_init)
    .def("foreach", py_foreach)
    .def("isMember", &Set::isMember)
    ;
}


