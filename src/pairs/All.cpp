#include "All.hpp"

#include <exception>
#include <cmath>

#include "particles/Computer.hpp"
#include "python.hpp"

using namespace espresso;
using namespace espresso::bc;
using namespace espresso::pairs;
using namespace espresso::storage;

using namespace boost;

All::All(particles::Set::SelfPtr _set)
  : set(_set) {}

namespace {
  struct CheckingComputer: public pairs::Computer {
    CheckingComputer(Computer &_computer, particles::Set &_set):
      computer(_computer), set(_set) {}

    void prepare(Storage::SelfPtr s1, Storage::SelfPtr s2)
    { computer.prepare(s1, s2); }
    bool apply(const Real3D &d,
	       ParticleHandle p1,
	       ParticleHandle p2) {
      if (set.contains(p1) && set.contains(p2)) {
	return computer.apply(d, p1, p2);
      }
      return true;
    }

    pairs::Computer &computer;
    particles::Set &set;
  };
}

bool All::foreachPairApply(Computer &computer) {
  CheckingComputer checkingComputer(computer, *set);
  return set->getStorage()->foreachPair(checkingComputer);
}

Storage::SelfPtr All::getLeftStorage()
{ return set->getStorage(); }

Storage::SelfPtr All::getRightStorage()
{ return set->getStorage(); }


//////////////////////////////////////////////////
// REGISTRATION WITH PYTHON
//////////////////////////////////////////////////

void 
All::registerPython() {
  using namespace espresso::python;
  class_< All, bases< Set > >
    ("pairs_All", 
     init< particles::Set::SelfPtr >())
    ;
}
