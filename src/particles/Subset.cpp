#include "python.hpp"
#include "boost/bind.hpp"

#include "particles/Subset.hpp"
#include "particles/Computer.hpp"

using namespace espresso::particles;
using namespace espresso::storage;

Subset::Subset(const Set::SelfPtr superset)
  : storage(superset->getStorage()) {}

Storage::SelfPtr
Subset::getStorage() { return storage; }

//////////////////////////////////////////////////
// REGISTRATION WITH PYTHON
//////////////////////////////////////////////////
void
Subset::registerPython() {
  using namespace espresso::python;
  
  class_< Subset, bases < Set >, boost::noncopyable >
    ("particles_Subset", no_init)
    ;
}


