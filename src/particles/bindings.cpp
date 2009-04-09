#include "bindings.hpp"
#include "types.hpp"

#include <boost/python.hpp>

#include <particles/Storage.hpp>
#include <particles/Computer.hpp>
#include <particles/PythonComputer.hpp>

using namespace boost::python;

using namespace espresso::particles;
 
void espresso::particles::registerPython() {
  class_<Computer, boost::noncopyable>("particles_Computer", no_init);

  PythonComputer::registerPython();

  Storage::registerPython();
}
