#include "bindings.hpp"
#include "types.hpp"

#include <boost/python.hpp>

#include <particles/Storage.hpp>
#include <particles/Computer.hpp>
#include <particles/PythonComputer.hpp>

using namespace boost::python;

using namespace espresso::particles;
 
void espresso::particles::registerPython() {
  /* class PropertyId has to be exported otherwise addProperty
     cannot return a Python object
     no_init: we do not need a constructor in Python
     copyable: otherwise we cannot assign it in Python */
  class_<PropertyId>("particles_PropertyId", no_init);

  /* Computer has to be exported otherwise Python will not accept
     something as argument for foreach in Storage.
     noncopyable: must be specified as it is an abstract class */
  class_<Computer, boost::noncopyable>("particles_Computer", no_init);

  PythonComputer::registerPython();

  Storage::registerPython();
}
