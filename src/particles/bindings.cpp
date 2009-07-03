#include "bindings.hpp"
#include "types.hpp"

#include <boost/python.hpp>

#include <particles/Storage.hpp>
#include <particles/Computer.hpp>
#include <particles/All.hpp>

using namespace boost::python;

using namespace espresso::particles;
 
void espresso::particles::registerPython() {
  Computer::registerPython();
  Storage::registerPython();
  Set::registerPython();
  All::registerPython();
}
