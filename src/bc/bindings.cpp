#include "bindings.hpp"

#include "BC.hpp"
#include "PBC.hpp"

using namespace espresso::bc;

void 
espresso::bc::registerPython() {
   BC::registerPython();
   PBC::registerPython();
}
