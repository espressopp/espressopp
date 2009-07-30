#include "bindings.hpp"
#include "types.hpp"

#include "Computer.hpp"
#include "Set.hpp"

using namespace espresso::particles;
 
void espresso::particles::registerPython() {
  Set::registerPython();
  Computer::registerPython();
}
