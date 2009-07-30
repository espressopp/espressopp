#include "bindings.hpp"
#include "types.hpp"

#include "Storage.hpp"

using namespace espresso::storage;
 
void espresso::storage::registerPython() {
  Storage::registerPython();
}
