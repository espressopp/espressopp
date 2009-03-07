#include "bindings.hpp"

#include <hello/bindings.hpp>
#include <interaction/bindings.hpp>

void espresso::registerPython() {
  espresso::hello::registerPython();
  espresso::interaction::registerPython();
}
