#include "bindings.hpp"

#include <hello/bindings.hpp>
#include <interaction/bindings.hpp>
#include <particles/bindings.hpp>
#include <esutil/bindings.hpp>

void espresso::registerPython() {
  espresso::hello::registerPython();
  espresso::particles::registerPython();
  espresso::interaction::registerPython();
  espresso::esutil::registerPython();
}
