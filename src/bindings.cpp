#include "bindings.hpp"
#include "Real3D.hpp"
#include <hello/bindings.hpp>
#include <particles/bindings.hpp>
#include <interaction/bindings.hpp>
#include <integrator/bindings.hpp>

void espresso::registerPython() {
  espresso::registerPythonReal3D();

  espresso::hello::registerPython();
  espresso::interaction::registerPython();
  espresso::particles::registerPython();
  espresso::integrator::registerPython();
}
