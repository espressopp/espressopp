#include "bindings.hpp"
#include "Real3D.hpp"
#include "Property.hpp"
#include <hello/bindings.hpp>
#include <particles/bindings.hpp>
#include <interaction/bindings.hpp>
#include <integrator/bindings.hpp>
#include <esutil/PyLogger.hpp>

void espresso::registerPython() {

  espresso::registerPythonReal3D();
  espresso::registerPythonParticle();
  espresso::registerPythonProperties();

  espresso::hello::registerPython();
  espresso::interaction::registerPython();
  espresso::particles::registerPython();
  espresso::integrator::registerPython();

  log4espp::PyLogger::registerPython();
}
