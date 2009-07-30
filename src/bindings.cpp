#include "bindings.hpp"
#include "Real3D.hpp"
#include "Property.hpp"
#include <esutil/Collectives.hpp>
#include <storage/bindings.hpp>
#include <particles/bindings.hpp>
#include <bc/bindings.hpp>
#include <pairs/bindings.hpp>
#include <potential/bindings.hpp>
#include <integrator/bindings.hpp>
#include <thermostat/bindings.hpp>
#include <esutil/PyLogger.hpp>

void espresso::registerPython() {

  espresso::registerPythonReal3D();
  espresso::registerPythonParticle();
  espresso::registerPythonProperty();
  espresso::esutil::Collectives::registerPython();
  espresso::particles::registerPython();
  espresso::storage::registerPython();
  espresso::bc::registerPython();
  espresso::pairs::registerPython();
  espresso::potential::registerPython();
  espresso::integrator::registerPython();
  espresso::thermostat::registerPython();

  log4espp::PyLogger::registerPython();
}
