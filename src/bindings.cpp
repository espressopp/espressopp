#include "bindings.hpp"
#include <esutil/Collectives.hpp>
#include <esutil/PyLogger.hpp>

void espresso::registerPython() {
  espresso::esutil::Collectives::registerPython();
  log4espp::PyLogger::registerPython();
}
