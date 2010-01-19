/* has to be included before any system headers according to Python API
   to avoid redefining _POSIX_C_SOURCE by Python.h */
#include <python.hpp>
#include <esutil/Collectives.hpp>
#include <BC.hpp>
#include <esutil/PyLogger.hpp>

#include "bindings.hpp"

void espresso::registerPython() {
  espresso::esutil::Collectives::registerPython();
//  espresso::bc::registerPython();
  espresso::BC::registerPython();
  log4espp::PyLogger::registerPython();
}
