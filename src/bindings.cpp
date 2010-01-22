/* has to be included before any system headers according to Python API
   to avoid redefining _POSIX_C_SOURCE by Python.h */
#include <python.hpp>
#include <esutil/bindings.hpp>
#include <bc/bindings.hpp>
#include <Real3D.hpp>
#include <BC.hpp>
#include <esutil/PyLogger.hpp>

#include "bindings.hpp"

void espresso::registerPython() {

  espresso::Real3D::registerPython();
  espresso::BC::registerPython();

  espresso::esutil::registerPython();
  espresso::bc::registerPython();

  log4espp::PyLogger::registerPython();
}
