/* python.hpp has to be included before any system headers according
   to Python API to avoid redefining _POSIX_C_SOURCE by Python.h */
#include <python.hpp>
#include "bindings.hpp"

#include <Particle.hpp>
#include <ParticleGroup.hpp>
#include <System.hpp>
#include <VerletList.hpp>
#include <VerletListAdress.hpp>
#include <FixedPairList.hpp>
#include <FixedTripleList.hpp>
#include <FixedQuadrupleList.hpp>
#include <FixedTupleList.hpp>
#include <Real3D.hpp>
#include <Tensor.hpp>
#include <Int3D.hpp>
#include <esutil/PyLogger.hpp>

#include <esutil/bindings.hpp>
#include <bc/bindings.hpp>
#include <storage/bindings.hpp>
#include <integrator/bindings.hpp>
#include <interaction/bindings.hpp>
#include <analysis/bindings.hpp>

void espresso::registerPython() {
  espresso::Particle::registerPython();
  espresso::ParticleGroup::registerPython();
  espresso::System::registerPython();
  espresso::VerletList::registerPython();
  espresso::VerletListAdress::registerPython();
  espresso::FixedPairList::registerPython();
  espresso::FixedTripleList::registerPython();
  espresso::FixedQuadrupleList::registerPython();
  espresso::FixedTupleList::registerPython();
  espresso::Real3D::registerPython();
  espresso::Tensor::registerPython();
  espresso::Int3D::registerPython();

  espresso::esutil::registerPython();
  espresso::bc::registerPython();
  espresso::storage::registerPython();
  espresso::integrator::registerPython();
  espresso::interaction::registerPython();
  espresso::analysis::registerPython();

  log4espp::PyLogger::registerPython();
}
