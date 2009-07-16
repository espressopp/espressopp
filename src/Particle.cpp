#include "python.hpp"
#include "Particle.hpp"

using namespace espresso;
using namespace espresso::python;

/// espresso::python ResultConverter for ParticleId
struct ParticleIdToInteger {
  static PyObject* convert(const ParticleId &id) {
    return incref(object(size_t(id)).ptr());
  }
};

void espresso::registerPythonParticle()
{
  /* rather than exporting the ParticleId class, we make sure it gets
     converted to a Python integer */
  implicitly_convertible<size_t, ParticleId>();
  to_python_converter<ParticleId, ParticleIdToInteger, false>();
}
