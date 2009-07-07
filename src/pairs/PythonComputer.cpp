#include "pairs/PythonComputer.hpp"
#include "Particle.hpp"
#include "python.hpp"

using namespace espresso::pairs;
using namespace espresso::particles;

void PythonComputer::operator()(const Real3D &dist,
                                const ParticleHandle p1,
                                const ParticleHandle p2)
{
  ParticleId id1 = set->getStorage()->getParticleId(p1);
  ParticleId id2 = set->getStorage()->getParticleId(p2);
  get_override("each")(id1, id2);
}

void PythonComputer::registerPython() {
  using namespace espresso::python;
  
  class_< PythonComputer, bases< Computer > >
    ("pairs_PythonComputer", init< pairs::Set::SelfPtr > ());
}
