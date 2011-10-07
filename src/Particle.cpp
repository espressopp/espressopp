#include <python.hpp>
#include "Particle.hpp"

namespace espresso {
  void 
  Particle::
  registerPython() {
    using namespace python;
    class_< Particle >("_TmpParticle", no_init)
      .add_property("id", &Particle::getId)
      .add_property("type", &Particle::getType, &Particle::setType)
      .add_property("mass", &Particle::getMass, &Particle::setMass)
      .add_property("pos", &Particle::getPos, &Particle::setPos)
      .add_property("image", &Particle::getImage)
      .add_property("f", &Particle::getF, &Particle::setF)
      .add_property("v", &Particle::getV, &Particle::setV)
      .add_property("q", &Particle::getQ, &Particle::setQ)
      .add_property("imageBox", &Particle::getImageBox, &Particle::setImageBox)
      .add_property("isGhost", &Particle::getGhostStatus, &Particle::setGhostStatus)
      ;
  }
}
