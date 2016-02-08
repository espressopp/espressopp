/*
  Copyright (C) 2012,2013
      Max Planck Institute for Polymer Research
  Copyright (C) 2008,2009,2010,2011
      Max-Planck-Institute for Polymer Research & Fraunhofer SCAI
  
  This file is part of ESPResSo++.
  
  ESPResSo++ is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.
  
  ESPResSo++ is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.
  
  You should have received a copy of the GNU General Public License
  along with this program.  If not, see <http://www.gnu.org/licenses/>. 
*/

#include <python.hpp>
#include "Particle.hpp"

namespace espressopp {
  void 
  Particle::
  registerPython() {
    using namespace python;
    class_< Particle >("_TmpParticle", no_init)
      .add_property("id", &Particle::getId)
      .add_property("type", &Particle::getType, &Particle::setType)
      .add_property("mass", &Particle::getMass, &Particle::setMass)
      .add_property("pos", &Particle::getPos, &Particle::setPos)
      .add_property("f", &Particle::getF, &Particle::setF)
      .add_property("v", &Particle::getV, &Particle::setV)
      .add_property("q", &Particle::getQ, &Particle::setQ)
      .add_property("radius", &Particle::getRadius, &Particle::setRadius)
      .add_property("fradius", &Particle::getFRadius, &Particle::setFRadius)
      .add_property("vradius", &Particle::getVRadius, &Particle::setVRadius)
      .add_property("imageBox", &Particle::getImageBox, &Particle::setImageBox)
      .add_property("isGhost", &Particle::getGhostStatus, &Particle::setGhostStatus)
      .add_property("lambda_adr", &Particle::getLambda, &Particle::setLambda)
      .add_property("lambda_adrd", &Particle::getLambdaDeriv, &Particle::setLambdaDeriv)
      .add_property("state", &Particle::getState, &Particle::setState)
      .add_property("res_id", &Particle::getResId, &Particle::setResId)
      .add_property("extVar", &Particle::getExtVar, &Particle::setExtVar)
      .add_property("drift_f", &Particle::getDrift, &Particle::setDrift)
      ;
  }

  void ParticleProperties::registerPython() {
    using namespace python;
    class_<ParticleProperties,
        boost::shared_ptr<ParticleProperties> >("_ParticleProperties")
      .add_property(
          "type",
          make_getter(&ParticleProperties::type),
          make_setter(&ParticleProperties::type))
      .add_property(
          "mass",
          make_getter(&ParticleProperties::mass),
          make_setter(&ParticleProperties::mass))
      .add_property(
          "q",
          make_getter(&ParticleProperties::q),
          make_setter(&ParticleProperties::q))
      .add_property(
          "state",
          make_getter(&ParticleProperties::state),
          make_setter(&ParticleProperties::state))
      .add_property(
          "res_id",
          make_getter(&ParticleProperties::res_id),
          make_setter(&ParticleProperties::res_id))
      .def(
          "init",
          &ParticleProperties::init
      );
  }
}
