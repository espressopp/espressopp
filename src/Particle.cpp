/*
  Copyright (C) 2016-2018
    Jakub Krajniak (jkrajniak at gmail.com)
  Copyright (C) 2012,2013
      Max Planck Institute for Polymer Research
  Copyright (C) 2008-2011
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
      .add_property("pib", &Particle::getPib, &Particle::setPib)
      .add_property("mass", &Particle::getMass, &Particle::setMass)
      .add_property("pos", &Particle::getPos, &Particle::setPos)
      .add_property("f", &Particle::getF, &Particle::setF)
      .add_property("fm", &Particle::getFm, &Particle::setFm)
      .add_property("v", &Particle::getV, &Particle::setV)
      .add_property("modepos", &Particle::getModepos, &Particle::setModepos)
      .add_property("modemom", &Particle::getModemom, &Particle::setModemom)
      .add_property("q", &Particle::getQ, &Particle::setQ)
      .add_property("radius", &Particle::getRadius, &Particle::setRadius)
      .add_property("fradius", &Particle::getFRadius, &Particle::setFRadius)
      .add_property("vradius", &Particle::getVRadius, &Particle::setVRadius)
      .add_property("imageBox", &Particle::getImageBox, &Particle::setImageBox)
      .add_property("isGhost", &Particle::getGhostStatus, &Particle::setGhostStatus)
      .add_property("lambda_adr", &Particle::getLambda, &Particle::setLambda)
      .add_property("lambda_adrd", &Particle::getLambdaDeriv, &Particle::setLambdaDeriv)
      .add_property("varmass", &Particle::getVarmass, &Particle::setVarmass)
      .add_property("state", &Particle::getState, &Particle::setState)
      .add_property("res_id", &Particle::getResId, &Particle::setResId)
      .add_property("extVar", &Particle::getExtVar, &Particle::setExtVar)
      .add_property("drift_f", &Particle::getDrift, &Particle::setDrift)
      ;
  }
}
