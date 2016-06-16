/*
  Copyright (C) 2016
      Jakub Krajniak (jkrajniak at gmail.com)

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

#include "python.hpp"
#include "MinimizeEnergy.hpp"
#include "storage/Storage.hpp"


namespace espressopp {
namespace integrator {

LOG4ESPP_LOGGER(MinimizeEnergy::theLogger, "MinimizeEnergy");

MinimizeEnergy::MinimizeEnergy(shared_ptr<System> system)
    : MDIntegrator(system) {
  LOG4ESPP_INFO(theLogger, "construct MinimizeEnergy");
}

MinimizeEnergy::~MinimizeEnergy() {
  LOG4ESPP_INFO(theLogger, "free MinimizeEnergy");
}

void MinimizeEnergy::run(int nsteps) {
  int nResorts = 0;
  System &system = getSystemRef();
  storage::Storage &storage = *system.storage;
  real skin_half = 0.5 * system.getSkin();

  // signal
  runInit();

  // Before start make sure that particles are on the right processor
  if (resort_flag_) {
    storage.decompose();
    max_dist_ = 0.0;
    resort_flag_ = false;
  }

  LOG4ESPP_INFO(theLogger, "recalc forces before starting main integration loop");

  // signal before recalculate forces
  recalc1();
  updateForces();
  // signal after recalculate forces
  recalc2();

  LOG4ESPP_INFO(theLogger, "starting main integration loop (nsteps=" << nsteps << ")");

  for (int i = 0; i < nsteps; i++) {
    LOG4ESPP_INFO(theLogger, "Next step " << i << " of " << nsteps << " starts");
    // signal before new position
    befIntP();

    max_dist_ += steepestDescentStep();

    // signal after new position
    aftIntP();

    LOG4ESPP_INFO(theLogger, "maxDist = " << max_dist_ << ", skin/2 = " << skin_half);

    resort_flag_ = max_dist_ > skin_half;

    if (resort_flag_) {
      LOG4ESPP_INFO(theLogger, "step " << i << ": resort particles");
      storage.decompose();
      max_dist_ = 0.0;
      resort_flag_ = false;
      nResorts++;
    }

    updateForces();

    // two signals, before and after calculate velocities but we don't
    // calculate velocities.
    befIntV();
    aftIntV();
  }
  LOG4ESPP_INFO(theLogger, "finished run");
}
void MinimizeEnergy::updateForces() {
  LOG4ESPP_INFO(theLogger,
                "update ghosts, calculate forces and collect ghost forces")
  real time;
  storage::Storage& storage = *getSystemRef().storage;
  storage.updateGhosts();

  calcForces();
  storage.collectGhostForces();

  // signal
  aftCalcF();
}

void MinimizeEnergy::calcForces() {

}

void VelocityVerlet::initForces() {
  // forces are initialized for real + ghost particles

  System& system = getSystemRef();
  CellList localCells = system.storage->getLocalCells();

  LOG4ESPP_INFO(theLogger, "init forces for real + ghost particles");

  for(CellListIterator cit(localCells); !cit.isDone(); ++cit) {
    cit->force() = 0.0;
    cit->drift() = 0.0;   // Can in principle be commented, when drift is not used.
  }
}

}  // end namespace integrator
}  // end namespace espressopp