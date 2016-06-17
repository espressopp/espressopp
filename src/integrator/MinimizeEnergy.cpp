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

namespace espressopp {
namespace integrator {

using namespace interaction;
using namespace iterator;
using namespace esutil;

LOG4ESPP_LOGGER(MinimizeEnergy::theLogger, "MinimizeEnergy");

MinimizeEnergy::MinimizeEnergy(shared_ptr<System> system,
                               real gamma,
                               real ftol,
                               real max_displacement)
    : SystemAccess(system), gamma_(gamma), ftol_(ftol),
      max_displacement_(max_displacement) {
  LOG4ESPP_INFO(theLogger, "construct MinimizeEnergy");

  nstep_ = 0;
}

MinimizeEnergy::~MinimizeEnergy() {
  LOG4ESPP_INFO(theLogger, "free MinimizeEnergy");
}

void MinimizeEnergy::run(int max_steps, bool verbose) {
  System &system = getSystemRef();
  storage::Storage &storage = *system.storage;
  real skin_half = 0.5 * system.getSkin();
  dp_sqr_max_ = 0.0;
  f_max_ = std::numeric_limits<real>::max();

  // Before start make sure that particles are on the right processor
  if (resort_flag_) {
    storage.decompose();
    max_dist_ = 0.0;
    resort_flag_ = false;
  }

  LOG4ESPP_INFO(theLogger, "recalc forces before starting main integration loop");

  updateForces();

  LOG4ESPP_INFO(theLogger,
                "starting energy minimalization loop (iters=" << max_steps << ")");

  if (verbose) {
    std::cout << "Minimize energy" << std::endl;
    std::cout << "  current force_max = " << f_max_ << std::endl;
    std::cout << "  f_tol = " << ftol_ << std::endl;
    std::cout << "  max_steps = " << max_steps << std::endl;
    std::cout << "  max displacement = " << max_displacement_ << std::endl;
  }
  int iters = 0;
  for (; iters < max_steps && f_max_ > ftol_; iters++) {
    steepestDescentStep();

    LOG4ESPP_DEBUG(theLogger, "step " << iters << " max_force=" << f_max_ << " displacement=" << dp_sqr_max_);

    resort_flag_ = sqrt(dp_sqr_max_) > skin_half;

    if (resort_flag_) {
      storage.decompose();
      max_dist_ = 0.0;
      resort_flag_ = false;
    }

    updateForces();

    if (verbose)
      std::cout << nstep_ << ": f_max=" << f_max_ << " max_dp=" << dp_sqr_max_ << std::endl;

    nstep_++;
  }

  if (verbose) {
    std::cout << "Minimize energy finished" << std::endl;
    std::cout << "  current force_max = " << f_max_ << std::endl;
    std::cout << "  run for steps = " << iters << std::endl;
    std::cout << "   max displacement = " << dp_sqr_max_ << std::endl;
  }

  LOG4ESPP_INFO(theLogger,
                "finished run, f_max_=" << f_max_ << " max_displ=" << dp_sqr_max_);
}


void MinimizeEnergy::updateForces() {
  LOG4ESPP_INFO(theLogger,
                "update ghosts, calculate forces and collect ghost forces");
  System& system = getSystemRef();

  system.storage->updateGhosts();

  LOG4ESPP_INFO(theLogger, "calculate forces");
  // First initialize
  CellList localCells = system.storage->getLocalCells();
  for(CellListIterator cit(localCells); !cit.isDone(); ++cit) {
    cit->force() = 0.0;
  }

  // Calculate force component from all interactions.
  const InteractionList& srIL = system.shortRangeInteractions;

  for (size_t i = 0; i < srIL.size(); i++) {
    LOG4ESPP_INFO(theLogger, "compute forces for srIL " << i << " of " << srIL.size());
    srIL[i]->addForces();
  }
  // Collect forces from ghost particles.
  system.storage->collectGhostForces();

  // Get max force in the system.
  real f_max = -std::numeric_limits<real>::max();
  CellList realCells = system.storage->getRealCells();
  for(CellListIterator cit(realCells); !cit.isDone(); ++cit) {
    f_max = std::max(f_max, cit->force().sqr());
  }
  mpi::all_reduce(*system.comm, f_max, f_max_, boost::mpi::maximum<real>());
}

template <typename T> int sgn(T val) {
  return (T(0) < val) - (val < T(0));
}

void MinimizeEnergy::steepestDescentStep() {
  LOG4ESPP_INFO(theLogger, "steepestDescent step");
  System& system = getSystemRef();

  real f_sqr, dp, dp_sqr;
  real f_max = std::numeric_limits<real>::min();
  real dp_sqr_max = std::numeric_limits<real>::min();

  // Iterate over only real particles.
  CellList realCells = system.storage->getRealCells();
  for (CellListIterator cit(realCells); !cit.isDone(); ++cit) {
    // Magnitiude of the force.
    f_sqr = cit->force().sqr();
    for (int i = 0; i < 3; i++) {   // Perhaps it can be done better.
      dp = gamma_ * cit->force()[i];
      if (fabs(dp) > max_displacement_)
        dp = sgn<real>(dp)*max_displacement_;
      dp_sqr += dp*dp;

      // Update position component by dp.
      cit->position()[i] += dp;
    }
    f_max = std::max(f_max, f_sqr);
    dp_sqr_max = std::max(dp_sqr_max, dp_sqr);
  }

  mpi::all_reduce(*system.comm, dp_sqr_max, dp_sqr_max_, boost::mpi::maximum<real>());
}

void MinimizeEnergy::registerPython() {
  using namespace espressopp::python;

  // Note: use noncopyable and no_init for abstract classes
  class_<MinimizeEnergy, boost::noncopyable>
    ("integrator_MinimizeEnergy", init<shared_ptr<System>, real, real, real>())
      .add_property("f_max", make_getter(&MinimizeEnergy::f_max_))
      .add_property("displacement", make_getter(&MinimizeEnergy::dp_sqr_max_))
      .add_property("step", make_getter(&MinimizeEnergy::nstep_), make_setter(&MinimizeEnergy::nstep_))
      .def("run", &MinimizeEnergy::run);
}

}  // end namespace integrator
}  // end namespace espressopp