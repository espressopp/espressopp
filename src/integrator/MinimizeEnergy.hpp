/*
  Copyright (C) 2016
      Max Planck Institute for Polymer Research
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

// ESPP_CLASS

#ifndef _INTEGRATOR_MINIMIZEENERGY_HPP
#define _INTEGRATOR_MINIMIZEENERGY_HPP

#include "logging.hpp"
#include "Real3D.hpp"
#include "iterator/CellListIterator.hpp"
#include "SystemAccess.hpp"
#include "storage/Storage.hpp"
#include "interaction/Interaction.hpp"
#include "interaction/Potential.hpp"

namespace espressopp {
namespace integrator {

// Adopted from espressomd: src/core/minimize_energy.cpp

class MinimizeEnergy : public SystemAccess {
 public:
  MinimizeEnergy(shared_ptr<class espressopp::System> system,
                 real gamma,
                 real ftol,
                 real max_displacement,
		 bool variable_step_flag);
  virtual ~MinimizeEnergy();

  bool run(int max_steps, bool verbose);

  /** Register this class so it can be used from Python. */
  static void registerPython();
 private:
  void steepestDescentStep();
  void updateForces();

  // Getters
  real getFMax() {
    return sqrt(f_max_sqr_);
  }

  real getDpMax() {
    return sqrt(dp_sqr_max_);
  }

  // Params
  real gamma_;
  real max_displacement_;  // Maximum displacement on particle.
  real ftol_sqr_;  // Force limit, when maximum force is lower then stop.

  real f_max_sqr_;  // Maximum force on particles.
  real dp_sqr_max_;   // Maximum particle displacement.
  real dp_MAX; // Summation of maximum particle displacement.

  bool variable_step_flag_;  //!< true implies that gamma is adjusted to the force strength. 

  bool resort_flag_;  //!< true implies need for resort of particles

  longint nstep_;

  static LOG4ESPP_DECL_LOGGER(theLogger);
};


}  // end namespace integrator
}  // end namespace espressopp
#endif
