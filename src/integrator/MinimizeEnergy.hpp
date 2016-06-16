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

// ESPP_CLASS
#ifndef _INTEGRATOR_MINIMIZEENERGY_HPP
#define _INTEGRATOR_MINIMIZEENERGY_HPP

#include "logging.hpp"
#include "Extension.hpp"
#include "boost/signals2.hpp"
#include "esutil/Timer.hpp"
#include "Real3D.hpp"
#include "Int3D.hpp"

namespace espressopp {
namespace integrator {

class MinimizeEnergy : public MDIntegrator {
 public:
  MinimizeEnergy(shared_ptr<class espressopp::System> system);
  virtual ~MinimizeEnergy();

  void run(int nsteps);

  /** Register this class so it can be used from Python. */
  static void registerPython();
 private:
  real steepestDescentStep();

  void updateForces();

  void calcForces();

  // Getters and setters.

  real f_max_;
  real gamma_;
  longint max_steps_;
  real max_displacement_;

  bool resort_flag_;  //!< true implies need for resort of particles
  real max_dist_;

  static LOG4ESPP_DECL_LOGGER(theLogger);
};


}  // end namespace integrator
}  // end namespace espressopp