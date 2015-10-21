/*
  Copyright (C) 2015
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
#ifndef _ANALYSIS_POTENTIALENERGY_HPP
#define _ANALYSIS_POTENTIALENERGY_HPP

#include "types.hpp"
#include "Observable.hpp"
#include "interaction/Interaction.hpp"

namespace espressopp {
namespace analysis {

class PotentialEnergy : public Observable {
 public:
  PotentialEnergy(shared_ptr<System> system, shared_ptr<interaction::Interaction> interaction,
                  bool compute_at):
      Observable(system), interaction_(interaction), compute_at_(compute_at) {
    result_type = real_scalar;
    compute_global_ = false;
  }
  PotentialEnergy(shared_ptr<System> system, shared_ptr<interaction::Interaction> interaction)
      : Observable(system), interaction_(interaction) {
    result_type = real_scalar;
    compute_global_ = true;
  }
  ~PotentialEnergy() {}
  real compute_real() const;

  static void registerPython();
 private:
  shared_ptr<interaction::Interaction> interaction_;
  bool compute_at_;  // set to true then computeEnergyAA, otherwise computeEnergyCG
  bool compute_global_;
};
}  // end namespace analysis
}  // end namespace espressopp

#endif
