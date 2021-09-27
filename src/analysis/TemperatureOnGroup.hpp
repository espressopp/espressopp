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
#ifndef _ANALYSIS_TemperatureOnGroup_HPP
#define _ANALYSIS_TemperatureOnGroup_HPP

#include "types.hpp"
#include "Observable.hpp"
#include "storage/Storage.hpp"
#include "ParticleGroup.hpp"

namespace espressopp
{
namespace analysis
{
/** Class to compute the Temperature in ParticleGroup. */
class TemperatureOnGroup : public Observable
{
public:
    static void registerPython();

    TemperatureOnGroup(std::shared_ptr<System> system, std::shared_ptr<ParticleGroup> pg)
        : Observable(system), particle_group_(pg), eKin_(0.0)
    {
        result_type = real_scalar;
    }

    virtual ~TemperatureOnGroup() {}

    real compute_real() const
    {
        int myN, systemN;
        real sumT = 0.0;
        real v2sum = 0.0;
        myN = 0;
        systemN = 0;

        for (ParticleGroup::iterator it = particle_group_->begin(); it != particle_group_->end();
             it++)
        {
            Real3D vel = it->velocity();
            v2sum += it->mass() * (vel * vel);
            myN += 1;
        }

        mpi::all_reduce(*getSystem()->comm, v2sum, sumT, std::plus<real>());
        mpi::all_reduce(*getSystem()->comm, myN, systemN, std::plus<int>());

        eKin_ = 0.5 * sumT;
        return sumT / (3.0 * systemN);
    }

    real getEkin() const { return eKin_; }

private:
    std::shared_ptr<ParticleGroup> particle_group_;
    mutable real eKin_;
};

}  // end namespace analysis
}  // end namespace espressopp
#endif
