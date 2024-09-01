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
#ifndef _INTEGRATOR_VelocityVerletHybrid_HPP
#define _INTEGRATOR_VelocityVerletHybrid_HPP

#include "types.hpp"
#include "MDIntegrator.hpp"
#include "esutil/Timer.hpp"
#include "FixedVSList.hpp"
#include <boost/signals2.hpp>

namespace espressopp
{
namespace integrator
{
/** Velocity Verlet Integrator */
class VelocityVerletHybrid : public MDIntegrator
{
public:
    VelocityVerletHybrid(shared_ptr<class espressopp::System> system,
                         shared_ptr<FixedVSList> vs_list);

    virtual ~VelocityVerletHybrid();

    void run(int nsteps);

    /** Load timings in array to export to Python as a tuple. */
    void loadTimers(std::vector<real> &return_vector, std::vector<std::string> &return_labels);

    /** Clean up all timers.*/
    void resetTimers();

    /** Register this class so it can be used from Python. */
    static void registerPython();

protected:
    bool resortFlag;  //!< true implies need for resort of particles
    real maxDist;

    real maxCut;

    real integrate1();

    void integrate2();

    void initForces();

    void updateForces();

    void calcForces();

    esutil::WallTimer timeIntegrate;  //!< used for timing

    // variables that keep time information about different phases
    real timeRun;
    real timeLost;
    real timeForce;
    std::vector<real> timeForceComp;
    real timeComm1;
    real timeComm2;
    real timeInt1;
    real timeInt2;
    real timeResort;

    // Signal timers
    real timeRunInitS;
    real timeRecalc1S;
    real timeRecalc2S;
    real timeBefIntPS;
    real timeAftIntPS;
    real timeAftInitFS;
    real timeAftCalcFS;
    real timeBefIntVS;
    real timeAftIntVS;

    real timeVS;
    real timeVSvel;
    real timeVSdistrF;

    static LOG4ESPP_DECL_LOGGER(theLogger);

private:
    void distributeVSforces();
    real updateVS();
    void updateVS_vel();

    real skinHalf;
    shared_ptr<FixedVSList> vs_list_;
    void checkConsistence(int postfix = 0);
};
}  // namespace integrator
}  // namespace espressopp

#endif
