/*
  Copyright (C) 2020-2022
      Max Planck Institute for Polymer Research & JGU Mainz

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
#ifndef HPX4ESPP_INTEGRATOR_VELOCITYVERLET_HPP
#define HPX4ESPP_INTEGRATOR_VELOCITYVERLET_HPP

#include "hpx4espp/storage/StorageHPX.hpp"
#include "hpx4espp/integrator/MDIntegratorHPX.hpp"
#include "types.hpp"
#include "integrator/MDIntegrator.hpp"
#include "esutil/Timer.hpp"
#include <boost/signals2.hpp>

namespace espressopp
{
namespace hpx4espp
{
namespace integrator
{
/// Velocity Verlet Integrator
class VelocityVerlet : public MDIntegratorHPX
{
public:
    typedef hpx4espp::integrator::MDIntegratorHPX baseClass;

    VelocityVerlet(shared_ptr<class SystemHPX> system,
                   shared_ptr<class storage::StorageHPX> storageHPX);

    virtual void run(int nsteps);

    real getTimeStep() { return baseClass::getTimeStep(); }

    /** Load timings in array to export to Python as a tuple. */
    void loadTimers(real t[10]);

    void resetTimers();

    /** Returns the number of resorts done during a single call to integrator.run().
        Its value is reset to zero at the beginning of each run. */
    int getNumResorts() const;

    /** Register this class so it can be used from Python. */
    static void registerPython();

protected:
    bool resortFlag;  //!< true implies need for resort of particles
    int nResorts;
    real maxDist;
    real maxCut;

    void run_(int nsteps);

    real integrate1();
    std::vector<size_t> vsidx;

    void integrate2();

    void calcForces();

    void updateForces();

    void updateForcesBlock();

    void initForcesPlist();

    void initForcesParray();

    espressopp::esutil::WallTimer timeIntegrate;  //!< used for timing

    // variables that keep time information about different phases
    real timeRun;
    real timeLost;
    real timeForce;
    real timeForceComp[100];
    real timeComm1;
    real timeComm2;
    real timeInt1;
    real timeInt2;
    real timeResort;

    static LOG4ESPP_DECL_LOGGER(theLogger);

    //--------------------------------------------------------//
    real timeOtherInitResort;
    real timeOtherInitForcesPlist;
    real timeOtherLoadCells;
    real timeOtherRecalcForces;
    real timeOtherInitForcesParray;
    real timeOtherUnloadCells;
    real timeOtherAftCalcF;

public:
    void loadOtherTimers(real *t);
    //--------------------------------------------------------//
};

}  // namespace integrator
}  // namespace hpx4espp
}  // namespace espressopp

#endif  // HPX4ESPP_INTEGRATOR_VELOCITYVERLET_HPP
