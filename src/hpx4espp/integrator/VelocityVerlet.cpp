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

#include <hpx/config.hpp>
#include <hpx/include/parallel_for_loop.hpp>

#include "hpx4espp/include/logging.hpp"
#include "hpx4espp/include/errors.hpp"
#include "hpx4espp/utils/algorithms/for_loop.hpp"
#include "hpx4espp/utils/algorithms/transform_reduce.hpp"
#include "hpx4espp/utils/multithreading.hpp"
#include "hpx4espp/integrator/VelocityVerlet.hpp"

#include "vec/integrator/VelocityVerlet.hpp"
#include "iterator/CellListIterator.hpp"
#include "interaction/Interaction.hpp"
#include "interaction/Potential.hpp"
#include "System.hpp"
#include "storage/Storage.hpp"
#include "mpi.hpp"
#include "python.hpp"

#include <iomanip>

#ifdef VTRACE
#include "vampirtrace/vt_user.h"
#else
#define VT_TRACER(name)
#endif

namespace espressopp
{
namespace hpx4espp
{
namespace integrator
{
using namespace interaction;
using namespace iterator;
using namespace esutil;

LOG4ESPP_LOGGER(VelocityVerlet::theLogger, "VelocityVerlet");

VelocityVerlet::VelocityVerlet(shared_ptr<SystemHPX> system,
                               shared_ptr<storage::StorageHPX> storageHPX)
    : baseClass(system, storageHPX)
{
    LOG4ESPP_INFO(theLogger, "construct VelocityVerlet");
    resortFlag = true;
    maxDist = 0.0;
    nResorts = 0;
}

void VelocityVerlet::run(int nsteps)
{
    HPX4ESPP_DEBUG_MSG("VelocityVerlet::run");
    hpx4espp::utils::runAsHPXThread([this, nsteps] { this->run_(nsteps); });
}

void VelocityVerlet::run_(int nsteps)
{
    nResorts = 0;
    real time;
    timeIntegrate.reset();
    resetTimers();
    System& system = getSystemRef();
    espressopp::storage::Storage& storage = *system.storage;
    real skinHalf = 0.5 * system.getSkin();

    // signal
    baseClass::runInit();

    // Before start make sure that particles are on the right processor
    if (resortFlag)
    {
        const real time = timeIntegrate.getElapsedTime();

        LOG4ESPP_INFO(theLogger, "resort particles");
        storageHPX->decomposeHPX();
        maxDist = 0.0;
        resortFlag = false;

        timeResort += timeIntegrate.getElapsedTime() - time;
    }
    {
        const real time = timeIntegrate.getElapsedTime();
        storageHPX->loadCells();
        timeOtherLoadCells += timeIntegrate.getElapsedTime() - time;
    }
    bool recalcForces = true;  // TODO: more intelligent
    if (recalcForces)
    {
        LOG4ESPP_INFO(theLogger, "recalc forces before starting main integration loop");

        // signal
        baseClass::recalc1();

        updateForces();

        // signal
        baseClass::recalc2();
    }

    for (int i = 0; i < nsteps; i++)
    {
        {
            const real time = timeIntegrate.getElapsedTime();

            const real maxSqDist = integrate1();

            // collective call to allreduce for dmax
            real maxAllSqDist = 0.0;
            mpi::all_reduce(*system.comm, maxSqDist, maxAllSqDist, boost::mpi::maximum<real>());
            maxDist += std::sqrt(maxAllSqDist);
            // std::cout << " step: " << i << " maxDist: " << maxDist << std::endl;

            timeInt1 += timeIntegrate.getElapsedTime() - time;
        }

        if (maxDist > skinHalf) resortFlag = true;

        if (resortFlag)
        {
            const real time = timeIntegrate.getElapsedTime();

            storageHPX->unloadCells();

            // storage.decompose();
            storageHPX->decomposeHPX();

            storageHPX->loadCells();

            maxDist = 0.0;
            resortFlag = false;
            nResorts++;

            timeResort += timeIntegrate.getElapsedTime() - time;
        }

        // update forces
        {
            updateForces();
        }

        {
            const real time = timeIntegrate.getElapsedTime();
            // second-half integration
            integrate2();
            step++;
            timeInt2 += timeIntegrate.getElapsedTime() - time;
        }
    }
    {
        const real time = timeIntegrate.getElapsedTime();
        storageHPX->unloadCells();
        timeOtherUnloadCells += timeIntegrate.getElapsedTime() - time;
    }

    timeRun = timeIntegrate.getElapsedTime();
    timeLost = timeRun - (timeForceComp[0] + timeForceComp[1] + timeForceComp[2] + timeComm1 +
                          timeComm2 + timeInt1 + timeInt2 + timeResort);
}

real VelocityVerlet::integrate1()
{
    auto& vs = this->storageHPX->virtualStorage;
    auto f_integrate1_vs = [this, &vs](size_t const& ivs) {
        auto& particles = vs[ivs].particles;
        return vec::integrator::VelocityVerlet::integrate1(particles, dt);
    };

    /// workaround to const& requirement for argument to convert operation
    /// see: https://github.com/STEllAR-GROUP/hpx/issues/3651
    if (vsidx.size() != vs.size())
    {
        vsidx.resize(vs.size());
        for (size_t i = 0; i < vsidx.size(); i++) vsidx[i] = i;
    }

    return utils::parallelTransformReduce(
        vsidx.begin(), vsidx.end(), 0.0,
        [](real const& a, real const& b) { return std::max(a, b); }, f_integrate1_vs);
}

void VelocityVerlet::integrate2()
{
    auto& vs = this->storageHPX->virtualStorage;
    auto f_integrate2_vs = [this, &vs](size_t const& ivs) {
        vec::integrator::VelocityVerlet::integrate2(vs[ivs].particles, dt);
    };

    utils::parallelForLoop(size_t(0), vs.size(), f_integrate2_vs);
}

void VelocityVerlet::initForcesParray()
{
    real time = timeIntegrate.getElapsedTime();

    auto& vs = this->storageHPX->virtualStorage;
    auto f_initf_vs = [this, &vs](size_t const& ivs) { vs[ivs].particles.zeroForces(); };

    utils::parallelForLoop(size_t(0), vs.size(), f_initf_vs);

    timeOtherInitForcesParray += timeIntegrate.getElapsedTime() - time;
}

void VelocityVerlet::calcForces()
{
    initForcesParray();
    {
        // TODO: Might need to place interaction list in HPXRuntime
        System& sys = getSystemRef();
        const espressopp::interaction::InteractionList& srIL = sys.shortRangeInteractions;

        for (size_t i = 0; i < srIL.size(); i++)
        {
            LOG4ESPP_INFO(theLogger, "compute forces for srIL " << i << " of " << srIL.size());
            real time;
            time = timeIntegrate.getElapsedTime();
            srIL[i]->addForces();
            timeForceComp[i] += timeIntegrate.getElapsedTime() - time;
        }
        // aftCalcFLocal();
    }
}

void VelocityVerlet::updateForces()
{
    // Implement force update here
    // Initial implementation: blocking update following original

    real time;

    time = timeIntegrate.getElapsedTime();
    storageHPX->updateGhostsBlocking();
    timeComm1 += timeIntegrate.getElapsedTime() - time;

    time = timeIntegrate.getElapsedTime();
    calcForces();
    timeForce += timeIntegrate.getElapsedTime() - time;

    time = timeIntegrate.getElapsedTime();
    storageHPX->collectGhostForcesBlocking();
    timeComm2 += timeIntegrate.getElapsedTime() - time;

    time = timeIntegrate.getElapsedTime();
    baseClass::aftCalcF();
    timeOtherAftCalcF += timeIntegrate.getElapsedTime() - time;
}

void VelocityVerlet::initForcesPlist()
{
    // forces are initialized for real + ghost particles

    System& system = getSystemRef();
    CellList localCells = system.storage->getLocalCells();

    LOG4ESPP_INFO(theLogger, "init forces for real + ghost particles");

    for (CellListIterator cit(localCells); !cit.isDone(); ++cit)
    {
        cit->force() = 0.0;
        cit->drift() = 0.0;  // Can in principle be commented, when drift is not used.
    }
}

void VelocityVerlet::resetTimers()
{
    timeForce = 0.0;
    for (int i = 0; i < 100; i++) timeForceComp[i] = 0.0;
    timeComm1 = 0.0;
    timeComm2 = 0.0;
    timeInt1 = 0.0;
    timeInt2 = 0.0;
    timeResort = 0.0;

    //--------------------------------------------------------//
    timeOtherInitResort = 0.0;
    timeOtherInitForcesPlist = 0.0;
    timeOtherLoadCells = 0.0;
    timeOtherRecalcForces = 0.0;
    timeOtherInitForcesParray = 0.0;
    timeOtherUnloadCells = 0.0;
    timeOtherAftCalcF = 0.0;
    //--------------------------------------------------------//
}

using namespace boost::python;

static object wrapGetTimers(class VelocityVerlet* obj)
{
    real tms[10];
    obj->loadTimers(tms);
    return boost::python::make_tuple(tms[0], tms[1], tms[2], tms[3], tms[4], tms[5], tms[6], tms[7],
                                     tms[8], tms[9]);
}

void VelocityVerlet::loadTimers(real t[10])
{
    t[0] = timeRun;
    t[1] = timeForceComp[0];
    t[2] = timeForceComp[1];
    t[3] = timeForceComp[2];
    t[4] = timeComm1;
    t[5] = timeComm2;
    t[6] = timeInt1;
    t[7] = timeInt2;
    t[8] = timeResort;
    t[9] = timeLost;
}

int VelocityVerlet::getNumResorts() const { return nResorts; }

void VelocityVerlet::loadOtherTimers(real* t)
{
    t[0] = timeOtherInitResort;
    t[1] = timeOtherInitForcesPlist;
    t[2] = timeOtherLoadCells;
    t[3] = timeOtherRecalcForces;
    t[4] = timeOtherInitForcesParray;
    t[5] = timeOtherUnloadCells;
    t[6] = timeOtherAftCalcF;
}

static object wrapGetOtherTimers(class VelocityVerlet* obj)
{
    real tms[7];
    obj->loadOtherTimers(tms);
    return boost::python::make_tuple(tms[0], tms[1], tms[2], tms[3], tms[4], tms[5], tms[6]);
}

/****************************************************
** REGISTRATION WITH PYTHON
****************************************************/

void VelocityVerlet::registerPython()
{
    using namespace espressopp::python;

    // Note: use noncopyable and no_init for abstract classes
    class_<hpx4espp::integrator::VelocityVerlet,
           bases<espressopp::integrator::MDIntegrator, MDIntegratorHPX>, boost::noncopyable>(
        "hpx4espp_integrator_VelocityVerlet",
        init<shared_ptr<SystemHPX>, shared_ptr<storage::StorageHPX> >())
        .def("run", &hpx4espp::integrator::VelocityVerlet::run)
        .def("getTimers", &wrapGetTimers)
        .def("getOtherTimers", &wrapGetOtherTimers)
        .def("resetTimers", &VelocityVerlet::resetTimers)
        .def("getNumResorts", &VelocityVerlet::getNumResorts);
}
}  // namespace integrator
}  // namespace hpx4espp
}  // namespace espressopp
