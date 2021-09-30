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
#include "VelocityVerletHybrid.hpp"
#include <iomanip>
#include "iterator/CellListIterator.hpp"
#include "interaction/Potential.hpp"
#include "storage/Storage.hpp"
#include "bc/BC.hpp"

namespace espressopp
{
using namespace std;
namespace integrator
{
using namespace interaction;
using namespace iterator;
using namespace esutil;

LOG4ESPP_LOGGER(VelocityVerletHybrid::theLogger, "VelocityVerletHybrid");

VelocityVerletHybrid::VelocityVerletHybrid(shared_ptr<System> system,
                                           shared_ptr<FixedVSList> vs_list)
    : MDIntegrator(system), vs_list_(vs_list)
{
    LOG4ESPP_INFO(theLogger, "construct VelocityVerletHybrid");
    resortFlag = true;
    maxDist = 0.0;
    timeIntegrate.reset();
    resetTimers();
}

VelocityVerletHybrid::~VelocityVerletHybrid()
{
    LOG4ESPP_INFO(theLogger, "free VelocityVerletHybrid");
}

real VelocityVerletHybrid::updateVS()
{
    LOG4ESPP_INFO(theLogger, "update position and velocity of VS");

    System &system = getSystemRef();
    storage::Storage &storage = *system.storage;
    const bc::BC &bc = *getSystemRef().bc;

    real maxSqDist = 0.0;

    FixedVSList::GlobalTuples vs = vs_list_->globalTuples;
    FixedVSList::GlobalTuples::iterator it = vs.begin();
    for (; it != vs.end(); ++it)
    {
        Particle *vp = storage.lookupLocalParticle(it->first);  // update local ghost particle.
        if (vp)
        {
            if (vp->id() != it->first) throw std::runtime_error("something is really wrong!");
            Real3D cmp(0.0, 0.0, 0.0);
            // Real3D cmv(0.0, 0.0, 0.0);
            for (auto & itp : it->second)
            {
                Particle *atp =
                    storage.lookupLocalParticle(itp);  // based on real or ghost position.
                if (atp)
                {
                    if (atp->id() != itp) throw std::runtime_error("something is really wrong!");
                    const Real3D &pos = atp->position();
                    real m = atp->mass();
                    LOG4ESPP_DEBUG(theLogger, "vp-" << vp->id() << " atp:" << atp->id() << " "
                                                    << atp->position());
                    Real3D vec12;
                    bc.getMinimumImageVectorBox(vec12, pos, vp->position());
                    cmp[0] += m * vec12[0];
                    cmp[1] += m * vec12[1];
                    cmp[2] += m * vec12[2];
                }
                else
                {
                    std::cout << " AT particle (" << itp << ") of VP " << vp->id() << "-"
                              << vp->ghost() << " not found in tuples ";
                    std::cout << " (" << vp->position() << ")" << std::endl;
                    exit(1);
                }
            }
            cmp /= vp->mass();
            // cmv /= vp->mass();
            LOG4ESPP_DEBUG(theLogger, "vp-" << vp->id() << " cmp=" << cmp);
            Real3D old_p = vp->position();
            vp->position() += cmp;
            // vp->velocity() = cmv;
            //  Check how fare we move new cg position.
            Real3D d_p = old_p - vp->position();
            maxSqDist = std::max(maxSqDist, d_p.sqr());
        }
    }

    real maxAllSqDist;
    mpi::all_reduce(*system.comm, maxSqDist, maxAllSqDist, boost::mpi::maximum<real>());

    LOG4ESPP_INFO(theLogger, " particles in updateVS"
                                 << ", max move local = " << sqrt(maxSqDist)
                                 << ", global = " << sqrt(maxAllSqDist));

    return sqrt(maxAllSqDist);
}

void VelocityVerletHybrid::updateVS_vel()
{
    LOG4ESPP_INFO(theLogger, "update velocity of VS");

    FixedVSList::GlobalTuples vs = vs_list_->globalTuples;
    FixedVSList::GlobalTuples::iterator it = vs.begin();

    System &system = getSystemRef();
    storage::Storage &storage = *system.storage;

    for (; it != vs.end(); ++it)
    {
        Particle *vp = storage.lookupRealParticle(it->first);
        if (vp)
        {
            Real3D cmv(0.0, 0.0, 0.0);
            for (FixedVSList::tuple::iterator itp = it->second.begin(); itp != it->second.end();
                 ++itp)
            {
                Particle *at = storage.lookupLocalParticle(*itp);
                if (at)
                {
                    cmv += at->mass() * at->velocity();
                }
                else
                {
                    std::cout << " AT particle (" << *itp << ") of VP " << vp->id() << "-"
                              << vp->ghost() << " not found in tuples ";
                    std::cout << " (" << vp->position() << ")" << std::endl;
                    exit(1);
                }
            }
            // Updates velocity.
            cmv /= vp->mass();
            vp->velocity() = cmv;
        }
    }
}

void VelocityVerletHybrid::run(int nsteps)
{
    int nResorts = 0;
    real time;
    timeIntegrate.reset();
    System &system = getSystemRef();
    storage::Storage &storage = *system.storage;
    skinHalf = 0.5 * system.getSkin();

    // Prepare the force comp timers if the size is not valid.
    const InteractionList &srIL = system.shortRangeInteractions;
    if (timeForceComp.size() < srIL.size())
    {
        LOG4ESPP_DEBUG(theLogger, "Prepare timeForceComp");
        timeForceComp.clear();
        for (size_t i = 0; i < srIL.size(); i++)
        {
            timeForceComp.push_back(0.0);
        }
    }

    time = timeIntegrate.getElapsedTime();
    // signal
    runInit();
    timeRunInitS += timeIntegrate.getElapsedTime() - time;

    // Before start make sure that particles are on the right processor
    if (resortFlag)
    {
        LOG4ESPP_INFO(theLogger, "resort particles");
        storage.decompose();
        maxDist = 0.0;
        resortFlag = false;
    }

    bool recalcForces = true;  // TODO: more intelligent

    if (recalcForces)
    {
        LOG4ESPP_INFO(theLogger, "recalc forces before starting main integration loop");

        time = timeIntegrate.getElapsedTime();
        // signal
        recalc1();
        timeRecalc1S += timeIntegrate.getElapsedTime() - time;

        updateForces();

        time = timeIntegrate.getElapsedTime();
        // signal
        recalc2();
        timeRecalc2S += timeIntegrate.getElapsedTime() - time;
    }

    LOG4ESPP_INFO(theLogger, "starting main integration loop (nsteps=" << nsteps << ")");

    for (int i = 0; i < nsteps; i++)
    {
        LOG4ESPP_INFO(theLogger, "Next step " << i << " of " << nsteps << " starts");

        time = timeIntegrate.getElapsedTime();
        // signal
        befIntP();
        timeBefIntPS += timeIntegrate.getElapsedTime() - time;

        LOG4ESPP_INFO(theLogger, "updating positions and velocities")
        maxDist += integrate1();
        timeInt1 += timeIntegrate.getElapsedTime() - time;

        time = timeIntegrate.getElapsedTime();
        // signal
        aftIntP();
        timeAftIntPS += timeIntegrate.getElapsedTime() - time;

        LOG4ESPP_INFO(theLogger, "maxDist = " << maxDist << ", skin/2 = " << skinHalf);

        if (maxDist > skinHalf) resortFlag = true;

        if (resortFlag)
        {
            time = timeIntegrate.getElapsedTime();
            LOG4ESPP_INFO(theLogger, "step " << i << ": resort particles");
            storage.decompose();
            maxDist = 0.0;
            resortFlag = false;
            nResorts++;
            timeResort += timeIntegrate.getElapsedTime() - time;
        }

        LOG4ESPP_INFO(theLogger, "updating forces")
        updateForces();

        timeIntegrate.startMeasure();
        // signal
        befIntV();
        timeBefIntVS += timeIntegrate.stopMeasure();

        time = timeIntegrate.getElapsedTime();
        integrate2();
        timeInt2 += timeIntegrate.getElapsedTime() - time;

        timeIntegrate.startMeasure();
        // signal
        aftIntV();
        timeAftIntVS += timeIntegrate.stopMeasure();
    }

    timeRun = timeIntegrate.getElapsedTime();

    LOG4ESPP_INFO(theLogger, "finished run");
}

void VelocityVerletHybrid::resetTimers()
{
    timeForce = 0.0;

    timeComm1 = 0.0;
    timeComm2 = 0.0;
    timeInt1 = 0.0;
    timeInt2 = 0.0;
    timeResort = 0.0;

    // Reset signal timers.
    timeRunInitS = 0.0;
    timeRecalc1S = 0.0;
    timeRecalc2S = 0.0;
    timeBefIntPS = 0.0;
    timeAftIntPS = 0.0;
    timeAftInitFS = 0.0;
    timeAftCalcFS = 0.0;
    timeBefIntVS = 0.0;
    timeAftIntVS = 0.0;

    timeVS = 0.0;
    timeVSvel = 0.0;
    timeVSdistrF = 0.0;
}

void VelocityVerletHybrid::loadTimers(std::vector<real> &return_vector,
                                      std::vector<std::string> &return_labels)
{
    return_vector.push_back(timeRun);
    return_labels.push_back("timeRun");
    for (size_t i = 0; i < timeForceComp.size(); i++)
    {
        return_vector.push_back(timeForceComp[i]);
        std::stringstream ss;
        ss << "f" << i;
        return_labels.push_back(ss.str());
    }

    // signal timers.
    return_vector.push_back(timeRunInitS);
    return_vector.push_back(timeRecalc1S);
    return_vector.push_back(timeRecalc2S);
    return_vector.push_back(timeBefIntPS);
    return_vector.push_back(timeAftIntPS);
    return_vector.push_back(timeAftCalcFS);
    return_vector.push_back(timeBefIntVS);
    return_vector.push_back(timeAftIntVS);

    return_vector.push_back(timeVS);
    return_vector.push_back(timeVSdistrF);
    return_vector.push_back(timeVSvel);
    return_vector.push_back(timeComm1);
    return_vector.push_back(timeComm2);
    return_vector.push_back(timeInt1);
    return_vector.push_back(timeInt2);
    return_vector.push_back(timeResort);

    return_labels.push_back("timeRunInitS");
    return_labels.push_back("timeRecalc1S");
    return_labels.push_back("timeRecalc2S");
    return_labels.push_back("timeBefIntPS");
    return_labels.push_back("timeAftIntPS");
    return_labels.push_back("timeAftCalcFS");
    return_labels.push_back("timeBefIntVS");
    return_labels.push_back("timeAftIntVS");
    return_labels.push_back("timeVS");
    return_labels.push_back("timeVSdistrF");
    return_labels.push_back("timeVSvel");
    return_labels.push_back("timeComm1");
    return_labels.push_back("timeComm2");
    return_labels.push_back("timeInt1");
    return_labels.push_back("timeInt2");
    return_labels.push_back("timeResort");
}

real VelocityVerletHybrid::integrate1()
{
    System &system = getSystemRef();
    CellList realCells = system.storage->getRealCells();

    // loop over all particles of the local cells
    int count = 0;
    real maxSqDist = 0.0;  // maximal square distance a particle moves
    for (CellListIterator cit(realCells); !cit.isDone(); ++cit)
    {
        if (cit->vp())  // propagate only real particles, skip virtual sites.
            continue;
        real sqDist = 0.0;
        LOG4ESPP_INFO(theLogger,
                      "updating first half step of velocities and full step of positions")
        LOG4ESPP_DEBUG(theLogger, "Particle " << cit->id() << ", pos = " << cit->position()
                                              << ", v = " << cit->velocity()
                                              << ", f = " << cit->force());

        real dtfm = 0.5 * dt / cit->mass();

        // Propagate velocities: v(t+0.5*dt) = v(t) + 0.5*dt * f(t)
        cit->velocity() += dtfm * cit->force();

        // Propagate positions (only NVT): p(t + dt) = p(t) + dt * v(t+0.5*dt)
        Real3D deltaP = cit->velocity();

        deltaP *= dt;
        cit->position() += deltaP;
        sqDist += deltaP * deltaP;

        count++;

        maxSqDist = std::max(maxSqDist, sqDist);
    }

    // signal
    inIntP(maxSqDist);

    real maxAllSqDist;
    mpi::all_reduce(*system.comm, maxSqDist, maxAllSqDist, boost::mpi::maximum<real>());

    LOG4ESPP_INFO(theLogger, "moved " << count << " particles in integrate1"
                                      << ", max move local = " << sqrt(maxSqDist)
                                      << ", global = " << sqrt(maxAllSqDist));

    return sqrt(maxAllSqDist);
}

void VelocityVerletHybrid::integrate2()
{
    LOG4ESPP_INFO(theLogger, "updating second half step of velocities")
    System &system = getSystemRef();
    CellList realCells = system.storage->getRealCells();

    // loop over all particles of the local cells
    real half_dt = 0.5 * dt;
    for (CellListIterator cit(realCells); !cit.isDone(); ++cit)
    {
        if (cit->vp()) continue;
        real dtfm = half_dt / cit->mass();
        /* Propagate velocities: v(t+0.5*dt) = v(t) + 0.5*dt * f(t) */
        cit->velocity() += dtfm * cit->force();
    }

    // Update velocity of VS based on the AT velocities.
    timeIntegrate.startMeasure();
    // updateVS_vel();
    timeVSvel += timeIntegrate.stopMeasure();

    step++;
}

void VelocityVerletHybrid::calcForces()
{
    LOG4ESPP_INFO(theLogger, "calculate forces");

    initForces();

    timeIntegrate.startMeasure();
    // signal
    aftInitF();
    timeAftInitFS += timeIntegrate.stopMeasure();

    System &sys = getSystemRef();
    const InteractionList &srIL = sys.shortRangeInteractions;
    real time;
    for (size_t i = 0; i < srIL.size(); i++)
    {
        LOG4ESPP_INFO(theLogger, "compute forces for srIL " << i << " of " << srIL.size());
        time = timeIntegrate.getElapsedTime();
        srIL[i]->addForces();
        timeForceComp[i] += timeIntegrate.getElapsedTime() - time;
    }
}

void VelocityVerletHybrid::updateForces()
{
    LOG4ESPP_INFO(theLogger, "update ghosts, calculate forces and collect ghost forces")
    real time;
    storage::Storage &storage = *getSystemRef().storage;

    // Make sure that positions and velocity of VS sites is correct with respect to the atoms.
    timeIntegrate.startMeasure();
    real maxDist = updateVS();
    timeVS += timeIntegrate.stopMeasure();

    if (maxDist > skinHalf)
    {
        LOG4ESPP_TRACE(theLogger,
                       "maxDist=" << maxDist << " >" << skinHalf << " storage.decompose()");
        time = timeIntegrate.getElapsedTime();
        storage.decompose();
        timeComm1 += timeIntegrate.getElapsedTime() - time;
    }
    else
    {
        time = timeIntegrate.getElapsedTime();
        storage.updateGhosts();
        timeComm1 += timeIntegrate.getElapsedTime() - time;
    }

    time = timeIntegrate.getElapsedTime();
    calcForces();
    timeForce += timeIntegrate.getElapsedTime() - time;

    time = timeIntegrate.getElapsedTime();
    storage.collectGhostForces();
    timeComm2 += timeIntegrate.getElapsedTime() - time;

    // Distribute forces from VS to AT.
    timeIntegrate.startMeasure();
    distributeVSforces();
    timeVSdistrF += timeIntegrate.stopMeasure();

    timeIntegrate.startMeasure();
    // signal
    aftCalcF();
    timeAftCalcFS += timeIntegrate.stopMeasure();
}

void VelocityVerletHybrid::initForces()
{
    // forces are initialized for real + ghost particles

    System &system = getSystemRef();
    CellList localCells = system.storage->getLocalCells();

    LOG4ESPP_INFO(theLogger, "init forces for real + ghost particles");

    for (CellListIterator cit(localCells); !cit.isDone(); ++cit)
    {
        cit->force() = 0.0;
        cit->drift() = 0.0;  // Can in principle be commented, when drift is not used.
    }
}

void VelocityVerletHybrid::distributeVSforces()
{
    // Zeros forces on ghost particles.
    System &system = getSystemRef();
    storage::Storage &storage = *system.storage;
    CellList ghostCells = storage.getGhostCells();
    LOG4ESPP_INFO(theLogger, "zeros forces on ghost particles");

    for (CellListIterator cit(ghostCells); !cit.isDone(); ++cit)
    {
        cit->force() = 0.0;
        cit->drift() = 0.0;
    }

    // Distribute forces to AT particles from VS.
    FixedVSList::GlobalTuples vs = vs_list_->globalTuples;
    FixedVSList::GlobalTuples::iterator it = vs.begin();

    for (; it != vs.end(); ++it)
    {
        Particle *vp = storage.lookupRealParticle(it->first);
        if (vp)
        {
            for (FixedVSList::tuple::iterator itp = it->second.begin(); itp != it->second.end();
                 ++itp)
            {
                Particle *at = storage.lookupLocalParticle(*itp);
                if (at)
                {
                    at->force() += (at->mass() / vp->mass()) * vp->force();
                }
                else
                {
                    std::cout << " AT particle (" << *itp << ") of VP " << vp->id() << "-"
                              << vp->ghost() << " not found in tuples ";
                    std::cout << " (" << vp->position() << ")" << std::endl;
                    exit(1);
                }
            }
        }
    }
    real time = timeIntegrate.getElapsedTime();
    storage.collectGhostForces();
    timeComm2 += timeIntegrate.getElapsedTime() - time;
}

static boost::python::object wrapGetTimers(class VelocityVerletHybrid *obj)
{
    std::vector<real> timers;
    std::vector<std::string> labels;
    obj->loadTimers(timers, labels);

    boost::python::list return_list;
    for (size_t i = 0; i < timers.size(); i++)
    {
        return_list.append(boost::python::make_tuple(labels[i], timers[i]));
    }
    return return_list;
}

/****************************************************
** REGISTRATION WITH PYTHON
****************************************************/

void VelocityVerletHybrid::registerPython()
{
    using namespace espressopp::python;

    // Note: use noncopyable and no_init for abstract classes
    class_<VelocityVerletHybrid, bases<MDIntegrator>, boost::noncopyable>(
        "integrator_VelocityVerletHybrid", init<shared_ptr<System>, shared_ptr<FixedVSList> >())
        .def("getTimers", &wrapGetTimers)
        .def("resetTimers", &VelocityVerletHybrid::resetTimers);
}
}  // namespace integrator
}  // namespace espressopp
