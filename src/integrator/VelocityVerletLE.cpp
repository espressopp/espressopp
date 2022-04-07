/*
  Copyright (C) 2017
      Jakub Krajniak (jkrajniak at gmail.com)
  Copyright (C) 2012,2013
      Max Planck Institute for Polymer Research
  Copyright (C) 2008,2009,2010,2011
      Max-Planck-Institute for Polymer Research & Fraunhofer SCAI
  Copyright (C) 2022
      Data Center, Johannes Gutenberg University Mainz

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

//#include <iomanip>
#include "python.hpp"
#include "VelocityVerletLE.hpp"
#include <iomanip>
#include "iterator/CellListIterator.hpp"
#include "interaction/Interaction.hpp"
#include "interaction/Potential.hpp"
#include "System.hpp"
#include "iostream"
#include "bc/BC.hpp"
#include "storage/Storage.hpp"
#include "mpi.hpp"
//#include <cstdlib>

#ifdef VTRACE
#include "vampirtrace/vt_user.h"
#else
#define VT_TRACER(name)
#endif

namespace espressopp
{
int shift_count = 0;

using namespace std;
namespace integrator
{
using namespace interaction;
using namespace iterator;
using namespace esutil;

LOG4ESPP_LOGGER(VelocityVerletLE::theLogger, "VelocityVerletLE");

VelocityVerletLE::VelocityVerletLE(shared_ptr<System> system, real _shearRate)
    : MDIntegrator(system), shearRate(_shearRate)
{
    LOG4ESPP_INFO(theLogger, "construct VelocityVerletLE");
    resortFlag = true;
    maxDist = 0.0;
    nResorts = 0;
}

VelocityVerletLE::~VelocityVerletLE() { LOG4ESPP_INFO(theLogger, "free VelocityVerletLE"); }

void VelocityVerletLE::run(int nsteps)
{
    VT_TRACER("run");
    nResorts = 0;
    real time;
    timeIntegrate.reset();
    resetTimers();
    System& system = getSystemRef();
    storage::Storage& storage = *system.storage;
    real skinHalf = 0.5 * system.getSkin();
    // storage::Storage& storage2 = *getSystemRef().storage;
    real Lx = system.bc->getBoxL()[0];
    real Lz = system.bc->getBoxL()[2];
    int ngrid =
        storage.getInt3DCellGrid()[0] * system.NGridSize[0];  // storage.getInt3DNodeGrid()[2];
    if (shearRate != .0)
    {
        system.shearRate = shearRate;
        system.ifShear = true;
	}
	else
	{
		throw std::runtime_error("VelocityVerletLE error: LE integrator is called \
		but the shear rate is set to zero \n");
	}

    // signal
    runInit();

    // Before start make sure that particles are on the right processor
    if (resortFlag)
    {
        VT_TRACER("resort");
        // time = timeIntegrate.getElapsedTime();
        LOG4ESPP_INFO(theLogger, "resort particles");
        storage.decompose();
        maxDist = 0.0;
        resortFlag = false;
        // timeResort += timeIntegrate.getElapsedTime();
    }

    bool recalcForces = true;  // TODO: more intelligent

    if (recalcForces)
    {
        LOG4ESPP_INFO(theLogger, "recalc forces before starting main integration loop");

        // signal
        recalc1();

        updateForces();
        if (LOG4ESPP_DEBUG_ON(theLogger))
        {
            // printForces(false);   // forces are reduced to real particles
        }

        // signal
        recalc2();
    }

    LOG4ESPP_INFO(theLogger, "starting main integration loop (nsteps=" << nsteps << ")");

    // if (rename("FLAG_P","FLAG_P")==0 && getenv("IRANK")!=NULL)
    // system.irank=atoi(getenv("IRANK"));
    if (getenv("LMODE") != NULL) system.lebcMode = atoi(getenv("LMODE"));

    for (int i = 0; i < nsteps; i++)
    {
        LOG4ESPP_INFO(theLogger, "Next step " << i << " of " << nsteps << " starts");

        // saveOldPos(); // save particle positions needed for constraints

        // signal
        befIntP();

        time = timeIntegrate.getElapsedTime();
        LOG4ESPP_INFO(theLogger, "updating positions and velocities")
        // if (rename("FLAG_P","FLAG_P")==0 && system.comm->rank()==system.irank){
        // std::cout<<" INT01> "<<" \n";}
        maxDist += integrate1();
        timeInt1 += timeIntegrate.getElapsedTime() - time;

        /*
        real cellsize = 1.4411685442;
        if (maxDist > 1.4411685442){
          cout<<"WARNING!!!!!! huge jump: "<<maxDist<<endl;
          exit(1);
        }*/

        // signal
        aftIntP();
        // if (rename("FLAG_P","FLAG_P")==0 && system.comm->rank()==system.irank){
        // std::cout<<" aftIntP> ("<<system.comm->rank()<<") \n";}

        // remap neighbour cells

        LOG4ESPP_INFO(theLogger, "maxDist = " << maxDist << ", skin/2 = " << skinHalf);

        int ctmp = static_cast<int>(floor(
            shearRate * static_cast<real>(getStep()) * getTimeStep() * ngrid * Lz / Lx + 0.5));
        int cshift = static_cast<int>(floor(
            shearRate * static_cast<real>(getStep() + 1) * getTimeStep() * ngrid * Lz / Lx + 0.5));

        if (cshift > ctmp)
        {
            shift_count++;
            storage.remapNeighbourCells(shift_count);
            system.ghostShift = shift_count;
            resortFlag = true;
        }
        else if (cshift < ctmp)
        {
            shift_count++;
            storage.remapNeighbourCells(-shift_count);
            system.ghostShift = -shift_count;
            resortFlag = true;
        }

        // resortFlag = true;
        if (ctmp > shift_count)
        {
            std::cout << " ERR> SHIFT: " << ctmp << " / " << shift_count << " \n";
            throw std::runtime_error("VelocityVerletLE error: numeric error leading to no cell shifts \n");
        }

        if (maxDist > skinHalf) resortFlag = true;

        if (resortFlag)
        {
            VT_TRACER("resort1");
            time = timeIntegrate.getElapsedTime();
            LOG4ESPP_INFO(theLogger, "step " << i << ": resort particles");

            storage.decompose();

            maxDist = 0.0;
            resortFlag = false;
            nResorts++;
            timeResort += timeIntegrate.getElapsedTime() - time;
        }

        // Analysis to get stress tensors
        if (system.ifViscosity)
        {
            system.dyadicP_xz = .0;
            system.dyadicP_zx = .0;
        }

        LOG4ESPP_INFO(theLogger, "updating forces")
        updateForces();

        // if (rename("FLAG_P","FLAG_P")==0 && system.comm->rank()==system.irank){
        // std::cout<<" UPDFC> "<<" \n";}
        // signal
        befIntV();

        time = timeIntegrate.getElapsedTime();
        /*if (i==nsteps-1 && rename("FLAG_V","FLAG_V")==0){
          integrate2_extra();
        }else*/
        integrate2();
        timeInt2 += timeIntegrate.getElapsedTime() - time;
        // if (rename("FLAG_P","FLAG_P")==0 && system.comm->rank()==system.irank){
        // std::cout<<" INT02> "<<" \n";}

        // signal
        aftIntV();
        // if (rename("FLAG_P","FLAG_P")==0 && system.comm->rank()==system.irank){
        // std::cout<<" aftIntV> "<<" \n";}
    }

    timeRun = timeIntegrate.getElapsedTime();
    timeLost = timeRun - (timeForceComp[0] + timeForceComp[1] + timeForceComp[2] + timeComm1 +
                          timeComm2 + timeInt1 + timeInt2 + timeResort);

    LOG4ESPP_INFO(theLogger, "finished run");
}

void VelocityVerletLE::resetTimers()
{
    timeForce = 0.0;
    for (int i = 0; i < 100; i++) timeForceComp[i] = 0.0;
    timeComm1 = 0.0;
    timeComm2 = 0.0;
    timeInt1 = 0.0;
    timeInt2 = 0.0;
    timeResort = 0.0;
}

using namespace boost::python;

static object wrapGetTimers(class VelocityVerletLE* obj)
{
    real tms[10];
    obj->loadTimers(tms);
    return boost::python::make_tuple(tms[0], tms[1], tms[2], tms[3], tms[4], tms[5], tms[6], tms[7],
                                     tms[8], tms[9]);
}

void VelocityVerletLE::loadTimers(real t[10])
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

void VelocityVerletLE::printTimers()
{
    using namespace std;
    real pct;

    cout << endl;
    cout << "run = " << setiosflags(ios::fixed) << setprecision(3) << timeRun << endl;
    pct = 100.0 * (timeForceComp[0] / timeRun);
    cout << "pair (%) = " << timeForceComp[0] << " (" << pct << ")" << endl;
    pct = 100.0 * (timeForceComp[1] / timeRun);
    cout << "FENE (%) = " << timeForceComp[1] << " (" << pct << ")" << endl;
    pct = 100.0 * (timeForceComp[2] / timeRun);
    cout << "angle (%) = " << timeForceComp[2] << " (" << pct << ")" << endl;
    pct = 100.0 * (timeComm1 / timeRun);
    cout << "comm1 (%) = " << timeComm1 << " (" << pct << ")" << endl;
    pct = 100.0 * (timeComm2 / timeRun);
    cout << "comm2 (%) = " << timeComm2 << " (" << pct << ")" << endl;
    pct = 100.0 * (timeInt1 / timeRun);
    cout << "int1 (%) = " << timeInt1 << " (" << pct << ")" << endl;
    pct = 100.0 * (timeInt2 / timeRun);
    cout << "int2 (%) = " << timeInt2 << " (" << pct << ")" << endl;
    pct = 100.0 * (timeResort / timeRun);
    cout << "resort (%) = " << timeResort << " (" << pct << ")" << endl;
    pct = 100.0 * (timeLost / timeRun);
    cout << "other (%) = " << timeLost << " (" << pct << ")" << endl;
    cout << endl;
}

int VelocityVerletLE::getNumResorts() const { return nResorts; }

real VelocityVerletLE::integrate1()
{
    System& system = getSystemRef();
    CellList realCells = system.storage->getRealCells();
    // if (rename("FLAG_P","FLAG_P")==0 && getenv("VAR1")!=NULL &&
    // system.comm->rank()==system.irank){ std:cout<<"===========\nRUN01> TimeStep= "<<getStep()<<" /
    // "<<system.shearOffset<<"\n------\n";}

    // loop over all particles of the local cells
    int count = 0;
    bool flag_inBoxShear = false;
    real maxSqDist = 0.0;  // maximal square distance a particle moves
    real Lx = system.bc->getBoxL()[0];
    real Lz = system.bc->getBoxL()[2];
    real halfL = Lz / 2.0;
    real cutoff = system.maxCutoff;

    for (CellListIterator cit(realCells); !cit.isDone(); ++cit)
    {
        real sqDist = 0.0;
        LOG4ESPP_INFO(theLogger,
                      "updating first half step of velocities and full step of positions")
        LOG4ESPP_DEBUG(theLogger, "Particle " << cit->id() << ", pos = " << cit->position()
                                              << ", v = " << cit->velocity()
                                              << ", f = " << cit->force());

        /* more precise for DEBUG:
        printf("Particle %d, pos = %16.12f %16.12f %16.12f, v = %16.12f, %16.12f %16.12f, f =
        %16.12f %16.12f %16.12f\n", cit->p.id, cit->r.p[0], cit->r.p[1], cit->r.p[2], cit->m.v[0],
        cit->m.v[1], cit->m.v[2], cit->f.f[0], cit->f.f[1], cit->f.f[2]);
        */

        real dtfm = 0.5 * dt / cit->mass();

        // Propagate velocities for X dim.
        cit->velocity()[0] += dtfm * cit->force()[0];
        real vshear = shearRate * (cit->position()[2] + 0.5 * cit->velocity()[2] * dt +
                                   dtfm * cit->force()[2] * dt / 3.0 - halfL);
        // Propagate velocities for Y-Z dim.: v(t+0.5*dt) = v(t) + 0.5*dt * f(t)
        cit->velocity()[2] += dtfm * cit->force()[2];
        cit->velocity()[1] += dtfm * cit->force()[1];

        // Propagate positions (only NVT): p(t + dt) = p(t) + dt * v(t+0.5*dt)
        Real3D deltaP = {.0, .0, .0};
        deltaP = cit->velocity();

        // Add shear speed into X dim.
        deltaP[0] += vshear;
        deltaP *= dt;
        cit->position() += deltaP;
        sqDist += deltaP * deltaP;
        // Map to ghost postion if z-position is within the top/bottom layer
        // real le_skin=system.getSkin();

        // if (rename("FLAG_P","FLAG_P")==0 && getenv("VAR1")!=NULL)
        // if (cit->id()==atoi(getenv("VAR1")))
        // std::cout<<"INT1-"<<system.comm->rank()<<"> "<<getStep()<<" "<<cit->id()<<"
        // ["<<cit->position()<<"] ("
        //<<deltaP<<") "<<static_cast<int>(floor(cit->position()[2]/Lz))<<" \n";
        int itmp = static_cast<int>(floor(cit->position()[2] / Lz));
        if (abs(itmp) > 1)
        {
            std::cout << "ERR> " << cit->id() << " [" << cit->position()[0] << ","
                      << cit->position()[1] << "," << cit->position()[2] << "] (" << deltaP[0]
                      << "," << deltaP[1] << "," << deltaP[2] << ") \n"
                      << cit->velocity()[2] << " | " << cit->force()[2] << " \n";
            throw std::runtime_error("VelocityVerletLE error: particle crosses over two images");
        }

        count++;

        maxSqDist = std::max(maxSqDist, sqDist);
    }

    // set boundary offset of a shear flow
    real offs;
    offs = shearRate * Lz * (getStep() + 1.0) * getTimeStep();
    int xtmp = static_cast<int>(floor(offs / Lx));
    system.shearOffset = offs - (xtmp + .0) * Lx;
    // if (rename("FLAG_P","FLAG_P")==0 && getenv("VAR1")!=NULL && system.comm->rank()==0)
    // std::cout<<"SHEAR> "<<system.shearOffset<<" \n";

    // signal
    inIntP(maxSqDist);

    real maxAllSqDist;
    mpi::all_reduce(*system.comm, maxSqDist, maxAllSqDist, boost::mpi::maximum<real>());

    LOG4ESPP_INFO(theLogger, "moved " << count << " particles in integrate1"
                                      << ", max move local = " << sqrt(maxSqDist)
                                      << ", global = " << sqrt(maxAllSqDist));

    return sqrt(maxAllSqDist);
}

void VelocityVerletLE::integrate2()
{
    LOG4ESPP_INFO(theLogger, "updating second half step of velocities")
    System& system = getSystemRef();
    CellList realCells = system.storage->getRealCells();

    // loop over all particles of the local cells
    real half_dt = 0.5 * dt;

    if (system.ifViscosity)
    {
        real mv2 = .0;
        for (CellListIterator cit(realCells); !cit.isDone(); ++cit)
        {
            real dtfm = half_dt / cit->mass();
            /* Propagate velocities: v(t+0.5*dt) = v(t) + 0.5*dt * f(t) */
            cit->velocity() += dtfm * cit->force();
            // Need to add propagation of shear speed if necessary
            // Collect xz-&zx- components from stress Tensor
            mv2 += cit->mass() * cit->velocity()[0] * cit->velocity()[2];
        }
    
        real allDyadicP_xz=.0;
        real P_xz=system.dyadicP_xz;
        //mpi::all_reduce(*system.comm, P_xz, allDyadicP_xz, boost::mpi::maximum<real>());
        boost::mpi::all_reduce(*system.comm, P_xz, allDyadicP_xz, std::plus<real>());
        // print the off-diagonal (XZ) component of stress Tensor
        real vol = system.bc->getBoxL()[2] * system.bc->getBoxL()[1] * system.bc->getBoxL()[0];
        if (system.comm->rank()==0) std::cout << "SIGXZ> " << getStep() + 1 << " "
                  << -1.0 / vol * (mv2 + allDyadicP_xz) / shearRate
                  //<<" "<<-1.0/vol*(mv2)/shearRate
                  << " \n";
        system.dyadicP_xz = .0;
        // std::cout<<"SIGZX> "<<getStep()+1<<" "<<-1.0/vol*(mv2+system.dyadicP_zx)/shearRate<<" \n";
        //system.dyadicP_zx = .0;
    }
    else
    {
        for (CellListIterator cit(realCells); !cit.isDone(); ++cit)
        {
            real dtfm = half_dt / cit->mass();
            /* Propagate velocities: v(t+0.5*dt) = v(t) + 0.5*dt * f(t) */
            cit->velocity() += dtfm * cit->force();
        }
    }

    step++;
}

void VelocityVerletLE::integrate2_extra()
{
    LOG4ESPP_INFO(theLogger, "updating second half step of velocities")
    System& system = getSystemRef();
    CellList realCells = system.storage->getRealCells();

    // loop over all particles of the local cells
    real half_dt = 0.5 * dt;

    int bin_size = 50, ntot = 0, *cnt;
    real mo[4] = {.0, .0, .0, .0};
    real *vx, *vz, *sx, *ty, *tz;
    Real3D* ori;
    int nchain = 10, clength = 20;

    if (getenv("NCHAIN") != NULL)
        if (atoi(getenv("NCHAIN")) > 0) nchain = atoi(getenv("NCHAIN"));
    if (getenv("CLENGTH") != NULL)
        if (atoi(getenv("CLENGTH")) > 0) clength = atoi(getenv("CLENGTH"));

    real mv2 = .0;

    cnt = new int[bin_size];
    vx = new real[bin_size];
    vz = new real[bin_size];
    sx = new real[bin_size];
    ty = new real[bin_size];
    tz = new real[bin_size];
    ori = new Real3D[nchain];

    // Init
    for (int i = 0; i < bin_size; i++)
    {
        cnt[i] = 0;
        vx[i] = .0;
        vz[i] = .0;
        sx[i] = .0;
        ty[i] = .0;
        tz[i] = .0;
    }
    real Lz = system.bc->getBoxL()[2];
    for (int i = 0; i < nchain; i++)
    {
        ori[i] = {.0, .0, .0};
    }

    for (CellListIterator cit(realCells); !cit.isDone(); ++cit)
    {
        real dtfm = half_dt / cit->mass();
        /* Propagate velocities: v(t+0.5*dt) = v(t) + 0.5*dt * f(t) */
        cit->velocity() += dtfm * cit->force();

        // Need to add propagation of shear speed if necessary
        real spdShear = shearRate * (cit->position()[2] - Lz / 2.0);
        mo[0] += cit->velocity()[0];
        mo[1] += cit->velocity()[1];
        mo[2] += cit->velocity()[2];
        mo[3] += cit->velocity()[0] + spdShear;

        int itmp = static_cast<int>(floor(cit->position()[2] / Lz * (bin_size + .0)));
        if (itmp < 0)
            itmp = 0;
        else if (itmp >= bin_size)
            itmp = bin_size - 1;

        ntot++;
        cnt[itmp]++;
        vx[itmp] += cit->velocity()[0];
        vz[itmp] += cit->velocity()[2];
        sx[itmp] += cit->velocity()[0] + spdShear;

        // Real3D vtmp={spdShear,.0,.0};
        // vtmp=vtmp+cit->velocity();
        ty[itmp] += (cit->velocity()[1] * cit->velocity()[1]) * cit->mass() / 1.0;
        tz[itmp] += (cit->velocity()[2] * cit->velocity()[2]) * cit->mass() / 1.0;

        // Collect orientation info. (always first N chains if not specific)
        int mono_idx = (cit->id() - 1) % clength;
        int chain_idx = (cit->id() - 1 - mono_idx) / clength;
        if (chain_idx < nchain)
        {
            if (mono_idx == 0)
            {
                ori[chain_idx] -= cit->position();
            }
            else if (mono_idx == clength - 1)
            {
                ori[chain_idx] += cit->position();
            }
        }

        // Collect xz-&zx- components from stress Tensor
        mv2 += cit->mass() * cit->velocity()[0] * cit->velocity()[2];
    }

    // Normalization
    for (int i = 0; i < 4; i++) mo[i] = mo[i] / (ntot + .0);
    for (int i = 0; i < bin_size; i++)
    {
        vx[i] = vx[i] / (cnt[i] + .0);
        vz[i] = vz[i] / (cnt[i] + .0);
        sx[i] = sx[i] / (cnt[i] + .0);
        ty[i] = ty[i] / (cnt[i] + .0);
        tz[i] = tz[i] / (cnt[i] + .0);
    }
    for (int i = 0; i < nchain; i++)
    {
        real val_tmp = ori[i].abs();
        ori[i] = ori[i] / val_tmp;
    }

    // print out
    if (rename("FLAG_VEL", "FLAG_VEL") == 0)
    {
        std::cout << "MOMENTUM> " << getStep() + 1 << " " << mo[0] << " " << mo[1] << " " << mo[2]
                  << " " << mo[3] << " \n";
        // std::cout<<"DISTR> "<<getStep()+1<<" 0.0 0.0 \n";
        for (int i = 0; i < bin_size; i++)
        {
            real z = (i + 0.5) / (bin_size + .0);
            std::cout << "DISTR> " << getStep() + 1 << " " << z << " "
                      << (cnt[i] + .0) / (ntot + .0) * (bin_size + .0) << " \n";
            std::cout << "VX> " << getStep() + 1 << " " << z << " " << vx[i] << " \n";
            std::cout << "VZ> " << getStep() + 1 << " " << z << " " << vz[i] << " \n";
            std::cout << "VSHEAR> " << getStep() + 1 << " " << z << " " << sx[i] / (shearRate)
                      << " \n";
            std::cout << "TP> " << getStep() + 1 << " " << z << " " << ty[i] << " \n";
            std::cout << "TZ> " << getStep() + 1 << " " << z << " " << tz[i] << " \n";
        }
    }

    if (rename("FLAG_OAF", "FLAG_OAF") == 0)
    {
        for (int i = 0; i < nchain; i++)
        {
            std::cout << "CHAIN-" << i + 1 << "> " << getStep() + 1 << " " << ori[i] << " \n";
        }
    }

    if (rename("FLAG_VIS", "FLAG_VIS") == 0)
    {
        real vol = system.bc->getBoxL()[2] * system.bc->getBoxL()[1] * system.bc->getBoxL()[0];
        std::cout << "SIGXZ> " << getStep() + 1 << " "
                  << -1.0 / vol * (mv2 + system.dyadicP_xz) / shearRate
                  //<<" "<<-1.0/vol*(mv2)/shearRate
                  << " \n";
        // std::cout<<"SIGZX> "<<getStep()+1<<" "<<-1.0/vol*(mv2+system.dyadicP_zx)/shearRate<<"
        // \n";
    }
    system.dyadicP_xz = .0;
    system.dyadicP_zx = .0;

    delete cnt;
    delete vx;
    delete vz;
    delete sx;
    delete ty;
    delete tz;
    delete ori;

    step++;
}

void VelocityVerletLE::calcForces()
{
    VT_TRACER("forces");

    LOG4ESPP_INFO(theLogger, "calculate forces");

    initForces();

    // signal
    aftInitF();

    System& sys = getSystemRef();
    const InteractionList& srIL = sys.shortRangeInteractions;

    for (size_t i = 0; i < srIL.size(); i++)
    {
        LOG4ESPP_INFO(theLogger, "compute forces for srIL " << i << " of " << srIL.size());
        real time;
        time = timeIntegrate.getElapsedTime();
        srIL[i]->addForces();
        timeForceComp[i] += timeIntegrate.getElapsedTime() - time;
    }
}

void VelocityVerletLE::updateForces()
{
    LOG4ESPP_INFO(theLogger, "update ghosts, calculate forces and collect ghost forces")
    real time;
    storage::Storage& storage = *getSystemRef().storage;
    time = timeIntegrate.getElapsedTime();

    {
        VT_TRACER("commF");
        storage.updateGhosts();
    }
    timeComm1 += timeIntegrate.getElapsedTime() - time;
    time = timeIntegrate.getElapsedTime();
    calcForces();
    timeForce += timeIntegrate.getElapsedTime() - time;
    time = timeIntegrate.getElapsedTime();
    {
        VT_TRACER("commR");
        storage.collectGhostForces();
    }
    timeComm2 += timeIntegrate.getElapsedTime() - time;

    // signal
    aftCalcF();
}

void VelocityVerletLE::initForces()
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

void VelocityVerletLE::printForces(bool withGhosts)
{
    // print forces of real + ghost particles

    System& system = getSystemRef();
    CellList cells;

    if (withGhosts)
    {
        cells = system.storage->getLocalCells();
        LOG4ESPP_DEBUG(theLogger, "local forces");
    }
    else
    {
        cells = system.storage->getRealCells();
        LOG4ESPP_DEBUG(theLogger, "real forces");
    }

    for (CellListIterator cit(cells); !cit.isDone(); ++cit)
    {
        LOG4ESPP_DEBUG(theLogger, "Particle " << cit->id() << ", force = " << cit->force());
    }
}

void VelocityVerletLE::printPositions(bool withGhosts)
{
    // print positions of real + ghost particles

    System& system = getSystemRef();
    CellList cells;

    if (withGhosts)
    {
        cells = system.storage->getLocalCells();
        LOG4ESPP_DEBUG(theLogger, "local positions");
    }
    else
    {
        cells = system.storage->getRealCells();
        LOG4ESPP_DEBUG(theLogger, "real positions");
    }

    for (CellListIterator cit(cells); !cit.isDone(); ++cit)
    {
        LOG4ESPP_DEBUG(theLogger, "Particle " << cit->id() << ", position = " << cit->position());
    }
}

/****************************************************
** REGISTRATION WITH PYTHON
****************************************************/

void VelocityVerletLE::registerPython()
{
    using namespace espressopp::python;

    // Note: use noncopyable and no_init for abstract classes
    class_<VelocityVerletLE, bases<MDIntegrator>, boost::noncopyable>(
        "integrator_VelocityVerletLE", init<shared_ptr<System>, real>())
        .def("getTimers", &wrapGetTimers)
        .def("resetTimers", &VelocityVerletLE::resetTimers)
        .def("getNumResorts", &VelocityVerletLE::getNumResorts)
        .add_property("shear", &VelocityVerletLE::getShearRate, &VelocityVerletLE::setShearRate);
}
}  // namespace integrator
}  // namespace espressopp
