/*
  Copyright (C) 2021
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

#include "vec/Vectorization.hpp"
#include "vec/storage/StorageVec.hpp"
#include "vec/iterator/ParticleArrayIterator.hpp"
#include "VelocityVerlet.hpp"

#include "python.hpp"
#include "iterator/CellListIterator.hpp"
#include "System.hpp"
#include "storage/Storage.hpp"

#include <iomanip>
#include <numeric>

#ifdef VTRACE
#include "vampirtrace/vt_user.h"
#else
#define VT_TRACER( name)
#endif

namespace espressopp { namespace vec {
  namespace integrator {

    using namespace interaction;
    using namespace iterator;
    using namespace esutil;

    using espressopp::System;
    using espressopp::storage::Storage;
    using espressopp::vec::storage::StorageVec;

    ///////////////////////////////////////////////////////////////////////////////////////////////
    // VelocityVerletBase

    LOG4ESPP_LOGGER(VelocityVerletBase::theLogger, "VelocityVerletBase");

    VelocityVerletBase::VelocityVerletBase(std::shared_ptr<System> system)
      : MDIntegratorVec(system)
    {
      LOG4ESPP_INFO(theLogger, "construct VelocityVerletBase");
      resortFlag = true;
      maxDist    = 0.0;
      nResorts   = 0;
    }

    void VelocityVerletBase::run(int nsteps)
    {
      if(!(getSystem()->vectorization->storageVec)) {
        throw std::runtime_error("Vectorization has no storageVec");
      }

      nResorts = 0;
      real time;
      timeIntegrate.reset();
      resetTimers();

      System& system = getSystemRef();
      Storage& storage = *system.storage;
      StorageVec& storageVec = *getSystem()->vectorization->storageVec;
      const real skinHalf = 0.5 * system.getSkin();

      // signal
      MDIntegratorVec::runInit();

      // Before start make sure that particles are on the right processor
      if (resortFlag) {
        time = timeIntegrate.getElapsedTime();
        LOG4ESPP_INFO(theLogger, "resort particles");
        storage.decompose();
        maxDist = 0.0;
        resortFlag = false;
        timeResort += timeIntegrate.getElapsedTime();
      }

      bool recalcForces = true;  // TODO: more intelligent
      if (recalcForces) {
        LOG4ESPP_INFO(theLogger, "recalc forces before starting main integration loop");

        initForcesPlist();

        // signal
        MDIntegratorVec::recalc1();

        updateForces();

        // signal
        MDIntegratorVec::recalc2();
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

          timeInt1 += timeIntegrate.getElapsedTime() - time;
        }

        if (maxDist > skinHalf) resortFlag = true;

        if (resortFlag)
        {
          const real time = timeIntegrate.getElapsedTime();

          storageVec.unloadCells();
          storage.decompose();

          maxDist  = 0.0;
          resortFlag = false;
          nResorts ++;

          timeResort += timeIntegrate.getElapsedTime() - time;
        }

        {
          updateForces();
        }

        {
          const real time = timeIntegrate.getElapsedTime();

          integrate2();

          timeInt2 += timeIntegrate.getElapsedTime() - time;
        }
        step++;
      }

      {
        // since load is counted in timeResort, unload should also be counted there
        const real time = timeIntegrate.getElapsedTime();
        storageVec.unloadCells();
        timeResort += timeIntegrate.getElapsedTime() - time;
      }

      timeRun = timeIntegrate.getElapsedTime();
      timeLost = timeRun - (timeForceComp[0] + timeForceComp[1] + timeForceComp[2] +
                 timeComm1 + timeComm2 + timeInt1 + timeInt2 + timeResort);
    }

    real VelocityVerletBase::integrate1()
    {
      real maxSqDist = 0.0;
      {
        auto& particles = getSystem()->vectorization->particles;
        for(iterator::ParticleArrayIterator pit(particles, true); pit.isValid(); ++pit)
        {
          const real dtfm = 0.5 * dt / pit.mass();
          pit.v_x() += dtfm * pit.f_x();
          pit.v_y() += dtfm * pit.f_y();
          pit.v_z() += dtfm * pit.f_z();
          const real dp_x = pit.v_x() * dt;
          const real dp_y = pit.v_y() * dt;
          const real dp_z = pit.v_z() * dt;
          pit.p_x() += dp_x;
          pit.p_y() += dp_y;
          pit.p_z() += dp_z;
          real sqDist = (dp_x*dp_x) + (dp_y*dp_y) + (dp_z*dp_z);
          maxSqDist = std::max(maxSqDist, sqDist);
        }
      }

      return maxSqDist;
    }

    void VelocityVerletBase::integrate2()
    {
      auto& particles = getSystem()->vectorization->particles;
      {
        for(iterator::ParticleArrayIterator pit(particles, true); pit.isValid(); ++pit)
        {
          const real dtfm = 0.5 * dt / pit.mass();
          pit.v_x() += dtfm * pit.f_x();
          pit.v_y() += dtfm * pit.f_y();
          pit.v_z() += dtfm * pit.f_z();
        }
      }
    }

    void VelocityVerletBase::calcForces()
    {
      initForcesParray();
      {
        // TODO: Might need to place interaction list in VecRuntime
        System& sys = getSystemRef();
        const espressopp::interaction::InteractionList& srIL = sys.shortRangeInteractions;

        for (size_t i = 0; i < srIL.size(); i++) {
        LOG4ESPP_INFO(theLogger, "compute forces for srIL " << i << " of " << srIL.size());
          real time;
          time = timeIntegrate.getElapsedTime();
          srIL[i]->addForces();
          timeForceComp[i] += timeIntegrate.getElapsedTime() - time;
        }
        // aftCalcFLocal();
      }
    }

    void VelocityVerletBase::updateForces()
    {
      auto storageVec = getSystem()->vectorization->storageVec;
      real time;

      time = timeIntegrate.getElapsedTime();
      storageVec->updateGhostsVec();
      timeComm1 += timeIntegrate.getElapsedTime() - time;

      time = timeIntegrate.getElapsedTime();
      calcForces();
      timeForce += timeIntegrate.getElapsedTime() - time;

      time = timeIntegrate.getElapsedTime();
      storageVec->collectGhostForcesVec();
      timeComm2 += timeIntegrate.getElapsedTime() - time;

      MDIntegratorVec::aftCalcF();
    }

    void VelocityVerletBase::initForcesPlist()
    {
      // forces are initialized for real + ghost particles

      System& system = getSystemRef();
      CellList localCells = system.storage->getLocalCells();

      LOG4ESPP_INFO(theLogger, "init forces for real + ghost particles");

      for(espressopp::iterator::CellListIterator cit(localCells); !cit.isDone(); ++cit) {
        cit->force() = 0.0;
        cit->drift() = 0.0;   // Can in principle be commented, when drift is not used.
      }
    }

    void VelocityVerletBase::initForcesParray()
    {
      getSystem()->vectorization->zeroForces();
    }

    void VelocityVerletBase::resetTimers() {
      timeForce  = 0.0;
      for(int i = 0; i < 100; i++)
        timeForceComp[i] = 0.0;
      timeComm1  = 0.0;
      timeComm2  = 0.0;
      timeInt1   = 0.0;
      timeInt2   = 0.0;
      timeResort = 0.0;
    }

    using namespace boost::python;

    static object wrapGetTimers(class VelocityVerletBase* obj) {
      real tms[14];
      obj->loadTimers(tms);
      return boost::python::make_tuple(
          tms[0]
        , tms[1]
        , tms[2]
        , tms[3]
        , tms[4]
        , tms[5]
        , tms[6]
        , tms[7]
        , tms[8]
        , tms[9]
        );
    }

    void VelocityVerletBase::loadTimers(real t[10]) {
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

    int VelocityVerletBase::getNumResorts() const
    {
      return nResorts;
    }

    /****************************************************
    ** REGISTRATION WITH PYTHON
    ****************************************************/

    void VelocityVerletBase::registerPython() {

      using namespace espressopp::python;

      // Note: use noncopyable and no_init for abstract classes
      class_< vec::integrator::VelocityVerletBase,
              bases<espressopp::integrator::MDIntegrator, MDIntegratorVec>,
              boost::noncopyable >
        ("vec_integrator_VelocityVerletBase", init< std::shared_ptr<System> >())
        .def("run", &vec::integrator::VelocityVerletBase::run)
        .def("getTimers", &wrapGetTimers)
        .def("resetTimers", &VelocityVerletBase::resetTimers)
        .def("getNumResorts", &VelocityVerletBase::getNumResorts)
        ;
    }

    ///////////////////////////////////////////////////////////////////////////////////////////////
    // VelocityVerlet

    LOG4ESPP_LOGGER(VelocityVerlet::theLogger, "VelocityVerlet");

    real VelocityVerlet::integrate1()
    {
      real maxSqDist = 0.0;
      {
        auto& particles                    = getSystem()->vectorization->particles;
        const auto& realCells              = particles.realCells();
        const size_t* __restrict cellRange = particles.cellRange().data();
        const size_t* __restrict sizes     = particles.sizes().data();
        for(const auto& rcell: realCells)
        {
          const size_t start = cellRange[rcell];
          const size_t size  = sizes[rcell];

          real* __restrict p_x = &(particles.p_x[start]);
          real* __restrict p_y = &(particles.p_y[start]);
          real* __restrict p_z = &(particles.p_z[start]);
          real* __restrict v_x = &(particles.v_x[start]);
          real* __restrict v_y = &(particles.v_y[start]);
          real* __restrict v_z = &(particles.v_z[start]);
          const real* __restrict f_x = &(particles.f_x[start]);
          const real* __restrict f_y = &(particles.f_y[start]);
          const real* __restrict f_z = &(particles.f_z[start]);
          const real* __restrict mass = &(particles.mass[start]);

          #pragma vector always
          #pragma vector aligned
          #pragma ivdep
          for(size_t ip=0; ip<size; ip++)
          {
            const real dtfm = 0.5 * dt / mass[ip];
            v_x[ip] += dtfm * f_x[ip];
            v_y[ip] += dtfm * f_y[ip];
            v_z[ip] += dtfm * f_z[ip];
            const real dp_x = v_x[ip] * dt;
            const real dp_y = v_y[ip] * dt;
            const real dp_z = v_z[ip] * dt;
            p_x[ip] += dp_x;
            p_y[ip] += dp_y;
            p_z[ip] += dp_z;
            real sqDist = (dp_x*dp_x) + (dp_y*dp_y) + (dp_z*dp_z);
            maxSqDist = std::max(maxSqDist, sqDist);
          }
        }
      }
      return maxSqDist;
    }

    void VelocityVerlet::integrate2()
    {
      auto& particles = getSystem()->vectorization->particles;
      {
        const auto& realCells              = particles.realCells();
        const size_t* __restrict cellRange = particles.cellRange().data();
        const size_t* __restrict sizes     = particles.sizes().data();

        for(const auto& rcell: realCells)
        {
          const size_t start = cellRange[rcell];
          const size_t size  = cellRange[rcell+1]-start;

          real* __restrict v_x = &(particles.v_x[start]);
          real* __restrict v_y = &(particles.v_y[start]);
          real* __restrict v_z = &(particles.v_z[start]);
          const real* __restrict f_x = &(particles.f_x[start]);
          const real* __restrict f_y = &(particles.f_y[start]);
          const real* __restrict f_z = &(particles.f_z[start]);
          const real* __restrict mass = &(particles.mass[start]);

          #pragma vector always
          #pragma vector aligned
          #pragma ivdep
          for(size_t ip=0; ip<size; ip++)
          {
            const real dtfm = 0.5 * dt / mass[ip];
            v_x[ip] += dtfm * f_x[ip];
            v_y[ip] += dtfm * f_y[ip];
            v_z[ip] += dtfm * f_z[ip];
          }
        }
      }
    }

    /****************************************************
    ** REGISTRATION WITH PYTHON
    ****************************************************/

    void VelocityVerlet::registerPython() {
      using namespace espressopp::python;

      VelocityVerletBase::registerPython();

      class_< vec::integrator::VelocityVerlet,
              bases<VelocityVerletBase>,
              boost::noncopyable >
        ("vec_integrator_VelocityVerlet", init< std::shared_ptr<System> >())
        ;
    }

  }
}}
