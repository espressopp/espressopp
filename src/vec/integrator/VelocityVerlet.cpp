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

// #include "vec/include/logging.hpp"
// #include "vec/include/errors.hpp"
// #include "vec/utils/algorithms/for_loop.hpp"
// #include "vec/utils/algorithms/transform_reduce.hpp"
// #include "vec/utils/multithreading.hpp"

#include "vec/Vectorization.hpp"

#include "python.hpp"
#include "VelocityVerlet.hpp"
#include "iterator/CellListIterator.hpp"
// #include "interaction/Interaction.hpp"
// #include "interaction/Potential.hpp"
#include "System.hpp"
#include "storage/Storage.hpp"
// #include "mpi.hpp"

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

    LOG4ESPP_LOGGER(VelocityVerlet::theLogger, "VelocityVerlet");

    VelocityVerlet::VelocityVerlet(shared_ptr<Vectorization> vectorization)
      : MDIntegratorVec(vectorization)
    {
      LOG4ESPP_INFO(theLogger, "construct VelocityVerlet");
      resortFlag = true;
      maxDist    = 0.0;
      nResorts   = 0;
    }

    void VelocityVerlet::run(int nsteps)
    {
      if(!(vectorization->storageVec)) {
        throw std::runtime_error("Vectorization has no storageVec");
      }

      nResorts = 0;
      real time;
      timeIntegrate.reset();
      resetTimers();

      System& system = getSystemRef();
      Storage& storage = *system.storage;
      StorageVec& storageVec = *vectorization->storageVec;
      const real skinHalf = 0.5 * system.getSkin();

      // signal
      MDIntegratorVec::runInit();

      // Before start make sure that particles are on the right processor
      if (resortFlag) {
        LOG4ESPP_INFO(theLogger, "resort particles");
        storage.decompose();
        maxDist = 0.0;
        resortFlag = false;
      }

      {
        storageVec.loadCells();
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

          storageVec.loadCells();

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
        storageVec.unloadCells();
      }

      timeRun = timeIntegrate.getElapsedTime();
      timeLost = timeRun - (timeForceComp[0] + timeForceComp[1] + timeForceComp[2] +
                 timeComm1 + timeComm2 + timeInt1 + timeInt2 + timeResort);
    }

    real VelocityVerlet::integrate1()
    {
      auto& particles                    = vectorization->particles;
      const auto& realCells              = particles.realCells();
      const size_t* __restrict cellRange = particles.cellRange().data();
      const size_t* __restrict sizes     = particles.sizes().data();

      real maxSqDist = 0.0;

      // first-half integration
      // apply integration scheme on every cell and obtain the maximum displacement value
      if(particles.mode_aos())
      {
        for(const auto& rcell: realCells)
        {
          size_t start = cellRange[rcell];
          size_t size  = sizes[rcell];

          using espressopp::vec::Real3DInt;
          using espressopp::vec::Real4D;
          Real3DInt*    __restrict p    = particles.position.data() + start;
          Real4D*       __restrict v    = particles.velocity.data() + start;
          const Real4D* __restrict f    = particles.force.data()    + start;
          const real*   __restrict mass = particles.mass.data()     + start;

          #pragma vector always
          #pragma vector aligned
          #pragma ivdep
          for(size_t ip=0; ip<size; ip++)
          {
            /// TODO: transform division by mass to multiplication by reciprocal or just store
            /// dtfm as an internal array that gets updated by loadCells()
            const real dtfm = 0.5 * dt / mass[ip];

            v[ip].x += dtfm * f[ip].x;
            v[ip].y += dtfm * f[ip].y;
            v[ip].z += dtfm * f[ip].z;

            const real dp_x = v[ip].x * dt;
            const real dp_y = v[ip].y * dt;
            const real dp_z = v[ip].z * dt;

            p[ip].x += dp_x;
            p[ip].y += dp_y;
            p[ip].z += dp_z;

            real sqDist = (dp_x*dp_x) + (dp_y*dp_y) + (dp_z*dp_z);
            maxSqDist = std::max(maxSqDist, sqDist);
          }
        }
      }
      else
      {
        for(const auto& rcell: realCells)
        {
          size_t start = cellRange[rcell];
          size_t size  = sizes[rcell];

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
            /// TODO: transform division by mass to multiplication by reciprocal or just store
            /// dtfm as an internal array that gets updated by loadCells()
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
      auto& particles                    = vectorization->particles;
      const auto& realCells              = particles.realCells();
      const size_t* __restrict cellRange = particles.cellRange().data();
      const size_t* __restrict sizes     = particles.sizes().data();

      if(particles.mode_aos())
      {
        for(const auto& rcell: realCells)
        {
          const size_t start = cellRange[rcell];
          const size_t size  = sizes[rcell];

          using espressopp::vec::Real4D;
          Real4D*       __restrict v    = particles.velocity.data() + start;
          const Real4D* __restrict f    = particles.force.data()    + start;
          const real*   __restrict mass = particles.mass.data()     + start;

          #pragma vector always
          #pragma vector aligned
          #pragma ivdep
          for(size_t ip=0; ip<size; ip++)
          {
            /// TODO: transform division by mass to multiplication by reciprocal or fixed dtfm array
            const real dtfm = 0.5 * dt / mass[ip];
            v[ip].x += dtfm * f[ip].x;
            v[ip].y += dtfm * f[ip].y;
            v[ip].z += dtfm * f[ip].z;
          }
        }
      }
      else
      {
        for(const auto& rcell: realCells)
        {
          const size_t start = cellRange[rcell];
          const size_t size  = sizes[rcell];

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
            /// TODO: transform division by mass to multiplication by reciprocal or fixed dtfm array
            const real dtfm = 0.5 * dt / mass[ip];
            v_x[ip] += dtfm * f_x[ip];
            v_y[ip] += dtfm * f_y[ip];
            v_z[ip] += dtfm * f_z[ip];
          }
        }
      }
    }

    void VelocityVerlet::calcForces()
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

    void VelocityVerlet::updateForces()
    {
      auto storageVec = vectorization->storageVec;
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

    void VelocityVerlet::initForcesPlist()
    {
      // forces are initialized for real + ghost particles

      System& system = getSystemRef();
      CellList localCells = system.storage->getLocalCells();

      LOG4ESPP_INFO(theLogger, "init forces for real + ghost particles");

      for(CellListIterator cit(localCells); !cit.isDone(); ++cit) {
        cit->force() = 0.0;
        cit->drift() = 0.0;   // Can in principle be commented, when drift is not used.
      }
    }

    void VelocityVerlet::initForcesParray()
    {
      vectorization->zeroForces();
    }

    void VelocityVerlet::resetTimers() {
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

    static object wrapGetTimers(class VelocityVerlet* obj) {
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

    void VelocityVerlet::loadTimers(real t[10]) {
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

    int VelocityVerlet::getNumResorts() const
    {
      return nResorts;
    }

    /****************************************************
    ** REGISTRATION WITH PYTHON
    ****************************************************/

    void VelocityVerlet::registerPython() {

      using namespace espressopp::python;

      // Note: use noncopyable and no_init for abstract classes
      class_< vec::integrator::VelocityVerlet,
              bases<espressopp::integrator::MDIntegrator, MDIntegratorVec>,
              boost::noncopyable >
        ("vec_integrator_VelocityVerlet", init< shared_ptr<Vectorization> >())
        .def("run", &vec::integrator::VelocityVerlet::run)
        .def("getTimers", &wrapGetTimers)
        .def("resetTimers", &VelocityVerlet::resetTimers)
        .def("getNumResorts", &VelocityVerlet::getNumResorts)
        ;
    }
  }
}}
