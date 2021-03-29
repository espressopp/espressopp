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

#include "python.hpp"
#include "VelocityVerlet.hpp"
// #include "iterator/CellListIterator.hpp"
// #include "interaction/Interaction.hpp"
// #include "interaction/Potential.hpp"
// #include "System.hpp"
// #include "storage/Storage.hpp"
// #include "mpi.hpp"

#include <iomanip>

#ifdef VTRACE
#include "vampirtrace/vt_user.h"
#else
#define VT_TRACER( name)
#endif

namespace espressopp { namespace vec {
  namespace integrator {
    // using namespace interaction;
    // using namespace iterator;
    // using namespace esutil;

    LOG4ESPP_LOGGER(VelocityVerlet::theLogger, "VelocityVerlet");

    VelocityVerlet::VelocityVerlet(
      shared_ptr<System> system,
      shared_ptr<vec::storage::StorageVec> storageVec
      ) : MDIntegratorVec(system, storageVec)
    {
      // LOG4ESPP_INFO(theLogger, "construct VelocityVerlet");
      // resortFlag = true;
      // maxDist    = 0.0;
      // nResorts   = 0;
      std::cout << "vec::integrator::VelocityVerlet::" << __FUNCTION__ << std::endl;
    }

    void VelocityVerlet::run(int nsteps)
    {
      // VEC_DEBUG_MSG("VelocityVerlet::run");
      // vec::utils::runAsVecThread([this,nsteps]{this->run_(nsteps);});
    }

#if 0
    void VelocityVerlet::run_(int nsteps)
    {
      // VEC_DEBUG_MSG_THREAD("VelocityVerlet::hpx_run");

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
      if (resortFlag) {
        const real time = timeIntegrate.getElapsedTime();

        LOG4ESPP_INFO(theLogger, "resort particles");
        storageVec->decomposeVec();
        maxDist = 0.0;
        resortFlag = false;

        timeOtherInitResort += timeIntegrate.getElapsedTime() - time;
      }
      {
        const real time = timeIntegrate.getElapsedTime();
        storageVec->loadCells();
        timeOtherLoadCells += timeIntegrate.getElapsedTime() - time;
      }
      bool recalcForces = true;  // TODO: more intelligent
      if (recalcForces) {
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

            storageVec->unloadCells();

            // storage.decompose();
            storageVec->decomposeVec();

            storageVec->loadCells();

            maxDist  = 0.0;
            resortFlag = false;
            nResorts ++;

            timeResort += timeIntegrate.getElapsedTime() - time;
          }

          #if 0
          else
          {
            /// FIXME: Simplify by updating only the outgoing commCells
            /// copy position/velocities from array to list
            /// if resort was done in the previous step, data in array should already be updated
            auto& virtualStorage = storageVec->virtualStorage;
            auto f = [&virtualStorage](size_t i){
              auto& vs = virtualStorage[i];
              vs.particles.updateToPositionVelocity(vs.localCells, true); /// integrate1 modified only real cells
            };
            const size_t nvs = virtualStorage.size();
            utils::parallelForLoop(0, nvs, f);
            // for(size_t i=0; i<nvs; i++) f(i);
          }
          #endif

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
      }
      {
        const real time = timeIntegrate.getElapsedTime();
        storageVec->unloadCells();
        timeOtherUnloadCells += timeIntegrate.getElapsedTime() - time;
      }

      timeRun = timeIntegrate.getElapsedTime();
      timeLost = timeRun - (timeForceComp[0] + timeForceComp[1] + timeForceComp[2] +
                 timeComm1 + timeComm2 + timeInt1 + timeInt2 + timeResort);
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
      // Implement force update here
      // Initial implementation: blocking update following original

      real time;

      time = timeIntegrate.getElapsedTime();
      storageVec->updateGhostsBlocking();
      timeComm1 += timeIntegrate.getElapsedTime() - time;

      time = timeIntegrate.getElapsedTime();
      calcForces();
      timeForce += timeIntegrate.getElapsedTime() - time;

      time = timeIntegrate.getElapsedTime();
      storageVec->collectGhostForcesBlocking();
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

      for(CellListIterator cit(localCells); !cit.isDone(); ++cit) {
        cit->force() = 0.0;
        cit->drift() = 0.0;   // Can in principle be commented, when drift is not used.
      }
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
      timeForceIntra = 0.0;
      for(int i = 0; i < 100; i++)
        timeForceIntraComp[i] = 0.0;
      timeOverlap= 0.0;

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

    static object wrapGetTimers(class VelocityVerlet* obj) {
      real tms[14];
      obj->loadTimers(tms);
      return boost::python::make_tuple(
          tms[0],
          tms[1],
          tms[2],
          tms[3],
          tms[4],
          tms[5],
          tms[6],
          tms[7],
          tms[8],
          tms[9],
          tms[10],
          tms[11],
          tms[12],
          tms[13]);
    }

    void VelocityVerlet::loadTimers(real t[14]) {
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
      t[10] = timeOverlap;
      t[11] = timeForceIntraComp[0];
      t[12] = timeForceIntraComp[1];
      t[13] = timeForceIntraComp[2];
    }

    int VelocityVerlet::getNumResorts() const
    {
      return nResorts;
    }

    void VelocityVerlet::loadOtherTimers(real *t) {
      t[0] = timeOtherInitResort;
      t[1] = timeOtherInitForcesPlist;
      t[2] = timeOtherLoadCells;
      t[3] = timeOtherRecalcForces;
      t[4] = timeOtherInitForcesParray;
      t[5] = timeOtherUnloadCells;
      t[6] = timeOtherAftCalcF;
    }

    static object wrapGetOtherTimers(class VelocityVerlet* obj) {
      real tms[7];
      obj->loadOtherTimers(tms);
      return boost::python::make_tuple(
          tms[0],
          tms[1],
          tms[2],
          tms[3],
          tms[4],
          tms[5],
          tms[6]);
    }
#endif

    /****************************************************
    ** REGISTRATION WITH PYTHON
    ****************************************************/

    void VelocityVerlet::registerPython() {

      using namespace espressopp::python;

      // Note: use noncopyable and no_init for abstract classes
      class_< vec::integrator::VelocityVerlet,
              bases<espressopp::integrator::MDIntegrator, MDIntegratorVec>,
              boost::noncopyable >
        ("vec_integrator_VelocityVerlet", init< shared_ptr<System>, shared_ptr< storage::StorageVec > >())
        // .def("run", &vec::integrator::VelocityVerlet::run)
        // .def("getTimers", &wrapGetTimers)
        // .def("getOtherTimers", &wrapGetOtherTimers)
        // .def("resetTimers", &VelocityVerlet::resetTimers)
        // .def("getNumResorts", &VelocityVerlet::getNumResorts)
        ;
    }
  }
}}
