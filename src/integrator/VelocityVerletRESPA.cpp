/*
  Copyright (C) 2018
      Max Planck Institute for Polymer Research

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
#include "VelocityVerletRESPA.hpp"
// #include <iomanip>
#include "iterator/CellListIterator.hpp"
#include "interaction/Interaction.hpp"
#include "interaction/Potential.hpp"
#include "System.hpp"
#include "storage/Storage.hpp"
#include "mpi.hpp"

namespace espressopp {
  using namespace std;
  namespace integrator {
    using namespace interaction;
    using namespace iterator;
    using namespace esutil;
    using namespace boost::python;

    VelocityVerletRESPA::VelocityVerletRESPA(shared_ptr< System > system) : MDIntegrator(system)
    {
      resortFlag = true;
      maxDist    = 0.0;
      multistep = 1;
      dtlong = dt*multistep;
    }

    VelocityVerletRESPA::~VelocityVerletRESPA() {}

    void VelocityVerletRESPA::run(int nsteps)
    {
      int nResorts = 0;
      System& system = getSystemRef();
      storage::Storage& storage = *system.storage;
      real skinHalf = 0.5 * system.getSkin();
      dtlong = dt*multistep;

      runInit(); // signal

      // Before start make sure that particles are on the right processor
      if (resortFlag) {
        storage.decompose();
        maxDist = 0.0;
        resortFlag = false;
      }

      updateForces(true);
      aftCalcSlow(); // signal

      for (int i = 0; i < nsteps; i++) {
        integrate2(true);
        aftIntSlow(); // signal

        recalc1(); // signal
        updateForces(false);
        recalc2(); // signal

        for (int j = 0; j < multistep; j++) {
          befIntP(); // signal

          maxDist += integrate1();
          aftIntP(); // signal

          if (maxDist > skinHalf) resortFlag = true;
          if (resortFlag) {
            storage.decompose();
            maxDist  = 0.0;
            resortFlag = false;
            nResorts ++;
          }

          updateForces(false);
          befIntV(); // signal

          integrate2(false);
          aftIntV(); // signal
        }

        updateForces(true);
        aftCalcSlow(); // signal

        integrate2(true);
        aftIntSlow(); // signal
      }
    }

    real VelocityVerletRESPA::integrate1()
    {
      System& system = getSystemRef();
      CellList realCells = system.storage->getRealCells();

      // loop over particles
      real maxSqDist = 0.0;
      for(CellListIterator cit(realCells); !cit.isDone(); ++cit) {
        real sqDist = 0.0;
        real dtfm = 0.5 * dt / cit->mass();

        // propagate velocities
        cit->velocity() += dtfm * cit->force();

        // propagate positions
        Real3D deltaP = cit->velocity();

        deltaP *= dt;
        cit->position() += deltaP;
        sqDist += deltaP * deltaP;
        maxSqDist = std::max(maxSqDist, sqDist);
      }

      inIntP(maxSqDist); // signal

      real maxAllSqDist;
      mpi::all_reduce(*system.comm, maxSqDist, maxAllSqDist, boost::mpi::maximum<real>());
      return sqrt(maxAllSqDist);
    }

    void VelocityVerletRESPA::integrate2(bool slow)
    {
      System& system = getSystemRef();
      CellList realCells = system.storage->getRealCells();

      // loop over particles
      real half_dt = 0.0;
      if(slow) {
        half_dt = 0.5 * dtlong;
      }
      else {
        half_dt = 0.5 * dt;
      }
      for(CellListIterator cit(realCells); !cit.isDone(); ++cit) {
        real dtfm = half_dt / cit->mass();
        // propagate velocities
        cit->velocity() += dtfm * cit->force();
      }

      step++;
    }

    void VelocityVerletRESPA::calcForces(bool slow)
    {
      initForces();
      aftInitF(); // signal

      System& sys = getSystemRef();
      const InteractionList& srIL = sys.shortRangeInteractions;

      if(slow == true) {
        for (size_t i = 0; i < srIL.size(); i++) {
          if(srIL[i]->bondType() == NonbondedSlow) {
            srIL[i]->addForces();
          }
        }
      }
      else {
        for (size_t i = 0; i < srIL.size(); i++) {
          if(srIL[i]->bondType() != NonbondedSlow) {
            srIL[i]->addForces();
          }
        }
      }
    }

    void VelocityVerletRESPA::updateForces(bool slow)
    {
      storage::Storage& storage = *getSystemRef().storage;
      storage.updateGhosts();
      calcForces(slow);
      storage.collectGhostForces();

      if (slow == false) {
        aftCalcF(); // signal
      }
    }

    void VelocityVerletRESPA::initForces()
    {
      // forces are initialized for real + ghost particles
      System& system = getSystemRef();
      CellList localCells = system.storage->getLocalCells();

      for(CellListIterator cit(localCells); !cit.isDone(); ++cit) {
        cit->force() = 0.0;
        cit->drift() = 0.0;
      }
    }

    void VelocityVerletRESPA::setmultistep(int _multistep)
    {
      if (_multistep <= 0) {
        throw std::invalid_argument("multistep must be larger than zero!");
      }
      multistep = _multistep;
      dtlong = dt*multistep;
    }

    void VelocityVerletRESPA::setTimeStep(real _dt)
    {
      if (_dt == 0.0) {
        throw std::invalid_argument("Timestep 'dt' must be non-zero!");
      }
      dt = _dt;
      dtlong = dt*multistep;
    }

    /****************************************************
    ** REGISTRATION WITH PYTHON
    ****************************************************/

    void VelocityVerletRESPA::registerPython() {

      using namespace espressopp::python;

      // Note: use noncopyable and no_init for abstract classes
      class_<VelocityVerletRESPA, bases<MDIntegrator>, boost::noncopyable >
      ("integrator_VelocityVerletRESPA", init< shared_ptr<System> >())
      .def("setmultistep", &VelocityVerletRESPA::setmultistep)
      .def("getmultistep", &VelocityVerletRESPA::getmultistep)
      .add_property("multistep", &VelocityVerletRESPA::getmultistep, &VelocityVerletRESPA::setmultistep)
      ;
    }
  }
}
