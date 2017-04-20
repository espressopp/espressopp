/*
  Copyright (C) 2017
      Jakub Krajniak (jkrajniak at gmail.com)
  Copyright (C) 2012,2013
      Max Planck Institute for Polymer Research
  Copyright (C) 2008,2009,2010,2011
      Max-Planck-Institute for Polymer Research & Fraunhofer SCAI
  
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
#include "VelocityVerlet.hpp"
#include <iomanip>
#include "iterator/CellListIterator.hpp"
#include "interaction/Interaction.hpp"
#include "interaction/Potential.hpp"
#include "System.hpp"
#include "storage/Storage.hpp"
#include "mpi.hpp"

#ifdef VTRACE
#include "vampirtrace/vt_user.h"
#else
#define VT_TRACER( name)
#endif

namespace espressopp {
  using namespace std;
  namespace integrator {
    using namespace interaction;
    using namespace iterator;
    using namespace esutil;

    LOG4ESPP_LOGGER(VelocityVerlet::theLogger, "VelocityVerlet");

    VelocityVerlet::VelocityVerlet(shared_ptr< System > system) : MDIntegrator(system)
    {
      LOG4ESPP_INFO(theLogger, "construct VelocityVerlet");
      resortFlag = true;
      maxDist    = 0.0;
    }

    VelocityVerlet::~VelocityVerlet()
    {
      LOG4ESPP_INFO(theLogger, "free VelocityVerlet");
    }

    void VelocityVerlet::run(int nsteps)
    {
      VT_TRACER("run");
      int nResorts = 0;
      real time;
      timeIntegrate.reset();
      resetTimers();
      System& system = getSystemRef();
      storage::Storage& storage = *system.storage;
      real skinHalf = 0.5 * system.getSkin();

      // signal
      runInit();

      // Before start make sure that particles are on the right processor
      if (resortFlag) {
        VT_TRACER("resort");
        // time = timeIntegrate.getElapsedTime();
        LOG4ESPP_INFO(theLogger, "resort particles");
        storage.decompose();
        maxDist = 0.0;
        resortFlag = false;
        // timeResort += timeIntegrate.getElapsedTime();
      }

      bool recalcForces = true;  // TODO: more intelligent

      if (recalcForces) {
        LOG4ESPP_INFO(theLogger, "recalc forces before starting main integration loop");

        // signal
        recalc1();

        updateForces();
        if (LOG4ESPP_DEBUG_ON(theLogger)) {
            // printForces(false);   // forces are reduced to real particles
        }

        // signal
        recalc2();
      }

      LOG4ESPP_INFO(theLogger, "starting main integration loop (nsteps=" << nsteps << ")");
  
      for (int i = 0; i < nsteps; i++) {
        LOG4ESPP_INFO(theLogger, "Next step " << i << " of " << nsteps << " starts");

        //saveOldPos(); // save particle positions needed for constraints

        // signal
        befIntP();

        time = timeIntegrate.getElapsedTime();
        LOG4ESPP_INFO(theLogger, "updating positions and velocities")
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

        LOG4ESPP_INFO(theLogger, "maxDist = " << maxDist << ", skin/2 = " << skinHalf);

        if (maxDist > skinHalf) resortFlag = true;
        
        if (resortFlag) {
            VT_TRACER("resort1");
            time = timeIntegrate.getElapsedTime();
            LOG4ESPP_INFO(theLogger, "step " << i << ": resort particles");
            storage.decompose();
            maxDist  = 0.0;
            resortFlag = false;
            nResorts ++;
            timeResort += timeIntegrate.getElapsedTime() - time;
        }

        LOG4ESPP_INFO(theLogger, "updating forces")
        updateForces();

        // signal
        befIntV();

        time = timeIntegrate.getElapsedTime();
        integrate2();
        timeInt2 += timeIntegrate.getElapsedTime() - time;

        // signal
        aftIntV();
      }

      timeRun = timeIntegrate.getElapsedTime();
      timeLost = timeRun - (timeForceComp[0] + timeForceComp[1] + timeForceComp[2] +
                 timeComm1 + timeComm2 + timeInt1 + timeInt2 + timeResort);

      LOG4ESPP_INFO(theLogger, "finished run");
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
      real tms[10];
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
          tms[9]);
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

    void VelocityVerlet::printTimers() {

      using namespace std;
      real pct;

      cout << endl;
      cout << "run = " << setiosflags(ios::fixed) << setprecision(1) << timeRun << endl;
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

    real VelocityVerlet::integrate1()
    {
      System& system = getSystemRef();
      CellList realCells = system.storage->getRealCells();

      // loop over all particles of the local cells
      int count = 0;
      real maxSqDist = 0.0; // maximal square distance a particle moves
      for(CellListIterator cit(realCells); !cit.isDone(); ++cit) {
        real sqDist = 0.0;
        LOG4ESPP_INFO(theLogger, "updating first half step of velocities and full step of positions")
        LOG4ESPP_DEBUG(theLogger, "Particle " << cit->id() << 
                ", pos = " << cit->position() <<
                ", v = " << cit->velocity() << 
                ", f = " << cit->force());

        /* more precise for DEBUG:
        printf("Particle %d, pos = %16.12f %16.12f %16.12f, v = %16.12f, %16.12f %16.12f, f = %16.12f %16.12f %16.12f\n",
            cit->p.id, cit->r.p[0], cit->r.p[1], cit->r.p[2],
                cit->m.v[0], cit->m.v[1], cit->m.v[2],
            cit->f.f[0], cit->f.f[1], cit->f.f[2]);
        */

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

      LOG4ESPP_INFO(theLogger, "moved " << count << " particles in integrate1" <<
		    ", max move local = " << sqrt(maxSqDist) <<
		    ", global = " << sqrt(maxAllSqDist));
      
      return sqrt(maxAllSqDist);
    }

    void VelocityVerlet::integrate2()
    {
      LOG4ESPP_INFO(theLogger, "updating second half step of velocities")
      System& system = getSystemRef();
      CellList realCells = system.storage->getRealCells();

      // loop over all particles of the local cells
      real half_dt = 0.5 * dt; 
      for(CellListIterator cit(realCells); !cit.isDone(); ++cit) {
        real dtfm = half_dt / cit->mass();
        /* Propagate velocities: v(t+0.5*dt) = v(t) + 0.5*dt * f(t) */
        cit->velocity() += dtfm * cit->force();
      }
      
      step++;
    }

    void VelocityVerlet::calcForces()
    {
      VT_TRACER("forces");

      LOG4ESPP_INFO(theLogger, "calculate forces");

      initForces();

      // signal
      aftInitF();

      System& sys = getSystemRef();
      const InteractionList& srIL = sys.shortRangeInteractions;

      for (size_t i = 0; i < srIL.size(); i++) {
	    LOG4ESPP_INFO(theLogger, "compute forces for srIL " << i << " of " << srIL.size());
        real time;
        time = timeIntegrate.getElapsedTime();
        srIL[i]->addForces();
        timeForceComp[i] += timeIntegrate.getElapsedTime() - time;
      }
    }

    void VelocityVerlet::updateForces()
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

    void VelocityVerlet::initForces()
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

    void VelocityVerlet::printForces(bool withGhosts)
    {
      // print forces of real + ghost particles

      System& system = getSystemRef();
      CellList cells;

      if (withGhosts) {
	    cells = system.storage->getLocalCells();
	    LOG4ESPP_DEBUG(theLogger, "local forces");
      } else {
	    cells = system.storage->getRealCells();
	    LOG4ESPP_DEBUG(theLogger, "real forces");
      }
  
      for(CellListIterator cit(cells); !cit.isDone(); ++cit) {
	    LOG4ESPP_DEBUG(theLogger, "Particle " << cit->id() << ", force = " << cit->force());
      }
    }

    void VelocityVerlet::printPositions(bool withGhosts)
    {
      // print positions of real + ghost particles

      System& system = getSystemRef();
      CellList cells;

      if (withGhosts) {
	    cells = system.storage->getLocalCells();
	    LOG4ESPP_DEBUG(theLogger, "local positions");
      } else {
	    cells = system.storage->getRealCells();
	    LOG4ESPP_DEBUG(theLogger, "real positions");
      }
  
      for(CellListIterator cit(cells); !cit.isDone(); ++cit) {
	    LOG4ESPP_DEBUG(theLogger, "Particle " << cit->id() << ", position = " << cit->position());
      }
    }

    /****************************************************
    ** REGISTRATION WITH PYTHON
    ****************************************************/

    void VelocityVerlet::registerPython() {

      using namespace espressopp::python;

      // Note: use noncopyable and no_init for abstract classes
      class_<VelocityVerlet, bases<MDIntegrator>, boost::noncopyable >
        ("integrator_VelocityVerlet", init< shared_ptr<System> >())
        .def("getTimers", &wrapGetTimers)
        .def("resetTimers", &VelocityVerlet::resetTimers)
        ;
    }
  }
}
