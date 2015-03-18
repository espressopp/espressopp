/*
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

#include "python.hpp"

#include "VelocityVerletOnGroup.hpp"

#include "LangevinThermostat.hpp"

#include "iterator/CellListIterator.hpp"
#include "interaction/Interaction.hpp"
#include "interaction/Potential.hpp"
#include "System.hpp"
#include "storage/Storage.hpp"
#include "mpi.hpp"

namespace espressopp {
  namespace integrator {

    using namespace interaction;
    using namespace iterator;
    using namespace esutil;

    VelocityVerletOnGroup::VelocityVerletOnGroup(shared_ptr< System > system,
            shared_ptr<class espressopp::ParticleGroup> group_) : MDIntegrator(system), group(group_)
    {
      LOG4ESPP_INFO(theLogger, "construct VelocityVerletOnGroup");

      resortFlag = true;
      maxDist = 0.0;
    }

    /*****************************************************************************/

    VelocityVerletOnGroup::~VelocityVerletOnGroup()
    {
      LOG4ESPP_INFO(theLogger, "free VelocityVerletOnGroup");
    }

    /*****************************************************************************/

    void VelocityVerletOnGroup::setLangevin(shared_ptr< LangevinThermostat > _langevin)
    {
      LOG4ESPP_INFO(theLogger, "set Langevin thermostat");
      langevin = _langevin;
    }

    /*****************************************************************************/

    void VelocityVerletOnGroup::run(int nsteps)
    {
      int nResorts = 0;

      real      time;

      timeIntegrate.reset();

      resetTimers();

      System& system = getSystemRef();

      storage::Storage& storage = *system.storage;

      //if (langevin) langevin->initialize(dt);

      // no more needed: setUp();

      // Before start make sure that particles are on the right processor

      if (resortFlag) {
        // time = timeIntegrate.getElapsedTime();
	LOG4ESPP_INFO(theLogger, "resort particles");
	storage.decompose();
	LOG4ESPP_INFO(theLogger, "particles resort");
	maxDist = 0.0;
	resortFlag = false;
        // timeResort += timeIntegrate.getElapsedTime();
      }

      bool recalcForces = true;  // TODO: more intelligent

      if (recalcForces) {

	LOG4ESPP_INFO(theLogger, "recalc Forces");

	if (langevin) langevin->heatUp();

        updateForces();

	if (LOG4ESPP_DEBUG_ON(theLogger)) {
	  // printForces(false);   // forces are reduced to real particles
	}

	if (langevin) langevin->coolDown();
      }

      LOG4ESPP_INFO(theLogger, "run " << nsteps << " iterations");
  
      real skinHalf = 0.5 * system.getSkin();

      for (int i = 0; i < nsteps; i++) {

	LOG4ESPP_INFO(theLogger, "Next step " << i << " of " << nsteps << " starts");
 
        time = timeIntegrate.getElapsedTime();
        maxDist += integrate1();
        timeInt1 += timeIntegrate.getElapsedTime() - time;

	LOG4ESPP_INFO(theLogger, "maxDist = " << maxDist << ", skin/2 = " << skinHalf);

	if (maxDist > skinHalf) resortFlag = true;

	if (resortFlag) {
          time = timeIntegrate.getElapsedTime();
	  LOG4ESPP_INFO(theLogger, "step " << i << ": resort particles");
	  storage.decompose();
	  maxDist  = 0.0;
	  resortFlag = false;
          nResorts ++;
          timeResort += timeIntegrate.getElapsedTime() - time;
	}

        updateForces();

        time = timeIntegrate.getElapsedTime();
        integrate2();
        timeInt2 += timeIntegrate.getElapsedTime() - time;

      }

      LOG4ESPP_INFO(theLogger, "finished run");

      // ToDo: print Timers only if INFO is enabled

      printTimers();
    }

    void VelocityVerletOnGroup::resetTimers()
    {
      timeResort = 0.0;
      timeForce  = 0.0;
      for(int i = 0; i < 100; i++)
        timeForceComp[i] = 0.0;
      timeComm1  = 0.0;
      timeComm2  = 0.0;
      timeInt1   = 0.0;
      timeInt2   = 0.0;
    }

    void VelocityVerletOnGroup::printTimers()
    {
      std::cout << "time: run = " << timeIntegrate <<
                   ", pair = " << timeForceComp[0] <<
                   ", FENE = " << timeForceComp[1] <<
                   ", angle = " << timeForceComp[2] <<
                   ", comm1 = " << timeComm1 <<
                   ", comm2 = " << timeComm2 <<
                   ", int1 = " << timeInt1 <<
                   ", int2 = " << timeInt2 <<
                   ", resort = " << timeResort << std::endl; 
    }

    /*****************************************************************************/

    real VelocityVerletOnGroup::integrate1()
    {
      System& system = getSystemRef();

      // loop over all particles of the local cells

      int count = 0;

      real maxSqDist = 0.0; // maximal square distance a particle moves

      for(ParticleGroup::iterator cit(group->begin()); cit != group->end(); ++cit) {

	real sqDist = 0.0;

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
        Real3D deltaP = dt * cit->velocity();
        cit->position() += deltaP;
        sqDist += deltaP * deltaP;

	count++;

	maxSqDist = std::max(maxSqDist, sqDist);
      }

      real maxAllSqDist;

      mpi::all_reduce(*system.comm, maxSqDist, maxAllSqDist, 
                      boost::mpi::maximum<real>());

      LOG4ESPP_INFO(theLogger, "moved " << count << " particles in integrate1" <<
		    ", max move local = " << sqrt(maxSqDist) <<
		    ", global = " << sqrt(maxAllSqDist));

      return sqrt(maxAllSqDist);
    }

    /*****************************************************************************/

    void VelocityVerletOnGroup::integrate2()
    {
      System& system = getSystemRef();

      // loop over all particles of the local cells

      for(ParticleGroup::iterator cit(group->begin()); cit != group->end(); ++cit) {

        real dtfm = 0.5 * dt / cit->mass();

	/* Propagate velocities: v(t+0.5*dt) = v(t) + 0.5*dt * f(t) */

	cit->velocity() += dtfm * cit->force();
      }
    }

    /*****************************************************************************/

    void VelocityVerletOnGroup::setUp()
    {
      System& system = getSystemRef();

      const InteractionList& srIL = system.shortRangeInteractions;

      maxCut = 0.0;

      for (size_t j = 0; j < srIL.size(); j++) {

	real cut = srIL[j]->getMaxCutoff();

	maxCut = std::max(maxCut, cut);
      }

      LOG4ESPP_INFO(theLogger, "maximal cutoff = " << maxCut);
    }

    /*****************************************************************************/

    void VelocityVerletOnGroup::calcForces()
    {
      LOG4ESPP_INFO(theLogger, "calculate forces");

      initForces();

      System& sys = getSystemRef();

      const InteractionList& srIL = sys.shortRangeInteractions;

      for (size_t i = 0; i < srIL.size(); i++) {

	LOG4ESPP_INFO(theLogger, "compute forces for srIL " << i 
                                  << " of " << srIL.size());

        real time;
        time = timeIntegrate.getElapsedTime();
	srIL[i]->addForces();
        timeForceComp[i] += timeIntegrate.getElapsedTime() - time;
      }
    }

    void VelocityVerletOnGroup::updateForces()
    {
      real time;
      storage::Storage& storage = *getSystemRef().storage;

      time = timeIntegrate.getElapsedTime();
      storage.updateGhosts();
      timeComm1 += timeIntegrate.getElapsedTime() - time;

      time = timeIntegrate.getElapsedTime();
      calcForces();
      timeForce += timeIntegrate.getElapsedTime() - time;

      time = timeIntegrate.getElapsedTime();
      storage.collectGhostForces();
      timeComm2 += timeIntegrate.getElapsedTime() - time;

      if (langevin) langevin->thermalize();
    }

    /*****************************************************************************/

    void VelocityVerletOnGroup::initForces()
    {
      // forces are initialized for real + ghost particles

      System& system = getSystemRef();

      CellList localCells = system.storage->getLocalCells();

      LOG4ESPP_INFO(theLogger, "init forces for real + ghost particles");

      for(CellListIterator cit(localCells); !cit.isDone(); ++cit) {
        cit->force() = 0.0;
      }

    }

    /*****************************************************************************/

    void VelocityVerletOnGroup::printForces(bool withGhosts)
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

        LOG4ESPP_DEBUG(theLogger,
                       "Particle " << cit->id()
                       << ", force = " << cit->force());
      }
    }

    /*****************************************************************************/

    void VelocityVerletOnGroup::printPositions(bool withGhosts)
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
  
      for(ParticleGroup::iterator cit(group->begin()); cit != group->end(); ++cit) {

	LOG4ESPP_DEBUG(theLogger, 
		       "Particle " << cit->id()
		       << ", position = " << cit->position());
      }
    }

    /****************************************************
    ** REGISTRATION WITH PYTHON
    ****************************************************/

    void VelocityVerletOnGroup::registerPython() {

      using namespace espressopp::python;

      // Note: use noncopyable and no_init for abstract classes

      class_<VelocityVerletOnGroup, bases<MDIntegrator>, boost::noncopyable >
        ("integrator_VelocityVerletOnGroup", init< shared_ptr<System>,  shared_ptr<ParticleGroup> >())

        .add_property("langevin", &VelocityVerletOnGroup::getLangevin, &VelocityVerletOnGroup::setLangevin)
        ;
    }
  }
}
