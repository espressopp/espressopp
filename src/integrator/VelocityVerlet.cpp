#include "VelocityVerlet.hpp"

#include "Langevin.hpp"

#include "iterator/CellListIterator.hpp"
#include "interaction/Interaction.hpp"
#include "System.hpp"
#include "storage/Storage.hpp"
#include "mpi.hpp"

using namespace espresso;
using namespace interaction;
using namespace iterator;

namespace espresso {
  namespace integrator {

    VelocityVerlet::VelocityVerlet(shared_ptr< System > system) : MDIntegrator(system)
    {
      LOG4ESPP_INFO(theLogger, "construct VelocityVerlet");

      resortFlag = true;
    }

    /*****************************************************************************/

    VelocityVerlet::~VelocityVerlet()
    {
      LOG4ESPP_INFO(theLogger, "free VelocityVerlet");
    }

    /*****************************************************************************/

    void VelocityVerlet::setLangevin(shared_ptr< Langevin > _langevin)
    {
      LOG4ESPP_INFO(theLogger, "set Langevin thermostat");
      langevin = _langevin;
    }

    /*****************************************************************************/

    void VelocityVerlet::run(int nsteps)
    {
      storage::Storage &storage = *(getSystem()->storage);

      if (langevin) langevin->init(dt);

      // no more needed: setUp();

      // Before start make sure that particles are on the right processor

      real maxDist;

      if (resortFlag) {
	LOG4ESPP_INFO(theLogger, "resort particles");
	storage.resortParticles();
	LOG4ESPP_INFO(theLogger, "particles resort");
	maxDist = 0.0;
	resortFlag = false;
      }

      bool recalcForces = true;  // TODO: more intelligent

      if (recalcForces) {

	LOG4ESPP_INFO(theLogger, "recalc Forces");

	if (langevin) langevin->heatUp();

	if (LOG4ESPP_DEBUG_ON(theLogger)) {
	  printPositions(false);
	}

	storage.updateGhosts();

	if (LOG4ESPP_DEBUG_ON(theLogger)) {
	  printPositions(true);
	}

	calcForces();

	if (LOG4ESPP_DEBUG_ON(theLogger)) {
	  printForces(true);    // check forces in real + ghost particles
	}

	storage.collectGhostForces();

	if (LOG4ESPP_DEBUG_ON(theLogger)) {
	  printForces(false);   // forces are reduced to real particles
	}

	storage.collectGhostForces();

	if (langevin) langevin->coolDown();
      }

      LOG4ESPP_INFO(theLogger, "run " << nsteps << " iterations");
  
      real skinHalf = 0.5 * getSystem()->skin;

      for (int i = 0; i < nsteps; i++) {

	LOG4ESPP_INFO(theLogger, "Next step " << i << " of " << nsteps << " starts");

	maxDist += integrate1();

	LOG4ESPP_INFO(theLogger, "maxDist = " << maxDist);

	if (maxDist > skinHalf) resortFlag = true;

	if (resortFlag) {
	  LOG4ESPP_INFO(theLogger, "step " << i << ": resort particles");
	  storage.resortParticles();
	  maxDist  = 0.0;
	  resortFlag = false;
	}

	storage.updateGhosts();

	calcForces();

	storage.collectGhostForces();

	if (langevin) langevin->thermalize();
   
	integrate2();
      }

      LOG4ESPP_INFO(theLogger, "proc " << mpiWorld.rank() << ": finished run");
    }

    /*****************************************************************************/

    real VelocityVerlet::integrate1()
    {
      CellList realCells = getSystem()->storage->getRealCells();

      // loop over all particles of the local cells

      int count = 0;

      real maxSqDist = 0.0; // maximal square distance a particle moves

      for(CellListIterator cit(realCells); !cit.isDone(); ++cit) {

	real sqDist = 0.0;

	LOG4ESPP_DEBUG(theLogger, "Particle " << cit->p.id << 
		       ", pos = " << cit->r.p[0] << " " 
		       << cit->r.p[1] << " " <<  cit->r.p[2] <<
		       ", v = " << cit->m.v[0] << " " 
		       << cit->m.v[1] << " " <<  cit->m.v[2] <<
		       ", f = " << cit->f.f[0] << " " 
		       << cit->f.f[1] << " " <<  cit->f.f[2]);

	for (int j = 0; j < 3; j++) {
	  /* Propagate velocities: v(t+0.5*dt) = v(t) + 0.5*dt * f(t) */
	  cit->m.v[j] += 0.5 * dt * cit->f.f[j];
	  /* Propagate positions (only NVT): p(t + dt)   = p(t) + dt * v(t+0.5*dt) */
	  real deltaP = dt * cit->m.v[j];
	  cit->r.p[j] += deltaP;
	  sqDist += deltaP * deltaP;
	}

	count++;

	maxSqDist = std::max(maxSqDist, sqDist);
      }

      // ToDo: here or outside: mpi::all_reduce(maxval)

      real maxAllSqDist;

      mpi::all_reduce(mpiWorld, maxSqDist, maxAllSqDist, boost::mpi::maximum<double>());

      LOG4ESPP_INFO(theLogger, "moved " << count << " particles in integrate1" <<
		    ", max move local = " << sqrt(maxSqDist) <<
		    ", global = " << sqrt(maxAllSqDist));

      return sqrt(maxAllSqDist);
    }

    /*****************************************************************************/

    void VelocityVerlet::integrate2()
    {
      CellList realCells = getSystem()->storage->getRealCells();

      // loop over all particles of the local cells

      for(CellListIterator cit(realCells); !cit.isDone(); ++cit) {

	for (int j = 0; j < 3; j++) {

	  /* Propagate velocities: v(t+0.5*dt) = v(t) + 0.5*dt * f(t) */
	  cit->m.v[j] += 0.5 * dt * cit->f.f[j];
	}
      }
    }

    /*****************************************************************************/

    void VelocityVerlet::setUp()
    {
      System& system = getSystemRef();

      const InteractionList& srIL = system.shortRangeInteractions;

      maxCut = 0.0;

      for (int j = 0; j < srIL.size(); j++) {

	real cut = srIL[j]->getMaxCutoff();

	maxCut = std::max(maxCut, cut);
      }

      LOG4ESPP_INFO(theLogger, "maximal cutoff = " << maxCut);
    }

    /*****************************************************************************/

    void VelocityVerlet::calcForces()
    {
      LOG4ESPP_INFO(theLogger, "calculate forces");

      initForces();

      System& sys = getSystemRef();

      const InteractionList& srIL = sys.shortRangeInteractions;

      real energy = 0.0;

      for (size_t i = 0; i < srIL.size(); i++) {

	LOG4ESPP_INFO(theLogger, "compute forces for srIL " << i << " of " << srIL.size());

	srIL[i]->addForces();
      }

      // Just for control now: compute + print energy

      // LOG4ESPP_INFO(theLogger, "energy  = " << energy);
    }

    /*****************************************************************************/

    void VelocityVerlet::initForces()
    {
      // forces are initialized for real + ghost particles

      // ToDo: make one loop when getLocalCells() works

      CellList localCells = getSystem()->storage->getLocalCells();

      LOG4ESPP_INFO(theLogger, "init forces for real + ghost particles");

      for(CellListIterator cit(localCells); !cit.isDone(); ++cit) {
	for (int j = 0; j < 3; j++) {
	  cit->f.f[j] = 0.0;
	}
      }
    }

    /*****************************************************************************/

    void VelocityVerlet::printForces(bool withGhosts)
    {
      // print forces of real + ghost particles

      CellList cells;

      if (withGhosts) {
	cells = getSystem()->storage->getLocalCells();
	LOG4ESPP_DEBUG(theLogger, "Proc " << mpiWorld.rank() << ": local forces");
      } else {
	cells = getSystem()->storage->getRealCells();
	LOG4ESPP_DEBUG(theLogger, "Proc " << mpiWorld.rank() << ": real forces");
      }
  
      for(CellListIterator cit(cells); !cit.isDone(); ++cit) {

	LOG4ESPP_DEBUG(theLogger, "Proc " << mpiWorld.rank()
		       << ": Particle " << cit->p.id 
		       << ", force = " << cit->f.f[0] << " "
		       << cit->f.f[1] << " " <<  cit->f.f[2]);
      }
    }

    /*****************************************************************************/

    void VelocityVerlet::printPositions(bool withGhosts)
    {
      // print positions of real + ghost particles

      CellList cells;

      if (withGhosts) {
	cells = getSystem()->storage->getLocalCells();
	LOG4ESPP_DEBUG(theLogger, "Proc " << mpiWorld.rank() << ": local positions");
      } else {
	cells = getSystem()->storage->getRealCells();
	LOG4ESPP_DEBUG(theLogger, "Proc " << mpiWorld.rank() << ": real positions");
      }
  
      for(CellListIterator cit(cells); !cit.isDone(); ++cit) {

	LOG4ESPP_DEBUG(theLogger, "Proc " << mpiWorld.rank()
		       << ": Particle " << cit->p.id
		       << ", position = " << cit->r.p[0] << " "
		       << cit->r.p[1] << " " <<  cit->r.p[2]);
      }
    }
  }
}
