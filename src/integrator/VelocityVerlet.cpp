#include "VelocityVerlet.hpp"

#include "VerletList.hpp"
#include "iterator/CellListIterator.hpp"
#include "Interaction.hpp"
#include "Langevin.hpp"
#include "System.hpp"

using namespace espresso;
using namespace integrator;
using namespace interaction;
using namespace iterator;

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

void VelocityVerlet::setLangevin(shared_ptr<Langevin> _langevin)
{
  LOG4ESPP_INFO(theLogger, "set Langevin thermostat");
  langevin = _langevin;
}

/*****************************************************************************/

void VelocityVerlet::run(int nsteps)
{
  System* pSystem = system.lock().get();

  if (langevin) langevin->init(dt);

  setUp();

  // Before start make sure that particles are on the right processor

  real maxDist;

  if (resortFlag) {
     LOG4ESPP_INFO(theLogger, "resort particles");
     pSystem->storage->resortParticles();
     LOG4ESPP_INFO(theLogger, "particles rosort, build VL");
     vl = make_shared<VerletList>(system.lock(), maxCut + pSystem->skin);
     maxDist = 0.0;
     resortFlag = false;
  }

  bool recalcForces = true;  // TODO: more intelligent

  if (recalcForces) {

     LOG4ESPP_INFO(theLogger, "recalc Forces");

     if (langevin) langevin->heatUp();

     pSystem->storage->updateGhosts();
     calcForces();
     pSystem->storage->collectGhostForces();

     if (langevin) langevin->coolDown();

  }

  LOG4ESPP_INFO(theLogger, "run " << nsteps << " iterations");
  
  real skinHalf = 0.5 * pSystem->skin;

  for (int i = 0; i < nsteps; i++) {

    LOG4ESPP_INFO(theLogger, "step " << i << " of " << nsteps << " iterations");

    maxDist += integrate1();

    LOG4ESPP_INFO(theLogger, "maxDist = " << maxDist);

    if (maxDist > skinHalf) resortFlag = true;

    if (resortFlag) {
      LOG4ESPP_INFO(theLogger, "resort particles + rebuild VL");
      pSystem->storage->resortParticles();
      vl = make_shared<VerletList>(system.lock(), maxCut + pSystem->skin);
      maxDist  = 0.0;
      resortFlag = false;
    }

    pSystem->storage->updateGhosts();
    calcForces();
    pSystem->storage->collectGhostForces();

    if (langevin) langevin->thermalize();
   
    integrate2();
  }
}

/*****************************************************************************/

real VelocityVerlet::integrate1()
{
  CellList realCells = system.lock().get()->storage->getRealCells();

  // loop over all particles of the local cells

  int count = 0;

  real maxSqDist = 0.0; // maximal square distance a particle moves

  for(CellListIterator cit(realCells); !cit.isDone(); ++cit) {

    real sqDist = 0.0;

    /*
    printf("Particle %d, pos = %f %f %f, vel = %f %f %f, f = %f %f %f\n", 
            cit->p.id, cit->r.p[0], cit->r.p[1], cit->r.p[2], 
            cit->m.v[0], cit->m.v[1], cit->m.v[2],
            cit->f.f[0], cit->f.f[1], cit->f.f[2]);
    */

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

  // ToDo: here or outside: mpi->sum 

  LOG4ESPP_INFO(theLogger, "moved " << count << " particles in integrate1" <<
                            ", max sqr move = " << maxSqDist);

  return sqrt(maxSqDist);
}

/*****************************************************************************/

void VelocityVerlet::integrate2()
{
  CellList realCells = system.lock().get()->storage->getRealCells();

  // loop over all particles of the local cells

  for(CellListIterator cit(realCells); !cit.isDone(); ++cit) {

    /*
    printf("Particle %d, pos = %f %f %f, vel = %f %f %f, f = %f %f %f\n", 
            cit->p.id, cit->r.p[0], cit->r.p[1], cit->r.p[2], 
            cit->m.v[0], cit->m.v[1], cit->m.v[2],
            cit->f.f[0], cit->f.f[1], cit->f.f[2]);
    */

    for (int j = 0; j < 3; j++) {
      /* Propagate velocities: v(t+0.5*dt) = v(t) + 0.5*dt * f(t) */
      cit->m.v[j] += 0.5 * dt * cit->f.f[j];
    }
  }
}

/*****************************************************************************/

void VelocityVerlet::setUp()
{
  System* pSystem = system.lock().get();

  const InteractionList& srIL = pSystem->shortRangeInteractions;

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
  initForces();

  System& sys = *(system.lock().get());

  const InteractionList& srIL = sys.shortRangeInteractions;

  real energy = 0.0;

  for (size_t i = 0; i < srIL.size(); i++) {

     LOG4ESPP_INFO(theLogger, "compute forces for srIL " << i << " of " << srIL.size());

     srIL[i]->addVerletListForces(vl);

     // energy += srIL[i]->computeVerletListEnergy(vl);
  }

  // Just for control now: compute + print energy

  // LOG4ESPP_INFO(theLogger, "energy  = " << energy);
}

/*****************************************************************************/

void VelocityVerlet::initForces()
{
  // forces are initialized for real + ghost particles

  // ToDo: make one loop when getLocalCells() works

  CellList realCells = system.lock().get()->storage->getRealCells();

  LOG4ESPP_INFO(theLogger, "init forces for real + ghost particles");

  for(CellListIterator cit(realCells); !cit.isDone(); ++cit) {
    for (int j = 0; j < 3; j++) {
      cit->f.f[j] = 0.0;
    }
  }

  CellList ghostCells = system.lock().get()->storage->getGhostCells();

  for(CellListIterator cit(ghostCells); !cit.isDone(); ++cit) {
    for (int j = 0; j < 3; j++) {
      cit->f.f[j] = 0.0;
    }
  }
}

