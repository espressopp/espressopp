#include "VelocityVerlet.hpp"

#include "VerletList.hpp"
#include "Interaction.hpp"
#include "Langevin.hpp"
#include "System.hpp"

using namespace espresso;
using namespace integrator;
using namespace interaction;

VelocityVerlet::VelocityVerlet(shared_ptr< System > system) : MDIntegrator(system)
{
  LOG4ESPP_INFO(theLogger, "construct VelocityVerlet");

  resortFlag = true;
}

VelocityVerlet::~VelocityVerlet()
{
  LOG4ESPP_INFO(theLogger, "free VelocityVerlet");
}

/*****************************************************************************/

void VelocityVerlet::setLangevin(shared_ptr<Langevin> _langevin)
{
  langevin = _langevin;
}

/*****************************************************************************/

void VelocityVerlet::run(int nsteps)
{

  System* pSystem = system.lock().get();

  if (langevin) langevin->init(dt);

  setUp();

  // Before start make sure that particles are on the right processor

  real maxSqDist;

  if (resortFlag) {
     LOG4ESPP_INFO(theLogger, "resort particles + rebuild VL");
     pSystem->storage->resortParticles();
     vl = make_shared<VerletList>(system.lock(), maxCut + pSystem->skin);
     maxSqDist = 0.0;
     resortFlag = false;
  }

  bool recalcForces = true;  // TODO: more intelligent

  if (recalcForces) {

     if (langevin) langevin->heatUp();

     pSystem->storage->updateGhosts();
     calcForces();
     pSystem->storage->collectGhostForces();

     if (langevin) langevin->coolDown();

  }

  LOG4ESPP_INFO(theLogger, "run " << nsteps << " iterations");
  
  real skinHalfSq = 0.25 * (pSystem->skin * pSystem->skin);

  for (int i = 0; i < nsteps; i++) {

    LOG4ESPP_INFO(theLogger, "step " << i << " of " << nsteps << " iterations");

    maxSqDist += integrate1();

    LOG4ESPP_INFO(theLogger, "maxSqDist = " << maxSqDist);

    if (maxSqDist > skinHalfSq) resortFlag = true;

    if (resortFlag) {
      LOG4ESPP_INFO(theLogger, "resort particles + rebuild VL");
      pSystem->storage->resortParticles();
      vl = make_shared<VerletList>(system.lock(), maxCut + pSystem->skin);
      maxSqDist  = 0.0;
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
  std::vector<Cell*>& realCells = system.lock().get()->storage->getRealCells();

  // loop over all particles of the local cells

  int count = 0;

  real maxSqDist = 0.0; // maximal square distance a particle moves

  for (size_t c = 0; c < realCells.size(); c++) {
    Cell* localCell = realCells[c];
    for (size_t index = 0; index < localCell->particles.size(); index++) {
      Particle& pt  = localCell->particles[index];
      real sqDist = 0.0;

      printf("Particle %d, pos = %f %f %f, vel = %f %f %f, f = %f %f %f\n", 
              pt.p.id, pt.r.p[0], pt.r.p[1], pt.r.p[2], 
              pt.m.v[0], pt.m.v[1], pt.m.v[2],
              pt.f.f[0], pt.f.f[1], pt.f.f[2]);

      for (int j = 0; j < 3; j++) {
        /* Propagate velocities: v(t+0.5*dt) = v(t) + 0.5*dt * f(t) */
        pt.m.v[j] += 0.5 * dt * pt.f.f[j];
        /* Propagate positions (only NVT): p(t + dt)   = p(t) + dt * v(t+0.5*dt) */
        real deltaP = dt * pt.m.v[j];
        pt.r.p[j] += deltaP;
        sqDist += deltaP * deltaP;
      }
      count++;
      maxSqDist = std::max(maxSqDist, sqDist);
    }
  }

  LOG4ESPP_DEBUG(theLogger, "moved " << count << " particles in integrate1" <<
                            ", max sqr move = " << maxSqDist);
  return maxSqDist;
}

/*****************************************************************************/

void VelocityVerlet::integrate2()
{
  std::vector<Cell*>& realCells = system.lock().get()->storage->getRealCells();

  // loop over all particles of the local cells

  for (size_t c = 0; c < realCells.size(); c++) {
    Cell* realCell = realCells[c];
    for (size_t index = 0; index < realCell->particles.size(); index++) {
      Particle& pt  = realCell->particles[index];
      printf("Particle %d, pos = %f %f %f, vel = %f %f %f, f = %f %f %f\n", 
              pt.p.id, pt.r.p[0], pt.r.p[1], pt.r.p[2], 
              pt.m.v[0], pt.m.v[1], pt.m.v[2],
              pt.f.f[0], pt.f.f[1], pt.f.f[2]);
      for (int j = 0; j < 3; j++) {
        /* Propagate velocities: v(t+0.5*dt) = v(t) + 0.5*dt * f(t) */
        pt.m.v[j] += 0.5 * dt * pt.f.f[j];
      }
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

     energy += srIL[i]->computeVerletListEnergy(vl);
  }

  // Just for control now: compute + print energy

  LOG4ESPP_INFO(theLogger, "energy  = " << energy);
}

/*****************************************************************************/

void VelocityVerlet::initForces()
{
  // forces are initialized for real + ghost particles

  std::vector<Cell>& localCells = system.lock().get()->storage->getLocalCells();

  for (size_t c = 0; c < localCells.size(); c++) {
    Cell* localCell = &localCells[c];
    for (size_t index = 0; index < localCell->particles.size(); index++) {
      Particle* particle  = &localCell->particles[index];
      for (int j = 0; j < 3; j++) {
        particle->f.f[j] = 0.0;
      }
    }
  }
}

