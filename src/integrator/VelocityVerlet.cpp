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

  bool recalcForces = true;  // TODO: more intelligent

  if (recalcForces) {

     if (langevin) langevin->heatUp();

     calcForces();

     if (langevin) langevin->coolDown();

  }

  LOG4ESPP_INFO(theLogger, "run " << nsteps << " iterations");
  
  // Before start make sure that particles are on the right processor

  pSystem->storage->resortParticles();

  real maxSqDist = 0.0;
  real skinHalfSq = 0.25 * (pSystem->skin * pSystem->skin);

  for (int i = 0; i < nsteps; i++) {
    LOG4ESPP_DEBUG(theLogger, "step " << i << " of " << nsteps << " iterations");
    maxSqDist += integrate1();
    if (maxSqDist > skinHalfSq) {
      LOG4ESPP_DEBUG(theLogger, "resort particles");
      pSystem->storage->resortParticles();
      maxSqDist = 0.0;
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
  std::vector<Cell>& localCells = system.lock().get()->storage->getLocalCells();

  // loop over all particles of the local cells

  real dt2 = dt * dt;
  int count = 0;

  real maxSqDist = 0.0; // maximal square distance a particle moves

  for (size_t c = 0; c < localCells.size(); c++) {
    Cell* localCell = &localCells[c];
    for (size_t index = 0; index < localCell->particles.size(); index++) {
      Particle* particle  = &localCell->particles[index];
      real sqDist = 0.0;
      for (int j = 0; j < 3; j++) {
        /* Propagate velocities: v(t+0.5*dt) = v(t) + 0.5*dt * f(t) */
        particle->m.v[j] += 0.5 * dt * particle->f.f[j];
        /* Propagate positions (only NVT): p(t + dt)   = p(t) + dt * v(t+0.5*dt) */
        real deltaP = 0.5 * dt2 * particle->m.v[j];
        particle->r.p[j] += deltaP;
        sqDist += deltaP * deltaP;
        count++;
      }
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
  std::vector<Cell>& localCells = system.lock().get()->storage->getLocalCells();

  // loop over all particles of the local cells

  for (size_t c = 0; c < localCells.size(); c++) {
    Cell* localCell = &localCells[c];
    for (size_t index = 0; index < localCell->particles.size(); index++) {
      Particle* particle  = &localCell->particles[index];
      for (int j = 0; j < 3; j++) {
        /* Propagate velocities: v(t+0.5*dt) = v(t) + 0.5*dt * f(t) */
        particle->m.v[j] += 0.5 * dt * particle->f.f[j];
      }
    }
  }
}

/*****************************************************************************/

void VelocityVerlet::calcForces()
{
  // build VerletList (currently each step, issue for optimization)

  System& sys = *(system.lock().get());

  InteractionList& interactions = sys.shortRangeInteractions;

  real cutForce = 0.0;

  for (int j = 0; j < interactions.size(); j++) {

     real cut = interactions[j]->getMaxCutoff();

     cutForce = std::max(cutForce, cut);
  }

  printf("maximal cut for all forces: %f\n", cutForce);

  shared_ptr< VerletList > vl = make_shared<VerletList>(system.lock(), 
                                      cutForce + sys.skin);

  VerletList::PairList pairs = vl->getPairs();

  LOG4ESPP_DEBUG(theLogger, "# of verlet list pairs =  " << pairs.size());

  real energy = 0.0;

  const InteractionList& srIL = sys.shortRangeInteractions;

  for (size_t i = 0; i < srIL.size(); i++) {

     srIL[i]->addVerletListForces(vl);
     // energy += srIL[i]->computeVerletListEnergy(vl);
  }

  // Just for control now: compute + print energy

  LOG4ESPP_INFO(theLogger, "energy  = " << energy);
}

