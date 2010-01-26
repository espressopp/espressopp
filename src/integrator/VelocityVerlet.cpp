#include "VelocityVerlet.hpp"

#include "VerletList.hpp"
#include "interaction/LennardJones.hpp"
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

  if (langevin) langevin->init(dt);

  bool recalcForces = true;  // TODO: more intelligent

  if (recalcForces) {

     if (langevin) langevin->heatUp();

     calcForces();

     if (langevin) langevin->coolDown();

  }

  LOG4ESPP_INFO(theLogger, "run " << nsteps << " iterations");
  
  // Before start make sure that particles are on the right processor

  system.lock().get()->storage->resortParticles();

  for (int i = 0; i < nsteps; i++) {
    LOG4ESPP_DEBUG(theLogger, "step " << i << " of " << nsteps << " iterations");
    integrate1();
    if (rebuild()) {
      system.lock().get()->storage->resortParticles();
    }
    system.lock().get()->storage->updateGhosts();
    calcForces();
    system.lock().get()->storage->collectGhostForces();

    if (langevin) langevin->thermalize();
   
    integrate2();
  }
}

/*****************************************************************************/

void VelocityVerlet::integrate1()
{
  std::vector<Cell>& localCells = system.lock().get()->storage->getLocalCells();

  // loop over all particles of the local cells

  real dt2 = dt * dt;
  int count = 0;

  for (size_t c = 0; c < localCells.size(); c++) {
    Cell* localCell = &localCells[c];
    for (size_t index = 0; index < localCell->particles.size(); index++) {
      Particle* particle  = &localCell->particles[index];
      for (int j = 0; j < 3; j++) {
        /* Propagate velocities: v(t+0.5*dt) = v(t) + 0.5*dt * f(t) */
        particle->m.v[j] += 0.5 * dt * particle->f.f[j];
        /* Propagate positions (only NVT): p(t + dt)   = p(t) + dt * v(t+0.5*dt) */
        particle->r.p[j] += 0.5 * dt2 * particle->m.v[j];
        count++;
      }
    }
  }

  LOG4ESPP_DEBUG(theLogger, "moved " << count << " particles in integrate1");
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

  real cut = 1.5;

  shared_ptr< VerletList > vl = 
    make_shared<VerletList>(system.lock(), cut);

  VerletList::PairList pairs = vl->getPairs();
  LOG4ESPP_DEBUG(theLogger, "# of verlet list pairs =  " << pairs.size());
  LennardJones lj = LennardJones();
  lj.setParameters(0, 0, 1.0, 1.0, 1.3);
  lj.addVerletListForces(vl);

  // Just for control now: compute energy
  real e2 = lj.computeVerletListEnergy(vl);
  LOG4ESPP_INFO(theLogger, "energy  = " << e2);
}

/*****************************************************************************/

bool VelocityVerlet::rebuild()
{
  // check for maximal movement
  return true;
}

