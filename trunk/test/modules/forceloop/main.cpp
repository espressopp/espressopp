
// Example C++ program doing the same as the force loop Python script

#include <types.hpp>
#include <logging.hpp>
#include <pmi.hpp>
#include <espresso_common.hpp>
#include <particles/Storage.hpp>
#include <particles/All.hpp>
#include <bc/PBC.hpp>
#include <interaction/LennardJones.hpp>
#include <interaction/FENE.hpp>
#include <pairs/All.hpp>
#include <pairs/ForceComputer.hpp>

#include "ParticleWriter.hpp"
#include "PairWriteComputer.hpp"

#include "integrator/VelocityVerlet.hpp"

#include <cstdio>
#include <vector>

using namespace espresso;
using namespace espresso::bc;
using namespace espresso::particles;
using namespace espresso::pairs;
using namespace espresso::interaction;
using namespace espresso::integrator;

/** N stands for number particles in each dimensions.

    N   SIZE 
  10    15.0

*/

#define N 3
#define SIZE 4.0

/** Main routine of a test program:

    - generate N * N * N particle in a box of size SIZE * SIZE * SIZE 
    - print out particle data
    - define periodic boundary conditions
    - define Lennard Jones interaction and apply it to all pairs
    - print out particle data
*/

void forceloop() {
  // Create a new particle storage

  Storage particleStorage;
  PropertyId position = particleStorage.addProperty<Real3D>();
  PropertyId velocity = particleStorage.addProperty<Real3D>();
  PropertyId force = particleStorage.addProperty<Real3D>();

  // generate particles in the particle storage

  for (int i = 0; i < N; i++) 
  for (int j = 0; j < N; j++) 
  for (int k = 0; k < N; k++) {
      
       real r;
       r = 0.4 + 0.2 * rand() / RAND_MAX;
       real x = (i + r) / N * SIZE;
       real y = (j + r) / N * SIZE; 
       real z = (k + r) / N * SIZE;

       ParticleReference ref = particleStorage.addParticle();
       PropertyReference<Real3D> positionRef = 
	 particleStorage.getPropertyReference<Real3D>(position);
       PropertyReference<Real3D> velocityRef = 
	 particleStorage.getPropertyReference<Real3D>(velocity);
       PropertyReference<Real3D> forceRef    = 
	 particleStorage.getPropertyReference<Real3D>(force);

       positionRef[ref] = Real3D(x, y, z);
       velocityRef[ref] = Real3D(x, y, z);
       forceRef[ref] = 0.0;
  }

  // For test only: ParticleWriter prints each particle

  ParticleWriter pWriter(particleStorage, position, force);

  // call pWriter(ref) for each particle reference ref of particle storage

  particleStorage.foreach(pWriter);

  // define periodic boundary conditions

  PBC pbc(SIZE);

  // define a set of all particles

  particles::All allSet(&particleStorage);

  // define allpairs with (x, y) for all x, y in allSet

  pairs::All allpairs(pbc, allSet, position);

  // For test only: PairWriter prints each particle pair 

  PairWriteComputer pairWriter(&particleStorage, position);

  // call pairWriter(ref1, ref2) for each particle ref pair of allSet

  // allpairs.foreach(pairWriter);

  // define LennardJones interaction

  LennardJones ljint;

  ljint.set(1.0, 1.0, 2.5);

  // make a FENE interaction

  FENE fene;

  fene.set(1.5, 1.0, 2.5);

  // force will be the vector of all forces in the particle storage
  // and force[ref] returns the force (as RealArrayRef) of particle reference ref

  PropertyReference<Real3D> forceRef = particleStorage.getPropertyReference<Real3D>(force);

  // Define a pair computer that computes the forces for particle pairs
  // ljint provides the routine computeForce for a particle pair
  // force (pointer to all forces of particles) tells us where the computed forces are added

  ForceComputer *forcecompute = ljint.createForceComputer(ForceComputer(forceRef));

  // call forcecompute(ref1, ref2) for each particle ref pair of allSet

  allpairs.foreach(*forcecompute);

  delete forcecompute;

  // print out all particle data to see that it calculates some forces

  particleStorage.foreach(pWriter);

  // create integrator and set properties

  VelocityVerlet integrator(&allSet, position, velocity, force);
  integrator.setTimeStep(0.005);
  integrator.addForce(&ljint, &allpairs);

  // integrate for a certain number of timesteps

  integrator.run(100);

  // check to see that particles have new positions

  particleStorage.foreach(pWriter);

}

int main() 

{

  LOG4ESPP_CONFIGURE();

  // Initialization of PMI and Python

  // initPythonEspresso();

#ifdef HAVE_MPI
  initMPI();

  std::cout << "Worker " << pmi::getWorkerId() << std::endl ;

  if (pmi::isController()) {

      std::cout << "Controller starts forceloop" << std::endl;

      forceloop();

      std::cout << "Controller ends forceloop" << std::endl;

      pmi::endWorkers();

      std::cout << "Controller has stopped workers" << std::endl;

  } else {

      std::cout << "Worker " << pmi::getWorkerId() << " starts mainLoop" << std::endl ;

      pmi::mainLoop();

      std::cout << "Worker " << pmi::getWorkerId() << " ends mainLoop" << std::endl ;
  }
  finalizeMPI();
#else
  forceloop();
#endif

}
