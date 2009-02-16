
// Example C++ program doing the same as the force loop Python script

#include "types.hpp"
#include "bc/PBC.hpp"
#include "interaction/LennardJones.hpp"
#include "interaction/FENE.hpp"
#include "pairs/All.hpp"
#include "pairs/PairForceComputer.hpp"
#include "particleset/All.hpp"
#include "particlestorage/ParticleStorage.hpp"

#include "ParticleWriter.hpp"
#include "PairWriteComputer.hpp"

#include "integrator/VelocityVerlet.hpp"

#include <cstdio>
#include <vector>

using namespace std;
using namespace espresso;
using namespace espresso::bc;
using namespace espresso::pairs;
using namespace espresso::particleset;
using namespace espresso::particlestorage;
using namespace espresso::interaction;
using namespace espresso::integrator;

/** N stands for number particles in each dimensions.
*/

#define N 10
#define SIZE 15.0

typedef espresso::pairs::All AllPairs;
typedef espresso::particleset::All AllParticles;

/** Main routine of a test program:

    - generate N * N * N particle in a box of size SIZE * SIZE * SIZE 
    - print out particle data
    - define periodic boundary conditions
    - define Lennard Jones interaction and apply it to all pairs
    - print out particle data
*/

int main() 

{
  // Create a new particle storage

  ParticleStorage particleStorage;
  size_t position = particleStorage.addProperty<Real3D>();
  size_t velocity = particleStorage.addProperty<Real3D>();
  size_t force = particleStorage.addProperty<Real3D>();

  // generate particles in the particle storage

  for (int i = 0; i < N; i++) 
  for (int j = 0; j < N; j++) 
  for (int k = 0; k < N; k++) {
      
       real r;
       r = 0.4 + 0.2 * rand() / RAND_MAX;
       real x = (i + r) / N * SIZE;
       real y = (j + r) / N * SIZE; 
       real z = (k + r) / N * SIZE;

       ParticleStorage::reference ref = particleStorage.addParticle();
       ParticleStorage::PropertyTraits<Real3D>::Reference positionRef = 
	 particleStorage.getProperty<Real3D>(position);
       ParticleStorage::PropertyTraits<Real3D>::Reference velocityRef = 
	 particleStorage.getProperty<Real3D>(velocity);
       ParticleStorage::PropertyTraits<Real3D>::Reference forceRef    = 
	 particleStorage.getProperty<Real3D>(force);

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

  AllParticles allSet(&particleStorage);

  // define allpairs with (x, y) for all x, y in allSet

  AllPairs allpairs(pbc, allSet, position);

  // For test only: PairWriter prints each particle pair 

  PairWriteComputer pairWriter(&particleStorage, position);

  // call pairWriter(ref1, ref2) for each particle ref pair of allSet

  // allpairs.foreach(pairWriter);

  // define LennardJones interaction

  LennardJones ljint;

  ljint.setCutoff(2.5);
  ljint.setEpsilon(1.0);
  ljint.setSigma(1.0);

  // make a FENE interaction

  FENE fene;

  // force will be the vector of all forces in the particle storage
  // and force[ref] returns the force (as RealArrayRef) of particle reference ref

  ParticleStorage::PropertyTraits<Real3D>::Reference forceRef = particleStorage.getProperty<Real3D>(force);

  // Define a pair computer that computes the forces for particle pairs
  // ljint provides the routine computeForce for a particle pair
  // force (pointer to all forces of particles) tells us where the computed forces are added

  PairForceComputer forcecompute(forceRef, ljint);

  // call forcecompute(ref1, ref2) for each particle ref pair of allSet

  allpairs.foreach(forcecompute);

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
