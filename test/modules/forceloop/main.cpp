
// Example C++ program doing the same as the force loop Python script

#include <types.hpp>
#include <logging.hpp>

#include <boost/shared_ptr.hpp>

#include <espresso_common.hpp>
#include <particles/Storage.hpp>
#include <particles/All.hpp>
#include <bc/PBC.hpp>
#include <interaction/LennardJones.hpp>
#include <interaction/FENE.hpp>
#include <pairs/All.hpp>
#include <pairs/ForceComputer.hpp>
#include <thermostat/Langevin.hpp>

#include "ParticleWriter.hpp"

// #include "PairWriteComputer.hpp"

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

  boost::shared_ptr<Storage> particleStorage = 
         boost::shared_ptr<Storage>(new Storage());

  boost::shared_ptr<Property<Real3D> > position = 
         boost::shared_ptr<Property<Real3D> >(new Property<Real3D>(particleStorage));
  boost::shared_ptr<Property<Real3D> > velocity = 
         boost::shared_ptr<Property<Real3D> >(new Property<Real3D>(particleStorage));
  boost::shared_ptr<Property<Real3D> > force    = 
         boost::shared_ptr<Property<Real3D> >(new Property<Real3D>(particleStorage));

  // generate particles in the particle storage

  for (int i = 0; i < N; i++) 
  for (int j = 0; j < N; j++) 
  for (int k = 0; k < N; k++) {
      
       real r;
       r = 0.4 + 0.2 * rand() / RAND_MAX;
       real x = (i + r) / N * SIZE;
       real y = (j + r) / N * SIZE; 
       real z = (k + r) / N * SIZE;

       ParticleId id = particleStorage->addParticle();
/*
       PropertyHandle<Real3D> positionHandle = 
	 particleStorage->getPropertyHandle<Real3D>(position);
       PropertyReference<Real3D> velocityRef = 
	 particleStorage.getPropertyReference<Real3D>(velocity);
       PropertyReference<Real3D> forceRef    = 
	 particleStorage.getPropertyReference<Real3D>(force);
*/

       (*position.get())[id] = Real3D(x, y, z);
       (*velocity.get())[id] = Real3D(x, y, z);
       (*force.get())[id] = 0.0;

/*
       positionRef[id] = Real3D(x, y, z);
       velocityRef[id] = Real3D(x, y, z);
       forceRef[id] = 0.0;
*/
  }

  // For test only: ParticleWriter prints each particle

  ParticleWriter pWriter(*particleStorage.get(), *position.get(), *force.get());

  // call pWriter(id) for each particle reference ref of particle storage

  particleStorage->foreach(pWriter);

  // define periodic boundary conditions

  boost::shared_ptr<PBC> pbc = boost::shared_ptr<PBC>(new PBC(SIZE));

  // define a set of all particles

  boost::shared_ptr<particles::All> allSet  = 
         boost::shared_ptr<particles::All>(new particles::All(particleStorage));

  // define allpairs with (x, y) for all x, y in allSet

  boost::shared_ptr<pairs::All> allpairs =
         boost::shared_ptr<pairs::All>(new pairs::All(pbc, allSet, position));

  // For test only: PairWriter prints each particle pair 

  // PairWriteComputer pairWriter(&particleStorage, position);

  // call pairWriter(ref1, ref2) for each particle ref pair of allSet

  // allpairs.foreach(pairWriter);

  // define LennardJones interaction

  boost::shared_ptr<LennardJones> ljint = boost::shared_ptr<LennardJones>(new LennardJones());

  ljint->set(1.0, 1.0, 2.5);

  // make a FENE interaction

  FENE fene;

  fene.set(1.5, 1.0, 2.5);

  // force will be the vector of all forces in the particle storage
  // and force[ref] returns the force (as RealArrayRef) of particle reference ref

  PropertyHandle<Real3D> forceRef = *force.get();

  // Define a pair computer that computes the forces for particle pairs
  // ljint provides the routine computeForce for a particle pair
  // force (pointer to all forces of particles) tells us where the computed forces are added

  // ForceComputer *forcecompute = ljint->createForceComputer(ForceComputer(forceRef));

  // call forcecompute(ref1, ref2) for each particle ref pair of allSet

  // allpairs->foreach(*forcecompute);

  // delete forcecompute;

  // print out all particle data to see that it calculates some forces

  // particleStorage->foreach(pWriter);

  // create integrator and set properties

  boost::shared_ptr<VelocityVerlet> integrator =
        boost::shared_ptr<VelocityVerlet>(new VelocityVerlet(allSet, position, velocity, force));

  integrator->setTimeStep(0.005);
  integrator->addForce(ljint, allpairs);

  // integrate for a certain number of timesteps

  integrator->run(100);

  // check to see that particles have new positions

  particleStorage->foreach(pWriter);

  // make a shared pointer to a Lanevin thermostat with T=298.15 and gamma=0.5
  boost::shared_ptr<thermostat::Langevin> lv =
    boost::shared_ptr<thermostat::Langevin>(new thermostat::Langevin(298.15, 0.5));
  // write out the gamma value
  std::cout << lv->getGamma() << "\t" << lv->getTemperature() << std::endl;
  // set the thermostat for the integrator
  integrator->setThermostat(lv);
  // get the integrator and call its thermalizeA method virtually
  integrator->getThermostat()->thermalizeA();

}

int main() 

{

  printf("configure forceloop\n");

  LOG4ESPP_CONFIGURE();

  forceloop();

}
