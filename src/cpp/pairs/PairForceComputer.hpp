#include "types.hpp"

#include "interaction/Interaction.hpp"
#include "particlestorage/ParticleStorage.hpp"

#include "ParticlePairComputer.hpp"


namespace espresso {

   namespace pairs {

     class PairForceComputer: public ParticlePairComputer {

     public:

       typedef particlestorage::ParticleStorage::PropertyTraits<Real3D>::Reference RealArrayRef;

     private: 

       typedef espresso::interaction::Interaction Interaction;

       RealArrayRef force;

       const Interaction& interaction;
    
       Real3D pressure;

       bool computesPressure;
       
     public:

       PairForceComputer(RealArrayRef _force, const Interaction& _interaction) 

         : force(_force),
           pressure(0.0),
           interaction(_interaction) {}
    
       virtual void operator()(const Real3D &dist,
                               const espresso::particleset::ParticleSet::reference p1,
                               const espresso::particleset::ParticleSet::reference p2)

       {
         Real3D f = interaction.computeForce(dist, p1, p2);

         force[p1] += f;
         force[p2] -= f;

         if (computesPressure) pressure = pressure + f * dist;
       }
     };
  }
}
