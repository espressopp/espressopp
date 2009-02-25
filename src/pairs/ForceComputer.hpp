#ifndef _PAIRS_FORCECOMPUTER_HPP
#define _PAIRS_FORCECOMPUTER_HPP

#include "types.hpp"

#include "interaction/Interaction.hpp"
#include "particles/Storage.hpp"

#include "Computer.hpp"

namespace espresso {
   namespace pairs {

     /** Class that provides a functional for the force calculation. This class
         implements the abstract class ParticlePairComputer.
     */

     class ForceComputer: public Computer {
     private: 
       typedef particles::Storage::PropertyTraits<Real3D>::Reference RealArrayRef;
       typedef espresso::interaction::Interaction Interaction;

       RealArrayRef force;
       const Interaction& interaction;
       Real3D pressure;
       bool computesPressure;
       
     public:
      /** Constructor for force computations needs the force property and the interaction.

         \item _force is the particle property that stands for the force
         \itme _interaction provides the function that computes the force between two particles.
      */
       ForceComputer(RealArrayRef _force, const Interaction& _interaction) 
         : force(_force),
           interaction(_interaction),
           pressure(0.0)
           {}
       
       // TODO: why do I have to care what kind of reference to a particle I have?
       virtual void operator()(const Real3D &dist,
                               const espresso::particles::Set::reference p1,
                               const espresso::particles::Set::reference p2)
       {
         Real3D f = interaction.computeForce(dist, p1, p2);

         force[p1] += f;
         force[p2] -= f;

         if (computesPressure) pressure = pressure + f * dist;
       }
     };
  }
}

#endif
