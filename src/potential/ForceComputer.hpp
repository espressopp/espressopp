#ifndef _POTENTIAL_FORCECOMPUTER_HPP
#define _POTENTIAL_FORCECOMPUTER_HPP

#include "types.hpp"
#include "Property.hpp"
#include "particles/Storage.hpp"
#include "pairs/Computer.hpp"

namespace espresso {
  namespace potential {
    template < class Potential >
    class ForceComputer : public pairs::Computer {
      // keep a pointer so that the object isn't destroyed
      typename Potential::SelfPtr potentialptr;
      // keep a reference for fast access
      Potential &potential;

      Property< Real3D >::SelfPtr forceProperty1;
      Property< Real3D >::SelfPtr forceProperty2;

      particles::PropertyHandle< Real3D > force1;
      particles::PropertyHandle< Real3D > force2;

      bool computesVirial;
      Real3D virial;
            
    public:
      ForceComputer(typename Potential::SelfPtr _potential,
		    Property< Real3D >::SelfPtr _forceProperty1,
		    Property< Real3D >::SelfPtr _forceProperty2,
		    bool _computesVirial = false) 
	: potentialptr(_potential), potential(*_potential), 
	  forceProperty1(_forceProperty1), forceProperty2(_forceProperty2),
	  computesVirial(_computesVirial)
      {}

      virtual ~ForceComputer() {};

      virtual void prepare(particles::Storage::SelfPtr storage1, 
			   particles::Storage::SelfPtr storage2) {
	force1 = forceProperty1->getHandle(storage1);
	force2 = forceProperty2->getHandle(storage2);
	virial = 0.0;
      }

      virtual void apply(const Real3D dist, 
			 const particles::ParticleHandle p1, 
			 const particles::ParticleHandle p2) {
	Real3D f = potential.computeForce(dist);
	force1[p1] += f;
	force2[p2] -= f;
        if (computesVirial) 
	  virial = virial + f * dist;
      }
    };


//     /** Class that provides a functional form for the force
//         calculation. Every Potential needs to provide such a force
//         computer. For standard forms of potential, there are
//         ForceCalculator facades automatically implementing a
//         ForceComputer. The facades are:
//         <ul>Storage
//         <li> VectorForceComputerFacade
//         </ul>

//         This class actually just serves as a container for the data
//         required by the force computer of _any_ interaction. By
//         passing and copying the data of an object of this type, we
//         can bypass that the interaction needs to know the signature
//         of the constructor.
//     */
//     class ForceComputer: public pairs::Computer {
//     public:
//       /** Constructor for force computations needs the force property
//           @param _force is the particle property that stands for the force
//           @param _computesPressure whether the total virial is also computed
//       */
//       ForceComputer(particles::PropertyHandle<Real3D> _force,
//                     bool _computesVirial = false) 
//         : force(_force), virial(0.0) {}
//       virtual ~ForceComputer() {};

//       virtual void operator()(const Real3D &dist,
//                               const particles::ParticleHandle p1,
//                               const particles::ParticleHandle p2) {};
       
//     protected: 
//       particles::PropertyHandle<Real3D> force;
//       Real3D virial;
//       bool computesVirial;

//       void addContribution(const Real3D &f,
//                            const Real3D &dist,
//                            const particles::ParticleHandle p1,
//                            const particles::ParticleHandle p2) {
//         force[p1] += f;
//         force[p2] -= f;
//         if (computesVirial) virial = virial + f * dist;
//       }
//     };

//     /** Class that provides a functional form for the force
//         calculation of interactions that only depend on the minimum
//         image vector between the two particles. InterBasicComputer
//         provides the basic calculation and constants for the
//         potential, and should have the following form:

//         class InterBasicComputer {
//         Real3D computeForce(const Real3D &dist);
//         };
//     */
//     template <class InterBasicComputer>
//     class VectorForceComputerFacade: public ForceComputer {
//     public:
//       VectorForceComputerFacade(const ForceComputer &_forceComputer,
//                                 const InterBasicComputer &_computer)
//         : ForceComputer(_forceComputer), computer(_computer) {}
//       virtual ~VectorForceComputerFacade() {};
       
//       virtual void operator()(const Real3D &dist,
//                               const particles::ParticleHandle p1,
//                               const particles::ParticleHandle p2) {
//         addContribution(computer.computeForce(dist), dist, p1, p2);
//       }
       
//     private:
//       InterBasicComputer computer;
//     };
  }
}

#endif
