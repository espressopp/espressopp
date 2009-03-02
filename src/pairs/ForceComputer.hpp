#ifndef _PAIRS_FORCECOMPUTER_HPP
#define _PAIRS_FORCECOMPUTER_HPP

#include "types.hpp"

#include "particles/Storage.hpp"

#include "Computer.hpp"

namespace espresso {
  namespace pairs {

    /** Class that provides a functional form for the force
        calculation. Every Potential needs to provide such a force
        computer. For standard forms of potential, there are
        ForceCalculator facades automatically implementing a
        ForceComputer. The facades are:
        <ul>
        <li> VectorForceComputerFacade
        </ul>

        This class actually just serves as a container for the data
        required by the force computer of _any_ interaction. By
        passing and copying the data of an object of this type, we
        can bypass that the interaction needs to know the signature
        of the constructor.
    */
    class ForceComputer: public pairs::Computer {
      typedef particles::Storage::PropertyTraits<Real3D>::Reference RealArrayRef;

    public:
      /** Constructor for force computations needs the force property
          @param _force is the particle property that stands for the force
          @param _computesPressure whether the total virial is also computed
      */
      ForceComputer(RealArrayRef _force, bool _computesVirial = false) 
        : force(_force), virial(0.0) {}
      virtual ~ForceComputer() {};

      virtual void operator()(const Real3D &dist,
                              const particles::Set::reference p1,
                              const particles::Set::reference p2) {};
       
    protected: 
      RealArrayRef force;
      Real3D virial;
      bool computesVirial;

      void addContribution(const Real3D &f,
                           const Real3D &dist,
                           const particles::Set::reference p1,
                           const particles::Set::reference p2) {
        force[p1] += f;
        force[p2] -= f;
        if (computesVirial) virial = virial + f * dist;
      }
    };

    /** Class that provides a functional form for the force
        calculation of interactions that only depend on the minimum
        image vector between the two particles. InterBasicComputer
        provides the basic calculation and constants for the
        potential, and should have the following form:

        class InterBasicComputer {
        Real3D computeForce(const Real3D &dist);
        };
    */
    template <class InterBasicComputer>
    class VectorForceComputerFacade: public ForceComputer {
    public:
      VectorForceComputerFacade(const ForceComputer &_forceComputer,
                                const InterBasicComputer &_computer)
        : ForceComputer(_forceComputer), computer(_computer) {}
      virtual ~VectorForceComputerFacade() {};
       
      virtual void operator()(const Real3D &dist,
                              const particles::Set::reference p1,
                              const particles::Set::reference p2) {
        addContribution(computer.computeForce(dist), dist, p1, p2);
      }
       
    private:
      InterBasicComputer computer;
    };
  }
}

#endif
