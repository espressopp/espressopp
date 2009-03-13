#ifndef _PAIRS_ENERGYCOMPUTER_HPP
#define _PAIRS_ENERGYCOMPUTER_HPP

#include "types.hpp"

#include "particles/Storage.hpp"

#include "Computer.hpp"

namespace espresso {
  namespace pairs {

    /** Class that provides a functional form for the energy
        calculation. Every short ranged potential needs to provide
        such an energy computer. For convenience, several facades are
        provided for this class; you only need to implement a
        derivate of this class directly, if you need to access
        additional information of particles, i.e. you need all
        information provided to the ()-operator. The facades are:
        <ul>
        <li> SquareDistEnergyComputerFacade
        </ul>


        This class actually just serves as a container for the data
        required by the energy computer of _any_ interaction. By
        passing and copying the data of an object of this type, we
        can bypass that the interaction needs to know the signature
        of the constructor.
    */
    class EnergyComputer: public pairs::Computer {
    public:
      /** Constructor for energy computations needs the energy property.

          @param _energy is the particle property that contains the energy
      */
      EnergyComputer(particles::PropertyReference<esutil::real> _energy) 
        : energy(_energy), totalEnergy(0.0) {}
      virtual ~EnergyComputer() {};

      virtual void operator()(const esutil::Real3D &dist,
                              const particles::ParticleReference p1,
                              const particles::ParticleReference p2) {};

      esutil::real getAccumulatedEnergy() const { return totalEnergy; }

    protected:
      particles::PropertyReference<esutil::real> energy;
      esutil::real totalEnergy;
      bool computesVirial;

      void addContribution(esutil::real e,
                           const particles::ParticleReference p1,
                           const particles::ParticleReference p2) {
        energy[p1] += 0.5*e;
        energy[p2] += 0.5*e;
        totalEnergy += e;
      }
    };

    /** Class that provides a functional form for the energy
        calculation for a potential depending only on the square of
        the particle distance. InterBasicComputer provides the very
        basic calculation of an interaction. It should provide at
        least the following interface:
         
        class InterBasicComputer {
        real computeEnergySqr(const real distanceSquared);
        };

        For other facades, see EnergyComputer.
    */
    template <class InterBasicComputer>
    class SquareDistEnergyComputerFacade: public EnergyComputer {
    public:
      SquareDistEnergyComputerFacade(const EnergyComputer &_energyComputer,
                                     const InterBasicComputer &_computer)
        : EnergyComputer(_energyComputer), computer(_computer)
      {}
      virtual ~SquareDistEnergyComputerFacade() {};

      virtual void operator()(const esutil::Real3D &dist,
                              const particles::ParticleReference p1,
                              const particles::ParticleReference p2) {
        addContribution(computer.computeEnergySqr(dist.sqr()), p1, p2);
      }

    private:
      InterBasicComputer computer;
    };
  }
}

#endif
