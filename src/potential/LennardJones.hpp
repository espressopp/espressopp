#ifndef _POTENTIAL_LENNARDJONES_HPP
#define _POTENTIAL_LENNARDJONES_HPP

#include <logging.hpp>
#include <potential/CentralPotential.hpp>

namespace espresso {
  namespace potential {
    /** This class provides routines to compute forces and energies
	of the Lennard Jones potential.

	\f[ V(r) = 4 \varepsilon \left[ \left( \frac{\sigma}{r} \right)^{12} -
	\left( \frac{\sigma}{r} \right)^{6} \right]
	\f]

    */
    class LennardJones: public CentralPotential {
    public:
      typedef boost::shared_ptr< LennardJones > SelfPtr;

      class BasicComputer {
        friend class LennardJones;

        real epsilon;
        real sigma;
        real cutoffSqr;

      public:
        Real3D computeForce(const Real3D &) const;
        real   computeEnergySqr(const real) const;
      };
      friend class BasicComputer;

    private:
      BasicComputer computer;
      real cutoff;

      static LOG4ESPP_DECL_LOGGER(theLogger);

    public:
      static void registerPython();

      /** Default constructor. */
      LennardJones();
      /** Destructor. */
      virtual ~LennardJones();

      // Setter and getter
      virtual void set(real _epsilon, real _sigma, real _cutoff);

      real getEpsilon() const;
      real getSigma() const;
      real getCutoff() const;

      virtual real computeEnergySqr(const real distSqr) const;
      virtual Real3D computeForce(const Real3D &dist) const;

      virtual pairs::EnergyComputer *createEnergyComputer(const pairs::EnergyComputer &) const;
      virtual pairs::ForceComputer  *createForceComputer (const pairs::ForceComputer &)  const;

      virtual real getCutoffSqr() const;
    };

  }
}

#endif
