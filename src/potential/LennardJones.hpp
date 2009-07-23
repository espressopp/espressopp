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
    class LennardJones 
      : public CentralPotentialTemplate< LennardJones >
    {
    public:
      typedef shared_ptr< LennardJones > SelfPtr;

    private:
      real epsilon;
      real sigma;
      real cutoff;
      real cutoffSqr;

      static LOG4ESPP_DECL_LOGGER(theLogger);

    public:
      static void registerPython();

      LennardJones() {}
      virtual ~LennardJones() {}

      // Setter and getter
      void setEpsilon(real _epsilon) { epsilon = _epsilon; }
      real getEpsilon() const { return epsilon; }

      void setSigma(real _sigma) { sigma = _sigma; }
      real getSigma() const { return sigma; }

      void setCutoff(real _cutoff) { cutoff = _cutoff; cutoffSqr = cutoff*cutoff; }
      real getCutoff() const { return cutoff; }

      virtual real computeEnergySqr(const real distSqr) const;
      virtual Real3D computeForce(const Real3D dist) const;

      virtual real getCutoffSqr() const { return cutoffSqr; }
    };

  }
}

#endif
