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
    class _LennardJones {
    private:
      real epsilon;
      real sigma;
      real cutoff;
      real cutoffSqr;

      static LOG4ESPP_DECL_LOGGER(theLogger);

    public:
      static void registerPython();

      _LennardJones(real _epsilon, real _sigma, real _cutoff) {
	setEpsilon(_epsilon);
	setSigma(_sigma);
	setCutoff(_cutoff);
      }

      // Setter and getter
      void setEpsilon(real _epsilon) { epsilon = _epsilon; }
      real getEpsilon() const { return epsilon; }

      void setSigma(real _sigma) { sigma = _sigma; }
      real getSigma() const { return sigma; }

      void setCutoff(real _cutoff) { cutoff = _cutoff; cutoffSqr = cutoff*cutoff; }
      real getCutoff() const { return cutoff; }

      real _getCutoffSqr() const { return cutoffSqr; }

      real _computeEnergySqr(const real distSqr) const;
      Real3D _computeForce(const Real3D dist) const;
    };

    typedef CentralPotentialWrapper< _LennardJones > LennardJones;

  }
}

#endif
