#ifndef _INTERACTION_LENNARDJONES_HPP
#define _INTERACTION_LENNARDJONES_HPP

#include <pmi/pmi.hpp>
#include <logging.hpp>
#include <interaction/Interaction.hpp>

namespace espresso {
  namespace interaction {

    /** This class provides routines to compute forces and energies
	of the Lennard Jones potential.

	\f[ V(r) = 4 \varepsilon \left[ \left( \frac{\sigma}{r} \right)^{12} -
	\left( \frac{\sigma}{r} \right)^{6} \right]
	\f]

    */

    class LennardJones: public Interaction {

    private:

      real sigma;       
      real epsilon;
      real cutoff;
      real cutoffSqr;

      typedef espresso::particleset::ParticleSet::const_reference const_reference;

      static LOG4ESPP_DECL_LOGGER(theLogger);

    public:

      /** Default constructor. Member variables are accessed via setters and getters. */

      LennardJones();
       
      virtual ~LennardJones();

      virtual real computeEnergy (const Real3D &dist,
				  const const_reference p1,
				  const const_reference p2) const;
      virtual real computeEnergy (const Real3D &dist) const;  
      virtual real computeEnergy(const real dist) const;
      virtual real computeEnergySqr (const real distSqr) const;
      virtual Real3D computeForce (const Real3D &dist,
				   const const_reference p1,
				   const const_reference p2) const;
      virtual Real3D computeForce (const Real3D &dist) const;
      virtual real getCutoff() const;
      virtual real getCutoffSqr() const;
      virtual void setCutoff(real _cutoff);

      virtual void setEpsilon(real _epsilon);
      virtual real getEpsilon() const;
      virtual void setSigma(real _sigma);
      virtual real getSigma() const;
    };
  }
}

#endif
