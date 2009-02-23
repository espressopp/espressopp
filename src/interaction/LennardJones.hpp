#ifndef _INTERACTION_LENNARDJONES_HPP
#define _INTERACTION_LENNARDJONES_HPP

#include <pmi.hpp>
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

      IF_MPI(pmi::ParallelClass<LennardJones> pmiObject;)

      real sigma;       
      real epsilon;
      real cutoff;
      real cutoffSqr;

      typedef espresso::particleset::ParticleSet::const_reference const_reference;

      static LOG4ESPP_DECL_LOGGER(theLogger);

    public:
      IF_PYTHON(static void registerPython();)

      /** Default constructor. */
      LennardJones();
      /** Destructor. */
      virtual ~LennardJones();

      // Setter and getter
      PMI_DECL_SETTER(setCutoff, real, _cutoff);
      virtual real getCutoff() const;

      PMI_DECL_SETTER(setEpsilon, real, _epsilon);
      virtual real getEpsilon() const;

      PMI_DECL_SETTER(setSigma, real, _sigma);
      virtual real getSigma() const;

      // PMI and Python visible
      virtual real computeEnergy(const Real3D &dist) const;  
      virtual real computeEnergy(const real dist) const;
      virtual Real3D computeForce(const Real3D &dist) const;

      // NOT visible on PMI/Python:
      virtual real computeEnergy(const Real3D &dist,
				 const const_reference p1,
				 const const_reference p2) const;
      virtual real computeEnergySqr(const real distSqr) const;

      virtual Real3D computeForce(const Real3D &dist,
				  const const_reference p1,
				  const const_reference p2) const;

      virtual real getCutoffSqr() const;

    };
  }
}

#endif
