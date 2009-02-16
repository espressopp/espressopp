#include <cmath>
#include "logging.hpp"
#include "interaction/Interaction.hpp"

namespace espresso {
  namespace interaction {

    /** This class provides routines to computer forces and energies
	based on the FENE potential.

	\f[ V(r) = -\frac{1}{2} \Delta r_\text{max}^2 K \log \left[ 1 - \left(
	\frac{r-r_0}{\Delta r_\text{max}} \right)^2 \right]

    */

    class FENE: public Interaction {

    private:

      real K;       
      real r0;
      real rMax;

      typedef espresso::particleset::ParticleSet::const_reference const_reference;

      static LOG4ESPP_DECL_LOGGER(theLogger);

    public:

      /** Default constructor. Member variables are accessed by setter and getter. */

      FENE() {}
       
      virtual ~FENE() {}

      virtual real computeEnergy (const Real3D &dist,
				  const const_reference p1,
				  const const_reference p2) const {
	         	
	return computeEnergy(dist);
      }

      virtual real computeEnergy (const Real3D &dist) const {
	return computeEnergySqr(dist.sqr());
      }
      
      virtual real computeEnergy (const real dist) const {
	real energy = -0.5 * pow(rMax, 2) * K * log(1 - pow((dist - r0) / rMax, 2));
	return energy;
      }

      virtual real computeEnergySqr (real distSqr) const {
	return computeEnergy(sqrt(distSqr));
      }
      
      virtual Real3D computeForce (const Real3D &dist,
				   const const_reference p1,
				   const const_reference p2) const {
	return computeForce(dist);
      }

      virtual Real3D computeForce (const Real3D &dist) const {
	Real3D f = 0.0;
	real r = sqrt(dist.sqr());

	real ffactor = -K * (r - r0) / (1 - pow((r - r0) / rMax, 2)) / r;
	f = ffactor * dist;

	return f;
      }
      
      
      virtual void setK(real _K) { K = _K; }
      virtual void setr0(real _r0) { r0 = _r0; }
      virtual void setrMax(real _rMax) { rMax = _rMax; }
      
      /* FENE should probably derived from a two-body interaction
	 without a cutoff to avoid the following */
      virtual real getCutoff() const { return -1; }
      virtual real getCutoffSqr() const { return -1; }
      
    };
  }
}
