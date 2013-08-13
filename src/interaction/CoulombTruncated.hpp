// ESPP_CLASS
#ifndef _INTERACTION_COULOMBTRUNCATED_HPP
#define _INTERACTION_COULOMBTRUNCATED_HPP

#include "Potential.hpp"
#include "FixedPairListInteractionTemplate.hpp"

namespace espresso {
  namespace interaction {
    /* This class provides methods to compute forces and energies of
	the truncated Coulomb potential.
    */
    class CoulombTruncated : public PotentialTemplate< CoulombTruncated > {
    private:
      real qq;

    public:
      static void registerPython();

      CoulombTruncated()
	: qq(0.0) {
	setShift(0.0);
	setCutoff(infinity);
      }

      CoulombTruncated(real _qq,
		   real _cutoff, real _shift)
	: qq(_qq) {
	setShift(_shift);
	setCutoff(_cutoff);
      }

      CoulombTruncated(real _qq,
		   real _cutoff)
	: qq(_qq)
      {
	autoShift = false;
	setCutoff(_cutoff);
	setAutoShift();
      }

      // Setter and getter
      void setQQ(real _qq) {
	qq = _qq;
	updateAutoShift();
      }
      real getQQ() const { return qq; }

      real _computeEnergySqrRaw(real distSqr) const {
	real energy = qq / sqrt(distSqr);
	return energy;
      }

      bool _computeForceRaw(Real3D& force,
                            const Real3D& dist,
                            real distSqr) const {

        real ffactor;
	ffactor = qq / pow(sqrt(distSqr), 3);
        force = dist * ffactor;
        return true;
      }

    };
    // provide pickle support
    struct CoulombTruncated_pickle : boost::python::pickle_suite
    {
      static
      boost::python::tuple
      getinitargs(CoulombTruncated const& pot)
      {
    	  real q2;
          real rc;
          real sh;
          q2 =pot.getQQ();
          rc =pot.getCutoff();
          sh =pot.getShift();
          return boost::python::make_tuple(q2, rc, sh);
      }
    };
  }
}

#endif
