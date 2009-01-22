#include "interaction/Interaction.hpp"

namespace espresso {
    namespace interaction {
	class LennardJones: public Interaction {
	private:
	    real sigma;
	    real epsilon;
	    real cutoff;
	    real cutoffSqr;
	public:
	    LennardJones() {}
	    virtual ~LennardJones() {}

	    virtual real computeEnergy (real distSqr,
					const espresso::particleset::ParticleSet::const_reference p1,
					const espresso::particleset::ParticleSet::const_reference p2) const {
		real frac2;
		real frac6;
	
                real energy = 0.0;

		if (distSqr < cutoffSqr) {
		    frac2 = sigma / distSqr;
		    frac6 = frac2 * frac2 * frac2;
		    energy = 4.0 * epsilon * (frac6 * frac6 - frac6);
		} 
		    
                return energy;
	    }

	    virtual Real3D computeForce (Real3D dist,
					const espresso::particleset::ParticleSet::const_reference p1,
					const espresso::particleset::ParticleSet::const_reference p2) const {
                Real3D f = 0.0;
		real   frac2;
		real   frac6;
	
                real distSqr = dist.sqr();
                
		if (distSqr < cutoffSqr) {
		    frac2 = sigma / distSqr;
		    frac6 = frac2 * frac2 * frac2;
		    real ffactor = 48.0 * epsilon * (frac6*frac6 - 0.5 * frac6) * frac2;
                    printf ("computeForce, distSqr = %f, ffactor = %f\n", distSqr, ffactor);
                    f = dist * ffactor;
		} 

                return f;
	    }

	    virtual real getCutoff() const { return cutoff; }
	    virtual real getCutoffSqr() const { return cutoffSqr; }
	    virtual void setCutoff(real _cutoff) { 
		cutoff = _cutoff; 
                cutoffSqr = cutoff * cutoff;
	    }

	    virtual void setEpsilon(real _epsilon) { epsilon = _epsilon; }
	    virtual void setSigma(real _sigma) { sigma = _sigma; }
	};
    }
}
