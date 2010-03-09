#ifndef _INTERACTION_LENNARDJONES_HPP
#define _INTERACTION_LENNARDJONES_HPP

#include "Interaction.hpp"
#include "VerletListInteractionTemplate.hpp"

namespace espresso {
  namespace interaction {
    class LennardJonesFunction
      : public InteractionFunction {
    private:
      real epsilon;
      real sigma;
      real shift;
      real ff1, ff2;
      real ef1, ef2;

      void preset();
      
    public:

      void setParameters(real _epsilon, real _sigma) {
	epsilon = _epsilon; 
	sigma = _sigma; 
	preset();
      }
      void setEpsilon(real _epsilon) { 
	epsilon = _epsilon; 
	preset(); 
      }
      real getEpsilon() const { return epsilon; }

      void setSigma(real _sigma) { 
	sigma = _sigma; 
	preset();
      }
      real getSigma() const { return sigma; }

      void setShift(real _shift) { shift = _shift; }
      real getShift() const { return shift; }

      /** sets the shift to the value of the function at _x. */
      void setAutoShift();
      real getEnergy(Particle &p1, Particle &p2,
		     const Real3DRef dist,
		     const real distSqr) const;
      
      void getForce(Particle &p1, Particle &p2,
		    const Real3DRef dist,
		    const real distSqr, Real3DRef force) const;
    };

    class VerletListLennardJones
      : public VerletListInteractionTemplate< LennardJonesFunction > 
    {
      VerletListLennardJones(shared_ptr< VerletList > _verletList) 
	: VerletListInteractionTemplate < LennardJonesFunction >(_verletList) 
      { ntypes = 0; }
      
      ~VerletListLennardJones() {}

      void setAutoShift(int type1, int type2);
    };

    //////////////////////////////////////////////////
    // INLINE IMPLEMENTATION
    //////////////////////////////////////////////////
    inline real 
    LennardJonesFunction::
    getEnergy(Particle &p1, Particle &p2,
		  const Real3DRef dist,
		  const real distSqr) const
    {
      real frac2  = 1.0 / distSqr;
      real frac6  = frac2 * frac2 * frac2;
      real energy = frac6 * (ef1 * frac6 - ef2) - shift;
      // printf ("Particle %d - %d, dist = %6.3g, energy = %7.3g\n", 
      //         p1.p.id, p2.p.id, distSqr, energy);
      return energy;
    }
    
    inline void 
    LennardJonesFunction::
    getForce(Particle &p1, Particle &p2,
	     const Real3DRef dist,
	     const real distSqr, Real3DRef force) const
    {
      real frac2   = 1.0 / distSqr;
      real frac6   = frac2 * frac2 * frac2;
      real ffactor = frac6 * (ff1 * frac6 - ff2) * frac2;
      
      // printf ("force of distsq %f  = %f, ff1 = %f, ff2 = %f\n", distSqr, ffactor, ff1, ff2);
      
      for (int i = 0; i < 3; i++) {
	force[i] = dist[i] * ffactor;
      }
    }


  }
}

#endif
