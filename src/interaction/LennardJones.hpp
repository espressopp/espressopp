#ifndef _INTERACTION_LENNARDJONES_HPP
#define _INTERACTION_LENNARDJONES_HPP

#include "InteractionTemplate.hpp"

namespace espresso {
  class LennardJones : public InteractionTemplate< LennardJones > {
  public:
    LennardJones() {}
    
    void setParameters(int type1, int type2, real ep, real sg, real rc) {
      Parameters& params = createParameters(type1, type2);
      params.setEpsilon(ep);
      params.setSigma(sg);
      params.setCutoff(rc);
    }

  private:
    friend class InteractionTemplate<>;
    /* nested class to store coefficients */
    class Parameters : public ParametersBase {
      real epsilon;
      real sigma;

    public:
      void setEpsilon(real _epsilon) { epsilon = _epsilon; }
      real getEpsilon() const { return epsilon; }

      void setSigma(real _sigma) { sigma = _sigma; }
      real getSigma() const { return sigma; }
    };

    real computeEnergy(Particle &p1, Particle &p2, 
		       Parameters &params, const real dist[3], 
		       const real distSqr) const;
    
    void computeForce(Particle &p1, Particle &p2, 
		      Parameters &params, const real dist[3],
		      const real distSqr,
		      real force[3]) const;
  };
}

#endif
