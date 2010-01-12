#ifndef _INTERACTION_LENNARDJONES_HPP
#define _INTERACTION_LENNARDJONES_HPP

#include "InteractionTemplate.hpp"

namespace espresso {
  namespace interaction {

    class LennardJones : public InteractionTemplate< LennardJones > {

    public:

      LennardJones() {}
    
      void setParameters(int type1, int type2, real ep, real sg, real rc) {
        Parameters& params = parameterArray[type1][type2];
        params.setEpsilon(ep);
        params.setSigma(sg);
        params.setCutoff(rc);
      }

    public:

      // friend class InteractionTemplate<>;

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

      /** Inline method for efficient reimplementations */

      inline double getCutoffSqr(int type1, int type2) { 
        const Parameters& params = parameterArray[type1][type2];
        return params.getCutoffSqr(); 
      }

      /** Inline method for reimplementation of methods to compute energies */

      inline real computeEnergy(Particle &p1, Particle &p2,
                                int type1, int type2, const real dist[3],
                                const real distSqr) const {
        const Parameters& params = parameterArray[type1][type2];
        real frac2 = params.getSigma() * params.getSigma() / distSqr;
        real frac6 = frac2 * frac2 * frac2;
        real energy = 4.0 * params.getEpsilon() * (frac6 * frac6 - frac6);
        printf ("Particle %d - %d, dist = %6.3g, energy = %7.3g\n", 
                 p1.p.id, p2.p.id, distSqr, energy);
        return energy;
      }

      /** Inline method for reimplementation of methods to compute forces */

      inline void computeForce(Particle &p1, Particle &p2,
                               int type1, int type2, const real dist[3],
                               const real distSqr, real force[3]) const 
      {
        const Parameters& params = parameterArray[type1][type2];
        real distSqrInv = 1.0 / distSqr;
        real frac2 = params.getSigma() * params.getSigma() * distSqrInv;
        real frac6 = frac2 * frac2 * frac2;
        real ffactor = 48.0 * params.getEpsilon() * (frac6 * frac6 - 0.5 * frac6) * distSqrInv;
    
        for (int i = 0; i < 3; i++)
          force[i] = dist[i] * ffactor;
      }

     private:

      /** ToDo: dynamic / flexible parameter array */

      Parameters parameterArray[1][1];

    };
  }
}

#endif
