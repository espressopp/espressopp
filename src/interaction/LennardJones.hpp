#ifndef _INTERACTION_LENNARDJONES_HPP
#define _INTERACTION_LENNARDJONES_HPP

#include "InteractionTemplate.hpp"

namespace espresso {

  namespace interaction {

    /** This class contains member variables for the Lennard Jones
        parameters and routines to compute the energy and the force
        between two particles.
    */

    class LennardJonesParameters : public Interaction::ParametersBase {

      private:

      real epsilon;
      real sigma;

      public:

      void setEpsilon(real _epsilon) { epsilon = _epsilon; }
      real getEpsilon() const { return epsilon; }

      void setSigma(real _sigma) { sigma = _sigma; }
      real getSigma() const { return sigma; }

      /** Inline method for reimplementation of methods to compute energies */

      real computeEnergy(Particle &p1, Particle &p2,
                         const real dist[3],
                         const real distSqr) const 
      {
        real frac2 = sigma * sigma / distSqr;
        real frac6 = frac2 * frac2 * frac2;
        real energy = 4.0 * epsilon * (frac6 * frac6 - frac6);
        printf ("Particle %d - %d, dist = %6.3g, energy = %7.3g\n", 
                 p1.p.id, p2.p.id, distSqr, energy);
        return energy;
      }

      /** Compute force */

      void computeForce(Particle &p1, Particle &p2,
                        const real dist[3],
                        const real distSqr, real force[3]) const 
      {
        real distSqrInv = 1.0 / distSqr;

        real frac2   = sigma * sigma * distSqrInv;
        real frac6   = frac2 * frac2 * frac2;
        real ffactor = 48.0 * epsilon * (frac6 * frac6 - 0.5 * frac6) * distSqrInv;
    
        for (int i = 0; i < 3; i++)
          force[i] = dist[i] * ffactor;
      }

    };

    class LennardJones : public InteractionTemplate< LennardJonesParameters > {

    public:

      virtual void addVerletListForces(shared_ptr<VerletList> vl);

      virtual real computeVerletListEnergy(shared_ptr<VerletList> vl);

      virtual real computeCellEnergy(ParticleList &pl);

      virtual real computeCellEnergy(ParticleList &pl1, ParticleList &pl2);

      LennardJones();

      ~LennardJones();
    
      void setParameters(int type1, int type2, real ep, real sg, real rc);

    };
  }
}

#endif
