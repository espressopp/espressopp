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
      real shift;
      real ff1, ff2;
      real ef1, ef2;

      public:

      void setEpsilon(real _epsilon) { epsilon = _epsilon; }
      real getEpsilon() const { return epsilon; }

      void setSigma(real _sigma) { sigma = _sigma; }
      real getSigma() const { return sigma; }

      void setShift(real _shift) { shift = _shift; }
      real getShift() const { return shift; }

      void enableShift() {
        real ratio = sigma / getCutoff();
        real rat2 = ratio * ratio;
        real rat6 = rat2 * rat2 * rat2;
        shift = 4.0 * epsilon * rat6 * (rat6 - 1.0);
      }

      void preset() {
        real sig2 = sigma * sigma;
        real sig6 = sig2 * sig2 * sig2;
        ff1 = 48.0 * epsilon * sig6 * sig6;
        ff2 = 24.0 * epsilon * sig6;
        ef1 =  4.0 * epsilon * sig6 * sig6;
        ef2 =  4.0 * epsilon * sig6;
      }

      /** Inline method for reimplementation of methods to compute energies */

      real computeEnergy(Particle &p1, Particle &p2,
                         const real dist[3],
                         const real distSqr) const 
      {
        real frac2  = 1.0 / distSqr;
        real frac6  = frac2 * frac2 * frac2;
        real energy = frac6 * (ef1 * frac6 - ef2) - shift;
        // printf ("Particle %d - %d, dist = %6.3g, energy = %7.3g\n", 
        //         p1.p.id, p2.p.id, distSqr, energy);
        return energy;
      }

      /** Compute force */

      void computeForce(Particle &p1, Particle &p2,
                        const real dist[3],
                        const real distSqr, real force[3]) const 
      {
        real frac2   = 1.0 / distSqr;
        real frac6   = frac2 * frac2 * frac2;
        real ffactor = frac6 * (ff1 * frac6 - ff2) * frac2;

        // printf ("force of distsq %f  = %f, ff1 = %f, ff2 = %f\n", distSqr, ffactor, ff1, ff2);

        for (int i = 0; i < 3; i++) {
          force[i] = dist[i] * ffactor;
        }
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

      void enableShift(int type1, int type2);
    };
  }
}

#endif
