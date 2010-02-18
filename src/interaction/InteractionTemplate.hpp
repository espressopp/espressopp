#ifndef _INTERACTION_INTERACTIONTEMPLATE_HPP
#define _INTERACTION_INTERACTIONTEMPLATE_HPP

#include "Interaction.hpp"
#include "BC.hpp"
#include "Storage.hpp"
#include "VerletList.hpp"
#include "esutil/Utils.hpp"

namespace espresso {
  namespace interaction {
    /** InteractionTemplate provides loop templates to compute
	forces and energies of the various interactions. */

    template< typename _Parameters >

    class InteractionTemplate : public Interaction {

    public:

      typedef _Parameters Parameters;

      // ENERGY COMPUTATION
      // full square over two particle lists

      virtual real computeCellEnergyImpl(ParticleList &pl1, ParticleList &pl2)
      {
        LOG4ESPP_INFO(theLogger, "compute energy for a cell pair");

	real e = 0.0;
	for (int i = 0, endi = pl1.size(); i < endi; i++) {
	  Particle &p1 = pl1[i];
          int type1 = p1.p.type;
	  for (int j = 0, endj = pl2.size(); j < endj; j++) {
	    Particle &p2 = pl2[j];
	    real dist[3];
	    real distSqr = esutil::distance2vec(p1.r.p, p2.r.p, dist);
            int type2 = p2.p.type;
            const Parameters& params = parameterArray[type1][type2];
            if (distSqr < params.getCutoffSqr()) {
	       e += params.computeEnergy(p1, p2, dist, distSqr);
            }
	  }
        }
	return e;
      }

      // half square over a list of a single cell

      virtual real computeCellEnergyImpl(ParticleList &pl)
      {
        LOG4ESPP_INFO(theLogger, "compute energy for a single cell");

        real e = 0.0;
        int size = pl.size();
        for (int i = 0; i < size; i++) {
          Particle &p1 = pl[i];
          int type1 = p1.p.type;
          for (int j = i+1; j < size; j++) {
            Particle &p2 = pl[j];
            real dist[3];
            real distSqr = esutil::distance2vec(p1.r.p, p2.r.p, dist);
            int type2 = p2.p.type;
            const Parameters& params = parameterArray[type1][type2];
            if (distSqr < params.getCutoffSqr()) {
               e += params.computeEnergy(p1, p2, dist, distSqr);
            }
          }
        }
        return e;
      }

      /** Reimplementation of pure routine for the derived class; the
          methods getCutoffSqr and computeEnergy will be inlined for efficiency.
      */

      virtual real computeVerletListEnergyImpl(shared_ptr<VerletList> vl) {

        LOG4ESPP_INFO(theLogger, "compute energy by the Verlet List");

        const VerletList::PairList pairList = vl->getPairs();
        int npairs = pairList.size();
        real e = 0.0;
        for (int i = 0; i < npairs; i++) {
          Particle& p1 = *pairList[i].first;
          Particle& p2 = *pairList[i].second;
          int type1 = p1.p.type;
          real dist[3];
          real distSqr = esutil::distance2vec(p1.r.p, p2.r.p, dist);
          int type2 = p2.p.type;
          const Parameters& params = parameterArray[type1][type2];
          if (distSqr < params.getCutoffSqr()) {
             e += params.computeEnergy(p1, p2, dist, distSqr);
          }
        }
        return e;
      }

      // FORCE COMPUTATION
      // full square over two particle lists

      virtual void addForcesStorage() {
      }

      virtual void addVerletListForcesImpl(shared_ptr<VerletList> vl) {

        LOG4ESPP_INFO(theLogger, "add forces computed by the Verlet List");

        const VerletList::PairList pairList = vl->getPairs();

        int npairs = pairList.size();
        real e = 0.0;
        for (int i = 0; i < npairs; i++) {
          Particle& p1 = *pairList[i].first;
          Particle& p2 = *pairList[i].second;
          int type1 = p1.p.type;
          real dist[3];
          real distSqr = esutil::distance2vec(p1.r.p, p2.r.p, dist);
          int type2 = p2.p.type;
          const Parameters& params = parameterArray[type1][type2];

          if (distSqr < params.getCutoffSqr()) {
             real force[3];
             params.computeForce(p1, p2, dist, distSqr, force);
             for (int k = 0; k < 3; k++) {
               p1.f.f[k] += force[k];
               p2.f.f[k] -= force[k];
             }
#if 0
             printf("dist(%d,%d), dist = %f -> %f %f %f\n", p1.p.id, p2.p.id, distSqr, dist[0], dist[1], dist[2]);
             printf("force(%d,%d), dist = %f -> %f %f %f\n", p1.p.id, p2.p.id, distSqr, force[0], force[1], force[2]);
             if (p1.p.id == 0) {
                printf("sum add force Particle 0 = %f %f %f\n", p1.f.f[0], p1.f.f[1], p1.f.f[0]);
             }
             if (p2.p.id == 0) {
                printf("sum sub force Particle 0 = %f %f %f\n", p2.f.f[0], p2.f.f[1], p2.f.f[0]);
             }
#endif
          }
        }
      }

      real getMaxCutoff() {

        real cutoff = 0.0;

        for (int i = 0; i < ntypes; i++) {
          for (int j = 0; i < ntypes; i++) {
             cutoff = std::max(cutoff, parameterArray[i][j].getCutoff());
          }
        }
        return cutoff;cutoff;
      }

     protected:

      int ntypes;

      Parameters parameterArray[1][1];

    };

  }
}
#endif
