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

    template< typename Derived >
    class InteractionTemplate : public Interaction {
    public:

      // ENERGY COMPUTATION
      // full square over two particle lists

      virtual real computeCellEnergy(ParticleList &pl1, ParticleList &pl2)
      {
        LOG4ESPP_INFO(theLogger, "compute energy for a cell pair");

        Derived* interaction = static_cast<Derived*>(this);
	real e = 0.0;
	for (int i = 0, endi = pl1.size(); i < endi; i++) {
	  Particle &p1 = pl1[i];
          int type1 = p1.p.type;
	  for (int j = 0, endj = pl2.size(); j < endj; j++) {
	    Particle &p2 = pl2[j];
	    real dist[3];
	    real distSqr = distance2vec(p1.r.p, p2.r.p, dist);
            int type2 = p2.p.type;
            if (distSqr < interaction->getCutoffSqr(type1, type2)) {
	       e += interaction->computeEnergy(p1, p2, type1, type2, dist, distSqr);
            }
	  }
        }
	return e;
      }

      /** Reimplementation of pure routine for the derived class; the
          methods getCutoffSqr and computeEnergy will be inlined for efficiency.
      */

      virtual real computeVerletListEnergy(shared_ptr<VerletList> vl) {

        LOG4ESPP_INFO(theLogger, "compute energy by the Verlet List");

        Derived* interaction = static_cast<Derived*>(this);

        const VerletList::PairList pairList = vl->getPairs();
        int npairs = pairList.size();
        real e = 0.0;
        for (int i = 0; i < npairs; i++) {
          Particle& p1 = *pairList[i].first;
          Particle& p2 = *pairList[i].second;
          int type1 = p1.p.type;
          real dist[3];
          real distSqr = distance2vec(p1.r.p, p2.r.p, dist);
          int type2 = p2.p.type;
          if (distSqr < interaction->getCutoffSqr(type1, type2)) {
             e += interaction->computeEnergy(p1, p2, type1, type2, dist, distSqr);
          }
        }
        return e;
      }

      // FORCE COMPUTATION
      // full square over two particle lists

      virtual void addForcesStorage() {
      }

      virtual void addVerletListForces(shared_ptr<VerletList> vl) {

       LOG4ESPP_INFO(theLogger, "add forces computed by the Verlet List");

       Derived* interaction = static_cast<Derived*>(this);

        const VerletList::PairList pairList = vl->getPairs();
        int npairs = pairList.size();
        real e = 0.0;
        for (int i = 0; i < npairs; i++) {
          Particle& p1 = *pairList[i].first;
          Particle& p2 = *pairList[i].second;
          int type1 = p1.p.type;
          real dist[3];
          real distSqr = distance2vec(p1.r.p, p2.r.p, dist);
          int type2 = p2.p.type;
          if (distSqr < interaction->getCutoffSqr(type1, type2)) {
             real force[3];
             interaction->computeForce(p1, p2, type1, type2, dist, distSqr, force);
             for (int k = 0; k < 3; k++) {
               p1.f.f[k] +=  force[k];
               p2.f.f[k] += -force[k];
             }
          }
        }
      }

    };

  }
}
#endif
