#ifndef _INTERACTION_VERLETLISTINTERACTIONTEMPLATE_HPP
#define _INTERACTION_VERLETLISTINTERACTIONTEMPLATE_HPP

#include "Interaction.hpp"
#include "Real3D.hpp"
#include "Particle.hpp"
#include "VerletList.hpp"

namespace espresso {
  namespace interaction {
    template < typename _InteractionFunction >
    class VerletListInteractionTemplate
        : public Interaction {
    protected:
      typedef _InteractionFunction InteractionFunction;
    public:
      VerletListInteractionTemplate
      (shared_ptr < VerletList > _verletList)
        : verletList(_verletList) {
      }

      void
      setVerletList(shared_ptr < VerletList > _verletList) {
        verletList = _verletList;
      }

      shared_ptr < VerletList > getVerletList() {
        return verletList;
      }

      void
      setFunction(int type1, int type2, const InteractionFunction &func) {
        // TODO: automatically resize array
        functionArray[type1][type2] = func;
      }

      InteractionFunction &getFunction(int type1, int type2) {
        // TODO: automatically resize array
        return functionArray[type1][type2];
      }

      virtual void addForces();
      virtual real computeEnergy();
      virtual real getMaxCutoff();

    protected:
      int ntypes;
      shared_ptr < VerletList > verletList;
      InteractionFunction functionArray[1][1];
    };

//////////////////////////////////////////////////
// INLINE IMPLEMENTATION
//////////////////////////////////////////////////
    template < typename _InteractionFunction > inline void
    VerletListInteractionTemplate < _InteractionFunction >::addForces() {
      LOG4ESPP_INFO(theLogger, "add forces computed by the Verlet List");
      const PairList pairList = verletList->getPairs();
      Real3D dist;
      real distSqr;

      int npairs = pairList.size();
      for(int i = 0; i < npairs; i++) {
        Particle &p1 = *pairList[i].first;
        Particle &p2 = *pairList[i].second;
        int type1 = p1.p.type;
        int type2 = p2.p.type;
        dist = p1.r.p - p2.r.p;
        distSqr = dist.sqr();
        const InteractionFunction &func = getFunction(type1, type2);

        if(distSqr < func.getCutoffSqr()) {
          Real3D force;
          func.getForce(p1, p2, dist, distSqr, force);
          for(int k = 0; k < 3; k++) {
            p1.f.f[k] += force[k];
            p2.f.f[k] -= force[k];
          }
#if 0
          printf
          ("dist(%d,%d), dist = %f -> %f %f %f\n",
           p1.p.id, p2.p.id, distSqr, dist[0], dist[1], dist[2]);
          printf
          ("force(%d,%d), dist = %f -> %f %f %f\n",
           p1.p.id, p2.p.id, distSqr,
           force[0], force[1], force[2]);
          if(p1.p.id == 0) {
            printf
            ("sum add force Particle 0 = %f %f %f\n",
             p1.f.f[0], p1.f.f[1], p1.f.f[0]);
          }
          if(p2.p.id == 0) {
            printf
            ("sum sub force Particle 0 = %f %f %f\n",
             p2.f.f[0], p2.f.f[1], p2.f.f[0]);
          }
#endif
        }
      }
    }

    template < typename _InteractionFunction >
    inline real
    VerletListInteractionTemplate < _InteractionFunction >::
    computeEnergy() {
      LOG4ESPP_INFO(theLogger, "compute energy by the Verlet List");

      const PairList pairList = verletList->getPairs();
      Real3D dist;
      real distSqr;

      int npairs = pairList.size();
      real e = 0.0;
      for(int i = 0; i < npairs; i++) {
        Particle &p1 = *pairList[i].first;
        Particle &p2 = *pairList[i].second;
        int type1 = p1.p.type;
        int type2 = p2.p.type;
        dist = p1.r.p - p2.r.p;
        distSqr = dist.sqr();
        const InteractionFunction &func = getFunction(type1, type2);
        if(distSqr < func.getCutoffSqr()) {
          e += func.getEnergy(p1, p2, dist, distSqr);
        }
      }
      return e;
    }

    template < typename _InteractionFunction >
    inline real
    VerletListInteractionTemplate <
    _InteractionFunction >::getMaxCutoff() {
      real cutoff = 0.0;

      for(int i = 0; i < ntypes; i++) {
        for(int j = 0; i < ntypes; i++) {
          cutoff = std::max(cutoff, getFunction(i, j).getCutoff());
        }
      }
      return cutoff;
    }
  }
}
#endif
