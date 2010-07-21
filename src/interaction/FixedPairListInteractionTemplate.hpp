// ESPP_CLASS
#ifndef _INTERACTION_FIXEDPAIRLISTINTERACTIONTEMPLATE_HPP
#define _INTERACTION_FIXEDPAIRLISTINTERACTIONTEMPLATE_HPP

#include "types.hpp"
#include "Interaction.hpp"
#include "Real3D.hpp"
#include "Particle.hpp"
#include "FixedPairList.hpp"
#include "esutil/Array2D.hpp"

namespace espresso {
  namespace interaction {
    template < typename _Potential >
    class FixedPairListInteractionTemplate
        : public Interaction {
    protected:
      typedef _Potential Potential;
    public:
      FixedPairListInteractionTemplate
      (shared_ptr < FixedPairList > _fixedpairList)
        : fixedpairList(_fixedpairList) 
      {
        potentialArray = esutil::Array2D<Potential, esutil::enlarge>(0, 0, Potential());
      }

      void
      setFixedPairList(shared_ptr < FixedPairList > _fixedpairList) {
        fixedpairList = _fixedpairList;
      }

      shared_ptr < FixedPairList > getFixedPairList() {
        return fixedpairList;
      }

      void
      setPotential(int type1, int type2, const Potential &potential) {
	potentialArray.at(type1, type2) = potential;
      }

      Potential &getPotential(int type1, int type2) {
        return potentialArray(type1, type2);
      }

      virtual void addForces();
      virtual real computeEnergy();
      virtual real computeVirial();
      virtual void computeVirialTensor(real* wij_);
      virtual real getMaxCutoff();

    protected:
      int ntypes;
      shared_ptr < FixedPairList > fixedpairList;
      esutil::Array2D<Potential, esutil::enlarge> potentialArray;
    };

    //////////////////////////////////////////////////
    // INLINE IMPLEMENTATION
    //////////////////////////////////////////////////
    template < typename _Potential > inline void
    FixedPairListInteractionTemplate < _Potential >::addForces() {
      LOG4ESPP_INFO(theLogger, "add forces computed by the FixedPair List");
      for (FixedPairList::Iterator it(*fixedpairList); it.isValid(); ++it) {
        Particle &p1 = *it->first;
        Particle &p2 = *it->second;
        int type1 = p1.p.type;
        int type2 = p2.p.type;
        const Potential &potential = getPotential(type1, type2);

	Real3D force;
	if (potential._computeForce(force, p1, p2))
	  for(int k = 0; k < 3; k++) {
	    p1.f.f[k] += force[k];
	    p2.f.f[k] -= force[k];
	  }
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
    
    template < typename _Potential >
    inline real
    FixedPairListInteractionTemplate < _Potential >::
    computeEnergy() {
      LOG4ESPP_INFO(theLogger, "compute energy of the FixedPair list pairs");

      real e = 0.0;
      for (FixedPairList::Iterator it(*fixedpairList); 
	   it.isValid(); ++it) {
        Particle &p1 = *it->first;
        Particle &p2 = *it->second;
        int type1 = p1.p.type;
        int type2 = p2.p.type;
        const Potential &potential = getPotential(type1, type2);
	e += potential._computeEnergy(p1, p2);
      }
      return e;
    }

    template < typename _Potential > inline real
    FixedPairListInteractionTemplate < _Potential >::computeVirial() {
      LOG4ESPP_INFO(theLogger, "compute the virial for the FixedPair List");
      
      real w = 0.0;
      for (FixedPairList::Iterator it(*fixedpairList);                
           it.isValid(); ++it) {                                         
        Particle &p1 = *it->first;                                       
        Particle &p2 = *it->second;                                      
        int type1 = p1.p.type;                                           
        int type2 = p2.p.type;
        const Potential &potential = getPotential(type1, type2);

        Real3D force;
        if (potential._computeForce(force, p1, p2)) {
          Real3D dist = Real3DRef(p1.r.p) - Real3DRef(p2.r.p);
          w = w + dist * force;
        }
      }
      return w; 
    }

    template < typename _Potential > inline void
    FixedPairListInteractionTemplate < _Potential >::computeVirialTensor(real* wij_) {
      LOG4ESPP_INFO(theLogger, "compute the virial tensor for the FixedPair List");

      wij_[0] = 0.0;
      wij_[1] = 0.0;
      wij_[2] = 0.0;
      wij_[3] = 0.0;
      wij_[4] = 0.0;
      wij_[5] = 0.0;
      for (FixedPairList::Iterator it(*fixedpairList);
           it.isValid(); ++it) {
        Particle &p1 = *it->first;
        Particle &p2 = *it->second;
        int type1 = p1.p.type;
        int type2 = p2.p.type;
        const Potential &potential = getPotential(type1, type2);

        Real3D force;
        if (potential._computeForce(force, p1, p2)) {
          Real3D dist = Real3DRef(p1.r.p) - Real3DRef(p2.r.p);
          wij_[0] += dist[0] * force[0];
          wij_[1] += dist[1] * force[1];
          wij_[2] += dist[2] * force[2];
          wij_[3] += dist[0] * force[1];
          wij_[4] += dist[0] * force[2];
          wij_[5] += dist[1] * force[2];
        }
      }
    }
 
    template < typename _Potential >
    inline real
    FixedPairListInteractionTemplate< _Potential >::
    getMaxCutoff() {
      real cutoff = 0.0;

      for (int i = 0; i < ntypes; i++) {
        for (int j = 0; j < ntypes; j++) {
          cutoff = std::max(cutoff, getPotential(i, j).getCutoff());
        }
      }
      return cutoff;
    }
  }
}
#endif
