// ESPP_CLASS
#ifndef _INTERACTION_FIXEDTRIPLELISTINTERACTIONTEMPLATE_HPP
#define _INTERACTION_FIXEDTRIPLELISTINTERACTIONTEMPLATE_HPP

#include "mpi.hpp"
#include "types.hpp"
#include "Interaction.hpp"
#include "Real3D.hpp"
#include "Particle.hpp"
#include "FixedTripleList.hpp"
#include "esutil/Array2D.hpp"
#include "bc/BC.hpp"
#include "SystemAccess.hpp"

namespace espresso {
  namespace interaction {
    template < typename _AngularPotential >
    class FixedTripleListInteractionTemplate : public Interaction, SystemAccess {
    protected:
      typedef _AngularPotential Potential;
    public:
      FixedTripleListInteractionTemplate
      (shared_ptr < System > _system,
       shared_ptr < FixedTripleList > _fixedtripleList)
        : SystemAccess(_system), fixedtripleList(_fixedtripleList) 
      {
        potentialArray = esutil::Array2D<Potential, esutil::enlarge>(0, 0, Potential());
      }

      void
      setFixedTripleList(shared_ptr < FixedTripleList > _fixedtripleList) {
        fixedtripleList = _fixedtripleList;
      }

      shared_ptr < FixedTripleList > getFixedTripleList() {
        return fixedtripleList;
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
      shared_ptr < FixedTripleList > fixedtripleList;
      esutil::Array2D<Potential, esutil::enlarge> potentialArray;
    };

    //////////////////////////////////////////////////
    // INLINE IMPLEMENTATION
    //////////////////////////////////////////////////
    template < typename _AngularPotential > inline void
    FixedTripleListInteractionTemplate < _AngularPotential >::addForces() {
      LOG4ESPP_INFO(theLogger, "add forces computed by FixedTripleList");
      for (FixedTripleList::Iterator it(*fixedtripleList); it.isValid(); ++it) {
        Particle &p1 = *it->first;
        Particle &p2 = *it->second;
        Particle &p3 = *it->third;
        int type1 = p1.p.type;
        int type2 = p2.p.type;
        const Potential &potential = getPotential(type1, type2);

	real force12[3];
	real force32[3];
        Real3D dist12 = getSystemRef().bc->getMinimumImageVector(p1.r.p, p2.r.p);
        Real3D dist32 = getSystemRef().bc->getMinimumImageVector(p3.r.p, p2.r.p);
	potential._computeForce(force12, force32, dist12.get(), dist32.get());
	for(int k = 0; k < 3; k++) {
	  p1.f.f[k] += force12[k];
	  p2.f.f[k] -= force12[k] + force32[k];
	  p3.f.f[k] += force32[k];
	}
      }
    }

    template < typename _AngularPotential >
    inline real
    FixedTripleListInteractionTemplate < _AngularPotential >::
    computeEnergy() {
      LOG4ESPP_INFO(theLogger, "compute energy of the triples");

      real e = 0.0;
      for (FixedTripleList::Iterator it(*fixedtripleList); it.isValid(); ++it) {
        Particle &p1 = *it->first;
        Particle &p2 = *it->second;
        Particle &p3 = *it->third;
        int type1 = p1.p.type;
        int type2 = p2.p.type;
        const Potential &potential = getPotential(type1, type2);
        Real3D dist12 = getSystemRef().bc->getMinimumImageVector(p1.r.p, p2.r.p);
        Real3D dist32 = getSystemRef().bc->getMinimumImageVector(p3.r.p, p2.r.p);
        e += potential._computeEnergy(dist12, dist32);
      }
      real esum;
      boost::mpi::reduce(*mpiWorld, e, esum, std::plus<real>(), 0);
      return esum;
    }

    template < typename _AngularPotential >
    inline real
    FixedTripleListInteractionTemplate < _AngularPotential >::
    computeVirial() {
      LOG4ESPP_INFO(theLogger, "compute scalar virial of the triples");

      real w = 0.0;
      for (FixedTripleList::Iterator it(*fixedtripleList); it.isValid(); ++it) {
        Particle &p1 = *it->first;
        Particle &p2 = *it->second;
        Particle &p3 = *it->third;
        int type1 = p1.p.type;
        int type2 = p2.p.type;
        const Potential &potential = getPotential(type1, type2);
        real force12[3];
        real force32[3];
        Real3D dist12 = getSystemRef().bc->getMinimumImageVector(p1.r.p, p2.r.p);
        Real3D dist32 = getSystemRef().bc->getMinimumImageVector(p3.r.p, p2.r.p);
        potential._computeForce(force12, force32, dist12.get(), dist32.get());
        w += dist12 * force12 + dist32 * force32;
      }
      return w;
    }

    template < typename _AngularPotential >
    inline void
    FixedTripleListInteractionTemplate < _AngularPotential >::
    computeVirialTensor(real* wij_) {
      LOG4ESPP_INFO(theLogger, "compute the virial tensor of the triples");

      /* 
      wij_[0] = 0.0;
      wij_[1] = 0.0;
      wij_[2] = 0.0;
      wij_[3] = 0.0;
      wij_[4] = 0.0;
      wij_[5] = 0.0;
      */
     
      for (FixedTripleList::Iterator it(*fixedtripleList); it.isValid(); ++it) {
        Particle &p1 = *it->first;
        Particle &p2 = *it->second;
        Particle &p3 = *it->third;
        int type1 = p1.p.type;
        int type2 = p2.p.type;
        const Potential &potential = getPotential(type1, type2);
        real force12[3];
        real force32[3];
        Real3D dist12 = getSystemRef().bc->getMinimumImageVector(p1.r.p, p2.r.p);
        Real3D dist32 = getSystemRef().bc->getMinimumImageVector(p3.r.p, p2.r.p);
        potential._computeForce(force12, force32, dist12.get(), dist32.get());
        wij_[0] += dist12[0] * force12[0] + dist32[0] * force32[0];
        wij_[1] += dist12[1] * force12[1] + dist32[1] * force32[1];
        wij_[2] += dist12[2] * force12[2] + dist32[2] * force32[2];
        wij_[3] += dist12[0] * force12[1] + dist32[0] * force32[1];
        wij_[4] += dist12[0] * force12[2] + dist32[0] * force32[2];
        wij_[5] += dist12[1] * force12[2] + dist32[1] * force32[2];
      }
    }

    template < typename _AngularPotential >
    inline real
    FixedTripleListInteractionTemplate< _AngularPotential >::
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
