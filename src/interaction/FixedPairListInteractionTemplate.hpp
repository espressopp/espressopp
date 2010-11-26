// ESPP_CLASS
#ifndef _INTERACTION_FIXEDPAIRLISTINTERACTIONTEMPLATE_HPP
#define _INTERACTION_FIXEDPAIRLISTINTERACTIONTEMPLATE_HPP

#include "mpi.hpp"
#include "types.hpp"
#include "Interaction.hpp"
#include "Real3D.hpp"
#include "Particle.hpp"
#include "FixedPairList.hpp"
#include "esutil/Array2D.hpp"
#include "bc/BC.hpp"
#include "SystemAccess.hpp"

namespace espresso {
  namespace interaction {
    template < typename _Potential >
    class FixedPairListInteractionTemplate
       : public Interaction, SystemAccess {
    protected:
      typedef _Potential Potential;
    public:
      FixedPairListInteractionTemplate
      (shared_ptr < System > system,
       shared_ptr < FixedPairList > _fixedpairList)

        : SystemAccess(system), fixedpairList(_fixedpairList) 

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
        const Potential &potential = getPotential(p1.type(), p2.type());

	Real3D force(0.0, 0.0, 0.0);
        Real3D dist = getSystemRef().bc->getMinimumImageVector(p1.position(), 
                                                               p2.position());
	if(potential._computeForce(force, dist)) {
          p1.force() += force;
          p2.force() -= force;
        }
      }
    }
    
    template < typename _Potential >
    inline real
    FixedPairListInteractionTemplate < _Potential >::
    computeEnergy() {

      LOG4ESPP_INFO(theLogger, "compute energy of the FixedPair list pairs");

      real e = 0.0;
      for (FixedPairList::Iterator it(*fixedpairList); 
	   it.isValid(); ++it) {
        const Particle &p1 = *it->first;
        const Particle &p2 = *it->second;
        const Potential &potential = getPotential(p1.type(), p2.type());
        Real3D dist = getSystemRef().bc->getMinimumImageVector(p1.position(), p2.position());
	e += potential._computeEnergy(dist);
      }
      real esum;
      boost::mpi::reduce(*mpiWorld, e, esum, std::plus<real>(), 0);
      return esum;
    }

    template < typename _Potential > inline real
    FixedPairListInteractionTemplate < _Potential >::computeVirial() {
      LOG4ESPP_INFO(theLogger, "compute the virial for the FixedPair List");
      
      real w = 0.0;
      for (FixedPairList::Iterator it(*fixedpairList);                
           it.isValid(); ++it) {                                         
        const Particle &p1 = *it->first;                                       
        const Particle &p2 = *it->second;                                      
        const Potential &potential = getPotential(p1.type(), p2.type());

        Real3D force(0.0, 0.0, 0.0);
        Real3D dist = getSystemRef().bc->getMinimumImageVector(p1.position(), p2.position());
        if(potential._computeForce(force, dist)) {
          w += dist * force;
        }
      }
      return w; 
    }

    template < typename _Potential > inline void
    FixedPairListInteractionTemplate < _Potential >::computeVirialTensor(real* wij_) {
      LOG4ESPP_INFO(theLogger, "compute the virial tensor for the FixedPair List");

      for (FixedPairList::Iterator it(*fixedpairList);
           it.isValid(); ++it) {
        const Particle &p1 = *it->first;
        const Particle &p2 = *it->second;
        const Potential &potential = getPotential(p1.type(), p2.type());
        Real3D force(0.0, 0.0, 0.0);
        Real3D dist = getSystemRef().bc->getMinimumImageVector(p1.position(), p2.position());
        if(potential._computeForce(force, dist)) { 
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
    FixedPairListInteractionTemplate< _Potential >::getMaxCutoff() {
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
