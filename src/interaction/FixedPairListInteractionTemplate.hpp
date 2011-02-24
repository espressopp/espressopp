// ESPP_CLASS
#ifndef _INTERACTION_FIXEDPAIRLISTINTERACTIONTEMPLATE_HPP
#define _INTERACTION_FIXEDPAIRLISTINTERACTIONTEMPLATE_HPP

#include "mpi.hpp"
#include "types.hpp"
#include "Interaction.hpp"
#include "Real3D.hpp"
#include "Tensor.hpp"
#include "Particle.hpp"
#include "FixedPairList.hpp"
#include "esutil/Array2D.hpp"
#include "bc/BC.hpp"
#include "SystemAccess.hpp"

namespace espresso {
  namespace interaction {
    template < typename _Potential >
    class FixedPairListInteractionTemplate: public Interaction, SystemAccess {
        
    protected:
      typedef _Potential Potential;
      
    public:
      FixedPairListInteractionTemplate
      (shared_ptr < System > system,
       shared_ptr < FixedPairList > _fixedpairList,
       shared_ptr < Potential > _potential)
        : SystemAccess(system), fixedpairList(_fixedpairList),
          potential(_potential)
      {
        if (! potential) {
          LOG4ESPP_ERROR(theLogger, "NULL potential");
        }
      }

      void
      setFixedPairList(shared_ptr < FixedPairList > _fixedpairList) {
        fixedpairList = _fixedpairList;
      }

      shared_ptr < FixedPairList > getFixedPairList() {
        return fixedpairList;
      }

      void
      setPotential(shared_ptr < Potential> _potential) {
        if (_potential) {
          potential = _potential;
        } else {
          LOG4ESPP_ERROR(theLogger, "NULL potential");
        }
      }

      shared_ptr < Potential > getPotential() {
        return potential;
      }

      virtual void addForces();
      virtual real computeEnergy();
      virtual real computeVirial();
      virtual void computeVirialTensor(Tensor& w);
      virtual real getMaxCutoff();

    protected:
      int ntypes;
      shared_ptr < FixedPairList > fixedpairList;
      shared_ptr < Potential > potential;
    };

    //////////////////////////////////////////////////
    // INLINE IMPLEMENTATION
    //////////////////////////////////////////////////
    template < typename _Potential > inline void
    FixedPairListInteractionTemplate < _Potential >::addForces() {
      LOG4ESPP_INFO(theLogger, "add forces computed by the FixedPair List");
      const bc::BC& bc = *getSystemRef().bc;  // boundary conditions
      for (FixedPairList::Iterator it(*fixedpairList); it.isValid(); ++it) {
        Particle &p1 = *it->first;
        Particle &p2 = *it->second;
        Real3D dist;
        bc.getMinimumImageVectorBox(dist, p1.position(), p2.position());
        Real3D force;
        if(potential->_computeForce(force, dist)) {
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
      const bc::BC& bc = *getSystemRef().bc;  // boundary conditions
      for (FixedPairList::Iterator it(*fixedpairList); 
	   it.isValid(); ++it) {
        const Particle &p1 = *it->first;
        const Particle &p2 = *it->second;
        Real3D dist;
        bc.getMinimumImageVectorBox(dist, p1.position(), p2.position());
        e += potential->_computeEnergy(dist);
      }
      real esum;
      boost::mpi::reduce(*mpiWorld, e, esum, std::plus<real>(), 0);
      return esum;
    }

    template < typename _Potential > inline real
    FixedPairListInteractionTemplate < _Potential >::
    computeVirial() {
      LOG4ESPP_INFO(theLogger, "compute the virial for the FixedPair List");
      
      real w = 0.0;
      const bc::BC& bc = *getSystemRef().bc;  // boundary conditions
      for (FixedPairList::Iterator it(*fixedpairList);                
           it.isValid(); ++it) {                                         
        const Particle &p1 = *it->first;                                       
        const Particle &p2 = *it->second;                                      

        Real3D dist;
        bc.getMinimumImageVectorBox(dist, p1.position(), p2.position());
        Real3D force;
        if(potential->_computeForce(force, dist)) {
          w += dist * force;
        }
      }
      return w; 
    }

    template < typename _Potential > inline void
    FixedPairListInteractionTemplate < _Potential >::computeVirialTensor(Tensor& w) {
      LOG4ESPP_INFO(theLogger, "compute the virial tensor for the FixedPair List");
      const bc::BC& bc = *getSystemRef().bc;  // boundary conditions
      for (FixedPairList::Iterator it(*fixedpairList);
           it.isValid(); ++it) {
        const Particle &p1 = *it->first;
        const Particle &p2 = *it->second;
        Real3D dist;
        bc.getMinimumImageVectorBox(dist, p1.position(), p2.position());
        Real3D force;
        if(potential->_computeForce(force, dist)) { 
          w += Tensor(dist, force);
        }
      }
    }
 
    template < typename _Potential >
    inline real
    FixedPairListInteractionTemplate< _Potential >::getMaxCutoff() {
      return potential->getCutoff();
    }
  }
}
#endif
