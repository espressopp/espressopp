// ESPP_CLASS
#ifndef _INTERACTION_FIXEDSINGLELISTINTERACTIONTEMPLATE_HPP
#define _INTERACTION_FIXEDSINGLELISTINTERACTIONTEMPLATE_HPP

#include "mpi.hpp"
#include "Interaction.hpp"
#include "Real3D.hpp"
#include "Tensor.hpp"
#include "Particle.hpp"
#include "FixedSingleList.hpp"
#include "esutil/Array2D.hpp"
#include "bc/BC.hpp"
#include "SystemAccess.hpp"
#include "Interaction.hpp"
#include "types.hpp"

namespace espresso {
  namespace interaction {
    template < typename _Potential >
    class FixedSingleListInteractionTemplate: public Interaction, SystemAccess {
        
    protected:
      typedef _Potential Potential;
      
    public:
      FixedSingleListInteractionTemplate
      (shared_ptr < System > system,
       shared_ptr < FixedSingleList > _fixedsingleList,
       shared_ptr < Potential > _potential)
        : SystemAccess(system), fixedsingleList(_fixedsingleList),
          potential(_potential)
      {
        if (! potential) {
          LOG4ESPP_ERROR(theLogger, "NULL potential");
        }
      }

      virtual ~FixedSingleListInteractionTemplate() {};

      void
      setFixedSingleList(shared_ptr < FixedSingleList > _fixedsingleList) {
        fixedsingleList = _fixedsingleList;
      }

      shared_ptr < FixedSingleList > getFixedSingleList() {
        return fixedsingleList;
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
      virtual void computeVirialTensor(Tensor& w, real z);
      virtual void computeVirialTensor(Tensor *w, int n);
      virtual real getMaxCutoff();
      virtual int bondType() { return Single; }

    protected:
      int ntypes;
      shared_ptr < FixedSingleList > fixedsingleList;
      shared_ptr < Potential > potential;
    };

    //////////////////////////////////////////////////
    // INLINE IMPLEMENTATION
    //////////////////////////////////////////////////
    template < typename _Potential > inline void
    FixedSingleListInteractionTemplate < _Potential >::addForces() {
      LOG4ESPP_INFO(theLogger, "add forces computed by the FixedSingle List");
      const bc::BC& bc = *getSystemRef().bc;  // boundary conditions
      real dist_dummy = 1.0;
      Real3D force;
      for (FixedSingleList::SingleList::Iterator it(*fixedsingleList); it.isValid(); ++it) {
        Particle &p1 = *it->first;
        if(potential->_computeForce(force, dist_dummy)) {
          p1.force() += force;
        }
      }
    }
    
    template < typename _Potential > inline real
    FixedSingleListInteractionTemplate < _Potential >::
    computeEnergy() {

      LOG4ESPP_INFO(theLogger, "compute energy of the FixedSingle list");

      real e = 0.0;
      const bc::BC& bc = *getSystemRef().bc;  // boundary conditions
      for (FixedSingleList::SingleList::Iterator it(*fixedsingleList);
	   it.isValid(); ++it) {
        const Particle &p1 = *it->first;
        Real3D r21_dummy = 1.0;
        e += potential->_computeEnergy(r21_dummy);
      }
      real esum;
      boost::mpi::all_reduce(*mpiWorld, e, esum, std::plus<real>());
      return esum;
    }

    template < typename _Potential > inline real
    FixedSingleListInteractionTemplate < _Potential >::
    computeVirial() {
      return 0.0;
    }

    template < typename _Potential > inline void
    FixedSingleListInteractionTemplate < _Potential >::computeVirialTensor(Tensor& w){
      LOG4ESPP_INFO(theLogger, "compute the virial tensor for the FixedPair List");
    }

    template < typename _Potential > inline void
    FixedSingleListInteractionTemplate < _Potential >::
    computeVirialTensor(Tensor& w, real z){
      LOG4ESPP_INFO(theLogger, "compute the virial tensor for the FixedPair List");
    }
    
    template < typename _Potential > inline void
    FixedSingleListInteractionTemplate < _Potential >::
    computeVirialTensor(Tensor *w, int n){
      LOG4ESPP_INFO(theLogger, "compute the virial tensor for the FixedPair List");
    }
    
    template < typename _Potential >
    inline real
    FixedSingleListInteractionTemplate< _Potential >::
    getMaxCutoff() {
      return potential->getCutoff();
    }
  }
}
#endif
