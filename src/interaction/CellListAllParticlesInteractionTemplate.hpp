// ESPP_CLASS
#ifndef _INTERACTION_CELLLISTALLPARTICLESINTERACTIONTEMPLATE_HPP
#define _INTERACTION_CELLLISTALLPARTICLESINTERACTIONTEMPLATE_HPP

#include "types.hpp"
#include "Tensor.hpp"
#include "Interaction.hpp"
#include "storage/Storage.hpp"
#include "iterator/CellListIterator.hpp"

namespace espresso {
  namespace interaction {
    template < typename _Potential >
    class CellListAllParticlesInteractionTemplate
        : public Interaction {
    protected:
      typedef _Potential Potential;
    public:
      CellListAllParticlesInteractionTemplate
      (shared_ptr < storage::Storage > _storage, shared_ptr < Potential > _potential)
        : storage(_storage), potential(_potential) {}

      shared_ptr< Potential > getPotential() {
        return potential;
      }

      virtual void addForces();
      virtual real computeEnergy();
      virtual real computeVirial();
      virtual void computeVirialTensor(Tensor& wij);
      virtual real getMaxCutoff() { return 0.0; }
      virtual int bondType() { return Nonbonded; }

    protected:
      int ntypes;
      shared_ptr< storage::Storage > storage;
      shared_ptr< Potential > potential;
    };

    //////////////////////////////////////////////////
    // INLINE IMPLEMENTATION
    //////////////////////////////////////////////////
    template < typename _Potential > inline void
    CellListAllParticlesInteractionTemplate < _Potential >::
    addForces() {
      LOG4ESPP_INFO(theLogger, "add forces computed for all particles in the cell lists");

      potential->_computeForce(storage->getRealCells());
    }

    template < typename _Potential >
    inline real
    CellListAllParticlesInteractionTemplate < _Potential >::
    computeEnergy() {
      LOG4ESPP_INFO(theLogger, "compute energy for all particles in cell list");
      
      // for the long range interaction the energy is already reduced in _computeEnergy
      return potential->_computeEnergy(storage->getRealCells());
    }
    
    template < typename _Potential > inline real 
    CellListAllParticlesInteractionTemplate < _Potential >::
    computeVirial() {
      LOG4ESPP_INFO(theLogger, "computed virial for all particles in the cell lists");
     
      // for the long range interaction the virial is already reduced in _computeVirial
      return potential -> _computeVirial(storage->getRealCells());
    }

    template < typename _Potential > inline void
    CellListAllParticlesInteractionTemplate < _Potential >::
    computeVirialTensor(Tensor& wij) {
      LOG4ESPP_INFO(theLogger, "computed virial tensor for all particles in the cell lists");
      
      // for the long range interaction the virialTensor is already reduced in _computeVirialTensor
      wij += potential -> _computeVirialTensor(storage->getRealCells());
    }

  }
}


#endif
