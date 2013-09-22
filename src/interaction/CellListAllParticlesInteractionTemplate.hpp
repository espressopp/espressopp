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
    class CellListAllParticlesInteractionTemplate: public Interaction {
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
      virtual real computeEnergyAA();
      virtual real computeEnergyCG();
      virtual void computeVirialX(std::vector<real> &p_xx_total, int bins); 
      virtual real computeVirial();
      virtual void computeVirialTensor(Tensor& wij);
      virtual void computeVirialTensor(Tensor& w, real z);
      virtual void computeVirialTensor(Tensor *w, int n);
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
    computeEnergyAA() {
      std::cout << "Warning! At the moment computeEnergyAA() in CellListAllParticlesInteractionTemplate does not work." << std::endl;
      return 0.0;
    }
    
    template < typename _Potential > inline real
    CellListAllParticlesInteractionTemplate < _Potential >::
    computeEnergyCG() {
      std::cout << "Warning! At the moment computeEnergyCG() in CellListAllParticlesInteractionTemplate does not work." << std::endl;
      return 0.0;
    }
    
    template < typename _Potential >
    inline void
    CellListAllParticlesInteractionTemplate < _Potential >::
    computeVirialX(std::vector<real> &p_xx_total, int bins) {
        std::cout << "Warning! At the moment computeVirialX in CellListAllParticlesInteractionTemplate does not work." << std::endl << "Therefore, the corresponding interactions won't be included in calculation." << std::endl;
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
    
    
    template < typename _Potential > inline void
    CellListAllParticlesInteractionTemplate < _Potential >::
    computeVirialTensor(Tensor& wij, real z) {
      std::cout<<"Warning! Calculating virial layerwise is not supported for "
              "long range interactions."<<std::endl;
    }
    
    template < typename _Potential > inline void
    CellListAllParticlesInteractionTemplate < _Potential >::
    computeVirialTensor(Tensor *wij, int n) {
      std::cout<<"Warning! Calculating virial layerwise is not supported for "
              "long range interactions."<<std::endl;
    }

  }
}


#endif
