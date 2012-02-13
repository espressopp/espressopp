// ESPP_CLASS
#ifndef _INTERACTION_CELLLISTALLPAIRSINTERACTIONTEMPLATE_HPP
#define _INTERACTION_CELLLISTALLPAIRSINTERACTIONTEMPLATE_HPP

#include "types.hpp"
#include "Tensor.hpp"
#include "Interaction.hpp"
#include "storage/Storage.hpp"
#include "esutil/Array2D.hpp"
#include "iterator/CellListAllPairsIterator.hpp"

namespace espresso {
  namespace interaction {
    template < typename _Potential >
    class CellListAllPairsInteractionTemplate
        : public Interaction {
    protected:
      typedef _Potential Potential;
    public:
      CellListAllPairsInteractionTemplate
      (shared_ptr < storage::Storage > _storage)
        : storage(_storage) 
      {
        potentialArray = esutil::Array2D<Potential, esutil::enlarge>(0, 0, Potential());
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
      virtual void computeVirialTensor(Tensor& wij);
      virtual real getMaxCutoff();
      virtual int bondType() { return Nonbonded; }

    protected:
      int ntypes;
      esutil::Array2D< Potential, esutil::enlarge > potentialArray;
      shared_ptr< storage::Storage > storage;
    };

    //////////////////////////////////////////////////
    // INLINE IMPLEMENTATION
    //////////////////////////////////////////////////
    template < typename _Potential > inline void
    CellListAllPairsInteractionTemplate < _Potential >::
    addForces() {
      LOG4ESPP_INFO(theLogger, "add forces computed for all pairs in the cell lists");

      for (iterator::CellListAllPairsIterator it(storage->getRealCells()); it.isValid(); ++it) {
        Particle &p1 = *it->first;
        Particle &p2 = *it->second;
        const Potential &potential = getPotential(p1.type(), p2.type());

        Real3D force(0.0, 0.0, 0.0);
        if(potential._computeForce(force, p1, p2)) {
          p1.force() += force;
          p2.force() -= force;
        }
      }
    }

    template < typename _Potential > inline real 
    CellListAllPairsInteractionTemplate < _Potential >::
    computeVirial() {
      LOG4ESPP_INFO(theLogger, "computed virial for all pairs in the cell lists");
     
      real w = 0.0;
      for (iterator::CellListAllPairsIterator it(storage->getRealCells());
           it.isValid(); ++it) {
        Particle &p1 = *it->first;
        Particle &p2 = *it->second;
        int type1 = p1.type();
        int type2 = p2.type();
        const Potential &potential = getPotential(type1, type2);

        Real3D force(0.0, 0.0, 0.0);
        if(potential._computeForce(force, p1, p2)) { 
          Real3D dist = p1.position() - p2.position();
          w = w + dist * force;
        }
      }
      return w;
    }

    template < typename _Potential > inline void
    CellListAllPairsInteractionTemplate < _Potential >::
    computeVirialTensor(Tensor& wij) {
      LOG4ESPP_INFO(theLogger, "computed virial tensor for all pairs in the cell lists");

      for (iterator::CellListAllPairsIterator it(storage->getRealCells());
           it.isValid(); ++it) {
        const Particle &p1 = *it->first;
        const Particle &p2 = *it->second;
        const Potential &potential = getPotential(p1.type(), p2.type());

        Real3D force(0.0, 0.0, 0.0);
        if(potential._computeForce(force, p1, p2)) {
          Real3D dist = p1.position() - p2.position();
          wij += Tensor(dist, force);
        }
      }
    }

    template < typename _Potential >
    inline real
    CellListAllPairsInteractionTemplate < _Potential >::
    computeEnergy() {
      LOG4ESPP_INFO(theLogger, "compute energy by the Verlet List");

      real e = 0.0;
      for (iterator::CellListAllPairsIterator it(storage->getRealCells());
	   it.isValid(); ++it) {
        const Particle &p1 = *it->first;
        const Particle &p2 = *it->second;
        const Potential &potential = getPotential(p1.type(), p2.type());
	e += potential._computeEnergy(p1, p2);
      }
      return e;
    }

    template < typename _Potential >
    inline real
    CellListAllPairsInteractionTemplate < _Potential >::getMaxCutoff() {
      real cutoff = 0.0;
      for(int i = 0; i < ntypes; i++) {
        for(int j = 0; j < ntypes; j++) {
          cutoff = std::max(cutoff, getPotential(i, j).getCutoff());
        }
      }
      return cutoff;
    }
  }
}


#endif
