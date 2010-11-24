// ESPP_CLASS
#ifndef _INTERACTION_CELLLISTALLPAIRSINTERACTIONTEMPLATE_HPP
#define _INTERACTION_CELLLISTALLPAIRSINTERACTIONTEMPLATE_HPP

#include "types.hpp"
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
      virtual void computeVirialTensor(real* wij_);
      virtual real getMaxCutoff();

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

      for (iterator::CellListAllPairsIterator it(storage->getRealCells());
	   it.isValid(); ++it) {
        Particle &p1 = *it->first;
        Particle &p2 = *it->second;
        int type1 = p1.p.type;
        int type2 = p2.p.type;
        const Potential &potential = getPotential(type1, type2);

	real force[3];
	potential._computeForce(force, p1, p2);
	  for (int k = 0; k < 3; k++) {
	    p1.f.f[k] += force[k];
	    p2.f.f[k] -= force[k];
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
        int type1 = p1.p.type;
        int type2 = p2.p.type;
        const Potential &potential = getPotential(type1, type2);

        real force[3];
        potential._computeForce(force, p1, p2); 
          Real3D dist = Real3DRef(p1.r.p) - Real3DRef(p2.r.p);
          w = w + dist * Real3DRef(force);
      }
      return w;
    }

    template < typename _Potential > inline void
    CellListAllPairsInteractionTemplate < _Potential >::
    computeVirialTensor(real* wij_) {
      LOG4ESPP_INFO(theLogger, "computed virial tensor for all pairs in the cell lists");

      /*
      wij_[0] = 0.0;
      wij_[1] = 0.0;
      wij_[2] = 0.0;
      wij_[3] = 0.0;
      wij_[4] = 0.0;
      wij_[5] = 0.0;
      */
      for (iterator::CellListAllPairsIterator it(storage->getRealCells());
           it.isValid(); ++it) {
        Particle &p1 = *it->first;
        Particle &p2 = *it->second;
        int type1 = p1.p.type;
        int type2 = p2.p.type;
        const Potential &potential = getPotential(type1, type2);

        real force[3];
        potential._computeForce(force, p1, p2); 
          Real3D dist = Real3DRef(p1.r.p) - Real3DRef(p2.r.p);
          wij_[0] += dist[0] * force[0];
          wij_[1] += dist[1] * force[1];
          wij_[2] += dist[2] * force[2];
          wij_[3] += dist[0] * force[1];
          wij_[4] += dist[0] * force[2];
          wij_[5] += dist[1] * force[2];
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
        Particle &p1 = *it->first;
        Particle &p2 = *it->second;
        int type1 = p1.p.type;
        int type2 = p2.p.type;
        const Potential &potential = getPotential(type1, type2);
	e += potential._computeEnergy(p1, p2);
      }
      return e;
    }

    template < typename _Potential >
    inline real
    CellListAllPairsInteractionTemplate < _Potential >::
    getMaxCutoff() {
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
