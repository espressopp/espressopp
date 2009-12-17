#ifndef _INTERACTION_INTERACTIONTEMPLATE_HPP
#define _INTERACTION_INTERACTIONTEMPLATE_HPP

#include "Interaction.hpp"

namespace espresso {
  namespace interaction {
    /** InteractionTemplate provides loop templates to compute
	forces and energies of the various interactions. */

    template< typename Derived >
    class InteractionTemplate : public Interaction {
    public:
      //////////////////////
      // STORAGE LOOP
      //////////////////////
      // ENERGY COMPUTATION
      // full square over two particle lists
      template < class ParticleList >
      real 
      computeCellEnergy(ParticleList &pl1, ParticleList &pl2) {
	real e = 0.0;
	for (int i = 0, endi = pl1.size(); i < endi; i++)
	  for (int j = 0, endj = pl2.size(); j < endj; j++) {
	    Particle &p1 = pl1[i];
	    Particle &p2 = pl2[j];
	    real dist[3];
	    real distSqr = esutil.getDist(dist, p1.r.p, p2.r.p);
	    e += computeEnergy(p1, p2, dist, distSqr);
	  }
	return e;
      }

      // half square over a single particle list
      template < class ParticleList >
      real 
      computeCellEnergy(ParticleList &pl) {
	real e = 0.0;
	for (int i = 0; i < pl.size(); i++)
	  for (int j = 0; j < i; j++) {
	    Particle &p1 = pl[i];
	    Particle &p2 = pl[j];
	    real dist[3];
	    real distSqr = esutil.getDist(dist, p1.r.p, p2.r.p);
	    e += computeEnergy(p1, p2, dist, distSqr);
	  }
	return e;
      }

      // full loop over a storage
      virtual real 
      computeEnergy(Storage::SelfPtr storage) {
	real e = 0.0;

	CellList &realCells = storage.getRealCells();
	for (CellList::const_iterator 
	       realCell = realCells.begin(),
	       realCellEnd = realCells.end();
	     realCell != realCellEnd;
	     realCell++) {
	  CellList &neighborCells = (*realCell)->neighborCells

	for (int i = 0; i < realCells.size(); i++) {
	  CellList &
	  for (int j = 0;
	}
      }

      // FORCE COMPUTATION
      // full square over two cells
      virtual void 
      computeCellForces(BoundaryConditions &bc,
			Cell &cell1, Cell &cell2) {
	for (int i = 0, endi = cell1.size(); i < endi; i++)
	  for (int j = 0, endj = cell2.size(); j < endj; j++) {
	    Particle &p1 = cell1[i];
	    Particle &p2 = cell2[j];
	    real dist[3];
	    real distSqr;
	    real force[3] = {0.0, 0.0, 0.0};
	    bc.getMinimumImageVector(dist, &distSqr, p1.r.p, p2.r.p);
	    computeForce(p1, p2, dist, distSqr, force);

	    for (int k = 0; k < 3; i++) {
	      p1.f.f[k] +=  force[k];
	      p2.f.f[k] += -force[k];
	    }
	  }
      }
   
      // half square over a single cell
      virtual void
      computeCellForces(BoundaryConditions &bc,
			Cell &cell) {
	for (int i = 0, endi = cell.size(); i < endi; i++)
	  for (int j = 0; j < i; j++) {
	    Particle &p1 = cell[i];
	    Particle &p2 = cell[j];
	    real dist[3];
	    real distSqr;
	    real force[3] = {0.0, 0.0, 0.0};
	    bc.getMinimumImageVector(dist, &distSqr, p1.r.p, p2.r.p);
	    computeForce(p1, p2, dist, distSqr, force);

	    for (int k = 0; k < 3; i++) {
	      p1.f.f[k] +=  force[k];
	      p2.f.f[k] += -force[k];
	    }
	  }
      }

      virtual void addForcesStorage() {
      
      }

      virtual void addForcesVerletLists() {
      }

      virtual real computeEnergyVerletLists() {
      }
    };

  }
}
#endif
