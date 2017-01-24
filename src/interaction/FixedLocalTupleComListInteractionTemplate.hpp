/*
  Copyright (C) 2016
      Max Planck Institute for Polymer Research
  
  This file is part of ESPResSo++.
  
  ESPResSo++ is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.
  
  ESPResSo++ is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.
  
  You should have received a copy of the GNU General Public License
  along with this program.  If not, see <http://www.gnu.org/licenses/>. 
*/

// ESPP_CLASS
#ifndef _INTERACTION_FIXEDLOCALTUPLECOMLISTINTERACTIONTEMPLATE_HPP
#define _INTERACTION_FIXEDLOCALTUPLECOMLISTINTERACTIONTEMPLATE_HPP

#include "mpi.hpp"
#include "Interaction.hpp"
#include "Real3D.hpp"
#include "Tensor.hpp"
#include "Particle.hpp"
#include "FixedLocalTupleList.hpp"
#include "esutil/Array2D.hpp"
#include "bc/BC.hpp"
#include "SystemAccess.hpp"
#include "storage/Storage.hpp"
#include "Interaction.hpp"
#include "types.hpp"

namespace espressopp {
    namespace interaction {
	template < typename _Potential >
	class FixedLocalTupleComListInteractionTemplate: public Interaction, SystemAccess {
	    
	protected:
	    typedef _Potential Potential;
	    
	private:
	    long unsigned int num_of_part; // total number of particles
	    long unsigned int num_of_subchain; // total number of subchains
	    std::map<long unsigned int, Real3D> com_origin;
	    
	    int N_Constrain; // number of particles constraining the interaction
	    int nParticles;  // number of particles in this node
	    
	public:
	    FixedLocalTupleComListInteractionTemplate
	    (shared_ptr < System > _system,
	     shared_ptr < FixedLocalTupleList > _fixedtupleList,
	     shared_ptr < Potential > _potential)
	      : SystemAccess(_system), fixedtupleList(_fixedtupleList), potential(_potential) {
		if (! potential) {
		    LOG4ESPP_ERROR(theLogger, "NULL potential");
		}
		System& system = getSystemRef();
		
		long unsigned int localN = system.storage -> getNRealParticles();
		std::cout << "The number of real particles " << localN << endl;
		boost::mpi::all_reduce(*system.comm, localN, num_of_part, std::plus<int>());

		FixedLocalTupleList::TupleList::Iterator it(*fixedtupleList);
		std::vector<Particle*> pList = it->second;

		N_Constrain = pList.size() + 1;
		
		num_of_subchain = num_of_part/N_Constrain;
		std::cout << "The number of subchain " << num_of_subchain << endl;

		Real3D* totCOM;
		totCOM = new Real3D[num_of_subchain];
		Real3D* COM;
		COM = new Real3D[num_of_subchain];
		
		Real3D Li = system.bc->getBoxL();
		
		for(long unsigned int i = 0; i < num_of_subchain; i++) {
		    COM[i] = 0.;
		    totCOM[i] = 0.;
		}

		for (FixedLocalTupleList::TupleList::Iterator it(*fixedtupleList); it.isValid(); ++it) {
		    Particle *p = it->first;
		    Real3D pos_i = p->position();
		    Real3D tmp_com = pos_i;
		    std::vector<Particle*> pList = it->second;

		    for (long unsigned int j = 1; j < N_Constrain; j++) {
			p = pList[j - 1];
			Real3D pos_j = p->position();
			Real3D dist;
			system.bc->getMinimumImageVectorBox(dist, pos_j, pos_i);
			pos_j = pos_i + dist;
			tmp_com += pos_j;
			pos_i = pos_j;
		    }
		    COM[(p->id() - 1)/N_Constrain] = tmp_com/N_Constrain;
		}

		boost::mpi::all_reduce(*system.comm, COM, num_of_subchain, totCOM, std::plus<Real3D>() );
		
		for (long unsigned int i = 0; i < num_of_subchain; i++) {
		    com_origin[i] = totCOM[i];
		    /*
		    std::cout << "COM "<< i << ": "
			 << com_origin[i][0] << " "
			 << com_origin[i][1] << " "
			 << com_origin[i][2] << endl;
		    */
		}
		
		delete [] COM;
		COM = NULL;
		delete [] totCOM;
		totCOM = NULL;
	    }
	    
	    virtual ~FixedLocalTupleComListInteractionTemplate() {};
	    
	    // set the center of mass of subchains
	    void setCom(longint id, const Real3D& pos) {
		com_origin[id] = pos;
	    }

	    void
	    setFixedLocalTupleList(shared_ptr < FixedLocalTupleList > _fixedtupleList) {
		fixedtupleList = _fixedtupleList;
	    }
	    
	    shared_ptr < FixedLocalTupleList > getFixedLocalTupleList() {
		return fixedtupleList;
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
	    virtual real computeEnergyDeriv();
	    virtual real computeEnergyAA();
	    virtual real computeEnergyCG();      
	    virtual void computeVirialX(std::vector<real> &p_xx_total, int bins); 
	    virtual real computeVirial();
	    virtual void computeVirialTensor(Tensor& w);
	    virtual void computeVirialTensor(Tensor& w, real z);
	    virtual void computeVirialTensor(Tensor *w, int n);
	    virtual real getMaxCutoff() { return 0.0; }
	    virtual int bondType() { return Nonbonded; }
	    
	protected:
	    std::map<long unsigned int, Real3D> computeCOM() {

		std::map<long unsigned int, Real3D> com;

		System& system = getSystemRef();

		for (FixedLocalTupleList::TupleList::Iterator it(*fixedtupleList); it.isValid(); ++it) {
		    Particle *p = it->first;
		    Real3D pos_i = p->position();
		    Real3D tmp_com = pos_i;
		    std::vector<Particle*> pList = it->second;

		    for (long unsigned int j = 1; j < N_Constrain; j++) {
			p = pList[j - 1];
			Real3D pos_j = p->position();
			Real3D dist;
			system.bc->getMinimumImageVectorBox(dist, pos_j, pos_i);
			pos_j = pos_i + dist;
			tmp_com += pos_j;
			pos_i = pos_j;
		    }
		    com[(p->id() - 1)/N_Constrain] = tmp_com/N_Constrain;
		}
		
		return com;
	    }

	    int ntypes;
	    shared_ptr < FixedLocalTupleList > fixedtupleList;
	    shared_ptr < Potential > potential;
	};
	
	//////////////////////////////////////////////////
	// INLINE IMPLEMENTATION
	//////////////////////////////////////////////////
	template < typename _Potential > inline void
	FixedLocalTupleComListInteractionTemplate < _Potential >::addForces() {
	    LOG4ESPP_INFO(_Potential::theLogger, "adding forces of FixedLocalTupleComList");

	    std::map<long unsigned int, Real3D> com;
	    com = computeCOM();

	    for (FixedLocalTupleList::TupleList::Iterator it(*fixedtupleList); it.isValid(); ++it) {
		Particle *p = it->first;
		std::vector<Particle*> pList = it->second;
		Real3D diff;
		long unsigned int pid = p->id() - 1;
		getSystemRef().bc->getMinimumImageVectorBox(diff,
							    com_origin[pid/N_Constrain],
							    com[pid/N_Constrain]);

		p->force() += potential->_computeForce(diff, N_Constrain);
		
		for (int i = 0; i < N_Constrain - 1; i++) {
		    p = pList[i];
		    p->force() += potential->_computeForce(diff, N_Constrain);
		}
	    }
	}
	
	template < typename _Potential > inline real
	FixedLocalTupleComListInteractionTemplate < _Potential >::
	computeEnergy() {
	    
	    LOG4ESPP_INFO(theLogger, "compute energy of the FixedLocalTupleComList pairs");
	    
	    real e = 0.0;
	    std::map<long unsigned int, Real3D> com;
	    com = computeCOM();
	    
	    for (FixedLocalTupleList::TupleList::Iterator it(*fixedtupleList); it.isValid(); ++it) {
		Particle *p = it->first;
		Real3D diff;
		long unsigned int pid = p->id() - 1;
		getSystemRef().bc->getMinimumImageVectorBox(diff,
							    com_origin[pid/N_Constrain],
							    com[pid/N_Constrain]);
		
		e += N_Constrain*potential->_computeEnergy(diff);
	    }

	    real esum;
	    boost::mpi::all_reduce(*mpiWorld, e, esum, std::plus<real>());
	    return esum;
	}

	template < typename _Potential > inline real
	FixedLocalTupleComListInteractionTemplate < _Potential >::
	computeEnergyDeriv() {
	    std::cout << "Warning! At the moment computeEnergyDeriv() in FixedLocalTupleComListInteractionTemplate does not work." << std::endl;
	    return 0.0;
	}
	
	template < typename _Potential > inline real
	FixedLocalTupleComListInteractionTemplate < _Potential >::
	computeEnergyAA() {
	    std::cout << "Warning! At the moment computeEnergyAA() in FixedLocalTupleComListInteractionTemplate does not work." << std::endl;
	    return 0.0;
	}
	
	template < typename _Potential > inline real
	FixedLocalTupleComListInteractionTemplate < _Potential >::
	computeEnergyCG() {
	    std::cout << "Warning! At the moment computeEnergyCG() in FixedLocalTupleComListInteractionTemplate does not work." << std::endl;
	    return 0.0;
	}   
	
	template < typename _Potential >
	inline void
	FixedLocalTupleComListInteractionTemplate < _Potential >::
	computeVirialX(std::vector<real> &p_xx_total, int bins) {
	    LOG4ESPP_INFO(theLogger, "compute virial p_xx of the pressure tensor slabwise");
	}
	
	
	template < typename _Potential > inline real
	FixedLocalTupleComListInteractionTemplate < _Potential >::
	computeVirial() {
	    LOG4ESPP_INFO(theLogger, "compute the virial for the FixedLocalTupleCom List");
	    
	    real w = 0.0;
	    const bc::BC& bc = *getSystemRef().bc;  // boundary conditions
	    for (FixedLocalTupleList::TupleList::Iterator it(*fixedtupleList); it.isValid(); ++it) {
		const Particle &p1 = *it->first;
	    }
	    real wsum;
	    boost::mpi::all_reduce(*mpiWorld, w, wsum, std::plus<real>());
	    return wsum;
	}
	
	template < typename _Potential > inline void
	FixedLocalTupleComListInteractionTemplate < _Potential >::computeVirialTensor(Tensor& w){
	    LOG4ESPP_INFO(theLogger, "compute the virial tensor for the FixedLocalTupleCom List");
	    
	    Tensor wlocal(0.0);
	    const bc::BC& bc = *getSystemRef().bc;  // boundary conditions
	    for (FixedLocalTupleList::TupleList::Iterator it(*fixedtupleList); it.isValid(); ++it) {
		const Particle &p1 = *it->first;
	    }
	    
	    // reduce over all CPUs
	    Tensor wsum(0.0);
	    boost::mpi::all_reduce(*mpiWorld, (double*)&wlocal, 6, (double*)&wsum, std::plus<double>());
	    w += wsum;
	}
	
	template < typename _Potential > inline void
	FixedLocalTupleComListInteractionTemplate < _Potential >::
	computeVirialTensor(Tensor& w, real z){
	    std::cout<<"Warning! Calculating virial layerwise is not supported for "
		"long range interactions."<<std::endl;
	}
	
	template < typename _Potential > inline void
	FixedLocalTupleComListInteractionTemplate < _Potential >::
	computeVirialTensor(Tensor *w, int n){
	    std::cout<<"Warning! Calculating virial layerwise is not supported for "
		"long range interactions."<<std::endl;
	}
    }
}
#endif
