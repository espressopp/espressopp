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
#ifndef _INTERACTION_FIXEDLOCALTUPLERGLISTINTERACTIONTEMPLATE_HPP
#define _INTERACTION_FIXEDLOCALTUPLERGLISTINTERACTIONTEMPLATE_HPP

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
	class FixedLocalTupleRgListInteractionTemplate: public Interaction, SystemAccess {
	    
	protected:
	    typedef _Potential Potential;
	    
	private:
	    long unsigned int num_of_part; // total number of particles
	    long unsigned int num_of_subchain; // total number of subchains
	    //std::map<long unsigned int, Real3D> com_origin;
	    std::map<long unsigned int, real> rg_origin;
	    
	    int N_Constrain; // number of particles constraining the interaction
	    
	public:
	    FixedLocalTupleRgListInteractionTemplate
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
		
		//for (long unsigned int i = 0; i < num_of_subchain; i++) {
		//    com_origin[i] = totCOM[i];
		//}

		real* totRG;
		totRG = new real[num_of_subchain];
		real* RG;
		RG = new real[num_of_subchain];
		
		for(long unsigned int i = 0; i < num_of_subchain; i++) {
		    RG[i] = 0.;
		    totRG[i] = 0.;
		}

		for (FixedLocalTupleList::TupleList::Iterator it(*fixedtupleList); it.isValid(); ++it) {
		    Particle *p = it->first;
		    Real3D pos_i = p->position();
		    std::vector<Particle*> pList = it->second;
		    long unsigned int id = (p->id() - 1)/N_Constrain;

		    Real3D dist;
		    system.bc->getMinimumImageVectorBox(dist, pos_i, totCOM[id]);
		    RG[id] = dist.sqr();
		    for (long unsigned int j = 1; j < N_Constrain; j++) {
			p = pList[j - 1];
			Real3D pos_j = p->position();
			system.bc->getMinimumImageVectorBox(dist, pos_j, pos_i);
			pos_j = pos_i + dist;
			system.bc->getMinimumImageVectorBox(dist, pos_j, totCOM[id]);
			RG[id] += dist.sqr();
			pos_i = pos_j;
		    }
		    RG[id] /= N_Constrain;
		}
		boost::mpi::all_reduce(*system.comm, RG, num_of_subchain, totRG, std::plus<real>() );
		
		if (system.comm->rank() == 0) {
		    for(long unsigned int i = 0; i < num_of_subchain; i++) {
			rg_origin[i] = totRG[i];
		    }
		    
		}

		delete [] RG;
		RG = NULL;
		delete [] totRG;
		totRG = NULL;
		
		delete [] COM;
		COM = NULL;
		delete [] totCOM;
		totCOM = NULL;
	    }
	    
	    virtual ~FixedLocalTupleRgListInteractionTemplate() {};
	    
	    // set the center of mass of subchains
	    void setRG(longint id, const real rg) {
		rg_origin[id] = rg;
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

	    // calcuclate the center of mass of subchain in the local cell list.
	    std::map<long unsigned int, real>
	    computeRG(std::map<long unsigned int, Real3D> com) {
		std::map<long unsigned int, real> rg;
		
		System& system = getSystemRef();


		for (FixedLocalTupleList::TupleList::Iterator it(*fixedtupleList); it.isValid(); ++it) {
		    Particle *p = it->first;
		    Real3D pos_i = p->position();
		    std::vector<Particle*> pList = it->second;
		    long unsigned int id = (p->id() - 1)/N_Constrain;

		    Real3D dist;
		    system.bc->getMinimumImageVectorBox(dist, pos_i, com[id]);
		    real tmp_rg = dist.sqr();
		    for (long unsigned int j = 1; j < N_Constrain; j++) {
			p = pList[j - 1];
			Real3D pos_j = p->position();
			system.bc->getMinimumImageVectorBox(dist, pos_j, pos_i);
			pos_j = pos_i + dist;
			system.bc->getMinimumImageVectorBox(dist, pos_j, com[id]);
			tmp_rg += dist.sqr();
			pos_i = pos_j;
		    }
		    rg[id] = tmp_rg/N_Constrain;
		}

		return rg;
	    }

	    int ntypes;
	    shared_ptr < FixedLocalTupleList > fixedtupleList;
	    shared_ptr < Potential > potential;
	};
	
	//////////////////////////////////////////////////
	// INLINE IMPLEMENTATION
	//////////////////////////////////////////////////
	template < typename _Potential > inline void
	FixedLocalTupleRgListInteractionTemplate < _Potential >::addForces() {
	    LOG4ESPP_INFO(_Potential::theLogger, "adding forces of FixedLocalTupleRgList");

	    System& system = getSystemRef();
	    
	    std::map<long unsigned int, Real3D> com;
	    com = computeCOM();

	    std::map<long unsigned int, real> rg_sq;
	    rg_sq = computeRG(com);

	    for (FixedLocalTupleList::TupleList::Iterator it(*fixedtupleList); it.isValid(); ++it) {
		Particle *p = it->first;
		std::vector<Particle*> pList = it->second;
		Real3D diff;
		long unsigned int id = (p->id() - 1)/N_Constrain;

		system.bc->getMinimumImageVectorBox(diff,
						    p->position(),
						    com[id]);

		real diff_rg = rg_origin[id] - rg_sq[id];
		
		p->force() += potential->_computeForce(diff, rg_sq[id], diff_rg, N_Constrain);

		for (long unsigned int i = 0; i < N_Constrain - 1; i++) {
		    Particle* pt = pList[i];
		    system.bc->getMinimumImageVectorBox(diff,
							pt->position(),
							com[id]);
		    pt->force() +=
			potential->_computeForce(diff, rg_sq[id], diff_rg, N_Constrain);
		}

	    }
	}
	
	template < typename _Potential > inline real
	FixedLocalTupleRgListInteractionTemplate < _Potential >::
	computeEnergy() {
	    
	    LOG4ESPP_INFO(theLogger, "compute energy of the FixedLocalTupleRgList pairs");
	    
	    real e = 0.0;
	    std::map<long unsigned int, Real3D> com;
	    com = computeCOM();

	    std::map<long unsigned int, real> rg;
	    rg = computeRG(com);
	    
	    for (FixedLocalTupleList::TupleList::Iterator it(*fixedtupleList); it.isValid(); ++it) {
		Particle *p = it->first;

		long unsigned int id = (p->id() - 1)/N_Constrain;
		real diff = rg_origin[id] - rg[id];
		
		e += N_Constrain*potential->_computeEnergy(diff);
	    }

	    real esum;
	    boost::mpi::all_reduce(*mpiWorld, e, esum, std::plus<real>());
	    return esum;
	}

	template < typename _Potential > inline real
	FixedLocalTupleRgListInteractionTemplate < _Potential >::
	computeEnergyDeriv() {
	    std::cout << "Warning! At the moment computeEnergyDeriv() in FixedLocalTupleRgListInteractionTemplate does not work." << std::endl;
	    return 0.0;
	}
	
	template < typename _Potential > inline real
	FixedLocalTupleRgListInteractionTemplate < _Potential >::
	computeEnergyAA() {
	    std::cout << "Warning! At the moment computeEnergyAA() in FixedLocalTupleRgListInteractionTemplate does not work." << std::endl;
	    return 0.0;
	}
	
	template < typename _Potential > inline real
	FixedLocalTupleRgListInteractionTemplate < _Potential >::
	computeEnergyCG() {
	    std::cout << "Warning! At the moment computeEnergyCG() in FixedLocalTupleRgListInteractionTemplate does not work." << std::endl;
	    return 0.0;
	}   
	
	template < typename _Potential >
	inline void
	FixedLocalTupleRgListInteractionTemplate < _Potential >::
	computeVirialX(std::vector<real> &p_xx_total, int bins) {
	    LOG4ESPP_INFO(theLogger, "compute virial p_xx of the pressure tensor slabwise");
	}
	
	
	template < typename _Potential > inline real
	FixedLocalTupleRgListInteractionTemplate < _Potential >::
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
	FixedLocalTupleRgListInteractionTemplate < _Potential >::computeVirialTensor(Tensor& w){
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
	FixedLocalTupleRgListInteractionTemplate < _Potential >::
	computeVirialTensor(Tensor& w, real z){
	    std::cout<<"Warning! Calculating virial layerwise is not supported for "
		"long range interactions."<<std::endl;
	}
	
	template < typename _Potential > inline void
	FixedLocalTupleRgListInteractionTemplate < _Potential >::
	computeVirialTensor(Tensor *w, int n){
	    std::cout<<"Warning! Calculating virial layerwise is not supported for "
		"long range interactions."<<std::endl;
	}
    }
}
#endif
