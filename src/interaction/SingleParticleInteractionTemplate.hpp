/*
  Copyright (C) 2014
      Pierre de Buyl
  Copyright (C) 2012,2013,2015,2016,2017,2018
      Max Planck Institute for Polymer Research
  Copyright (C) 2008,2009,2010,2011
      Max-Planck-Institute for Polymer Research & Fraunhofer SCAI

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
#ifndef _INTERACTION_SINGLEPARTICLEINTERACTIONTEMPLATE_HPP
#define _INTERACTION_SINGLEPARTICLEINTERACTIONTEMPLATE_HPP

#include "mpi.hpp"
#include "Interaction.hpp"
#include "Real3D.hpp"
#include "Tensor.hpp"
#include "Particle.hpp"
#include "bc/BC.hpp"
#include "SystemAccess.hpp"
#include "storage/Storage.hpp"
#include "Interaction.hpp"
#include "types.hpp"
#include "SingleParticlePotential.hpp"
#include "iterator/CellListIterator.hpp"

namespace espressopp {
  namespace interaction {

    /** This class provides a template for a single-particle interaction,
        typically used for external forces on the system.

        The force is applied to all the particles in the system.
     */
    template < typename _Potential >
    class SingleParticleInteractionTemplate: public Interaction, SystemAccess {

    protected:
      typedef _Potential Potential;

    public:
      SingleParticleInteractionTemplate
      (shared_ptr < System > system,
       shared_ptr < _Potential > _potential)
        : SystemAccess(system),
          potential(_potential)
      {
        if (! potential) {
          LOG4ESPP_ERROR(theLogger, "NULL potential");
        }
      }

      virtual ~SingleParticleInteractionTemplate() {};

      void
      setPotential(shared_ptr < _Potential> _potential) {
        if (_potential) {
          potential = _potential;
        } else {
          LOG4ESPP_ERROR(theLogger, "NULL potential");
        }
      }

      shared_ptr < _Potential > getPotential() {
        return potential;
      }

      virtual void addForces();
      virtual real computeEnergy();
      virtual real computeEnergyDeriv();
      virtual real computeEnergyAA();
      virtual real computeEnergyCG();
      virtual real computeEnergyAA(int atomtype);
      virtual real computeEnergyCG(int atomtype);
      virtual void computeVirialX(std::vector<real> &p_xx_total, int bins);
      virtual real computeVirial();
      virtual void computeVirialTensor(Tensor& w);
      virtual void computeVirialTensor(Tensor& w, real z);
      virtual void computeVirialTensor(Tensor *w, int n);
      virtual real getMaxCutoff();
      virtual int bondType() { return Single; }

    protected:
      int ntypes;
      shared_ptr < Potential > potential;
    };

    //////////////////////////////////////////////////
    // INLINE IMPLEMENTATION
    //////////////////////////////////////////////////
    template < typename _Potential > inline void
    SingleParticleInteractionTemplate < _Potential >::addForces() {
      LOG4ESPP_INFO(theLogger, "add forces computed for all particles");

      System& system = getSystemRef();
      CellList realCells = system.storage->getRealCells();
      const bc::BC& bc = *system.bc;

      for(iterator::CellListIterator cit(realCells); !cit.isDone(); ++cit) {
        Real3D force;
        if(potential->_computeForce(force, *cit, bc)) {
          cit->force() += force;
        }
      }
    }

    template < typename _Potential > inline real
    SingleParticleInteractionTemplate < _Potential >::
    computeEnergy() {
      LOG4ESPP_INFO(theLogger, "compute energy for all particles");

      System& system = getSystemRef();
      CellList realCells = system.storage->getRealCells();
      const bc::BC& bc = *system.bc;

      real e=0.0;
      for(iterator::CellListIterator cit(realCells); !cit.isDone(); ++cit) {
        e += potential->_computeEnergy(*cit, bc);
      }

      real esum;
      boost::mpi::all_reduce(*mpiWorld, e, esum, std::plus<real>());
      return esum;
    }

    template < typename _Potential > inline real
    SingleParticleInteractionTemplate < _Potential >::
    computeEnergyDeriv() {
      std::cout << "Warning! At the moment computeEnergyDeriv() in SingleParticleInteractionTemplate does not work." << std::endl;
      return 0.0;
    }

    template < typename _Potential > inline real
    SingleParticleInteractionTemplate < _Potential >::
    computeEnergyAA() {
      LOG4ESPP_INFO(theLogger, "compute force of the VSphere Self potential");
      LOG4ESPP_INFO(theLogger, "Warning! At the moment computeEnergyAA() in SingleParticleInteractionTemplate does not work.");
      return 0.0;
    }

    template < typename _Potential > inline real
    SingleParticleInteractionTemplate < _Potential >::
    computeEnergyAA(int atomtype) {
      LOG4ESPP_INFO(theLogger, "compute force of the VSphere Self potential");
      LOG4ESPP_INFO(theLogger, "Warning! At the moment computeEnergyAA(int atomtype) in SingleParticleInteractionTemplate does not work.");
      return 0.0;
    }

    template < typename _Potential > inline real
    SingleParticleInteractionTemplate < _Potential >::
    computeEnergyCG() {
      LOG4ESPP_INFO(theLogger, "Warning! At the moment computeEnergyCG() in SingleParticleInteractionTemplate does not work.");
      return 0.0;
    }

    template < typename _Potential > inline real
    SingleParticleInteractionTemplate < _Potential >::
    computeEnergyCG(int atomtype) {
      LOG4ESPP_INFO(theLogger, "Warning! At the moment computeEnergyCG(int atomtype) in SingleParticleInteractionTemplate does not work.");
      return 0.0;
    }

    template < typename _Potential >
    inline void
    SingleParticleInteractionTemplate < _Potential >::
    computeVirialX(std::vector<real> &p_xx_total, int bins) {
              LOG4ESPP_INFO(theLogger, "compute virial p_xx of the pressure tensor slabwise");
              LOG4ESPP_INFO(theLogger, "Warning! At the moment computeVirialX() in SingleParticleInteractionTemplate does not work.");
    }

    template < typename _Potential > inline real
    SingleParticleInteractionTemplate < _Potential >::
    computeVirial() {
      LOG4ESPP_INFO(theLogger, "compute the virial for all particles");
      LOG4ESPP_INFO(theLogger, "Warning! At the moment computeVirialX() in SingleParticleInteractionTemplate does not work.");
      return 0.0;
    }

    template < typename _Potential > inline void
    SingleParticleInteractionTemplate < _Potential >::computeVirialTensor(Tensor& w){
      LOG4ESPP_INFO(theLogger, "compute the virial tensor for all particles");
      LOG4ESPP_INFO(theLogger, "Warning! At the moment computeVirialTensor() in SingleParticleInteractionTemplate does not work.");
    }

    template < typename _Potential > inline void
    SingleParticleInteractionTemplate < _Potential >::
    computeVirialTensor(Tensor& w, real z){
      LOG4ESPP_INFO(theLogger, "compute the virial tensor for all particles");
      LOG4ESPP_INFO(theLogger, "Warning! At the moment computeVirialTensor() in SingleParticleInteractionTemplate does not work.");
    }

    template < typename _Potential > inline void
    SingleParticleInteractionTemplate < _Potential >::
    computeVirialTensor(Tensor *w, int n){
      LOG4ESPP_INFO(theLogger, "compute the virial tensor for all particles");
      LOG4ESPP_INFO(theLogger, "Warning! At the moment computeVirialTensor() in SingleParticleInteractionTemplate does not work.");
    }

    template < typename _Potential >
    inline real
    SingleParticleInteractionTemplate< _Potential >::
    getMaxCutoff() {
      return potential->getMaxCutoff();
    }

  }
}
#endif
