/*
  Copyright (C) 2012,2013
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

#include "python.hpp"
#include "CapForce.hpp"

#include "types.hpp"
#include "System.hpp"
#include "storage/Storage.hpp"
#include "iterator/CellListIterator.hpp"

namespace espressopp {

  using namespace iterator;

  namespace integrator {

    LOG4ESPP_LOGGER(CapForce::theLogger, "CapForce");

    CapForce::CapForce(shared_ptr<System> system, const Real3D& _capForce)
    : Extension(system), capForce(_capForce)
    {
      LOG4ESPP_INFO(theLogger, "Force capping for all particles constructed");
      allParticles = true;
      absCapping   = false;
      adress = false;
    }

    CapForce::CapForce(shared_ptr<System> system, real _absCapForce)
    : Extension(system), absCapForce(_absCapForce)
    {
      LOG4ESPP_INFO(theLogger, "Force capping for all particles constructed");
      allParticles = true;
      absCapping   = true;
      adress = false;
    }

    CapForce::CapForce(shared_ptr<System> system, const Real3D& _capForce, shared_ptr< ParticleGroup > _particleGroup)
    : Extension(system), capForce(_capForce), particleGroup(_particleGroup)
    {
      LOG4ESPP_INFO(theLogger, "Force capping for particle group constructed");
      allParticles = false;
      absCapping   = false;
      adress = false;
    }

    CapForce::CapForce(shared_ptr<System> system, real _absCapForce, shared_ptr< ParticleGroup > _particleGroup)
    : Extension(system), absCapForce(_absCapForce), particleGroup(_particleGroup)
    {
      LOG4ESPP_INFO(theLogger, "Force capping for particle group constructed");
      allParticles = false;
      absCapping   = true;
      adress = false;
    }

    void CapForce::disconnect(){
      _aftCalcF.disconnect();
    }

    void CapForce::connect(){
      // connection to initialisation
      if (!allParticles) {
        _aftCalcF  = integrator->aftCalcF.connect( boost::bind(&CapForce::applyForceCappingToGroup, this), boost::signals2::at_back);
      } else {
    	_aftCalcF  = integrator->aftCalcF.connect( boost::bind(&CapForce::applyForceCappingToAll, this), boost::signals2::at_back);
      }
    }

    void CapForce::setCapForce(Real3D& _capForce) {
    	capForce = _capForce;
    }

    void CapForce::setAbsCapForce(real _absCapForce) {
    	absCapForce = _absCapForce;
    }

    Real3D& CapForce::getCapForce() {
    	return capForce;
    }

    real CapForce::getAbsCapForce() {
    	return absCapForce;
    }
    
    void CapForce::setAdress(bool _adress){
        adress = _adress;
    }

    bool CapForce::getAdress(){
        return adress;
    }

    void CapForce::setParticleGroup(shared_ptr< ParticleGroup > _particleGroup) {
        particleGroup = _particleGroup;
        if (allParticles) {
     	  disconnect();
     	  allParticles = false;
     	  connect();
        }
    }

     shared_ptr< ParticleGroup > CapForce::getParticleGroup() {
    	//if (!allParticles) {
     	  return particleGroup;
    	//}
     }

     void CapForce::applyForceCappingToGroup() {
       LOG4ESPP_DEBUG(theLogger, "applying force capping to particle group of size " << particleGroup->size());

       if (absCapping) {
    	 real capfsq = absCapForce * absCapForce;
         for (ParticleGroup::iterator it=particleGroup->begin(); it != particleGroup->end(); it++ ) {
       	   LOG4ESPP_DEBUG(theLogger, "applying scalar force capping to particle " << it->getId());
           real fsq = it->force().sqr();
           if (fsq > capfsq) {
             real scaling = sqrt(capfsq / fsq);
             Real3D& f=it->force();
      	     for (int dir=0; dir<3; dir++) {
      	    	f[dir] *= scaling;
      	     }
           }
         }
       } else {
         for (ParticleGroup::iterator it=particleGroup->begin(); it != particleGroup->end(); it++ ) {
           LOG4ESPP_DEBUG(theLogger, "applying vector force capping to particle " << it->getId());
           Real3D& f=it->force();
      	   for (int dir=0; dir<3; dir++) {
             if (f[dir]>0 && f[dir]>capForce[dir]) {
               f[dir] = capForce[dir];
             }
             else if (f[dir]<0 && f[dir] < -capForce[dir]) {
               f[dir] = -capForce[dir];
             }
      	   }
         }
       }
     }

     void CapForce::applyForceCappingToAll() {
       LOG4ESPP_DEBUG(theLogger, "applying force capping to all particles");

       System& system = getSystemRef();
       CellList realCells = system.storage->getRealCells();
       if (absCapping) {
    	 real capfsq = absCapForce * absCapForce;
         for(CellListIterator cit(realCells); !cit.isDone(); ++cit) {
           real fsq = cit->force().sqr();
           if (fsq > capfsq) {
             real scaling = sqrt(capfsq / fsq);
             Real3D& f=cit->force();
      	     for (int dir=0; dir<3; dir++) f[dir] *= scaling;
             //std::cout << "Force Capping applied on particle " << cit->getId() << "\n";
           }
         }
       } else {
           for(CellListIterator cit(realCells); !cit.isDone(); ++cit) {
        	 Real3D& f=cit->force();
      	     for (int dir=0; dir<3; dir++) {
               if (f[dir]>0 && f[dir]>capForce[dir]) {
              	 f[dir] = capForce[dir];
               }
               else if (f[dir]<0 && f[dir] < -capForce[dir]) {
            	 f[dir] = -capForce[dir];
               }
      	     }
           }
       }
       
       if (adress && absCapping){
         real capfsq2 = absCapForce * absCapForce;
         ParticleList& adrATparticles = system.storage->getAdrATParticles();
         for (std::vector<Particle>::iterator it = adrATparticles.begin();
                     it != adrATparticles.end(); it++) {             
             real fsq2 = it->force().sqr();
             if (fsq2 > capfsq2) {
                real scaling2 = sqrt(capfsq2 / fsq2);
                Real3D& f2=it->force();
                for (int dir=0; dir<3; dir++) f2[dir] *= scaling2;
                //std::cout << "AdResS Force Capping applied on particle " << it->getId() << "\n";
             }                     
         }
       }       
     }

    /****************************************************
    ** REGISTRATION WITH PYTHON
    ****************************************************/

    void CapForce::registerPython() {

      using namespace espressopp::python;

      class_<CapForce, shared_ptr<CapForce>, bases<Extension> >

        ("integrator_CapForce", init< shared_ptr< System >, const Real3D& >())
        .def(init< shared_ptr< System >, real >())
        .def(init< shared_ptr< System >, const Real3D&, shared_ptr< ParticleGroup > >())
        .def(init< shared_ptr< System >, real, shared_ptr< ParticleGroup > >())
        .add_property("particleGroup", &CapForce::getParticleGroup, &CapForce::setParticleGroup)
        .add_property("adress", &CapForce::getAdress, &CapForce::setAdress)
        .def("getCapForce", &CapForce::getCapForce, return_value_policy< reference_existing_object >())
        .def("getAbsCapForce", &CapForce::getAbsCapForce)
        .def("setCapForce", &CapForce::setCapForce )
        .def("setAbsCapForce", &CapForce::setAbsCapForce )
        .def("connect", &CapForce::connect)
        .def("disconnect", &CapForce::disconnect)
        ;
    }

  }
}
