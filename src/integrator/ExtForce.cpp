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
#include "ExtForce.hpp"

#include "types.hpp"
#include "System.hpp"
#include "storage/Storage.hpp"
#include "iterator/CellListIterator.hpp"

namespace espressopp {

  using namespace iterator;

  namespace integrator {

    LOG4ESPP_LOGGER(ExtForce::theLogger, "ExtForce");

    ExtForce::ExtForce(shared_ptr<System> system, const Real3D& _extForce)
    : Extension(system), extForce(_extForce)
    {
      LOG4ESPP_INFO(theLogger, "External Force for all particles constructed");
      allParticles = true;
    }

    ExtForce::ExtForce(shared_ptr<System> system, const Real3D& _extForce, shared_ptr< ParticleGroup > _particleGroup)
    : Extension(system), extForce(_extForce), particleGroup(_particleGroup)
    {
      LOG4ESPP_INFO(theLogger, "External Force for particle group constructed");
      allParticles = false;
    }

    void ExtForce::disconnect(){
      _aftInitF.disconnect();
    }

    void ExtForce::connect(){
      // connection to initialisation
      if (!allParticles) {
        _aftInitF  = integrator->aftInitF.connect( boost::bind(&ExtForce::applyForceToGroup, this));
      } else {
    	_aftInitF  = integrator->aftInitF.connect( boost::bind(&ExtForce::applyForceToAll, this));
      }
    }

    void ExtForce::setExtForce(Real3D& _extForce) {
    	extForce = _extForce;
    }

    Real3D& ExtForce::getExtForce() {
    	return extForce;
    }

    void ExtForce::setParticleGroup(shared_ptr< ParticleGroup > _particleGroup) {
        particleGroup = _particleGroup;
        if (allParticles) {
     	  disconnect();
     	  allParticles = false;
     	  connect();
        }
    }

     shared_ptr< ParticleGroup > ExtForce::getParticleGroup() {
    	//if (!allParticles) {
     	  return particleGroup;
    	//}
     }

     void ExtForce::applyForceToGroup() {
       LOG4ESPP_DEBUG(theLogger, "applying external force to particle group of size " << particleGroup->size());
       for (ParticleGroup::iterator it=particleGroup->begin(); it != particleGroup->end(); it++ ) {
    	 LOG4ESPP_DEBUG(theLogger, "applying external force to particle " << it->getId());
         it->force() += extForce;
       }
     }

     void ExtForce::applyForceToAll() {
       System& system = getSystemRef();
       CellList realCells = system.storage->getRealCells();
       for(CellListIterator cit(realCells); !cit.isDone(); ++cit) {
           cit->force() += extForce;
       }
     }

    /****************************************************
    ** REGISTRATION WITH PYTHON
    ****************************************************/

    void ExtForce::registerPython() {

      using namespace espressopp::python;

      class_<ExtForce, shared_ptr<ExtForce>, bases<Extension> >

        ("integrator_ExtForce", init< shared_ptr< System >, const Real3D& >())
        .def(init< shared_ptr< System >, const Real3D&, shared_ptr< ParticleGroup > >())
        .add_property("particleGroup", &ExtForce::getParticleGroup, &ExtForce::setParticleGroup)
        .def("getExtForce", &ExtForce::getExtForce, return_value_policy< reference_existing_object >())
        .def("setExtForce", &ExtForce::setExtForce )
        .def("connect", &ExtForce::connect)
        .def("disconnect", &ExtForce::disconnect)
        ;
    }

  }
}
