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
#include "FixPositions.hpp"

#include "types.hpp"
#include "System.hpp"
#include "storage/Storage.hpp"
#include "iterator/CellListIterator.hpp"

namespace espressopp {

  using namespace iterator;

  namespace integrator {

    LOG4ESPP_LOGGER(FixPositions::theLogger, "FixPositions");

    FixPositions::FixPositions(shared_ptr<System> system, shared_ptr< ParticleGroup > _particleGroup, const Int3D& _fixMask)
          : Extension(system), particleGroup(_particleGroup), fixMask(_fixMask)
    {
      LOG4ESPP_INFO(theLogger, "Isokinetic constructed");
    }

    void FixPositions::disconnect(){
      _befIntP.disconnect();
      _aftIntP.disconnect();
    }

    void FixPositions::connect(){
      // connection to initialisation
      _befIntP  = integrator->befIntP.connect( boost::bind(&FixPositions::savePositions, this));
      _aftIntP  = integrator->aftIntP.connect( boost::bind(&FixPositions::restorePositions, this));
    }

    void FixPositions::setParticleGroup(shared_ptr< ParticleGroup > _particleGroup) {
     	particleGroup = _particleGroup;
     }

     shared_ptr< ParticleGroup > FixPositions::getParticleGroup() {
     	return particleGroup;
     }

     void FixPositions::setFixMask(Int3D& _fixMask) {
     	fixMask = _fixMask;
     }

     Int3D& FixPositions::getFixMask() {
     	return fixMask;
     }

     void FixPositions::savePositions() {
    	 savePos.clear();
    	 for (ParticleGroup::iterator it=particleGroup->begin(); it != particleGroup->end(); it++ ) {
             savePos.push_back(std::pair<Particle *, Real3D>(*it, it->getPos()));
    	 }
     }

     void FixPositions::restorePositions() {
    	 for (std::list< std::pair<Particle *, Real3D> >::iterator it=savePos.begin(); it!=savePos.end(); it++) {
    		 Real3D savpos = it->second;
    		 Real3D newpos = it->first->getPos();
                 Real3D velo   = it->first->getV();
    		 for (int i=0; i<3; i++) {
    			 if (fixMask[i] != 0) {
    				 newpos[i] = savpos[i];
                                 velo[i] = 0.0;
    			 }
    		 }
   	         it->first->setPos(newpos);
   	         it->first->setV( velo );
   	         //it->first->setF( Real3D(0,0,0) );
    	 }
     }

    /****************************************************
    ** REGISTRATION WITH PYTHON
    ****************************************************/

    void FixPositions::registerPython() {

      using namespace espressopp::python;

      class_<FixPositions, shared_ptr<FixPositions>, bases<Extension> >

        ("integrator_FixPositions", init< shared_ptr< System >, shared_ptr< ParticleGroup >, const Int3D& >())
        .add_property("particleGroup", &FixPositions::getParticleGroup, &FixPositions::setParticleGroup)
        .def("getFixMask", &FixPositions::getFixMask, return_value_policy< reference_existing_object >() )
        .def("setFixMask", &FixPositions::setFixMask )
        .def("connect", &FixPositions::connect)
        .def("disconnect", &FixPositions::disconnect)
        ;
    }

  }
}
