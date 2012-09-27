#include "python.hpp"
#include "FixPositions.hpp"

#include "types.hpp"
#include "System.hpp"
#include "storage/Storage.hpp"
#include "iterator/CellListIterator.hpp"

namespace espresso {

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
    		 for (int i=0; i<3; i++) {
    			 if (fixMask[i] != 0) {
    				 newpos[i] = savpos[i];
    			 }
    		 }
   	         it->first->setPos(newpos);
   	         it->first->setV( Real3D(0.0, 0.0, 0.0) );
   	         //it->first->setF( Real3D(0,0,0) );
    	 }
     }

    /****************************************************
    ** REGISTRATION WITH PYTHON
    ****************************************************/

    void FixPositions::registerPython() {

      using namespace espresso::python;

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
