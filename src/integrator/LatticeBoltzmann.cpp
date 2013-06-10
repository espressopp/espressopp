#include "python.hpp"
#include "LatticeBoltzmann.hpp"

#include "types.hpp"
#include "System.hpp"
#include "storage/Storage.hpp"
#include "iterator/CellListIterator.hpp"

namespace espresso {

  using namespace iterator;

  namespace integrator {

    LOG4ESPP_LOGGER(LatticeBoltzmann::theLogger, "LatticeBoltzmann");

    LatticeBoltzmann::LatticeBoltzmann(shared_ptr<System> system, shared_ptr< ParticleGroup > _particleGroup, const Int3D& _fixMask)
          : Extension(system), particleGroup(_particleGroup), fixMask(_fixMask)
    {
      LOG4ESPP_INFO(theLogger, "Isokinetic constructed");
    }

    void LatticeBoltzmann::disconnect(){
      _befIntP.disconnect();
      _aftIntP.disconnect();
    }

    void LatticeBoltzmann::connect(){
      // connection to initialisation
      _befIntP  = integrator->befIntP.connect( boost::bind(&LatticeBoltzmann::savePositions, this));
      _aftIntP  = integrator->aftIntP.connect( boost::bind(&LatticeBoltzmann::restorePositions, this));
    }

    void LatticeBoltzmann::setParticleGroup(shared_ptr< ParticleGroup > _particleGroup) {
     	particleGroup = _particleGroup;
     }

     shared_ptr< ParticleGroup > LatticeBoltzmann::getParticleGroup() {
     	return particleGroup;
     }

     void LatticeBoltzmann::setFixMask(Int3D& _fixMask) {
     	fixMask = _fixMask;
     }

     Int3D& LatticeBoltzmann::getFixMask() {
     	return fixMask;
     }

     void LatticeBoltzmann::savePositions() {
    	 savePos.clear();
    	 for (ParticleGroup::iterator it=particleGroup->begin(); it != particleGroup->end(); it++ ) {
             savePos.push_back(std::pair<Particle *, Real3D>(*it, it->getPos()));
    	 }
     }

     void LatticeBoltzmann::restorePositions() {
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

    void LatticeBoltzmann::registerPython() {

      using namespace espresso::python;

      class_<LatticeBoltzmann, shared_ptr<LatticeBoltzmann>, bases<Extension> >

        ("integrator_LatticeBoltzmann", init< shared_ptr< System >, shared_ptr< ParticleGroup >, const Int3D& >())
        .add_property("particleGroup", &LatticeBoltzmann::getParticleGroup, &LatticeBoltzmann::setParticleGroup)
        .def("getFixMask", &LatticeBoltzmann::getFixMask, return_value_policy< reference_existing_object >() )
        .def("setFixMask", &LatticeBoltzmann::setFixMask )
        .def("connect", &LatticeBoltzmann::connect)
        .def("disconnect", &LatticeBoltzmann::disconnect)
        ;
    }

  }
}
