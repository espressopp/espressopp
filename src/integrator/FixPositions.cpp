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
          : SystemAccess(system), particleGroup(_particleGroup), fixMask(_fixMask)
    {
      LOG4ESPP_INFO(theLogger, "Isokinetic constructed");
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

    void FixPositions::apply(longint pid, Real3D& vel, Real3D& dp) {
    	if (particleGroup->has(pid)) {
    	    for (int i=0; i<3; i++) {
    		     dp[i] *= fixMask[i];
    		    vel[i] *= fixMask[i];
    	    }
    	}
    }

    /****************************************************
    ** REGISTRATION WITH PYTHON
    ****************************************************/

    void FixPositions::registerPython() {

      using namespace espresso::python;

      class_<FixPositions, shared_ptr<FixPositions> >

        ("integrator_FixPositions", init< shared_ptr< System >, shared_ptr< ParticleGroup >, const Int3D& >())
        .add_property("particleGroup", &FixPositions::getParticleGroup, &FixPositions::setParticleGroup)
        .def("getFixMask", &FixPositions::getFixMask, return_value_policy< reference_existing_object >() )
        .def("setFixMask", &FixPositions::setFixMask )
        ;
    }

  }
}
