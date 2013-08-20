#include "python.hpp"
#include "VelocityVerletOnRadius.hpp"

#include "types.hpp"
#include "System.hpp"
#include "storage/Storage.hpp"
#include "iterator/CellListIterator.hpp"

namespace espresso {

  using namespace iterator;

  namespace integrator {

    LOG4ESPP_LOGGER(VelocityVerletOnRadius::theLogger, "VelocityVerletOnRadius");

    VelocityVerletOnRadius::VelocityVerletOnRadius(shared_ptr<System> system) : Extension(system)
    {
      LOG4ESPP_INFO(theLogger, "VelocityVerletOnRadius constructed");
    }

    void VelocityVerletOnRadius::disconnect(){
      _befIntP.disconnect();
      _aftIntP.disconnect();
    }

    void VelocityVerletOnRadius::connect(){
      // connection to initialisation
      //_befIntP  = integrator->befIntP.connect( boost::bind(&FixPositions::savePositions, this));
      //_aftIntP  = integrator->aftIntP.connect( boost::bind(&FixPositions::restorePositions, this));
    }


    // aftInitF(); // signal after initForces {cit->force = 0.0}

    /****************************************************
    ** REGISTRATION WITH PYTHON
    ****************************************************/

    void VelocityVerletOnRadius::registerPython() {

      using namespace espresso::python;

      class_<VelocityVerletOnRadius, shared_ptr<VelocityVerletOnRadius>, bases<Extension> >
        ("integrator_VelocityVerletOnRadius", init< shared_ptr< System > >())
        .def("connect", &VelocityVerletOnRadius::connect)
        .def("disconnect", &VelocityVerletOnRadius::disconnect)
        ;
    }

  }
}
