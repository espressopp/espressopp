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
      _aftIntP.disconnect();
      _aftIntV.disconnect();
      _aftInitF.disconnect();
    }

    void VelocityVerletOnRadius::connect(){
      // connection to initialisation
      _aftIntP  = integrator->aftIntP.connect( boost::bind(&VelocityVerletOnRadius::integrate1, this));
      _aftIntV  = integrator->aftIntV.connect( boost::bind(&VelocityVerletOnRadius::integrate2, this));
      _aftInitF  = integrator->aftInitF.connect( boost::bind(&VelocityVerletOnRadius::initForces, this));
    }

    void VelocityVerletOnRadius::integrate1() {
      System& system = getSystemRef();
      CellList realCells = system.storage->getRealCells();
      // loop over all particles of the local cells
      for(CellListIterator cit(realCells); !cit.isDone(); ++cit) {
        LOG4ESPP_DEBUG(theLogger, "Particle " << cit->id() <<  ", radius = " << cit->radius());
        real dt = integrator->getTimeStep();
        real dtfm = 0.5 * dt / cit->mass();
        // Propagate velocities: v(t+0.5*dt) = v(t) + 0.5*dt * f(t)
        cit->vradius() += dtfm * cit->fradius();
        // Propagate radius (only NVT): p(t + dt) = p(t) + dt * v(t+0.5*dt)
        real deltaP = cit->vradius();
        deltaP *= dt;
        cit->radius() += deltaP;
      }
    }

    void VelocityVerletOnRadius::integrate2() {
      System& system = getSystemRef();
      CellList realCells = system.storage->getRealCells();
      // loop over all particles of the local cells
      real dt = integrator->getTimeStep();
      real half_dt = 0.5 * dt;
      for(CellListIterator cit(realCells); !cit.isDone(); ++cit) {
        real dtfm = half_dt / cit->mass();
        /* Propagate radius velocities: v(t+0.5*dt) = v(t) + 0.5*dt * f(t) */
        cit->vradius() += dtfm * cit->fradius();
      }
    }

    void VelocityVerletOnRadius::initForces() {
      // fradius is initialized for real + ghost particles
      System& system = getSystemRef();
      CellList localCells = system.storage->getLocalCells();
      LOG4ESPP_INFO(theLogger, "init fradius for real + ghost particles");
      for(CellListIterator cit(localCells); !cit.isDone(); ++cit) {
        cit->fradius() = 0.0;
      }
    }


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
