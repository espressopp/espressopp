#include "VelocityVerletFixedParticles.hpp"
#include "iterator/CellListIterator.hpp"
#include "interaction/Interaction.hpp"
#include "interaction/Potential.hpp"
#include "System.hpp"
#include "storage/Storage.hpp"
#include "mpi.hpp"

#ifdef VTRACE
#include "vampirtrace/vt_user.h"
#else
#define VT_TRACER( name)
#endif

namespace espresso {

  using namespace std;
  
  namespace integrator {

    using namespace interaction;
    using namespace iterator;
    using namespace esutil;

    VelocityVerletFixedParticles::VelocityVerletFixedParticles(shared_ptr< System > system, shared_ptr< ParticleGroup > _fixedParticles, Int3D _fixMask) : VelocityVerlet(system) {
      LOG4ESPP_INFO(theLogger, "construct VelocityVerletFixedParticles");
      fixedParticles = _fixedParticles;
      fixMask        = _fixMask;
      resortFlag     = true;
      maxDist        = 0.0;
    }

    VelocityVerletFixedParticles::~VelocityVerletFixedParticles() {
      LOG4ESPP_INFO(theLogger, "free VelocityVerlet");
    }

    real VelocityVerletFixedParticles::integrate1() {
      System& system     = getSystemRef();
      CellList realCells = system.storage->getRealCells();

      int count = 0;
      real maxSqDist = 0.0; // maximal square distance a particle moves
      // loop over all particles of the local cells
      for(CellListIterator cit(realCells); !cit.isDone(); ++cit) {
        LOG4ESPP_DEBUG(theLogger, "Particle " << cit->id() << ", pos = " << cit->position() << ", v = " << cit->velocity() << ", f = " << cit->force());
        real sqDist = 0.0;
        real dtfm   = 0.5 * dt / cit->mass();
        // Propagate velocities: v(t+0.5*dt) = v(t) + 0.5*dt * f(t) 
        cit->velocity() += dtfm * cit->force();
        // Propagate positions (only NVT): p(t + dt) = p(t) + dt * v(t+0.5*dt) 
        Real3D deltaP = cit->velocity();

        deltaP          *= dt;

        // if particle is in fixedParticles list
        // do not change position and set velocity to 0 (according to fixMask)
        if (fixedParticles->has(cit->id())) {
      	  for (int i=0; i<3; i++) {
      	    deltaP[i]          *= fixMask[i];
      	    cit->velocity()[i] *= fixMask[i];
          }
        }

        cit->position() += deltaP;
        sqDist          += deltaP * deltaP;
        maxSqDist        = std::max(maxSqDist, sqDist);
        count++;
      }
      
      // signal
      inIntP(maxSqDist);

      real maxAllSqDist;
      mpi::all_reduce(*system.comm, maxSqDist, maxAllSqDist, boost::mpi::maximum<real>());
      LOG4ESPP_INFO(theLogger, "moved " << count << " particles in integrate1" << ", max move local = " << sqrt(maxSqDist) << ", global = " << sqrt(maxAllSqDist));
      
      return sqrt(maxAllSqDist);
    }

    void VelocityVerletFixedParticles::setFixedParticles(shared_ptr< ParticleGroup > _fixedParticles) {
      fixedParticles = _fixedParticles;
    }

    shared_ptr< ParticleGroup > VelocityVerletFixedParticles::getFixedParticles() {
      return fixedParticles;
    }

    void VelocityVerletFixedParticles::setFixMask(Int3D& _fixMask) {
      fixMask = _fixMask;
    }

    Int3D& VelocityVerletFixedParticles::getFixMask() {
      return fixMask;
    }

    /****************************************************
    ** REGISTRATION WITH PYTHON
    ****************************************************/
    void VelocityVerletFixedParticles::registerPython() {

      using namespace espresso::python;

      // Note: use noncopyable and no_init for abstract classes
      class_<VelocityVerletFixedParticles, bases<VelocityVerlet>, boost::noncopyable >
        ("integrator_VelocityVerletFixedParticles", init< shared_ptr<System>, shared_ptr< ParticleGroup >, const Int3D& >())
        .add_property("fixedParticles", &VelocityVerletFixedParticles::getFixedParticles, &VelocityVerletFixedParticles::setFixedParticles)
        .def("getFixMask", &VelocityVerletFixedParticles::getFixMask, return_value_policy< reference_existing_object >() )
        .def("setFixMask", &VelocityVerletFixedParticles::setFixMask )
        ;
    }
  }
}
