#include "python.hpp"
#include "StochasticVelocityRescaling.hpp"

#include "types.hpp"
#include "System.hpp"
#include "storage/Storage.hpp"
#include "iterator/CellListIterator.hpp"
#include "esutil/RNG.hpp"

namespace espresso {

  using namespace iterator;

  namespace integrator {

    LOG4ESPP_LOGGER(StochasticVelocityRescaling::theLogger, "StochasticVelocityRescaling");

    StochasticVelocityRescaling::StochasticVelocityRescaling(shared_ptr<System> system) : SystemAccess(system)
    {
      temperature = 0.0;
      coupling    = 1; // couple to thermostat in every md step
      couplecount = 0;

      // if (!system->rng) {
      //  throw std::runtime_error("system has no RNG");
      // }

      // rng = system->rng;

      LOG4ESPP_INFO(theLogger, "StochasticVelocityRescaling constructed");
    }

    void StochasticVelocityRescaling::setTemperature(real _temperature)
    {
      temperature = _temperature;
    }

    real StochasticVelocityRescaling::getTemperature()
    {
      return temperature;
    }

    void StochasticVelocityRescaling::setCoupling(int _coupling)
    {
      coupling = _coupling;
    }

    int StochasticVelocityRescaling::getCoupling()
    {
      return coupling;
    }

    StochasticVelocityRescaling::~StochasticVelocityRescaling()
    {
    }
    
    void StochasticVelocityRescaling::rescaleVelocities() {
      LOG4ESPP_DEBUG(theLogger, "rescaleVelocities");

      couplecount++;
      if (couplecount < coupling) {
    	  return;
      } else {
    	  couplecount = 0;
      }

      int NPart_local, NPart;
      real EKin = 0.0;
      real EKin_local = 0.0;
      real DegreesOfFreedom;
      real currentTemperature;
      real ScalingFactor;

      System& system = getSystemRef();
      CellList realCells = system.storage->getRealCells();

      for (CellListIterator cit(realCells); !cit.isDone(); ++cit) {
        Real3D vel = cit->velocity();
        EKin_local += 0.5 * cit->mass() * (vel * vel);
      }

      NPart_local = system.storage->getNRealParticles();

      boost::mpi::all_reduce(*getSystem()->comm, EKin_local, EKin, std::plus<real>());
      boost::mpi::all_reduce(*getSystem()->comm, NPart_local, NPart, std::plus<int>());

      DegreesOfFreedom   = 3.0/2.0 * NPart;
      currentTemperature = EKin / DegreesOfFreedom;
      ScalingFactor      = sqrt(temperature / currentTemperature);

      for (CellListIterator cit(realCells); !cit.isDone(); ++cit) {
        cit->velocity() *= ScalingFactor;
      }

    }

    /****************************************************
    ** REGISTRATION WITH PYTHON
    ****************************************************/

    void StochasticVelocityRescaling::registerPython() {

      using namespace espresso::python;

      class_<StochasticVelocityRescaling, shared_ptr<StochasticVelocityRescaling> >

        ("integrator_StochasticVelocityRescaling", init< shared_ptr<System> >())

        .add_property("temperature", &StochasticVelocityRescaling::getTemperature, &StochasticVelocityRescaling::setTemperature)
        .add_property("coupling", &StochasticVelocityRescaling::getCoupling, &StochasticVelocityRescaling::setCoupling)
        ;
    }

  }
}
