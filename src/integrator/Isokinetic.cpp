#include "python.hpp"
#include "Isokinetic.hpp"

#include "types.hpp"
#include "System.hpp"
#include "storage/Storage.hpp"
#include "iterator/CellListIterator.hpp"
#include "esutil/RNG.hpp"

namespace espresso {

  using namespace iterator;

  namespace integrator {

    LOG4ESPP_LOGGER(Isokinetic::theLogger, "Isokinetic");

    Isokinetic::Isokinetic(shared_ptr<System> system) : SystemAccess(system)
    {
      temperature = 0.0;
      coupling    = 1; // couple to thermostat in every md step
      couplecount = 0;

      // if (!system->rng) {
      //  throw std::runtime_error("system has no RNG");
      // }

      // rng = system->rng;

      LOG4ESPP_INFO(theLogger, "Isokinetic constructed");
    }

    void Isokinetic::setTemperature(real _temperature)
    {
      temperature = _temperature;
    }

    real Isokinetic::getTemperature()
    {
      return temperature;
    }

    void Isokinetic::setCoupling(int _coupling)
    {
      coupling = _coupling;
    }

    int Isokinetic::getCoupling()
    {
      return coupling;
    }

    Isokinetic::~Isokinetic()
    {
    }
    
    void Isokinetic::rescaleVelocities() {
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
      real SkalingFaktor;

      System& system = getSystemRef();
      CellList realCells = system.storage->getRealCells();

      NPart_local = system.storage->getNRealParticles();
      printf("CPU %i: NPart=%i NPart_local=%i\n", getSystem()->comm->rank(), NPart, NPart_local);

      for (CellListIterator cit(realCells); !cit.isDone(); ++cit) {
        Real3D vel = cit->velocity();
        EKin_local += 0.5 * cit->mass() * (vel * vel);
      }

      // boost::mpi::reduce(*getSystem()->comm, EKin_local, EKin, std::plus<real>(), 0);
      boost::mpi::all_reduce(*getSystem()->comm, NPart_local, NPart, std::plus<int>());
      //printf("CPU %i: Ekin=%f   NPart=%i NPart_local=%i\n", getSystem()->comm->rank(), EKin, NPart, NPart_local);

      DegreesOfFreedom   = 3.0/2.0 * NPart;
      currentTemperature = EKin / DegreesOfFreedom;
      SkalingFaktor      = sqrt(temperature / currentTemperature);

      for (CellListIterator cit(realCells); !cit.isDone(); ++cit) {
        cit->velocity() *= SkalingFaktor;
      }

    }

    /****************************************************
    ** REGISTRATION WITH PYTHON
    ****************************************************/

    void Isokinetic::registerPython() {

      using namespace espresso::python;

      class_<Isokinetic, shared_ptr<Isokinetic> >

        ("integrator_Isokinetic", init< shared_ptr<System> >())

        .add_property("temperature", &Isokinetic::getTemperature, &Isokinetic::setTemperature)
        .add_property("coupling", &Isokinetic::getCoupling, &Isokinetic::setCoupling)
        ;
    }

  }
}
