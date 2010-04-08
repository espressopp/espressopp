#include "python.hpp"
#include "Langevin.hpp"

#include "types.hpp"
#include "System.hpp"
#include "storage/Storage.hpp"
#include "iterator/CellListIterator.hpp"
#include "esutil/RNG.hpp"

namespace espresso {

  using namespace iterator;

  namespace integrator {

    LOG4ESPP_LOGGER(Langevin::theLogger, "Langevin");

    Langevin::Langevin(shared_ptr<System> system) : SystemAccess(system)
    {
      gamma  = 0.0;
      temperature = 0.0;

      if (!system->rng) {
        throw std::runtime_error("system has no RNG");
      }

      rng = system->rng;

      LOG4ESPP_INFO(theLogger, "Langevin constructed");
    }

    void Langevin::setGamma(real _gamma)
    {
      gamma = _gamma;
    }

    real Langevin::getGamma()
    {
      return gamma;
    }

    void Langevin::setTemperature(real _temperature)
    {
      temperature = _temperature;
    }

    real Langevin::getTemperature()
    {
      return temperature;
    }

    Langevin::~Langevin()
    {
    }

    void Langevin::thermalize()
    {
      LOG4ESPP_DEBUG(theLogger, "thermalize");

      System& system = getSystemRef();

      CellList cells = system.storage->getRealCells();

      for(CellListIterator cit(cells); !cit.isDone(); ++cit) {
        frictionThermo(*cit);
      }
    }

    void Langevin::frictionThermo(Particle& p)
    {
      // double massf = sqrt(PMASS(*p));
 
      double massf = 1.0;

      for (int j = 0 ; j < 3 ; j++) {
	p.f.f[j] += pref1 * p.m.v[j] + pref2*((*rng)()-0.5)*massf;
      }

      LOG4ESPP_TRACE(theLogger, "new force of p = " << p.f.f[0] << " "
		                 << p.f.f[1] << " " << p.f.f[2]);
    }

    void Langevin::initialize(real timestep)

    { // calculate the prefactors

      LOG4ESPP_INFO(theLogger, "init, timestep = " << timestep <<
		    ", gamma = " << gamma << 
		    ", temperature = " << temperature);

      pref1 = -gamma * timestep;
      pref2 = sqrt(24.0 * temperature * gamma / timestep);

    }

    /** very nasty: if we recalculate force when leaving/reentering the integrator,
	a(t) and a((t-dt)+dt) are NOT equal in the vv algorithm. The random
	numbers are drawn twice, resulting in a different variance of the random force.
	This is corrected by additional heat when restarting the integrator here.
	Currently only works for the Langevin thermostat, although probably also others
	are affected.
    */

    void Langevin::heatUp()
    {
      LOG4ESPP_INFO(theLogger, "heatUp");

      pref2buffer = pref2;
      pref2       *= sqrt(3.0);
    }

    /** Opposite to heatUp */

    void Langevin::coolDown()
    {
      LOG4ESPP_INFO(theLogger, "coolDown");

      pref2 = pref2buffer;
    }

    /****************************************************
    ** REGISTRATION WITH PYTHON
    ****************************************************/

    void Langevin::registerPython() {

      using namespace espresso::python;

      class_<Langevin, shared_ptr<Langevin> >

        ("integrator_Langevin", init< shared_ptr<System> >())

        .add_property("gamma", &Langevin::getGamma, &Langevin::setGamma)
        .add_property("temperature", &Langevin::getTemperature, &Langevin::setTemperature)
        ;
    }

  }
}
