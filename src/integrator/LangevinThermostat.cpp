#include "python.hpp"
#include "LangevinThermostat.hpp"

#include "types.hpp"
#include "System.hpp"
#include "storage/Storage.hpp"
#include "iterator/CellListIterator.hpp"
#include "esutil/RNG.hpp"

namespace espresso {

  namespace integrator {

    using namespace espresso::iterator;


    LangevinThermostat::LangevinThermostat(shared_ptr<System> system)
    :Extension(system) {


      gamma  = 0.0;
      temperature = 0.0;

      if (!system->rng) {
        throw std::runtime_error("system has no RNG");
      }

      rng = system->rng;

      LOG4ESPP_INFO(theLogger, "Langevin constructed");


    }

    void LangevinThermostat::setGamma(real _gamma)
    {
      gamma = _gamma;
    }

    real LangevinThermostat::getGamma()
    {
      return gamma;
    }

    void LangevinThermostat::setAdress(bool _adress){
        adress = _adress;
    }

    bool LangevinThermostat::getAdress(){
        return adress;
    }

    void LangevinThermostat::setTemperature(real _temperature)
    {
      temperature = _temperature;
    }

    real LangevinThermostat::getTemperature()
    {
      return temperature;
    }

    LangevinThermostat::~LangevinThermostat() {
        disconnect();
    }


    void LangevinThermostat::disconnect() {

        _initialize.disconnect();
        _heatUp.disconnect();
        _coolDown.disconnect();
        _thermalize.disconnect();
        _thermalizeAdr.disconnect();

    }

    void LangevinThermostat::connect() {

        // connect to initialization inside run()
        _initialize = integrator->runInit.connect(
                boost::bind(&LangevinThermostat::initialize, this));

        _heatUp = integrator->recalc1.connect(
                boost::bind(&LangevinThermostat::heatUp, this));

        _coolDown = integrator->recalc2.connect(
                boost::bind(&LangevinThermostat::coolDown, this));

        if (adress) {
            _thermalizeAdr = integrator->aftCalcF.connect(
                boost::bind(&LangevinThermostat::thermalizeAdr, this));
        }
        else {
            _thermalize = integrator->aftCalcF.connect(
                boost::bind(&LangevinThermostat::thermalize, this));
        }
    }


    void LangevinThermostat::thermalize()
    {
      LOG4ESPP_DEBUG(theLogger, "thermalize");

      System& system = getSystemRef();
      
      CellList cells = system.storage->getRealCells();

      for(CellListIterator cit(cells); !cit.isDone(); ++cit) {
        frictionThermo(*cit);
      }
    }

    // for AdResS
    void LangevinThermostat::thermalizeAdr()
    {
      LOG4ESPP_DEBUG(theLogger, "thermalize");

      System& system = getSystemRef();

      // thermalize CG particles
      CellList cells = system.storage->getRealCells();
      for(CellListIterator cit(cells); !cit.isDone(); ++cit) {
        frictionThermo(*cit);
      }

      // thermalize AT particles
      ParticleList& adrATparticles = system.storage->getAdrATParticles();
      for (std::vector<Particle>::iterator it = adrATparticles.begin();
              it != adrATparticles.end(); it++) {
        frictionThermo(*it);
      }
      
    }

    void LangevinThermostat::frictionThermo(Particle& p)
    {
      real massf = sqrt(p.mass());

      // get a random value for each vector component

      Real3D ranval((*rng)() - 0.5, (*rng)() - 0.5, (*rng)() - 0.5);

      p.force() += pref1 * p.velocity() * p.mass() +
                   pref2 * ranval * massf;

      LOG4ESPP_TRACE(theLogger, "new force of p = " << p.force());
    }

    void LangevinThermostat::initialize()
    { // calculate the prefactors

        real timestep = integrator->getTimeStep();

      LOG4ESPP_INFO(theLogger, "init, timestep = " << timestep <<
		    ", gamma = " << gamma << 
		    ", temperature = " << temperature);

      pref1 = -gamma;
      pref2 = sqrt(24.0 * temperature * gamma / timestep);

    }

    /** very nasty: if we recalculate force when leaving/reentering the integrator,
	a(t) and a((t-dt)+dt) are NOT equal in the vv algorithm. The random
	numbers are drawn twice, resulting in a different variance of the random force.
	This is corrected by additional heat when restarting the integrator here.
	Currently only works for the Langevin thermostat, although probably also others
	are affected.
    */

    void LangevinThermostat::heatUp()
    {
      LOG4ESPP_INFO(theLogger, "heatUp");

      pref2buffer = pref2;
      pref2       *= sqrt(3.0);
    }

    /** Opposite to heatUp */

    void LangevinThermostat::coolDown()
    {
      LOG4ESPP_INFO(theLogger, "coolDown");

      pref2 = pref2buffer;
    }

    /****************************************************
    ** REGISTRATION WITH PYTHON
    ****************************************************/

    void LangevinThermostat::registerPython() {


      using namespace espresso::python;


      class_<LangevinThermostat, shared_ptr<LangevinThermostat>, bases<Extension> >
        ("integrator_LangevinThermostat", init<shared_ptr<System> >())
        .def("connect", &LangevinThermostat::connect)
        .def("disconnect", &LangevinThermostat::disconnect)
        .add_property("adress", &LangevinThermostat::getAdress, &LangevinThermostat::setAdress)
        .add_property("gamma", &LangevinThermostat::getGamma, &LangevinThermostat::setGamma)
        .add_property("temperature", &LangevinThermostat::getTemperature, &LangevinThermostat::setTemperature)
        ;


    }

  }
}

