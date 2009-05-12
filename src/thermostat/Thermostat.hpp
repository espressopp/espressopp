#ifndef _THERMOSTAT_HPP
#define _THERMOSTAT_HPP

#include "types.hpp"
#include "particles/Set.hpp"

namespace espresso {
  namespace thermostat {

    class Thermostat {

    protected:
      
      boost::shared_ptr<particles::Set> particles;
      real temperature;

    public:

      Thermostat(boost::shared_ptr<particles::Set> _particles, real _temperature):
        particles(_particles), temperature(_temperature) { }

      virtual ~Thermostat() { }
      
      virtual void setTemperature(real _temperature) { temperature = _temperature; }
      
      virtual real getTemperature() const { return temperature; }

      virtual void thermalizeA() = 0;

      virtual void thermalizeB() = 0;
      
      /** Abstract class needs also registration in Python */
      static void registerPython();
    };

  }
}

#endif
