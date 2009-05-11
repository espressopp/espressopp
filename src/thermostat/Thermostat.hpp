#ifndef _THERMOSTAT_HPP
#define _THERMOSTAT_HPP

#include "types.hpp"
#include "particles/Set.hpp"

namespace espresso {
  namespace thermostat {

    class Thermostat {

    protected:
      
      boost::shared_ptr<particles::Set> particles;
      //Integrator i;
      real temperature;

    public:

      virtual ~Thermostat() {}
      
      virtual void setTemperature(real _temperature) { temperature = _temperature; }
      
      virtual real getTemperature() const { return temperature; }
      
      /** Abstract class needs also registration in Python */
      static void registerPython();
    };

  }
}

#endif
