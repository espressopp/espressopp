#ifndef _THERMOSTAT_HPP
#define _THERMOSTAT_HPP

#include "types.hpp"
#include "error.hpp"
#include "particles/Set.hpp"

namespace espresso {
  namespace thermostat {

    /** Base class for a thermostat. A thermostat consists at least
        of a temperature and a set of particles to which the thermostat is
        applied.
	
    */
    
    class Thermostat {
    public:
      typedef shared_ptr< Thermostat > SelfPtr;
      
    protected:
      static LOG4ESPP_DECL_LOGGER(theLogger);
      particles::Set::SelfPtr set;  //!< Particles to which the thermostat is applied
      real temperature;   //!< target temperature for the thermostat
      
    public:
      
      /** A constructor of a thermostat needs at least a particle set and a temperature. */
      
      Thermostat(real _temperature) { 
	setTemperature(_temperature); 
      }
      
      /** Setter routine for the temperature
	  
          \param temperature is the new target temperature (must not be negative)
      */
      
      virtual void setTemperature(real _temperature);

      virtual real getTemperature() const { return temperature; }

      virtual void setSet(particles::Set::SelfPtr _set) { set = _set; }

      virtual particles::Set::SelfPtr getSet() { return set; }

      /** Abstract class needs also registration in Python */

      static void registerPython();
    };

  }
}

#endif
