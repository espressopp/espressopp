#ifndef _LANGEVIN
#define _LANGEVIN

#include "types.hpp"
#include "Thermostat.hpp"

namespace espresso {
  namespace thermostat {

    class Langevin: public Thermostat {
    private:

      real gamma;

    public:

      /** Method to register the Python bindings. */
      static void registerPython();

      /** Default constructor. */
      Langevin();

      Langevin(real _gamma);

      Langevin(real _temperature, real _gamma);

      /** Method to get gamma. */
      real getGamma() const;

      /** Method to set gamma. */
      void setGamma(real _gamma);

      virtual ~Langevin();

    };

  }
}

#endif
