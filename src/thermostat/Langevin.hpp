#ifndef _LANGEVIN
#define _LANGEVIN

#include "types.hpp"

namespace espresso {
  namespace thermostat {

    class Langevin {
    private:

      real gamma;

    public:

      static void registerPython();

      Langevin();

      Langevin(real _gamma);

      real getGamma() const;

      void setGamma(real _gamma);

      virtual ~Langevin();

    };

  }
}

#endif
