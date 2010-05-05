// ESPP_CLASS
#ifndef _ANALYSIS_PRESSURE_HPP
#define _ANALYSIS_PRESSURE_HPP

#include "types.hpp"
#include "Observable.hpp"

namespace espresso {
  namespace analysis {
    /** Class to compute the pressure. */
    class Pressure : public Observable {
    public:
      Pressure(shared_ptr< System > system) : Observable(system) {}
      ~Pressure() {}
      virtual real compute() const;

      static void registerPython();
    };
  }
}

#endif
