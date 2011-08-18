// ESPP_CLASS
#ifndef _ANALYSIS_TEMPERATURE_HPP
#define _ANALYSIS_TEMPERATURE_HPP

#include "types.hpp"
#include "Observable.hpp"

namespace espresso {
  namespace analysis {
    /** Class to compute the temperature. */
    class Temperature : public Observable {
    public:
      Temperature(shared_ptr< System > system) : Observable(system) {}
      virtual ~Temperature() {}
      virtual real compute() const;
      //virtual real time_average() const;

      static void registerPython();
    };
  }
}

#endif
