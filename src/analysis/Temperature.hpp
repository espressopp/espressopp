#ifndef _ANALYSIS_TEMPERATURE_HPP
#define _ANALYSIS_TEMPERATURE_HPP

#include "types.hpp"
#include "Observable.hpp"

namespace espresso {
  namespace analysis {
    /** Class to compute the temperature. */
    class Temperature : public Observable {
    private:
      virtual real compute() const;
      //      virtual real average();
    };
  }
}

#endif
