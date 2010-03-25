#ifndef _ANALYSIS_OBSERVABLE_HPP
#define _ANALYSIS_OBSERVABLE_HPP

#include "types.hpp"

namespace espresso {
  namespace analysis {
    /** All quantities to be measured derive from this abstract base class. */
    class Observable {
    public:
      virtual real compute() const = 0;
      virtual real average() = 0;
    };
  }
}

#endif
