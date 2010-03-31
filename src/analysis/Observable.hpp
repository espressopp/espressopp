#ifndef _ANALYSIS_OBSERVABLE_HPP
#define _ANALYSIS_OBSERVABLE_HPP

#include "types.hpp"

// TODO: (1) need to access storage through system
//       (2) compute pressure
//       (3) export to Python

namespace espresso {
  namespace analysis {
    /** All quantities to be measured derive from this abstract base class. */
    class Observable {
    public:
      Observable() {}
      Observable(shared_ptr< storage::Storage > _storage) : storage(_storage) {}
      ~Observable() {}

    public:
      virtual real compute() const = 0;
      //virtual real time_average() const = 0;
      //virtual real spatial_average() const = 0;

      shared_ptr< storage::Storage > storage;
    };
  }
}

#endif
