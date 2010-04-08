#ifndef _ANALYSIS_OBSERVABLE_HPP
#define _ANALYSIS_OBSERVABLE_HPP

#include "types.hpp"
#include "SystemAccess.hpp"

// TODO: compute pressure

namespace espresso {
  namespace analysis {
    /** All quantities to be measured derive from this abstract base class. */
    class Observable : public SystemAccess {
    public:
      Observable(shared_ptr< System > system) : SystemAccess(system) {}
      ~Observable() {}

    public:
      virtual real compute() const = 0;
      //virtual real time_average() const = 0;
      //virtual real spatial_average() const = 0;

      static void registerPython();

     protected:
      static LOG4ESPP_DECL_LOGGER(logger);

    };
  }
}

#endif
