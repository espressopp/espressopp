// ESPP_CLASS
#ifndef _ANALYSIS_XPRESSURE_HPP
#define _ANALYSIS_XPRESSURE_HPP

#include "types.hpp"
#include "Observable.hpp"

#include "python.hpp"

namespace espresso {
  namespace analysis {
    // Class to compute the pressure profile along slabs in the x-direction of the system.
    class XPressure : public Observable {
    public:
      XPressure(shared_ptr< System > system) : Observable(system) {}
      ~XPressure() {}
      virtual real compute() const;
      virtual python::list computeArray(int) const;

      static void registerPython();
    };
  }
}

#endif
