// ESPP_CLASS
#ifndef _ANALYSIS_XDENSITY_HPP
#define _ANALYSIS_XDENSITY_HPP

#include "types.hpp"
#include "Observable.hpp"

#include "python.hpp"

namespace espresso {
  namespace analysis {
    // Class to compute the density profile along slabs in the x-direction of the system.
    class XDensity : public Observable {
    public:
      XDensity(shared_ptr< System > system) : Observable(system) {}
      ~XDensity() {}
      virtual real compute() const;
      virtual python::list computeArray(int) const;

      static void registerPython();
    };
  }
}

#endif
