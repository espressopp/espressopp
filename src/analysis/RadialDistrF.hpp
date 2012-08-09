// ESPP_CLASS
#ifndef _ANALYSIS_RADIALDISTRF_HPP
#define _ANALYSIS_RADIALDISTRF_HPP

#include "types.hpp"
#include "Observable.hpp"

#include "python.hpp"

namespace espresso {
  namespace analysis {
    /** Class to compute the radial distribution function of the system. */
    class RadialDistrF : public Observable {
    public:
      RadialDistrF(shared_ptr< System > system) : Observable(system) {}
      ~RadialDistrF() {}
      virtual real compute() const;
      virtual python::list computeArray(int) const;

      static void registerPython();
    };
  }
}

#endif
