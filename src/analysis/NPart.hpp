// ESPP_CLASS
#ifndef _ANALYSIS_NPART_HPP
#define _ANALYSIS_NPART_HPP

#include "types.hpp"
#include "Observable.hpp"

namespace espresso {
  namespace analysis {
    /** Class to get the number of particles in the system. */
    class NPart : public Observable {
    public:
      NPart(shared_ptr< System > system) : Observable(system) {result_type=real_scalar;}
      virtual ~NPart() {}
      virtual real compute_real() const;

      static void registerPython();
    };
  }
}

#endif
