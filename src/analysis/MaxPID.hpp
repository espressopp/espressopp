// ESPP_CLASS
#ifndef _ANALYSIS_MAXPID_HPP
#define _ANALYSIS_MAXPID_HPP

#include "types.hpp"
#include "Observable.hpp"

namespace espresso {
  namespace analysis {
    /** Class to get the number of particles in the system. */
    class MaxPID : public Observable {
    public:
      MaxPID(shared_ptr< System > system) : Observable(system) {}
      virtual ~MaxPID() {}
      virtual real compute() const;

      static void registerPython();
    };
  }
}

#endif
