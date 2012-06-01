// ESPP_CLASS
#ifndef _ANALYSIS_NEIGHBORFLUCTUATION_HPP
#define _ANALYSIS_NEIGHBORFLUCTUATION_HPP

#include "types.hpp"
#include "Observable.hpp"

namespace espresso {
  namespace analysis {
    /** Class to get the number of particles in the system. */
    class NeighborFluctuation : public Observable {
    public:
      NeighborFluctuation(shared_ptr< System > system, real _radius) : Observable(system), radius(_radius){}
      virtual ~NeighborFluctuation() {}
      virtual real compute() const;

      static void registerPython();

    protected:
      real radius;
    };
  }
}

#endif
