// ESPP_CLASS
#ifndef _ANALYSIS_CENTEROFMASS_HPP
#define _ANALYSIS_CENTEROFMASS_HPP

#include "types.hpp"
#include "Observable.hpp"
#include "Real3D.hpp"

namespace espresso {
  namespace analysis {
    /** Class to compute the center-of-mass of the system. */
    class CenterOfMass : public Observable {
    public:
      CenterOfMass(shared_ptr< System > system) : Observable(system) {}
      ~CenterOfMass() {}
      virtual real compute() const;
      virtual Real3D computeVector() const;

      static void registerPython();
    };
  }
}

#endif
