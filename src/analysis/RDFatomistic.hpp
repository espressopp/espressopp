// ESPP_CLASS
#ifndef _ANALYSIS_RDFATOMISTIC_HPP
#define _ANALYSIS_RDFATOMISTIC_HPP

#include "types.hpp"
#include "Observable.hpp"

#include "python.hpp"

namespace espresso {
  namespace analysis {
    /** Class to compute the radial distribution function of the system. */
    class RDFatomistic : public Observable {
    public:
      RDFatomistic(shared_ptr< System > system) : Observable(system) {}
      ~RDFatomistic() {}
      virtual real compute() const;
      virtual python::list computeArray(int) const;

      static void registerPython();
    };
  }
}

#endif
