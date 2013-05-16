// ESPP_CLASS
#ifndef _ANALYSIS_TEST_HPP
#define _ANALYSIS_TEST_HPP

#include "types.hpp"
#include "AnalysisBase.hpp"

namespace espresso {
  namespace analysis {
    /** Class to test AnalysisBase */
    class Test : public AnalysisBaseTemplate< int > {
    public:
      Test(shared_ptr< System > system) : AnalysisBaseTemplate< int >(system) {}
      int computeRaw();
      static void registerPython();
    };
  }
}

#endif
