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
      static void registerPython();

      Test(shared_ptr< System > system) : AnalysisBaseTemplate< int >(system) {};
      virtual ~Test() {};

      int computeRaw() {
        return 99;
      }

      python::list getInstantValue() {
        python::list ret;
        int res = computeRaw();
        ret.append(res);
        return ret;
      }

      python::list getAverageValue() {
        python::list ret;
        int res;
        res = nMeasurements>0 ? newAverage : 0;
        ret.append(res);
        res = nMeasurements>0 ? newVariance : 0;
        ret.append(res);
        return ret;
      }

      void resetAverage() {
        newAverage = 0.0;
        newVariance = 0.0;
      }

      void updateAverage(int res) {
        // compare e. g. Knuth TAOCP vol 2, 3rd edition, page 232
        if (nMeasurements == 1) {
            newAverage     = res;
            lastAverage    = newAverage;
        } else {
            newAverage  = lastAverage  + (res - lastAverage) / nMeasurements;
            newVariance = lastVariance + (res - lastAverage) * (res - newAverage);
            lastAverage = newAverage;
            lastVariance = newVariance;
        }
        return;
      }
    };
  }
}

#endif
