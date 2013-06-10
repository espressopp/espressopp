// ESPP_CLASS
#ifndef _INTEGRATOR_EXTANALYZE_HPP
#define _INTEGRATOR_EXTANALYZE_HPP

#include "types.hpp"
#include "logging.hpp"
#include "Extension.hpp"
#include "boost/signals2.hpp"
#include "analysis/AnalysisBase.hpp"

namespace espresso {
  using namespace analysis;
  namespace integrator {

    /** ExtAnalyze */
    class ExtAnalyze : public Extension {
      public:
        ExtAnalyze(shared_ptr< AnalysisBase > _analysis, int _interval);
        virtual ~ExtAnalyze() {};
        /** Register this class so it can be used from Python. */
        static void registerPython();

      private:
        boost::signals2::connection _aftIntV;
        void connect();
        void disconnect();
        void performMeasurement();

        shared_ptr< AnalysisBase > analysis;
        int interval;
        int counter;

        /** Logger */
        static LOG4ESPP_DECL_LOGGER(theLogger);
    };
  }
}

#endif
