// ESPP_CLASS
#ifndef _INTEGRATOR_EXTANALYZE_HPP
#define _INTEGRATOR_EXTANALYZE_HPP

#include "types.hpp"
#include "logging.hpp"
#include "Extension.hpp"
#include "boost/signals2.hpp"
#include "analysis/Observable.hpp"

namespace espresso {
  using namespace analysis;
  namespace integrator {

    /** ExtAnalyze */
    class ExtAnalyze : public Extension {
      public:
        ExtAnalyze(shared_ptr< Observable > _observable, int _interval);
        virtual ~ExtAnalyze() {};
        /** Register this class so it can be used from Python. */
        static void registerPython();

      private:
        boost::signals2::connection _aftIntV;
        void connect();
        void disconnect();
        void compute();

        /** return the averaged observable */
        real getAverage();
        /** return the standard deviation of the averaged observable*/
        real getVariance();
        /** return the number of measurements that have been taken for this observable*/
        int getN();
        /** reset measurement */
        void reset() {
          obs_ave        = 0.0;
          obs_ave_old    = 0.0;
          obs_var        = 0.0;
          obs_var_old    = 0.0;
          n_measurements = 0;
        }

        shared_ptr< Observable > observable;
        int interval;
        real obs_ave, obs_ave_old;
        real obs_var, obs_var_old;
        int n_measurements;

        /** Logger */
        static LOG4ESPP_DECL_LOGGER(theLogger);
    };
  }
}

#endif
