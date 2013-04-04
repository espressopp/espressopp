#include "types.hpp"
#include "python.hpp"
#include "ExtAnalyze.hpp"
#include "SystemAccess.hpp"
#include "analysis/Observable.hpp"

namespace espresso {
  using namespace analysis;
  namespace integrator {

    LOG4ESPP_LOGGER(ExtAnalyze::theLogger, "ExtAnalyze");

    ExtAnalyze::ExtAnalyze(shared_ptr< Observable > _observable, int _interval) : Extension(_observable->getSystem()), interval(_interval)
    {
      LOG4ESPP_INFO(theLogger, "Analyze observable in integrator");
      observable     = _observable;
      reset();
      type = Extension::ExtAnalysis;
    }

    void ExtAnalyze::disconnect(){
      _aftIntV.disconnect();
    }

    void ExtAnalyze::connect(){
      // connection to end of integrator
      _aftIntV  = integrator->aftIntV.connect( boost::bind(&ExtAnalyze::compute, this));
    }

    /** return the averaged observable */
    real ExtAnalyze::getAverage() {
      return obs_ave;
    }

    /** return the standard deviation of the averaged observable*/
    real ExtAnalyze::getVariance() {
      if (n_measurements > 1)
        return obs_var / (n_measurements - 1);
      else
    	return 0.0;
    }

    /** return the number of measurements that have been taken for this observable*/
    int ExtAnalyze::getN() {
      return n_measurements;
    }

    void ExtAnalyze::compute() {
      LOG4ESPP_INFO(theLogger, "computing observable in integrator");
      real res;

      res = observable->compute();
      n_measurements++;

      // compare e. g. Knuth TAOCP vol 2, 3rd edition, page 232
      if (n_measurements == 1) {
          obs_ave     = res;
          obs_ave_old = obs_ave;
          obs_var     = 0.0;
      } else {
          obs_ave = obs_ave_old + (res - obs_ave_old) / n_measurements;
          obs_var = obs_var_old + (res - obs_ave_old) * (res - obs_ave);
          obs_ave_old = obs_ave;
          obs_var_old = obs_var;
      }
    }

    /****************************************************
    ** REGISTRATION WITH PYTHON
    ****************************************************/
    void ExtAnalyze::registerPython() {
      using namespace espresso::python;
      class_<ExtAnalyze, shared_ptr<ExtAnalyze>, bases<Extension> >
        ("integrator_ExtAnalyze", init< shared_ptr< Observable > , int >())
        .def("connect", &ExtAnalyze::connect)
        .def("disconnect", &ExtAnalyze::disconnect)
        .def("getAverage", &ExtAnalyze::getAverage)
        .def("getVariance", &ExtAnalyze::getVariance)
        .def("getN", &ExtAnalyze::getN)
        .def("reset", &ExtAnalyze::reset)
        ;
    }
  }
}
