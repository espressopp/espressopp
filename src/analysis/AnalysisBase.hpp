// ESPP_CLASS
#ifndef _ANALYSIS_ANALYSISBASE_HPP
#define _ANALYSIS_ANALYSISBASE_HPP
#include "types.hpp"
#include "python.hpp"
#include "SystemAccess.hpp"

namespace espresso {
  namespace analysis {

    /** All quantities to be measured derive from this abstract base class and the
    * corresponding template.
    * The following functions have to be implemented by the derived class:
    * virtual ResultType computeRaw()
    * virtual python::list compute()
    * virtual python::list getAverageValue()
    * virtual void resetAverage()
    * virtual void updateAverage(ResultType res) = 0;
    * */
    class AnalysisBase : public SystemAccess {
    public:
      AnalysisBase(shared_ptr< System > system) : SystemAccess(system) {}
      virtual ~AnalysisBase() {}
      virtual void performMeasurement() = 0;
      virtual void reset() = 0;
      virtual python::list compute() = 0;
      virtual python::list getAverageValue() = 0;
      virtual int getNumberOfMeasurements() = 0;
      static void registerPython();
      static LOG4ESPP_DECL_LOGGER(logger);
    };

    template < class ResultType >
    class AnalysisBaseTemplate : public AnalysisBase {
    public:
      AnalysisBaseTemplate(shared_ptr< System > system) : AnalysisBase(system) {
    	  //this->reset();
      };
      virtual ~AnalysisBaseTemplate() {};
      virtual ResultType computeRaw() = 0;
      virtual void performMeasurement();
      virtual void reset();
      virtual void resetAverage() = 0;
      virtual void updateAverage(ResultType res) = 0;
      virtual int getNumberOfMeasurements();
    protected:
      ResultType newAverage, lastAverage;
      ResultType newVariance, lastVariance;
      int nMeasurements;
    };

    /**
    *  These 2 functions have to be implemented by the derived class and
    *  represent the interface to python.
    *  This is an example that only works if ResultType is scalar
    *
    template < class ResultType >
    inline python::list
    AnalysisBaseTemplate< ResultType >::compute() {
        python::list ret;
        ResultType res = computeRaw();
        // BOOST_FOREACH(real value, result_real_vector) ret.append(value);
        ret.append(res);
        return ret;
    }
    template < class ResultType >
    inline python::list
    AnalysisBaseTemplate< ResultType >::getAverageValue() {
        python::list ret;
        ResultType res = newAverage;
        // BOOST_FOREACH(real value, result_real_vector) ret.append(value);
        ret.append(res);
        ResultType res = newVariance;
        // BOOST_FOREACH(real value, result_real_vector) ret.append(value);
        ret.append(res);
        return ret;
    }
    */

    /** This method will measure the observable and can be called from within the integrator
     * by using the ExtAnalyse extension. It can also be called from python.
     * It will also update the average value (typically also the variance) of the observable **/
    template < class ResultType >
    inline void
    AnalysisBaseTemplate< ResultType >::performMeasurement() {
    	ResultType res = computeRaw();
    	nMeasurements++;
    	updateAverage(res);
    }

    /** This method resets/initializes the measurement **/
    template < class ResultType >
    inline void
    AnalysisBaseTemplate< ResultType >::reset() {
    	nMeasurements = 0;
    	resetAverage();
    };

    /**
     *  This function has to be implemented by the derived analysis class
     *  this is an example that only works if ResultType is scalar
     template < class ResultType >
     inline void
     AnalysisBaseTemplate< ResultType >::updateAverage(ResultType res) {
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
    };
    */

    /** return the number of measurements since the last reset **/
    template < class ResultType >
    inline int
    AnalysisBaseTemplate< ResultType >::getNumberOfMeasurements() {
    	return nMeasurements;
    };
  }
}
#endif
